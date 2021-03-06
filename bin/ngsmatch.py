#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Mar 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: adjusted paths relative to app directory
@last_edit: 2022-02-16
@edit_comment: Adds a lock (PID) file to prevent restarting during execution

Read a NGS Genotyping Stage3.csv file
Read an NGS Genotyping reference file containing assay target sequences.
Process raw FASTQ files and check the sequences they contain against all known target sequences.
Raw MiSeq reads shold be in the 'raw' directory.
It cleans reads into files in the 'clean' directory, then merges them into an 'mclean'
directory ... along with a log of the cleaning and merging that tells us how many
reads there are at each stage.
This code uses BBtools to clean and merge paired-end read files.
Code is now multiprocessing at the "well" level.
Match and mismatch caches are employed to massively speed up matching. ~15x?
Old Bio.pairwise2 algorithm replaced with much faster (~8-9x) Bio.Align

Note logging was removed because it is incompatible with the multiprocessing method employed here.

Needs
"""

from logging import exception
import os
import sys
import argparse
import datetime
import re
import csv
import concurrent.futures
import multiprocessing as mp
from queue import Empty as QueueEmpty
import threading
import collections
from collections import Counter, OrderedDict
import itertools
import subprocess
import gzip
import glob
import json
import time
from functools import partial

import Bio.pairwise2 as pw2
import Bio.SeqIO
import Bio.Align as Align

from bin.echo_primer import padwell, unpadwell
from bin.util import TemplateFiles, output_error

bclen = 8 # length of pseudo barcode
matchcutoff = 1.5 # should control via command line - must be >1.x

basetrans = str.maketrans("ACGT-", "TGCA-") # include '-' to allow for translating matched sequences


class Target(collections.namedtuple("TgtSeqTuple", "hdr seq")):
    """ Target sequence: name/hdr and R1 and R2 versions """
    def __new__(cls, idx, seq):
        """ construct from FASTA record - tuple of header and sequence """
        return super(Target, cls).__new__(cls, idx, seq)

    
    def offmatch(self, s, off):
        """ exact match with target - allow for an initial offset """
        ts = self.seq if not off else self.seq[off:]
        return ts.startswith(s) if len(s)<=len(ts) else s.startswith(ts)


def myopen(fn):
    """ handles gz files seamlessly """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt")
    return open(fn, errors="ignore")


def domatch(mcnt, s, n, cs):
    """
    exact matching for a (merged) read with all known target sequences

    Parameters
    ----------
    mcnt : TYPE
        DESCRIPTION. list counting number of matches to match candidates
    s : TYPE
        DESCRIPTION. sequence to match (merged read)
    n : TYPE
        DESCRIPTION. number of reads that have this (merged) sequence - if there 
        is a match add this to the appropriate value in mcnt
    cs : TYPE
        DESCRIPTION. Match candidates - list of sequences for the assay

    Returns
    -------
    True on success, else False

    """
    try:
        # it would be good to count the offsets below to see what they look like
        for j, cx in enumerate(cs):
            c =cx.seq
            # match with small offset - len(a) <= len(b)
            a, b = s, c
            if len(c)<len(s):
                a, b = c, s
                # note: 'in' operator does exact substring matching
            if len(b)<=len(a)+5 and a in b:
                mcnt[j] += n
                return True
        return False
    except Exception as exc:
        output_error(exc, msg='Error in ngsmatch.domatch')


def exact_match(match_cnt, s, n, targets, mc):
    """
    exact matching for a (merged) read with all known target sequences.
    Length of seq to target must be within 5 bases.

    Parameters
    ----------
    match_cnt : Counter
        DESCRIPTION. list counting number of matches to match candidates
    s : sequence (string)
        DESCRIPTION. sequence to match (merged read)
    n : int
        DESCRIPTION. number of reads that have this (merged) sequence - if there 
        is a match add this to the appropriate value in mcnt
    targets : list(Target)
        DESCRIPTION. Match candidates - list of sequences for the assay

    Returns
    -------
    True on success, else False

    """
    try:
        # it would be good to count the offsets below to see what they look like
        for j, cx in enumerate(targets):
            c = cx.seq
            # match with small offset - len(a) <= len(b)
            a, b = s, c
            if len(c)<len(s):
                a, b = c, s
                # note: 'in' operator does exact substring matching
            if len(b)<=len(a)+5 and a in b:
                match_cnt[cx.hdr] += n
                mc[s] = cx.hdr
                #print('Match cache', match_cache)
                return True
        return False
    except Exception as exc:
        output_error(exc, msg='Error in ngsmatch.exact_match')


def best_match(match_cnt, s, n, targets, mc, sum_parent=True):
    """
    Sequence name - inexact match of sequence to targets - returns the best matching target name
    s is a sequence, ts -s a list (set) or targets to base the name on
    Pairwise match - base name on matching name pluse difference
    
    The name is the most similar sequence with a description of the sequence 
    difference appended ... or just the sequence if there is no sufficiently
    similar sequence.

    Parameters
    ----------
    match_cnt : Counter
    s : string
        DNA sequence.
    n : int 
        counts of this sequence
    targets : list of Targets
        namedtuple with hdr and sequence.
    mc: match_cache
    sum_parent: [T/F] if True add counts to original sequence as well as individual variation

    Returns
    -------
    True if an inexact within 75% was found, False if not.

    """
    try:
        if not s or not targets:
            return False

        if s in mc:
            match_cnt[mc[s]] += n
            return True

        def bestalign(naxs):
            score, best_align, bn = max((x[0].score, x[0], n) for n, x in naxs)
            return bn, best_align
    
        def mkstr2(query_seq, target_seq, start):
            d1 = {i:[cq,ct] for i,(cq,ct) in enumerate(zip(query_seq,target_seq), start=start) if cq!=ct}
            #print(f"mkstr2 {d1}", file=sys.stdout)
            d2 = {}
            prev_i = -5
            for i in d1:
                if prev_i not in d2:
                    d2[i] = d1[i]
                    prev_i = i
                else:
                    if i == prev_i+len(d2[prev_i][0]):  # append change if they are consecutive
                        d2[prev_i][0] += d1[i][0]
                        d2[prev_i][1] += d1[i][1]
                    else:  # otherwise make a new entry
                        d2[i] = d1[i]
                        prev_i = i
            for i in sorted(d2):
                yield str(i+1) + d2[i][1] + '>' + d2[i][0]  # target -> query
            return

        def mkstr(s1, s2, start):
            pos = -10
            for p, (ca, cb) in enumerate(zip(s1, s2), start=start):
                if ca!=cb:
                    # only provide a number for first of a sequence
                    yield ('' if pos+1==p else str(p+1))+cb+'>'+ca
                    pos = p
            return
        varsep = " // "

        def best_new_align(named_aligns):
            best_score, best_query, best_target, best_name = max((algn.score, \
                    str(algn[0]).split('\n')[0], str(algn[0]).split('\n')[2], name) for algn, name in named_aligns)
            #print(best_score, best_target, best_query, best_name, file=sys.stdout)
            return best_name, best_query, best_target, best_score

        use_new_aligner =  True
        if use_new_aligner:
            #print('Attempt inexact matching', file=sys.stdout)
            aligner = Align.PairwiseAligner()
            aligner.mode = 'local'
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -1
            aligner.target_end_gap_score = -1.0
            aligner.query_end_gap_score = -1.0
            aligner.match_score = 1
            aligner.mismatch_score = -1
            all_alignments = [(aligner.align(s, t.seq),t.hdr) for t in targets if varsep not in t.hdr]
            #print(f"Retrieved {len(all_alignments)} alignments with max score {max([a[0].score for a in all_alignments])}", file=sys.stdout)
            best_name, best_query, best_target, best_score = best_new_align(all_alignments)
            if best_score > len(best_query)*0.75:
                #print(f"{os.getpid()} best hit: {best_name} score: {best_score} cutoff: {len(best_query)*0.75} target: {best_target} query: {best_query}", file=sys.stdout)
                left_trimmed_best_query = best_query.lstrip('-')
                start_pos = len(s) - len(left_trimmed_best_query)
                dx = ''.join(mkstr2(best_query, best_target, start_pos))
                #print(f"{os.getpid()} dx: {dx}", file=sys.stdout)
                match_name = varsep.join((best_name, dx))
                #print(f"{os.getpid()} match_name {match_name}", file=sys.stdout)
                if sum_parent:
                    match_cnt[best_name] += n
                match_cnt[match_name] += n
                mc[s] = match_name
            
                #print('Returning True', file=sys.stdout)
                return True
            return False
        else:
            zs = [(t.hdr, pw2.align.localms(s, t.seq, 1, -1, -1, -1,one_alignment_only=True)) for t in targets if varsep not in t.hdr]
            #print(f"{os.getpid()} Bestmatch {len(zs)} {zs}", file=sys.stdout)
            if zs:
            #    print(f"{os.getpid()} Bestmatch, calling bestalign()", file=sys.stdout)
                best_name, ax = bestalign(zs)
            #    print(f"{os.getpid()} Bestmatch, returned from bestalign() {best_name}:{ax}", file=sys.stdout)
                pos = ax.start
                sa, sb = [s[pos:ax.end] for s in (ax.seqA, ax.seqB)]
            #    print(f"{os.getpid()} Bestmatch sa:{sa} \t sb:{sb}", file=sys.stdout)
            #    print(f"{os.getpid()} Bestmatch score: {ax.score} ?>? len(sa)*0.75={len(sa)*0.75}", file=sys.stdout)
                if ax.score>len(sa)*0.75: # only use a name if we have a close-ish (75%) match
                    dx = ''.join(mkstr(sa, sb, pos))
            #        print(f"{os.getpid()} Bestmatch dx:{dx}", file=sys.stdout)
                    match_name = varsep.join((best_name, dx))
                    if sum_parent:
                        match_cnt[best_name] += n
                    match_cnt[match_name] += n
                    mc[s] = match_name
            #        print(f"{os.getpid()} Bestmatch, done... returning True to process_well()", file=sys.stdout)
                    return True
            #    print(f"{os.getpid()} Bestmatch returning False", file=sys.stdout)
            return False
    except Exception as exc:
        output_error(exc, msg='Error in ngsmatch.best_match')


def process_well(work_block, wrs_od_list, targets, match_cache, miss_cache, results, log, lock_c, lock_mc, 
                 lock_r, lock_l, exhaustive, rundir, sum_parent=True, debug=False):
    """ 
    Worker process that runs bbduk, bbmerge, and matching
    work_block: work block number
    wrs_od_list: list of well records ordered dict. Required for multiprocessing
    targets: target sequence dictionary
    match_cache: dict of sequence to sequence identity
    miss_cache: dict of unknown sequences, used as a set
    results: shared list of results + log entries (might separate in future)
    lock_c: lock object for the match cache
    lock_ms: lock object for the miss cache
    lock_r: lock object for the results                        
    exhaustive: [T/F] don't skip any sequences, not matter how low their count
    sum_parent: [T/F] best match adds counts to parent as well as specific variation
    debug: display debugging info
    """
    try:
        varsep = " // "
        PID = os.getpid()
        if debug:
            msg = f"Executing work block {work_block*len(wrs_od_list)} on process {PID}"
            print(msg, file=sys.stdout)
        #cache_hits = 0  # a hit in either the match_cache or miss_cache
        #cache_misses = 0
        result_list = []
        for wr in wrs_od_list:
            #print(f"wr is type: {type(wr)}", file=sys.stdout)
            # find source files
            if debug:
                print(f"{PID} Working on: {wr['pcrplate']} {wr['pcrwell']}", file=sys.stdout)
            rfn= "*{}-{}_*_R1_*.fastq.gz".format(wr['pcrplate'], padwell(wr['pcrwell']))
            fn1s = glob.glob(os.path.join("raw", rfn))
            if not fn1s:
                with lock_l:
                    log.append(f"no data for {fn1s}")
                with lock_r:
                    results.append(tuple(wr.values()))
                return
            
            if len(fn1s)>1:
                with lock_l:
                    log.append(f"too many reads files for pattern{rfn}")
                    log.append(f"   R1 files= {fn1s}")
                with lock_r:
                    results.append(tuple(wr.values()))
                return
            
            fnr1 = fn1s[0]
            fn1 = os.path.basename(os.path.join(rundir,fnr1)) # just the one file
            fn2 = fn1.replace('_R1_', '_R2_')
            fnr2 = os.path.join(rundir, 'raw', fn2) 
            
            # find the data file
            if not os.path.isfile(os.path.join(rundir,fnr2)):
                with lock_l:
                    log.append(f"missing file: {fnr2}")
                with lock_r:
                    results.append(tuple(wr.values()))
                return
        
            # use Windows file separators
            fncs = tuple(os.path.join("clean", fn) for fn in (fn1, fn2))
            bbmapd = os.path.join('bbmap', 'current')
            # unpick the file name a bit
            fnparts = fn1.split('_', 2)
            fnmfmt = "{}-{}_{}{{}}.{{}}".format(wr['pcrplate'], padwell(wr['pcrwell']), fnparts[1])
            fnms = tuple(os.path.join(rundir, "mclean", fnmfmt.format('_'+tx, "fastq.gz")) for tx in ('M', 'U1', 'U2'))
            fnlog = os.path.join(rundir, "mclean", fnmfmt.format('', 'log'))
            if not all(os.path.isfile(os.path.join(rundir,fn)) for fn in (fnms[0], fnlog)):
                # java and bbmap need to be properly installed
                # run BBDuk to clean reads - note: bash isn't available on Windows
                cmd = r"java -ea -Xmx1g -cp {} jgi.BBDuk in1={} in2={} out1={} out2={} t=2 qtrim=rl trimq=20 minlen=50 k=23 ktrim=r mink=11 hdist=1 overwrite=true ref=..\bbmap\resources\adapters.fa".format(bbmapd, fnr1, fnr2, fncs[0], fncs[1])
                pres1 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
            
                # run BBMerge to join paired-end reads
                # cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} verystrict=t".format(*((bbmapd,)+fncs+fnms)).split()
                cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} overwrite=true pfilter=1".format(*((bbmapd,)+fncs+fnms)).split()
                pres2 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
                # could delete cleaned data once it's been merged.
                
                # keep the log output as record counts get used
                with open(fnlog, "wt") as logfile:
                    for s in (pres1.stdout, pres1.stderr):
                        if s:
                            print(s, file=logfile)
                    print('======', file=logfile)
                    for s in (pres2.stdout, pres2.stderr):
                        if s:
                            print(s, file=logfile)
                log1, log2 = ['\n'.join((p.stdout, p.stderr)) for p in (pres1, pres2)]

            else:
                fnm_fmt = fnmfmt.replace('{}.{}','')
                if debug:
                    with lock_l:
                        log.append(f"Skipping cleaning and merging of reads for {fnm_fmt}")
                with open(fnlog) as logfile:
                    logdata = logfile.read()
                log1, log2 = logdata.split("======\n", 1)
            
            # get counts from BBDuk & BBmerge output
            readCount, cleanCount, cleanCount2, joinCount = -1, -1, -1, -1
            m = re.search(r'Input:\s+(\d+) reads', log1)
            if m:
                readCount = int(m.group(1))//2
            m = re.search(r'Result:\s+(\d+) reads', log1) 
            if m:
                cleanCount = int(m.group(1))//2
            m = re.search(r'Pairs:\s+(\d+)\s', log2)
            if m:
                cleanCount2 = int(m.group(1))
            m = re.search(r'Joined:\s+(\d+)\s', log2)
            if m:
                joinCount = int(m.group(1))
                
            if cleanCount!=cleanCount2:
                with lock_l:
                    log.append(f"Pair count mismatch. {cleanCount} != {cleanCount2}")
                with lock_r:
                    results.append(tuple(wr.values()))
                return
                
            fn = fnms[0] # merged reads file
            with myopen(os.path.join(rundir,fn)) as src:
                # read a FASTQ file - sequence in in second line of each 4 line
                gs = (r.strip() for i, r in enumerate(src) if i%4==1)
                seqcnt = collections.Counter(gs)
            mergeCount = sum(seqcnt.values())
            if mergeCount != joinCount:
                with lock_l:
                    log.append(f"Merged counts mismatch. {mergeCount} != {joinCount}") 
                with lock_r:
                    results.append(tuple(wr.values()))
                return
            #assert mergeCount==joinCount # check data is "good"
            
            # we are looking for sequences that match any target in fadict
            # what's the best matching algorithm?
            assayfam = wr['primer'].split('_',1)[0] # assay family name for this well/sample
            #assayfams = [p.split('_',1)[0] for p in wr['primer']] # multiple assay families in same well!!!
            prime_targets = [] #[t for t in targets if any([af.lower() in t[0].lower() for af in assayfams])]
            other_targets = []
            #print(f"Target names: {[t[0] for t in targets]}")
            #print(f"Assay names: {[af for af in assayfams]}")
            for t in targets:
                if assayfam.lower() in t[0].lower():
                    prime_targets.append(t)
                else:
                    other_targets.append(t)
            #print(f"{PID} Prime targets: {len(prime_targets)} {[pt[0] for pt in prime_targets]}")
            #print(f"{PID} Other targets: {len(other_targets)} {[ot[0] for ot in other_targets]}")
            #targets = [fadict[k] for k in sorted(fadict)]
            match_cnt = Counter()
            low_cov_cutoff = max(mergeCount/20, 50)
            other_count = 0
            mc = dict(match_cache)
            original_cache_size = len(mc)
            msc = {}
            if miss_cache is not None:
                msc = dict(miss_cache)
                original_miss_cache_size = len(msc)
                #print(PID, 'Miss cache', type(msc), msc)
            for seq, num in seqcnt.most_common():
                if num >= low_cov_cutoff or exhaustive:
                    if debug:
                        with lock_l:
                            log.append(f"Processing {wr['pcrplate']} {wr['pcrwell']} on process {PID} with counts {num}")
                    if seq in mc:
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Cache hit {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        if sum_parent:
                            if varsep in mc[seq]:
                                match_cnt[mc[seq].split(varsep)[0]] += num
                        #print(f"Cache hit of {num} counts")
                        match_cnt[mc[seq]] += num
                        #cache_hits += 1
                    elif seq in msc:
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Miss cache hit {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        match_cnt[msc[seq]] += num
                        #cache_hits += 1
                    elif exact_match(match_cnt, seq, num, prime_targets, mc):
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Exact match against prime targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        pass
                    elif exact_match(match_cnt, seq, num, other_targets, mc):
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Exact match against other targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        #print(f"{PID}: Exact match of {num} counts")
                        #cache_misses += 1
                        pass
                    elif best_match(match_cnt, seq, num, prime_targets, mc):
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Inexact match against prime targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        pass
                    elif best_match(match_cnt, seq, num, other_targets, mc):
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Inexact match against other_targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        #print(f"Inexact match of {num} counts")
                        #cache_misses += 1
                        pass
                    else:
                        if debug:
                            with lock_l:
                                log.append(f"{PID}: Missed {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        if miss_cache is not None:
                            msc[seq] = seq
                        match_cnt[seq] += num
                        #cache_misses += 1
                else:
                    #if debug:
                    #    print(f"Skipped {wr['pcrplate']} {wr['pcrwell']} on process {PID} with counts {num}", file=sys.stdout)
                    match_cnt['other'] += num
                    other_count += 1
            if len(mc) > original_cache_size:
                #print(f"{PID} Getting cache lock")
                with lock_c:
                    match_cache.update(mc)
                #print(f"{PID} leaving cache lock")
            if miss_cache is not None:
                if len(msc) > original_miss_cache_size:
                    #print(f"{PID} Size of miss cache: {len(msc)}")
                    #print(f"{PID} Getting miss cache lock")
                    with lock_mc:
                        miss_cache.update(msc)
                    #print(f"{PID} Leaving miss cache lock")
            # rename 'other' to include number of separate sequence variations
            match_cnt['other ('+str(other_count)+')'] = match_cnt['other']
            match_cnt['other'] = 0
            res1 = readCount, cleanCount, mergeCount
            if not exhaustive:
                try:
                    resx = sorted(((match_cnt[key], key) for key in match_cnt if match_cnt[key] >= low_cov_cutoff), reverse=True)
                except TypeError:
                    with lock_l:
                        log.append(f"!!! match_cnt: {match_cnt}", file=sys.stdout)
                
            else:
                resx = sorted(((match_cnt[key], key) for key in match_cnt), reverse=True)
            res2 = tuple(x for p in resx for x in p)
            retval = tuple((str(x) for x in tuple(wr.values()) + res1 + res2))
            with lock_r:
                results.append(retval)
            return
        with lock_r:
            results.append(tuple(wr.values()))
        return
    except Exception as exc:
        output_error(exc, msg='Error in ngsmatch.process_well')


def main():
    """
    Read background data: target reference sequences file
    Then processes merged pairs files producing NGS report.
    
    Raises
    ------
    ProgramFail - fail silently...
        if the program fails for any (programmed) reason.

        We catch this and still output our log
    """
    parser = argparse.ArgumentParser(description="NGS Reporting Program")
    parser.add_argument("-d", "--debug", action="store_true", help="more reporting/output for debugging purposes")
    parser.add_argument("-r", "--rundir", required=True, help="Path to experiment/run folder")
    #parser.add_argument("-d", "--datadir", default="mclean", help="data directory")
    #parser.add_argument("-t", "--targets", default=None, help="file of targets in FASTA format")
    # parser.add_argument("-a", "--assays", default=None, help="file linking mice/samples to assay names")
    parser.add_argument('-o','--outfn', default='Results.csv', help='Path to output file name (CSV format)')
    parser.add_argument('-c','--clear_caches',action='store_true', help='Clear match and miss caches in run folder')
    parser.add_argument('-k','--chunk_size', type=int, default=1, help='Number of unique sequences per work unit, default 1')
    parser.add_argument('-n','--ncpus', type=int, default=os.cpu_count(), help='Number of processes to run simultaneously, default=number of CPUs in system')

    parser.add_argument('-x','--exhaustive',action='store_true',help='Try to match every sequence, no matter how few counts')
    parser.add_argument('-q','--quiet',action='store_true',help='Run with minimal messages')
    parser.add_argument('--custom', action='store_true',help='Musterer data columns not present')
    parser.add_argument('-s','--stagefile', default="Stage3.csv", help="NGS genotyping Stage 3 file (default=Stage3.csv")
    args = parser.parse_args()

    print('Run with the following command line options:', file=sys.stdout)
    for arg in vars(args):
        print(arg, getattr(args, arg), file=sys.stdout)
    print('', file=sys.stdout)

    templates = TemplateFiles()

    # read the wells data
    start_time = datetime.datetime.now()
    print(f"Started: {start_time}", file=sys.stdout)
    try:
        with open(os.path.join(args.rundir,args.stagefile)) as srcfd:
            src = csv.reader(srcfd, dialect="unix")
            hdr = next(src)
            WRec = collections.namedtuple("WRec", hdr)
            wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrplate, x.pcrwell[0], int(x.pcrwell[1:]), x.primer))
        if args.debug:
            print(wdata, file=sys.stdout)
        print(f"{len(wdata)} sample wells to process.", file=sys.stdout)
        
        ## get a set of assay family names - family name ends with first underscore char  
        #assays = frozenset(r.primer.split('_',1)[0] for r in wdata)
        
        # read the target sequence file into a dictionary of target sequences
        # make sure there are no lowercase chars in target sequences
         
        # by default use the mouse reference            
        with myopen(os.path.join(args.rundir,templates.files['mouse_ref'])) as src:
            fadict = dict((x.id, Target(x.id, str(x.seq).upper())) for x in Bio.SeqIO.parse(src, "fasta") if len(x.seq) != 0 and len(x.id) != 0)
        # fagen = ((h, s.upper()) for h,s in fasta.fasta_reader(args.targets))
        # fadict = dict(sorted((x.hdr, x) for x in map(Target, fagen)))
        print(f"read {len(fadict)} target sequences.\n", file=sys.stdout)
    
        if not os.path.isdir(os.path.join(args.rundir,"raw")):
            msg = "raw directory is absent - please transfer MiSeq data."
            try:
                raise FileNotFoundError(msg)
            except Exception as exc:
                output_error(exc)
    
        for d in ("clean", "mclean"):
            if not os.path.isdir(os.path.join(args.rundir,d)):
                os.mkdir(os.path.join(args.rundir,d))
                
        # process data for each sample well
        with open(os.path.join(args.rundir,args.outfn), "wt", buffering=1) as dstfd:
            print(f"Opening {args.outfn} for results", file=sys.stdout)
            dst = csv.writer(dstfd, dialect="unix")
            hdrres1 = ("readCount", "cleanCount", "mergeCount")
            hdrres2 = ("seqCount", "seqName")
            complete_row_hdr = tuple((x for xs in (hdr, hdrres1, hdrres2) for x in xs))
            #dst.writerow(x for xs in (hdr, hdrres1, hdrres2) for x in xs)
            dst.writerow(complete_row_hdr)
            #print(','.join(complete_row_hdr), file=dst)
            # process samples wells - code handles bad (manually prepared) plate picklists with >1 assay in a well
            # wdata is already sorted.

            # parallel execute over wdata and collect results. These will then need to be sorted by (x.pcrplate, x.pcrwell)
            grouped_wrs = itertools.groupby(wdata, key=lambda x:(x.pcrplate, x.pcrwell))
            #wrs = list((key,list(group)) for key,group in grouped_wrs)
            wrs = []
            for key, group in grouped_wrs:
                for g in group:
                    wrs.append(OrderedDict(g._asdict()))
        
            with mp.Manager() as manager:
                # Use a server process to manage shared data structures. Heavy option, but works on Windows
                match_cache = manager.dict()
                miss_cache = manager.dict()  # just use this as a set
                results = manager.list()  # at the moment it holds both the main results and the log, but maybe we should split these out?
                log = manager.list()
                lock_mc = manager.Lock()  # match_cache locking
                lock_msc = manager.Lock()  # miss_cache locking
                lock_r = manager.Lock()  # result list locking
                lock_l = manager.Lock()  # log list locking
                #print(type(match_cache))
                if args.clear_caches:
                    if os.path.exists(os.path.join(args.rundir,"match_cache.jsn")):
                        os.remove(os.path.join(args.rundir,"match_cache.jsn"))
                    if os.path.exists(os.path.join(args.rundir,"miss_cache.jsn")):
                        os.remove(os.path.join(args.rundir,"miss_cache.jsn"))
                else:
                    if os.path.exists(os.path.join(args.rundir,"match_cache.jsn")):
                        try:
                            with open(os.path.join(args.rundir,"match_cache.jsn", 'rt')) as json_f:
                                c = json.load(os.path.join(args.rundir,json_f))
                                match_cache.update(c)
                        except:  # it's okay to fail silently, we'll make this file as we close the program
                            print(f"Couldn't open match cache file match_cache.jsn", file=sys.stdout)
                    if os.path.exists(os.path.join(args.rundir,"miss_cache.jsn")):
                        try:
                            with open(os.path.join(args.rundir,"miss_cache.jsn", 'rt')) as json_f:
                                c = json.load(os.path.join(args.rundir,json_f))
                                miss_cache.update(c)
                        except:  # it's okay to fail silently, we'll make this file as we close the program
                            print(f"Couldn't open miss cache file miss_cache.jsn", file=sys.stdout)
        
                targets = [fadict[k] for k in sorted(fadict)]
                # Single process
                #for i,wr in enumerate(wrs):
                #    res = process_well(wr, targets, match_cache)
                #    results[(wr['pcrplate'],wr['pcrwell'])] = res[0]
                #    for entry in res[1]:
                #        logging.info(f"{i} {entry}")

            
                # multiprocessing
                NUMPROCS = args.ncpus
           
                pool = mp.Pool(NUMPROCS)
                chunksize = args.chunk_size
                reports = []
                previous_match_cache_size = len(match_cache)
                previous_miss_cache_size = len(miss_cache)
                for i, chunk in enumerate([wrs[i:i+chunksize] for i in range(0,len(wrs),chunksize)]):
                    while True:
                        reports_waiting = [r for r in reports if not r[0].ready()]
                        if len(reports_waiting) >= NUMPROCS:
                            time.sleep(0.05 + chunksize/20)
                            continue
                        else:
                            break
                    # update caches before we try to launch any more processes
                    with lock_mc:
                        matchc = match_cache
                    with lock_msc:
                        missc = miss_cache
                
                    if missc is not None and not args.quiet:
                        print(f"Adding work to pool... chunk {i*chunksize}, cache size {len(matchc)}, miss cache size {len(missc)}", file=sys.stdout)
                    elif not args.quiet:
                        print(f"Adding work to pool... chunk {i*chunksize}, cache size {len(matchc)}, miss cache disabled", file=sys.stdout)
                    r1 = pool.apply_async(process_well, (i, chunk, targets, matchc, missc, results, log, lock_mc, 
                            lock_msc, lock_r, lock_l, args.exhaustive, args.rundir, args.debug))
                    reports.append((r1,chunk))
                    #if i%4 == 0 and i != 0:
                    #    time.sleep(2 + chunksize/10)

                print('Waiting on all processes to return', file=sys.stdout)
                while True:
                    if any([not r[0].ready() for r in reports]):
                        print('Waiting 1 minute for jobs to complete... you may see many of these messages', file=sys.stdout)
                        time.sleep(60)
                    else:
                        break

                for r in reports:
                    if r[0].ready():
                        if not r[0].successful():
                            print(f"Failure in {r[0]} {r[1]}", file=sys.stdout)
                    try:
                        r[0].get(60)
                        if args.debug:
                            print(f"Final report {r[0]} {r[1]}", file=sys.stdout)
                    except mp.TimeoutError:
                        print(f"Process timeout {r[0]} {r[1]}", file=sys.stdout)
                        r[0].get(1)
            
                print("All processes completed", file=sys.stdout)
                # rewrite cache files
                if args.debug:
                    print(f"Size of match cache {len(match_cache)}", file=sys.stdout)
                with lock_mc:
                    matchc = match_cache
                with open(os.path.join(args.rundir,"match_cache.jsn", 'wt')) as json_f:
                    json_f.write(json.dumps(dict(matchc)))
                
                if args.debug:
                    print(f"Size of miss cache {len(miss_cache)}", file=sys.stdout)
                with lock_msc:
                    missc = miss_cache
                with open(os.path.join(args.rundir,"miss_cache.jsn", 'wt')) as json_f:
                    json_f.write(json.dumps(dict(missc)))

                print('Writing results to' + args.outfn, file=sys.stdout)

                complete_results = {}
                for k,r in enumerate(results):
                    # order results by pcrplate and pcrwell
                    if args.custom:
                        try:
                           complete_results[(r[9],int(r[10][1:]),r[10][0])] = r
                        except IndexError as e:
                           print(e, r, file=sys.stdout)
                           exit(1)
                    else:
                        complete_results[(r[15],int(r[16][1:]),r[16][0])] = r
                
                for i,k in enumerate(sorted(complete_results.keys())):
                    dst.writerow(complete_results[k])
                dstfd.flush()

                pool.close()
                pool.join()
    except Exception as exc:
        output_error(exc, msg='Error in ngsmatch.main')
    end_time = datetime.datetime.now()
    print(f"Finished: {end_time} took: {end_time - start_time}", file=sys.stdout)
            
   
if __name__=="__main__":
    if not os.path.exists(os.path.join(args.rundir,"ngsgeno_lock")):
        with open(os.path.join(args.rundir,"ngsgeno_lock","wt")):
            try:
                main()
            except Exception as exc:
                output_error(exc, msg='Error in ngsmatch')
    else:
        print("Analysis already running", file=sys.stdout)
        exit(2)
    os.remove("ngsgeno_lock")
    
