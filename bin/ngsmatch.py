#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Mar 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021


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

Note: logging module was removed because it is incompatible with the multiprocessing method employed here.
"""

#from ast import Str
from ast import Bytes, Str
import os
from pathlib import Path
from pickletools import bytes1
import sys
import argparse
import datetime
import re
import csv
import difflib
#import concurrent.futures
import multiprocessing as mp
from queue import Empty as QueueEmpty
#import threading
import collections
from collections import Counter, OrderedDict
import itertools
import subprocess
import gzip
import glob
import json
import time
from functools import partial
from math import ceil
try:
    import bin.util as util
except ModuleNotFoundError:
    import util

import Bio.SeqIO
import Bio.Align as Align

# Add/remove leading 0's from well ids - replicated from util.py to allay path issues
def padwell(w):
    "convert a Nimbus well ID to an EP plate (proper/sortable) well ID"
    return '0'.join([w[0], w[1]]) if len(w)==2 else w

def unpadwell(w):
    "Echo software doesn't like well IDs like A01, it wants A1"
    if len(w) == 2:
        return w
    try:
        well = w[0]+w[2] if w[1]=='0' else w
    except:
        print(f"unpadwell failed with {w=}")
    return well

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


def exact_match(match_cnt, query, n, targets, mtc):
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
    mtc: match_cache

    Returns
    -------
    True on success, else False

    """
    for target in targets:
        if len(query) < len(target.seq):
            if query in target.seq:
                match_cnt[target.hdr] += n
                mtc[query] = target.hdr
                return True
        else:
            if target.seq in query:
                match_cnt[target.hdr] += n
                mtc[query] = target.hdr
                return True
    return False


def best_match(match_cnt, query, n, targets, mtc, sum_parent=True):
    """
    Inexact match of sequence to targets - returns the best matching target name
    
    The name is the most similar sequence with at least 75% identity, with a description 
    of the sequence difference appended ... or just the sequence if there is no sufficiently
    similar sequence.

    Parameters
    ----------
    match_cnt : Counter
    query : DNA sequence (string)
    n : counts of this sequence (int)
    targets : list of Targets (namedtuple with hdr and seq)
    mtc: match_cache (dictionary)
    sum_parent: [T/F] if True add counts to original sequence as well as individual variation

    Returns
    -------
    True if an inexact within 75% was found, False if not.

    """
    if not query or not targets:
        return False

    # check if already seen, and add counts
    if query in mtc:
        match_cnt[mtc[query]] += n
        return True

    def compare_seq_strings(best_query, best_target):
        """
        using difflib, make a modification report string
        """
        try:
            #print('In compare_seq_strings', file=sys.stderr)
            left_trimmed_best_query_seq = best_query.lstrip('-')
            start_pos = len(query) - len(left_trimmed_best_query_seq)
            mod_list = [(i+1,li[0],li[2]) for i,li in enumerate(difflib.ndiff(best_target, left_trimmed_best_query_seq)) if li[0] != ' ']
            #print(f'{best_query=} {best_target=} {left_trimmed_best_query_seq=}', file=sys.stderr)
            #print(f'{mod_list=}', file=sys.stderr)
            merged_list = []
            for i,s,c in mod_list:  # integer, sign, character
                if len(merged_list) == 0 and s == '-':
                    merged_list.append([i,c,''])  # [1]: deleted chars, [2]: inserted chars
                elif len(merged_list) == 0 and s == '+':
                    merged_list.append([i,'',c])
                else:
                    if merged_list[-1][0] != i-(len(merged_list[-1][1])+len(merged_list[-1][2])):  # not consecutive difference events
                        if s == '-':
                            merged_list.append([i,c,''])  # deletion (or substitution)
                        elif s == '+':
                            merged_list.append([i,'',c])  # insertion
                    else:  # consecutive with stored sequence
                        if len(merged_list[-1][1]) > 0 and len(merged_list[-1][2]) == 0 and s == '-':  # longer deletion event
                            merged_list[-1][1] += c
                        elif len(merged_list[-1][1]) > 0 and len(merged_list[-1][2]) > 0 and s == '-':  # separate event
                            merged_list.append([i,c,''])
                        elif len(merged_list[-1][1]) > 0 and s == '+':  # substitution
                            merged_list[-1][2] += c
                        elif len(merged_list[-1][1]) == 0 and s == '+':  # insertion
                            merged_list[-1][2] += c
            #print(f'{merged_list=}', file=sys.stderr)     
            final_string_list = []
            for i, del_seq, ins_seq in merged_list:
                if len(del_seq) > 0 and len(ins_seq) > 0:
                    final_string_list.append(f'{i+start_pos}-{del_seq}/+{ins_seq}')
                elif len(del_seq) > 0:
                    final_string_list.append(f'{i+start_pos}-{del_seq}')
                elif len(ins_seq) > 0:
                    final_string_list.append(f'{i+start_pos}+{ins_seq}')

            diff_str = ''.join(final_string_list)
            return diff_str
        except Exception as exc:
            print(f'{exc} failure in compare_seq_strings', file=sys.stderr)
            return 'fail'


    def best_new_align(named_aligns):
        best_score, best_target, best_query, best_name = max([(algn.score, \
                algn.sequences[0], algn.sequences[1], name) for algn, name in named_aligns])
        #print(best_score, best_target, best_query, best_name, file=sys.stderr)
        return best_name, best_target, best_query, best_score

    varsep = " // "
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    aligner.target_end_gap_score = -1.0
    aligner.query_end_gap_score = -1.0
    aligner.match_score = 1
    aligner.mismatch_score = -1
    all_alignments = [(aligner.align(t.seq, query),t.hdr) for t in targets if varsep not in t.hdr]
    #print('alignment', all_alignments[0][0].sequences, all_alignments[0][0].score, file=sys.stderr)
    #print(f"Retrieved {len(all_alignments)} alignments with max score {max([(a[0].score,a[0].sequences[1]) for a in all_alignments])}", file=sys.stderr)
    best_name, best_target, best_query, best_score = best_new_align(all_alignments)
    #print (f'{best_name=} {best_target=} {best_query=} {best_score=}', file=sys.stderr)
    if best_score > len(best_query)*0.75:
        #print(f"{os.getpid()} best hit: {best_name} score: {best_score} cutoff: {len(best_query)*0.75} target: {best_target} query: {best_query}", file=sys.stderr)
        #print('best query', best_query, file=sys.stderr)
        #print('best target', best_target, file=sys.stderr)
        diff_str = compare_seq_strings(best_query, best_target)
        
        #print(f"{os.getpid()} dx: {dx}", file=sys.stderr)
        match_name = varsep.join((best_name, diff_str))
        #print(f"{os.getpid()} match_name {match_name}", file=sys.stderr)
        if sum_parent:
            match_cnt[best_name] += n
        match_cnt[match_name] += n
        mtc[query] = match_name  
        return True
    return False


def preprocess_seqs(wr, rundir, results, log, lock_r, lock_l, debug=False):
    """
    Call bbduk for cleaning of reads
    Call bbmerge to merge read pairs
    If number of merged reads is less than half the number of unmerged reads then call this a failure
    Return seqcnt (Counter of unique sequences) and passing (T/F)
    """
    rfn= "*{}-{}_*_R1_*.fastq.gz".format(util.unguard(wr['pcrplate'],silent=True), padwell(wr['pcrwell']))
    fn1s = glob.glob(os.path.join(rundir, "raw", rfn))
    if not fn1s:
        with lock_l:
            log.append(f"no data for {fn1s}")
            if debug:
                wdb(f"no data for {fn1s}", rundir)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
            
    if len(fn1s)>1:
        with lock_l:
            log.append(f"too many reads files for pattern{rfn}")
            log.append(f"   R1 files= {fn1s}")
            if debug:
                wdb(f"too many reads files for pattern{rfn}", rundir)
                wdb(f"   R1 files= {fn1s}", rundir)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
            
    fnr1 = fn1s[0]
    fn1 = os.path.basename(fnr1) # just the one file
    fn2 = fn1.replace('_R1_', '_R2_')
    fnr2 = os.path.join(rundir, 'raw', fn2) 
            
    # find the data file
    if not os.path.isfile(fnr2):
        with lock_l:
            log.append(f"missing file: {fnr2}")
            if debug:
                wdb(f"missing file: {fnr2}", rundir)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
        
    # use Windows file separators
    fncs = tuple(os.path.join(rundir,"cleaned", fn) for fn in (fn1, fn2))
    bbmapd = os.path.join('bbmap','current')
    # unpick the file name a bit
    fnparts = fn1.split('_', 2)
    fnmfmt = "{}-{}_{}{{}}.{{}}".format(wr['pcrplate'], padwell(wr['pcrwell']), fnparts[1])
    fnms = tuple(os.path.join(rundir,"merged", fnmfmt.format('_'+tx, "fastq.gz")) for tx in ('M', 'U1', 'U2'))
    fnlog = os.path.join(rundir,"merged", fnmfmt.format('', 'log'))
    if not all(os.path.isfile(fn) for fn in (fnms[0], fnlog)):
        # java and bbmap need to be properly installed
        cmd = r"java -ea -Xmx1g -cp {} jgi.BBDuk in1={} in2={} out1={} out2={} t=2 qtrim=rl trimq=20 minlen=50 k=23 ktrim=r mink=11 hdist=1 overwrite=true ref=bbmap\resources\adapters.fa".format(bbmapd, fnr1, fnr2, fncs[0], fncs[1])
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
        with lock_l:
            log.append(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}")
            if debug:
                wdb(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}")
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
            msg = f"Pair count mismatch. {cleanCount} != {cleanCount2}"
            log.append(msg)
            if debug:
                wdb(msg, rundir)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
                
    fn = fnms[0] # merged reads file
    with myopen(fn) as src:
        # read a FASTQ file - sequence in in second line of each 4 line
        gs = (r.strip() for i, r in enumerate(src) if i%4==1)
        seqcnt = collections.Counter(gs)
    mergeCount = sum(seqcnt.values())
    if mergeCount != joinCount:
        with lock_l:
            msg = f"Merged counts mismatch. {mergeCount} != {joinCount}"
            log.append(msg)
            if debug:
                wdb(msg, rundir)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
    #assert mergeCount==joinCount # check data is "good"
    return seqcnt, True, readCount, cleanCount, mergeCount

def wdb(msg, rundir):
    """ wdb = write debug. Call from within a lock """
    debugfn = os.path.join(rundir, 'debug.log')
    with open(debugfn, 'at') as df:
        print(msg, file=df)

def process_well(work_block, wrs_od_list, rundir, targets, match_cache, miss_cache, results, log, lock_c, lock_mc, 
                 lock_r, lock_l, mincov, minprop, exhaustive, sum_parent=True, debug=False):
    """ 
    Worker process that runs bbduk, bbmerge, and matching
    work_block: work block number
    wrs_od_list: list of well records ordered dict. Required for multiprocessing
    targets: list of Target(hdr, seq)
    match_cache: dict of sequence to sequence identity
    miss_cache: dict of unknown sequences, used as a set
    results: shared list of results + log entries (might separate in future)
    lock_c: lock object for the match cache
    lock_ms: lock object for the miss cache
    lock_r: lock object for the results                        
    exhaustive: [T/F] don't skip any sequences, not matter how low their count
    sum_parent: [T/F] best match adds counts to parent as well as specific variation
    debug: write messages direct to file (slow I/O)
    """
    varsep = " // "
    PID = os.getpid()
    msg = f"Debug: Executing work block {work_block*len(wrs_od_list)} on process {PID}"
    print(msg, file=sys.stderr)
    with lock_l:
        log.append(msg)
        if debug:
            wdb(msg, rundir)

    for wr in wrs_od_list: # work record in chunk
        if debug:
            with lock_l:
                msg = f"Debug:{PID} Working on: {wr['pcrplate']} {wr['pcrwell']}"
                log.append(msg)
                wdb(msg, rundir)

        # clean and merge FASTQs
        seqcnt, success, readCount, cleanCount, mergeCount = preprocess_seqs(wr, rundir, results, log, lock_r, lock_l, debug=debug)
        #print(f'{log=} {success=} {wr=}', file=sys.stderr)
        if not success:
            continue        
            
        # we are looking for sequences that match any target in fadict
        # what's the best matching algorithm?
        assayfam = wr['primer'].split('_',1)[0] # assay family name for this well/sample
        #assayfams = [p.split('_',1)[0] for p in wr['primer']] # multiple assay families in same well!!!
        prime_targets = [] #[t for t in targets if any([af.lower() in t[0].lower() for af in assayfams])]
        other_targets = []
        #print(f"Target names: {[t[0] for t in targets]}")
        #print(f"Assay names: {[af for af in assayfams]}")
        #print(f'{targets=}', file=sys.stderr)
        for t in targets:
            if assayfam.lower() in t.hdr.lower():
                prime_targets.append(t)
            else:
                other_targets.append(t)
        #print(f"{PID} Prime targets: {len(prime_targets)} {[pt[0] for pt in prime_targets]}")
        #print(f"{PID} Other targets: {len(other_targets)} {[ot[0] for ot in other_targets]}")
        #targets = [fadict[k] for k in sorted(fadict)]
        match_cnt = Counter()
        if mincov < 1:
            mincov = 1
        if minprop > 1.0:
            minprop = 1.0
        elif minprop < 0.0:
            minprop = 0.0
        low_cov_cutoff = max(mergeCount*minprop, mincov)
        other_count = 0
        mc = dict(match_cache)
        original_match_cache_size = len(mc)
        msc = dict(miss_cache)
        original_miss_cache_size = len(msc)
            #print(PID, 'Miss cache', type(msc), msc)
        for seq, num in seqcnt.most_common():
            if num >= low_cov_cutoff or exhaustive:
                if debug:
                    with lock_l:
                        msg = f"Processing {wr['pcrplate']} {wr['pcrwell']} on process {PID} with counts {num}"
                        log.append(msg)
                        wdb(msg, rundir)
                if seq in mc:
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Cache hit {wr['pcrplate']} {wr['pcrwell']} {num} {seq}"
                            log.append(msg)
                            wdb(msg, rundir)
                    if sum_parent:
                        if varsep in mc[seq]:
                            match_cnt[mc[seq].split(varsep)[0]] += num
                    #print(f"Cache hit of {num} counts")
                    match_cnt[mc[seq]] += num
                    #cache_hits += 1
                elif seq in msc:
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Miss cache hit {wr['pcrplate']} {wr['pcrwell']} {num} {seq}"
                            log.append(msg)
                            wdb(msg, rundir)
                    match_cnt[msc[seq]] += num
                    #cache_hits += 1
                elif exact_match(match_cnt, seq, num, prime_targets, mc):
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Exact match against prime targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}"
                            log.append(msg)
                            wdb(msg, rundir)
                    pass
                elif exact_match(match_cnt, seq, num, other_targets, mc):
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Exact match against other targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}"
                            log.append(msg)
                            wdb(msg, rundir)
                    #print(f"{PID}: Exact match of {num} counts")
                    #cache_misses += 1
                    pass
                elif best_match(match_cnt, seq, num, prime_targets, mc):
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Inexact match against prime targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}"
                            log.append(msg)
                            wdb(msg, rundir)
                    pass
                elif best_match(match_cnt, seq, num, other_targets, mc):
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Inexact match against other_targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}" 
                            log.append(msg)
                            wdb(msg, rundir)
                    #print(f"Inexact match of {num} counts")
                    #cache_misses += 1
                    pass
                else:
                    if debug:
                        with lock_l:
                            msg = f"{PID}: Missed {wr['pcrplate']} {wr['pcrwell']} {num} {seq}"
                            log.append(msg)
                            wdb(msg, rundir)
                    if miss_cache is not None:
                        msc[seq] = seq
                    match_cnt[seq] += num
            else:
                match_cnt['other'] += num
                other_count += 1
        if len(mc) > original_match_cache_size:
            #print(f"{PID} Getting cache lock")
            with lock_c:
                match_cache.update(mc) 
        if len(msc) > original_miss_cache_size:
            with lock_mc:
                miss_cache.update(msc)
        # rename 'other' to include number of separate sequence variations
        match_cnt['other ('+str(other_count)+')'] = match_cnt['other']
        match_cnt['other'] = 0
        res1 = [readCount, cleanCount, mergeCount]
        if not exhaustive:
            resx = sorted(((match_cnt[key], key) for key in match_cnt if match_cnt[key] >= low_cov_cutoff), reverse=True)
        else:
            resx = sorted(((match_cnt[key], key) for key in match_cnt), reverse=True)
        seqCounts = []
        seqNames = []                                        
        otherCounts = []
        otherNames = []
        try:
            for count, key in resx:
                if varsep not in key and assayfam.lower() in key.lower():
                    seqCounts.append(str(count))
                    seqNames.append(str(key))
                else:
                    otherCounts.append(str(count))
                    otherNames.append(str(key))
        except Exception as exc:
            msg = f'Exception {exc} for {count=} {key=} in result {resx=}'
            print(msg, file=sys.stderr) 
            wdb(msg, rundir)
        #res2 = tuple(x for p in resx for x in p)
        try:
            res2 = [';'.join(seqCounts), ';'.join(seqNames), ';'.join(otherCounts), ';'.join(otherNames)]
        except Exception as exc:
            msg = f'Exception {exc} in joining seqCounts and seqNames {seqCounts=} {seqNames=} {otherCounts=} {otherNames=}'
            print(msg, file=sys.stderr)
            wdb(msg, rundir)
        #print(res2, file=sys.stderr)
        try:
            retval = [str(wr[x]) for x in wr] + res1 + res2
            with lock_r:
                results.append(retval)
        except Exception as exc:
            msg = f'Exception {exc} in forming return value from matching {wr=} {res1=} {res2=}'
            print(msg, file=sys.stderr)       
            wdb(msg, rundir)


def write_log(log, logfn):
    print(f"Writing log to {logfn}", file=sys.stderr)
    with open(logfn, 'wt') as logf:
        for l in log:
            print(l, file=logf)

 
def get_raw_fastq_pairs(dirpath):  
    """ return a sorted list of tuple(R1_path, R2_path) to raw FASTQ files """
    valid_pairs = []
    r1s = [dirpath/Path(f) for f in os.listdir(dirpath) if f.endswith('.fastq.gz') and '_R1_' in f]
    for r1 in r1s:
        r2 = Path(str(r1).replace('_R1_','_R2_'))
        if r1.is_file() and r2.is_file():
            valid_pairs.append((r1,r2))
    return sorted(valid_pairs)


def report_progress(rundir, launch_progress, match_progress):
    """
    Clear all previous progress files and touch a new file with the launch and match progress values in the name
    Only do the operation is progress is a multiple of 5 to save on disk writes
    """
    progress_files = list(Path(rundir).glob('match_progress_*'))
    if len(progress_files) == 1:
        lp = int(str(progress_files[0]).split('_')[-2])
        mp = int(str(progress_files[0]).split('_')[-1])
        if lp == launch_progress and mp == match_progress:
            # no change
            return
    for fn in list(Path(rundir).glob('match_progress_*')):
        os.remove(fn)
    progress_fn = os.path.join(args.rundir,'match_progress_'+str(launch_progress)+'_'+str(match_progress))
    with open(progress_fn, 'wt') as fp:
        pass  # just touch the file


def main(args):
    """
    Read background data: target reference sequences file
    Then processes merged pairs files producing NGS report.
    """
    log = []
    if args.debug:
        db_log_fn = os.path.join(args.rundir, 'debug.log')
        if Path(db_log_fn).exists():
            try:
                os.remove(db_log_fn)
            except Exception as exc:
                print('Could not clear debug.log, perhaps you have it open?', file=sys.stderr)
                return
        print(f'Writing debug information to {db_log_fn}', file=sys.stderr)
    
    log.append('Info: Run with the following command line options:')
    for arg in vars(args):
        log.append(f'Info: {arg} {getattr(args, arg)}')
                                                                                                                                            
    # read the wells data
    start_time = datetime.datetime.now()
    log.append(f"Begin: {start_time}")
    raw_pair_list = get_raw_fastq_pairs(os.path.join(args.rundir, 'raw'))
    raw_file_identifiers = set([str(f[0].name).split('_')[0] for f in raw_pair_list])
    #print(sorted(raw_file_identifiers), file=sys.stderr)
    with open(os.path.join(args.rundir,args.stagefile)) as srcfd:
        src = csv.reader(srcfd, dialect="unix")
        hdr = next(src)
        WRec = collections.namedtuple("WRec", hdr)
        wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrplate, x.pcrwell[0], int(x.pcrwell[1:]), x.primer))
    #print(wdata, file=sys.stderr)
    wdata = [rec for rec in wdata if util.unguard(rec.pcrplate, silent=True) +'-'+ util.padwell(rec.pcrwell) in raw_file_identifiers]
    #print(f'After filtering by available files {wdata=}', file=sys.stderr)
    log.append(f"Info: {len(wdata)} sample wells to process.")
    #print(log, file=sys.stderr)
    ## get a set of assay family names - family name ends with first underscore char  
    #assays = frozenset(r.primer.split('_',1)[0] for r in wdata)
        
    # read the target sequence file into a dictionary of target sequences
                    
    fadict = {}
    with myopen(os.path.join(args.rundir,args.targets)) as src:
        for id_seq in Bio.SeqIO.parse(src, "fasta"):
            if len(id_seq.id) > 0 and len(id_seq.seq) > 0:
                try:
                    seq_id = bytes(str(id_seq.id).strip(), 'ascii', errors='strict').decode()
                except Exception as exc:
                    msg = f'Warning: {exc} non-ascii character in reference sequence id {id_seq.id}'
                    log.append(msg)
                    if args.debug:
                        wdb(msg, args.rundir)
                    seq_id = bytes(str(id_seq.id).strip(), 'ascii',errors='ignore').decode()

                try:
                    seq = bytes(str(id_seq.seq).strip(), 'ascii', errors='strict').decode().upper()
                except Exception as exc:
                    msg = f'Warning: {exc} non-ascii character in reference sequence {id_seq.seq}'
                    log.append(msg)
                    if args.debug:
                        wdb(msg, args.rundir)
                    seq = bytes(str(id_seq.seq).strip(), 'ascii', errors='ignore').decode().upper()

                if seq_id in fadict:
                    msg = f"Warning: sequence for {seq_id} already appears in reference target list "+\
                            f"with sequence {fadict[seq_id]=}, alternative {seq=}"
                    log.append(msg)
                    if args.debug:
                        wdb(msg, args.rundir)

                fadict[seq_id] = Target(seq_id, seq)
        #fadict = dict((x.id, Target(x.id, str(x.seq).upper())) for x in Bio.SeqIO.parse(src, "fasta") if len(x.seq) != 0 and len(x.id) != 0)
    
    log.append(f"Info: read {len(fadict)} target sequences.\n")
    
    if not os.path.isdir(os.path.join(args.rundir,"raw")):
        log.append("Error: raw FASTQ data folder is absent - please transfer MiSeq data.")
        write_log(log, os.path.join(args.rundir,args.logfn))
        return
    
    for d in ("cleaned", "merged"):
        dp = os.path.join(args.rundir, d)
        if not os.path.isdir(dp):
            os.mkdir(dp)
                
    # process data for each sample well
    
    if True:
        # process samples wells - code handles bad (manually prepared) plate picklists with >1 assay in a well
        # wdata is already sorted.

        # parallel execute over wdata and collect results
        grouped_wrs = itertools.groupby(wdata, key=lambda x:(x.pcrplate, x.pcrwell))
        wrs = []
        for key, group in grouped_wrs:
            for g in group:
                wrs.append(OrderedDict(g._asdict()))
        
        with mp.Manager() as manager:
            # Use a server process to manage shared data structures. Heavy option, but works on Windows
            match_cache = manager.dict()
            miss_cache = manager.dict()  # just use this as a set
            results = manager.list()  # at the moment it holds both the main results and the log, but maybe we should split these out?
            logm = manager.list()
            lock_mtc = manager.Lock()  # match_cache locking
            lock_msc = manager.Lock()  # miss_cache locking
            lock_r = manager.Lock()  # result list locking
            lock_l = manager.Lock()  # log list locking
      
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
            launch_progress = 0
            match_progress = 0
            for i, chunk in enumerate([wrs[i:i+chunksize] for i in range(0,len(wrs),chunksize)]):
                while True:
                    reports_waiting = [r for r in reports if not r[0].ready()]
                    if len(reports_waiting) >= NUMPROCS:
                        time.sleep(0.05 + chunksize/40)
                        launch_progress = int(100*(i+1)/ceil(len(wrs)/chunksize))
                        match_progress = int(100*len([r for r in reports if r[0].ready()])/(len(wrs)/chunksize))
                        report_progress(args.rundir, launch_progress, match_progress)
                        continue
                    else:
                        # update caches before we try to launch any more processes
                        with lock_mtc:
                            matchc = match_cache
                        with lock_msc:
                            missc = miss_cache
                        msg = f"Adding work to pool... chunk {i} of size {chunksize}, cache size {len(matchc)}, "+\
                                f"miss cache size {len(missc)}, outstanding reports {len(reports_waiting)}"
                        print(msg, file=sys.stderr)
                        log.append(msg)
                        if args.debug:
                            wdb(msg, args.rundir)                          
                
                        r1 = pool.apply_async(process_well, (i, chunk, args.rundir, targets, matchc, missc, results, logm, lock_mtc, 
                                lock_msc, lock_r, lock_l, args.mincov, args.minprop, args.exhaustive, True, args.debug))
                        reports.append((r1,chunk))
                        launch_progress = int(100*(i+1)/ceil(len(wrs)/chunksize))
                        match_progress = int(100*len([r for r in reports if r[0].ready()])/(len(wrs)/chunksize))
                        report_progress(args.rundir, launch_progress, match_progress)
                        break
                
                
                #if i%4 == 0 and i != 0:
                #    time.sleep(2 + chunksize/10)

            print('Waiting on all processes to return', file=sys.stderr)
            while True:
                if any([not r[0].ready() for r in reports]):
                    match_progress = int(100*len([r for r in reports if r[0].ready()])/(len(wrs)/chunksize))
                    report_progress(args.rundir, launch_progress, match_progress)
                    print('Waiting 20 seconds for jobs to complete... you may see many of these messages', file=sys.stderr)
                    time.sleep(20)
                else:
                    match_progress = int(100*len([r for r in reports if r[0].ready()])/(len(wrs)/chunksize))
                    report_progress(args.rundir, launch_progress, match_progress)
                    break

            for r in reports:
                if r[0].ready():
                    if not r[0].successful():
                        log.append(f"Error: Failure in {r[0]} {r[1]}")
                try:
                    r[0].get(60)
                    if args.debug:
                        log.append(f"Debug: Final report {r[0]} {r[1]}")
                except mp.TimeoutError:
                    print(f"Process timeout {r[0]} {r[1]}", file=sys.stderr)
                    r[0].get(1)
            
            log.append("Info: All processes completed")
            with lock_l:
                for l in logm:
                    log.append(l)

            #print(f'{list(results)=}', file=sys.stderr)
            log.append(f'Info: Writing results to {args.outfn}')
            #print('hello0.4', file=sys.stderr)
            complete_results = {}  # order results by sampleNo
            for k,r in enumerate(results):
                complete_results[int(r[0])] = r

        with open(os.path.join(args.rundir,args.outfn), "wt", buffering=1) as dstfd:
            print(f"Opening {args.outfn} for results", file=sys.stderr)
            dst = csv.writer(dstfd, dialect="unix", quoting=csv.QUOTE_ALL)
            hdrres1 = ("readCount", "cleanCount", "mergeCount")
            hdrres2 = ("seqCount", "seqName", "otherCount", "otherName")
            complete_row_hdr = tuple((x for xs in (hdr, hdrres1, hdrres2) for x in xs))
            dst.writerow(complete_row_hdr)
            for k in sorted(complete_results):
                if args.debug:
                    print(f'{complete_results[k]=}', file=sys.stderr)
                dst.writerow(complete_results[k])
            #for i,k in enumerate(sorted(complete_results.keys())):
            #    dst.writerow(complete_results[k])
            dstfd.flush()
            #print('hello0.5', file=sys.stderr)
            pool.close()
            #print('hello0.7', file=sys.stderr)
            pool.join()
            #print('hello1', file=sys.stderr)
        #print('hello2', file=sys.stderr)
    #print('hello3', file=sys.stderr)
    end_time = datetime.datetime.now()
    log.append(f"End: {end_time} took: {end_time - start_time}")
    write_log(log, os.path.join(args.rundir,args.logfn))
            
   
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="NGS Reporting Program")
    parser.add_argument("-d", "--debug", action="store_true", help="more reporting/output for debugging purposes")
    parser.add_argument("-t", "--targets", default="targets.fa", help="file of targets in FASTA format (default=targets.fa)")
    # parser.add_argument("-a", "--assays", default=None, help="file linking mice/samples to assay names")
    parser.add_argument('-o','--outfn', default='results.csv', help='Name of output file (CSV format)')
    parser.add_argument('-k','--chunk_size', type=int, default=1, help='Number of unique sequences per work unit, (default=1)')
    parser.add_argument('-n','--ncpus', type=int, default=os.cpu_count()-1, help='Number of processes to run simultaneously, default=number of CPUs in system - 1')
    parser.add_argument('-l','--logfn', default='match.log', help='Name of logging file (default=match.log)')
    parser.add_argument('-r','--rundir', required=True, help='Path to experiment folder')
    parser.add_argument('-m','--mincov', type=int, default=50, help='Do not match unique sequences with less than this many reads coverage, default 50')
    parser.add_argument('-p','--minprop', type=float, default=0.2, help='Do not match unique sequences '+\
            'with less than this proportion of the total number of reads, default 0.2. Must be between 0.0 and 1.0')
    parser.add_argument('-x','--exhaustive',action='store_true',help='Try to match every sequence, '+\
            'no matter how few counts. Ignores --minseqs and --minprop')
    #parser.add_argument('--custom', action='store_true',help='Musterer data columns not present')
    parser.add_argument('-s','--stagefile', default="Stage3.csv", help="Name of the NGS genotyping Stage 3 file (default=Stage3.csv)")
    args = parser.parse_args()
    in_error = False
    print(f"{args=}", file=sys.stderr)
    lock_path = os.path.join(args.rundir,"ngsgeno_lock") 
    if not os.path.exists(lock_path):
        try:
            with open(lock_path,"wt"):
                report_progress(args.rundir, 0, 0)  # set this up asap
                main(args)
                os.remove(lock_path)
        except:
            os.remove(lock_path)
    else:
        print("Analysis already running", file=sys.stderr)
        exit(2)
    
    
    
