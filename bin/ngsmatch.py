#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Mar 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: Major overhaul to take Experiment object and work with the GUI
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
import re
import csv
import concurrent.futures
import multiprocessing as mp
from queue import Empty as QueueEmpty
import collections
from collections import Counter, OrderedDict
import itertools
import subprocess
import gzip
import glob
import time
from functools import partial

import Bio.pairwise2 as pw2
import Bio.SeqIO
import Bio.Align as Align

from bin.util import padwell, unpadwell

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


#def domatch(mcnt, s, n, cs):
#    """
#    DEPRECATED
#    exact matching for a (merged) read with all known target sequences

#    Parameters
#    ----------
#    mcnt : TYPE
#        DESCRIPTION. list counting number of matches to match candidates
#    s : TYPE
#        DESCRIPTION. sequence to match (merged read)
#    n : TYPE
#        DESCRIPTION. number of reads that have this (merged) sequence - if there 
#        is a match add this to the appropriate value in mcnt
#    cs : TYPE
#        DESCRIPTION. Match candidates - list of sequences for the assay

#    Returns
#    -------
#    True on success, else False

#    """
#    try:
#        # it would be good to count the offsets below to see what they look like
#        for j, cx in enumerate(cs):
#            c =cx.seq
#            # match with small offset - len(a) <= len(b)
#            a, b = s, c
#            if len(c)<len(s):
#                a, b = c, s
#                # note: 'in' operator does exact substring matching
#            if len(b)<=len(a)+5 and a in b:
#                mcnt[j] += n
#                return True
#        return False
#    except Exception as exc:
#        output_error(exc, msg='Error in ngsmatch.domatch')


def exact_match(match_cnt, s, n, targets, mtc):
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
    mtc: match cache

    Returns
    -------
    True on success, else False

    """
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
            mtc[s] = cx.hdr
            return True
    return False


def best_match(match_cnt, s, n, targets, mtc, sum_parent=True):
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
    mtc: match_cache
    sum_parent: [T/F] if True add counts to original sequence as well as individual variation

    Returns
    -------
    True if an inexact match within 75% was found, False if not.

    """
   
    if not s or not targets:
        return False

    if s in mtc:
        match_cnt[mtc[s]] += n
        return True
    
    def mkstr2(query_seq, target_seq, start):
        """ combine the query sequence, target sequence and starting position into one string to indicate what has been mutated """
        d1 = {i:[cq,ct] for i,(cq,ct) in enumerate(zip(query_seq,target_seq), start=start) if cq!=ct}
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

    varsep = " // "

    def best_new_align(named_aligns):
        best_score, best_query, best_target, best_name = max((algn.score, \
                str(algn[0]).split('\n')[0], str(algn[0]).split('\n')[2], name) for algn, name in named_aligns)
        #print(best_score, best_target, best_query, best_name, file=sys.stdout)
        return best_name, best_query, best_target, best_score
      
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
        mtc[s] = match_name
            
        #print('Returning True', file=sys.stdout)
        return True
    return False
        


def process_well(work_block, wrs_od_list, targets, match_cache, miss_cache, results, logs, lock_mtc, lock_msc, 
                 lock_r, lock_l, exhaustive, datadir, cleandir, mergeddir, sum_parent=True, debug=False):
    """ 
    Worker process that runs bbduk, bbmerge, and matching
    work_block: work block number
    wrs_od_list: list of well records ordered dict. Required for multiprocessing
    targets: target sequence dictionary
    match_cache: shared dict of sequence to sequence identity
    miss_cache: shared dict of unknown sequences, used as a set
    results: shared lockable list of results 
    logs: shared lockable list of log entries
    lock_mtc: lock object for the match cache
    lock_msc: lock object for the miss cache
    lock_r: lock object for the results
    lock_l: lock object for the log
    exhaustive: [T/F] don't skip any sequences, not matter how low their count
    datadir: dirpath for raw fastq files
    cleandir: dirpath for cleaned fastq files
    mergeddir: dirpath for merged fastq paired ends
    sum_parent: [T/F] best match adds counts to parent as well as specific variation
    debug: display debugging info
    """
    log_entries = []
    try:
        varsep = " // "
        PID = os.getpid()
        log_entries.append(f"Info: Executing work block {work_block*len(wrs_od_list)} on process {PID}")          
        #cache_hits = 0  # a hit in either the match_cache or miss_cache
        #cache_misses = 0
        
        for wr in wrs_od_list:
            # find source files
            log_entries.append(f"Info: {PID} Working on: {wr['pcrplate']} {wr['pcrwell']}")
            rfn= "*{}-{}_*_R1_*.fastq.gz".format(wr['pcrplate'], padwell(wr['pcrwell']))
            fn1s = glob.glob(os.path.join(datadir, rfn))
            if not fn1s:
                log_entries.append(f"Warning: no data for {fn1s}")
                with lock_l:
                    for l in log_entries:
                        logs.append(l)
                with lock_r:
                    results.append(tuple(wr.values()))
                return
            
            if len(fn1s)>1:
                log_entries.append(f"Error: too many reads files for pattern{rfn}")
                log_entries.append(f"Info: R1 files {fn1s=}")
                with lock_l:
                    for l in log_entries:
                        logs.append(l)
                with lock_r:
                    results.append(tuple(wr.values()))
                return
            
            fnr1 = fn1s[0]
            fn1 = os.path.basename(os.path.join(datadir,fnr1)) # just the one file
            fn2 = fn1.replace('_R1_', '_R2_')
            fnr2 = os.path.join(datadir, 'raw', fn2) 
            
            # find the data file
            if not os.path.isfile(os.path.join(datadir,fnr2)):
                log_entries.append(f"Error: missing file: {fnr2}")
                with lock_l:
                    for l in log_entries:
                        logs.append(l)
                with lock_r:
                    results.append(tuple(wr.values()))
                return
        
            # use Windows file separators
            fncs = tuple(os.path.join(cleandir, fn) for fn in (fn1, fn2))
            bbmapd = os.path.join('bbmap', 'current')
            # unpick the file name a bit
            fnparts = fn1.split('_', 2)
            fnmfmt = "{}-{}_{}{{}}.{{}}".format(wr['pcrplate'], padwell(wr['pcrwell']), fnparts[1])
            fnms = tuple(os.path.join(mergeddir, fnmfmt.format('_'+tx, "fastq.gz")) for tx in ('M', 'U1', 'U2'))
            fnlog = os.path.join(mergeddir, fnmfmt.format('', 'log'))
            if not all(os.path.isfile(os.path.join(datadir,fn)) for fn in (fnms[0], fnlog)):
                # java and bbmap need to be properly installed
                # run BBDuk to clean reads - note: bash isn't available on Windows
                cmd = f"java -ea -Xmx1g -cp {bbmapd} jgi.BBDuk in1={fnr1} in2={fnr2} out1={fncs[0]} out2={fncs[1]} "+\
                        f"t=2 qtrim=rl trimq=20 minlen=50 k=23 ktrim=r mink=11 hdist=1 overwrite=true ref=..\bbmap\resources\adapters.fa"
                pres1 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
            
                # run BBMerge to join paired-end reads
                cmd = f"java -ea -Xmx1g -cp {bbmapd} jgi.BBMerge in1={fncs[0]} in2={fncs[1]} out={fnms[0]} "+\
                        f"outu1={fnms[1]} outu2={fnms[2]} overwrite=true pfilter=1"
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
                    log_entries.append(f"Info: Skipping cleaning and merging of reads for {fnm_fmt}")
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
                log_entries.append(f"Error: Pair count mismatch. {cleanCount} != {cleanCount2}")
                with lock_l:
                    for l in log_entries:
                        logs.append(l)
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
                log_entries.append(f"Error: Merged counts mismatch. {mergeCount} != {joinCount}") 
                with lock_l:
                    for l in log_entries:
                        logs.append(l)
                with lock_r:
                    results.append(tuple(wr.values()))
                return
            
            assayfam = wr['primer'].split('_',1)[0] # assay family name for this well/sample
            #assayfams = [p.split('_',1)[0] for p in wr['primer']] # multiple assay families in same well!!!
            prime_targets = [] #[t for t in targets if any([af.lower() in t[0].lower() for af in assayfams])]
            other_targets = []
            for t in targets:
                if assayfam.lower() in t[0].lower():
                    prime_targets.append(t)
                else:
                    other_targets.append(t)
            
            match_cnt = Counter()
            low_cov_cutoff = max(mergeCount/20, 50)
            other_count = 0
            mtc = dict(match_cache)
            original_match_cache_size = len(mtc)
            msc = dict(miss_cache)
            original_miss_cache_size = len(msc)
            for seq, num in seqcnt.most_common():
                if num >= low_cov_cutoff or exhaustive:
                    if debug:
                        log_entries.append(f"Debug: Processing {wr['pcrplate']} {wr['pcrwell']} on process {PID} with counts {num}")
                    if seq in mtc:
                        if debug:
                            log_entries.append(f"Debug: {PID} Match cache hit {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        if sum_parent and varsep in mtc[seq]:
                            match_cnt[mtc[seq].split(varsep)[0]] += num
                        match_cnt[mtc[seq]] += num
                    elif seq in msc:
                        if debug:
                            log_entries.append(f"Debug: {PID} Miss cache hit {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        if sum_parent and varsep in msc[seq]:
                            match_cnt[msc[seq].split(varsep)[0]] += num
                        match_cnt[msc[seq]] += num
                    elif exact_match(match_cnt, seq, num, prime_targets, mtc):
                        if debug:
                            log_entries.append(f"Debug: {PID} Exact match against prime targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        pass
                    elif exact_match(match_cnt, seq, num, other_targets, mtc):
                        if debug:
                            log_entries.append(f"Debug: {PID} Exact match against other targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        pass
                    elif best_match(match_cnt, seq, num, prime_targets, mtc):
                        if debug:
                            log_entries.append(f"Debug: {PID} Inexact match against prime targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        pass
                    elif best_match(match_cnt, seq, num, other_targets, mtc):
                        if debug:
                            log_entries.append(f"Debug: {PID} Inexact match against other_targets {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        pass
                    else:
                        if debug:
                            log_entries.append(f"Debug: {PID} Missed {wr['pcrplate']} {wr['pcrwell']} {num} {seq}")
                        msc[seq] = seq
                        match_cnt[seq] += num
                else:
                    if debug:
                        log_entries(f"Warning: Skipped {wr['pcrplate']} {wr['pcrwell']} on process {PID} with counts {num}")
                    match_cnt['other'] += num
                    other_count += 1
            if len(mtc) > original_match_cache_size: 
                with lock_mtc:
                    match_cache.update(mtc)
            if len(msc) > original_miss_cache_size:
                with lock_msc:
                    miss_cache.update(msc)
                    
            # rename 'other' to include number of separate sequence variations
            match_cnt['other ('+str(other_count)+')'] = match_cnt['other']
            match_cnt['other'] = 0
            res1 = readCount, cleanCount, mergeCount
            if not exhaustive:
                try:
                    resx = sorted(((match_cnt[key], key) for key in match_cnt if match_cnt[key] >= low_cov_cutoff), reverse=True)
                except TypeError:
                    log_entries.append(f"Critical: match_cnt: {match_cnt}")        
            else:
                resx = sorted(((match_cnt[key], key) for key in match_cnt), reverse=True)
            res2 = tuple(x for p in resx for x in p)
            retval = tuple((str(x) for x in tuple(wr.values()) + res1 + res2))
            with lock_r:
                results.append(retval)            
        
    except Exception as exc:
        log_entries.append(f"Error: failed to run matching process {exc}")
    with lock_l:
        for l in log_entries:
            logs.append(l)
    return


def match_alleles_main(exp, ncpus, chunk_size, exhaustive, stagefile):
    try:
        with open(stagefile, 'rt') as stage_fd:
            src = csv.reader(stage_fd, dialect="unix")
            hdr = next(src)
            WRec = collections.namedtuple("WRec", hdr)
            wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrplate, x.pcrwell[0], int(x.pcrwell[1:]), x.primer))
        
        exp.log(f"Info: {len(wdata)} sample wells to process")
        hdrres1 = ("readCount", "cleanCount", "mergeCount")
        hdrres2 = ("seqCount", "seqName")
        complete_row_hdr = tuple((x for xs in (hdr, hdrres1, hdrres2) for x in xs))
        ## get a set of assay family names - family name ends with first underscore char  
        #assays = frozenset(r.primer.split('_',1)[0] for r in wdata)
        
        # read the target sequence file into a dictionary of target sequences
        # make sure there are no lowercase chars in target sequences
        fadict = {}
        for group in exp.reference_sequences:
            for id in exp.reference_sequences[group]:
                fadict[id] = Target(id, exp.reference_sequences[group][id])
            
        exp.log(f"Info: Read {len(fadict)} target reference sequences")

        results_rodentity = {}  # dict so we can order the outputs
        results_custom = {}
        # process data for each sample well    

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
            NUMPROCS = ncpus
           
            pool = mp.Pool(NUMPROCS)
            chunksize = chunk_size
            reports = []
            
            for i, chunk in enumerate([wrs[i:i+chunksize] for i in range(0,len(wrs),chunksize)]):
                while True:
                    reports_waiting = [r for r in reports if not r[0].ready()]
                    if len(reports_waiting) >= NUMPROCS:
                        time.sleep(0.05 + chunksize/20)
                        continue
                    else:
                        break
                # update caches before we try to launch any more processes
                with lock_mtc:
                    matchc = match_cache
                with lock_msc:
                    missc = miss_cache
                     
                exp.log(f"Debug:Adding work to pool... chunk {i*chunksize}, cache size {len(matchc)}, miss cache disabled", file=sys.stdout)
                r1 = pool.apply_async(process_well, (i, chunk, targets, matchc, missc, results, log, lock_mtc, 
                        lock_msc, lock_r, lock_l, exhaustive, exp.get_raw_dirpath(), exp.get_clean_dirpath(), exp.get_merged_dirpath()))
                reports.append((r1,chunk))
                #if i%4 == 0 and i != 0:
                #    time.sleep(2 + chunksize/10)

            exp.log('Debug: Waiting on all processes to return')
            while True:
                if any([not r[0].ready() for r in reports]):
                    exp.log('Debug: Waiting 1 minute for jobs to complete... you may see many of these messages')
                    time.sleep(60)
                else:
                    break

            for r in reports:
                if r[0].ready():
                    if not r[0].successful():
                        exp.log(f"Failure: with record {r[0]} {r[1]}")
                try:
                    r[0].get(60)
                    exp.log(f"Debug: Final report {r[0]} {r[1]}")
                except mp.TimeoutError:
                    exp.log(f"Error: Process timeout {r[0]} {r[1]}")
                    r[0].get(1)
            
            exp.log("Info: All processes completed")
            
            complete_results = {}
            for k,r in enumerate(results):
                print(f"{k=} {r=}")
                return
            #    # order results by pcrplate and pcrwell
            #    if args.custom:
            #        try:
            #            complete_results[(r[9],int(r[10][1:]),r[10][0])] = r
            #        except IndexError as e:
            #            print(e, r, file=sys.stdout)
            #            exit(1)
            #    else:
            #        complete_results[(r[15],int(r[16][1:]),r[16][0])] = r
                
            #for i,k in enumerate(sorted(complete_results.keys())):
            #    dst.writerow(complete_results[k])

            pool.close()
            pool.join()
        exp.log("Info: writing Rodentity records to results_rodentity.csv")
        with open(exp.get_exp_fp('results_rodentity.csv')) as f_rod:
            print(complete_row_hdr, file=f_rod)
            for res in results_rodentity:
                print(res, file=f_rod)
        exp.log("Info: writing custom records to results_custom.csv")
        with open(exp.get_exp_fp('results_custom.csv')) as f_cust:
            print(complete_row_hdr, file=f_cust)
            for res in results_custom:
                print(res, file=f_cust)

    except Exception as exc:
        exp.log(f'Failure: could not complete processing of sequences from FASTQ {exc}')
        return
    exp.log('Success: completed calculating allele counts from sequencing data')


def match_alleles(exp, ncpus=os.cpu_count(), chunk_size=1, exhaustive=False, stagefile=None, force_restart=False):
    """
    Entry point that handles locking of multiprocessing job.
    Process FASTQ files by cleaning and merging pairs
    Then match merged sequences against known targets and produce the NGS report.
    """
    if stagefile is None:
        stagefile = exp.get_exp_fp('Stage3.csv')

    lock_fn = exp.get_exp_fp('ngsgeno_lock')
    if force_restart:
        if os.path.exists(lock_fn):
            try:
                os.remove(lock_fn)
            except Exception as exc:
                exp.log('Critical: Could not remove lock file! Check process list in Task Manager')
                return
            exp.log('Info: processing lock file removed')
    if not os.path.exists(lock_fn):
        with open(lock_fn,"wt"):
            try:
                exp.log(f'Begin: match_alleles() {ncpus=} {chunk_size=} {exhaustive=} {stagefile=}')
                exp.log(f'Info: creating process lock file')
                match_alleles_main(exp, ncpus, chunk_size, exhaustive, stagefile)
            except Exception as exc:
                exp.log(f"Error: {exc}")
    else:
        exp.log("Warning: Analysis already running - lock file exists")
        return

    try:
        os.remove(lock_fn)
    except Exception as exc:
        exp.log(f'Error: could not remove processing lock file {exc}')
        return
    exp.log(f'Info: processing lock file successfully removed')
    
   
if __name__=="__main__":
    print('Call as a library only')
    
    
