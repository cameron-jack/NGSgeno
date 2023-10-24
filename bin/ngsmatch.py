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
from timeit import default_timer as timer
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

def myopen(fn):
    """ handles gz files seamlessly """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt")
    return open(fn, errors="ignore")


def exact_match(seq, refseqs, margin=70):
    """
    Exact matching a merged read with a set of reference target sequences
    Return True/False, matching reference sequence/None
    If the length difference is greater than margin, don't try to match 
    """
    for rs in refseqs:
        if abs(len(rs) - len(seq)) > margin:
            continue  # skip matching when length is too mismatched
        if len(seq) < len(rs):
            if seq in rs:
                return True, rs
        else:
            if rs in seq:
                return True, rs
    return False, None


def join_diffs(last_diffs):
    """ 
    Create a new diff string from a contiguous set of same-type differences
    Positions should be zero-based
    """
    pos = last_diffs[0][0] + 1
    ref_diff = ''.join([rc for pos,rc,ac in last_diffs])
    alt_diff = ''.join([ac for pos,rc,ac in last_diffs])
    if ref_diff.startswith('-'):  # insert in alternative
        return ''.join([str(pos),'+',alt_diff])
    elif alt_diff.startswith('-'):  # deletion in alternative
        return ''.join([str(pos),'-',ref_diff])
    # substitution
    return ''.join([str(pos),ref_diff,'/',alt_diff])


def annotate_seq(alnref, alnalt):
    """
    Return a string of PosRef/Alt e.g. 8-ATG11+TGC14A/G
    Ignore leading and trailing unmatched positions (ref or alt ---), only for global alignments
    """
    #print(f'Ref: {alnref}', flush=True)
    #print(f'Alt: {alnalt}', flush=True)
    lmargin = max(len(alnref)-len(alnref.lstrip('-')),len(alnalt)-len(alnalt.lstrip('-')))
    rmargin = max(len(alnref)-len(alnref.rstrip('-')),len(alnalt)-len(alnalt.rstrip('-')))
    #print(f'{lmargin=} {rmargin=}', flush=True)
    all_diffs = []
    last_diffs = []  # contiguous diffs of the same type: substitution, deletion, insertion
    for i,(rc, ac) in enumerate(zip(alnref, alnalt)):
        if i < lmargin:
            continue
        if i == len(alnref)-rmargin:
            break
        if rc != ac:
            if rc == '-':  # insertion in alnalt
                if last_diffs:
                    if last_diffs[-1][1] == '-':
                        last_diffs.append((i,rc,ac))
                    else:  # clear previously seen diffs
                        all_diffs.append(join_diffs(last_diffs))
                        last_diffs = [(i,rc,ac)]
                else:
                    last_diffs = [(i,rc,ac)]
            elif ac == '-':  # deletion in alnalt
                if last_diffs:
                    if last_diffs[-1][2] == '-':
                        last_diffs.append((i,rc,ac))
                    else:  # clear previously seen diffs
                        all_diffs.append(join_diffs(last_diffs))
                        last_diffs = [(i,rc,ac)]
                else:
                    last_diffs = [(i,rc,ac)]
            else:  # substitution
                if last_diffs:
                    if last_diffs[-1][1] != '-' and last_diffs[-1][2] != '-':
                        last_diffs.append((i,rc,ac))
                    else:
                        all_diffs.append(join_diffs(last_diffs))
                        last_diffs = [(i,rc,ac)]
                else:
                    last_diffs = [(i,rc,ac)]
        else:    
            if last_diffs:
                all_diffs.append(join_diffs(last_diffs))
                last_diffs = []
    if last_diffs:
        all_diffs.append(join_diffs(last_diffs))

    return ''.join(all_diffs)


def test_annotation():
    """
    Some test cases to ensure that annotation of two aligned sequences is working correctly
    """
    refseq = 'ATTACGGTCT'
    altseq = 'ATTACGGTCT'
    ds = annotate_seq(refseq, altseq)
    if ds != '':
        return False
    altseq = '--TACGGTCT'
    ds = annotate_seq(refseq, altseq)
    if ds != '':
        return False
    altseq = '--TACGGACT'
    ds = annotate_seq(refseq, altseq)
    if ds != '8T/A':
        return False
    altseq = '--TACGG-CT'
    ds = annotate_seq(refseq, altseq)
    if ds != '8-T':
        return False
    refseq = 'ATTACGGT-CT'
    altseq = '--TACGGTTCT'
    ds = annotate_seq(refseq, altseq)
    if ds != '9+T':
        return False
    refseq = '--TACGGTCT--'
    altseq = 'ATTACCTTCTAA'
    ds = annotate_seq(refseq, altseq)
    if ds != '6GG/CT':
        return False
    return True


def inexact_match(seq, refseqs, rundir, debug, lock_d, identity=0.75, margin=70):
    """
    Inexact matching of a merged sequence with a set of reference target sequences
    Declare a match if we have 75% sequence identity over the matched length or better
    Identity must be >= 0 and <= 1.0
    Return True/False, and matching reference sequence/None
    If the length difference is greater than margin, don't try to match
    rundir, debug, lock_d are used for logging debug messages
    """
    if identity < 0 or identity > 1.0:
        return False, None, None
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -2.5
    aligner.extend_gap_score = -1.0
    aligner.end_gap_score = 0
    #aligner.match_score = 1.0
    aligner.mismatch_score = -1.5
    try:
        # do alignment against a reference sequence if they aren't too different in length
        all_alignments = [(aligner.align(rs, seq)[0], rs) for rs in refseqs if abs(len(rs)-len(seq))<margin]
    except Exception as exc:
        msg = f'alignment failed with exception {exc}'
        print(msg, flush=True)
        wdb(msg, rundir, debug, lock_d)
        return False, None, None

    if len(all_alignments) == 0:
        wdb(f'No alignments found {seq=}', rundir, debug, lock_d)
        return False, None, None
    best_score = 0  # best scoring reference sequence according to alignment rules, favours substitutions
    aligned_rs = None  # the aligned reference sequence
    aligned_seq = None  # the aligned query sequence
    aligned_bars = 0  # the number of exact matches - used to determine sequence identity
    refseq = None  # the actual best matching reference sequence
    for algn, rs in all_alignments:
        if algn.score > best_score:
            best_score = algn.score
            target_sections = []
            query_sections = []
            bar_sections = []
            for segment in str(algn).split('\n'):
                tokens = [s for s in segment.split(' ') if s != '']
                if len(tokens) < 3:  # skip rows without actual sequence or bars
                    continue
                if tokens[-1][-1] in '0123456789':
                    if tokens[0] == 'target':
                        target_sections.append(tokens[-2])
                    elif tokens[0] == 'query':
                        query_sections.append(tokens[-2])
                    else:
                        bar_sections.append(tokens[-2])
                else:
                    if tokens[0] == 'target':
                        target_sections.append(tokens[-1])
                    elif tokens[0] == 'query':
                        query_sections.append(tokens[-1])
                    else:
                        bar_sections.append(tokens[-1])
            aligned_rs = ''.join(target_sections)
            aligned_seq = ''.join(query_sections)
            aligned_bars = ''.join(bar_sections).count('|')
            refseq = rs

    # print(f'{best_score=} {aligned_bars=} {aligned_rs=} {aligned_seq=} {refseq=}', flush=True)
    if aligned_bars >= min(len(seq.strip('-')), len(aligned_rs.strip('-'))) * identity:
        # Annotate match
        seq_anno = annotate_seq(aligned_rs, aligned_seq)
        return True, refseq, seq_anno
    msg = f'Best alignment identity worse than 75% {aligned_rs=} {aligned_seq=} {refseq=} {seq=}'
    wdb(msg, rundir, debug, lock_d)
    return False, None, None


def preprocess_seqs(wr, rundir, results, log, lock_r, lock_l, lock_d, debug=False):
    """
    Call bbduk for cleaning of reads
    Call bbmerge to merge read pairs
    If number of merged reads is less than half the number of unmerged reads then call this a failure
    Return seqcnt (Counter of unique sequences) and passing (T/F)
    """
    rfn= "*{}-{}_*_R1_*.fastq.gz".format(util.unguard(wr['pcrPlate'],silent=True), padwell(wr['pcrWell']))
    fn1s = glob.glob(os.path.join(rundir, "raw", rfn))
    lrecs = []
    if not fn1s:
        with lock_l:
            log.append(f"no data for {fn1s}")
        wdb(f"no data for {fn1s}", rundir, debug, lock_d)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
            
    if len(fn1s)>1:
        with lock_l:
            log.append(f"too many reads files for pattern{rfn}")
            log.append(f"   R1 files= {fn1s}")
        wdb(f"too many reads files for pattern{rfn}\n    R1 files= {fn1s}", rundir, debug, lock_d)
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
        wdb(f"missing file: {fnr2}", rundir, debug, lock_d)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False
        
    # use Windows file separators
    fncs = tuple(os.path.join(rundir,"cleaned", fn) for fn in (fn1, fn2))
    bbmapd = os.path.join('bbmap','current')
    # unpick the file name a bit
    fnparts = fn1.split('_', 2)
    fnmfmt = "{}-{}_{}{{}}.{{}}".format(wr['pcrPlate'], padwell(wr['pcrWell']), fnparts[1])
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
        lrecs.append(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}")
        wdb(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}", rundir, debug, lock_d)
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
        msg = f"Pair count mismatch. {cleanCount} != {cleanCount2}"
        with lock_l:
            log.append(msg)
        wdb(msg, rundir, debug, lock_d)
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
            wdb(msg, rundir, debug, lock_d)
        with lock_r:
            results.append(tuple(wr.values()))
        return {}, False

    with lock_l:
        for l in lrecs:
            log.append(l)
    
    return seqcnt, True, readCount, cleanCount, mergeCount


def wdb(msg, rundir, debug, lock_d=None):
    """ 
    wdb = write debug.
    debug [T/F] - do this here to reduce shared logic
    lock_d: Call from within a lock 
    """
    #print (f'In wdb {msg=}')
    if debug:
        debugfn = os.path.join(rundir, 'debug.log')
        if lock_d:
            with lock_d:
                try:
                    with open(debugfn, 'at') as df:
                        print(msg, file=df, flush=True)
                except:
                    print('Could not write to debug')
                    print(f'{msg}', flush=True)
        else:
            try:
                with open(debugfn, 'at') as df:
                    print(msg, file=df, flush=True)
            except:
                print('Could not write to debug', flush=True)
                print(f'{msg}')
    #print('leaving wdb')


def archetypes(seq_counts, rundir, debug, lock_d):
    """ 
    Collapse substrings and superstrings
    inputs: dictionary of sequences and counts
    output: collections.Counter
    """
    start = timer()
    final_seqs = {}
    unique_seqs = set(seq_counts.keys())
    msg = f'Start of archetypes. Number of sequences {len(seq_counts)} total counts {sum([seq_counts[s] for s in seq_counts])}'
    wdb(msg, rundir, debug, lock_d)
    
    unique_seqs = sorted(unique_seqs, key=len)
    skip_indices = set()
    kept_seqs = []
    for i, u1 in enumerate(unique_seqs):
        if i in skip_indices:
            continue
        kept_seqs.append(u1)
        skip_indices.add(i)
        for j, u2 in enumerate(unique_seqs):
            if j <= i:
                continue
            if u2 in u1:
                skip_indices.add(j)
                seq_counts[u1] += seq_counts[u2]
    for ks in kept_seqs:
        final_seqs[ks] = seq_counts[ks]

    #print(f'Final {collections.Counter(final_seqs).most_common()=}')
    msg = f'End of archetypes. Number of sequences {len(final_seqs)} total counts {sum([final_seqs[s] for s in final_seqs])}'
    wdb(msg, rundir, debug, lock_d)
   
    end = timer()
    wdb(f'Archetypes took {end-start}', rundir, debug, lock_d)
    
    return collections.Counter(final_seqs)


def process_well(work_block, wr, rundir, seq_ids, id_seq, primer_assayfam, assayfam_primers, 
        match_cache, anno_cache, miss_cache, results, reps, log, lock_mtc, lock_msc, lock_r, 
        lock_l, lock_d, mincov, minprop, exhaustive, debug=False):
    """ 
    Worker process that runs bbduk, bbmerge, and matching
    work_block: work block number
    wr: well record
    seq_ids:  {seq:set([ids]}
    id_seq: {seq_id:seq}
    primer_assayfam: {primer:assayfam}
    assayfam_primers: {assayfam:[primers]}
    match_cache: dict of observed sequence to known sequence
    anno_cache: dict of annotations for variant sequences
    miss_cache: dict of unknown sequences, used as a set
    results: shared list of results + log entries (might separate in future)
    lock_mtc: lock object for the match cache
    lock_msc: lock object for the miss cache
    lock_r: lock object for the results    
    lock_l: lock for logging
    lock_d: lock for writing to debug
    exhaustive: [T/F] don't skip any sequences, not matter how low their count
    debug: write messages direct to file (slow I/O)

    Parallel, independent processing of a single well.
    Clean and merge reads using bbduk
    Group sequences into archetypes and order by decreasing count
    Check whether any sequences exactly match the expected assay family and base cutoffs on this
    Check all remaining sequences until counts falls below cutoff
      -> Exact matches against the match cache (includes all assay families)
      -> Exact matches against the miss cache
      -> Inexact matches against the match cache -> add variants to match cache
      -> Inexact matches against the miss cache -> add variants to miss cache
      -> Add any remaining entries to the miss cache

    Match cache ["ATC..."]="ATC..." <- a sequence in seq_ids and id_seq
    Miss cache
    """
    PID = os.getpid()
    lrecs = [] # log records, lock and write on exit/return

    msg = f"Info:{PID} Working on: {wr['pcrPlate']} {wr['pcrWell']}"
    lrecs.append(msg)
    wdb(msg, rundir, debug, lock_d)
    print(msg, flush=True)

    # clean and merge FASTQs
    seqcnt, success, readCount, cleanCount, mergeCount = preprocess_seqs(wr, rundir, results, log, 
            lock_r, lock_l, lock_d, debug=debug)
    lrecs.append(f'{log=} {success=} {wr=}')
    if not success:
        msg = f"Failed: to completed preprocessing for {wr['pcrPlate']} {wr['pcrWell']}"
        lrecs.append(msg)
        wdb(msg, rundir, debug, lock_d)
        retval = [str(wr[x]) for x in wr]
        with lock_r:
            reps[int(wr['sampleNo'])] = retval
        wdb(f'{retval=}', rundir, debug, lock_d)
        with lock_l:
            for l in lrecs:
                log.append(l)
        return False

    msg = f"Info: Completed preprocessing for {wr['pcrPlate']} {wr['pcrWell']} "+\
            f"{readCount=} {cleanCount=} {mergeCount}"
    lrecs.append(msg)
    wdb(msg, rundir, debug, lock_d)
                
    primer = wr['primer']
    if primer not in primer_assayfam:
        msg = f'Warning: primer {primer} unknown, skipping. {primer_assayfam=}'
        lrecs.append(msg)
        wdb(msg, rundir, debug, lock_d)
        retval = [str(wr[x]) for x in wr]
        with lock_r:
            reps[int(wr['sampleNo'])] = retval
        wdb(f'{retval=}', rundir, debug, lock_d)
        with lock_l:
            for l in lrecs:
                log.append(l)
        return False
        
    print(f"Processing {wr['pcrPlate']} {wr['pcrWell']} with {mergeCount} reads and primer {primer}", flush=True)

    # decide which assays are on-target vs off-target
    on_target_ids = set()
    off_target_ids = set()
    on_target_seqs = set()
    off_target_seqs = set()
    for name in id_seq:
        if name.lower().startswith(primer.split('_')[0].lower()):
            on_target_ids.add(name)
            on_target_seqs.add(id_seq[name])
        else:
            off_target_ids.add(name)
    # avoid looking at the same sequence under different names - used by best_match
    off_target_seqs = set(list(seq_ids.keys())).difference(on_target_seqs)
    #print(f"{wr['pcrWell']} {on_target_ids} {on_target_seqs}", flush=True)
    msg = f"Debug: {on_target_ids=} {len(off_target_seqs)=}"
    wdb(msg, rundir, debug, lock_d)
    match_cnt = Counter()
    if mincov < 1:
        mincov = 1
    if minprop > 1.0:
        minprop = 1.0
    elif minprop < 0.0:
        minprop = 0.0
        
    # unique sequences only (substrings collapsed), from most to least common
    if exhaustive:
        print('Aggregating archetype sequences', flush=True)
        seqcnt = archetypes(seqcnt, rundir, debug, lock_d)
    # calculate min count proportion from exact matches to our expected targets
    family_exact_counts = 0
    for on_target_id in on_target_ids:
        on_target_seq = id_seq[on_target_id]
        if on_target_seq in seqcnt:
            family_exact_counts += seqcnt[on_target_seq]
    low_cov_cutoff = max(family_exact_counts*minprop, mincov)

    other_count = 0
    mtc = dict(match_cache)  # force copy of shared object
    original_match_cache_size = len(mtc)
    anc = dict(anno_cache)  # force copy of shared object
    original_anno_cache_size = len(anc)
    msc = dict(miss_cache)  # force copy of shared object
    original_miss_cache_size = len(msc)

    for seq, num in seqcnt.most_common():
        msg=f"{wr['pcrPlate']} {wr['pcrWell']} with {seq=} and counts {num}"
        #wdb(msg, rundir, debug, lock_d)
        if (num >= low_cov_cutoff) or exhaustive:
            msg = f"Processing {wr['pcrPlate']} {wr['pcrWell']} on process {PID} with counts {num} and {seq}"
            wdb(msg, rundir, debug, lock_d)      
                
            if seq in mtc:
                msg = f"{PID}: Match cache hit {num} {seq}"
                wdb(msg, rundir, debug, lock_d)
                if seq in anc:
                    anno = anc[seq]
                    match_cnt[seq+'//'+anno] += num
                else:
                    match_cnt[seq] += num
                # top up the original reference counts if needed
                ref_seq = mtc[seq]
                if ref_seq != seq: 
                    match_cnt[ref_seq] += num
                continue
                
            # miss cache - doesn't match any known sequence
            if seq in msc:
                msg = f"{PID}: Miss cache hit {num} {seq}"
                wdb(msg, rundir, debug, lock_d)
                match_cnt[seq] += num
                continue

            # substring or superstring
            is_match, ref_seq = exact_match(seq, on_target_seqs)
            if is_match:
                msg = f"{PID}: Exact match against {seq_ids[ref_seq]} with {seq} and counts {num}"
                wdb(msg, rundir, debug, lock_d)
                match_cnt[ref_seq] += num
                continue

            is_match, ref_seq = exact_match(seq, off_target_seqs)
            if is_match:
                msg = f"{PID}: Exact match against {seq_ids[ref_seq]} with {seq} and counts {num}"
                wdb(msg, rundir, debug, lock_d)
                match_cnt[ref_seq] += num
                continue

            is_match, ref_seq, seq_anno = inexact_match(seq, on_target_seqs, rundir, debug, lock_d)
            if is_match:
                msg = f"{PID}: Inexact match against {seq_ids[ref_seq]} {wr['pcrPlate']} {wr['pcrWell']} {num} {seq}\n"
                wdb(msg, rundir, debug, lock_d)
                mtc[seq] = ref_seq
                match_cnt[ref_seq] += num
                if seq_anno:
                    anc[seq] = seq_anno
                    match_cnt[ref_seq+'//'+seq_anno] += num
                else:
                    wdb(f'No annotation. Ref: {ref_seq} Alt: {seq}', rundir, debug, lock_d)
                continue

            is_match, ref_seq, seq_anno = inexact_match(seq, off_target_seqs, rundir, debug, lock_d)
            if is_match:
                msg = f"{PID}: Inexact match against {seq_ids[ref_seq]} {wr['pcrPlate']} {wr['pcrWell']} {num} {seq}\n" 
                wdb(msg, rundir, debug, lock_d)
                mtc[seq] = ref_seq
                match_cnt[ref_seq] += num
                if seq_anno:
                    anc[seq] = seq_anno
                    match_cnt[ref_seq+'//'+seq_anno] += num
                else:
                    wdb(f'No annotation. Ref: {ref_seq} Alt: {seq}', rundir, debug, lock_d)
                continue

            # add to miss cache
            msg = f"{PID}: No match for {wr['pcrPlate']} {wr['pcrWell']} {num} {seq}\n"
            wdb(msg, rundir, debug, lock_d)
            msc[seq] = None
            match_cnt[seq] += num
                    
        else:
            #msg = f"{PID}: too few reads {num} for matching {seq}"
            #wdb(msg, rundir, debug, lock_d)
            match_cnt['other'] += num
            other_count += 1
        # try another sequence

    if len(mtc) > original_match_cache_size:
        with lock_mtc:
            match_cache.update(mtc) 
            anno_cache.update(anc)
        wdb('Info: Updating match and annotation caches', rundir, debug, lock_d)
           
    if len(msc) > original_miss_cache_size:
        with lock_msc:
            miss_cache.update(msc)
        wdb('Updating miss cache', rundir, debug, lock_d)

    msg = f"Completed matching for {PID=} {wr['pcrPlate']=} {wr['pcrWell']=}"
    lrecs.append(msg)
    wdb(msg, rundir, debug, lock_d)
    #print(msg, flush=True)
    # rename 'other' to include number of separate sequence variations
    total_others = match_cnt['other']
    match_cnt['other ('+str(other_count)+')'] = total_others
    del match_cnt['other']
    res1 = [readCount, cleanCount, mergeCount]
    # name the outputs and counts
    seqCounts = []
    seqNames = []
    otherCounts = []
    otherNames = []
    for seq, count in match_cnt.most_common():
        if seq.startswith('other') or seq in msc:
            otherNames.append(seq)
            otherCounts.append(count)
            continue

        if '//' not in seq:
            if seq in on_target_seqs:
                seqCounts.append(count)
                seqNames.append([name for name in seq_ids[seq] if name in on_target_ids][0])
            else:
                for name in seq_ids[seq]:
                    otherCounts.append(count)
                    otherNames.append(name)
            continue
        ref_seq = seq.split('//')[0]
        anno = seq.split('//')[1]
        for name in seq_ids[ref_seq]:
            otherCounts.append(count)
            otherNames.append(name+'//'+anno)

    try:
        res2 = [';'.join(map(str,seqCounts)), ';'.join(seqNames), ';'.join(map(str,otherCounts)), ';'.join(otherNames)]
    except Exception as exc:
        wdb(f"{exc=}", rundir, debug, lock_d)
        res2 = ['','','','']
    retval = [str(wr[x]) for x in wr] + res1 + res2
    with lock_r:
        reps[int(wr['sampleNo'])] = retval
    wdb(f'{retval=}', rundir, debug, lock_d)
    with lock_l:
        for l in lrecs:
            log.append(l)
    msg = f"Exiting process_well() {PID} {wr['pcrPlate']} {wr['pcrWell']}"
    wdb(msg, rundir, debug, lock_d)
    print(msg, flush=True)


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


def report_progress(rundir, match_progress):
    """
    Clear all previous progress files and touch a new file with the launch and match progress values in the name
    Only do the operation is progress is a multiple of 5 to save on disk writes
    """
    launch_progress = match_progress
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


def parse_targets(rundir, targets, log, debug):
    """
    Return a dictionary of sequences to a set of ids, and
    a dictionary of ids to sequences
    """
    seq_ids = {}
    id_seq = {}
    with myopen(os.path.join(rundir,targets)) as src:
        for entry in Bio.SeqIO.parse(src, "fasta"):
            if len(entry.id) > 0 and len(entry.seq) > 0:
                try:
                    name = bytes(str(entry.id).strip(), 'ascii', errors='strict').decode()
                except Exception as exc:
                    msg = f'Warning: {exc} non-ascii character in reference sequence id {entry.id}'
                    log.append(msg)
                    wdb(msg, args.rundir, debug)
                    seq_id = bytes(str(entry.id).strip(), 'ascii',errors='ignore').decode()

                try:
                    seq = bytes(str(entry.seq).strip(), 'ascii', errors='strict').decode().upper()
                except Exception as exc:
                    msg = f'Warning: {exc} non-ascii character in reference sequence {entry.seq}'
                    log.append(msg)
                    wdb(msg, args.rundir, debug)
                    seq = bytes(str(entry.seq).strip(), 'ascii', errors='ignore').decode().upper()
                if seq not in seq_ids:
                    seq_ids[seq] = set()
                seq_ids[seq].add(name)
                if name in id_seq:
                    msg = f'Warning: {name} duplicated in reference'
                    log.append(msg)
                    wdb(msg, args.rundir, debug)
                id_seq[name] = seq
    return seq_ids, id_seq


def parse_primer_assayfams(rundir, paf):
    """
    Two column file (tab separated) of primer and assay families
    Return a mapping of primer to set of assay family, and mapping of assay family to set of primers
    """
    primer_assayfam = {}  # [primer] = []
    assayfam_primers = {}  # [assayfam] = []
    pmr_fn = os.path.join(rundir, paf)
    try:
        with open(pmr_fn, 'rt', errors='ignore') as f:
            for line in f:
                cols = line.strip().split(',')
                if len(cols) != 2:
                    print(f'skipping line {line=}')
                    continue
                primer = cols[0]
                assayfam = cols[1]
                if not primer or not assayfam:
                    continue
                if primer not in primer_assayfam:
                    primer_assayfam[primer] = set()
                primer_assayfam[primer].add(assayfam)
                if assayfam not in assayfam_primers:
                    assayfam_primers[assayfam] = set()
                assayfam_primers[assayfam].add(primer)
        return primer_assayfam, assayfam_primers
    except Exception as exc:
        msg = f'Critical: failed to open {pmr_fn} for primer assay family info {exc=}'
        print(msg)
    return primer_assayfam, assayfam_primers

def parse_primer_file(paf):
    """
    Parse the user provided assay list file. Use this if testing without a primer_assayfam file
    """
    primer_assayfam = {}
    assayfam_primers = {}  # [assayfam] = []
    with open(paf, 'rt') as f:
        for line in f:
            cols = [c.strip() for c in line.split(',')]
            assayfam = cols[1]
            pmrs = [c for c in cols[2:] if len(c) != 0]
            if assayfam not in assayfam_primers:
                assayfam_primers[assayfam] = set()
            for pmr in pmrs:
                if pmr not in primer_assayfam:
                    primer_assayfam[pmr] = set()
                primer_assayfam[pmr].add(assayfam)
                assayfam_primers[assayfam].add(pmr)
    return primer_assayfam, assayfam_primers


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
        wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrPlate, x.pcrWell[0], int(x.pcrWell[1:]), x.primer))
    #print(wdata, file=sys.stderr)
    wdata = [rec for rec in wdata if util.unguard(rec.pcrPlate, silent=True) +'-'+ util.padwell(rec.pcrWell) in raw_file_identifiers]
    #print(f'After filtering by available files {wdata=}', file=sys.stderr)
    log.append(f"Info: {len(wdata)} sample wells to process.")
    ## get a set of assay family names - family name ends with first underscore char  
    #assays = frozenset(r.primer.split('_',1)[0] for r in wdata)
        
    # read the target sequence file into dictionaries [seq] = set([ids]), and [id] = seq             
    seq_ids, id_seq = parse_targets(args.rundir, args.targets, log, args.debug)
    print('Parsed targets', flush=True)

    # get relationship of primers to assay family
    #print(f'{args.primer_assayfam=}')
    #primer_assayfam, assayfam_primers = parse_primer_assayfams(args.rundir, args.primer_assayfam)
    primer_assayfam, assayfam_primers = parse_primer_file(os.path.join(args.rundir, args.primer_assayfam))
    print(f'Parsed primers/assays', flush=True) # {primer_assayfam=} {assayfam_primers=}')
    
    if not os.path.isdir(os.path.join(args.rundir,"raw")):
        log.append("Error: raw FASTQ data folder is absent - please transfer MiSeq data.")
        write_log(log, os.path.join(args.rundir,args.logfn))
        return
    
    for d in ("cleaned", "merged"):
        dp = os.path.join(args.rundir, d)
        if not os.path.isdir(dp):
            os.mkdir(dp)
                
    # parallel execute over wdata and collect results
    grouped_wrs = itertools.groupby(wdata, key=lambda x:(x.pcrPlate, x.pcrWell))
    wrs = []
    for key, group in grouped_wrs:
        for g in group:
            wrs.append(OrderedDict(g._asdict()))
        
    with mp.Manager() as manager:
        # Use a server process to manage shared data structures. Heavy option, but works on Windows
        match_cache = manager.dict()  # [seq] = known_seq
        anno_cache = manager.dict()  # [seq] = annotation
        miss_cache = manager.dict()  # just use this as a set
        results = manager.list()  # at the moment it holds both the main results and the log, but maybe we should split these out?
        reps = manager.dict()  # [sampleNo] = all columns
        logm = manager.list()
        lock_mtc = manager.Lock()  # match_cache locking
        lock_msc = manager.Lock()  # miss_cache locking
        lock_r = manager.Lock()  # result list locking
        lock_l = manager.Lock()  # log list locking
        lock_d = manager.Lock()  # debug file locking
      
        with lock_mtc:
            for seq in seq_ids:
                match_cache[seq] = seq  # variants of these will also go in here 

        #print('Before launching jobs', file=sys.stderr, flush=True)
        # multiprocessing
        NUMPROCS = args.ncpus
           
        pool = manager.Pool(NUMPROCS)
        reports = []
        launch_progress = 0
        match_progress = 0
        print('launching jobs', flush=True)
        total_jobs = len(wrs)
        # multiprocessing pool counter from https://superfastpython.com/multiprocessing-pool-asyncresult/
        reports = [pool.apply_async(process_well, args=(i, wr, args.rundir, seq_ids, id_seq, primer_assayfam, 
                assayfam_primers, match_cache, anno_cache, miss_cache, results, reps, logm, lock_mtc, 
                lock_msc, lock_r, lock_l, lock_d, args.mincov, args.minprop, args.exhaustive, args.debug))\
                for i, wr in enumerate(wrs)]

        count = len(results)
        while count:
            # check all tasks and count the number that are not done
            count = sum([not r.ready() for r in results])
            completed = sum([r.ready() for r in results])
            # report the number of remaining tasks
            print(f'{count}/{len(results)} tasks remain')
            match_progress = int(100*completed/total_jobs)
            report_progress(args.rundir, match_progress)
            # wait a moment
            time.sleep(0.5)       

        pool.close()
        pool.join()
        log.append("Info: All processes completed")
        with lock_l:
            for l in logm:
                log.append(l)

        log.append(f'Info: Writing results to {args.outfn}')
        
        with open(os.path.join(args.rundir,args.outfn), "wt", buffering=1) as dstfd:
            print(f"Opening {args.outfn} for results", file=sys.stderr)
            dst = csv.writer(dstfd, dialect="unix", quoting=csv.QUOTE_ALL)
            hdrres1 = ("readCount", "cleanCount", "mergeCount")
            hdrres2 = ("seqCount", "seqName", "otherCount", "otherName")
            complete_row_hdr = tuple((x for xs in (hdr, hdrres1, hdrres2) for x in xs))
            dst.writerow(complete_row_hdr)
            with lock_r:
                for key in sorted(reps.keys()):
                    dst.writerow(reps[key])
            dstfd.flush()

    end_time = datetime.datetime.now()
    log.append(f"End: {end_time} took: {end_time - start_time}")
    write_log(log, os.path.join(args.rundir,args.logfn))
   
    
def run_matches():
    """
    Entry point when used as a library. Calls main()
    """
    pass

   
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="NGS Reporting Program")
    parser.add_argument("-d", "--debug", action="store_true", help="more reporting/output for debugging purposes")
    parser.add_argument("-t", "--targets", default="targets.fa", help="file of targets in FASTA format (default=targets.fa)")
    parser.add_argument('-P', '--primer_assayfam', default="primers.csv", help='file of primer to assay family mappings')
    parser.add_argument('-o','--outfn', default='results.csv', help='Name of output file (CSV format)')
    parser.add_argument('-n','--ncpus', type=int, default=os.cpu_count()-1, help='Number of processes to run simultaneously, default=number of CPUs in system - 1')
    parser.add_argument('-l','--logfn', default='match.log', help='Name of logging file (default=match.log)')
    parser.add_argument('-r','--rundir', required=True, help='Path to experiment folder')
    parser.add_argument('-m','--mincov', type=int, default=50, help='Do not match unique sequences with less than this many reads coverage, default 50')
    parser.add_argument('-p','--minprop', type=float, default=0.2, help='Do not match unique sequences '+\
            'with less than this proportion of the total number of exact matched on-target reads, default 0.2. Must be between 0.0 and 1.0')
    parser.add_argument('-x','--exhaustive',action='store_true',help='Try to match every sequence, '+\
            'no matter how few counts. Ignores --minseqs and --minprop')
    parser.add_argument('-s','--stagefile', default="Stage3.csv", help="Name of the NGS genotyping Stage 3 file (default=Stage3.csv)")
    args = parser.parse_args()
    in_error = False
    print(f"{args=}", file=sys.stderr)
    lock_path = os.path.join(args.rundir,"ngsgeno_lock") 
    if not os.path.exists(lock_path):
        try:
            with open(lock_path,"wt"):
                report_progress(args.rundir, 0)  # set this up asap
                main(args)
                os.remove(lock_path)
        except:
            os.remove(lock_path)
    else:
        print("Analysis already running", file=sys.stderr)
        exit(2)
    
    
    
