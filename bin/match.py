#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Jun 2025
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@license: GNU General Public License v3.0
@description:
Uses bbduk and bbmerge to clean and merge reads from a MiSeq run in the 'raw' directory
Performs exact matching of reads to a set of reference sequences
Uses inexact matching to find variants in the reads
@depends: bbduk, bbmerge, BioPython, Python 3.6+
"""

#from ast import Str
from timeit import default_timer as timer
import os
from pathlib import Path
import sys
import argparse
import datetime
import re
import csv
import queue
import threading
#import concurrent.futures
import multiprocessing as mp
import warnings

import collections
from collections import Counter, OrderedDict
import itertools
import subprocess
import gzip
import glob
import time
from math import ceil, floor

import Bio.SeqIO
import Bio.Align as Align

from cogent3 import get_app, make_unaligned_seqs, make_tree
from cogent3.align.progressive import tree_align

try:
    from bin.util import unguard, padwell
except ModuleNotFoundError:
    from util import unguard, padwell

# bclen = 8 # length of pseudo barcode
# matchcutoff = 1.5 # should control via command line - must be >1.x

# def myopen(fn):
#     """ handles gz files seamlessly """
#     if fn.endswith('.gz') :
#         return gzip.open(fn, "rt")
#     return open(fn, errors="ignore")


# def exact_match(seq, refseqs, margin=0.9):
#     """
#     Exact matching a merged read with a set of reference target sequences
    
#     args:
#     seq: the read sequence to match
#     refseqs: list of reference target sequences
#     margin: seq and refseq must be within this factor of each other in length

#     returns:
#     True/False, matching reference sequence/None

#     If the length difference is greater than margin, don't try to match 
#     If there are parentheses in the sequence, these are masked regions
#     and we should split the sequence on these and match the unmasked regions.

#     Simple algorithm for a single masked region and known amplicons:
#     1. Cut out the two reference flanks at brackets and compare these to the sequence
#     2. If they match, then we have a match

#     Two bracket algorithm for a pair of masked regions:
#     1. Cut out the two reference flanks at parentheses
#     2. remove any square brackets, and compare these to the sequence
#     3. If they match, then we have a match
#     4. Assuming it's not a match repeat from 1, but finding the flanks of the square brackets
#     5. Remove any parentheses and compare these to the sequence
#     5. If they match, then we have a match

#     General algorithm for masked regions (more than 1):
#     1. split the reference sequence on masked regions
#     2. for each unmasked reference region, match this to the sequence
#     3. if all unmasked regions match, then we have a match
#     """
#     for rs in refseqs:
#         found_match = False
#         # compare length of sequenced region to reference seq
#         rs = rs.replace('[','(').replace(']',')')
#         if '(' not in rs:
#             if len(seq) > len(rs):
#                 if len(seq)*margin > len(rs):
#                     continue  # skip matching when length is too mismatched
#             else:
#                 if len(rs)*margin > len(seq):
#                     continue  # skip matching when length is too mismatched
#         # unmasked case
#         if rs.count('(') == 0 and rs.count(')') == 0:
#             if len(seq) < len(rs):
#                 if seq in rs:
#                     return True, rs
#             else:
#                 if rs in seq:
#                     return True, rs
#         # two bracket type - DEPRECATED
#         elif rs.count('(') == 1 and rs.count('[') == 1 and \
#                 rs.count(')') == 1 and rs.count(']') == 1:
#             left_flank = rs.split('(')[0].replace('[','').replace('[','')
#             right_flank = rs.split(')')[1].replace('[','').replace('[','')
#             if left_flank in seq:
#                 centre_pos = seq.find(left_flank) + len(left_flank)
#                 if right_flank in seq[centre_pos:]:
#                     return True, rs
#             else:
#                 left_flank = rs.split('[')[0].replace('(','').replace(')','')
#                 right_flank = rs.split(']')[1].replace('(','').replace(')','')
#                 if left_flank in seq:
#                     centre_pos = seq.find(left_flank) + len(left_flank)
#                     if right_flank in seq[centre_pos:]:
#                         return True, rs
#         # simple masked case
#         elif rs.count('(') == 1 and rs.count(')') == 1:
#             left_flank = rs.split('(')[0]
#             right_flank = rs.split(')')[1]
#             if left_flank in seq:
#                 centre_pos = seq.find(left_flank) + len(left_flank)
#                 if right_flank in seq[centre_pos:]:
#                     return True, rs
#         # general masked case with multiple masked regions    
#         elif rs.count('(') > 1 and rs.count(')') == rs.count('('):
#             seq_index = 0
#             failed_rs = False
#             for section in rs.split('('):
#                 if ')' in section:
#                     non_variable_seq = section.split(')')[1]
#                 else:
#                     non_variable_seq = section
#                 seq_index = seq.find(non_variable_seq, seq_index)
#                 if seq_index == -1:
#                     failed_rs = True
#                     break
#                 seq_index += len(non_variable_seq)  # move past the non-variable sequence
#             if not failed_rs:
#                 return True, rs        
#     return False, None


# def test_exact_match():
#     """
#     Some test cases to ensure that exact matching of a sequence to a set of reference sequences is working correctly
#     """
#     refseqs = ['ATGCGTGTTCAAGTACACCCAAGTTGACAGTGCA']
#     seq = 'ATGCGTGTTCAAGTACACCCAAGTTGACAGTGCA'
#     success, rs = exact_match(seq, refseqs)
#     if not success:
#         print('Test exact match failed')
#         return False
    
#     refseqs = ['ATGCGTG(TTC)AAGTACACCCAAGTTGACAGTGCA']
#     seq = 'ATGCGTGTCAAGTACACCCAAGTTGACAGTGCA'
#     success, rs = exact_match(seq, refseqs)
#     if not success:
#         print('Test exact match with round brackets failed')
#         return False
    
#     refseqs = ['ATGCGTG(TTC)AAGTACAC(CC)AAGTTGAC(A)GTGCA']
#     seq = 'ATGCGTGAAGTACACAAGTTGACAGTGCA'
#     success, rs = exact_match(seq, refseqs)
#     if not success:
#         print('Test exact match with multiple round brackets failed')
#         return False
    
#     refseqs = ['ATGCGTG(TTC)AAGTACACCCAAGTTGAC[A]GTGCA']
#     seq1 = 'ATGCGTGTTCAAGTACACCCAAGTTGACAGTGCA'
#     seq2 = 'ATGCGTGAAGTACACCCAAGTTGACAGTGCA'
#     seq3 = 'ATGCGTGTTCAAGTACACCCAAGTTGACGTGCA'
#     success1, rs1 = exact_match(seq1, refseqs)
#     success2, rs2 = exact_match(seq2, refseqs)
#     success3, rs3 = exact_match(seq3, refseqs)
#     if not success1:
#         print(f'Test exact match with square brackets 1 failed for {seq1}')
#         return False
#     if not success2:
#         print(f'Test exact match with square brackets 2 failed for {seq2}')
#         return False
#     if not success3:
#         print(f'Test exact match with square brackets 3 failed for {seq3}')
#         return False
#     return True


# def join_diffs(last_diffs):
#     """ 
#     Create a new diff string from a contiguous set of same-type differences
#     Positions should be zero-based
#     Each diff is of the form (i,ref_c,alt_c)
#     """
#     pos = last_diffs[0][0] + 1
#     ref_diff = ''.join([rc for pos,rc,ac in last_diffs])
#     alt_diff = ''.join([ac for pos,rc,ac in last_diffs])
#     if ref_diff.startswith('-'):  # insert in alternative
#         return ''.join([str(pos),'+',alt_diff])
#     elif alt_diff.startswith('-'):  # deletion in alternative
#         return ''.join([str(pos),'-',ref_diff])
#     # substitution
#     return ''.join([str(pos),ref_diff,'/',alt_diff])


# def annotate_seq(alnref, alnalt):
#     """
#     Return a string of PosRef/Alt e.g. 8-ATG11+TGC14A/G
#     Ignore leading and trailing unmatched positions (ref or alt ---), only for global alignments
#     TODO: handle brackets from variable sections in reference sequences
#     """
#     #print(f'Ref: {alnref}', flush=True)
#     #print(f'Alt: {alnalt}', flush=True)
#     lmargin = max(len(alnref)-len(alnref.lstrip('-')),len(alnalt)-len(alnalt.lstrip('-')))
#     rmargin = max(len(alnref)-len(alnref.rstrip('-')),len(alnalt)-len(alnalt.rstrip('-')))
#     #print(f'{lmargin=} {rmargin=}', flush=True)
#     all_diffs = []
#     last_diffs = []  # contiguous diffs of the same type: substitution, deletion, insertion
#     for i,(rc, ac) in enumerate(zip(alnref, alnalt)):
#         if i < lmargin:
#             continue
#         if i == len(alnref)-rmargin:
#             break
#         if rc != ac:
#             if rc == '-':  # insertion in alnalt
#                 if last_diffs:
#                     if last_diffs[-1][1] == '-':
#                         last_diffs.append((i,rc,ac))
#                     else:  # clear previously seen diffs
#                         all_diffs.append(join_diffs(last_diffs))
#                         last_diffs = [(i,rc,ac)]
#                 else:
#                     last_diffs = [(i,rc,ac)]
#             elif ac == '-':  # deletion in alnalt
#                 if last_diffs:
#                     if last_diffs[-1][2] == '-':
#                         last_diffs.append((i,rc,ac))
#                     else:  # clear previously seen diffs
#                         all_diffs.append(join_diffs(last_diffs))
#                         last_diffs = [(i,rc,ac)]
#                 else:
#                     last_diffs = [(i,rc,ac)]
#             else:  # substitution
#                 if last_diffs:
#                     if last_diffs[-1][1] != '-' and last_diffs[-1][2] != '-':
#                         last_diffs.append((i,rc,ac))
#                     else:
#                         all_diffs.append(join_diffs(last_diffs))
#                         last_diffs = [(i,rc,ac)]
#                 else:
#                     last_diffs = [(i,rc,ac)]
#         else:    
#             if last_diffs:
#                 all_diffs.append(join_diffs(last_diffs))
#                 last_diffs = []
#     if last_diffs:
#         all_diffs.append(join_diffs(last_diffs))

#     return ''.join(all_diffs)


# def test_annotation():
#     """
#     Some test cases to ensure that annotation of two aligned sequences is working correctly
#     """
#     refseq = 'ATTACGGTCT'
#     altseq = 'ATTACGGTCT'
#     ds = annotate_seq(refseq, altseq)
#     if ds != '':
#         return False
#     altseq = '--TACGGTCT'
#     ds = annotate_seq(refseq, altseq)
#     if ds != '':
#         return False
#     altseq = '--TACGGACT'
#     ds = annotate_seq(refseq, altseq)
#     if ds != '8T/A':
#         return False
#     altseq = '--TACGG-CT'
#     ds = annotate_seq(refseq, altseq)
#     if ds != '8-T':
#         return False
#     refseq = 'ATTACGGT-CT'
#     altseq = '--TACGGTTCT'
#     ds = annotate_seq(refseq, altseq)
#     if ds != '9+T':
#         return False
#     refseq = '--TACGGTCT--'
#     altseq = 'ATTACCTTCTAA'
#     ds = annotate_seq(refseq, altseq)
#     if ds != '6GG/CT':
#         return False
#     return True


# def inexact_match(seq, refseqs, run_dn, debug, lock_d, identity=0.90):
#     """
#     Inexact matching of a merged sequence with a set of reference target sequences
#     Declare a match if we have 90% sequence identity over the matched length or better
#     Identity must be >= 0 and <= 1.0
#     We must check for the best match against all targets or we will incorrectly call variant sequences!
#     Return True/False, and matching reference sequence/None
#     run_dn, debug, lock_d are used for logging debug messages
#     """
#     msg = f"debug: inexact match {seq=}"
#     wdb(msg, run_dn, debug, lock_d)
#     #print(msg, flush=True)
#     if len(refseqs) == 0:
#         wdb('Warning: No reference sequences to align to', run_dn, debug, lock_d)
#         return False, None, None

#     if identity < 0 or identity > 1.0:
#         return False, None, None
#     aligner = Align.PairwiseAligner()
#     aligner.mode = 'global'
#     aligner.open_gap_score = -2.0
#     aligner.extend_gap_score = -0.7
#     aligner.end_gap_score = -0.5
#     aligner.match_score = 3.5
#     aligner.mismatch_score = -1.0
#     #aligner.mode = 'local'
#     #aligner.open_gap_score = -1.0
#     #aligner.extend_gap_score = -1.0
#     #aligner.end_gap_score = -1.0
#     #aligner.match_score = 1
#     #aligner.mismatch_score = -1.0
#     #print('beginning alignment', flush=True)
#     try:
#         # do alignment against a reference sequence if they aren't too different in length
#         # order is align(Target, Query)
#         all_alignments = [(aligner.align(rs, seq)[0], rs) for rs in refseqs]
#     except Exception as exc:
#         msg = f'alignment failed with exception {exc}'
#         print(msg, flush=True)
#         wdb(msg, run_dn, debug, lock_d)
#         return False, None, None

#     if len(all_alignments) == 0:
#         wdb(f'No alignments found {seq=}', run_dn, debug, lock_d)
#         return False, None, None
#     #print(f"{len(all_alignments)} alignments completed", flush=True)
#     best_score = 0
#     best_algn = None
#     best_rs = None
#     for algn, rs in all_alignments:
#         if algn.score > best_score:
#             best_score = algn.score
#             best_algn = algn
#             best_rs = rs
#     if best_score < 1:
#         wdb(f'No alignments found {seq=}', run_dn, debug, lock_d)
#         return False, None, None
    
#     target_sections = []
#     query_sections = []
#     for segment in str(best_algn).split('\n'):
#         if segment.startswith('target'):
#             tokens = [s for s in segment.split(' ') if s != '']
#             if len(tokens) > 2:
#                 #print('target', tokens[2], flush=True)
#                 target_sections.append(tokens[2])
#         elif segment.startswith('query'):
#             tokens = [s for s in segment.split(' ') if s != '']
#             if len(tokens) > 2:
#                 #print('query', tokens[2], flush=True)
#                 query_sections.append(tokens[2])
#     aligned_rs = ''.join(target_sections)
#     aligned_seq = ''.join(query_sections)
#     adjusted_query_length = len(aligned_seq.lstrip('-').rstrip('-'))
#     adjusted_target_length = len(aligned_rs.lstrip('-').rstrip('-'))
#     match_length = min(adjusted_query_length, adjusted_target_length)

#     identities = best_algn.counts().identities
#     cutoff = match_length * identity
#     #print(f'{best_algn.score=} {best_rs=} {identities=} {cutoff=} {aligned_rs=} {aligned_seq=}', flush=True)
#     if identities > cutoff:
#         # Annotate match
#         #print(f'Matched! {best_rs=}', flush=True)
#         seq_anno = annotate_seq(aligned_rs, aligned_seq)
#         #print(f"{seq_anno=} {best_rs=}", flush=True)
#         return True, best_rs, seq_anno
#     msg = f'Best alignment identity worse than {identity*100}% {aligned_rs=} {aligned_seq=} {best_rs=} {seq=}'
#     wdb(msg, run_dn, debug, lock_d)
#     return False, None, None


# def clean_reads(run_dn, raw_dn, clean_dn, filename, bbduk_path):
#     """
#     Clean reads using bbduk
#     run_dn: directory where the raw reads are stored
#     filename: string name of raw reads file (R1)
#     Returns: True if successful, False otherwise
#     """
#     fn1 = Path(filename)
#     fn2 = Path(str(fn1).replace('_R1_001', '_R2_001'))
#     fn_raw1 = Path(run_dn/raw_dn/fn1)
#     fn_raw2 = Path(run_dn/raw_dn/fn2)                    
#     fn_cln1 = Path(run_dn/clean_dn/fn1)
#     fn_cln2 = Path(run_dn/clean_dn/fn2)
#     cmd = f"{bbduk_path} in1={fn_raw1} in2={fn_raw2} out1={fn_cln1} out2={fn_cln2} t=2 qtrim=rl trimq=20 minlen=50 k=23 ktrim=r mink=11 hdist=1 overwrite=true ref=bbmap/resources/adapters.fa"
    
#     try:
#         res1 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
#         if res1.returncode == 0:
#             return True
#         return False
#     except subprocess.CalledProcessError as e:
#         print(f"Error running bbduk: {e.stderr}", file=sys.stderr)
#         return False
    
        
# def merge_reads(run_dn, clean_dn, merge_dn, filename, bbmerge_path):
#     """
#     run BBMerge to join paired-end reads (from "cleaned")
#     run_dn: directory where the cleaned reads are stored
#     filename: string name of cleaned reads file (R1)
#     Returns: True if successful, False otherwise
#     """
#     fn1 = Path(filename)
#     fn2 = Path(str(fn1).replace('_R1_001', '_R2_001'))
#     fn_cln1 = Path(run_dn/clean_dn/fn1)
#     fn_cln2 = Path(run_dn/clean_dn/fn2)
#     fn_m = Path(run_dn/merge_dn/str(fn1).replace('_R1_001','_M'))
#     fn_u1 = Path(run_dn/merge_dn/str(fn1).replace('_R1_001', '_U1'))
#     fn_u2 = Path(run_dn/merge_dn/str(fn1).replace('_R1_001', '_U2'))
#     cmd = f"{bbmerge_path} in1={fn_cln1} in2={fn_cln2} out={fn_m} outu1={fn_u1} outu2={fn_u2} overwrite=true pfilter=1"
#     try:
#         res2 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
#         if res2.returncode != 0:
#             return False
#         return True
#     except Exception as exc:
#         print(f"Failed to run bbmerge: {exc}")
#         return False


# def run_preprocess(run_dn, raw_dn, clean_dn, merge_dn, filename, logq):
#     """
#     Run the preprocessing of reads in a separate thread
#     run_dn: directory where the run is stored
#     raw_dn: directory where the raw reads are stored
#     clean_dn: directory where the cleaned reads are stored
#     logq: queue for logging messages
#     """
#     seqcnt, passing, readCount, cleanCount, mergeCount, msg = preprocess_seqs(run_dn, raw_dn, clean_dn, logq)
#     if not passing:
#         logq.push(f"Preprocessing failed: {msg}")
#         return None
#     return seqcnt, readCount, cleanCount, mergeCount


# def launch_preprocess(run_dn, raw_dn, clean_dn, taskq, logq):
#     """
#     create a task for each sample FASTQ and place in the task queue
#     Call bbduk for cleaning of reads
#     Call bbmerge to merge read pairs
#     If number of merged reads is less than half the number of unmerged reads then call this a failure
#     args:
#     run_dn: directory where the run is stored
#     raw_dn: directory where the raw reads are stored
#     clean_dn: directory where the cleaned reads are stored
#     taskq: queue for tasks to be processed
#     logq: queue for logging messages
#     output:
#     Return seqcnt (Counter of unique sequences), passing (T/F), readCount, cleanCount, MergeCount
#     """
#     raw_dp = Path(run_dn)/raw_dn
#     fn1s = raw_dp.glob("*_R1_001.fastq.gz")
#     for fn1 in fn1s:
#         fn2 = str(fn1).replace('_R1_001', '_R2_001')
#         fn2 = Path(fn2)
#         if not fn2.is_file():
#             logq.push(f"Missing R2 file for {fn1}")
#             return {}, False, -1, -1, -1, "R2 file missing"
        
#         # Clean reads
#         if not clean_reads(run_dn, fn1.name, 'bbduk'):
#             with lock_l:
#                 log.append(f"Failed to clean reads for {fn1}")
#             wdb(f"Failed to clean reads for {fn1}", run_dn, debug, lock_d)
#             return {}, False, -1, -1, -1, "Cleaning failed"
        
#         # Merge reads
#         if not merge_reads(run_dn, fn1.name, 'bbmerge'):
#             with lock_l:
#                 log.append(f"Failed to merge reads for {fn1}")
#             wdb(f"Failed to merge reads for {fn1}", run_dn, debug, lock_d)
#             return {}, False, -1, -1, -1, "Merging failed"
#     readCount, cleanCount, cleanCount2, joinCount, mergeCount = -1, -1, -1, -1, -1
#     rfn= "*{}-{}_*_R1_*.fastq.gz".format(unguard(wr['pcrPlate'],silent=True), padwell(wr['pcrWell']))
#     fn1s = glob.glob(os.path.join(run_dn, "raw", rfn))
#     if not fn1s:
#         # we need to look for newer Illumina files with an underscore
#         rfn= "*{}_{}_*_R1_*.fastq.gz".format(unguard(wr['pcrPlate'],silent=True), padwell(wr['pcrWell']))
#         fn1s = glob.glob(os.path.join(run_dn, "raw", rfn))
#     fn1s = glob.glob(os.path.join(run_dn, "raw", rfn))
#     lrecs = []
#     if not fn1s:
#         with lock_l:
#             log.append(f"no data for {fn1s}")
#         wdb(f"no data for {fn1s}", run_dn, debug, lock_d)
#         return {}, False, readCount, cleanCount, mergeCount, "No files"
            
#     if len(fn1s)>1:
#         with lock_l:
#             log.append(f"too many reads files for pattern{rfn}")
#             log.append(f"   R1 files= {fn1s}")
#         wdb(f"too many reads files for pattern{rfn}\n    R1 files= {fn1s}", run_dn, debug, lock_d)
#         return {}, False, readCount, cleanCount, mergeCount, "Too many files"
            
#     fnr1 = fn1s[0]
#     fn1 = os.path.basename(fnr1) # just the one file
#     fn2 = fn1.replace('_R1_001', '_R2_001')
#     fnr2 = os.path.join(run_dn, 'raw', fn2) 
            
#     # find the data file
#     if not os.path.isfile(fnr2):
#         with lock_l:
#             log.append(f"missing file: {fnr2}")
#         wdb(f"missing file: {fnr2}", run_dn, debug, lock_d)
#         return dict(), False, readCount, cleanCount, mergeCount, "R2 file missing"
        
#     # use Windows file separators
#     fncs = tuple(os.path.join(run_dn,"cleaned", fn) for fn in (fn1, fn2))
#     bbmapd = os.path.join('bbmap','current')
#     # unpick the file name a bit
#     fnparts = fn1.split('_', 2)
#     fnmfmt = "{}-{}_{}{{}}.{{}}".format(wr['pcrPlate'], padwell(wr['pcrWell']), fnparts[1])
#     fnms = tuple(os.path.join(run_dn,"merged", fnmfmt.format('_'+tx, "fastq.gz")) for tx in ('M', 'U1', 'U2'))
#     fnlog = os.path.join(run_dn,"merged", fnmfmt.format('', 'log'))
#     if not all(os.path.isfile(fn) for fn in (fnms[0], fnlog)):
#         # java and bbmap need to be properly installed
#         cmd = r"java -ea -Xmx1g -cp {} jgi.BBDuk in1={} in2={} out1={} out2={} t=2 qtrim=rl trimq=20 minlen=50 k=23 ktrim=r mink=11 hdist=1 overwrite=true ref=bbmap\resources\adapters.fa".format(bbmapd, fnr1, fnr2, fncs[0], fncs[1])
#         try:
#             pres1 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
#         except Exception as exc:
#             msg = f"failed to run bbmerge {exc}".replace("'","").replace('"',"'").replace('\\','/')
#             return None, False, readCount, cleanCount, mergeCount, msg
#         if pres1.returncode != 0:
#             return None, False, readCount, cleanCount, mergeCount, "Failed: bbduk run failed"
#         # run BBMerge to join paired-end reads
#         # cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} verystrict=t".format(*((bbmapd,)+fncs+fnms)).split()
#         cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} overwrite=true pfilter=1".format(*((bbmapd,)+fncs+fnms)).split()
#         try:
#             pres2 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
#         except Exception as exc:
#             msg = f"failed to run bbmerge {exc}".replace("'","").replace('"',"'").replace('\\','/')
#             return None, False, readCount, cleanCount, mergeCount, msg
#         if pres2.returncode != 0:
#             return None, False, readCount, cleanCount, mergeCount, "Failed: bbmerge run failed"
#         # could delete cleaned data once it's been merged.
                
#         # keep the log output as record counts get used
#         with open(fnlog, "wt") as logfile:
#             for s in (pres1.stdout, pres1.stderr):
#                 if s:
#                     print(s, file=logfile)
#             print('======', file=logfile)
#             for s in (pres2.stdout, pres2.stderr):
#                 if s:
#                     print(s, file=logfile)
#         log1, log2 = ['\n'.join((p.stdout, p.stderr)) for p in (pres1, pres2)]

#     else:
#         fnm_fmt = fnmfmt.replace('{}.{}','')
#         lrecs.append(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}")
#         wdb(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}", run_dn, debug, lock_d)
#         with open(fnlog) as logfile:
#             logdata = logfile.read()
#         log1, log2 = logdata.split("======\n", 1)
            
#     # get counts from BBDuk & BBmerge output
    
#     m = re.search(r'Input:\s+(\d+) reads', log1)
#     if m:
#         readCount = int(m.group(1))//2
#     m = re.search(r'Result:\s+(\d+) reads', log1) 
#     if m:
#         cleanCount = int(m.group(1))//2
#     m = re.search(r'Pairs:\s+(\d+)\s', log2)
#     if m:
#         cleanCount2 = int(m.group(1))
#     m = re.search(r'Joined:\s+(\d+)\s', log2)
#     if m:
#         joinCount = int(m.group(1))
                
#     if cleanCount!=cleanCount2:
#         msg = f"Pair count mismatch. {cleanCount} != {cleanCount2}"
#         with lock_l:
#             log.append(msg)
#         wdb(msg, run_dn, debug, lock_d)
#         return {}, False, readCount, cleanCount, joinCount, msg
                
#     fn = fnms[0] # merged reads file
#     with myopen(fn) as src:
#         # read a FASTQ file - sequence in in second line of each 4 line
#         gs = (r.strip() for i, r in enumerate(src) if i%4==1)
#         seqcnt = collections.Counter(gs)
#     mergeCount = sum(seqcnt.values())
#     if mergeCount != joinCount:
#         with lock_l:
#             msg = f"Merged counts mismatch. {mergeCount} != {joinCount}"
#             log.append(msg)
#             wdb(msg, run_dn, debug, lock_d)
#         return {}, False, readCount, cleanCount, mergeCount, msg

#     with lock_l:
#         for l in lrecs:
#             log.append(l)
    
#     return seqcnt, True, readCount, cleanCount, mergeCount, None


# def wdb(msg, run_dn, debug, lock_d=None):
#     """ 
#     wdb = write debug.
#     debug [T/F] - do this here to reduce shared logic
#     lock_d: Call from within a lock 
#     """
#     #print (f'In wdb {msg=}')
#     if debug:
#         debugfn = os.path.join(run_dn, 'debug.log')
#         if lock_d:
#             with lock_d:
#                 try:
#                     with open(debugfn, 'at') as df:
#                         print(msg, file=df, flush=True)
#                 except:
#                     print('Could not write to debug')
#                     print(f'{msg}', flush=True)
#         else:
#             try:
#                 with open(debugfn, 'at') as df:
#                     print(msg, file=df, flush=True)
#             except:
#                 print('Could not write to debug', flush=True)
#                 print(f'{msg}')
#     #print('leaving wdb')


# def archetypes(seq_counts, run_dn, debug, lock_d):
#     """ 
#     Collapse substrings and superstrings
#     inputs: dictionary of sequences and counts
#     output: collections.Counter
#     """
#     start = timer()
#     final_seqs = {}
#     unique_seqs = set(seq_counts.keys())
#     msg = f'Start of archetypes. Number of sequences {len(seq_counts)} total counts {sum([seq_counts[s] for s in seq_counts])}'
#     wdb(msg, run_dn, debug, lock_d)
    
#     unique_seqs = sorted(unique_seqs, key=len)
#     skip_indices = set()
#     kept_seqs = []
#     for i, u1 in enumerate(unique_seqs):
#         if i in skip_indices:
#             continue
#         kept_seqs.append(u1)
#         skip_indices.add(i)
#         for j, u2 in enumerate(unique_seqs):
#             if j <= i:
#                 continue
#             if u2 in u1:
#                 skip_indices.add(j)
#                 seq_counts[u1] += seq_counts[u2]
#     for ks in kept_seqs:
#         final_seqs[ks] = seq_counts[ks]

#     #print(f'Final {collections.Counter(final_seqs).most_common()=}')
#     msg = f'End of archetypes. Number of sequences {len(final_seqs)} total counts {sum([final_seqs[s] for s in final_seqs])}'
#     wdb(msg, run_dn, debug, lock_d)
   
#     end = timer()
#     wdb(f'Archetypes took {end-start}', run_dn, debug, lock_d)
    
#     return collections.Counter(final_seqs)


# def process_well(work_block, wr, run_dn, seq_ids, id_seq, primer_assayfam, assayfam_primers, 
#         match_cache, anno_cache, miss_cache, reps, log, lock_mtc, lock_msc, lock_r, 
#         lock_l, lock_d, margin, identity, mincov, minprop, inexact, no_miss_cache,
#         exhaustive, debug=False):
#     """ 
#     Worker process that runs bbduk, bbmerge, and matching
#     work_block: work block number
#     wr: well record
#     seq_ids:  {seq:set([ids]}
#     id_seq: {seq_id:seq}
#     primer_assayfam: {primer:assayfam}
#     assayfam_primers: {assayfam:[primers]}
#     match_cache: dict of observed sequence to known sequence
#     anno_cache: dict of annotations for variant sequences
#     miss_cache: dict of unknown sequences, used as a set
#     reps: dict of sampleNo to all output fields as a list {sampleNo:list} 
#     lock_mtc: lock object for the match cache
#     lock_msc: lock object for the miss cache
#     lock_r: lock object for the results    
#     lock_l: lock for logging
#     lock_d: lock for writing to debug
#     margin: minimum proportion overlap for exact matching
#     identity: proportional identity score for inexact matching
#     inexact: enable inexact matching
#     no_offtarget_inexact: disable offtarget inexact matching
#     no_miss_cache: disable miss cache
#     exhaustive: [T/F] don't skip any sequences, not matter how low their count
#     debug: write messages direct to file (slow I/O)

#     Parallel, independent processing of a single well.
#     Clean and merge reads using bbduk
#     Group sequences into archetypes and order by decreasing count
#     Check whether any sequences exactly match the expected assay family and base cutoffs on this
#     Check all remaining sequences until counts falls below cutoff
#       -> Exact matches against prefered assays
#       -> Exact matches against all assay families
#       -> Exact matches against the miss cache
#       The miss cache is problematic... we may miss amplicons that are out of context
#       -> Inexact matches against the match cache -> add variants to match cache
#       -> Inexact matches against the miss cache -> add variants to miss cache
#       -> Add any remaining entries to the miss cache

#     Match cache ["ATC..."]="ATC..." <- a sequence in seq_ids and id_seq
#     Writes back result (seqCount,seqName,Efficiency,otherCount,otherName) to result thread
#     """
#     PID = os.getpid()
#     lrecs = [] # log records, lock and write on exit/return
#     sampleNo = wr.get('sampleNo', None)
#     pcrPlate = wr.get('pcrPlate', None)
#     pcrWell = wr.get('pcrWell', None)
#     if sampleNo is None or pcrPlate is None or pcrWell is None:
#         msg = f"Critical: {sampleNo=}, {pcrPlate=} or {pcrWell=} missing from Stage3.csv!"
#         wdb(msg, run_dn, debug, lock_d)
#         retval = [str(wr[x]) for x in wr] + [-1,-1,-1, '', '', -1, msg]
#         with lock_r:
#             reps[sampleNo] = retval
#         with lock_l:
#             for l in lrecs:
#                 log.append(l)
#         return
#     pcrPlate = unguard(pcrPlate, silent=True)
#     pcrWell = padwell(pcrWell)
#     msg = f"Info: {PID} Working on {sampleNo=} {pcrPlate=} {pcrWell=}"
#     lrecs.append(msg)
#     wdb(msg, run_dn, debug, lock_d)
#     if debug:
#         print(msg, flush=True)

#     # clean and merge FASTQs
#     seqcnt, success, readCount, cleanCount, mergeCount, fault_msg = preprocess_seqs(wr, run_dn, log, 
#             lock_l, lock_d, debug=debug)
#     if not success:
#         msg = f"Failed: preprocessing for {pcrPlate} {pcrWell} with {fault_msg}"
#         if debug:
#             print(msg, flush=True)
#             print(f'{readCount=} {cleanCount=} {mergeCount=} {fault_msg=}', flush=True)
#         lrecs.append(msg)
#         wdb(msg, run_dn, debug, lock_d)
        
#         retval = [str(wr[x]) for x in wr] + [readCount, cleanCount, mergeCount, '','',-1,'Failed to run bbduk']
        
#         wdb(f'Failed {retval=}', run_dn, debug, lock_d)
#         if debug:
#             print(f'Failed {retval=}', flush=True)
#         with lock_r:
#             reps[int(sampleNo)] = retval
#         with lock_l:
#             for l in lrecs:
#                 log.append(l)
#         return

#     msg = f"Info: Completed preprocessing for {pcrPlate} {pcrWell} "+\
#             f"{readCount=} {cleanCount=} {mergeCount=}"
#     lrecs.append(msg)
#     wdb(msg, run_dn, debug, lock_d)
                
#     # if we get an unknown primer then we should run archetypes and call it a day?
#     primer = wr.get('primer', 'No_primer')
#     if primer not in primer_assayfam:
#         msg = f'Warning: primer {primer} unknown'
#         lrecs.append(msg)
#         wdb(msg, run_dn, debug, lock_d)
        
#     # decide which assays are on-target vs off-target
#     on_target_ids = set()
#     off_target_ids = set()
#     on_target_seqs = set()
#     off_target_seqs = set()
#     for name in id_seq:
#         if name.lower().startswith(primer.split('_')[0].lower() + '_'):
#             on_target_ids.add(name)
#             on_target_seqs.add(id_seq[name])
#         else:
#             off_target_ids.add(name)
#     # avoid looking at the same sequence under different names - used by best_match
#     off_target_seqs = set(list(seq_ids.keys())).difference(on_target_seqs)
#     #print(f"{wr['pcrWell']} {on_target_ids} {on_target_seqs}", flush=True)
#     msg = f"Debug: {on_target_ids=} {len(off_target_seqs)=}"
#     wdb(msg, run_dn, debug, lock_d)
#     match_cnt = Counter()
#     if mincov < 1:
#         mincov = 1
#     if minprop > 1.0:
#         minprop = 1.0
#     elif minprop < 0.0:
#         minprop = 0.0
        
#     # unique sequences only (substrings collapsed), from most to least common
#     if exhaustive:
#         wdb('Info: Aggregating archetype sequences', run_dn, debug, lock_d)
#         seqcnt = archetypes(seqcnt, run_dn, debug, lock_d)
        
#     # calculate min count proportion from exact matches to our expected targets
#     family_exact_counts = 0
#     for on_target_id in on_target_ids:
#         on_target_seq = id_seq[on_target_id]
#         if on_target_seq in seqcnt:
#             family_exact_counts += seqcnt[on_target_seq]
#     low_cov_cutoff = max(family_exact_counts*minprop, mincov)

#     other_count = 0
#     mtc = dict(match_cache)  # force copy of shared object
#     original_match_cache_size = len(mtc)
#     anc = dict(anno_cache)  # force copy of shared object
#     original_anno_cache_size = len(anc)
#     msc = dict(miss_cache)  # force copy of shared object
#     original_miss_cache_size = len(msc)

#     for seq, num in seqcnt.most_common():
#         if (num >= low_cov_cutoff) or exhaustive:
#             msg = f"Info: Processing {pcrPlate} {pcrWell} on process {PID} with counts {num} and {seq}"
#             wdb(msg, run_dn, debug, lock_d)      
#             if seq in mtc:
#                 msg = f"Info: Match cache hit {num=} {seq=} {mtc[seq]=}"
#                 wdb(msg, run_dn, debug, lock_d)
#                 ref_seq = mtc[seq]
#                 match_cnt[ref_seq] += num
#                 if seq in anc:
#                     anno = anc[seq]
#                     match_cnt[ref_seq+'//'+anno] += num
#                 continue
#             # miss cache - doesn't match any known sequence
#             if not no_miss_cache and seq in msc:
#                 msg = f"Info: Miss cache hit {num} {seq}"
#                 wdb(msg, run_dn, debug, lock_d)
#                 match_cnt[seq] += num
#                 continue
#             # substring or superstring
#             is_match, ref_seq = exact_match(seq, on_target_seqs, margin=margin)
#             if is_match:
#                 msg = f"Info: Exact match against {seq_ids[ref_seq]} with {seq} and counts {num}"
#                 wdb(msg, run_dn, debug, lock_d)
#                 match_cnt[ref_seq] += num
#                 continue
#             is_match, ref_seq = exact_match(seq, off_target_seqs, margin=margin)
#             if is_match:
#                 msg = f"Info: Exact match against {seq_ids[ref_seq]} with {seq} and counts {num}"
#                 wdb(msg, run_dn, debug, lock_d)
#                 match_cnt[ref_seq] += num
#                 continue
            
#             if  inexact:
#                 # inexact matching must be done against all known sequences at once or it risks false association
#                 is_match, ref_seq, seq_anno = inexact_match(seq, on_target_seqs.union(off_target_seqs),run_dn, debug, 
#                         lock_d, identity=identity)

#                 if is_match:
#                     msg = f"Info: Inexact match against {seq_ids[ref_seq]} {wr['pcrPlate']}"+\
#                             f" {wr['pcrWell']} {num} {seq}\n"
#                     wdb(msg, run_dn, debug, lock_d)
#                     mtc[seq] = ref_seq
#                     match_cnt[ref_seq] += num
#                     if seq_anno:
#                         anc[seq] = seq_anno
#                         match_cnt[ref_seq+'//'+seq_anno] += num
#                     continue
            
#             # add to miss cache
#             msg = f"Info: No match for {pcrPlate=} {pcrWell=} {num=} {seq=}\n"
#             wdb(msg, run_dn, debug, lock_d)
#             if not no_miss_cache:
#                 msc[seq] = None
#             match_cnt[seq] += num
                    
#         else:
#             #msg = f"Debug: too few reads {num} for matching {seq}"
#             #wdb(msg, run_dn, debug, lock_d)
#             match_cnt['other'] += num
#             other_count += 1
#         # try another sequence

#     if len(mtc) > original_match_cache_size:
#         with lock_mtc:
#             match_cache.update(mtc) 
#             anno_cache.update(anc)
#         wdb('Info: Updating match and annotation caches', run_dn, debug, lock_d)
           
#     if len(msc) > original_miss_cache_size:
#         with lock_msc:
#             miss_cache.update(msc)
#         wdb('Info: Updating miss cache', run_dn, debug, lock_d)

#     msg = f"Info: Completed matching for {PID=} {pcrPlate=} {pcrWell=}"
#     lrecs.append(msg)
#     wdb(msg, run_dn, debug, lock_d)
#     # rename 'other' to include number of separate sequence variations
#     total_others = match_cnt['other']
#     other_name = 'other (' + str(other_count) +')'
#     match_cnt[other_name] = total_others
#     res1 = [readCount, cleanCount, mergeCount]
#     # name the outputs and counts
#     seqCounts = []
#     seqNames = []
#     otherCounts = []
#     otherNames = []

#     for seq, count in match_cnt.most_common():
#         if seq == 'other':
#             continue
#         if (seq.startswith('other') and seq != 'other') or seq in msc:
#             otherNames.append(seq)
#             otherCounts.append(count)
#             continue

#         if '//' not in seq:
#             if seq in mtc:
#                 refseq = mtc[seq]
#                 if refseq in on_target_seqs:
#                     seqCounts.append(count)
#                     refname = [name for name in seq_ids[refseq] if name in on_target_ids][0]
#                     seqNames.append(refname)
#                     if seq in anc:
#                         otherCounts.append(count)
#                         otherNames.append(refname+'//'+anc[seq])
#                 else:
#                     for name in seq_ids[refseq]:
#                         otherCounts.append(count)
#                         otherNames.append(name)
#                 continue
#             else:
#                 msg = f'Critical: {seq=} not found in match cache {mtc=}'
#                 print(msg, flush=True)
#                 wdb(msg, run_dn, debug, lock_d)
#         else:
#             ref_seq = seq.split('//')[0]
#             anno = seq.split('//')[1]
#             if ref_seq not in seq_ids:
#                 msg = f'Critical: {ref_seq=} not found in sequence to ID mapping {seq_ids=}'
#                 print(msg, flush=True)
#                 wdb(msg, run_dn, debug, lock_d)
#             else:
#                 for name in seq_ids[ref_seq]:
#                     otherCounts.append(count)
#                     otherNames.append(name+'//'+anno)
#     # calculate efficiency of PCR as the sum(seqCounts)/wr['mergeCount']
#     try:
#         efficiency = sum(map(int,seqCounts))/mergeCount
#     except Exception as exc:
#         wdb(f"Error: calculating efficiency {seqCounts=} {mergeCount=} {exc=}", run_dn, debug, lock_d)
#         effiency = -1.0
#     # combine outputs
#     try:
#         res2 = [';'.join(map(str,seqCounts)), ';'.join(map(str,seqNames)), f'{round(efficiency,3)}', 
#                 ';'.join(map(str,otherCounts)), ';'.join(map(str,otherNames))]
#     except Exception as exc:
#         wdb(f"Error: joining names and counts {exc=}", run_dn, debug, lock_d)
#         res2 = ['','','','','']
#     retval = [str(wr[x]) for x in wr] + res1 + res2
#     sn = None
#     try:
#         sn = int(sampleNo)
#     except Exception as exc:
#         msg = f"Critical: sampleNo not an integer {exc=}"
#         print(msg, flush=True)
#         wdb(msg, run_dn, debug, lock_d)
#     if sn:    
#         wdb(f'{retval=}', run_dn, debug, lock_d)    
#         with lock_r:
#             reps[sn] = retval
#     wdb(f'{retval=}', run_dn, debug, lock_d)
#     with lock_l:
#         for l in lrecs:
#             log.append(l)
#     msg = f"Info: Exiting process_well() {PID=} {sampleNo=} {pcrPlate=} {pcrWell=}"
#     wdb(msg, run_dn, debug, lock_d)
#     if debug:
#         print(msg)


# def write_log(log, logfn):
#     print(f"Writing log to {logfn}", file=sys.stderr)
#     with open(logfn, 'wt') as logf:
#         for l in log:
#             print(l, file=logf)

 
# def get_raw_fastq_pairs(dirpath):  
#     """ return a sorted list of tuple(R1_path, R2_path) to raw FASTQ files """
#     valid_pairs = []
#     r1s = [dirpath/Path(f) for f in os.listdir(dirpath) if f.endswith('.fastq.gz') and '_R1_' in f]
#     for r1 in r1s:
#         r2 = Path(str(r1).replace('_R1_001','_R2_001'))
#         if r1.is_file() and r2.is_file():
#             valid_pairs.append((r1,r2))
#     return sorted(valid_pairs)


# def report_progress(run_dn, launch_progress, match_progress):
#     """
#     Clear all previous progress files and touch a new file with the launch and match progress values in the name
#     Only do the operation is progress is a multiple of 5 to save on disk writes
#     """
#     progress_files = list(Path(run_dn).glob('match_progress_*'))
#     if len(progress_files) == 1:
#         lp = int(str(progress_files[0]).split('_')[-2])
#         mp = int(str(progress_files[0]).split('_')[-1])
#         if lp == launch_progress and mp == match_progress:
#             # no change
#             return
#     for fn in list(Path(run_dn).glob('match_progress_*')):
#         os.remove(fn)
#     progress_fn = os.path.join(args.run_dn,'match_progress_'+str(launch_progress)+'_'+str(match_progress))
#     with open(progress_fn, 'wt') as fp:
#         pass  # just touch the file


# def parse_targets(run_dn, targets, log, debug):
#     """
#     Return a dictionary of sequences to a set of ids, and
#     a dictionary of ids to sequences
#     """
#     seq_ids = {}
#     id_seq = {}
#     with myopen(os.path.join(run_dn,targets)) as src:
#         for entry in Bio.SeqIO.parse(src, "fasta"):
#             if len(entry.id) > 0 and len(entry.seq) > 0:
#                 try:
#                     name = bytes(str(entry.id).strip(), 'ascii', errors='strict').decode()
#                 except Exception as exc:
#                     msg = f'Warning: {exc} non-ascii character in reference sequence id {entry.id}'
#                     log.append(msg)
#                     wdb(msg, args.run_dn, debug)
#                     seq_id = bytes(str(entry.id).strip(), 'ascii',errors='ignore').decode()
#                 try:
#                     seq = bytes(str(entry.seq).strip(), 'ascii', errors='strict').decode().upper()
#                 except Exception as exc:
#                     msg = f'Warning: {exc} non-ascii character in reference sequence {entry.seq}'
#                     log.append(msg)
#                     wdb(msg, args.run_dn, debug)
#                     seq = bytes(str(entry.seq).strip(), 'ascii', errors='ignore').decode().upper()
#                 if seq not in seq_ids:
#                     seq_ids[seq] = set()
#                 seq_ids[seq].add(name)
#                 if name in id_seq:
#                     msg = f'Warning: {name} duplicated in reference'
#                     log.append(msg)
#                     wdb(msg, args.run_dn, debug)
#                 id_seq[name] = seq
#     return seq_ids, id_seq


# def parse_primer_assayfams(run_dn, paf):
#     """
#     Two column file (tab separated) of primer and assay families
#     Return a mapping of primer to set of assay family, and mapping of assay family to set of primers
#     """
#     primer_assayfam = {}  # [primer] = []
#     assayfam_primers = {}  # [assayfam] = []
#     pmr_fn = os.path.join(run_dn, paf)
#     try:
#         with open(pmr_fn, 'rt', errors='ignore') as f:
#             for line in f:
#                 cols = line.strip().split(',')
#                 if len(cols) != 2:
#                     print(f'skipping line {line=}')
#                     continue
#                 primer = cols[0]
#                 assayfam = cols[1]
#                 if not primer or not assayfam:
#                     continue
#                 if primer not in primer_assayfam:
#                     primer_assayfam[primer] = set()
#                 primer_assayfam[primer].add(assayfam)
#                 if assayfam not in assayfam_primers:
#                     assayfam_primers[assayfam] = set()
#                 assayfam_primers[assayfam].add(primer)
#         return primer_assayfam, assayfam_primers
#     except Exception as exc:
#         msg = f'Critical: failed to open {pmr_fn} for primer assay family info {exc=}'
#         print(msg)
#     return primer_assayfam, assayfam_primers


# def parse_primer_file(paf):
#     """
#     Parse the user provided assay list file. Use this if testing without a primer_assayfam file
#     """
#     primer_assayfam = {}
#     assayfam_primers = {}  # [assayfam] = []
#     with open(paf, 'rt') as f:
#         for line in f:
#             cols = [c.strip() for c in line.split(',')]
#             assayfam = cols[1]
#             pmrs = [c for c in cols[2:] if len(c) != 0]
#             if assayfam not in assayfam_primers:
#                 assayfam_primers[assayfam] = set()
#             for pmr in pmrs:
#                 if pmr not in primer_assayfam:
#                     primer_assayfam[pmr] = set()
#                 primer_assayfam[pmr].add(assayfam)
#                 assayfam_primers[assayfam].add(pmr)
#     return primer_assayfam, assayfam_primers


# def get_variant_seq(var_name, id_seq):
#     """
#     Reverse engineer the variant sequence from annotations in the var name
#     Positions are 1-based
#     """
#     if '//' not in var_name:
#         return var_name  # We should never see this happen, but to be safe...
#     primer = var_name.split('//')[0]
#     if primer not in id_seq:
#         print(f'Error: no known primer: {primer}')
#         return ''
#     original_seq = id_seq[primer]
#     parts = re.split(r'(\d+)', var_name.split('//')[1])
#     rev_parts = parts[::-1]
#     new_seq = original_seq
#     # need to go through changes in reverse order to avoid length changes from affecting position
#     for i, p in enumerate(rev_parts):
#         if i % 2 == 0:
#             if p == '':
#                 break  # we've reached the end
#             change = p
#             if len(parts) <= i+1:
#                 print(f'Error: no matching change for position: {pos} in {parts} from {var_name}', flush=True)
#                 return ''
#             try:
#                 pos = int(rev_parts[i+1]) -1  # 1-based
#             except Exception as e:
#                 print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
#                 return ''
#             #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
#             if '+' in change:
#                 new_seq = new_seq[:pos] + change[1:] + new_seq[pos:]
#             if '-' in change:
#                 new_seq = new_seq[:pos] + new_seq[pos+len(change)-1:]
#             if '/' in change:
#                 repl_bases = change.split('/')[1] 
#                 new_seq = new_seq[:pos] + repl_bases +new_seq[pos+len(repl_bases):]
            
#     #print(f'{original_seq=}', flush=True)
#     #print(f'{new_seq=}', flush=True)
#     return new_seq


# def test_variant_seq():
#     """
#     unit tests for get_variant_seq()
#     """
#     print('beginning unit tests for variant sequence recreation', flush=True)
#     test_id_seq = {'orig':'ACTGAACCTTGG'}
#     test1 = 'orig//4+C'
#     expected = 'ACTCGAACCTTGG'
#     new_seq = get_variant_seq(test1, test_id_seq)
#     if expected != new_seq:
#         print(f'Test 1: {expected} does not match {new_seq}!', flush=True)
#     else:
#         print('Test 1: pass', flush=True)
        
#     test2 = 'orig//4-G'
#     expected = 'ACTAACCTTGG'
#     new_seq = get_variant_seq(test2, test_id_seq)
#     if expected != new_seq:
#         print(f'Test 2: {expected} does not match {new_seq}!', flush=True)
#     else:
#         print('Test 2: pass', flush=True)
        
#     test3 = 'orig//4G/C'
#     expected = 'ACTCAACCTTGG'
#     new_seq = get_variant_seq(test3, test_id_seq)
#     if expected != new_seq:
#         print(f'Test 3: {expected} does not match {new_seq}!', flush=True)
#     else:
#         print('Test 3: pass', flush=True)
        

# # def main(run_dn, stagefile, logfn, outfn, targets, exhaustive=False, debugfn='debug.log', debug=False):
# def main(args):
#     """
#     Read background data: target reference sequences file
#     Then processes merged pairs files producing NGS report.
#     """
#     log = []
#     if args.debug:
#         db_log_fn = os.path.join(args.run_dn, 'debug.log')
#         if Path(db_log_fn).exists():
#             try:
#                 os.remove(db_log_fn)
#             except Exception as exc:
#                 print('Could not clear debug.log, perhaps you have it open?', file=sys.stderr)
#                 return
#         print(f'Writing debug information to {db_log_fn}', file=sys.stderr)
#     try:
#     #if True:  # helps with debugging
#         log.append('Info: Run with the following command line options:')
#         for arg in vars(args):
#             log.append(f'Info: {arg} {getattr(args, arg)}')
        
#         # read the wells data
#         start_time = datetime.datetime.now()
#         log.append(f"Begin: {start_time}")
#         raw_pair_list = get_raw_fastq_pairs(os.path.join(args.run_dn, 'raw'))
#         raw_file_identifiers = [str(f[0].name).split('_')[0] for f in raw_pair_list]
#         if len(raw_file_identifiers) == 0:
#             return
#         # New Illumina doesn't have a dash in the file name
#         if raw_file_identifiers[0][-4] != '-':
#             raw_file_identifiers = ['-'.join(str(f[0].name).split('_')[0:2]) for f in raw_pair_list]
#         raw_file_identifiers = set(raw_file_identifiers)
#         #print(sorted(raw_file_identifiers), file=sys.stderr)
#         with open(os.path.join(args.run_dn,args.stagefile)) as srcfd:
#             src = csv.reader(srcfd, dialect="unix")
#             hdr = next(src)
#             WRec = collections.namedtuple("WRec", hdr)
#             wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrPlate, x.pcrWell[0], int(x.pcrWell[1:]), x.primer))
#         #print(wdata, file=sys.stderr)
#         wdata = [rec for rec in wdata if unguard(rec.pcrPlate, silent=True) +'-'+ padwell(rec.pcrWell) in raw_file_identifiers]
#         #print(f'After filtering by available files {wdata=}', file=sys.stderr)
#         log.append(f"Info: {len(wdata)} sample wells to process.")
#         if len(wdata) == 0:
#             return
#         ## get a set of assay family names - family name ends with first underscore char  
#         #assays = frozenset(r.primer.split('_',1)[0] for r in wdata)
            
#         # read the target sequence file into dictionaries [seq] = set([ids]), and [id] = seq             
#         seq_ids, id_seq = parse_targets(args.run_dn, args.targets, log, args.debug)
#         print('Parsed targets', flush=True)

#         # get relationship of primers to assay family
#         #print(f'{args.primer_assayfam=}')
#         primer_assayfam, assayfam_primers = parse_primer_assayfams(args.run_dn, args.primer_assayfam)
#         #primer_assayfam, assayfam_primers = parse_primer_file(os.path.join(args.run_dn, args.primer_assayfam))
#         print(f'Parsed primers/assays', flush=True) # {primer_assayfam=} {assayfam_primers=}')
        
#         if not os.path.isdir(os.path.join(args.run_dn,"raw")):
#             log.append("Error: raw FASTQ data folder is absent - please transfer MiSeq data.")
#             write_log(log, os.path.join(args.run_dn,args.logfn))
#             return
        
#         for d in ("cleaned", "merged"):
#             dp = os.path.join(args.run_dn, d)
#             if not os.path.isdir(dp):
#                 os.mkdir(dp)
                    
#         # parallel execute over wdata and collect results
#         grouped_wrs = itertools.groupby(wdata, key=lambda x:(x.pcrPlate, x.pcrWell))
#         wrs = []
#         for key, group in grouped_wrs:
#             for g in group:
#                 wrs.append(OrderedDict(g._asdict()))
            
#         with mp.Manager() as manager:
#             # Use a server process to manage shared data structures. Heavy option, but works on Windows
#             match_cache = manager.dict()  # [seq] = known_seq
#             anno_cache = manager.dict()  # [seq] = annotation
#             miss_cache = manager.dict()  # just use this as a set
#             reps = manager.dict()  # [sampleNo] = all columns
#             logm = manager.list()
#             lock_mtc = manager.Lock()  # match_cache locking
#             lock_msc = manager.Lock()  # miss_cache locking
#             lock_r = manager.Lock()  # result list locking
#             lock_l = manager.Lock()  # log list locking
#             lock_d = manager.Lock()  # debug file locking
        
#             with lock_mtc:
#                 for seq in seq_ids:
#                     match_cache[seq] = seq  # variants of these will also go in here 

#             #print('Before launching jobs', file=sys.stderr, flush=True)
#             # multiprocessing
#             NUMPROCS = args.ncpus
            
#             pool = manager.Pool(NUMPROCS)
#             reports = []
#             launch_progress = 0
#             match_progress = 0
#             print('launching jobs', flush=True)
#             total_jobs = len(wrs)
#             # multiprocessing pool counter from https://superfastpython.com/multiprocessing-pool-asyncresult/
#             reports = []
#             for i, wr in enumerate(wrs):
#                 r = pool.apply_async(process_well, args=(i, wr, args.run_dn, seq_ids, id_seq, primer_assayfam, 
#                     assayfam_primers, match_cache, anno_cache, miss_cache, reps, logm, lock_mtc, 
#                     lock_msc, lock_r, lock_l, lock_d, args.margin, args.identity, args.mincov, args.minprop, 
#                     args.inexact, args.no_miss_cache, args.exhaustive, args.debug))
#                 reports.append(r)
#                 if i % 3 == 0:
#                     launch_progress = ceil(100*i/total_jobs)
#                     completed = sum([r.ready() for r in reports])
#                     match_progress = floor(100*completed/total_jobs)
#                     report_progress(args.run_dn, launch_progress, match_progress)
                    
#             while match_progress < 100:
#                 completed = sum([r.ready() for r in reports])
#                 # report the number of remaining tasks
#                 print(f'Match completion: {100*completed/total_jobs}')
#                 match_progress = 100*completed/total_jobs
#                 report_progress(args.run_dn, 100, floor(match_progress))
#                 # wait a moment
#                 time.sleep(2.5)       
#             report_progress(args.run_dn, 100, 100)
                            
#             pool.close()
#             pool.join()
#             print('All processes completed', flush=True)

#             with lock_l:
#                 for l in logm:
#                     log.append(l)

#             with lock_r:
#                 completed_jobs = [reps[key] for key in sorted(reps.keys())]

#             test_variant_seq()

#             with open(os.path.join(args.run_dn,args.outfn), "wt", buffering=1) as dstfd,\
#                     open(os.path.join(args.run_dn,args.variants), 'wt', buffering=1) as varfd:
#                 print(f"Opening {args.outfn} for results", flush=True)
#                 print(f"Opening {args.variants} for variant sequences", flush=True)

#             #with open(os.path.join(args.run_dn,args.outfn), "wt", buffering=1) as dstfd:
#             #    print(f"Opening {args.outfn} for results", file=sys.stderr)
#                 dst = csv.writer(dstfd, dialect="unix", quoting=csv.QUOTE_ALL)
#                 hdrres1 = ("readCount", "cleanCount", "mergeCount")
#                 hdrres2 = ("seqCount", "seqName", "efficiency", "otherCount", "otherName")
#                 complete_row_hdr = tuple((x for xs in (hdr, hdrres1, hdrres2) for x in xs))
#                 primer_col = [i for i,col in enumerate(complete_row_hdr) if col=='primer'][0]
#                 dst.writerow(complete_row_hdr)
#                 for i,job in enumerate(completed_jobs):
#                     try:
#                         dst.writerow(job)
#                     except Exception as exc:
#                         print(f'{exc=}', flush=True)
#                     #print(i, job, flush=True)
#                     var_count_entries = job[-2].split(';')
#                     var_name_entries = job[-1].split(';')
#                     primer_name = job[primer_col]
#                     #print(f'{primer_name=} {var_count_entries=} {var_name_entries=}')
#                     for var_count, var_name in zip(var_count_entries, var_name_entries):
#                         if var_name.startswith('other'):
#                             continue
#                         var_row_name = f'>Sample:{i+1};Primer:{primer_name}'
#                         if '//' in var_name:
#                             var_row_name += f';{var_name}'
#                         else:
#                             var_row_name += f';No_match'
#                         var_row_name += f';count:{var_count}'
#                         print(var_row_name, file=varfd)
#                         if '//' in var_name:
#                             print(get_variant_seq(var_name, id_seq), file=varfd)
#                         else:
#                             print(var_name, file=varfd)
                            
#                 dstfd.flush()
#                 varfd.flush()
#     except Exception as exc:
#         print(f'Error: {exc}', flush=True)
#         log.append(f'Error: {exc}')

#     end_time = datetime.datetime.now()
#     msg = f"End: {end_time} took: {end_time - start_time}"
#     print(msg, flush=True)
#     log.append(msg)
#     write_log(log, os.path.join(args.run_dn,args.logfn))
     
    
# def run_matches(exp, ncpus, mincov, minprop, exhaustive, debug):
#     """
#     Entry point when used as a library. Calls main()
#     cmd_str = f'python {matching_prog} --ncpus {num_cpus} --run_dn {run_dn} --mincov {mincov} --minprop {minprop}'
#                             if exhaustive_mode:
#                                 cmd_str += ' --exhaustive'
#                             if debug_mode:
#                                 cmd_str += ' --debug'
#     """
#     stagefile="Stage3.csv"
#     logfn = "match.log"

#     pass


def compare_var_to_ref(ref_seq, var_anno, display_width=120, colour_all=False, colour_changes=True):
    """
    Build the views for the reference sequences, sequence, variable and summary
    Args:
        ref_seq (str): dictionary of reference sequences
        var_anno (str): the sequence annotation chosen by the user
        display_width (int): width of the display in characters, wrap outputs to this width
        colour_all (bool): False, whether to colour every position in the display
        colour_changes (bool): True, whether to colour the changes of variable sites only
    Returns:
        outputs (list of tuples): each tuple contains four strings:
            ref_view, link_view, var_view, summary_view (str): the views for the reference sequences,
            the link between the reference and variable sequences, 
            the variable sequence, and a summary of the variable sequence
    Notes:
        This function builds on the variant sequence generated in ngsmatch.get_variant_seq()
        It also adjusts the reference sequence so that bases remain in matching positions
    """
    if ref_seq is None or var_anno is None:
        return [('', '', '', '')]

    parts = re.split(r'(\d+)', var_anno)
    rev_parts = parts[::-1]  # work from the end to the start
    var_seq = ref_seq  # variant sequence to be built up
    mod_seq = ref_seq  # modified sequence to be built up
    link_view = []  # link between reference and variable sequences
    summary_view = []  # summary of the variable sequence
    # need to go through changes in reverse order to avoid length changes from affecting position
    for i, p in enumerate(rev_parts):
        if i % 2 == 0:
            if p == '':
                break  # we've reached the end
            change = p
            if len(parts) <= i+1:
                print(f'Error: no matching change for position: {pos} in {parts} from {var_name}', flush=True)
                return ''
            try:
                pos = int(rev_parts[i+1]) -1  # 1-based
            except Exception as e:
                print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
                return ''
            #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
            if '+' in change:
                var_seq = var_seq[:pos] + change[1:] + var_seq[pos:]
                mod_seq = mod_seq[:pos] + '-'*len(change[1:]) + mod_seq[pos:]
            if '-' in change:
                var_seq = var_seq[:pos] + '-'*len(change[1:]) + var_seq[pos+len(change)-1:]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                var_seq = var_seq[:pos] + repl_bases +var_seq[pos+len(repl_bases):]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
    link_view = ['|' if c1 == c2 else ' ' for c1, c2 in zip(mod_seq, var_seq)]
    summary_view = []
    for m,v in zip(mod_seq, var_seq):
        if m == v:
            summary_view.append('*')
        elif m == '-':
            summary_view.append('+')
        elif v == '-':
            summary_view.append('-')
        elif v != m:
            summary_view.append(v)
        else:
            summary_view.append(' ')
    summary_view = ''.join(summary_view)
    outputs = []
    link_chrs = ''.join(link_view)
    summary_chrs = ''.join(summary_view)
    for i in range(0, len(var_seq), display_width):
        outputs.append((mod_seq[i:i+display_width],
                link_chrs[i:i+display_width],
                var_seq[i:i+display_width],
                summary_chrs[i:i+display_width]))
    return outputs


def reconstruct_sequence(var_anno, ref_seq):
    """
    Recreate the full variant sequence based on the variant annotations and the reference sequence
    args:
        var_anno (str): variant annotations
        ref_seq (str): the reference sequence variants are based on
    returns:
        variant_sequence (str): the aligned variant sequence
        mod_ref_seq (str): the matching aligned
    """
    if ref_seq is None or var_anno is None:
        return ''

    parts = re.split(r'(\d+)', var_anno)
    rev_parts = parts[::-1]  # work from the end to the start
    var_seq = ref_seq  # variant sequence to be built up
    mod_seq = ref_seq  # modified sequence to be built up
    # need to go through changes in reverse order to avoid length changes from affecting position
    for i, p in enumerate(rev_parts):
        print(f'reconstruct_sequence() {i=} {p=}', flush=True)
        if i % 2 == 0:
            if p == '':
                break  # we've reached the end
            change = p
            if len(parts) <= i+1:
                 print(f'Error: no matching change for position: {pos} in {parts}', flush=True)
                 return ''
            try:
                pos = int(rev_parts[i+1]) -1  # 1-based
            except Exception as e:
                print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
                return ''
            #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
            if '+' in change:
                var_seq = var_seq[:pos] + change[1:] + var_seq[pos:]
                mod_seq = mod_seq[:pos] + '-'*len(change[1:]) + mod_seq[pos:]
            if '-' in change:
                var_seq = var_seq[:pos] + '-'*len(change[1:]) + var_seq[pos+len(change)-1:]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                var_seq = var_seq[:pos] + repl_bases +var_seq[pos+len(repl_bases):]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
    return var_seq


def build_aligned_pair(var_anno, ref_seq):
    """
    Recreate the full variant sequence based on the variant annotations and the reference sequence
    args:
        var_anno (str): variant annotations
        ref_seq (str): the reference sequence variants are based on
    returns:
        var_seq (str): the aligned variant sequence
        mod_seq (str): the matching aligned reference sequence
    """
    if ref_seq is None or var_anno is None:
        return '',''

    parts = re.split(r'(\d+)', var_anno)
    rev_parts = parts[::-1]  # work from the end to the start
    var_seq = ref_seq  # variant sequence to be built up
    mod_seq = ref_seq  # modified sequence to be built up
    # need to go through changes in reverse order to avoid length changes from affecting position
    for i, p in enumerate(rev_parts):
        if i % 2 == 0:
            if p == '':
                break  # we've reached the end
            change = p
            if len(parts) <= i+1:
                print(f'Error: no matching change for position: {pos} in {parts} from {var_name}', flush=True)
                return ''
            try:
                pos = int(rev_parts[i+1]) -1  # 1-based
            except Exception as e:
                print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
                return ''
            #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
            if '+' in change:
                var_seq = var_seq[:pos] + change[1:] + var_seq[pos:]
                mod_seq = mod_seq[:pos] + '-'*len(change[1:]) + mod_seq[pos:]
            if '-' in change:
                var_seq = var_seq[:pos] + '-'*len(change[1:]) + var_seq[pos+len(change)-1:]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                var_seq = var_seq[:pos] + repl_bases +var_seq[pos+len(repl_bases):]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
    return var_seq, mod_seq



def run_msa(ref_id_seq, var_list):
    """
    Reconstruct full variant sequences and then perform multiple sequence alignment
    args:
        ref_id_seq (tuple(str,str)): tuple of reference id and reference seq
        var_list (list[str]): a list of variant annotations
    returns:
        aligned (cogent3 alignment object)
    """
    # need to reconstruct the full sequences from the annotations
    reconstructed_seqs = {ref_id_seq[0]:ref_id_seq[1].replace('-','')}  # start with the reference sequence
    for var_anno in var_list:
        reconstructed_seqs[var_anno] = reconstruct_sequence(var_anno, ref_id_seq[1]).replace('-','')  # remove gaps for alignment
    no_tree = make_tree(tip_names=list(reconstructed_seqs.keys()), underscore_unmunge=True)
    unaligned_seqs = make_unaligned_seqs(reconstructed_seqs, moltype='dna')
    with warnings.catch_warnings(action="ignore"):
        aln, tree = tree_align("HKY85", unaligned_seqs, tree=no_tree, show_progress=False)
    return aln


class PrepManager(object):
    """
    Producer consumer manager for preparing NGS data
    Try using multithreading or multiprocessing to run the preparation tasks
    """

    def __init__(self, run_dn, outfn, ncpus, debug):
        self.run_dn = run_dn
        self.outfn = outfn
        self.ncpus = ncpus
        self.debug = debug
        self.prep_tasks = queue.Queue()
        self.prep_results = queue.Queue()
        #self.lock = threading.Lock()


    def populate_prep_tasks(self, pids:list):
        """
        Populate the prep tasks queue with the necessary tasks for the matching process
        pids: list of plateIDs to use for populating the prep tasks
        """
        # Here you would add the logic to populate the prep_tasks queue
        # For example, you might add tasks related to preparing data for matching
        for pid in pids:
            
            # Assuming pid is a sample or plate ID, you would add a task for it
            sample = {'pid': pid}


    def add_prep_task(self, sample):
        """
        Add a preparation task to the queue
        """
        # Here you would implement the logic to add a task to the prep_tasks queue
        # For example, you might add a task related to preparing a sample for matching
        self.prep_tasks.put(sample)


    def remove_prep_task(self):
        """
        Remove a preparation task from the queue
        """
        # Here you would implement the logic to remove a task from the prep_tasks queue
        # For example, you might remove a task that has been completed
        if not self.prep_tasks.empty():
            return self.prep_tasks.get()
        return None


    def get_prep_tasks(self):
        """
        Get the list of preparation tasks
        """
        # Here you would implement the logic to return the list of preparation tasks
        # For example, you might return a list of all tasks in the prep_tasks queue
        return list(self.prep_tasks.queue)


    def run_prep_tasks(self):
        """
        Run the preparation tasks using multithreading or multiprocessing
        """
        # Here you would implement the logic to run the preparation tasks
        # For example, you might use a ThreadPool or ProcessPool to execute the tasks in parallel
        pass


    def report_prep_tasks(self):
        """
        Report the progress of the preparation tasks
        """
        # Here you would implement the logic to report the progress of the preparation tasks
        # For example, you might print the number of tasks completed or remaining
        return len(self.prep_tasks, self.prep_results)


class MatchManager(object):
    """
    Producer consumer manager for matching NGS data
    Try using multithreading or multiprocessing to run the matching tasks
    """

    def __init__(self, run_dn, targets, primer_assayfam, outfn, variants, ncpus, mincov, minprop, exhaustive, debug):
        self.run_dn = run_dn
        self.targets = targets
        self.primer_assayfam = primer_assayfam
        self.outfn = outfn
        self.variants = variants
        self.ncpus = ncpus
        self.mincov = mincov
        self.minprop = minprop
        self.exhaustive = exhaustive
        self.debug = debug
        self.prep_tasks = queue.Queue()
        self.match_tasks = queue.Queue()


    def populate_match_tasks(self):
        """
        Populate the match tasks queue with the necessary tasks for the matching process
        """
        # Here you would add the logic to populate the match_tasks queue
        # For example, you might add tasks related to performing the actual matching
        pass

    
    def run_match_tasks(self):
        """
        Run the matching tasks using multithreading or multiprocessing
        """
        # Here you would implement the logic to run the matching tasks
        # For example, you might use a ThreadPool or ProcessPool to execute the tasks in parallel
        pass


    def report_match_tasks(self):
        """
        Report the progress of the matching tasks
        """
        # Here you would implement the logic to report the progress of the matching tasks
        # For example, you might print the number of tasks completed or remaining
        pass

 
# producer task
def producer(queue):
    print('Producer: Running')
    # generate items
    for i in range(10):
        # generate a value
        value = random()
        # block, to simulate effort
        sleep(value)
        # create a tuple
        item = (i, value)
        # add to the queue
        queue.put(item)
        # report progress
        print(f'>producer added {item}')
    # signal that there are no further items
    queue.put(None)
    print('Producer: Done')
 
# consumer task
def consumer(queue):
    print('Consumer: Running')
    # consume items
    while True:
        # get a unit of work
        item = queue.get()
        # check for stop
        if item is None:
            break
        # block, to simulate effort
        sleep(item[1])
        # report
        print(f'>consumer got {item}')
    # all done
    print('Consumer: Done')
        

    def get_preprocess_progress(self):
        return 0
    
    def get_match_progress(self):
        return 0

   
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="NGS Reporting Program")
    parser.add_argument("-d", "--debug", action="store_true", help="more reporting/output for debugging purposes")
    parser.add_argument("-t", "--targets", default="targets.fa", help="file of targets in FASTA format (default=targets.fa)")
    parser.add_argument('-P', '--primer_assayfam', default="primers.csv", help='file of primer to assay family mappings')
    parser.add_argument('-o','--outfn', default='results.csv', help='Name of output file (CSV format)')
    parser.add_argument('-v','--variants', default='variant_seqs.fa', help='Name of variant sequences file (FASTA format)'),
    parser.add_argument('-n','--ncpus', type=int, default=os.cpu_count()-1, help='Number of processes to run simultaneously, default=number of CPUs in system - 1')
    parser.add_argument('-l','--logfn', default='match.log', help='Name of logging file (default=match.log)')
    parser.add_argument('-r','--run_dn', required=True, help='Path to experiment folder')
    parser.add_argument('-i','--identity',type=float,default=0.9,help='Proportional score for inexact matching')
    parser.add_argument('-m','--mincov', type=int, default=5, help='Do not match unique sequences with less than this many reads coverage, default 50')
    parser.add_argument('-p','--minprop', type=float, default=0.1, help='Do not match unique sequences '+\
            'with less than this proportion of the total number of exact matched on-target reads, default 0.2. Must be between 0.0 and 1.0')
    parser.add_argument('--inexact', action="store_true", help="enable inexact matching")
    parser.add_argument('-x','--exhaustive',action='store_true',help='Try to match every sequence, '+\
            'no matter how few counts. Ignores --minseqs and --minprop')
    parser.add_argument('-C','--no_miss_cache', action="store_true", help="disable miss cache")
    parser.add_argument('-M','--margin',type=float,default=0.9,help="Sequences must be this proportion of the reference seq length")
    parser.add_argument('-s','--stagefile', default="Stage3.csv", help="Name of the NGS genotyping Stage 3 file (default=Stage3.csv)")
    args = parser.parse_args()
    in_error = False
    print(f"{args=}", file=sys.stderr)
    lock_path = os.path.join(args.run_dn,"ngsgeno_lock") 
    if not os.path.exists(lock_path):
        try:
            with open(lock_path,"wt"):
                report_progress(args.run_dn, 0, 0)  # set this up asap
                main(args)
            print('Completed regular execution')
        except Exception as exc:
            print(f'Completed with exception {exc}', flush=True)
        if os.path.exists(lock_path):
            os.remove(lock_path)
    else:
        print("Analysis already running", file=sys.stderr)
        exit(2)
    
    
    
