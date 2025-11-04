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
import aiofiles
import asyncio
from timeit import default_timer as timer
#from ast import Bytes, Str
import os
from pathlib import Path
#from pickletools import bytes1
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
from math import ceil, floor


import Bio.SeqIO
import Bio.Align as Align

# guard/unguard duplicated from util

class UnguardedBarcodeError(Exception):
    """ Exception raised for barcodes with degraded or missing guards

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class EmptyBarcodeError(Exception):
    """ Exception raised when a barcode is an empty string

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message


class ExistingGuardError(Exception):
    """ Exception raised when a barcode appears to already have a valid set of guards

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message



# Confirm guards
def is_guarded_mbc(mbc):
    """ Return True if Musterer barcode is correctly guarded, else False """
    mbc = str(mbc)
    if mbc.startswith('m') and mbc.endswith('m'):
        return True
    return False

def is_guarded_rbc(rbc):
    """ Return True if Rodentity barcode is correctly guarded (rNNNNNr), else False """
    rbc = str(rbc)
    if rbc.startswith('r') and rbc.endswith('r'):
        return True
    return False

def is_guarded_cbc(cbc):
    """ Return True if custom barcode is correctly guarded, else False """
    cbc = str(cbc)
    if cbc.startswith('c') and cbc.endswith('c'):
        return True
    return False

def is_guarded_pbc(pbc):
    """ Return True if plate barcode is correctly guarded, else False """
    pbc = str(pbc)
    if pbc.startswith('p') and pbc.endswith('p'):
        return True
    return False

def is_guarded_abc(abc):
    """ Return True if amplicon sample barcode is correctly guarded, else False """
    abc = str(abc)
    if abc.startswith('a') and abc.endswith('a'):
        return True
    return False

def is_guarded(bc):
    """ generic check of guard status """
    if is_guarded_mbc(bc):
        return True
    elif is_guarded_rbc(bc):
        return True
    elif is_guarded_cbc(bc):
        return True
    elif is_guarded_pbc(bc):
        return True
    elif is_guarded_abc(bc):
        return True
    return False

def get_guard_type(bc):
    """ return the character of the guard in a barcode """
    if is_guarded(bc):
        return bc[0]
    return None


# Add guards to barcodes
def guard_mbc(mbc, silent=False):
    """ Add guards to Musterer barcode """
    mbc = str(mbc).strip()
    if mbc == '':
        msg = f"Musterer barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(mbc):
        if silent:
            return mbc
        msg = f"Musterer barcode appears to already have Musterer guards: {mbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(mbc):
        msg = f"Musterer barcode appears to already have Rodentity guards: {mbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(mbc):
        msg = f"Musterer barcode appears to already have custom guards: {mbc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(mbc):
        msg = f"Musterer barcode appears to already have plate guards: {mbc}"
        raise ExistingGuardError(msg)
    return 'm' + mbc + 'm'

def guard_rbc(rbc, silent=False):
    """ Add guards to Rodentity barcode """
    rbc = str(rbc).strip()
    if rbc == '':
        msg = f"Rodentity barcode is zero length string!" 
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(rbc):
        msg = f"Rodentity barcode appears to already have Musterer guards: {rbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(rbc):
        if silent:
            return rbc
        msg = f"Rodentity barcode appears to already have Rodentity guards: {rbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(rbc):
        msg = f"Rodentity barcode appears to already have custom guards: {rbc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(rbc):
        msg = f"Rodentity barcode appears to already have plate guards: {rbc}"
        raise ExistingGuardError(msg)
    return 'r' + rbc + 'r'  # Rodentity barcodes already start with 'M'

def guard_cbc(cbc, silent=False):
    """ Add guards to custom barcode """
    cbc = str(cbc).strip()
    if cbc == '':
        msg = f"Custom barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(cbc):
        msg = f"Custom barcode appears to already have Musterer guards: {cbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(cbc):
        msg = f"Custom barcode appears to already have Rodentity guards: {cbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(cbc):
        if silent:
            return cbc
        msg = f"Custom barcode appears to already have custom guards: {cbc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(cbc):
        msg = f"Custom barcode appears to already have plate guards: {cbc}"
        raise ExistingGuardError(msg)
    return 'c' + cbc + 'c'

def guard_pbc(pbc, silent=False):
    """ Add guards to plate barcode, silent ignores already correctly guarded case """
    pbc = str(pbc).strip()
    if pbc == '':
        msg = f"Plate barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(pbc):
        msg = f"Plate barcode appears to already have Musterer guards: {pbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(pbc):
        msg = f"Plate barcode appears to already have Rodentity guards: {pbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(pbc):
        msg = f"Plate barcode appears to already have custom guards: {pbc}"
        raise ExistingGuardError(msg)
    if is_guarded_abc(pbc):
        msg = f"Plate barcode appears to already have amplicon guards: {pbc}"
    if is_guarded_pbc(pbc):
        if silent:
            return pbc
        msg = f"Plate barcode appears to already have plate guards: {pbc}"
        raise ExistingGuardError(msg)
    return 'p' + pbc + 'p'

def guard_abc(abc, silent=False):
    """ Add guards to amplicon sample barcodes, silent ignores already correctly guarded case """
    abc = str(abc).strip()
    if abc == '':
        msg = f"Amplicon sample barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(abc):
        msg = f"Amplicon sample barcode appears to already have Musterer guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(abc):
        msg = f"Amplicon sample barcode appears to already have Rodentity guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(abc):
        msg = f"Amplicon sample barcode appears to already have custom guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(abc):
        msg = f"Amplicon sample barcode appears to already have plate guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(abc):
        if silent:
            return abc
        msg = f"Amplicon sample barcode appears to already have amplicon guards: {abc}"
        raise ExistingGuardError(msg)
    return 'a' + abc + 'a'


# Unguard barcodes
def unguard_mbc(mbc, silent=False):
    """ remove guards from a Musterer barcode """
    if type(mbc) != str:
        msg = f"Musterer barcode is not a string! {mbc} is type {type(mbc)}"
        raise AttributeError(msg)
    if not mbc.startswith('m') and not mbc.endswith('m') and silent:  # just return unguarded barcodes as themselves
        return mbc
    if not mbc.startswith('m') or not mbc.endswith('m'):
        msg = f"Musterer barcode guards degraded or missing in: {mbc}"
        raise UnguardedBarcodeError(msg)
    return mbc[1:-1]

def unguard_rbc(rbc, silent=False):
    """ remove guards from a Rodentity barcode """
    if type(rbc) != str:
        msg = f"Rodentity barcode is not a string! {rbc} is type {type(rbc)}"
        raise AttributeError(msg)
    if not rbc.startswith('r') and not rbc.endswith('r') and silent:  # just return unguarded barcodes as themselves
        return rbc
    if not rbc.startswith('r') or not rbc.endswith('r'):
        msg = f"Musterer barcode guards degraded or missing in: {rbc}"
        raise UnguardedBarcodeError(msg)
    return rbc[1:-1]  # Rodentity barcodes keep their 'M' prefix

def unguard_cbc(cbc, silent=False):
    """ remove guards from a custom barcode """
    if type(cbc) != str:
        msg = f"Custom barcode is not a string! {cbc} is type {type(cbc)}"
        raise AttributeError(msg)
    if not cbc.startswith('c') and not cbc.endswith('c') and silent:  # just return unguarded barcodes as themselves
        return cbc
    if not cbc.startswith('c') or not cbc.endswith('c'):
        msg = f"Custom barcode guards degraded or missing in: {cbc}"
        raise UnguardedBarcodeError(msg)
    return cbc[1:-1]

def unguard_pbc(pbc, silent=False):
    """ remove guards from a plate barcode """
    if type(pbc) != str:
        msg = f"Plate barcode is not a string! {pbc} is type {type(pbc)}"
        raise AttributeError(msg)
    if not pbc.startswith('p') and not pbc.endswith('p') and silent:  # just return unguarded barcodes as themselves
        #print(f'{pbc=}', file=sys.stderr)
        return pbc
    if not pbc.startswith('p') or not pbc.endswith('p'):
        msg = f"Plate barcode guards degraded or missing in: {pbc}"
        raise UnguardedBarcodeError(msg)
    return pbc[1:-1]

def unguard_abc(abc, silent=False):
    """ remove guards from an amplicon sample barcode """
    if type(abc) != str:
        msg = f"Amplicon sample barcode is not a string! {abc} is type {type(abc)}"
        raise AttributeError(msg)
    if not abc.startswith('a') and not abc.endswith('a') and silent:  # just return unguarded barcodes as themselves
        return abc
    if not abc.startswith('a') or not abc.endswith('a'):
        msg = f"Amplicon sample barcode guards degraded or missing in: {abc}"
        raise UnguardedBarcodeError(msg)
    return abc[1:-1]

def unguard(bc, silent=False):
    """ shortcut, for when we don't know what kind of guarding we have but need to remove it for output """
    if is_guarded_mbc(bc):
        return unguard_mbc(bc)
    elif is_guarded_rbc(bc):
        return unguard_rbc(bc)
    elif is_guarded_cbc(bc):
        return unguard_cbc(bc)
    elif is_guarded_pbc(bc):
        return unguard_pbc(bc)
    elif is_guarded_abc(bc):
        return unguard_abc(bc)
    else:
        if not silent:
            print('Possibly unguarded already:',bc, file=sys.stdout)
        return bc



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


def exact_match(seq, refseqs, rundir=None, debug=None, lock_d=None, margin=0.9):
    """
    Exact matching a merged read with a set of reference target sequences
    
    args:
    seq: the read sequence to match
    refseqs: list of reference target sequences
    rundir, debug, lock_d: used for logging debug messages (optional)
    margin: seq and refseq must be within this factor of each other in length

    returns:
    True/False, matching reference sequence/None

    If the length difference is greater than margin, don't try to match 
    If there are parentheses in the sequence, these are masked regions
    and we should split the sequence on these and match the unmasked regions.

    Simple algorithm for a single masked region and known amplicons:
    1. Cut out the two reference flanks at brackets and compare these to the sequence
    2. If they match, then we have a match

    Two bracket algorithm for a pair of masked regions:
    1. Cut out the two reference flanks at parentheses
    2. remove any square brackets, and compare these to the sequence
    3. If they match, then we have a match
    4. Assuming it's not a match repeat from 1, but finding the flanks of the square brackets
    5. Remove any parentheses and compare these to the sequence
    5. If they match, then we have a match

    General algorithm for masked regions (more than 1):
    1. split the reference sequence on masked regions
    2. for each unmasked reference region, match this to the sequence
    3. if all unmasked regions match, then we have a match
    """
    if rundir and debug and lock_d:
        msg = f"debug: exact_match() {seq=} {refseqs=}"
        asyncio.run(wdb(msg, rundir, debug, lock_d))
    for refseq in refseqs:
        if rundir and debug and lock_d:
            msg = f'debug: exact_match() Trying to match {seq=} to {refseq=}'
            asyncio.run(wdb(msg, rundir, debug, lock_d))
        
        found_match = False
        # compare length of sequenced region to reference seq
        rs = refseq.replace('[','(').replace(']',')')
        if '(' not in rs:
            msg = f"debug: exact_match() unmasked {rs=}"
            if len(seq) > len(rs):
                if len(seq)*margin > len(rs):
                    continue  # skip matching when length is too mismatched
            else:
                if len(rs)*margin > len(seq):
                    continue  # skip matching when length is too mismatched
        # unmasked case
        if rs.count('(') == 0 and rs.count(')') == 0:
            if len(seq) < len(rs):
                if seq in rs:
                    if rundir and debug and lock_d:
                        msg = f"debug: exact_match() no masking {seq=} found in {rs=}"
                        asyncio.run(wdb(msg, rundir, debug, lock_d))
                    return True, rs
            else:
                if rs in seq:
                    if rundir and debug and lock_d:
                        msg = f"debug: exact_match() no masking {rs=} found in {seq=}"
                        asyncio.run(wdb(msg, rundir, debug, lock_d))
                    return True, rs
        # two bracket type - DEPRECATED
        elif rs.count('(') == 1 and rs.count('[') == 1 and \
                rs.count(')') == 1 and rs.count(']') == 1:
            if rundir and debug and lock_d:
                msg = f"debug: exact_match() two types of masking (deprecated) {rs=}"
                asyncio.run(wdb(msg, rundir, debug, lock_d))
            left_flank = rs.split('(')[0].replace('[','').replace('[','')
            right_flank = rs.split(')')[1].replace('[','').replace('[','')
            if left_flank in seq:
                centre_pos = seq.find(left_flank) + len(left_flank)
                if right_flank in seq[centre_pos:]:
                    return True, rs
            else:
                left_flank = rs.split('[')[0].replace('(','').replace(')','')
                right_flank = rs.split(']')[1].replace('(','').replace(')','')
                if left_flank in seq:
                    centre_pos = seq.find(left_flank) + len(left_flank)
                    if right_flank in seq[centre_pos:]:
                        return True, rs
        # simple masked case
        elif rs.count('(') == 1 and rs.count(')') == 1:
            if rundir and debug and lock_d:
                msg = f"debug: exact_match() single type of masking (simple) {rs=}"
                asyncio.run(wdb(msg, rundir, debug, lock_d))
            left_flank = rs.split('(')[0]
            right_flank = rs.split(')')[1]
            if left_flank in seq:
                centre_pos = seq.find(left_flank) + len(left_flank)
                if right_flank in seq[centre_pos:]:
                    return True, rs
        # general masked case with multiple masked regions    
        elif rs.count('(') > 1 and rs.count(')') == rs.count('('):
            if rundir and debug and lock_d:
                msg = f"debug: exact_match() multiple masked regions {rs=}"
                asyncio.run(wdb(msg, rundir, debug, lock_d))
            seq_index = 0
            failed_rs = False
            for section in rs.split('('):
                if ')' in section:
                    non_variable_seq = section.split(')')[1]
                else:
                    non_variable_seq = section
                seq_index = seq.find(non_variable_seq, seq_index)
                if seq_index == -1:
                    failed_rs = True
                    break
                seq_index += len(non_variable_seq)  # move past the non-variable sequence
            if not failed_rs:
                return True, rs
    if rundir and debug and lock_d:
        msg = f"debug: exact_match() no match found for {seq=}"
        asyncio.run(wdb(msg, rundir, debug, lock_d))
    return False, None


def test_exact_match():
    """
    Some test cases to ensure that exact matching of a sequence to a set of reference sequences is working correctly
    """
    refseqs = ['ATGCGTGTTCAAGTACACCCAAGTTGACAGTGCA']
    seq = 'ATGCGTGTTCAAGTACACCCAAGTTGACAGTGCA'
    success, rs = exact_match(seq, refseqs)
    if not success:
        print('Test exact match failed')
        return False
    
    refseqs = ['ATGCGTG(TTC)AAGTACACCCAAGTTGACAGTGCA']
    seq = 'ATGCGTGTCAAGTACACCCAAGTTGACAGTGCA'
    success, rs = exact_match(seq, refseqs)
    if not success:
        print('Test exact match with round brackets failed')
        return False
    
    refseqs = ['ATGCGTG(TTC)AAGTACAC(CC)AAGTTGAC(A)GTGCA']
    seq = 'ATGCGTGAAGTACACAAGTTGACAGTGCA'
    success, rs = exact_match(seq, refseqs)
    if not success:
        print('Test exact match with multiple round brackets failed')
        return False
    
    refseqs = ['ATGCGTG(TTC)AAGTACACCCAAGTTGAC[A]GTGCA']
    seq1 = 'ATGCGTGTTCAAGTACACCCAAGTTGACAGTGCA'
    seq2 = 'ATGCGTGAAGTACACCCAAGTTGACAGTGCA'
    seq3 = 'ATGCGTGTTCAAGTACACCCAAGTTGACGTGCA'
    success1, rs1 = exact_match(seq1, refseqs)
    success2, rs2 = exact_match(seq2, refseqs)
    success3, rs3 = exact_match(seq3, refseqs)
    if not success1:
        print(f'Test exact match with square brackets 1 failed for {seq1}')
        return False
    if not success2:
        print(f'Test exact match with square brackets 2 failed for {seq2}')
        return False
    if not success3:
        print(f'Test exact match with square brackets 3 failed for {seq3}')
        return False
    return True


def join_diffs(last_diffs):
    """ 
    Create a new diff string from a contiguous set of same-type differences
    Positions should be zero-based
    Each diff is of the form (i,ref_c,alt_c)
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
    Must be free of brackets and parentheses
    """
    lmargin = max(len(alnref)-len(alnref.lstrip('-')),len(alnalt)-len(alnalt.lstrip('-')))
    rmargin = max(len(alnref)-len(alnref.rstrip('-')),len(alnalt)-len(alnalt.rstrip('-')))
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


def inexact_match(seq, refseqs, rundir=None, debug=None, lock_d=None, identity=0.90):
    """
    Inexact matching of a merged sequence with a set of reference target sequences
    args:
    seq (str): query sequence
    refseqs (list[str]): reference sequences

    Notes:
    Declare a match if we have 90% sequence identity over the matched length or better
    Removes all brackets in reference prior to matching
    Identity must be >= 0 and <= 1.0
    We must check for the best match against all targets or we will incorrectly call variant sequences!
    Return True/False, and matching reference sequence/None
    rundir, debug, lock_d are used for logging debug messages
    """
    if rundir and debug and lock_d:
        msg = f"debug: inexact match {seq=}"
        asyncio.run(wdb(msg, rundir, debug, lock_d))
    if len(refseqs) == 0:
        if rundir and debug and lock_d:
            asyncio.run(wdb('Warning: No reference sequences to align to', rundir, debug, lock_d))
        return False, None, None

    if identity < 0:
        identity = 0.0
        if rundir and debug and lock_d:
            asyncio.run(wdb('Warning: identity < 0, setting to 0.0', rundir, debug, lock_d))
    elif identity > 1.0:
        identity = 1.0
        if rundir and debug and lock_d:
            asyncio.run(wdb('Warning: identity > 1.0, setting to 1.0', rundir, debug, lock_d))

    # clean up brackets in reference sequences
    cleanseqs_refseqs = {rs.replace('(','').replace(')','').replace('[','').replace(']',''):rs for rs in refseqs}
    cleanseqs = cleanseqs_refseqs.keys()
        
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.7
    aligner.end_gap_score = -0.5
    aligner.match_score = 3.5
    aligner.mismatch_score = -1.0
    #aligner.mode = 'local'
    #aligner.open_gap_score = -1.0
    #aligner.extend_gap_score = -1.0
    #aligner.end_gap_score = -1.0
    #aligner.match_score = 1
    #aligner.mismatch_score = -1.0
    try:
        # do alignment against a reference sequence if they aren't too different in length
        # order is align(Target, Query)
        all_alignments = [(aligner.align(rs, seq)[0], rs) for rs in refseqs]
    except Exception as exc:
        if rundir and debug and lock_d:
            msg = f'critical: inexact_match() alignment failed with exception {exc}'
            print(msg, flush=True)
            asyncio.run(wdb(msg, rundir, debug, lock_d))
        return False, None, None

    if len(all_alignments) == 0:
        if rundir and debug and lock_d:
            asyncio.run(wdb(f'critical: inexact_match() No alignments found {seq=}', rundir, debug, lock_d))
        return False, None, None
    best_score = 0
    best_algn = None
    best_rs = None
    for algn, rs in all_alignments:
        if algn.score > best_score:
            best_score = algn.score
            best_algn = algn
            best_rs = rs
    if best_score < 1:
        if debug:
            asyncio.run(wdb(f'No alignments found {seq=}', rundir, debug, lock_d))
        return False, None, None
    
    target_sections = []
    query_sections = []
    for segment in str(best_algn).split('\n'):
        if segment.startswith('target'):
            tokens = [s for s in segment.split(' ') if s != '']
            if len(tokens) > 2:
                target_sections.append(tokens[2])
        elif segment.startswith('query'):
            tokens = [s for s in segment.split(' ') if s != '']
            if len(tokens) > 2:
                query_sections.append(tokens[2])
    aligned_rs = ''.join(target_sections)
    aligned_seq = ''.join(query_sections)
    adjusted_query_length = len(aligned_seq.lstrip('-').rstrip('-'))
    adjusted_target_length = len(aligned_rs.lstrip('-').rstrip('-'))
    match_length = min(adjusted_query_length, adjusted_target_length)

    identities = best_algn.counts().identities
    cutoff = match_length * identity
    if identities > cutoff:
        # Annotate match
        seq_anno = annotate_seq(aligned_rs, aligned_seq)
        best_ref = cleanseqs_refseqs[best_rs]
        return True, best_ref, seq_anno
    msg = f'Best alignment identity worse than {identity*100}% {aligned_rs=} {aligned_seq=} {best_rs=} {seq=}'
    if debug:
        asyncio.run(wdb(msg, rundir, debug, lock_d))
    return False, None, None


async def preprocess_seqs(wr, rundir, log, lock_l, lock_d, debug=False):
    """
    Call bbduk for cleaning of reads
    Call bbmerge to merge read pairs
    If number of merged reads is less than half the number of unmerged reads then call this a failure
    Return seqcnt (Counter of unique sequences), passing (T/F), readCount, cleanCount, MergeCount, msg
    """
    readCount, cleanCount, cleanCount2, joinCount, mergeCount = -1, -1, -1, -1, -1
    rfn= "*{}-{}_*_R1_*.fastq.gz".format(unguard(wr['pcrPlate'],silent=True), padwell(wr['pcrWell']))
    fn1s = glob.glob(os.path.join(rundir, "raw", rfn))
    if not fn1s:
        # we need to look for newer Illumina files with an underscore
        rfn= "*{}_{}_*_R1_*.fastq.gz".format(unguard(wr['pcrPlate'],silent=True), padwell(wr['pcrWell']))
        fn1s = glob.glob(os.path.join(rundir, "raw", rfn))
    fn1s = glob.glob(os.path.join(rundir, "raw", rfn))
    lrecs = []
    if not fn1s:
        with lock_l:
            log.append(f"no data for {fn1s}")
        await wdb(f"no data for {fn1s}", rundir, debug, lock_d)
        return {}, False, readCount, cleanCount, mergeCount, "No files"
            
    if len(fn1s)>1:
        with lock_l:
            log.append(f"too many reads files for pattern{rfn}")
            log.append(f"   R1 files= {fn1s}")
        await wdb(f"too many reads files for pattern{rfn}\n    R1 files= {fn1s}", rundir, debug, lock_d)
        return {}, False, readCount, cleanCount, mergeCount, "Too many files"
            
    fnr1 = fn1s[0]
    fn1 = os.path.basename(fnr1) # just the one file
    fn2 = fn1.replace('_R1_001', '_R2_001')
    fnr2 = os.path.join(rundir, 'raw', fn2) 
            
    # find the data file
    if not os.path.isfile(fnr2):
        with lock_l:
            log.append(f"missing file: {fnr2}")
        await wdb(f"missing file: {fnr2}", rundir, debug, lock_d)
        return dict(), False, readCount, cleanCount, mergeCount, "R2 file missing"
        
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
        try:
            pres1 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
        except Exception as exc:
            msg = f"failed to run bbmerge {exc}".replace("'","").replace('"',"'").replace('\\','/')
            return None, False, readCount, cleanCount, mergeCount, msg
        if pres1.returncode != 0:
            return None, False, readCount, cleanCount, mergeCount, "Failed: bbduk run failed"
        # run BBMerge to join paired-end reads
        # cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} verystrict=t".format(*((bbmapd,)+fncs+fnms)).split()
        cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} overwrite=true pfilter=1".format(*((bbmapd,)+fncs+fnms)).split()
        try:
            pres2 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
        except Exception as exc:
            msg = f"failed to run bbmerge {exc}".replace("'","").replace('"',"'").replace('\\','/')
            return None, False, readCount, cleanCount, mergeCount, msg
        if pres2.returncode != 0:
            return None, False, readCount, cleanCount, mergeCount, "Failed: bbmerge run failed"
        # could delete cleaned data once it's been merged.
                
        # keep the log output as record counts get used
        async with aiofiles.open(fnlog,'wt') as log_file:
            for s in (pres1.stdout, pres1.stderr):
                if s:
                    await log_file.write(s)
            await log_file.write('======\n')
            for s in (pres2.stdout, pres2.stderr):
                if s:
                    await log_file.write(s)
        log1, log2 = ['\n'.join((p.stdout, p.stderr)) for p in (pres1, pres2)]
    else:
        fnm_fmt = fnmfmt.replace('{}.{}','')
        lrecs.append(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}")
        if debug:
            await wdb(f"Debug: Skipping cleaning and merging of reads for {fnm_fmt}", rundir, debug, lock_d)
        async with aiofiles.open(fnlog, 'rt') as f:
            log_data = await f.read()
        log1, log2 = log_data.split("======\n", 1)
            

    await wdb(f'Got here {log1=} {log2=}', rundir, debug, lock_d)
    # get counts from BBDuk & BBmerge output
    
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
        if debug:
            await wdb(msg, rundir, debug, lock_d)
        return {}, False, readCount, cleanCount, joinCount, msg
                
    
    fn = fnms[0] # merged reads file
    seq_strings = await read_FASTQ_async(fn, rundir, debug, lock_d, seqs_only=True)
    seqcnt = collections.Counter(seq_strings)
    
    mergeCount = sum(seqcnt.values())
    if mergeCount != joinCount:
        with lock_l:
            msg = f"Merged counts mismatch. {mergeCount} != {joinCount}"
            log.append(msg)
        if debug:
            await wdb(msg, rundir, debug, lock_d)
        return {}, False, readCount, cleanCount, mergeCount, msg

    with lock_l:
        for l in lrecs:
            log.append(l)

    return seqcnt, True, readCount, cleanCount, mergeCount, None


async def wdb(msg, rundir, debug, lock_d=None):
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
                f = None
                try:
                    f = await aiofiles.open(debugfn, 'at')
                    await f.write(f'{msg}\n')
                except:
                    print('Could not write to debug')
                    print(f'{msg}', flush=True)
                finally:
                    if f:
                        await f.close()                  
        else:
            f = None
            try:
                f = await aiofiles.open(debugfn, 'at')
                await f.write(f'{msg}\n')
            except:
                print('Could not write to debug', flush=True)
                print(f'{msg}')
            finally:
                if f:
                    await f.close()
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
        match_cache, anno_cache, reps, log, lock_mtc, lock_r, 
        lock_l, lock_d, margin, identity, mincov, minprop, inexact,
        exhaustive, debug=False):
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
    reps: dict of sampleNo to all output fields as a list {sampleNo:list} 
    lock_mtc: lock object for the match cache
    lock_r: lock object for the results    
    lock_l: lock for logging
    lock_d: lock for writing to debug
    margin: minimum proportion overlap for exact matching
    identity: proportional identity score for inexact matching
    inexact: enable inexact matching (genotyping only)
    exhaustive: [T/F] don't skip any sequences, not matter how low their count
    debug: write messages direct to file (slow I/O)

    Parallel, independent processing of a single well.
    Clean and merge reads using bbduk
    Group sequences into archetypes and order by decreasing count
    Check whether any sequences exactly match the expected assay family and base cutoffs on this
    Check all remaining sequences until counts falls below cutoff
      -> Exact matches against prefered assays
      -> Exact matches against all assay families
      The miss cache is problematic... we may miss amplicons that are out of context (removed)
      -> Inexact matches against the match cache -> add variants to match cache

    Match cache ["ATC..."]="ATC..." <- a sequence in seq_ids and id_seq
    Writes back result (seqCount,seqName,Efficiency,otherCount,otherName) to result thread
    """
    PID = os.getpid()
    msg = f"Info: {PID} Starting work block {work_block} for {wr}"
    if debug:
        asyncio.run(wdb(msg, rundir, debug, lock_d))
    lrecs = [] # log records, lock and write on exit/return
    amplicon_run = False
    if not primer_assayfam or not assayfam_primers:
        amplicon_run = True
    sampleNo = wr.get('sampleNo', None)
    pcrPlate = wr.get('pcrPlate', None)
    pcrWell = wr.get('pcrWell', None)
    if sampleNo is None or pcrPlate is None or pcrWell is None:
        msg = f"Critical: {sampleNo=}, {pcrPlate=} or {pcrWell=} missing from Stage3.csv!"
        if debug:
            asyncio.run(wdb(msg, rundir, debug, lock_d))
        retval = [str(wr[x]) for x in wr] + [-1,-1,-1, '', '', -1, msg]
        with lock_r:
            reps[sampleNo] = retval
        with lock_l:
            for l in lrecs:
                log.append(l)
        return
    pcrPlate = unguard(pcrPlate, silent=True)
    pcrWell = padwell(pcrWell)
    msg = f"Info: {PID} Working on {sampleNo=} {pcrPlate=} {pcrWell=}"
    lrecs.append(msg)
    if debug:
        asyncio.run(wdb(msg, rundir, debug, lock_d))
        print(msg, flush=True)

    # clean and merge FASTQs
    seqcnt, success, readCount, cleanCount, mergeCount, fault_msg = asyncio.run(preprocess_seqs(wr, rundir, log, 
            lock_l, lock_d, debug=debug))
    if not success:
        msg = f"Failed: preprocessing for {pcrPlate} {pcrWell} with {fault_msg}"
        lrecs.append(msg)
        if debug:
            print(msg, flush=True)
            print(f'{readCount=} {cleanCount=} {mergeCount=} {fault_msg=}', flush=True)
            asyncio.run(wdb(msg, rundir, debug, lock_d))
        
        retval = [str(wr[x]) for x in wr] + [readCount, cleanCount, mergeCount, '','',-1,'Failed to run bbduk']
        if debug:
            asyncio.run(wdb(msg, rundir, debug, lock_d))
            print(f'Failed {retval=}', flush=True)
        with lock_r:
            reps[int(sampleNo)] = retval
        with lock_l:
            for l in lrecs:
                log.append(l)
        return

    msg = f"Info: Completed preprocessing for {pcrPlate} {pcrWell} "+\
            f"{readCount=} {cleanCount=} {mergeCount=}"
    lrecs.append(msg)
    if debug:
        asyncio.run(wdb(msg, rundir, debug, lock_d))
                
    # if we get an unknown primer then we should run archetypes and call it a day?
    primer = wr.get('primer', 'No_primer')
    if not amplicon_run and primer not in primer_assayfam:
        msg = f'Warning: primer {primer} unknown'
        lrecs.append(msg)
        if debug:
            asyncio.run(wdb(msg, rundir, debug, lock_d))
        
    # decide which assays are on-target vs off-target
    on_target_ids = set()
    off_target_ids = set()
    on_target_seqs = set()
    off_target_seqs = set()
   
    if amplicon_run:
        for name in id_seq:
            if primer.lower() in name.lower():
                # no worrying about splitting names here
                on_target_ids.add(name)
                on_target_seqs.add(id_seq[name])
            else:
                off_target_ids.add(name)
        # avoid looking at the same sequence under different names - used by best_match
        off_target_seqs = set(list(seq_ids.keys())).difference(on_target_seqs)
    else:  # genotyping
        for name in id_seq:
            if name.lower().startswith(primer.split('_')[0].lower() + '_'):
                on_target_ids.add(name)
                on_target_seqs.add(id_seq[name])
            else:
                off_target_ids.add(name)
        # avoid looking at the same sequence under different names - used by best_match
        off_target_seqs = set(list(seq_ids.keys())).difference(on_target_seqs)
    if debug:
        msg = f"Debug: {on_target_ids=} {len(off_target_seqs)=}"
        asyncio.run(wdb(msg, rundir, debug, lock_d))

    match_cnt = Counter()
    if mincov < 1:
        mincov = 1
    if minprop > 1.0:
        minprop = 1.0
    elif minprop < 0.0:
        minprop = 0.0
        
    # unique sequences only (substrings collapsed), from most to least common
    if exhaustive:
        if debug:
            asyncio.run(wdb('Info: Aggregating archetype sequences', rundir, debug, lock_d))
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
    msg = f"Info: Matching {len(seqcnt)} sequences for {pcrPlate} {pcrWell} "+\
            f"with {family_exact_counts=} {low_cov_cutoff=}"
    if debug:
        asyncio.run(wdb(msg, rundir, debug, lock_d))

    for seq, num in seqcnt.most_common():
        if (num >= low_cov_cutoff) or exhaustive:
            if debug:
                msg = f"debug: Matching {pcrPlate} {pcrWell} on process {PID} with counts {num} and {seq=}"
                asyncio.run(wdb(msg, rundir, debug, lock_d))
            
            # check if we already have this on record
            if seq in mtc:
                if debug:
                    msg = f"debug: hit to match cache with {seq=} and counts {num}"
                    asyncio.run(wdb(msg, rundir, debug, lock_d))
                ref_seq = mtc[seq]
                match_cnt[ref_seq] += num
                if seq in anc:
                    seq_anno = anc[seq]
                    match_cnt[ref_seq + '//' + seq_anno] += num
                continue

            # on target exact matching
            is_match, ref_seq = exact_match(seq, on_target_seqs, rundir, debug, lock_d, margin=margin)
            if is_match:
                if debug:
                    msg = f"debug: Exact match against {seq_ids[ref_seq]} with {seq} and counts {num}"
                    asyncio.run(wdb(msg, rundir, debug, lock_d))
                match_cnt[ref_seq] += num
                continue   

            # off target exact matching for everyone!
            is_match, ref_seq = exact_match(seq, off_target_seqs, rundir, debug, lock_d, margin=margin)
            if is_match:
                if debug:
                    msg = f"debug: Exact match against {seq_ids[ref_seq]} with {seq} and counts {num}"
                    asyncio.run(wdb(msg, rundir, debug, lock_d))
                match_cnt[ref_seq] += num
                continue

            # We want to only compare inexact matches against the primer that was used - but don't add counts!
            if amplicon_run:
                is_match, ref_seq, seq_anno = inexact_match(seq, on_target_seqs, rundir, debug, 
                        lock_d, identity=identity)
                if is_match:
                    if debug:
                        msg = f"debug: Inexact match against {seq_ids[ref_seq]} {wr['pcrPlate']}"+\
                                f" {wr['pcrWell']} {num} {seq}\n"
                        asyncio.run(wdb(msg, rundir, debug, lock_d))
                    mtc[seq] = ref_seq
                    # match_cnt[ref_seq] += num  # don't add to the original count!
                    if seq_anno:
                        anc[seq] = seq_anno
                        match_cnt[ref_seq+'//'+seq_anno] += num
                    continue
            else:  # only inexact for genotypes if requested, and match against all targets        
                if inexact:
                    # inexact matching must be done against all known sequences at once or it risks false association
                    is_match, ref_seq, seq_anno = inexact_match(seq, on_target_seqs.union(off_target_seqs),rundir, debug, 
                            lock_d, identity=identity)

                    if is_match:
                        if debug:
                            msg = f"debug: Inexact match against {seq_ids[ref_seq]} {wr['pcrPlate']}"+\
                                    f" {wr['pcrWell']} {num} {seq}\n"
                            asyncio.run(wdb(msg, rundir, debug, lock_d))
                        mtc[seq] = ref_seq
                        match_cnt[ref_seq] += num
                        if seq_anno:
                            anc[seq] = seq_anno
                            match_cnt[ref_seq+'//'+seq_anno] += num
                        continue
            # Lastly, no match, map to itself
            match_cnt[seq] += num
 
        else:
            if debug:
                msg = f"Debug: too few reads {num} for matching {seq}"
                asyncio.run(wdb(msg, rundir, debug, lock_d))
            match_cnt['other'] += num
            other_count += 1
        # try another sequence

    if len(mtc) > original_match_cache_size:
        with lock_mtc:
            match_cache.update(mtc) 
            anno_cache.update(anc)
        if debug:
            asyncio.run(wdb('Info: Updating match and annotation caches', rundir, debug, lock_d))
           
    msg = f"Info: Completed matching for {PID=} {pcrPlate=} {pcrWell=}"
    lrecs.append(msg)
    if debug:
        asyncio.run(wdb(msg, rundir, debug, lock_d))
    # rename 'other' to include number of separate sequence variations
    total_others = match_cnt['other']
    other_name = 'other (' + str(other_count) +')'
    match_cnt[other_name] = total_others
    res1 = [readCount, cleanCount, mergeCount]
    # name the outputs and counts
    seqCounts = []
    seqNames = []
    otherCounts = []
    otherNames = []

    for seq, count in match_cnt.most_common():
        if seq == 'other':
            continue
        if seq.startswith('other') and seq != 'other':
            otherNames.append(seq)
            otherCounts.append(count)
            continue

        if '//' not in seq:  # includes novel sequences
            if seq in mtc:
                refseq = mtc[seq]
                if refseq in on_target_seqs:
                    seqCounts.append(count)
                    refname = [name for name in seq_ids[refseq] if name in on_target_ids][0]
                    seqNames.append(refname)
                    if seq in anc:
                        otherCounts.append(count)
                        otherNames.append(refname+'//'+anc[seq])
                else:
                    for name in seq_ids[refseq]:
                        otherCounts.append(count)
                        otherNames.append(name)
                continue
            else:  # novel sequences
                otherCounts.append(count)
                otherNames.append(seq)

        else:
            ref_seq = seq.split('//')[0]
            anno = seq.split('//')[1]
            if ref_seq not in seq_ids:
                msg = f'Critical: {ref_seq=} not found in sequence to ID mapping'
                print(msg, flush=True)
                if debug:
                    asyncio.run(wdb(msg, rundir, debug, lock_d))
            else:
                for name in seq_ids[ref_seq]:
                    otherCounts.append(count)
                    otherNames.append(name+'//'+anno)
    # calculate efficiency of PCR as the sum(seqCounts)/wr['mergeCount']
    try:
        efficiency = sum(map(int,seqCounts))/mergeCount
    except Exception as exc:
        if debug:
            asyncio.run(wdb(f"Error: calculating efficiency {seqCounts=} {mergeCount=} {exc=}", rundir, debug, lock_d))
        effiency = -1.0
    # combine outputs
    try:
        res2 = [';'.join(map(str,seqCounts)), ';'.join(map(str,seqNames)), f'{round(efficiency,3)}', 
                ';'.join(map(str,otherCounts)), ';'.join(map(str,otherNames))]
    except Exception as exc:
        if debug:
            asyncio.run(wdb(f"Error: joining names and counts {exc=}", rundir, debug, lock_d))
        res2 = ['','','','','']
    retval = [str(wr[x]) for x in wr] + res1 + res2
    sn = None
    try:
        sn = int(sampleNo)
    except Exception as exc:
        msg = f"Critical: sampleNo not an integer {exc=}"
        print(msg, flush=True)
        if debug:
            asyncio.run(wdb(msg, rundir, debug, lock_d))
    if sn:
        if debug:    
            asyncio.run(wdb(f'{retval=}', rundir, debug, lock_d))
        with lock_r:
            reps[sn] = retval
    if debug:
        asyncio.run(wdb(f'{retval=}', rundir, debug, lock_d))
    with lock_l:
        for l in lrecs:
            log.append(l)
    if debug:        
        msg = f"Info: Exiting process_well() {PID=} {sampleNo=} {pcrPlate=} {pcrWell=}"
        asyncio.run(wdb(msg, rundir, debug, lock_d))
        print(msg)


async def write_log(log, logfn):
    print(f"Writing log to {logfn}", file=sys.stderr)
    async with aiofiles.open(logfn, 'wt') as f:
        for l in log:
            await f.write(l + '\n')


def get_raw_fastq_pairs(dirpath):
    """ return a sorted list of tuple(R1_path, R2_path) to raw FASTQ files """
    valid_pairs = []
    r1s = [dirpath/Path(f) for f in os.listdir(dirpath) if f.endswith('.fastq.gz') and '_R1_' in f]
    for r1 in r1s:
        r2 = Path(str(r1).replace('_R1_001','_R2_001'))
        if r1.is_file() and r2.is_file():
            valid_pairs.append((r1,r2))
    return sorted(valid_pairs)


async def write_report_progress(rundir, launch_progress, match_progress):
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
    progress_fn = os.path.join(rundir,'match_progress_'+str(launch_progress)+'_'+str(match_progress))
    async with aiofiles.open(progress_fn, 'wt') as f:
        await f.write('')


async def read_FASTA_async_gen(filepath, errors='ignore'):
    """
    Asynchronously read a FASTA file and yield entries
    """
    async with aiofiles.open(filepath, 'rt', errors=errors) as f:
        header = None
        seq = []
        async for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield (header, ''.join(seq))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        yield (header, ''.join(seq))


async def read_FASTQ_async(filepath, rundir, debug, lock_d, seqs_only=False, errors='ignore'):
    """
    Asynchronously read a FASTQ file and yield entries
    """
    try:
        if filepath.endswith('.gz'):
            try:
                async with aiofiles.open(filepath, mode='rb') as f:
                    decompressed_data = gzip.decompress(await f.read())
                    try:
                        lines = decompressed_data.decode('utf-8', errors=errors).splitlines()
                    except Exception as exc:
                        await wdb(f'could not decode FASTQ {exc}', rundir, debug, lock_d)
                    if not seqs_only:
                        return lines
                    seqs = []
                    for i, line in enumerate(lines):
                        if i % 4 == 1:
                            seqs.append(line)
                    return seqs
            except Exception as exc:
                await wdb(f'read_FASTQ() {exc}', rundir, debug, lock_d)
                return []
        else:
            async with aiofiles.open(filepath, mode='rt', errors=errors) as f:
                lines = await f.readlines()
                if not seqs_only:
                    return lines
                seqs = []
                for i, line in enumerate(lines):
                    if i % 4 == 1:
                        seqs.append(line)
                return seqs
    except Exception as exc:
        await wdb(f'read_FASTQ failed {exc}', rundir, debug, lock_d)
        return []


async def read_file_async_gen(filepath, errors='ignore'):
    """
    Asynchronously read a file and yield lines
    """
    async with aiofiles.open(filepath, mode='rt', errors=errors) as f:
        async for line in f:
            yield line.strip()


async def parse_targets(rundir, targets, log, debug):
    """
    Return a dictionary of sequences to a set of ids, and
    a dictionary of ids to sequences
    """
    seq_ids = {}
    id_seq = {}
    async for header,seq in read_FASTA_async_gen(os.path.join(rundir,targets)):
        if not header or not seq:
            continue
        name = header[1:] if header.startswith('>') else header
        seq = seq.strip().upper()
        if seq not in seq_ids:
            seq_ids[seq] = set()
        seq_ids[seq].add(name)
        if name in id_seq:
            msg = f'Warning: {name} duplicated in reference'
            log.append(msg)
            wdb(msg, rundir, debug)
        id_seq[name] = seq
    return seq_ids, id_seq


async def parse_primer_assayfams(rundir, paf):
    """
    Two column file (tab separated) of primer and assay families
    Return a mapping of primer to set of assay family, and mapping of assay family to set of primers
    """
    primer_assayfam = {}  # [primer] = []
    assayfam_primers = {}  # [assayfam] = []
    pmr_fn = os.path.join(rundir, paf)
    try:
        async for line in read_file_async_gen(pmr_fn):
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


async def parse_primer_file(paf):
    """
    Parse the user provided assay list file. Use this if testing without a primer_assayfam file
    """
    primer_assayfam = {}
    assayfam_primers = {}  # [assayfam] = []
    try:
        async for line in read_file_async_gen(paf):
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
    except Exception as exc:
        print(f'Critical: failed to open {paf} for primer assay family info {exc=}')
    return primer_assayfam, assayfam_primers


def get_variant_seq(var_name, id_seq):
    """
    Reverse engineer the variant sequence from annotations in the var name
    Positions are 1-based
    """
    if '//' not in var_name:
        return var_name  # We should never see this happen, but to be safe...
    primer = var_name.split('//')[0]
    if primer not in id_seq:
        print(f'Error: no known primer: {primer}')
        return ''
    original_seq = id_seq[primer].replace('(','').replace('(','').replace('[','').replace(']','')
    parts = re.split(r'(\d+)', var_name.split('//')[1])
    rev_parts = parts[::-1]
    new_seq = original_seq
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
            if '+' in change:
                new_seq = new_seq[:pos] + change[1:] + new_seq[pos:]
            if '-' in change:
                new_seq = new_seq[:pos] + new_seq[pos+len(change)-1:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                new_seq = new_seq[:pos] + repl_bases +new_seq[pos+len(repl_bases):]
    return new_seq


def test_variant_seq():
    """
    unit tests for get_variant_seq()
    """
    print('beginning unit tests for variant sequence recreation', flush=True)
    test_id_seq = {'orig':'ACTGAACCTTGG'}
    test1 = 'orig//4+C'
    expected = 'ACTCGAACCTTGG'
    new_seq = get_variant_seq(test1, test_id_seq)
    if expected != new_seq:
        print(f'Test 1: {expected} does not match {new_seq}!', flush=True)
    else:
        print('Test 1: pass', flush=True)
        
    test2 = 'orig//4-G'
    expected = 'ACTAACCTTGG'
    new_seq = get_variant_seq(test2, test_id_seq)
    if expected != new_seq:
        print(f'Test 2: {expected} does not match {new_seq}!', flush=True)
    else:
        print('Test 2: pass', flush=True)
        
    test3 = 'orig//4G/C'
    expected = 'ACTCAACCTTGG'
    new_seq = get_variant_seq(test3, test_id_seq)
    if expected != new_seq:
        print(f'Test 3: {expected} does not match {new_seq}!', flush=True)
    else:
        print('Test 3: pass', flush=True)
        

# def main(rundir, stagefile, logfn, outfn, targets, exhaustive=False, debugfn='debug.log', debug=False):
async def main(args):
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
    try:
    #if True:  # helps with debugging
        log.append('Info: Run with the following command line options:')
        for arg in vars(args):
            log.append(f'Info: {arg} {getattr(args, arg)}')
        
        # read the wells data
        start_time = datetime.datetime.now()
        log.append(f"Begin: {start_time}")
        raw_pair_list = get_raw_fastq_pairs(os.path.join(args.rundir, 'raw'))
        if args.amplicons:
            amplicon_plates = set([unguard_pbc(a, silent=True) for a in args.amplicons])
            raw_pair_list = [f for f in raw_pair_list if str(f[0].name).split('_')[0] in amplicon_plates]
        raw_file_identifiers = [str(f[0].name).split('_')[0] for f in raw_pair_list]
        if len(raw_file_identifiers) == 0:
            print('No raw FASTQ files found in the raw directory', file=sys.stderr)
            return
        # New Illumina doesn't have a dash in the file name
        if raw_file_identifiers[0][-4] != '-':
            raw_file_identifiers = ['-'.join(str(f[0].name).split('_')[0:2]) for f in raw_pair_list]
        raw_file_identifiers = set(raw_file_identifiers)
        #print(sorted(raw_file_identifiers), file=sys.stderr)
        with open(os.path.join(args.rundir,args.stagefile)) as srcfd:
            src = csv.reader(srcfd, dialect="unix")
            hdr = next(src)
            WRec = collections.namedtuple("WRec", hdr)
            wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrPlate, x.pcrWell[0], int(x.pcrWell[1:]), x.primer))
        #print(wdata, file=sys.stderr)
        if args.amplicons:
            # select only those wells that have amplicon guards on the sample barcodes
            wdata = [rec for rec in wdata if is_guarded_abc(rec.sampleBarcode)]
        else:
            # select only those wells that have no amplicon guards on the sample barcodes
            wdata = [rec for rec in wdata if not is_guarded_abc(rec.sampleBarcode)]
        wdata = [rec for rec in wdata if unguard(rec.pcrPlate, silent=True) +'-'+ padwell(rec.pcrWell) in raw_file_identifiers]
        #print(f'After filtering by available files {wdata=}', file=sys.stderr)
        log.append(f"Info: {len(wdata)} sample wells to process.")
        if len(wdata) == 0:
            return
                    
        # read the target sequence file into dictionaries [seq] = set([ids]), and [id] = seq
        seq_ids, id_seq = await parse_targets(args.rundir, args.targets, log, args.debug)
        print(f'Parsed {len(seq_ids)} targets from reference', flush=True)

        # get relationship of primers to assay family for genotyping
        #print(f'{args.primer_assayfam=}')
        if not args.amplicons:
            primer_assayfam, assayfam_primers = await parse_primer_assayfams(args.rundir, args.primer_assayfam)
            print(f'Parsed {len(primer_assayfam)} primers/assays', flush=True) # {primer_assayfam=} {assayfam_primers=}')
        else:
            primer_assayfam = {}
            assayfam_primers = {}

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
            reps = manager.dict()  # [sampleNo] = all columns
            logm = manager.list()
            lock_mtc = manager.Lock()  # match_cache locking
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
            reports = []
            for i, wr in enumerate(wrs):
                print(f'Launching process for sampleNo {wr.get("sampleNo", "")} {i+1} of {total_jobs}', flush=True)
                r = pool.apply_async(process_well, args=(i, wr, args.rundir, seq_ids, id_seq, primer_assayfam, 
                    assayfam_primers, match_cache, anno_cache, reps, logm, lock_mtc, 
                    lock_r, lock_l, lock_d, args.margin, args.identity, args.mincov, args.minprop, 
                    args.inexact, args.exhaustive, args.debug))
                reports.append(r)
                if i % 3 == 0:
                    launch_progress = ceil(100*i/total_jobs)
                    completed = sum([r.ready() for r in reports])
                    match_progress = floor(100*completed/total_jobs)
                    await write_report_progress(args.rundir, launch_progress, match_progress)
                    
            while match_progress < 100:
                completed = sum([r.ready() for r in reports])
                # report the number of remaining tasks
                print(f'Match completion: {100*completed/total_jobs}')
                match_progress = 100*completed/total_jobs
                await write_report_progress(args.rundir, 100, floor(match_progress))
                # wait a moment
                time.sleep(2.5)       
            await write_report_progress(args.rundir, 100, 100)
                            
            pool.close()
            pool.join()
            print('All processes completed', flush=True)

            with lock_l:
                for l in logm:
                    log.append(l)

            with lock_r:
                completed_jobs = [reps[key] for key in sorted(reps.keys())]

            test_variant_seq()

            if args.amplicons:
                args.outfn = 'amplicon_' + args.outfn
                args.variants = 'amplicon_' + args.variants
            
            with open(os.path.join(args.rundir,args.outfn), "wt", buffering=1) as dstfd,\
                    open(os.path.join(args.rundir,args.variants), 'wt', buffering=1) as varfd:
                print(f"Opening {args.outfn} for results", flush=True)
                print(f"Opening {args.variants} for variant sequences", flush=True)

                dst = csv.writer(dstfd, dialect="unix", quoting=csv.QUOTE_ALL)
                hdrres1 = ("readCount", "cleanCount", "mergeCount")
                hdrres2 = ("seqCount", "seqName", "efficiency", "otherCount", "otherName")
                complete_row_hdr = tuple((x for xs in (hdr, hdrres1, hdrres2) for x in xs))
                primer_col = [i for i,col in enumerate(complete_row_hdr) if col=='primer'][0]
                dst.writerow(complete_row_hdr)
                for i,job in enumerate(completed_jobs):
                    try:
                        dst.writerow(job)
                    except Exception as exc:
                        print(f'{exc=}', flush=True)
                    #print(i, job, flush=True)
                    var_count_entries = job[-2].split(';')
                    var_name_entries = job[-1].split(';')
                    primer_name = job[primer_col]
                    ### Genotyping only want to report the sample and the count
                    #print(f'{primer_name=} {var_count_entries=} {var_name_entries=}')
                    for var_count, var_name in zip(var_count_entries, var_name_entries):
                        if var_name.startswith('other'):
                            continue
                        # only print out actual variants, not off-target hits
                        if '//' in var_name:
                            var_row_name = f'>Sample:{i+1}'   # ;Primer:{primer_name}'
                            var_row_name += f';count:{var_count}'
                            print(var_row_name, file=varfd)
                            print(get_variant_seq(var_name, id_seq), file=varfd)
                            
                dstfd.flush()
                varfd.flush()
    except Exception as exc:
        print(f'Error: {exc}', flush=True)
        log.append(f'Error: {exc}')

    end_time = datetime.datetime.now()
    msg = f"End: {end_time} took: {end_time - start_time}"
    print(msg, flush=True)
    log.append(msg)
    await write_log(log, os.path.join(args.rundir,args.logfn))
     

   
if __name__=="__main__":
    """

    """
    parser = argparse.ArgumentParser(description="NGS Reporting Program")
    parser.add_argument("-d", "--debug", action="store_true", help="more reporting/output for debugging purposes")
    parser.add_argument("-t", "--targets", default="targets.fa", help="file of targets in FASTA format (default=targets.fa)")
    parser.add_argument('-a','--amplicons', nargs='*', help='Enable amplicon matching by including plateIDs, default is to match genotypes')
    parser.add_argument('-P', '--primer_assayfam', default="primers.csv", help='file of primer to assay family mappings')
    parser.add_argument('-o','--outfn', default='results.csv', help='Name of output file (CSV format)')
    parser.add_argument('-v','--variants', default='variant_seqs.fa', help='Name of variant sequences file (FASTA format)'),
    parser.add_argument('-n','--ncpus', type=int, default=os.cpu_count()-1, help='Number of processes to run simultaneously, default=number of CPUs in system - 1')
    parser.add_argument('-l','--logfn', default='match.log', help='Name of logging file (default=match.log)')
    parser.add_argument('-r','--rundir', required=True, help='Path to experiment folder')
    parser.add_argument('-i','--identity',type=float,default=0.9,help='Proportional score for inexact matching')
    parser.add_argument('-m','--mincov', type=int, default=5, help='Do not match unique sequences with less than this many reads coverage, default 50')
    parser.add_argument('-p','--minprop', type=float, default=0.1, help='Do not match unique sequences '+\
            'with less than this proportion of the total number of exact matched on-target reads, default 0.2. Must be between 0.0 and 1.0')
    parser.add_argument('--inexact', action="store_true", help="enable inexact matching")
    parser.add_argument('-x','--exhaustive',action='store_true',help='Try to match every sequence, '+\
            'no matter how few counts. Ignores --minseqs and --minprop')
    parser.add_argument('-M','--margin',type=float,default=0.9,help="Sequences must be this proportion of the reference seq length")
    parser.add_argument('-s','--stagefile', default="Stage3.csv", help="Name of the NGS genotyping Stage 3 file (default=Stage3.csv)")
    args = parser.parse_args()
    in_error = False
    print(f"{args=}", file=sys.stderr)
    lock_path = os.path.join(args.rundir,"ngsgeno_lock") 
    if not os.path.exists(lock_path):
        try:
            with open(lock_path,"wt"):
                asyncio.run(write_report_progress(args.rundir, 0, 0))  # set this up asap
                asyncio.run(main(args))
            print('Completed regular execution')
        except Exception as exc:
            print(f'Completed with exception {exc}', flush=True)
        finally:
            if os.path.exists(lock_path):
                os.remove(lock_path)
            for fn in list(Path(args.rundir).glob('match_progress_*')):
                os.remove(fn)
    else:
        print("Analysis already running", file=sys.stderr)