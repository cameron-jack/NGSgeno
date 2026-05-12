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

from util import unguard, padwell


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
        #print(f'reconstruct_sequence() {i=} {p=}', flush=True)
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
    
    
    
