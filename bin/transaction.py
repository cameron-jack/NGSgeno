#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from copy import deepcopy
from pathlib import Path
import os
import sys
import functools
import jsonpickle

"""
@created: Nov 2021
@author: Cameron Jack ANU Bioinformatics Consultancy, 2019-2021

transactions.py contains the transactional protection code for files and plates to ensure repeatibility.

Experiment.reproducible_steps is built around atomic file writes acting as checkpoints, which are initially
prefixed with 'pending_' until accepted (either because of no blockers, or by the user). This operation also prefixes
any plates that are being written with "pending_". When file changes are accepted this prefix is removed.

If a plate is modified during a file operation, then the change in state (e.g. volume) is held in the 
reproducible_steps dictionary entry.

The user shouldn't ever see this directly, it should handle everything silently. The interfaces should
always take the experiment and a queue for clashes that the interface can then deal with. No operation is allowed
to replace an existing file or plate that is prefixed as "pending". These entries must be cleared first.
"""

def transact(id, is_plate=False):
    """
    take a file path or plate_id e.g. \example\myfile.csv and convert to 
    pending name eg \example\pending_myfile.csv
    """
    #print("convert id to pending id", file=sys.stderr)
    #print(f"{id=}", file=sys.stderr)
    if is_plate:
        pending_id = 'pending_' + id
    else:
        p = Path(id)
        parent = p.parent
        file_name = str(p.name)
        pending_name = "pending_" + file_name
        pending_id = str(parent / pending_name)

    #print(f"transact:{id=} {p=} {file_name=} {pending_name=} {pending_id=}", file=sys.stderr)
    return pending_id


def untransact(pending_id, is_plate=False):
    """
    take a pending path \example\pending_myfile.csv and convert to final name
    eg \example\myfile.csv
    """
    #print(f"{id=}", file=sys.stderr)
    if is_plate:
        id = pending_id[len('pending_'):]
    else:
        p = Path(pending_id)
        parent = p.parent
        file_name = str(p.name)
        final_name = file_name[len('pending_'):]  # cut off the leading "pending_"
        id = str(parent / final_name)
        #print(f'untransact:{p=} {parent=} {file_name=} {id=}')
    #print(f"{id=}", file=sys.stderr)
    return id

 
def is_pending(exp):
    """
    True if either exp.pending_uploads or exp.pending_steps have entries (both are sets)
    """
    if exp.pending_uploads or exp.pending_steps:
        return True
    return False

def get_pids_affected_by_file(exp, fn):
    """
    Return a list of all plateIDs that require fn
    """
    pass

def get_pids_affected_by_pid(exp, pid):
    """
    Return a list of all plateIDs that require pid
    """
    pass

def get_files_affected_by_file(exp, fn):
    """
    Return a list of all files that require fn
    """
    pass

def get_pids_affected_by_file(exp, fn):
    """
    Return a list of all plateIDs that require fn
    """
    pass


def add_pending_transactions(exp, transactions):
        """                 
        Add to the current pending steps for a "step" of file generation. 
        It will fail (return False) if it attempts to overwrite an existing stored change.
        You can call this multiple times and the results will be combined into a single "step"
        """
        if not transactions:
            return True
        if exp.pending_steps is None:
            exp.pending_steps = {}
        for t in transactions:
            if t in exp.pending_steps:
                exp.log(f'Critical: file generation {t} already performed in this stage of the pipeline')
                return False
            exp.pending_steps[t] = deepcopy(transactions[t])
            exp.log(f'Info: Adding generated file {t} to pending pipeline stage history')
            #print('These are the pending transactions', exp.pending_steps, file=sys.stderr)
        return True

def convert_pending_to_final(exp, pending_name):
    """
    take a pending path \example\pending_myfile.csv and convert to final name
    eg \example\myfile.csv
    """
    #print(f"{pending_name=}", file=sys.stderr)
    p = Path(pending_name)
    parent = p.parent
    file_name = str(p.name)
    final_name = file_name[len('pending_'):]  # cut off the leading "pending_"
    final_path = str(parent / final_name)
    #print(f"{final_path=}", file=sys.stderr)
    return final_path

def convert_final_to_pending(exp, final_name):
    """
    take a final path \example\myfile.csv and convert to pending name
    eg \example\pending_myfile.csv
    """
    #print("convert final to pending", file=sys.stderr)
    #print(f"{final_name=}", file=sys.stderr)
    p = Path(final_name)
    parent = p.parent
    file_name = str(p.name)
    pending_name = "pending_" + file_name
    pending_path = str(parent / pending_name)
    #print(f"{pending_path=}", file=sys.stderr)
    return pending_path

def check_for_clashing_transactions(exp, filenames=None, pids=None):
    """
    If not filenames or pids are provided, go through self.pending_steps and check whether any of the 
    transacted filenames clash with existing filenames in self.reproducible_steps, noting the filenames in
    a set called self.pending_steps['clashing_filenames'] and pids in self.pending_steps['clashing_pids']
    If filenames is a list then find all filenames and pids that clash
    if pids is a list then find all filenames and pids that clash
    Clashes get added to temporary filename and pid sets and if these entries aren't in the existing sets
    we call this function recursively
    Once recursions are finished, all clashes are noted in self.pending_steps['clashing_filenames'] and
    self.pending_steps['clashing_pids']
    """
    if exp.pending_steps is None:
        exp.pending_steps = {}
    clashing_filenames = set()
    clashing_pids = set()

    # we'll check these again at the end before adding them to self.pending_steps['clashing_XXX']
    local_clashing_filenames = set()
    local_clashing_pids = set() 

    # go through all pending steps
    if not filenames and not pids:
        for fn in exp.pending_steps:
            for rs in exp.reproducible_steps:
                for rs_fn in rs:
                    if rs_fn == exp.convert_pending_to_final(fn):
                        local_clashing_filenames.add(fn)
                        for pid in rs['rs_fn']:
                            local_clashing_pids.add(pid)

    # check for clashing filenames, and also note the PIDs involved
    if filenames:
        for fn in filenames:
            for rs in exp.reproducible_steps:
                for rs_fn in rs:
                    if rs_fn == convert_pending_to_final(exp,fn):
                        local_clashing_filenames.add(fn)
                        for pid in rs['rs_fn']:
                            local_clashing_pids.add(pid)

    # check for clashing PIDs that either match or are required by other plates
    if pids:
        for PID in exp.plate_location_sample:
            for pid in pids:
                if pid == PID:
                    local_clashing_pids.add(pid)
                elif pid in exp.plate_location_sample[PID]['source']:
                    local_clashing_pids.add(PID)
        
    # we'll recursively call this function on these
    extra_clashing_filenames = set()
    extra_clashing_pids = set()

    # add them now to self.pending_steps to avoid infinite loop
    for fn in local_clashing_filenames:
        if fn not in clashing_filenames:
            extra_clashing_filenames.add(fn)
        clashing_filenames.add(fn)
            
    for pid in local_clashing_pids:
        if pid not in clashing_pids:
            extra_clashing_pids.add(pid)
        clashing_pids.add(pid)

    # call recursively if required
    if extra_clashing_filenames or extra_clashing_pids:
        returned_clashing_filenames, returned_clashing_pids =\
                check_for_clashing_transactions(exp, filenames=extra_clashing_filenames, pids=extra_clashing_pids)
        for fn in returned_clashing_filenames:
            clashing_filenames.add(fn)
        for pid in returned_clashing_pids:
            clashing_pids.add(pid)
    return clashing_filenames, clashing_pids

def clashing_pending_transactions(exp):
    """
    We need the user to know if accepting pending transactions will result in overwriting an existing
    file. We return the list of all clashing filepaths
    """
    clashes = []
    if exp.pending_steps is None:
        #print('no pending steps')
        return clashes
    for transaction in exp.pending_steps:
        final_path = convert_pending_to_final(exp,transaction)
        #print('Pending and final paths: ', str(transaction), str(final_path), file=sys.stderr)
        if Path(final_path).exists():
            clashes.append(final_path)
            exp.log(f'Warning: file path {final_path} already' +
                    'exists and will be overwritten if pending changes are accepted')
        else:
            for step in exp.reproducible_steps:
                if step is None or len(step) == 0:
                    continue
                # each step is dict['filenames'] = {PID: {well:change}}
                if final_path in step:
                    exp.log(f'Warning: file path {final_path} already' +
                            'exists and will be overwritten if pending changes are accepted')
                    clashes.append(final_path)                  
    return clashes
    
def clashing_pending_transaction(exp, file_upload):
    clashes = clashing_pending_transactions()
    file_path = exp.get_exp_fn(file_upload.name, trans=False)
    if file_path in clashes:
        return file_path
    return False

def clear_pending_transaction(exp, file_upload):
    """ provided file path shouldn't be prefixed 'pending_' """
    if exp.pending_steps is None:
        return True
    try:
        pf = exp.get_exp_fn(file_upload, trans=True)
        if pf not in exp.pending_steps:
            return False
        if Path(pf).exists():
            print(f'Clearing pending file: {pf}', file=sys.stderr)
            os.remove(pf)
        exp.pending_steps.pop(pf)
    except Exception as exc:
        exp.log(f'Critical: Could not clear pending transactions, possbile locked file. {exc}')
        return False
    return True


def clear_pending_transactions(exp):
    """
    Nukes all pending steps and files (presumably the user declined to keep them)
    """
    if exp.pending_steps is None:
        return True
    try:
        for transaction in exp.pending_steps:
            if Path(transaction).exists():
                os.remove(transaction)
        exp.pending_steps = None
        pending_files_hanging = Path(exp.get_exp_dn()).glob('pending_*')+\
                Path(exp.get_exp_dn(subdir='uploads')).glob('pending_*')
        for p in pending_files_hanging:
            os.remove(p)
            #print(f'removing pending file {p}', file=sys.stderr)
    except Exception as exc:
        exp.log(f'Critical: Could not clear pending transactions, possbile locked file. {exc}')
        return False
    return True
    
def accept_pending_transaction(exp, file_name):
    """
    Replace the existing file with the one that is pending
    """
    print(file_name)
    pending_file = exp.get_exp_fn(filename=file_name, trans=True)
    print(pending_file)
    if not Path(pending_file).exits():
        exp.log(f"Error: pending file {pending_file} does not exist, reverting to original")
        if pending_file in exp.pending_uploads:
            exp.pending_uploads.remove(pending_file)
        if pending_file in exp.pending_steps:
            del exp.pending_steps[pending_file]
        return False
    if not exp.pending_steps:
        exp.log("Warning: there are no pending transactions to record")
        return True  # It didn't actually fail
    elif pending_file in exp.pending_steps:
        #print(f'Pending steps: {exp.pending_steps}, reproducible steps: {exp.reproducible_steps}', 
            #       file=sys.stderr)
        clashes = exp.clashing_pending_transactions()
        #print(f'Clashes seen: {clashes}', file=sys.stderr)
        final_path = exp.convert_pending_to_final(pending_file)
        #print(f'File name: {final_path}', file=sys.stderr)

        #if file already exists, remove original from files and reproducible_steps
        if final_path in clashes:
            if Path(final_path).exists():
                print(f"Removing the old file: {str(final_path)}", file=sys.stderr)
                os.remove(final_path)

            #Need to remove entries from reproducible steps that contain the same plate id?
            for i, step in enumerate(exp.reproducible_steps):
                if final_path in step:
                    del exp.reproducible_steps[i]

        p = Path(pending_file)
        if not p.exists():
            exp.log('Warning: file does not exist')
            return False
        os.rename(pending_file, final_path)
        print(f'Accepted file: {final_path}', file=sys.stderr)

        this_step = {}
        record = exp.pending_steps[pending_file]
        this_step[final_path] = record
        #print(f'{this_step}', file=sys.stderr)
        exp.reproducible_steps.append(this_step)
        exp.pending_steps.pop(pending_file)

    return True

def accept_pending_transactions(exp):
    """ 
    Look up the keys from exp.pending_steps in exp.reproducible_steps then
        remove all entries matching and following this then
        append the pending steps, rename files, and clear pending steps

    Returns True on success
    """
    print(f"accept_pending_transactions for {exp.pending_steps.keys()=}", file=sys.stderr)
    if not exp.pending_steps:
        exp.log("Warning: there are no pending transactions to record")
        return True  # It didn't actually fail
    ps_keys = exp.pending_steps.copy()
    for ps in ps_keys:
        if not Path(ps).exists():
            exp.log(f"Error: pending file {ps} does not exist, reverting to original")
            if ps in exp.pending_uploads:
                exp.pending_uploads.remove(ps)
            if ps in exp.pending_steps:
                del exp.pending_steps[ps]
            if ps in exp.uploaded_files:
                exp.del_file_record(ps)

    if len(exp.pending_steps.keys()) == 0:
        return True

    clashes = clashing_pending_transactions(exp)
    #print(f"Clashes seen {clashes=}", file=sys.stderr)
    if len(clashes) == 0:
        for transaction in exp.pending_steps:
            p = Path(transaction)
            if not p.exists():
                exp.log(f'Warning: {str(p)} not found')
                continue
            final_path = convert_pending_to_final(exp,transaction)
            os.rename(p, final_path)
            #print(f"No clash {str(p)=} {str(final_path)=}", file=sys.stderr)
    else:
        MAX_STAGES=99999
        clashing_index = MAX_STAGES
        #print(f"{exp.reproducible_steps=}", file=sys.stderr)
        for i,step in enumerate(exp.reproducible_steps):
            for dest in step:
                dp = convert_final_to_pending(exp,dest)
                #print(f"{dp=}", file=sys.stderr)
                if dp in exp.pending_steps:
                    if i < clashing_index:
                        clashing_index = i
                        break
        #print(f"{clashing_index=}", file=sys.stderr)
        if clashing_index == MAX_STAGES:  # this should NEVER happen
            exp.log(
                    'Critical: pipeline detects clashing' +
                    f'transaction {clashing_index=} for {exp.pending_steps=}') 
            return False

        # keep everything prior to the clash, then add on the pending steps
        remove_these_steps = exp.reproducible_steps[clashing_index:]
        #print(f"{remove_these_steps=}", file=sys.stderr)
        for step in remove_these_steps:
            for fp in step:
                if fp in exp.uploaded_files:
                    del exp.uploaded_files[fp]
                if Path(fp).exists():
                    #print(f"removing the original file: {str(fp)}", file=sys.stderr)
                    os.remove(fp)
        # rename pending filepaths
        for transaction in exp.pending_steps:
            p = Path(transaction)
            if not p.exists():
                exp.log(f'Warning: {str(p)} not found')
                continue
            final_path = convert_pending_to_final(exp,transaction)
            #print(f"Renaming {str(p)=} to {str(final_path)=}", file=sys.stderr)
            os.rename(p, final_path)
            op = str(p)
            if op in exp.uploaded_files:
                exp.uploaded_files[final_path] = exp.uploaded_file[op].copy()
                del exp.uploaded_file[op]
        exp.reproducible_steps = exp.reproducible_steps[0:clashing_index]
    if exp.pending_steps is None:
        exp.save()
        return True
    this_step = {}
    for ps in exp.pending_steps:
        if ps is None:
            continue
        final_name = convert_pending_to_final(exp,ps)
        record = exp.pending_steps[ps]
        this_step[final_name] = record
        op = str(final_name).split('/')[-1]
        if exp.pending_steps[ps]:
            exp.uploaded_files[op] = {'plates':exp.pending_steps[ps].keys()}
        else:
            exp.uploaded_files[op] = {'plates':[]}
    exp.reproducible_steps.append(this_step)
    for fn in this_step:
        if fn in exp.pending_steps:
            exp.pending_steps.remove(fn)
    exp.save()
    return True


def enforce_file_consistency(exp):
    """ 
    Iterate over all expected files and ensure that they are present. If they aren't, remove these steps and any that follow.
    """
    pids_to_delete = set()
    files_to_delete = set()
    for i, step in enumerate(exp.reproducible_steps):
        for f in step:
            if not Path(f).exists():
                files_to_delete.add(f)
                if step[f] is None:
                    continue
                for pid in step[f]:
                    pids_to_delete.add(pid)
    # clear out files and look for pipeline breakages from missing sections
    pipe_broken_step = None
    for i, step in enumerate(exp.reproducible_steps):
        for f in files_to_delete:
            if f in exp.uploaded_files:
                del exp.uploaded_files[f]
            if f in step:
                del step[f]
            if len(step) == 0:
                pipe_broken_step = i
                break
        if pipe_broken_step is not None:
            break
    if pipe_broken_step is not None:
        exp.reproducible_steps = exp.reproducible_steps[:pipe_broken_step]

    # now clean out pids 
    for i, step in enumerate(exp.reproducible_steps):
        for f in step:
            for pid in pids_to_delete:
                if pid in step[f]:
                    del exp.reproducible_steps[i][f][pid]


def get_plate(exp, PID, transactions=None):
    """ 
    Look up a plate ID in exp.plate_location_sample and then in exp.reproducible_steps
    Apply all changes seen in exp.reproducible steps
    Return the modified plate
    PID should be a guarded plate ID 
    """
    if PID not in exp.plate_location_sample:
        exp.log(f'Error: {PID} not found in plate records')
        return None
    mod_plate = deepcopy(exp.plate_location_sample[PID])
    # exp.reproducible_steps is a list, so stage will be a dictionary
    for stage in exp.reproducible_steps:
        if stage is None:
            continue
        # each stage is dict['filenames'] = {PID: {well:change}}
        for fn in stage:
            if stage[fn] is None:
                continue
            if PID in stage[fn]:
                for well in stage[fn][PID]:
                    if 'volume' in mod_plate[well]:
                        mod_plate[well]['volume'] += stage[fn][PID][well]
                    #except:
                    #    print(f'{mod_plate[well]=} {stage[fn][PID][well]=} {fn=} {PID=} {well=}', file=sys.stderr)
                    #    exit()
    if transactions:
        for t in transactions:
            if PID in transactions[t]:
                for well in transactions[t][PID]:
                    if 'volume' in mod_plate[well]:
                        mod_plate[well] += transactions[t][PID][well]
    return mod_plate
