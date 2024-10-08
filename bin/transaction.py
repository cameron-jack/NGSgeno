#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from copy import deepcopy
from pathlib import Path
import os
import sys
import functools
import jsonpickle

from stutil import m

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
    True if either exp.uploaded_files['_upload_pending'] or exp.pending_steps have entries (both are sets)
    """
    if exp.uploaded_files['_upload_pending'] or exp.pending_steps:
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
                exp.log(f'file generation {t} already performed in this stage of the pipeline', level='critical')
                return False
            exp.pending_steps[t] = deepcopy(transactions[t])
            exp.log(f'Adding generated file {t} to pending pipeline stage history', level='info')
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
            #exp.log(f'file path {final_path} already ' +
            #        'exists and will be overwritten if pending changes are accepted', level='warning')
        else:
            for step in exp.reproducible_steps:
                if step is None or len(step) == 0:
                    continue
                # each step is dict['filenames'] = {PID: {well:change}}
                if final_path in step:
                    #exp.log(f'file path {final_path} already ' +
                    #        'exists and will be overwritten if pending changes are accepted', level='warning')
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
        exp.log(f'could not clear pending transactions, possbile locked file. {exc}', level='critical')
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
        pending_files_hanging = list(Path(exp.get_exp_dn()).glob('pending_*'))+\
                list(Path(exp.get_exp_dn(subdir='uploads')).glob('pending_*'))
        for p in pending_files_hanging:
            os.remove(p)
            #print(f'removing pending file {p}', file=sys.stderr)
    except Exception as exc:
        m(f'Could not clear pending transactions, possbile locked file. {exc}', level='error')
        return False
    return True

    
def accept_pending_transactions(exp, file_name=None, caller_id=None):
    """ 
    Look up the keys from exp.pending_steps in exp.reproducible_steps then
        remove all entries matching and following this then
        append the pending steps, rename files, and clear pending steps
    If file_name then use only this key and ignore others
        
    Returns True on success
    """
    #print(f"accept_pending_transactions for {exp.pending_steps.keys()=}", file=sys.stderr)
    if not exp.pending_steps:
        m("there are no pending transactions to record", level='debug', dest=('noGUI',))
        return True  # It didn't actually fail
    ps_keys = exp.pending_steps.keys()
    if file_name:
        ps_keys = [ps for ps in ps_keys if file_name == ps_keys]
    elif not ps_keys:
        m(f'expected {file_name} in {exp.pending_steps.keys()}', level='error', dest=('noGUI',), caller_id=caller_id)
        return False
    MAX_STAGES=len(exp.reproducible_steps)
    clashing_index = MAX_STAGES
    saw_a_clash = False
    #print(f"{exp.reproducible_steps=}", file=sys.stderr)
    for i,step in enumerate(exp.reproducible_steps):
        for dest in step:
            dp = convert_final_to_pending(exp,dest)
            #print(f"{dp=} {ps_keys=}", flush=True, file=sys.stderr)
            if dp in ps_keys:
                if i < clashing_index:
                    clashing_index = i
                    saw_a_clash = True
                    break
    if not saw_a_clash: # maybe no entries in reproducible_steps yet, but files exist for some reason
        exp.reproducible_steps.append({})
        for ps in ps_keys:
            exp.reproducible_steps[-1][ps] = exp.pending_steps[ps].copy()
            fs = convert_pending_to_final(exp, ps)
            # check if a file already exists with this name
            if os.path.exists(fs):
                try:
                    os.remove(fs)
                except Exception as exc:
                    m(f'Could not remove existing file {fs}, {exc}', level='error', caller_id=caller_id)
                m(f'removed existing file {fs}', level='info', dest=('noGUI',))
            try:
                os.rename(ps, fs)
            except Exception as exc:
                m(f'failed to rename pending file {ps} to {fs}, {exc}', level='error', dest=('noGUI',), caller_id=caller_id)
                return False
            m(f'renamed pending file {ps} to {fs}', level='success', dest=('noGUI',))
        exp.pending_steps = {}
        return True
    else:
        # find the clashing location, nuke everything from that index and add the pending_steps instead
        affected_files = set()
        affected_pids = set()  # currently we're leaving this, but if plates are generated then they need to be removed too
        for step in exp.reproducible_steps[clashing_index:]:
            for fn in step.keys():
                affected_files.add(fn)
                for pid in step[fn]:
                    affected_pids.add(pid)
        for af in affected_files:
            success = exp.del_file_record(af)
            if success:
                m(f'obsolete tracked pipeline file {af} removed', level='info', dest=('noGUI',), caller_id=caller_id)
            else:
                try:
                    os.remove(af)
                except Exception as exc:
                    m(f'could not delete obsolete file {af}, {exc}', level='error', dest=('noGUI',), caller_id=caller_id)
                    return False
                m(f'obsolete tracked pipeline file {af} removed', level='info', dest=('noGUI',), caller_id=caller_id)
            
        exp.reproducible_steps = exp.reproducible_steps[:clashing_index]
        exp.reproducible_steps.append({})
        for ps in ps_keys:
            exp.reproducible_steps[-1][ps] = exp.pending_steps[ps].copy()
            fs = convert_pending_to_final(exp, ps)
            try:
                os.rename(ps, fs)
            except Exception as exc:
                m(f'failed to rename pending file {ps} to {fs}, {exc}', level='error', caller_id=caller_id)
                return False
            m(f'renamed pending file {ps} to {fs}', level='info', dest=('noGUI',))
        exp.pending_steps = {}
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
        m(f'{PID} not found in plate records', level='error')
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
                    if well in mod_plate:
                        if 'volume' in mod_plate[well]:
                            mod_plate[well]['volume'] += stage[fn][PID][well]
                    #except:
                    #    print(f'{mod_plate[well]=} {stage[fn][PID][well]=} {fn=} {PID=} {well=}', file=sys.stderr)
                    #    exit()
    if transactions:
        for t in transactions:
            if PID in transactions[t]:
                for well in transactions[t][PID]:
                    if well in mod_plate:
                        if 'volume' in mod_plate[well]:
                            mod_plate[well] += transactions[t][PID][well]
    return mod_plate


### Manage file to plate and plate to file mappings

def find_affected(exp, files, plates):
    """
    Recurse through files and plates, catching all that rely on each other
    If a plate needs a file, then the file to plate mapping should reflect this and vice versa
    This should be a better way of handling transactions than simply winding back all files 
    produced so far, but removes the idea of a truly linear pipeline
    files and plates are sets with elements found in:
    exp.denied_assays[0] is for files
    exp.denied_primers[0] is for plates
    """
    mapped_files = exp.denied_assays[0]
    mapped_plates = exp.denied_plates[0]
    extra_files = set()
    extra_plates = set()
    for f in files:
        if f in mapped_files:
            for p in mapped_files[f]:
                if p not in plates:
                    extra_plates.add(p)
                    plates.add(p)
    for p in plates:
        if p in mapped_plates:
            for f in mapped_plates[p]:
                if f not in files:
                    extra_files.add(f)
                    files.add(f)
    if extra_files or extra_plates:
        returned_files, returned_plates = find_affected(exp, extra_files, extra_plates)
        for rf in returned_files:
            files.add(rf)
        for rp in returned_plates:
            plates.add(rp)
    return files, plates