#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: 22 Feb 2020
@author: Cameron Jack, Gabrielle Ryan, ANU Bioinformatics Consultancy, JCSMR, Australian National University 

Define an experiment - all parts and steps that go into a NGS Genotyping experiment will be held here (all state).
The GUI interacts with a single Experiment object at one time. Methods are called on this to activate pipeline
functionality. The Experiment then deals directly with the pipeline logic.

Stores information in experiment.ini - replaces the old TemplateFiles class
Stores both the inventory - all the plates and samples - but also the options used at each step to enable reproduction of a run
"""

import os
import sys
import csv
import jsonpickle
from itertools import combinations
from io import StringIO, BytesIO
import datetime
import json
from Bio import SeqIO
from copy import deepcopy
from math import ceil, floor
from pathlib import Path
import inspect  

import pandas as pd

## Ugly hack to find things in ../bin
#from pathlib import Path
#file = Path(__file__).resolve()
#sys.path.append(os.path.join(file.parents[1],'bin'))
try:
    import bin.db_io as db_io
except ModuleNotFoundError:
    import db_io
try:
    import bin.file_io as file_io
except ModuleNotFoundError:
    import file_io
try:
    import util
    
except ModuleNotFoundError:
    import bin.util as util


EXP_FN = 'experiment.json'      

# Defines content of experiment.ini which stores the current pipeline configuration
class Experiment():
    """

    All pipeline inputs go here. Logs are housed here. Everything the pipeline needs to operate can be found here.
    All plates are recorded in self.plate_location_sample which is defined by {barcode/name:{well:{info}}}
    Info must have 'purpose' in {'sample', 'primer', 'index', 'amplicon', 'dna', 'pcr', 'taq_water'}
    We need functions to link results back to sample
    samples are their own dict structure, separate from plates, which have the form plate_locations[pid] = {well:sample_id}

    wells are always unpadded. Echo needs unpadded wells, possibly Miseq?

    Do we read raw inputs and then offer the ability to manually correct them? Maybe only the manifests...
    
    How to track usage of wells and not accidentally drain a well twice if we regenerate earlier files?
    Transaction tracking:
        - "generate" methods, which create fixed outputs (and usually modify state) should create a "transaction".
        - self.reproducible_steps - [] is a list (stack) of transactions
        - Only one instance of a generate method's use can exist in the list for the same set of outputs
        - Because it is ordered, we can "rewind" to this location.
        - Transactions store all plate modifications (i.e. +/- volume changes)
        - These modifications are applied to a plate record before any operations are carried out with these plates

    Once sequence data is available, the earlier stages of the pipeline should be locked.
        
    """
    def __init__(self, name):
        """ set up defaults, an experiment name is required, should match folder name less then run_ prefix """
        self.name = name
        self.description = ''
        # Use "locked" for everything up to the allele calling stage. Once sequences are available you shouldn't
        # be allowed to change anything in the pipeline itself
        self.locked = False  # This is meant to prevent modification to the experiment when True
        self.unassigned_plates = {1:'', 2:'', 3:'', 4:''}  # plate_id:info - we can import these and then let users edit the results - they aren't checked until added to a plate set
        self.dest_sample_plates = {}  # {dest_pid:[4 sample plate ids]}
        self.plate_location_sample = {}  # pid:{well:sample_dict}  # use this for everything!
        self.deleted_plates = {}  # pid:[list of deleted plates with this pid]
        self.stage = 'new'  # one of 'new', 'nimbus', 'primers', 'barcodes', 'miseq', 'genotyping', 'review'  <-- probably unnecessary
        ###
        # sample_dict is of the form:
        # All: {"barcode":str, "sex":str, mouseId:int (this could be mouseKey or mouseNumber in Rodentity
        # Rodentiy: "alleles":[{alleleKey, name, symbol, options:[], assays:[{assayKey, name, method}]}]
        # Musterer: "strain", "assays":[{name, method, value(g/t), options:[]}] <- rename to must_assays
        # Custom: "assays": [str]
        # Standardised: "assays": [str], "assayFamilies": [str] <- musterer assay[name], rodentity allele[assay[name]]
        # The specific import function should be responsisble for creating the "assays" and "assayFamilies" fields.
        # Note the ugly camel case for assayFamilies - to match headers with this case type
        ###
        self.denied_assays = []
        self.denied_primers = []
        self.assay_synonyms = {}  # {source:{reference:alternative}}
        self.primer_assay = {}  # {source:{mapping of primer to assay}}
        self.reference_sequences = {}  # {source:{name:seq}} mapping of sequence name to sequence
        ###
        # All generate_XXX() functions produce one or more files
        # a transaction record in self.reproducible_steps consists of a dictionary of form:
        #   {sorted output filenames:{plate_id:[modifications to the plate]}
        # This is a set of modifications to our saved plates which are applied before we do anything further with these plates
        # On loading an experiment we should go through the steps and make sure all the expected output files are present.
        # In order to access a plate with all its modifications we need a function for looking up plates and applying 
        #   modifications during the pipeline. We'll call this Experiment.get_plate(PID)
        ###
        self.uploaded_files = {}
        self.reproducible_steps = []  # a strict subset of the log that only includes working instructions
        self.pending_steps = None # reproducible steps that are awaiting user approval to replace existing steps
        self.log_entries = []  # use self.log(message, level='') to add to this
        self.log_time = datetime.datetime.now()  # save the log after a certain time has elapsed, to avoid too much IO
        self.transfer_volumes = {'DNA_VOL':util.DNA_VOL, 'PRIMER_VOL':util.PRIMER_VOL, 'PRIMER_TAQ_VOL':util.PRIMER_TAQ_VOL,
                'PRIMER_WATER_VOL':util.PRIMER_WATER_VOL, 'INDEX_VOL':util.INDEX_VOL, 'INDEX_TAQ_VOL':util.INDEX_TAQ_VOL,
                'INDEX_WATER_VOL':util.INDEX_WATER_VOL}  # load the defaults from util.py
        

    def __setstate__(self, state):
        for k, v in state.items():
            setattr(self, k, v)

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __repr__(self) -> str:
        return str(self.__dict__)

    def lock(self):
        self.log('Info: Locking experiment. No modification allowed to plates while lock remains')
        #self.locked = True
        #self.save()

    def unlock(self):
        self.log('Info: Unlocking experiment. Modification is now possible. There should be good reason for this!')
        #self.locked = False
        #self.save()

    ### transactions/reproducible steps

    def add_pending_transactions(self, transactions):
        """                 
        Add to the current pending steps for a "step" of file generation. 
        It will fail (return False) if it attempts to overwrite an existing stored change.
        You can call this multiple times and the results will be combined into a single "step"
        """
        if not transactions:
            return True
        if self.pending_steps is None:
            self.pending_steps = {}
        for t in transactions:
            if t in self.pending_steps:
                self.log(f'Critical: file generation {t} already performed in this stage of the pipeline')
                return False
            self.pending_steps[t] = deepcopy(transactions[t])
            self.log(f'Info: Adding generated file {t} to pending pipeline stage history')
            #print('These are the pending transactions', self.pending_steps, file=sys.stderr)
        return True

    def convert_pending_to_final(self, pending_name):
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

    def convert_final_to_pending(self, final_name):
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


    def clashing_pending_transactions(self):
        """
        We need the user to know if accepting pending transactions will result in overwriting an existing
        file. We return the list of all clashing filepaths
        """
        clashes = []
        if self.pending_steps is None:
            #print('no pending steps')
            return clashes
        
        for transaction in self.pending_steps:
            final_path = self.convert_pending_to_final(transaction)
            #print('Pending and final paths: ', str(p), str(final_path), file=sys.stderr)
            if Path(final_path).exists():
                clashes.append(final_path)
                self.log(f'Warning: file path {final_path} already exists and will be overwritten if pending changes are accepted')
            else:
                for step in self.reproducible_steps:
                    if step is None or len(step) == 0:
                        continue
                    # each step is dict['filenames'] = {PID: {well:change}}
                    if final_path in step:
                        self.log(f'Warning: file path {final_path} already exists and will be overwritten if pending changes are accepted')
                        clashes.append(final_path)                  
        return clashes


    def clear_pending_transactions(self):
        """
        Nukes all pending steps and files (presumably the user declined to keep them)
        """
        if self.pending_steps is None:
            return True
        try:
            for transaction in self.pending_steps:
                if Path(transaction).exists():
                    os.remove(transaction)
            self.pending_steps = None
            pending_files_hanging = Path(self.get_exp_dir()).glob('pending_*')
            for p in pending_files_hanging:
                os.remove(p)
                #print(f'removing pending file {p}', file=sys.stderr)
        except Exception as exc:
            self.log(f'Critical: Could not clear pending transactions, possbile locked file. {exc}')
            return False
        return True

    def accept_pending_transactions(self):
        """ 
        Look up the keys from self.pending_steps in self.reproducible_steps then
            remove all entries matching and following this then
            append the pending steps, rename files, and clear pending steps

        Returns True on success
        """
        #print(f"accept_pending_transactions for {self.pending_steps.keys()=}", file=sys.stderr)
        if not self.pending_steps:
            self.log("Warning: there are no pending transactions to record")
            return True  # It didn't actually fail

        clashes = self.clashing_pending_transactions()
        #print(f"Clashes seen {clashes=}", file=sys.stderr)
        if len(clashes) == 0:
            for transaction in self.pending_steps:
                p = Path(transaction)
                if not p.exists():
                    self.log(f'Warning: {str(p)} not found')
                    continue
                final_path = self.convert_pending_to_final(transaction)
                os.rename(p, final_path)
                #print(f"No clash {str(p)=} {str(final_path)=}", file=sys.stderr)
        else:
            MAX_STAGES=99999
            clashing_index = MAX_STAGES
            #print(f"{self.reproducible_steps=}", file=sys.stderr)
            for i,step in enumerate(self.reproducible_steps):
                for dest in step:
                    dp = self.convert_final_to_pending(dest)
                    #print(f"{dp=} {self.pending_steps.keys()=}", file=sys.stderr)
                    if dp in self.pending_steps:
                        if i < clashing_index:
                            clashing_index = i
                            break
            #print(f"{clashing_index=}", file=sys.stderr)
            if clashing_index == MAX_STAGES:  # this should NEVER happen
                self.log(f"Critical: pipeline detects clashing transaction {clashing_index=} for {self.pending_steps=}") 
                return False

            # keep everything prior to the clash, then add on the pending steps
            remove_these_steps = self.reproducible_steps[clashing_index:]
            #print(f"{remove_these_steps=}", file=sys.stderr)
            for step in remove_these_steps:
                for fp in step:
                    if fp in self.uploaded_files:
                        del self.uploaded_files[fp]
                    if Path(fp).exists():
                        #print(f"removing the original file: {str(fp)}", file=sys.stderr)
                        os.remove(fp)
            # rename pending filepaths
            for transaction in self.pending_steps:
                p = Path(transaction)
                if not p.exists():
                    self.log(f'Warning: {str(p)} not found')
                    continue
                final_path = self.convert_pending_to_final(transaction)
                #print(f"Renaming {str(p)=} to {str(final_path)=}", file=sys.stderr)
                os.rename(p, final_path)
                op = str(p)
                if op in self.uploaded_files:
                    self.uploaded_files[final_path] = self.uploaded_file[op].copy()
                    del self.uploaded_file[op]
            self.reproducible_steps = self.reproducible_steps[0:clashing_index]
        if self.pending_steps is None:
            self.save()
            return True
        this_step = {}
        for ps in self.pending_steps:
            if ps is None:
                continue
            final_name = self.convert_pending_to_final(ps)
            record = self.pending_steps[ps]
            this_step[final_name] = record
            self.uploaded_files[final_name] = self.pending_steps[ps].keys()
        self.reproducible_steps.append(this_step)
        self.pending_steps = None
        self.save()
        return True

    def enforce_file_consistency(self):
        """ 
        Iterate over all expected files and ensure that they are present. If they aren't, remove these steps and any that follow.
        """
        pids_to_delete = set()
        files_to_delete = set()
        for i, step in enumerate(self.reproducible_steps):
            for f in step:
                if not Path(f).exists():
                    files_to_delete.add(f)
                    if step[f] is None:
                        continue
                    for pid in step[f]:
                        pids_to_delete.add(pid)
        # clear out files and look for pipeline breakages from missing sections
        pipe_broken_step = None
        for i, step in enumerate(self.reproducible_steps):
            for f in files_to_delete:
                if f in self.uploaded_files:
                    del self.uploaded_files[f]
                if f in step:
                    del step[f]
                if len(step) == 0:
                    pipe_broken_step = i
                    break
            if pipe_broken_step is not None:
                break
        if pipe_broken_step is not None:
            self.reproducible_steps = self.reproducible_steps[:pipe_broken_step]

        # now clean out pids 
        for i, step in enumerate(self.reproducible_steps):
            for f in step:
                for pid in pids_to_delete:
                    if pid in step[f]:
                        del self.reproducible_steps[i][f][pid]


    def get_plate(self, PID, transactions=None):
        """ 
        Look up a plate ID in self.plate_location_sample and then in self.reproducible_steps
        Apply all changes seen in self.reproducible steps
        Return the modified plate
        PID should be a guarded plate ID 
        """
        if PID not in self.plate_location_sample:
            self.log(f'Error: {PID} not found in plate records')
            return None
        mod_plate = deepcopy(self.plate_location_sample[PID])
        # self.reproducible_steps is a list, so stage will be a dictionary
        for stage in self.reproducible_steps:
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

    ### functions for returning locally held file paths

    def get_exp_dir(self):
        """ return the experiment directory name """
        return 'run_' + self.name

    def get_exp_fp(self, filename, transaction=False):
        """ 
        Return the expected experiment path to filename as a string
        if transaction is True, check for an existing match and append "_pending" if required        
        """
        dirname = self.get_exp_dir()
        if dirname not in filename:
            fp = Path(os.path.join(dirname, filename))
        else:
            fp = Path(filename)
        if transaction:
            parent_path = fp.parent
            file_name = 'pending_' + str(fp.name)
            fp = parent_path / file_name                   
        #print(str(fp), file=sys.stderr)
        return str(fp)

    def get_raw_dirpath(self):
        """ 
        Return the absolute path to the where raw fastq files should be stored
        Create it if it doesn't exist
        """
        dp = self.get_exp_fp('raw')
        if not os.path.exists(dp):
            os.mkdir(dp)
        return dp

    def get_clean_dirpath(self):
        """
        Return the absolute path to the where cleaned fastq files should be stored
        Create it if it doesn't exist
        """
        dp = self.get_exp_fp('cleaned')
        if not os.path.exists(dp):
            os.mkdir(dp)
        return dp

    def get_merged_dirpath(self):
        """
        Return the absolute path to the where cleaned, merged fastq files should be stored
        Create it if it doesn't exist
        """
        dp = self.get_exp_fp('merged')
        if not os.path.exists(dp):
            os.mkdir(dp)
        return dp

    def get_raw_fastq_pairs(self):
        """ return a sorted list of tuple(R1_path, R2_path) to raw FASTQ files """
        valid_pairs = []
        rdp = self.get_raw_dirpath()
        r1s = [rdp/Path(f) for f in os.listdir(rdp) if f.endswith('.fastq.gz') and '_R1_' in f]
        for r1 in r1s:
            r2 = Path(str(r1).replace('_R1_','_R2_'))
            if r1.is_file() and r2.is_file():
                valid_pairs.append((r1,r2))
            else:
                if not r1.is_file():
                    self.log(f'Warning: {r1} expected raw FASTQ file does not exist')
                elif not r2.is_file():
                    self.log(f'Warning: {r2} expected raw FASTQ file does not exist')
        return sorted(valid_pairs)

    ### functions for handling plates

    def add_uploaded_file(self, file_name, PIDs=[], purpose=None):
        """
        Add an uploaded file to the experiment with associated plates and purpose
        """
        if file_name in self.uploaded_files:
            self.log(f'Warning: file name {file_name} has already been uploaded, overwriting')

        self.uploaded_files[file_name] = {'plates': [util.guard_pbc(PID, silent=True) for PID in PIDs], 'purpose': purpose}

   
    def add_rodentity_plate_set(self, sample_plate_ids, dna_plate_id):
        """ 
        Add mouse sample information to the experiment for a set of up to 4 ear punch plates from Rodentity JSON files
        Return True on success, False on failure
        """
        #print(f"{sample_plate_ids=}{dna_plate_id=}")
        if self.locked:
            self.log('Error: Cannot add Rodentity plate set while lock is turned on')
            return False
        try:
            self.log('Adding Rodentity plate set', level='b')  # records the function starting
            dna_plate_id = util.guard_pbc(dna_plate_id, silent=True)
            sample_plate_ids = sorted([util.guard_pbc(spid, silent=True) for spid in sample_plate_ids if spid])
        
            # check required data is present
            if len(sample_plate_ids) == 0:
                self.log(f"Error: Sample plate ids not present {sample_plate_ids=}")
                return False
            if not dna_plate_id:  # '' or None
                self.log(f"Error: DNA plate id not present {dna_plate_id=}")
                return False
 
            # get info from JSON files
            sample_info = {}
            for spid in sample_plate_ids:
                fn = os.path.join('run_'+self.name, db_io._eppfn_r(spid))  # rodentity file name format
                if os.path.isfile(fn): # if cached in local file
                    self.log('Info: Reading cached data: ' + fn)
                    with open(fn) as src:
                        info = json.load(src)
                    if not info:
                        self.log(f"Error: No data returned for barcode {spid}")
                        return False
                    sample_info[spid] = info
                    if fn in self.uploaded_files:
                        self.log(f'Warning: {fn} already present in records, overwriting')
                    self.uploaded_files[fn] = {'plates':[util.guard_pbc(spid, silent=True)], 'purpose':'rodentity ear punch'}
                else:
                    self.log("Error: JSON file doesn't exist: " + fn)
                    return False

            # no failures so update the experiment
            self.dest_sample_plates[dna_plate_id] = sample_plate_ids
            self.plate_location_sample[dna_plate_id] = {'purpose':'dna', 'source':','.join(sample_plate_ids), 'wells':set()}
            for spid in sample_plate_ids:
                if spid not in self.plate_location_sample:
                    self.plate_location_sample[spid] = {'purpose':'sample','source':'rodentity', 'wells':set()}
                #print(f"{sample_info=}")
                for record in sample_info[spid]['wells']:  # "wellLocation", "mouse" 
                    pos = util.unpadwell(record['wellLocation'])
                    self.plate_location_sample[spid]['wells'].add(pos)
                    if pos not in self.plate_location_sample[spid]:
                        self.plate_location_sample[spid][pos] = {}
                    self.plate_location_sample[spid][pos]['barcode'] = util.guard_rbc(record['mouse']['barcode'], silent=True)
                    self.plate_location_sample[spid][pos]['strain'] = record['mouse']['mouselineName']
                    self.plate_location_sample[spid][pos]['sex'] = record['mouse']['sex']
                    self.plate_location_sample[spid][pos]['mouse'] = deepcopy(record['mouse'])
                    # "alleles":[{alleleKey, name, symbol, options:[], assays:[{assayKey, name, method}]}]
                    self.plate_location_sample[spid][pos]['assay_records'] = {}
                    assays = []
                    assayFamilies = set()
                    unknown_assays = []
                    unknown_assayFamilies = set()
                    
                    if 'alleles' in record['mouse']:
                        for allele in record['mouse']['alleles']:  # allele is a dict from a list
                            for assay in allele['assays']:  # assay is a dict from a list
                                if assay['name'] not in self.plate_location_sample[spid][pos]['assay_records']:
                                    self.plate_location_sample[spid][pos]['assay_records'][assay['name']] =\
                                            {'assayFamily':assay['name'].split('_')[0]}
                                self.plate_location_sample[spid][pos]['assay_records'][assay['name']]['alleleKey'] = str(allele['alleleKey'])
                                self.plate_location_sample[spid][pos]['assay_records'][assay['name']]['alleleSymbol'] = str(allele['symbol'])
                                self.plate_location_sample[spid][pos]['assay_records'][assay['name']]['assayKey'] = str(assay['assay_key'])
                                self.plate_location_sample[spid][pos]['assay_records'][assay['name']]['assayName'] = str(assay['name'])
                                self.plate_location_sample[spid][pos]['assay_records'][assay['name']]['assayMethod'] = str(assay['method'])
                                if assay['method'] == 'NGS' or assay['name'].startswith('NGS'):
                                    assays.append(assay['name'])
                                    assayFamilies.add(assay['name'].split('_')[0])
                                else:
                                    unknown_assays.append(assay['name'])
                                    unknown_assayFamilies.add(assay['name'].split('_')[0])
                            
                    self.plate_location_sample[spid][pos]['assays'] = assays.copy()
                    self.plate_location_sample[spid][pos]['assayFamilies'] = list(assayFamilies)
                    self.plate_location_sample[spid][pos]['unknown_assays'] = unknown_assays.copy()
                    self.plate_location_sample[spid][pos]['unknown_assayFamilies'] = list(unknown_assayFamilies) 
                
                self.log(f"Success: added sample plate {spid} with destination {dna_plate_id} ")

            #print(f"{self.name=} {self.plate_location_sample=} {self.sample_plates=}")
        except Exception as exc:
            print(f"Failed to load Rodentity plate set {exc=}")
            if dna_plate_id in self.dest_sample_plates:
                self.dest_sample_plates.pop(dna_plate_id)
            return False
         
        finally:
            self.save()
        return True


    def add_musterer_plate_set(self, sample_plate_ids, dna_plate_id):
        """ Add mouse sample information to the experiment for a set of up to 4 ear punch plates from Musterer JSON files
            Checks validity of inputs - duplicate epps is a warning, duplicate dnap is disallowed
            
            Input can come from either DB or JSON files
            DEPRECATED - keep as template for reading from DB
        """ 
        if self.locked:
            self.log('Error: Cannot add Musterer plate set while lock is turned on')
            return False
        try:
            self.log(f"Begin: Adding Musterer plate set")

            sample_plate_ids = sorted([spid for spid in sample_plate_ids if spid])
        
            # check required data is present
            if not dna_plate_id:  # '' or None
                self.log(f"Error: DNA plate id not present {dna_plate_id=}")
                return False
            # check if duplicates exist
            for subset in combinations(sample_plate_ids, 2):
                if subset[0] == subset[1] and subset[0] != '':
                    self.log(f"Error: Ear punch plate barcode {subset[0]} is duplicated in inputs {sample_plate_ids}")
                    return False
            # check if DNA plate id already used
            if dna_plate_id in self.plate_location_sample:
                self.log(f"Error: DNA plate {dna_plate_id} already used as an output plate id")
                return False
            # check that none of the sample plate ids match the DNA plate id, or any other DNA plate ids
            for epp in sample_plate_ids:
                if epp == dna_plate_id:
                    self.log(f"Error: Sample plate barcode matches DNA plate barcode {dna_plate_id}")
                    return False
                if epp in self.plate_location_sample:
                    self.log(f"Error: Sample plate {epp} already used as a DNA plate barocde")
                    return False

            # get info from JSON files
            sample_info = {}
            for spid in sample_plate_ids:
                fn = os.path.join('run_'+self.name, db_io._eppfn_m(spid))
                if os.path.isfile(fn): # if cached in local file
                    self.log('Info: Reading cached data: ' + fn)
                    with open(fn) as src:
                        info = json.load(src)
                    if not info:
                        self.log(f"Error: No data returned for barcode {spid}")
                        return False
                    sample_info[spid] = info['wells']  # list of well contents
                    if fn in self.uploaded_files:
                        self.log(f'Warning: {fn} already present in records, overwriting')
                    self.uploaded_files[fn] = {'plates':[util.guard_pbc(spid, silent=True)], 'purpose':'musterer ear punch'}
                else:
                    self.log("Error: JSON file doesn't exist: " + fn)
                    return False

            # no failures so update the experiment
            self.dest_sample_plates[dna_plate_id] = sorted(sample_plate_ids)
            for spid in sample_plate_ids:
                if spid not in self.plate_location_sample:  # should always be the case
                    self.plate_location_sample[spid] = {'purpose':'sample','source':'musterer','wells':set()}
                for i, record in enumerate(sample_info):  # 'wellLocation', 'mouse', 'mouseBarcode', 'mouseId'
                    pos = util.unpadwell(record['wellLocation'])
                    self.plate_location_sample[spid]['wells'].add(pos)
                    self.plate_location_sample[spid][pos] = record['mouse'].copy()  # dict
                    self.plate_location_sample[spid][pos]['mouseId'] = record['mouseId']
                    self.plate_location_sample[spid][pos]['barcode'] = util.guard_mbc(record['mouseBarcode'])
                    self.plate_location_sample[spid][pos]['must_assays'] = record['mouse']['assays'].copy()
                    self.plate_location_sample[spid][pos]['sampleNumber'] = i+1
                    assays = [assay['assayName'] for assay in record['mouse']['assays']]
                    assayFamilies = set([a.split('_')[0] for a in assays])
                    self.plate_location_sample[spid][pos]['assays'] = assays.copy()
                    self.plate_location_sample[spid][pos]['assayFamilies'] = list(assayFamilies)
                
                self.log(f"Success: added sample plate {spid} with destination {dna_plate_id} ")

            #print(f"{self.name=} {self.plate_location_sample=} {self.sample_plates=}")  
        finally:
            self.save()
        return True

    #def add_manifest(self, manifest_stream, default_manifest_type='c'):
    #    """ Because streamlit's uploader is derived from BytesIO we have a byte stream we have to deal with.
    #        The manifest is a CSV. Turn this into a list of dicts which can then handle sensibly.
    #        We also take an optional default_manifest_type ['c','r','m'] for custom, rodentity, or musterer respectively
    #        Returns True on success, False on failure
    #       DEPRECATED
    #    """
    #    if self.locked:
    #        self.log('Error: Cannot add manifest while lock is active')
    #        return False
    #    try:
    #        self.log(f"Begin: add custom manifest")
        
    #        #self.log(f"Debug: {self.name=} {manifest_strm=} {default_manifest_type=}")
    #        manifest_name = manifest_stream.name
    #        if manifest_name in self.uploaded_files:
    #            self.log(f'Warning: {manifest_name} already present in records, overwriting')
    #        self.uploaded_files[manifest_name] = {'plates':[], 'purpose':'custom manifest'}
    #        manifest_io = StringIO(manifest_stream.getvalue().decode("utf-8"))
        
    #        header = ''
    #        manifest_info = []
    #        in_csv = csv.reader(manifest_io, delimiter=',')
    #        for i, row in enumerate(in_csv):
    #            if 'barcode' in ''.join(row).lower():
    #                header = row
    #                continue  # header
    #            elif all([col.strip() for col in row]) == '': # blank
    #                continue
    #            elif any(col.strip() == '' and j < 5 for j,col in enumerate(row)):
    #                self.log(f"Warning: not adding blank manifest entry in row {i+1}: {row}")
    #                continue
    #            else:
    #                manifest_info.append({header[i]:r.strip() for i,r in enumerate(row) if i<len(header)})

    #        # gather up dest plates and sample plates so we can clear these if we already have them in the experiment
    #        dpids_spids = {}  # {dpid:set([spid,..])}
    #        # go through entries and guard barcodes as needed
    #        for i,entry in enumerate(manifest_info):
    #            for key in entry:
    #                if key == 'Dest barcode':
    #                    manifest_info[i][key] = util.guard_pbc(entry[key], silent=True)
    #                    dpids_spids[manifest_info[i][key]] = set()
    #                if key == 'Plate barcode':
    #                    manifest_info[i][key] = util.guard_pbc(entry[key], silent=True)
    #                if key == 'Sample barcode':
    #                    if not util.is_guarded(entry[key]):
    #                        if default_manifest_type == 'c':
    #                            manifest_info[i][key] = util.guard_cbc(entry[key])
    #                        elif default_manifest_type == 'r':
    #                            manifest_info[i][key] = util.guard_rbc(entry[key])
    #                        elif default_manifest_type == 'm':
    #                            manifest_info[i][key] = util.guard_mbc(entry[key])
    #        for i,entry in enumerate(manifest_info):
    #            dpids_spids[entry['Dest barcode']].add(entry['Plate barcode'])

    #        # now clear any existing plateIds
    #        for dpid in dpids_spids:
    #            if dpid in self.plate_location_sample:
    #                self.log(f"Info: Clearing existing plate details for {dpid=}")
    #                self.plate_location_sample.pop(dpid)
    #            for spid in dpids_spids[dpid]:
    #                if spid in self.plate_location_sample:
    #                    self.log(f"Info: Clearing existing plate details for {spid=}")
    #                    self.plate_location_sample.pop(spid)

    #        #print(f"{manifest_info=}")
    #        dest_pids_sample_pids = {}
    #        sample_pids = set()
    #        for entry in manifest_info:
    #            if entry['Dest barcode'] not in dest_pids_sample_pids:
    #                dest_pids_sample_pids[entry['Dest barcode']] = set()
    #            dest_pids_sample_pids[entry['Dest barcode']].add(entry['Plate barcode'])
    #            sample_pids.add(entry['Plate barcode'])
        
    #        #self.dest_sample_plates = {}  # {dest_pid:[4 sample plate ids]}
    #        #self.plate_location_sample = {}  # pid:{well:sample_id}
    #        #self.sample_info = {}  # sample_id: {strain:str, assays: [], possible_gts: [], 
    #        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] }

    #        # check for overused destination PIDs and clashes with sample plate barcodes. Note: we can't meaningfully check for duplicated destination PIDs
    #        dp_count = {dp:len(sps) for dp, sps in dest_pids_sample_pids.items()}     
    #        for dp in dp_count:
    #            if dp_count[dp] > 4:
    #                self.log(f"Error: 384-well DNA destination plate {dp} used by too many sample plates: {dp_count[dp]}")
    #                return False
    #            if dp in sample_pids:
    #                self.log(f"Error: 384-well DNA destination plate barcode also used as a sample plate barcode {dp}")
    #                return False
            

    #        # Warn for duplicate sample plate barcodes?
    #        for sp in sample_pids:
    #            dp_set = set()
    #            for dp in dest_pids_sample_pids:
    #                for entry in manifest_info:
    #                    if entry['Plate barcode'] == sp and entry['Dest barcode'] == dp:
    #                        dp_set.add(dp)
    #            if len(dp_set) > 1:
    #                self.log(f"Warning: Sample plate barcode {sp} used in multiple 384-well DNA destination plates: {dp_set}")
        
    #        # we need to collect destination, sample lists
    #        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], gts: [],
    #        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
    #        dp_samples = {}
    #        for i,entry in enumerate(manifest_info):
    #            if 'Sample No' in entry:
    #                sample_number = entry['Sample No']
    #            else:
    #                sample_number = i+1
    #            dest_pid = entry['Dest barcode']
    #            if dest_pid not in self.uploaded_files[manifest_name]['plates']:
    #                self.uploaded_files[manifest_name]['plates'].append(dest_pid)
    #            source_pid = entry['Plate barcode']
    #            if source_pid not in self.uploaded_files[manifest_name]['plates']:
    #                self.uploaded_files[manifest_name]['plates'].append(source_pid)
    #            well = util.unpadwell(entry['Well'])
    #            assays = [entry[key] for key in entry if 'assay' in key.lower() and entry[key].strip()!='']
    #            assayFamilies = set([a.split('_')[0] for a in assays])
                
    #            if dest_pid not in dp_samples:
    #                dp_samples[dest_pid] = set()
    #            dp_samples[dest_pid].add(source_pid)
    #            if source_pid not in self.plate_location_sample:
    #                self.plate_location_sample[source_pid] = {'purpose':'sample','source':'manifest', 'wells':set()}
    #            if well in self.plate_location_sample[source_pid] and self.plate_location_sample[source_pid][well] != {}:
    #                self.log(f"Error: duplicate {well=} in {source_pid=}. Continuing...")

    #            self.plate_location_sample[source_pid]['wells'].add(well)
    #            self.plate_location_sample[source_pid][well] = {}
    #            self.plate_location_sample[source_pid][well]['barcode'] = entry['Sample barcode']
    #            self.plate_location_sample[source_pid][well]['assays'] = assays.copy()
    #            self.plate_location_sample[source_pid][well]['assayFamilies'] = list(assayFamilies)
    #            self.plate_location_sample[source_pid][well]['sex'] = ''
    #            self.plate_location_sample[source_pid][well]['strain'] = ''
    #            self.plate_location_sample[source_pid][well]['other_id'] = ''
    #            self.plate_location_sample[source_pid][well]['sampleNumber'] = str(sample_number)
    #            for key in entry:
    #                if 'sex' in key.lower():
    #                    self.plate_location_sample[source_pid][well]['sex'] = entry[key]
    #                elif 'strain' in key.lower():
    #                    self.plate_location_sample[source_pid][well]['strain'] = entry[key]
    #                elif 'other_id' in key.lower():
    #                    self.plate_location_sample[source_pid][well]['other_id'] = entry[key]
    #            self.plate_location_sample[source_pid][well]['gts'] = []
    #            self.plate_location_sample[source_pid][well]['possible_gts'] = []
    #            self.plate_location_sample[source_pid][well]['parents'] = []      

    #        #print(f"\n{dp_samples=}")
    #        for dp in dp_samples:
    #            self.dest_sample_plates[dp] = list(dp_samples[dp])
    #            for sample_pid in dp_samples[dp]:
    #                self.log(f"Success: added sample plate {sample_pid} with destination {dp}")
    #                #print('')
    #                #print(f"{sample_pid=} {self.plate_location_sample[sample_pid]=}")
    #    finally:
    #        self.save()
    #    return True

    def read_custom_manifests(self, manifests):
        """
        Parse any number of custom manifest files and store them in self.unassigned_plates['custom'] =\
                {plateBarcode={Assay=[],...}}
        Also add an entry for each file into self.uploaded_files
        Example headers:
        sampleNum	plateBarcode	well	sampleBarcode	assay	assay	clientName	sampleName 	alleleSymbol
        Required columns: [plateBarcode, well, sampleBarcode, assay*, clientName] *Duplicates allowed
        Optional columns: [sampleNo, sampleName, alleleSymbol] and anything else
        """
        if self.locked:
            self.log('Error: Cannot add manifest while lock is active')
            return False
        #try:
        if True:
            self.log(f"Begin: read custom manifest")
            file_names, file_tables = file_io.read_csv_or_excel_from_stream(manifests)
            for file_name, file_table in zip(file_names, file_tables):
                if file_name in self.uploaded_files:
                    self.log(f'Warning: {file_name} already present in records')
                #manifest_name = manifest_stream.name
                #if manifest_name in self.uploaded_files:
                #    self.log(f'Warning: {manifest_name} already present in records')
                ##print(f'{manifest_name=}', file=sys.stderr)
                #if manifest_name.lower().endswith('xlsx'):
                #    workbook = openpyxl.load_workbook(BytesIO(manifest_stream.getvalue()))
                #    sheet = workbook.active
                #    rows = [','.join(map(str,cells)) for cells in sheet.iter_rows(values_only=True)]
                #    #print(rows, file=sys.stderr)
                #else:
                #    rows = StringIO(manifest_stream.getvalue().decode("utf-8"))

                plate_entries = {}  # looks like self.plate_location_sample, but temporary
                plate_barcode_col = None
                assay_cols = []
                for i, line in enumerate(file_table):
                    cols = [c.strip() if c is not None else '' for c in line.split(',')]  # lower case column names
                    #print(i, cols, file=sys.stderr)
                    if i==0:  # process header
                        cols_lower = [c.lower() for c in cols]
                        header_dict = {k:c for k,c in enumerate(cols_lower)}
                        #print(header_dict, file=sys.stderr)
                        matching_cols = [col_name in cols_lower for col_name in ['platebarcode', 'well', 'samplebarcode', 'assay', 'clientname']]
                        #print(matching_cols, file=sys.stderr)
                        if not all(matching_cols):
                            self.log(f'Error: manifest {file_name} requires at least columns plateBarcode, well, sampleBarcode, assay, clientName')
                            return False
                        else:
                            self.log(f'Info: parsing manifest {file_name}')
                        for k in header_dict:
                            if header_dict[k] == 'platebarcode':
                                plate_barcode_col = k
                            elif header_dict[k] == 'assay':
                                assay_cols.append(k)
                        continue
                    if len(cols) < 5:  # skip empty rows
                        continue
                    gpid = util.guard_pbc(cols[plate_barcode_col], silent=True)
                    if gpid not in plate_entries:
                        plate_entries[gpid] = {'purpose':'sample','source':'manifest', 'wells':set()}  # create plate_location_sample entries here
                    assays = [] 
                    # do a first pass to set up recording a sample in a well
                    for k,c in enumerate(cols):
                        if c.lower() == 'none':
                            continue
                        if k in assay_cols:
                            if c != '':  # ignore empty assay entries
                                assays.append(c)
                        if header_dict[k] == 'well':
                            well = util.unpadwell(c.upper())
                            if well in plate_entries[gpid]['wells'] and plate_entries[gpid][well] != {}:
                                self.log(f"Error: duplicate {c} in {gpid}. Skipping row {i+2} {cols=}")
                                break
                            plate_entries[gpid]['wells'].add(well)
                            plate_entries[gpid][well] = {}
                    # now do a second pass to collect everything together
                    for k,c in enumerate(cols):
                        if c.lower() == 'none':
                            c = ''
                        if header_dict[k] == 'well':
                            continue
                        if header_dict[k] == 'samplebarcode':
                            #if c.startswith('C'):
                            #    sid = util.guard_cbc(c, silent=True)
                            #elif c.startswith('M'):
                            #    sid = util.guard_rbc(c, silent=True)
                            #else:
                            sid = util.guard_cbc(c, silent=True) # fall back to custom?
                            plate_entries[gpid][well]['barcode'] = sid
                        elif header_dict[k] == 'platebarcode':
                            plate_entries[gpid][well]['platebarcode'] = gpid
                        elif header_dict[k] == 'sampleno':
                            plate_entries[gpid][well]['sampleNumber'] = str(c)
                        elif header_dict[k] == 'assay':
                            continue
                        else:
                            try:
                                plate_entries[gpid][well][header_dict[k]] = str(c)
                            except:
                                print(type(header_dict[k]), header_dict[k], str(c), file=sys.stderr)
                    plate_entries[gpid][well]['assays'] = assays
                    plate_entries[gpid][well]['assayFamilies'] = list(set([a.split('_')[0] for a in assays]))

                for gpid in plate_entries:
                    if gpid in self.unassigned_plates or gpid in self.plate_location_sample:
                        self.log(f'Warning: plate records exist for {util.unguard_pbc(gpid, silent=True)}, potentially overwriting')
                    self.unassigned_plates['custom'][gpid] = plate_entries[gpid]
                self.uploaded_files[file_name] = {'plates':plate_entries.keys(), 'purpose':'custom manifest'}
        return True
                    

    def accept_custom_manifests(self, dest_pid, sample_pids):
        """
        Move four custom plate records from self.unassigned_plates to self.plate_location_sample and self.dest_sample_plates
        """
        if self.locked:
            self.log('Error: Cannot add manifest while lock is active')
            return False
        #try:
        if True:
            self.log(f"Begin: accept custom manifest")
            # Do checks
            dest_pid = util.guard_pbc(dest_pid, silent=True)
            sample_pids = [util.guard_pbc(pid, silent=True) for pid in sample_pids if pid != 'None' and pid is not None]
            for pid in [dest_pid] + sample_pids:
                if pid in self.dest_sample_plates or pid in self.plate_location_sample:
                        self.log(f'Warning: Plate barcode {util.unguard_pbc(pid, silent=True)} already in use, overwriting')
            for pid in sample_pids:
                if pid not in self.unassigned_plates['custom']:
                    self.log(f'Critical: {util.unguard_pbc(pid, silent=True)} not found in {self.unassigned_plates=}')
                    return False
            # do assignment
            for pid in sample_pids:
                self.plate_location_sample[pid] = self.unassigned_plates['custom'][pid]
            self.dest_sample_plates[dest_pid] = sample_pids
            self.plate_location_sample[dest_pid] = {'purpose':'dna', 'source':','.join(sample_pids), 'wells':set()}
            self.unassigned_plates['custom'] = {'None':{}}
        return True


    def get_musterer_pids(self):
        """ return [pids] with samples that are sourced from Musterer. DEPRECATED """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        musterer_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if util.is_guarded_mbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    musterer_pids.add(pid)
        return musterer_pids


    def get_rodentity_pids(self):
        """ return [pids] with samples that are sourced from Rodentity """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        rodentity_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if util.is_guarded_rbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    rodentity_pids.add(pid)
        return rodentity_pids


    def get_custom_pids(self):
        """ return [pids] with samples that are sourced as custom - useful for mixed content plates """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        custom_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if util.is_guarded_cbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    custom_pids.add(pid)
        return custom_pids
        

    def summarise_inputs(self):
        """ return a list of fields and a list of headers, summarising the contents of plate sets """
        plate_set_summary = []
        plate_set_headers = ['DNA PID', 'Sample PID1', 'Sample PID2', 'Sample PID3', 'Sample PID4', 'Custom samples', 'Rodentity samples']
        total_wells = 0
        total_unique_samples = set()
        total_unique_assays = set()
        total_well_counts = {'c':0,'m':0,'r':0}
        for dna_pid in self.dest_sample_plates:
            plate_set_details = [util.unguard_pbc(dna_pid)]
            custom_wells = 0
            rodentity_wells = 0
            
            for i,sample_pid in enumerate(sorted(self.dest_sample_plates[dna_pid])):
                plate_set_details.append(util.unguard_pbc(sample_pid))
                for info in (self.plate_location_sample[sample_pid][well] for well in self.plate_location_sample[sample_pid]['wells']):
                    samp_barcode = info['barcode']     
                    if util.is_guarded_cbc(samp_barcode):
                        custom_wells += 1
                        total_well_counts['c'] += 1
                    #elif util.is_guarded_mbc(samp_barcode):
                    #    d[' counts']['m'] += 1
                    #    total_well_counts['m'] += 1  
                    elif util.is_guarded_rbc(samp_barcode):
                        rodentity_wells += 1
                        total_well_counts['r'] += 1   
                    #for assay in info['assays']:
                    #    d['Primers required'].add(assay)
                    #    total_unique_assays.add(assay)
                    
                    #d['Unique samples'].add(info['barcode'])
                    total_unique_samples.add(info['barcode'])
            for j in range(3-i):
                plate_set_details.append('')

            plate_set_details.append(str(custom_wells))
            plate_set_details.append(str(rodentity_wells))
            # print(f'{plate_set_details=}', file=sys.stderr)
            plate_set_summary.append(plate_set_details)
        if len(plate_set_summary) > 0:
            plate_set_summary.append(['Total','','','','',total_well_counts['c'], total_well_counts['r']])
                    # 'Total assays': len(total_unique_assays)})
        #print(f"{plate_set_summary=}")
        return plate_set_summary, plate_set_headers


    def inputs_as_dataframe(self):
        """ return the experiment contents as a pandas dataframe """
        rows, headers = self.summarise_inputs()
        if not rows:
            return None
        return pd.DataFrame(rows, columns=headers)


    def summarise_consumables(self):
        """
        Return a dictionary of all consumables: primers, indexes, taq+water plates
        The descriptions we provide here will likely vary as we find new things to display
        TODO: separate custom and NGS primers
        """
        d = {'taqwater_pids':[], 'taq_vol':0, 'water_vol':0, 'primer_pids':[], 'primer_count_ngs':0, 
                'primer_count_custom':0, 'unique_primers':set(), 'primer_well_count':0, 'assay_primer_mappings':0,
                'reference_files':[], 'unique_references':set(), 'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
        consumable_plate_purposes = set(['primer', 'taq_water', 'index'])
        for pid in self.plate_location_sample:
            if 'purpose' not in self.plate_location_sample[pid]:
                continue
            if self.plate_location_sample[pid]['purpose'] not in consumable_plate_purposes:
                continue
            plate = self.get_plate(pid)  # get the plate contents with all usage modifications applied
            if plate['purpose'] == 'taq_water':
                d['taqwater_pids'].append(pid)
                for well in util.TAQ_WELLS:
                    if well in plate and 'volume' in plate[well]:
                        d['taq_vol'] += plate[well]['volume']
                for well in util.WATER_WELLS:
                    if well in plate and 'volume' in plate[well]:
                        d['water_vol'] += plate[well]['volume']
            elif plate['purpose'] == 'primer':
                d['primer_pids'].append(pid)
                for well in plate['wells']:
                    if 'primer' in plate[well]:
                        d['unique_primers'].add(plate[well]['primer'])
                        d['primer_well_count'] += 1
            elif plate['purpose'] == 'index':
                d['index_pids'].append(pid)
                for well in plate['wells']:
                    if 'idt_name' in plate[well]:
                        if '_i7F_' in plate[well]['idt_name']:
                            d['unique_i7s'].add(plate[well]['idt_name'])
                        elif '_i5R_' in plate[well]['idt_name']:
                            d['unique_i5s'].add(plate[well]['idt_name'])
        for f in self.reference_sequences:
            d['reference_files'].append(f)
            for refname in self.reference_sequences[f]:
                d['unique_references'].add(refname)
        d['assay_primer_mappings'] = len(self.primer_assay)
        return d


    def add_nimbus_outputs(self, nim_outputs):
        """
        Copy Nimbus output files into the project folder
        """
        transactions={}
        #print(f'{type(nim_outputs)=} {len(nim_outputs)=}', file=sys.stderr)
        try:
            print('Starting add_nimbus_outputs', file=sys.stderr)
            for nim_output in nim_outputs:
                #print(f'{nim_output.name=}', file=sys.stderr) 
                fp = self.get_exp_fp(nim_output.name, transaction=True)
                self.log(f"Info: copying {fp} to experiment folder")
                plate_set = set()
                with open(fp, 'wt') as outf:
                    #print(nim_output.getvalue().decode("utf-8"))
                    nim_outstr = nim_output.getvalue().decode("utf-8").replace('\r\n','\n')
                    outf.write(nim_outstr)
                    for i, line in enumerate(nim_outstr.split('\n')):
                        if i == 0:
                            continue
                        cols = line.split('\t')
                        #print(cols, file=sys.stderr)
                        if len(cols) != 7:
                            continue
                        # RecordId	TRackBC	TLabwareId	TPositionId	SRackBC	SLabwareId	SPositionId
                        # 1	p32542p	Echo_384_COC_0001	A1	p83115p	ABg_96_PCR_NoSkirt_0001	A1
                        plate_set.add(util.guard_pbc(cols[1], silent=True))
                        plate_set.add(util.guard_pbc(cols[4], silent=True))
                transactions[fp] = {pid:{} for pid in plate_set}
                final_fp = self.convert_pending_to_final(fp)
                if final_fp in self.uploaded_files:
                    self.log(f'Warning: {final_fp} already recorded as uploaded, overwriting')
                self.uploaded_files[final_fp] = {'plates': list(plate_set), 'purpose':'DNA plate'} 
        except Exception as exc:
            self.log(f'Error: could not upload Hamilton Nimbus output files, {exc}')
            return False
        self.add_pending_transactions(transactions)
        return True

    def add_pcr_plates(self, pcr_plate_list=[]):
        """
        Add one or more empty 384-well PCR plates to self.plate_location_sample
        """
        for pid in pcr_plate_list:
            p = util.guard_pbc(pid, silent=True)
            if p in self.plate_location_sample:
                self.log(f'Error: {p} already in use, skipping')
                continue
            self.plate_location_sample[p] = {'purpose':'pcr', 'wells':set(), 'source':'user', 
                    'plate_type':util.PLATE_TYPES['PCR384'], 'barcode':p}
        self.save()
        return True

                                                                                            
    def get_pcr_pids(self):
        """ return a list of user supplied PCR plate ids """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'pcr']


    def get_assay_usage(self, dna_plate_list=[], filtered=True, included_guards=util.GUARD_TYPES):
        """ We will want ways to get assay usage from various subsets, but for now do it for everything that 
            has a Nimbus destination plate.
        """
        assay_usage = {}
        for dest in self.dest_sample_plates:
            if dna_plate_list and dest not in dna_plate_list:
                continue
            for sample_pid in self.dest_sample_plates[dest]:
                plate = self.get_plate(sample_pid)
                assay_usage.update(util.calc_plate_assay_usage(plate,filtered=filtered, included_guards=included_guards))
        return assay_usage


    def get_volumes_required(self, assay_usage=None, dna_plate_list=[], filtered=True):
        """
            Standard usage is to call get_assay_usage() first and pass it to this.
            if filtered, only apply those assays allowed by the assay/primer list.
        """
        if not assay_usage:
            assay_usage = self.get_assay_usage(dna_plate_list=dna_plate_list, filtered=filtered)  # all assays present
            #assay_list = set([self.primer_assay[prmr] for prmr in self.primer_assay])
            #assay_usage = {a:assay_usage[a] for a in assay_usage if a in assay_list}
        #print (f'{assay_usage=}', assay_usage.values(), file=sys.stderr)
        reactions = sum([v for v in assay_usage.values()])

        # convert assay name to primer name
        #for a in assay_usage:
        #    assay_primer_conversion_list = set([self.primer_assay[prmr] for prmr in self.primer_assay])
        #    if a in assay_primer_conversion_list:
        #        assay_usage(self.assay

        primer_taq_vol = reactions * self.transfer_volumes['PRIMER_TAQ_VOL']
        primer_water_vol = reactions * self.transfer_volumes['PRIMER_WATER_VOL']
        index_water_vol = reactions * self.transfer_volumes['INDEX_WATER_VOL']
        index_taq_vol = reactions * self.transfer_volumes['INDEX_TAQ_VOL']
        primer_vols = {a:assay_usage[a]*self.transfer_volumes['PRIMER_VOL'] for a in assay_usage}
        return primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol
    

    def get_index_avail(self, included_pids=None):
        """
        Returns {primer:count}, {primer:vol} from what's been loaded 
        If included_pids is an iterable, only include plates with these ids
        TODO: check the validity of line 1768ish fwd_idx and rev_idx. Assignment/comparison...
        """
        index_pids = []
        warning_idxs = ''
        for pid in self.plate_location_sample:
            if self.plate_location_sample[pid]['purpose'] == 'index':
                index_pids.append(pid)

        if included_pids:
            guarded_included_pids = [util.guard_pbc(pid, silent=True) for pid in included_pids]
            index_pids = [pid for pid in index_pids if pid in guarded_included_pids]
                
        fwd_idx = {}
        rev_idx = {}

        for idx_pid in index_pids:
            idx_plate = self.get_plate(idx_pid)

            #print(f'get_index_remaining_available_volume() {idx_plate=}', file=sys.stderr)
            
            for well in idx_plate['wells']:
                if 'idt_name' not in idx_plate[well]:
                    continue
                name = idx_plate[well]['idt_name']
                # if 'volume' in idx_plate[well]:
                #     #print(f"get_index_remaining_available_volume() {idx_plate[well]['volume']=}")
                if 'i7F' in name:
                    if name not in fwd_idx:
                        fwd_idx[name] = {'count':0, 'req_vol':[], 'avail_vol':[]}
                    fwd_idx[name]['count'] += 1
                    if 'volume' in idx_plate[well]:
                        fwd_idx[name]['avail_vol'].append(max(idx_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0))
                elif 'i5R' in name:
                    if name not in rev_idx:
                        rev_idx[name] = {'count':0, 'req_vol':[],'avail_vol':[]}
                    rev_idx[name]['count'] += 1
                    if 'volume' in idx_plate[well]:
                        rev_idx[name]['avail_vol'].append(max(idx_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0))
                else:
                    self.log('Unexpected index name:' + name, level='Warning')

        max_i7F = len(fwd_idx)
        max_i5R = len(rev_idx)
        for idx in fwd_idx.keys():
            fwd_idx[idx]['req_vol'].append(max_i5R*self.transfer_volumes['INDEX_VOL']/1000)
            if 'avail_vol' not in fwd_idx[idx] or max_i5R*self.transfer_volumes['INDEX_VOL']/1000 > sum(fwd_idx[idx]['avail_vol']):
                warning_idxs += idx + ', '
        for idx in rev_idx.keys():
            rev_idx[idx]['req_vol'].append(max_i7F*self.transfer_volumes['INDEX_VOL']/1000)
            if 'avail_vol' not in rev_idx[idx] or max_i7F*self.transfer_volumes['INDEX_VOL']/1000 > sum(rev_idx[idx]['avail_vol']):
                warning_idxs += idx + ', '

        return fwd_idx, rev_idx, warning_idxs
     
    
    def get_index_remaining_available_volume(self, assay_usage=None, fwd_idx=None, rev_idx=None):
        """
        Returns the barcode pairs remaining, max available barcode pairs, max barcode pairs allowed by volume
        """
        if not assay_usage:
            assay_usage = self.get_assay_usage()
        reactions = sum([v for v in assay_usage.values()])

        if not fwd_idx or not rev_idx:
            return 0, 0, False

        max_i7F = len(fwd_idx)
        max_i5R = len(rev_idx)

        fwd_idx_reactions = {}
        reaction_vol_capacity = 0
        for name in fwd_idx.keys():
            # get the number of possible reactions and keep the lower number from possible reactions or possible reaction partners
            max_reactions = sum([floor((vol-util.DEAD_VOLS[util.PLATE_TYPES['Echo384']])/self.transfer_volumes['INDEX_VOL']) \
                    for vol in fwd_idx[name]['avail_vol']])
            reaction_vol_capacity += min(max_reactions, max_i5R)

        max_idx_pairs = max_i7F * max_i5R
        
        return max_idx_pairs-reactions, max_idx_pairs, reaction_vol_capacity

    
        
    # def get_index_remaining_available_volume(self, assay_usage=None):
    #    """
    #    Returns the barcode pairs remaining, max available barcode pairs, max barcode pairs allowed by volume
    #    """
    #    index_pids = []
    #    for pid in self.plate_location_sample:
    #        if self.plate_location_sample[pid]['purpose'] == 'index':
    #            index_pids.append(pid)

    #    #print(f'get_index_remaining_available_volume() {index_pids=}', file=sys.stderr)
    #    if not index_pids:
    #        return 0, 0, False

    #    if not assay_usage:
    #        assay_usage = self.get_assay_usage()
    #    reactions = sum([v for v in assay_usage.values()])
        
    #    fwd_barcode_vols = {}  # name=[vol, vol, ...]
    #    rev_barcode_vols = {}
        
    #    for idx_pid in index_pids:
    #        idx_plate = self.get_plate(idx_pid)
    #        #print(f'get_index_remaining_available_volume() {idx_plate=}', file=sys.stderr)
            
    #        for well in idx_plate['wells']:
    #            if 'idt_name' not in idx_plate[well]:
    #                continue
    #            name = idx_plate[well]['idt_name']
    #            if 'volume' not in idx_plate[well]:
    #                continue
    #                #print(f"get_index_remaining_available_volume() {idx_plate[well]['volume']=}")
    #            if 'i7F' in name:
    #                if name not in fwd_barcode_vols:
    #                    fwd_barcode_vols[name] = []
    #                    fwd_barcode_vols[name].append(max(idx_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0))
    #            elif 'i5R' in name:
    #                if name not in rev_barcode_vols:
    #                    rev_barcode_vols[name] = []
    #                    rev_barcode_vols[name].append(max(idx_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0))
    #            else:
    #                self.log('Unexpected index name:' + name, level='Warning')
    #    max_i7F = len(fwd_barcode_vols)
    #    max_i5R = len(rev_barcode_vols)
    #    reaction_vol_capacity = 0
    #    for name in fwd_barcode_vols:
    #        # get the number of possible reactions and keep the lower number from possible reactions or possible reaction partners
    #        max_reactions = sum([floor((vol-util.DEAD_VOLS[util.PLATE_TYPES['Echo384']])/self.transfer_volumes['INDEX_VOL']) \
    #                for vol in fwd_barcode_vols[name]])
    #        reaction_vol_capacity += min(max_reactions, max_i5R)
    #    max_barcode_pairs = max_i7F * max_i5R

    #    return max_barcode_pairs-reactions, max_barcode_pairs, reaction_vol_capacity


    def get_primers_avail(self, included_pids=None):
        """ 
        Returns {primer:count}, {primer:vol} from what's been loaded 
        If included_pids is an iterable, only include plates with these ids
        """
        primer_counts = {}
        primer_vols = {}
        for pid in self.plate_location_sample:
            if self.plate_location_sample[pid]['purpose'] == 'primer':
                if included_pids:
                    if pid not in included_pids:
                        continue
                pmr_plate = self.get_plate(pid)
                for well in pmr_plate['wells']:
                    if 'primer' in pmr_plate[well]:
                        primer_name = pmr_plate[well]['primer']
                        if primer_name not in primer_counts:
                            primer_counts[primer_name] = 0
                        if primer_name not in primer_vols:
                            primer_vols[primer_name] = 0
                        primer_counts[primer_name] += 1
                        if 'volume' in pmr_plate[well]:
                            primer_vols[primer_name] += max(pmr_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0)
        return primer_counts, primer_vols


    def get_taqwater_avail(self, taqwater_bcs=None, transactions=None):
        """ 
        Returns (int) taq and (int) water volumes (in nanolitres) loaded as available 
        Will work from a list of barcodes if provided, or will calculate for all available plates
        transactions (Optional): a dictionary of plates and their unprocessed changes for accurate calculation
        """
        taq_avail = 0
        water_avail = 0
        if taqwater_bcs is not None:
            pids = [p for p in taqwater_bcs if self.plate_location_sample[p]['purpose'] == 'taq_water']
        else:
            pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'taq_water']

        for pid in pids:
            tw_plate = self.get_plate(pid) # get the plate records with adjusted usage
            for well in ['A1','A2','A3']:
                water_avail += tw_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo6']]
                if transactions is not None:
                    if pid in transactions:
                        if well in transactions[pid]:
                            water_avail += transactions[pid][well] # the transactions are recorded as change, so will be negative values

            for well in ['B1','B2','B3']:
                taq_avail += tw_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo6']]
                if transactions is not None:
                    if pid in transactions:
                        if well in transactions[pid]:
                            taq_avail += transactions[pid][well]
        return taq_avail, water_avail, pids

    def get_taqwater_pids(self):
        pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'taq_water']
        return pids


    def generate_nimbus_inputs(self):
        success = file_io.nimbus_gen(self)
        return success


    def get_nimbus_filepaths(self):
        """ Return the lists of nimbus input files, echo input file (nimbus outputs), 
            and barcodes that are only seen in nimbus """
        #print("In get_nimbus_filepaths")
        nimbus_input_filepaths, echo_input_paths, xbc = file_io.match_nimbus_to_echo_files(self)
        return nimbus_input_filepaths, echo_input_paths, xbc


    def remove_entries(self, selected_rows):
        """ remove JSON elements in selection from experiment
            [{"DNA PID":"p12202p"
            "Sample PID1":"p71561p"
            "Sample PID2":""
            "Sample PID3":""
            "Sample PID4":""
            "Source":"Musterer"
            "Well count":83
            "Unique samples":83
            "Unique assays":15
            }]
            Returns True on success
        """                                                                                              
        #print(f"In remove_entries. {selected_rows=}", file=sys.stderr)
        if type(selected_rows) is dict:
            selected_rows = [selected_rows]
        for row in selected_rows:
            dest_pid = util.guard_pbc(row['DNA PID'], silent=True)
            if dest_pid not in self.dest_sample_plates:
                self.log(f"Error: {row['DNA PID']=} doesn't actually exist in the experiment!")
                continue
            sample_pids = self.dest_sample_plates[dest_pid]
            delete_pids = sample_pids + [dest_pid]
            self.delete_plates(delete_pids)
            del self.dest_sample_plates[dest_pid]
            #print(f'remove_entries() {self.dest_sample_plates=}', file=sys.stderr)
        self.save()
        return True

    def add_references(self, uploaded_references):
        """
        read in reference (target) IDs and sequences.
        May be added when pipeline is locked.
        """
        #if self.locked:
        #    self.log('Error: Cannot add reference sequences while lock is active')
        #    self.save()
        #    return False
        partial_fail = False
        for uploaded_reference in uploaded_references:
            ref_name = uploaded_reference.name
            if ref_name in self.reference_sequences:
                self.log(f"Warning: Duplicate reference file name: {ref_name}. Overwriting")
            try:
                with StringIO(uploaded_reference.getvalue().decode()) as sio:
                    self.reference_sequences[ref_name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
            except UnicodeDecodeError:
                try:
                    with StringIO(uploaded_reference.getvalue().decode('utf-8')) as sio:
                        self.reference_sequences[ref_name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
                except UnicodeDecodeError:
                    try:
                        with StringIO(uploaded_reference.getvalue().decode('latin-1')) as sio:
                            self.reference_sequences[ref_name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
                    except UnicodeDecodeError:
                        try:
                            with StringIO(uploaded_reference.getvalue().decode('ISO-8859')) as sio:
                                self.reference_sequences[ref_name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
                        except Exception as exc:
                            self.log(f"Error: Couldn't parse reference file {ref_name} {exc}")
                            self.save()
                            partial_fail = True
                            continue                
            self.log(f'Success: uploaded {len(self.reference_sequences[ref_name])} reference sequences from {ref_name}')
            if ref_name in self.uploaded_files:
                self.log(f'Warning: {ref_name} already present in records. Overwriting in file list')
            self.uploaded_files[ref_name] = {'plates':[], 'purpose': 'reference sequences'}
        self.save()
        if partial_fail:
            return False
        
        return True


    def generate_targets(self):
        """ create target file based on loaded references """
        transactions = {}
        target_fn = self.get_exp_fp('targets.fa', transaction=True)
        transactions[target_fn] = {} # add plates and modifications to this
        counter = 0
        try:
            with open(target_fn, 'wt') as targetf:
                for group in self.reference_sequences:
                    for id in self.reference_sequences[group]:
                        print(f'>{id}', file=targetf)
                        print(f'{self.reference_sequences[group][id]}', file=targetf)
                        counter += 1    
        except Exception as exc:
            self.log(f'Critical: could not write reference sequences to {target_fn} {exc}')
            self.save()
            return False
        self.add_pending_transactions(transactions)
        self.accept_pending_transactions()
        self.log(f'Success: created reference sequences file {target_fn} containing {counter} sequences')
        self.save()
        return True
        

    def add_assaylists(self, uploaded_assaylists):
        """ mapping of assay family to primer family """
        if self.locked:
            self.log('Error: cannot add assay list while lock is active.')
            return False
        for uploaded_assaylist in uploaded_assaylists:
            file_name = uploaded_assaylist.name
            data = StringIO(uploaded_assaylist.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    continue  # header
                if row == '' or row[0] == '' or row[1] == '':
                    continue
                p = row[1]
                a = row[0]
                if p in self.primer_assay and self.primer_assay[p] != a:
                    self.log(f'Warning: Existing primer:assay pair {p}:{self.primer_assay[p]} being overwritten by {p}:{a}')
                self.primer_assay[row[1]] = row[0]
            if file_name in self.uploaded_files:
                self.log(f'Warning {file_name} already present in recorded files, overwriting')
            self.uploaded_files[file_name] = {'plates':[], 'purpose':'Assay list'}
        self.log(f"Success: added primer-assay lists: {', '.join([ual.name for ual in uploaded_assaylists])}")
        self.save()
        return True


    def add_primer_layouts(self, uploaded_primer_layouts):
        """ add primer plate definition with well and name columns """
        if self.locked:
            self.log('Error: cannot add primer plates while lock is active.')
            return False
        for uploaded_primer_layout in uploaded_primer_layouts:
            upl_name = uploaded_primer_layout.name
            PID = upl_name.split('_')[0]   # assumes first field is PID
            gPID = util.guard_pbc(PID, silent=True)   
                
            if gPID in self.plate_location_sample:
                if self.plate_location_sample[gPID]['purpose'] == 'primer':
                    self.log(f"Info: Primer file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Primer file with PID {PID} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[gPID]}")
                    return False
            else:  # create new plate entry and set purpose
                self.plate_location_sample[gPID] = {'purpose':'primer', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
                self.log(f"Info: Creating new primer plate record for {PID}")
            # register file and plates with self.uploaded_files
            if upl_name in self.uploaded_files:
                self.log(f"Warning: file {upl_name} has already been uploaded, overwriting")
            self.uploaded_files[upl_name] = {'plates': [gPID], 'purpose': "primer layout"}
              
            # load data into plate
            data = StringIO(uploaded_primer_layout.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    continue  # header
                if row == '' or row[0] == '' or row[1] == '':
                    continue
                well = util.unpadwell(row[0])
                if well not in self.plate_location_sample[gPID]:
                    self.plate_location_sample[gPID]['wells'].add(well)
                    self.plate_location_sample[gPID][well] = {}
                #print(row[0], row[1])
                self.plate_location_sample[gPID][well]['primer'] = row[1]
        self.log(f"Success: added primer layouts from {', '.join([upl.name for upl in uploaded_primer_layouts])}")
        self.save()
        return True


    def add_primer_volumes(self, uploaded_primer_volumes):
        """ add primer plate volumes with well and volume columns """
        if self.locked:
            self.log('Error: cannot add primer plate volumes while lock is active.')
            return False
        for uploaded_primer_volume in uploaded_primer_volumes:
            upv_name = uploaded_primer_volume.name
            PID = uploaded_primer_volume.name.split('_')[0]  # assumes first field is PID
            gPID = util.guard_pbc(PID, silent=True)
            if gPID in self.plate_location_sample:
                if self.plate_location_sample[gPID]['purpose'] == 'primer':
                    self.log(f"Info: Primer file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Primer file with PID {PID} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[gPID]}")
                    return False
            else:  # create new plate entry and set purpose
                self.plate_location_sample[gPID] = {'purpose':'primer', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
                self.log(f"Info: Creating new primer plate record for {PID}")
            # register file and plates with self.uploaded_files
            if upv_name in self.uploaded_files:
                self.log(f"Warning: file {upv_name} has already been uploaded, overwriting")
            self.uploaded_files[upv_name] = {'plates': [gPID], 'purpose': "primer volumes"}

            # load data into plate
            data = StringIO(uploaded_primer_volume.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    if row[0].lower().startswith('date'):
                        # plate style layout
                        layout = 'plate'
                    else:
                        layout = 'columns'
                    continue  # header
                
                if layout == 'plate':
                    if not row or row[0] == '' or row[0][0] not in {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'}:
                       continue
                    for j,col in enumerate(row):
                        if j==0:
                            continue  # row name cell
                        well = row[0] + str(j)
                        if well not in self.plate_location_sample[gPID]:
                            self.plate_location_sample[gPID]['wells'].add(well)
                            self.plate_location_sample[gPID][well] = {}
                        self.plate_location_sample[gPID][well]['volume'] = float(col)*1000
                else:
                    well = util.unpadwell(row[0])
                    if well not in self.plate_location_sample[gPID]:
                        self.plate_location_sample[gPID]['wells'].add(well)
                        self.plate_location_sample[gPID][well] = {}
                    self.plate_location_sample[gPID][well]['volume'] = float(row[1])*1000
        self.log(f"Success: added primer volumes from {', '.join([upv.name for upv in uploaded_primer_volumes])}")
        self.save()
        return True


    def add_amplicon_manifests(self, uploaded_amplicon_manifests):
        """
        Amplicon plate layouts (containing pre-amplified sequences) are added here via manifest files.
        Only column layout (comma separate) is supported. col1: plate barcode; col2: well position; col3: sample barcode; col4 (optional): volume.
        The first row must be a header row: plate, well, sample, (volume in uL, if provided).
        plateBarcode, well, sampleBarcode, amplicon* (*multiple columns with this name are allowed)
        Adding this will result in changing the Miseq and Stage3 output files, so needs to be wrapped as a transaction
        Needs to check whether there are sufficient indexes 
        """
        if self.locked:
            self.log('Error: cannot add amplicon plate layouts while lock is active.')
            return False

        file_names, file_tables = file_io.read_csv_or_excel_from_stream(uploaded_amplicon_manifests)
        for file_name, file_table in zip(file_names, file_tables):
            if file_name in self.uploaded_files:
                self.log(f'Warning: {file_name} already present in records')
                
            plate_entries = {}  # looks like self.plate_location_sample, but temporary
            plate_barcode_col = None
            well_col = None
            amplicon_cols = []  # combine these
            for i, line in enumerate(file_table): # csv_reader is super limiting
                cols = [c.strip() if c is not None else '' for c in line.split(',')]  # lower case column names
                #print(i, cols, file=sys.stderr)
                if i==0:  # process header
                    cols_lower = [c.lower() for c in cols]
                    header_dict = {k:c for k,c in enumerate(cols_lower)}
                    matching_cols = [col_name in cols_lower for col_name in ['platebarcode', 'well', 'samplebarcode']] 
                    if not all(matching_cols):
                        self.log(f'Error: amplicon manifest {file_name} requires at least columns plateBarcode, well, sampleBarcode')
                        return False
                    else:
                        self.log(f'Info: parsing amplicon manifest {file_name}')
                    for k in header_dict:
                        if header_dict[k] == 'platebarcode':
                            plate_barcode_col = k
                        elif header_dict[k] == 'well':
                            well_col = k
                        elif header_dict[k] == 'amplicon': # combine amplicon columns
                            amplicon_cols.append(k)
                    continue
                if len(cols) < 3:  # skip empty rows
                    continue
                # set up plate
                gpid = util.guard_pbc(cols[plate_barcode_col], silent=True)
                if gpid not in plate_entries:
                    # create plate_location_sample entries here and copy them across when we've got them all
                    plate_entries[gpid] = {'purpose':'amplicon','source':'manifest', 'wells':set(), 
                            'plate_type':'384PP_AQ_BP', 'barcode':gpid}
                # set up well
                well = util.unpadwell(cols[well_col].upper())
                if well in plate_entries[gpid]['wells'] and plate_entries[gpid][well] != {}:
                    self.log(f"Error: duplicate {c} in {gpid}. Skipping row {i+2} {cols=}")
                    continue
                plate_entries[gpid]['wells'].add(well)
                plate_entries[gpid][well] = {}
                amplicons = []
                # now collect everything together
                for k,c in enumerate(cols):
                    if c.lower() == 'none':
                        c = ''
                    if header_dict[k] == 'well':
                        continue  # we've already got this!
                    if header_dict[k] == 'samplebarcode':
                        # amplicon guards
                        sid = plate_entries[gpid][well][header_dict[k]] = util.guard_abc(c, silent=True)
                        plate_entries[gpid][well]['barcode'] = sid
                    elif header_dict[k] == 'platebarcode':
                        # don't really need it but whatever
                        plate_entries[gpid][well]['platebarcode'] = gpid
                    elif header_dict[k] == 'sampleno':
                        plate_entries[gpid][well]['sampleNumber'] = str(c)
                    elif header_dict[k] == 'amplicon':
                        if c != '':  # ignore empty amplicon entries
                            amplicons.append(c)
                    elif header_dict[k] == 'volume':
                        plate_entries[gpid][well]['volume'] = float(c)*1000  # save as nL
                    else:
                        try:
                            plate_entries[gpid][well][header_dict[k]] = str(c)
                        except:
                            print(type(header_dict[k]), header_dict[k], str(c), file=sys.stderr)
                plate_entries[gpid][well]['amplicons'] = amplicons

            # Check whether we have existing records and clean up if necessary
            # For amplicons which are added late in the pipeline, this may mean invalidating the final stages of the pipeline
            kept_pids = []
            clashing_unassigned_pids = []
            clashing_existing_pids = []
            for gpid in plate_entries:
                if gpid in self.unassigned_plates: 
                    self.log(f'Warning: uploaded records exist for {util.unguard_pbc(gpid, silent=True)}, overwriting existing plate')
                    clashing_unassigned_pids.append(gpid)
                    kept_pids.append(gpid)
                elif gpid in self.plate_location_sample:
                    existing_purpose = self.plate_location_sample[gpid]['purpose']
                    if existing_purpose == 'amplicon':
                        self.log(f'Warning: Amplicon plate records exist for {util.unguard_pbc(gpid, silent=True)}, overwriting existing plate')
                        clashing_existing_pids.append(gpid)
                        kept_pids.append(gpid)
                    else:
                        self.log(f'Error: Plate {util.unguard_pbc(gpid, silent=True)} already exists with purpose {existing_purpose}, skipping this plate')
                        continue
                else:
                    kept_pids.append(gpid)
            for cp in clashing_unassigned_pids:
                del self.unassigned_plates[cp]
            for cp in clashing_existing_pids:
                del self.plate_location_sample[cp]
            for kp in kept_pids:
                self.plate_location_sample[kp] = plate_entries[gpid]
            self.uploaded_files[file_name] = {'plates':kept_pids, 'purpose':'amplicon manifest'}
            self.log(f"Success: added amplicon plate info from {file_name} for plates '+\
                    f'{', '.join([util.unguard_pbc(kp, silent=True) for kp in kept_pids])}")
       
        self.save()
        return True


    def get_amplicon_pids(self):
        """ return a list of user supplied amplicon plate ids """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'amplicon']


    def generate_echo_primer_survey(self, primer_survey_filename='primer-svy.csv'):
        """ 
        Generate a primer survey file for use by the Echo. Replaces echovolume.py
        Not strictly necessary, but the old code reads this file in making the picklists.
        """
        if self.locked:
            self.log('Error: cannot generate primer survey file while lock is active.')
            return False
        primer_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'primer']
        header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'plate position on Echo 384 PP',
                'primer names pooled', 'volume']
        transactions = {}                                       
        primer_survey_fn = self.get_exp_fp(primer_survey_filename, transaction=True)
        transactions[primer_survey_fn] = {} # add plates and modifications to this
        try:
            with open(primer_survey_fn, 'wt') as fout:
                print(','.join(header), file=fout)
                for i,pid in enumerate(primer_pids):
                    transactions[primer_survey_fn][pid] = {}
                    plate=self.plate_location_sample[pid]
                    for well in util.row_ordered_384:
                        if well not in plate['wells']:
                            continue
                        if 'primer' not in plate[well]:
                            continue
                        if 'volume' not in plate[well]:
                            continue
                        outline = ','.join([f"Source[{i+1}]",util.unguard_pbc(pid,silent=True),plate['plate_type'],
                                well,plate[well]['primer'],f"{int(plate[well]['volume'])/1000}"])
                        print(outline, file=fout)
        except Exception as exc:
            self.log(f'Failure: could not write primer survey {exc}')
            self.save()
            return False
        self.add_pending_transactions(transactions)
        self.log(f"Success: written Echo primer survey to {primer_survey_fn}")
        return True


    def check_plate_presence(self, pids, purpose, messages):
        """ 
        Check given plate for presence in self.plate_location_sample and compare purpose
        return success (and messages by reference)
        Required by exp.check_ready_pcr1() and exp.check_ready_pcr2()
        """
        success = True  
        for pid in pids:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != purpose:
                    msg = f"Error: plate already exists with PID {pid} with purpose "+\
                            f"{self.plate_location_sample[pid]['purpose']}, expected {purpose}"
                    self.log(msg)
                    messages.append(msg)
                    success = False
            else:
                msg = f"Critical: No plate exists with barcode {pid}"
                self.log(msg)
                messages.append(msg)
                success = False
        return success


    def check_ready_pcr1(self, dna_plates, pcr_plates, taq_water_plates, messages):
        """
        Check that everything required to successfully generate PCR1 picklists is available
        TODO: Check that volumes and wells are sufficient!
        messages is an optional list, allowing us to pass useful info back to the UI
        Returns success
        """
        success = True
        if not dna_plates or not pcr_plates or not taq_water_plates:
            success = False

        dna_success= self.check_plate_presence(dna_plates, 'dna', messages)
        if not dna_success:
            success = False

        pcr_success = self.check_plate_presence(pcr_plates, 'pcr', messages)
        if not pcr_success:
            success = False

        taq_success = self.check_plate_presence(taq_water_plates, 'taq_water', messages)
        if not taq_success:
            success = False
      
        return success


    def check_ready_pcr2(self, pcr_plates, taq_water_plates, index_plates, amplicon_plates, messages):
        """
        Check the everything required to successfully generate PCR2 picklists is available
        TODO: Check that volumes and wells are sufficient!
        messages is a list, allowing us to pass useful info back to the UI
        Return success
        """
        success = True

        if not pcr_plates or not taq_water_plates or not index_plates:
            success = False

        pcr_success = self.check_plate_presence(pcr_plates, 'pcr', messages)
        if not pcr_success:
            success = False

        taq_success = self.check_plate_presence(taq_water_plates, 'taq_water', messages)
        if not taq_success:
            success = False

        index_success = self.check_plate_presence(index_plates, 'index', messages)
        if not index_success:
            success = False

        amplicon_success = self.check_plate_presence(amplicon_plates, 'amplicon', messages)
        if not amplicon_success:
            success = False
      
        return success


    def generate_echo_PCR1_picklists(self, dna_plates, pcr_plates, taq_water_plates):
        """
        Calls echo_primer.generate_echo_PCR1_picklist() to do the work, needs a set of accepted DNA_plates,
        the final assay list, and a set of destination PCR plate barcodes, taq+water plates, primer plates, primer volumes.
        Returns True on success
        TODO: We need to reduce available taq and water during this operation.
        """
        #print(f"Experiment.generate_echo_PCR1 {dna_plates=} {pcr_plates=} {taq_water_plates=}")
        
        #success = self.generate_echo_primer_survey()
        #if not success:
        #    self.log('Failure: failed to generate primer survey file')
        #    return False
        # do transaction handling in generate_echo_PCR1_picklist()
        success = file_io.generate_echo_PCR1_picklist(self, dna_plates, pcr_plates, taq_water_plates)
        if not success:
            self.log('Failure: could not generate PCR1 picklists correctly')
            return False
        self.log('Success: generated PCR1 picklists')
        return True

    def get_echo_PCR1_picklist_filepaths(self):
        """
        Return file paths for PCR1_dna-picklist_XXX.csv, PCR1_primer-picklist_XXX.csv, PCR1_taqwater-picklist_XXX.csv
        """
        all_files = os.listdir(self.get_exp_dir())
        dna_picklist_paths = []
        primer_picklist_paths = []
        taqwater_picklist_paths = []
        for f in all_files:
            if f.startswith('PCR1_dna-picklist_') and f.endswith('.csv'):
                dna_picklist_paths.append(self.get_exp_fp(f))
            elif f.startswith('PCR1_primer-picklist_') and f.endswith('.csv'):
                primer_picklist_paths.append(self.get_exp_fp(f))
            elif f.startswith('PCR1_taqwater-picklist_') and f.endswith('.csv'):
                taqwater_picklist_paths.append(self.get_exp_fp(f))
        return dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths


    def generate_echo_PCR2_picklists(self, pcr_plates, index_plates, taq_water_plates, amplicon_plates=None):
        """
        Calls echo_index.generate_echo_PCR2_picklist() to do the work, needs a set of destination PCR plate barcodes, 
        taq+water plates, index plates, index volumes, and optionally any amplicon plates (post-PCR).
        Returns True on success
        """
        for pid in pcr_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'pcr':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
        for pid in index_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'index':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
        for pid in taq_water_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'taq_water':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
        for pid in amplicon_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'amplicon':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
            
        success = self.generate_echo_index_survey(index_plates)
        if not success:
            self.log('Failure: failed to generate index survey file')
            return False
        self.log('Success: generated index survey file')
        success = file_io.generate_echo_PCR2_picklist(self, pcr_plates, index_plates, taq_water_plates, amplicon_plates)
        if not success:
            self.log('Failure: could not generate PCR2 (index) picklists correctly')
            return False
        self.log('Success: generated PCR2 (index) picklists')
        return success

    def get_echo_PCR2_picklist_filepaths(self):
        """
        Return file paths for PCR2_index-picklist_XXX.csv, PCR2_taqwater-picklist_XXX.csv, 
                (optionally) PCR1_amplicon-picklist_XXX.csv (per amplicon plate?)
        """
        all_files = os.listdir(self.get_exp_dir())
        index_picklist_paths = []
        amplicon_picklist_paths = []
        taqwater_picklist_paths = []
        for f in all_files:
            if f.startswith('PCR2_index-picklist_') and f.endswith('.csv'):
                index_picklist_paths.append(self.get_exp_fp(f))
            elif f.startswith('PCR2_amplicon-picklist_') and f.endswith('.csv'):
                amplicon_picklist_paths.append(self.get_exp_fp(f))
            elif f.startswith('PCR2_taqwater-picklist_') and f.endswith('.csv'):
                taqwater_picklist_paths.append(self.get_exp_fp(f))
        return index_picklist_paths, taqwater_picklist_paths, amplicon_picklist_paths

    def delete_plates(self, pids):
        """
        Soft-delete plates with the selected pids
        """
        for pid in pids:
            if pid in self.plate_location_sample:
                if pid not in self.deleted_plates:
                    self.deleted_plates[pid] = []
                self.deleted_plates[pid].append(deepcopy(self.plate_location_sample[pid]))
                # need to find this in reproducible steps and delete
                del self.plate_location_sample[pid]
                self.log(f'Warning: moved {pid} to deleted plate bin')
            else:
                self.log(f'Warning: {pid} has no definition loaded')
        self.save()

    def get_stages(self):
        """ get all information on reproducible steps and pending steps for display purposes """
        header = ['Stage order','Staged file', 'Affected plates', 'Status']
        stages = []
        counter = 1
        for steps_dict in self.reproducible_steps:
            for file_name in steps_dict:
                pids = [util.unguard_pbc(pid, silent=True) for pid in steps_dict[file_name].keys()]
                stages.append([str(counter), file_name, ', '.join(pids), 'committed'])
            counter += 1
        if self.pending_steps:
            for file_name in self.pending_steps:
                pids = [util.unguard_pbc(pid, silent=True) for pid in self.pending_steps[file_name].keys()]
                stages.append([str(counter), file_name, ', '.join(pids), 'pending'])
        return stages, header

    def get_stage2_pcr_plates(self):
        """
        Used by Indexing stage to find all the used PCR plate IDs
        """
        stage2_pcr_plates = set()
        stage2_fn = self.get_exp_fp('Stage2.csv')
        if not os.path.exists(stage2_fn):
            return stage2_pcr_plates

        with open(stage2_fn, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue  # header
                cols = line.split(',')
                stage2_pcr_plates.add(cols[8].strip())
        return stage2_pcr_plates

    def get_file_usage(self):
        file_usage = {}
        for filename in self.uploaded_files:
            file_usage[filename] = {'plates' :[], 'purpose':None}
            if 'purpose' in self.uploaded_files[filename]:
                file_usage[filename]['purpose'] = self.uploaded_files[filename]['purpose']
            if 'plates' in self.uploaded_files[filename]:
                for pid in self.uploaded_files[filename]['plates']:
                    file_usage[filename]['plates'].append(util.unguard_pbc(pid, silent=True))
            file_usage[filename]['plates'] = ', '.join(file_usage[filename]['plates'])
        return file_usage

    def get_plate_usage(self):
        """return list of plate information to display"""
        plate_usage = []
        for p in self.plate_location_sample:
            if 'purpose' not in self.plate_location_sample[p]:
                continue            
            pid = util.unguard_pbc(p, silent=True)
            plate_usage.append([pid, len(self.plate_location_sample[p]['wells']), self.plate_location_sample[p]['purpose']])
        return plate_usage

    def get_index_pids(self):
        """
        Used by Indexing stage
        """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'index']

  


    def add_index_layouts(self, uploaded_index_layouts):
        """ add index plates with well and index columns """
        if self.locked:
            self.log('Error: cannot add index plates while lock is active.')
            return False
        for uploaded_index_layout in uploaded_index_layouts:
            layout_name = uploaded_index_layout.name
            PID = layout_name.split('_')[0]  # assumes first field is PID
            gPID = util.guard_pbc(PID, silent=True)   
            if gPID in self.plate_location_sample:
                if self.plate_location_sample[gPID]['purpose'] == 'index':
                    self.log(f"Info: Index plate {PID} already exists, adding data")
                else:
                    self.log(f"Error: Index plate PID: {PID} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[gPID]['purpose']}")
                    return False
            else:  # create new plate entry and set purpose
                self.plate_location_sample[gPID] = {'purpose': 'index', 'source':'user', 'wells':set(), 
                        'plate_type':util.PLATE_TYPES['Echo384'], 'barcode':gPID}
                self.log(f"Info: Creating new index plate record for {PID}")
            # register file and plates with self.uploaded_files
            if layout_name not in self.uploaded_files:
                self.log(f'Warning: File name {layout_name} already exists, overwriting')
            self.uploaded_files[layout_name] = {'plates': [gPID], 'purpose': "index layout"} 
                
            # load data into plate
            data = StringIO(uploaded_index_layout.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0 or row == '':
                    continue  # header or blank
                well = util.unpadwell(row[0])
                if well not in self.plate_location_sample[gPID]:
                    self.plate_location_sample[gPID]['wells'].add(well)
                    self.plate_location_sample[gPID][well] = {}
                    idt_name = row[1]
                    index = row[2]
                    bc_name = row[3]
                    oligo = row[4]
                    self.plate_location_sample[gPID][well]['idt_name'] = idt_name
                    self.plate_location_sample[gPID][well]['index'] = index
                    self.plate_location_sample[gPID][well]['bc_name'] = bc_name
                    self.plate_location_sample[gPID][well]['oligo'] = oligo
                    
        self.log(f"Success: added index plate layouts from {', '.join([uil.name for uil in uploaded_index_layouts])}")        
        self.save()
        return True

    def add_index_volumes(self, uploaded_index_volumes):
        """ add volume information to a single barcode plate, with a matched, existing plate name.
        To make life fun there are two possible formats:
        Format 1: plate layout (generated through Echo test software)
        Format 2: column layout (generated through Echo main interface)
        Returns True on success
        """
        if self.locked:
            self.log('Error: cannot add index plate volumes while lock is active.')
            return False
        for uploaded_index_volume in uploaded_index_volumes:
            volume_name = uploaded_index_volume.name
            PID = uploaded_index_volume.name.split('_')[0]  # assumes first field is PID
            gPID = util.guard_pbc(PID, silent=True)   
                
            if gPID not in self.plate_location_sample:
                self.log(f"Info: Adding index plate {PID}")
                self.plate_location_sample[gPID] = {'purpose':'index', 'source':'user', 'wells':set(), 
                        'plate_type':util.PLATE_TYPES['Echo384'], 'barcode':gPID}
            else:
                if self.plate_location_sample[gPID]['purpose'] != 'index':
                    self.log(f"Error: {PID} plate purpose is "+\
                            f"{self.plate_location_sample[gPID]['purpose']}, expected 'index'")
                    return False
                self.log(f"Info: {PID} exists, appending index volumes")
            # register file and plate with self.uploaded_files
            if volume_name not in self.uploaded_files:
                self.log(f'Warning: File name {volume_name} already exists, overwriting')
            self.uploaded_files[volume_name] = {'plates': [gPID], 'purpose': "index volumes"}
        
            plate_format = False
            data = StringIO(uploaded_index_volume.getvalue().decode("utf-8"), newline='')
            hdr_seen = False
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    if row[0].lower().startswith('date'):
                        plate_format = True
                    continue
                if plate_format:
                    # format 1 - volume in r,c plate matrix
                    # matrix format - 2 initial lines, a blank line then a matrix
                    #next(src), next(src) # skip two more lines
                    if not row or row[0].strip() == '':
                        continue
                    elif row[0].strip() == 'A':
                        hdr_seen = True
                    elif not hdr_seen:
                        continue
                    
                    plate_row = row[0].strip()
                    for col, v in enumerate(row):
                        if col > 0 and v.strip() != '':
                            well = plate_row+str(col)
                            self.plate_location_sample[gPID]['wells'].add(well)
                            if well not in self.plate_location_sample[gPID]:
                                self.plate_location_sample[gPID][well] = {}
                            self.plate_location_sample[gPID][well]['volume'] = float(v)*1000
                else:                
                    #format 2 - one row per well - well ID in r[3], volume in r[5] or r[6]
                    #if i == 0 and line.startswith('Run ID'):
                    if row.startswith('[DETAILS]'):
                        continue
                    if ''.join(row).strip() == '':
                        continue
                    if not hdr_seen:
                        self.hdr = row
                        hdr_seen = True
                    else:
                        well = util.unpadwell(row[3])
                        self.plate_location_sample[gPID]['wells'].add(well)
                        if well not in self.plate_location_sample[gPID]:
                            self.plate_location_sample[gPID][well] = {}
                        self.plate_location_sample[gPID][well]['volume'] = int(float[5])*1000
                
        self.log(f"Success: added index plate volumes from {', '.join([uiv.name for uiv in uploaded_index_volumes])}")
        self.save()
        return True

    def generate_echo_index_survey(self, user_index_pids, index_survey_filename='index-svy.csv'):
        """ 
        Generate an index survey file for use by the Echo. Replaces echovolume.py 
        Not strictly necessary, but the old code reads this file in making the picklists.
        Called internally by generate_echo_PCR2_picklist_interface()
        Requires the set of user chosen user_index_pids
        """
        if self.locked:
            self.log('Error: cannot generate index survey while lock is active.')
            return False
        index_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'index' and p in user_index_pids]
        header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'Source well', 
                'Name for IDT ordering', 'index', 'Name', 'Oligo * -phosphothioate modif. against exonuclease', 'volume']
        fn = self.get_exp_fp(index_survey_filename)
        if os.path.exists(fn):
            self.log(f'Warning: overwriting index survey file {fn}')
        with open(fn, 'wt') as fout:
            print(','.join(header), file=fout)
            for i,pid in enumerate(index_pids):
                plate=self.plate_location_sample[pid]
                for well in util.row_ordered_384:
                    if well not in plate['wells']:
                        continue
                    if 'idt_name' not in plate[well]:
                        continue
                    if 'volume' not in plate[well]:
                        continue
                    try:
                        outline = ','.join([f'Source[{i+1}]',util.unguard_pbc(pid,silent=True),plate['plate_type'],
                                well,plate[well]['idt_name'],plate[well]['index'],plate[well]['bc_name'],
                                plate[well]['oligo'],str(plate[well]['volume']/1000)])
                    except:
                        print(f"Source[{i+1}]")
                        print(f"{util.unguard_pbc(pid,silent=True)=},{plate['plate_type']=}")
                        print(f"{well=},{plate[well]=}")
                        print(f"{plate[well]['idt_name']=}")
                        print(f"{plate[well]['index']=},{plate[well]['bc_name']=}"+\
                                f"{plate[well]['oligo']=},{plate[well]['volume']=}")
                    print(outline, file=fout)
        self.log(f'Success: index survey file written to {fn}')
        return True

    def add_standard_taqwater_plates(self, plate_barcodes):  # <- they're the same for primer and barcode stages
        """ These plates have a fixed layout with water in row A and taq in row B of a 6-well reservoir plate.
        See echo_primer.py or echo_barcode.py mytaq2()
        We may need to also store volumes at some stage, but for now that isn't necessary
        """
        if self.locked:
            self.log('Error: cannot add taq+water plates while lock is active.')
            return False
        try:
            for plate_barcode in plate_barcodes:
                pid = util.guard_pbc(plate_barcode, silent=True)
                if pid in self.plate_location_sample:
                    if self.plate_location_sample[pid]['purpose'] != 'taq_water':
                        self.log(f"Error: cannot add taq+water plate {util.unguard(pid)}. "+\
                            f"A plate with this barcode already exists with purpose {self.plate_location_sample[pid]['purpose']}")
                        return False
                    else:
                        self.log(f"Warning: overwriting existing taq/water plate entry with barcode {util.unguard(pid)}")
                else:
                    self.plate_location_sample[pid] = {}
                self.plate_location_sample[pid]['purpose'] = 'taq_water'
                self.plate_location_sample[pid]['plate_type'] = util.PLATE_TYPES['Echo6']
                self.plate_location_sample[pid]['capacity'] = util.CAP_VOLS[util.PLATE_TYPES['Echo6']]
                self.plate_location_sample[pid]['dead_vol'] = util.DEAD_VOLS[util.PLATE_TYPES['Echo6']]
                self.plate_location_sample[pid]['wells'] = ['A1','A2','A3','B1','B2','B3']
                self.plate_location_sample[pid]['water_wells'] = util.WATER_WELLS
                self.plate_location_sample[pid]['taq_wells'] = util.TAQ_WELLS
                self.plate_location_sample[pid]['A1'] = {'name': 'water', 'volume': util.CAP_VOLS[util.PLATE_TYPES['Echo6']]}
                self.plate_location_sample[pid]['A2'] = {'name': 'water', 'volume': util.CAP_VOLS[util.PLATE_TYPES['Echo6']]}
                self.plate_location_sample[pid]['A3'] = {'name': 'water', 'volume': util.CAP_VOLS[util.PLATE_TYPES['Echo6']]}
                self.plate_location_sample[pid]['B1'] = {'name': 'taq', 'volume': util.CAP_VOLS[util.PLATE_TYPES['Echo6']]}
                self.plate_location_sample[pid]['B2'] = {'name': 'taq', 'volume': util.CAP_VOLS[util.PLATE_TYPES['Echo6']]}
                self.plate_location_sample[pid]['B3'] = {'name': 'taq', 'volume': util.CAP_VOLS[util.PLATE_TYPES['Echo6']]}
        except Exception as exc:
            self.log(f"Error: adding taq+water plate barcode failed {plate_barcodes=} {exc}")
            return False
        self.save()
        return True

    def get_miseq_samplesheets(self):
        """ return the MiSeq-XXX.csv samplesheet, if it exists """
        miseq_fps = list(Path(self.get_exp_dir()).glob('MiSeq_*.csv'))
        return miseq_fps  


    def check_sequence_upload_ready(self, messages):
        """
        Prevent users from uploading sequence files without references being loaded, or without generating Stage3 or MiSeq files
        Requires a list of messages for the GUI which can be appended to (pass by reference)
        """
        success = True
        # check for MiSeq.csv and Stage3 files
        fn1 = self.get_exp_fp(f'MiSeq_{self.name}.csv')
        fn2 = self.get_exp_fp(f'Stage3.csv')
        fns = [fn1, fn2]
        for fn in fns:
            if not Path(fn).exists():
                success = False
                msg = f"Error: {fn} not present"
                self.log(msg)
                messages.append(msg)

        # check that at least one reference sequence has been uploaded
        if len(self.reference_sequences) == 0:
            success = False
            msg = 'No reference sequences have been uploaded yet, please add these'
            self.log(msg)
            messages.append(msg)

        self.save()
        return success


    def check_allele_calling_ready(self, messages):
        """ 
        Return True if everything needed for allele calling is present
            - Stage3.csv is present - done in self.check_sequence_upload_ready()
            - Miseq file is present - done in self.check_sequence_upload_ready()
            - check that reference sequences are uploaded - done in self.check_sequence_upload_ready()
            - the raw directory is present (can be empty)

        Requires a list "messages" to be given for returning feedback to the GUI
        """
        success = self.check_sequence_upload_ready(messages)

        # check whether the raw directory for FASTQs exists yet - implies at least one FASTQ has been uploaded        
        dn = self.get_exp_fp(f'raw')
        if not Path(dn).exists() or not Path(dn).is_dir():
            success = False
            msg = f"Error: {dn} does not exist"
            self.log(msg)
            messages.append(msg)

        self.save()
        return success


    def read_allele_results(self):
        """
        Read in allele_results
        [SampleNo, EPplate, EPwell, sampleBarcode, assays, assayFamilies, dnaplate, dnawell, 
        primer, pcrplate, pcrwell, i7bc, i7name, i7well, i5bc, i5name, i5well, allele_prop, 
        efficiency, readCount, cleanCount, mergeCount, seqCount*, seqName*]

        """
        pass


    def custom_output_report(exp):
        """
        """
        pass


    def rodentity_output_report(exp):
        """
        """
        pass

     
    def save(self):
        """ save experiment details to self.name/experiment.json, returns True on success and False on fail. Assumes the correct working directory """
        try:
            exp = jsonpickle.encode(self, indent=4, keys=True, warn=True)
            if not self.name:
                return False
            with open(os.path.join('run_'+self.name, EXP_FN), 'wt') as f:
                print(exp, file=f)
        except Exception as exc:
            print(f"Error saving {self.name=} {exc}", file=sys.stderr)
            return False
        return True

    def log(self, message, level=''):
        """ Always add the date/time, function where the log was run and the caller function to the log """
        now = datetime.datetime.now()
        t = f"{now:%Y-%m-%d %H:%M}"
        func = sys._getframe(1).f_code.co_name
        func_line = inspect.getframeinfo(sys._getframe(1)).lineno
        caller = sys._getframe(2).f_code.co_name
        caller_line = inspect.getframeinfo(sys._getframe(2)).lineno

        # levels are: Debug/Info/Warning/Error/Critical/Success/Begin/End/Success/Failure
        if not level:
            if message.lower().startswith('d:') or message.lower().startswith('debug:'):
                message = message.split(':',1)[1].strip()
                level = 'Debug'
            elif message.lower().startswith('i:') or message.lower().startswith('info:'):
                message = message.split(':',1)[1].strip()
                level = 'Info'
            elif message.lower().startswith('w:') or message.lower().startswith('warning:'):
                message = message.split(':',1)[1].strip()
                level = 'Warning'
            elif message.lower().startswith('e:') or message.lower().startswith('error:'):
                message = message.split(':',1)[1].strip()
                level = 'Error'
            elif message.lower().startswith('c:') or message.lower().startswith('critical:'):
                message = message.split(':',1)[1].strip()
                level = 'Critical'
            elif message.lower().startswith('b:') or message.lower().startswith('begin:'):
                message = message.split(':',1)[1].strip()
                level = 'Begin'
            elif message.lower().startswith('n:') or message.lower().startswith('end:'):
                message = message.split(':',1)[1].strip()
                level = 'End'
            elif message.lower().startswith('s:') or message.lower().startswith('success:'):
                message = message.split(':',1)[1].strip()
                level = 'Success'
            elif message.lower().startswith('f:') or message.lower().startswith('failure:'):
                message = message.split(':',1)[1].strip()
                level = 'Failure'

        if level.lower() == 'd':
            level = 'Debug'
        elif level.lower() == 'i':
            level = 'Info'
        elif level.lower() == 'w':
            level = 'Warning'
        elif level.lower() == 'e':
            level = 'Error'
        elif level.lower() == 'c':
            level = 'Critical'
        elif level.lower() == 'b':
            level = 'Begin'
        elif level.lower() == 'n':
            level = 'End'
        elif level.lower() == 's':
            level = 'Success'
        elif level.lower() == 'f':
            level = 'Failure'

        #if level == '':
        #    level = 'Debug'
        self.log_entries.append([t, func, func_line, caller, caller_line, level, message])
        if (now - self.log_time).seconds > 10 or level in ['Error', 'Critical', 'End', 'Success', 'Failure']:
            self.save()

    def get_log_header(self):
        return ['Time', 'Function name', 'Func line', 'Calling function', 'Call line', 'Level', 'Message']

    def get_log(self, num_entries=10):
        """ return a chunk of the log. -1 gives everything """
        if num_entries == -1:
            return self.log[::-1]
        return self.log_entries[:-(num_entries+1):-1]

    def clear_log(self, message_type=None):
        """ delete all the log entries with level='Debug' """
        if message_type is not None:
            if message_type not in ['Debug','Info','Warning','Error','Critical','Begin','End','Success','Failure']:
                self.log_entries.append(f"Error: didn't recognise message type {message_type}")
                self.save()
                return
        clean_log = []
        if message_type is not None:
            for entry in self.log_entries:
                if entry[-2] != message_type:
                    clean_log.append(entry)
        self.log_entries = clean_log
        self.save()

#def save_experiment(experiment, exp_path):
#    """ save experiment details to experiment.json """
#    exp_file_path = os.path.join(exp_path, EXP_FN)
#    exp = jsonpickle.encode(experiment, indent=4)
#    with open(exp_file_path, 'wt') as f:
#        print(exp, f)


def load_experiment(exp_path):
    """ load experiment details from experiment.json, or return None. exp_path is the folder path only """
    
    exp_file_path = os.path.join(exp_path, EXP_FN)
    print('Attempting to load from path:', exp_file_path)
    if os.path.exists(exp_file_path):
        with open(exp_file_path, 'rt') as f:
            expt = jsonpickle.decode(f.read(), keys=True) # this seems totally flaky, but it must come down to what the contents are?
        if isinstance(expt, Experiment):
            return expt
        elif isinstance(expt, str):
            expt = dict(expt)
        exp = Experiment(expt['name'])
        exp.__setstate__(expt)
        return exp
    else:
        return None