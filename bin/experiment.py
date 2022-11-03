#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: 22 Feb 2020
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New, replaces template_files.ini
@last_edit: 2022-02-22
@edit_comment: 

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
from collections import Counter
from io import StringIO
import datetime
import json
from Bio import SeqIO
from copy import deepcopy
from math import ceil, floor
from pathlib import Path
import inspect

import pandas as pd
from streamlit import _update_logger

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
    Tracks information on the following file paths (when available)
    mouse pipeline files:
    mouse_ref - defaults to library .fa or.txt (.fa takes precedence)
    
    assay_list - defaults to library
    conversions - defaults to library (this file may on borrowed time thanks to Rodentity)

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
        
    """
    def __init__(self, name):
        """ set up defaults, an experiment name is required, should match folder name less then run_ prefix """
        self.name = name
        self.description = ''
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
        self.log('Info: Locking experiment. Not modification should take place while lock remains')
        self.locked = True

    def unlock(self):
        self.log('Info: Unlocking experiment. Modification is now possible. There should be good reason for this!')
        self.locked = False

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
            print('These are the pending transactions', self.pending_steps, file=sys.stderr)
        return True


    def clashing_pending_transactions(self):
        """
        We need the user to know if accepting pending transactions will result in overwriting an existing
        file. We return the list of all clashing filepaths
        """
        clashes = []
        if self.pending_steps is None:
            print('no pending steps')
            return clashes
        
        for transaction in self.pending_steps:
            p = Path(transaction)
            parent = p.parent
            file_name = str(p.name)
            final_name = file_name[len('pending_'):]  # cut off the leading "pending_"
            final_path = str(parent / final_name)
            #print('Pending and final paths: ', str(p), str(final_path), file=sys.stderr)
            if Path(final_path).exists():
                clashes.append(final_path)
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
        if not self.pending_steps:
            self.log("Warning: there are no pending transactions to record")
            return True  # It didn't actually fail

        clashes = self.clashing_pending_transactions()
        if len(clashes) == 0:
            for transaction in self.pending_steps:
                p = Path(transaction)
                if not p.exists():
                    self.log(f'Warning: {str(p)} not found')
                    continue
                parent = p.parent
                file_name = str(p.name)
                final_name = file_name[len('pending_'):]
                final_path = parent / final_name
                os.rename(p, final_path)
            return True

        clashing_index = -1
        for i,step in enumerate(self.reproducible_steps):
            for dest in step:
                if dest in self.pending_steps:
                    clashing_index = i
                    break
        
        if clashing_index == -1:  # this should NEVER happen
            self.log(f"Critical: pipeline detects clashing transaction {clashing_index=} for {self.pending_steps=}") 
            return False

        # keep everything prior to the clash, then add on the pending steps
        remove_these_steps = self.reproducible_steps[clashing_index:]
        for step in remove_these_steps:
            for fp in step:
                if Path(fp).exists():
                    os.remove(fp)
        # rename pending filepaths
        for transaction in self.pending_steps:
            p = Path(transaction)
            if not p.exists():
                self.log(f'Warning: {str(p)} not found')
                continue
            parent = p.parent
            file_name = str(p.name)
            final_name = file_name[len('pending_'):]
            final_path = parent / final_name
            os.rename(p, final_path)
        self.reproducible_steps = self.reproducible_steps[0:clashing_index]
        self.reproducible_steps.append(self.pending_steps.copy())
        self.pending_steps = None
        return True


    def get_plate(self, PID):
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
                        mod_plate[well] += stage[fn][PID][well]
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
                else:
                    self.log("Error: JSON file doesn't exist: " + fn)
                    return False

            # no failures so update the experiment
            self.dest_sample_plates[dna_plate_id] = sample_plate_ids
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
                    assays = []
                    assayFamilies = set()
                    for allele in record['mouse']['alleles']:  # allele is a dict
                        for assay in allele['assays']:  # assay is a dict
                            assays.append(assay['name'])
                            assayFamilies.add(assay['name'].split('_')[0])
                    self.plate_location_sample[spid][pos]['assays'] = assays.copy()
                    self.plate_location_sample[spid][pos]['assayFamilies'] = list(assayFamilies)
                
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
                else:
                    self.log("Error: JSON file doesn't exist: " + fn)
                    return False

            # no failures so update the experiment
            self.dest_sample_plates[dna_plate_id] = sorted(sample_plate_ids)
            for spid in sample_plate_ids:
                if spid not in self.plate_location_sample:  # should always be the case
                    self.plate_location_sample[spid] = {'purpose':'sample','source':'musterer','wells':set()}
                for record in sample_info:  # 'wellLocation', 'mouse', 'mouseBarcode', 'mouseId'
                    pos = util.unpadwell(record['wellLocation'])
                    self.plate_location_sample[spid]['wells'].add(pos)
                    self.plate_location_sample[spid][pos] = record['mouse'].copy()  # dict
                    self.plate_location_sample[spid][pos]['mouseId'] = record['mouseId']
                    self.plate_location_sample[spid][pos]['barcode'] = util.guard_mbc(record['mouseBarcode'])
                    self.plate_location_sample[spid][pos]['must_assays'] = record['mouse']['assays'].copy()
                    assays = [assay['assayName'] for assay in record['mouse']['assays']]
                    assayFamilies = set([a.split('_')[0] for a in assays])
                    self.plate_location_sample[spid][pos]['assays'] = assays.copy()
                    self.plate_location_sample[spid][pos]['assayFamilies'] = list(assayFamilies)
                
                self.log(f"Success: added sample plate {spid} with destination {dna_plate_id} ")

            #print(f"{self.name=} {self.plate_location_sample=} {self.sample_plates=}")  
        finally:
            self.save()
        return True

    def add_manifest(self, manifest_strm, default_manifest_type='c'):
        """ Because streamlit's uploader is derived from BytesIO we have a byte stream we have to deal with.
            The manifest is a CSV. Turn this into a list of dicts which can then handle sensibly.
            We also take an optional default_manifest_type ['c','r','m'] for custom, rodentity, or musterer respectively
            Returns True on success, False on failure
        """
        if self.locked:
            self.log('Error: Cannot add manifest while lock is active')
            return False
        try:
            self.log(f"Begin: add custom manifest")
        
            self.log(f"DEBUG: {self.name=} {manifest_strm=} {default_manifest_type=}")
            manifest_io = StringIO(manifest_strm.getvalue().decode("utf-8"))
        
            header = ''
            manifest_info = []
            in_csv = csv.reader(manifest_io, delimiter=',')
            for i, row in enumerate(in_csv):
                if 'barcode' in ''.join(row).lower():
                    header = row
                    continue  # header
                elif all([col.strip() for col in row]) == '': # blank
                    continue
                elif any(col.strip() == '' and j < 5 for j,col in enumerate(row)):
                    self.log(f"Warning: not adding blank manifest entry in row {i+1}: {row}")
                    continue
                else:
                    manifest_info.append({header[i]:r.strip() for i,r in enumerate(row) if i<len(header)})

            # gather up dest plates and sample plates so we can clear these if we already have them in the experiment
            dpids_spids = {}  # {dpid:set([spid,..])}
            # go through entries and guard barcodes as needed
            for i,entry in enumerate(manifest_info):
                for key in entry:
                    if key == 'Dest barcode':
                        manifest_info[i][key] = util.guard_pbc(entry[key], silent=True)
                        dpids_spids[manifest_info[i][key]] = set()
                    if key == 'Plate barcode':
                        manifest_info[i][key] = util.guard_pbc(entry[key], silent=True)
                    if key == 'Sample barcode':
                        if not util.is_guarded(entry[key]):
                            if default_manifest_type == 'c':
                                manifest_info[i][key] = util.guard_cbc(entry[key])
                            elif default_manifest_type == 'r':
                                manifest_info[i][key] = util.guard_rbc(entry[key])
                            elif default_manifest_type == 'm':
                                manifest_info[i][key] = util.guard_mbc(entry[key])
            for i,entry in enumerate(manifest_info):
                dpids_spids[entry['Dest barcode']].add(entry['Plate barcode'])

            # now clear any existing plateIds
            for dpid in dpids_spids:
                if dpid in self.plate_location_sample:
                    self.log(f"Clearing existing plate details for {dpid=}")
                    self.plate_location_sample.pop(dpid)
                for spid in dpids_spids[dpid]:
                    if spid in self.plate_location_sample:
                        self.log(f"Clearing existing plate details for {spid=}")
                        self.plate_location_sample.pop(spid)

            #print(f"{manifest_info=}")
            dest_pids_sample_pids = {}
            sample_pids = set()
            for entry in manifest_info:
                if entry['Dest barcode'] not in dest_pids_sample_pids:
                    dest_pids_sample_pids[entry['Dest barcode']] = set()
                dest_pids_sample_pids[entry['Dest barcode']].add(entry['Plate barcode'])
                sample_pids.add(entry['Plate barcode'])
        
            #self.dest_sample_plates = {}  # {dest_pid:[4 sample plate ids]}
            #self.plate_location_sample = {}  # pid:{well:sample_id}
            #self.sample_info = {}  # sample_id: {strain:str, assays: [], possible_gts: [], 
            #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] }

            # check for overused destination PIDs and clashes with sample plate barcodes. Note: we can't meaningfully check for duplicated destination PIDs
            dp_count = {dp:len(sps) for dp, sps in dest_pids_sample_pids.items()}     
            for dp in dp_count:
                if dp_count[dp] > 4:
                    self.log(f"Error: 384-well DNA destination plate {dp} used by too many sample plates: {dp_count[dp]}")
                    return False
                if dp in sample_pids:
                    self.log(f"Error: 384-well DNA destination plate barcode also used as a sample plate barcode {dp}")
                    return False
            

            # Warn for duplicate sample plate barcodes?
            for sp in sample_pids:
                dp_set = set()
                for dp in dest_pids_sample_pids:
                    for entry in manifest_info:
                        if entry['Plate barcode'] == sp and entry['Dest barcode'] == dp:
                            dp_set.add(dp)
                if len(dp_set) > 1:
                    self.log(f"Warning: Sample plate barcode {sp} used in multiple 384-well DNA destination plates: {dp_set}")
        
            # we need to collect destination, sample lists
            # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], gts: [],
            #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
            dp_samples = {}
            for i,entry in enumerate(manifest_info):
                if 'Sample No' in entry:
                    sample_number = entry['Sample No']
                else:
                    sample_number = i+1
                dest_pid = entry['Dest barcode']
                source_pid = entry['Plate barcode']
                well = util.unpadwell(entry['Well'])
                assays = [entry[key] for key in entry if 'assay' in key.lower() and entry[key].strip()!='']
                assayFamilies = set([a.split('_')[0] for a in assays])
                
                if dest_pid not in dp_samples:
                    dp_samples[dest_pid] = set()
                dp_samples[dest_pid].add(source_pid)
                if source_pid not in self.plate_location_sample:
                    self.plate_location_sample[source_pid] = {'purpose':'sample','source':'manifest', 'wells':set()}
                if well in self.plate_location_sample[source_pid] and self.plate_location_sample[source_pid][well] != {}:
                    self.log(f"Error: duplicate {well=} in {source_pid=}. Continuing...")

                self.plate_location_sample[source_pid]['wells'].add(well)
                self.plate_location_sample[source_pid][well] = {}
                self.plate_location_sample[source_pid][well]['barcode'] = entry['Sample barcode']
                self.plate_location_sample[source_pid][well]['assays'] = assays.copy()
                self.plate_location_sample[source_pid][well]['assayFamilies'] = list(assayFamilies)
                self.plate_location_sample[source_pid][well]['sex'] = ''
                self.plate_location_sample[source_pid][well]['strain'] = ''
                self.plate_location_sample[source_pid][well]['other_id'] = ''
                self.plate_location_sample[source_pid][well]['sampleNumber'] = str(sample_number)
                for key in entry:
                    if 'sex' in key.lower():
                        self.plate_location_sample[source_pid][well]['sex'] = entry[key]
                    elif 'strain' in key.lower():
                        self.plate_location_sample[source_pid][well]['strain'] = entry[key]
                    elif 'other_id' in key.lower():
                        self.plate_location_sample[source_pid][well]['other_id'] = entry[key]
                self.plate_location_sample[source_pid][well]['gts'] = []
                self.plate_location_sample[source_pid][well]['possible_gts'] = []
                self.plate_location_sample[source_pid][well]['parents'] = []      

            print(f"\n{dp_samples=}")
            for dp in dp_samples:
                self.dest_sample_plates[dp] = list(dp_samples[dp])
                for sample_pid in dp_samples[dp]:
                    self.log(f"SUCCESS: added sample plate {sample_pid} with destination {dp}")
                    print('')
                    print(f"{sample_pid=} {self.plate_location_sample[sample_pid]=}")
        finally:
            self.save()
        return True

    def get_musterer_pids(self):
        """ return [pids] with samples that are sourced from Musterer """
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
        """ return a dictionary summarising the contents of plate sets """
        plate_set_summary = []
        total_wells = 0
        total_unique_samples = set()
        total_unique_assays = set()
        total_well_counts = {'c':0,'m':0,'r':0}
        for dna_pid in self.dest_sample_plates:
            d = {'DNA PID': util.unguard_pbc(dna_pid)}
            d['Sample PID1'] = ''  # we need to define a regular structure for Pandas dataframe, I know this is ugly
            d['Sample PID2'] = ''
            d['Sample PID3'] = ''
            d['Sample PID4'] = ''
         
            #d['Well counts'] = {'c':0,'m':0,'r':0}
            d['Custom samples'] = 0
            d['Rodentity samples'] = 0
            #d['All primers required'] = 0
            #d['Unique samples'] = set()
            #d['Unique assays'] = set()
            
            for i, sample_pid in enumerate(sorted(self.dest_sample_plates[dna_pid])):
                d['Sample PID'+str(i+1)] = util.unguard_pbc(sample_pid)
                for info in (self.plate_location_sample[sample_pid][well] for well in self.plate_location_sample[sample_pid]['wells']):
                    samp_barcode = info['barcode']     

                    if util.is_guarded_cbc(samp_barcode):
                        d['Custom samples'] += 1
                        total_well_counts['c'] += 1
                    #elif util.is_guarded_mbc(samp_barcode):
                    #    d[' counts']['m'] += 1
                    #    total_well_counts['m'] += 1  
                    elif util.is_guarded_rbc(samp_barcode):
                        d['Rodentity samples'] += 1
                        total_well_counts['r'] += 1   
                    #for assay in info['assays']:
                    #    d['Primers required'].add(assay)
                    #    total_unique_assays.add(assay)
                    
                    #d['Unique samples'].add(info['barcode'])
                    total_unique_samples.add(info['barcode'])
                    
            #d['Unique samples'] = len(d['Unique samples'])
            #d['Unique assays'] = len(d['Unique assays'])
            #total_wells += sum(k[v] for k,v in d['Well counts'].items())
            plate_set_summary.append(d)
        plate_set_summary.append({'DNA PID':'Total','Sample PID1':'','Sample PID2':'','Sample PID3':'',
                #'Sample PID4':'','Well count':total_wells,'Unique samples': len(total_unique_samples),
                'Sample PID4':'','Custom samples': total_well_counts['c'], 'Rodentity samples': total_well_counts['r']})
                # 'Total assays': len(total_unique_assays)})
        #print(f"{plate_set_summary=}")
        return plate_set_summary

    def inputs_as_dataframe(self):
        """ return the experiment contents as a pandas dataframe """
        inputs_dict = self.summarise_inputs()
        if not inputs_dict:
            return None
        # need to reshape this list of dicts as dict of lists
        inputs = {key: [i[key] for i in inputs_dict] for key in inputs_dict[0]}
        return pd.DataFrame(inputs)

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

    def get_assay_usage(self, dna_plate_list=[], included_guards=util.GUARD_TYPES):
        """ We will want ways to get assay usage from various subsets, but for now do it for everything that 
            has a Nimbus destination plate.
        """
        assay_usage = {}
        for dest in self.dest_sample_plates:
            if dna_plate_list and dest not in dna_plate_list:
                continue
            for sample_pid in self.dest_sample_plates[dest]:
                plate = self.get_plate(sample_pid)
                assay_usage.update(util.calc_plate_assay_usage(plate,included_guards=included_guards))
        return assay_usage

    def get_volumes_required(self, assay_usage=None, dna_plate_list=[]):
        """
            Standard usage is to call get_assay_usage() first and pass it to this.
        """
        if not assay_usage:
            assay_usage = self.get_assay_usage(dna_plate_list=dna_plate_list)
        #print (f'{assay_usage=}', assay_usage.values(), file=sys.stderr)
        reactions = sum([v for v in assay_usage.values()])

        primer_taq_vol = reactions * self.transfer_volumes['PRIMER_TAQ_VOL']
        primer_water_vol = reactions * self.transfer_volumes['PRIMER_WATER_VOL']
        index_water_vol = reactions * self.transfer_volumes['INDEX_WATER_VOL']
        index_taq_vol = reactions * self.transfer_volumes['INDEX_TAQ_VOL']

        primer_vols = {a:assay_usage[a]*self.transfer_volumes['PRIMER_VOL'] for a in assay_usage}
        return primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol
        
        
    def get_index_remaining_available_volume(self, assay_usage=None):
        """
        Returns the barcode pairs remaining, max available barcode pairs, max barcode pairs allowed by volume
        """
        index_pids = []
        for pid in self.plate_location_sample:
            if self.plate_location_sample[pid]['purpose'] == 'index':
                index_pids.append(pid)

        print(f'get_index_remaining_available_volume() {index_pids=}', file=sys.stderr)
        if not index_pids:
            return 0, 0, False

        if not assay_usage:
            assay_usage = self.get_assay_usage()
        reactions = sum([v for v in assay_usage.values()])
        
        fwd_barcode_vols = {}  # name=[vol, vol, ...]
        rev_barcode_vols = {}
        
        for idx_pid in index_pids:
            idx_plate = self.get_plate(idx_pid)
            #print(f'get_index_remaining_available_volume() {idx_plate=}', file=sys.stderr)
            
            for well in idx_plate['wells']:
                if 'idt_name' not in idx_plate[well]:
                    continue
                name = idx_plate[well]['idt_name']
                if 'volume' in idx_plate[well]:
                    print(f"get_index_remaining_available_volume() {idx_plate[well]['volume']=}")
                if 'i7F' in name:
                    if name not in fwd_barcode_vols:
                        fwd_barcode_vols[name] = []
                    if 'volume' in idx_plate[well]:
                        fwd_barcode_vols[name].append(max(idx_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0))
                elif 'i5R' in name:
                    if name not in rev_barcode_vols:
                        rev_barcode_vols[name] = []
                    if 'volume' in idx_plate[well]:
                        rev_barcode_vols[name].append(max(idx_plate[well]['volume'] - util.DEAD_VOLS[util.PLATE_TYPES['Echo384']],0))
                else:
                    self.log('Unexpected index name:' + name, level='Warning')
        max_i7F = len(fwd_barcode_vols)
        max_i5R = len(rev_barcode_vols)
        reaction_vol_capacity = 0
        for name in fwd_barcode_vols:
            # get the number of possible reactions and keep the lower number from possible reactions or possible reaction partners
            max_reactions = sum([floor((vol-util.DEAD_VOLS[util.PLATE_TYPES['Echo384']])/self.transfer_volumes['INDEX_VOL']) \
                    for vol in fwd_barcode_vols[name]])
            reaction_vol_capacity += min(max_reactions, max_i5R)
        max_barcode_pairs = max_i7F * max_i5R

        return max_barcode_pairs-reactions, max_barcode_pairs, reaction_vol_capacity


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
        print(f"In remove_entries. {selected_rows=}")
        for row in selected_rows:
            if row['DNA PID'] not in self.dest_sample_plates:
                self.log(f"Error: {row['DNA PID']=} doesn't actually exist in the experiment!")
                continue
            dest_pid = row['DNA PID']
            sample_pids = self.dest_sample_plates.pop(row['DNA PID'])
            delete_pids = sample_pids.append(dest_pid)
            self.delete_plates(delete_pids)
        self.save()
        return True

    def add_references(self, uploaded_references):
        """
        read in reference (target) IDs and sequences
        """
        if self.locked:
            self.log('Error: Cannot add reference sequences while lock is active')
            self.save()
            return False
        partial_fail = False
        for uploaded_reference in uploaded_references:
            if uploaded_reference.name in self.reference_sequences:
                self.log(f"Warning: Duplicate reference file name: {uploaded_reference.name=} in . Overwriting")
            try:
                with StringIO(uploaded_reference.getvalue().decode()) as sio:
                    self.reference_sequences[uploaded_reference.name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
            except UnicodeDecodeError:
                try:
                    with StringIO(uploaded_reference.getvalue().decode('utf-8')) as sio:
                        self.reference_sequences[uploaded_reference.name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
                except UnicodeDecodeError:
                    try:
                        with StringIO(uploaded_reference.getvalue().decode('latin-1')) as sio:
                            self.reference_sequences[uploaded_reference.name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
                    except UnicodeDecodeError:
                        try:
                            with StringIO(uploaded_reference.getvalue().decode('ISO-8859')) as sio:
                                self.reference_sequences[uploaded_reference.name] = {str(seqitr.id):str(seqitr.seq) for seqitr in SeqIO.parse(sio, "fasta")}
                        except Exception as exc:
                            self.log(f"Error: Couldn't parse reference file {uploaded_reference.name} {exc}")
                            self.save()
                            partial_fail = True
                            continue                
            self.log(f'Success: uploaded {len(self.reference_sequences[uploaded_reference.name])} reference sequences from {uploaded_reference.name}')
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
        self.add_reproducible_steps(transactions)
        self.log(f'Success: created reference sequences file {target_fn} containing {counter} sequences')
        self.save()
        return True
        

    def add_assaylists(self, uploaded_assaylists):
        """ mapping of assay to primer name - may not be required in future in this capacity.
        What we really need/want is a list of allowed assays
        """
        if self.locked:
            self.log('Error: cannot add assay list while lock is active.')
            return False
        for uploaded_assaylist in uploaded_assaylists:
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
        self.log(f"Success: added primer-assay lists: {', '.join([ual.name for ual in uploaded_assaylists])}")
        self.save()
        return True

    def add_primer_layouts(self, uploaded_primer_plates):
        """ add primer plate definition with well and name columnes """
        if self.locked:
            self.log('Error: cannot add primer plates while lock is active.')
            return False
        for uploaded_primer_plate in uploaded_primer_plates:
            PID = util.guard_pbc(uploaded_primer_plate.name.split('_')[0], silent=True)  # assumes first field is PID
            if PID in self.plate_location_sample:
                if self.plate_location_sample[PID]['purpose'] == 'primer':
                    self.log(f"Info: Primer file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Primer file PID: {uploaded_primer_plate=} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return
            else:  # create new plate entry and set purpose
                self.plate_location_sample[PID] = {'purpose':'primer', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
                self.log(f"Info: Creating new primer plate record for {PID}")
            # load data into plate
            data = StringIO(uploaded_primer_plate.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    continue  # header
                if row == '' or row[0] == '' or row[1] == '':
                    continue
                well = util.unpadwell(row[0])
                if well not in self.plate_location_sample[PID]:
                    self.plate_location_sample[PID]['wells'].add(well)
                    self.plate_location_sample[PID][well] = {}
                print(row[0], row[1])
                self.plate_location_sample[PID][well]['primer'] = row[1]
        self.log(f"Success: added primer layouts from {', '.join([upl.name for upl in uploaded_primer_plates])}")
        self.save()
        return True

    def add_primer_volumes(self, uploaded_primer_volumes):
        """ add primer plate volumes with well and volume columns """
        if self.locked:
            self.log('Error: cannot add primer plate volumes while lock is active.')
            return False
        for uploaded_primer_volume in uploaded_primer_volumes:
            PID = util.guard_pbc(uploaded_primer_volume.name.split('_')[0], silent=True)  # assumes first field is PID
            if PID in self.plate_location_sample:
                if self.plate_location_sample[PID]['purpose'] == 'primer':
                    self.log(f"Info: Primer file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Primer file PID: {uploaded_primer_volume=} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return
            else:  # create new plate entry and set purpose
                self.plate_location_sample[PID] = {'purpose':'primer', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
                self.log(f"Info: Creating new primer plate record for {PID}")
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
                        if well not in self.plate_location_sample[PID]:
                            self.plate_location_sample[PID]['wells'].add(well)
                            self.plate_location_sample[PID][well] = {}
                        self.plate_location_sample[PID][well]['volume'] = float(col)*1000
                else:
                    well = util.unpadwell(row[0])
                    if well not in self.plate_location_sample[PID]:
                        self.plate_location_sample[PID]['wells'].add(well)
                        self.plate_location_sample[PID][well] = {}
                    self.plate_location_sample[PID][well]['volume'] = float(row[1])*1000
        self.log(f"Success: added primer volumes from {', '.join([upv.name for upv in uploaded_primer_volumes])}")
        self.save()
        return True

    def add_amplicon_manifests(self, uploaded_amplicon_manifests):
        """
        Amplicon plate layouts (containing pre-amplified sequences) are added here via manifest files.
        Only column layout (comma separate) is supported. col1: plate barcode; col2: well position; col3: sample barcode; col4 (optional): volume.
        The first row must be a header row: plate, well, sample, (volume in uL, if provided).
        """
        if self.locked:
            self.log('Error: cannot add amplicon plate layouts while lock is active.')
            return False

        for uploaded_amplicon_manifest in uploaded_amplicon_manifests:
            data = StringIO(uploaded_amplicon_manifest.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    continue  # header
                cols = [c.strip() for c in row]
                PID = util.guard_pbc(cols[0], silent=True)
                if PID in self.plate_location_sample and self.plate_location_sample[PID]['purpose'] != 'amplicon':
                    self.log(f"Error: Amplicon plate barcode: {PID=} matches "+\
                            f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return False
                if PID not in self.plate_location_sample:
                    self.plate_location_sample[PID] = {'purpose':'amplicon', 'source':'user', 'wells':set(), 
                            'plate_type':'384PP_AQ_BP', 'barcode':PID}
                    self.log(f"Info: Creating new amplicon plate record for {PID}")
                well = util.unpadwell(cols[1])
                try:
                    sample = util.guard_abc(cols[2], silent=True)
                except Exception as exc:
                    self.log(f"Error: Amplicon sample barcode could not be guarded {exc=}")
                    self.delete_plates([PID])
                    return False
                if len(cols) == 4:
                    try:
                        vol = float(cols[3])*1000  # save as nanolitres
                    except ValueError:
                        self.log(f"Warning: could not interpret volume as numeric. Ignoring volume {cols[3]=}")
                        vol = None

                self.plate_location_sample[PID]['wells'].add(well)
                self.plate_location_sample[PID][well] = {'barcode':sample}
                if vol:
                    self.plate_location_sample[PID][well]['volume'] = vol
            
        self.log(f"Success: added amplicon plate info from {', '.join([uam.name for uam in uploaded_amplicon_manifests])}")
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
                        #print(outline, file=fout)
        except Exception as exc:
            self.log(f'Failure: could not write primer survey {exc}')
            self.save()
            return False
        self.add_reproducible_steps(transactions)
        self.log(f"Success: written Echo primer survey to {primer_survey_fn}")
        return True


    def check_ready_echo1(self, dna_plates, pcr_plates, taq_water_plates):
        """
        TODO: Should check that everything required to successfully generate PCR1 picklists is available
        """
        for pid in pcr_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'pcr':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
            self.log(f'Warning: existing entry for PCR plate {pid}. This will be overwritten!')

        for pid in dna_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'dna':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
        return True


    def check_ready_echo2(self, pcr_plates, taq_water_plates, included_index_plates):
        """
        TODO: Should check the everything required to successfully generate PCR2 picklists is available
        """
        return True


    def generate_echo_PCR1_picklists(self, dna_plates, pcr_plates, taq_water_plates):
        """
        Calls echo_primer.generate_echo_PCR1_picklist() to do the work, needs a set of accepted DNA_plates,
        the final assay list, and a set of destination PCR plate barcodes, taq+water plates, primer plates, primer volumes.
        Returns True on success
        TODO: We need to reduce available taq and water during this operation.
        """
        #print(f"Experiment.generate_echo_PCR1 {dna_plates=} {pcr_plates=} {taq_water_plates=}")
        
        success = self.generate_echo_primer_survey()
        if not success:
            self.log('Failure: failed to generate primer survey file')
            return False
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


    def generate_echo_PCR2_picklist_interface(self, pcr_plates, index_plates, taq_water_plates, amplicon_plates=None):
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
            
        success = self.generate_echo_index_survey()
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

    def get_index_pids(self):
        """
        Used by Indexing stage
        """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'index']

    def add_index_layouts(self, uploaded_index_plates):
        """ add index plates with well and index columns """
        if self.locked:
            self.log('Error: cannot add index plates while lock is active.')
            return False
        for uploaded_index_plate in uploaded_index_plates:
            PID = util.guard_pbc(uploaded_index_plate.name.split('_')[0], silent=True)  # assumes first field is PID
            if PID in self.plate_location_sample:
                if self.plate_location_sample[PID]['purpose'] == 'index':
                    self.log(f"Info: Index plate {PID} already exists, adding data")
                else:
                    self.log(f"Error: Index plate PID: {uploaded_index_plate=} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return
            else:  # create new plate entry and set purpose
                self.plate_location_sample[PID] = {'purpose': 'index', 'source':'user', 'wells':set(), 
                        'plate_type':util.PLATE_TYPES['Echo384'], 'barcode':PID}
                self.log(f"Info: Creating new index plate record for {PID}")
            # load data into plate
            data = StringIO(uploaded_index_plate.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0 or row == '':
                    continue  # header or blank
                well = util.unpadwell(row[0])
                if well not in self.plate_location_sample[PID]:
                    self.plate_location_sample[PID]['wells'].add(well)
                    self.plate_location_sample[PID][well] = {}
                    idt_name = row[1]
                    index = row[2]
                    bc_name = row[3]
                    oligo = row[4]
                    self.plate_location_sample[PID][well]['idt_name'] = idt_name
                    self.plate_location_sample[PID][well]['index'] = index
                    self.plate_location_sample[PID][well]['bc_name'] = bc_name
                    self.plate_location_sample[PID][well]['oligo'] = oligo
                    
        self.log(f"Success: added index plate layouts from {', '.join([uip.name for uip in uploaded_index_plates])}")        
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
            PID = util.guard_pbc(uploaded_index_volume.name.split('_')[0], silent=True)  # assumes first field is PID
            if PID not in self.plate_location_sample:
                self.log(f"Info: Adding index plate {PID}")
                self.plate_location_sample[PID] = {'purpose':'index', 'source':'user', 'wells':set(), 
                        'plate_type':util.PLATE_TYPES['Echo384'], 'barcode':PID}
            else:
                if self.plate_location_sample[PID]['purpose'] != 'index':
                    self.log(f"Error: {PID} plate purpose is "+\
                            f"{self.plate_location_sample[PID]['purpose']}, expected 'index'")
                    return False
                self.log(f"Info: {PID} exists, appending index volumes")
        
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
                            self.plate_location_sample[PID]['wells'].add(well)
                            if well not in self.plate_location_sample[PID]:
                                self.plate_location_sample[PID][well] = {}
                            self.plate_location_sample[PID][well]['volume'] = float(v)*1000
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
                        self.plate_location_sample[PID]['wells'].add(well)
                        if well not in self.plate_location_sample[PID]:
                            self.plate_location_sample[PID][well] = {}
                        self.plate_location_sample[PID][well]['volume'] = int(float[5])*1000
                
        self.log(f"Success: added index plate volumes from {', '.join([uiv.name for uiv in uploaded_index_volumes])}")
        self.save()
        return True

    def generate_echo_index_survey(self, index_survey_filename='index-svy.csv'):
        """ 
        Generate an index survey file for use by the Echo. Replaces echovolume.py 
        Not strictly necessary, but the old code reads this file in making the picklists.
        Called internally by generate_echo_PCR2_picklist_interface()
        """
        if self.locked:
            self.log('Error: cannot generate index survey while lock is active.')
            return False
        index_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'index']
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
            self.log(f"Error: adding taq+water plate barcode failed {plate_barcodes=}")
            return False
        self.save()
        return True

    def get_miseq_samplesheets(self):
        """ return the MiSeq-XXX.csv samplesheet, if it exists """
        miseq_fps = Path(self.get_exp_dir()).glob('MiSeq-*.csv')
        return miseq_fps

     
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
