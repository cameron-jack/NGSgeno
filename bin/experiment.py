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

from pandas.core.base import NoNewAttributesMixin
import jsonpickle
from itertools import combinations
from io import StringIO, BytesIO
import datetime
from shutil import copyfileobj
import json
from Bio import SeqIO
from copy import deepcopy
import weakref
from math import ceil, floor
from pathlib import Path
import inspect  
from collections import defaultdict

import pandas as pd

try:
    import bin.db_io as db_io
except ModuleNotFoundError:
    import db_io
try:
    import bin.generate as generate
except ModuleNotFoundError:
    import generate
try:
    import bin.parse as parse
except ModuleNotFoundError:
    import parse
try:
    import bin.util as util
except ModuleNotFoundError:
    import util
try:
    import bin.transaction as transaction
except ModuleNotFoundError:
    import transaction

from stutil import m 

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
        # Standardised: "assays": [str], rodentity allele[assay[name]]
        # The specific import function should be responsisble for creating the "assays" field
        ###
        self.denied_assays = []  # file to plate mapping
        self.denied_primers = []  # plate to file mapping
        self.assay_synonyms = {}  # {source:{reference:alternative}}
        self.primer_assayfam = {}  # {mapping of primer to assayfam} for reverse lookups
        self.assayfam_primers = {}  # {mapping of assay families to list of primers}
        self.assay_assayfam = {}  # mapping of assay to assay family
        self.assayfam_assays = {}  # mapping of assay family to list of assays for reverse lookups
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
        self.uploaded_files = {}  # {'_upload_queue':{fp:tuple}} defers all uploads, {'_upload_pending':{fp:tuple}} for pending uploads
        self.reproducible_steps = []  # a strict subset of the log that only includes working instructions
        self.pending_steps = set() # reproducible steps (files) that are awaiting user approval to replace existing steps
        self.pending_uploads = set() # DEPRECATED - replaced with uploaded_files['_upload_pending'] = {}
        self.log_entries = []  # use self.log(message, level='') to add to this
        self.log_time = datetime.datetime.now()  # save the log after a certain time has elapsed, to avoid too much IO
        self.transfer_volumes = {'DNA_VOL':util.DNA_VOL, 'PRIMER_VOL':util.PRIMER_VOL, 'PRIMER_TAQ_VOL':util.PRIMER_TAQ_VOL,
                'PRIMER_WATER_VOL':util.PRIMER_WATER_VOL, 'INDEX_VOL':util.INDEX_VOL, 'INDEX_TAQ_VOL':util.INDEX_TAQ_VOL,
                'INDEX_WATER_VOL':util.INDEX_WATER_VOL}  # load the defaults from util.py
        self._finalizer = weakref.finalize(self, self.save)
        

    def __setstate__(self, state):
        for k, v in state.items():
            setattr(self, k, v)

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __repr__(self) -> str:
        return str(self.__dict__)

    def lock(self):
        self.log('Info: Locking/unlocking of experiment is currently disabled')
        #self.log('Info: Locking experiment. No modification allowed to plates while lock remains')
        #self.locked = True

    def unlock(self):
        self.log('Info: Locking/unlocking of experiment is currently disabled')
        #self.log('Info: Unlocking experiment. Modification is now possible. There should be good reason for this!')
        #self.locked = False


    ### functions for returning locally held file paths
    
    def get_exp_dn(self, subdir=None, caller_id=None):
        """ 
        Return the experiment directory name
        If subdir is supplied return the path to this:    
        'raw' - fastq files should be stored
        'cleaned' - fastq files should be stored
        'merged' - fastq files should be stored
        'uploads' - uploaded file should be stored
        """
        dirpath = 'run_' + self.name 
        if not subdir:
            return dirpath
        allowed_subdirs = ['raw','cleaned','merged','uploads']
        if subdir in allowed_subdirs:
            dp2 = Path(dirpath)/subdir
            if not dp2.exists():
                dp2.mkdir()
            else:
                if not dp2.is_dir():
                    m(f'failed to make directory {str(dp2)}, which already exists as a file', level='critical', caller_id=caller_id)
            dirname = str(dp2)
            return str(dirname)
        else:                                                       
            m(f'{subdir=} not in {allowed_subdirs=}', level='critical', caller_id=caller_id)


    def get_exp_fn(self, filename, subdir=None, trans=False, caller_id=None):
        """ 
        Return the expected experiment path to filename as a string
        If trans is True, check for an existing match and append "_pending" if required
        If subdir is in allowed_subdirs, then include this in the returned filepath
        """
        dirname = self.get_exp_dn(subdir)
        if dirname not in filename:
            fp = Path(os.path.join(dirname, filename))
        else:
            fp = Path(filename)
        if trans:
            fp = transaction.transact(str(fp)) 
        return str(fp)

    
    def get_raw_fastq_pairs(self, caller_id=None):
        """ return a sorted list of tuple(R1_path, R2_path) to raw FASTQ files """
        valid_pairs = []
        rdp = Path(self.get_exp_dn(subdir='raw'))
        r1s = [rdp/Path(f) for f in os.listdir(rdp) if f.endswith('.fastq.gz') and '_R1_' in f]
        for r1 in r1s:
            r2 = Path(str(r1).replace('_R1_001','_R2_001'))
            if r1.is_file() and r2.is_file():
                valid_pairs.append((r1,r2))
            else:
                if not r1.is_file():
                    m(f'{r1} expected raw FASTQ file does not exist', level='warning', caller_id=caller_id)
                elif not r2.is_file():
                    m(f'{r2} expected raw FASTQ file does not exist', level='warning', caller_id=caller_id)
        return sorted(valid_pairs)
    
    ### functions for managing data file records

    def add_file_record(self, file_name, PIDs=None, purpose=None, caller_id=None):
        """
        Add an uploaded file to the experiment with associated plates and purpose
        """
        if file_name in self.uploaded_files:
            m(f'file name {file_name} has already been uploaded, overwriting', level='warning', caller_id=caller_id)

        if not PIDs:
            PIDs = []

        file_md5 = util.get_md5(file_name)
        ft = os.path.getmtime(file_name)
        ft = datetime.datetime.fromtimestamp(ft).strftime('%Y/%m/%d %H:%M:%S')
        self.uploaded_files[file_name] = {'plates': [util.guard_pbc(PID, silent=True) for PID in PIDs], 
                'purpose': purpose, 'md5':file_md5, 'date modified':ft}
        return True


    def mod_file_record(self, existing_name, new_name=None, new_purpose=None,
            extra_PIDs=None, PID_name_updates=None, caller_id=None):
        """
        Modify an existing file record (self.uploaded_files)
        Mostly used to change a pending file to normal file path, or to modify the list of plates
        new_name (str)
        new_purpose (str)
        extra_pids ([str])
        PID_name_updates ([(str,str)]) - list of (from,to) tuples.
        """
        if existing_name not in self.uploaded_files:
            m(f'{existing_name} does not exist in uploaded file records', level='critical', caller_id=caller_id)
            return False

        if new_purpose:
            self.uploaded_files[existing_name]['purpose'] = new_purpose

        if extra_PIDs:
            for ep in extra_PIDs:
                self.uploaded_files[existing_name]['plates'].append(ep)

        if PID_name_updates:
            for pnu in PID_name_updates:
                existing_pid, new_pid = pnu
                if existing_pid not in self.uploaded_files[existing_name]['plates']:
                    m(f'{existing_pid} not present in uploaded_files', level='critical', caller_id=caller_id)
                    return False
                self.uploaded_files[existing_name]['plates'].remove(existing_pid)
                self.uploaded_files[existing_name]['plates'].append(new_pid)
                                                                                    
        if new_name:
            
            self.uploaded_files[new_name] = self.uploaded_files[existing_name].copy()
            del self.uploaded_files[existing_name]
        
        return True
    

    def del_file_record(self, file_name, soft=True, caller_id=None):
        """
        Remove a file record, and its plate definitions (to deleted_plates) if required
        If 'soft' then rename the file to Dyyyymmddhhmmss_filename
        If not 'rename' then delete the file permanently
        """
        if file_name in self.uploaded_files:
            # get rid of actual file
            if Path(file_name).exists():
                util.delete_file(file_name, soft=soft, caller_id=caller_id)
            
            # clear any plates first
            for pid in self.uploaded_files[file_name].get('plates', []):
                if pid in self.plate_location_sample:
                    filepath_purpose = self.plate_location_sample[pid].get('filepath_purpose', None)
                    if isinstance(filepath_purpose, dict) and file_name not in filepath_purpose:
                        continue
                    else:
                        success = self.delete_plate(pid)
                        if success:
                            m(f'Removed plate entry {pid} from deleted file {file_name}',
                                    level='info', dest=('noGUI',), caller_id=caller_id)
                        else:
                            m(f'Failed to remove plate entry {pid} from deleted file {file_name}', 
                                    level='error', dest=('noGUI',), caller_id=caller_id)
                        
            # now remove the actual entry
            del self.uploaded_files[file_name]
            m(f'removed file record of {file_name}, attempting to remove file', level='success', dest=('noGUI',), caller_id=caller_id)
            return True
        
        if Path(file_name).exists():
            #m(f'No file record exists, but file found. Attempting to remove file', 
            #        level='warning', caller_id=caller_id)
            util.delete_file(file_name, soft=soft, caller_id=caller_id)
            return True
        m(f'{file_name} not present in file records, and file does not exist', level='warning', caller_id=caller_id)
        return False
    

    ### Plate related operations

    def add_plate(self, pid, dependent_files, caller_id=None):
        """
        Standardise the interface for a plate record. Useful due to the complexity of plates
        All plate records go in:
                self.plate_location_sample = {}  # pid:{well:sample_dict}  # use this for everything!
        It's important that we record a bunch of plate properties:
            'barcode' (str) plate ID
            'plate_type' an entry from util.PLATE_TYPE
            'wells' set()  # either occupied wells, or the set of all possible wells
            'timestamp'  # when this entity came into existence
            'purpose' [''pcr',')
            'source' ['user','custom','rodentity']
            'dependent_files' set(str)  # the collection of files created from this plate record, can't be set at creation
            [well] = {} the actual data about what's in a well. Don't set this here, but know that it exists!
            purpose ['sample','dna','pcr','primer','index',''] 
                                             'wells':set(), 
                                             'source':'user', 
                                            'plate_type':util.PLATE_TYPES['PCR384'], 
                                            'barcode':p}
        """
        pass

    def get_plate(self, pid, caller_id=None):
        """ 
        Return the actual plate as a nested dictionary object.
        This is required because plates may have had things added or removed from them
        This function uses transactions to build a record of an up-to-date plate dictionary
        """
        return transaction.get_plate(self, pid)
    

    def build_dna_plate_entry(self, sample_plate_ids, dna_plate_id, source=None, caller_id=None):
        """
        Replaces the rodentity- and custom-specific code for combining 96-well sample plates into
        384-well DNA plates.
        """
        if self.locked:
            m('Cannot add DNA plate set while lock is turned on', level='failure', caller_id=caller_id)
            return False
        
        dna_plate_id = util.guard_pbc(dna_plate_id, silent=True)

        if dna_plate_id in self.plate_location_sample:
            m(f'plate {dna_plate_id} already exists! Please delete this plate before trying again', level='failure', caller_id=caller_id)
            return False

        sample_plate_ids = sorted([util.guard_pbc(spid, silent=True) for spid in sample_plate_ids if spid])

        if source == 'rodentity':
            m(f'combining Rodentity plate set into 384-well DNA plate {dna_plate_id}', level='info', caller_id=caller_id)
            for spid in sample_plate_ids:
                purpose = self.plate_location_sample[spid]['purpose']
                source = self.plate_location_sample[spid]['source']
                if purpose != 'sample' or source != 'rodentity':
                    m(f'cannot combine {spid} with {purpose=} and {source=}', level='error', caller_id=caller_id)
                    return False
        elif source == 'custom':
            m(f'combining custom plate set into 384-well DNA plate {dna_plate_id}', level='begin', caller_id=caller_id)
            for spid in sample_plate_ids:
                purpose = self.plate_location_sample[spid]['purpose']
                source = self.plate_location_sample[spid]['source']
                if purpose != 'sample' or source != 'custom':
                    m(f'cannot combine {spid} with {purpose=} and {source=}', level='error', caller_id=caller_id)
                    return False
        else:
            m('must choose "rodentity" or "custom" as input source', level='failure', caller_id=caller_id)
            return False

        self.dest_sample_plates[dna_plate_id] = sample_plate_ids
        self.plate_location_sample[dna_plate_id] = {'purpose':'dna', 'source':','.join(sample_plate_ids), 
                'wells':set(), 'plate_type':util.PLATE_TYPES['Echo384']}
        if source == 'rodentity':
            self.unassigned_plates[1] = ''
            self.unassigned_plates[2] = ''
            self.unassigned_plates[3] = ''
            self.unassigned_plates[4] = ''
        elif source == 'custom':
            self.unassigned_plates['custom'] = {}
        return True


    def get_musterer_pids(self, caller_id=None):
        """ return [pids] with samples that are sourced from Musterer. DEPRECATED """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        musterer_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if util.is_guarded_mbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    musterer_pids.add(pid)
        return musterer_pids


    def get_rodentity_pids(self, caller_id=None):
        """ return [pids] with samples that are sourced from Rodentity """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        rodentity_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if util.is_guarded_rbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    rodentity_pids.add(pid)
        return rodentity_pids

    
    def get_custom_pids(self, caller_id=None):
        """ return [pids] with samples that are sourced as custom - useful for mixed content plates """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        custom_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if util.is_guarded_cbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    custom_pids.add(pid)
        return custom_pids
        
    
    def summarise_inputs(self, caller_id=None):
        """ return a list of fields and a list of headers, summarising the contents of sample, amplicon, and DNA plate sets """
        plate_set_summary = []
        plate_set_headers = ['DNA/amplicon PID', 'Sample PID1', 'Sample PID2', 'Sample PID3', 'Sample PID4', 'Custom samples', 'Rodentity samples']
        total_wells = 0
        total_unique_samples = set()
        total_unique_assays = set()
        total_well_counts = {'a':0,'c':0,'m':0,'r':0}
        for amp_pid in self.get_amplicon_pids():
            amp_wells = 0
            for well in self.plate_location_sample[amp_pid]['wells']:
                if 'barcode' in self.plate_location_sample[amp_pid][well]:
                    barcode = self.plate_location_sample[amp_pid][well]['barcode']
                    if util.is_guarded_abc(barcode):
                        amp_wells += 1
                        total_well_counts['a'] += 1
                    elif util.is_guarded_cbc(barcode):  # these should be guarded with 'a', but we'll note them anyway
                        amp_wells += 1
                        total_well_counts['a'] += 1
            plate_set_summary.append([util.unguard_pbc(amp_pid, silent=True), '', '', '', '', amp_wells, 0])
            
        for dna_pid in self.dest_sample_plates:
            plate_set_details = [util.unguard_pbc(dna_pid, silent=True)]
            custom_wells = 0
            rodentity_wells = 0
            pid_count = 0
            # check the sample plates still exist
            new_dest_scheme = []
            for i, spid in enumerate(self.dest_sample_plates[dna_pid]):
                if spid in self.plate_location_sample:
                    new_dest_scheme.append(spid)
                # drop any deleted plates
            self.dest_sample_plates[dna_pid] = new_dest_scheme
            #print(f'{dna_pid=} {self.dest_sample_plates[dna_pid]=}')

            for i,sample_pid in enumerate(sorted(self.dest_sample_plates[dna_pid])):
                plate_set_details.append(util.unguard_pbc(sample_pid, silent=True))
                for well in self.plate_location_sample[sample_pid]['wells']:
                    if 'barcode' in self.plate_location_sample[sample_pid][well]:
                        barcode = self.plate_location_sample[sample_pid][well]['barcode']
                        if util.is_guarded_cbc(barcode):
                            custom_wells += 1
                            total_well_counts['c'] += 1 
                        elif util.is_guarded_rbc(barcode):
                            rodentity_wells += 1
                            total_well_counts['r'] += 1   
                        total_unique_samples.add(barcode)
                #print(f'{dna_pid=} {sample_pid=} {custom_wells=} {rodentity_wells=}')
                pid_count += 1
            for j in range(4-pid_count):
                plate_set_details.append('')

            plate_set_details.append(str(custom_wells))
            plate_set_details.append(str(rodentity_wells))
            # print(f'{plate_set_details=}', file=sys.stderr)
            plate_set_summary.append(plate_set_details)

        viewed_plates = set([p[0] for p in plate_set_summary])
        # now get any individual 384-well DNA plates that were loaded separately
        for dna_pid in self.get_dna_pids(plate_type='Echo384'):
            dpid = util.unguard_pbc(dna_pid, silent=True)
            if dpid in viewed_plates:
                continue
            custom_wells = 0
            rodentity_wells = 0
            if dna_pid not in self.dest_sample_plates:
                for well in self.plate_location_sample[dna_pid]['wells']:
                    if 'barcode' in self.plate_location_sample[dna_pid][well]:
                        barcode = self.plate_location_sample[dna_pid][well]['barcode']
                        if util.is_guarded_cbc(barcode):
                            custom_wells += 1
                            total_well_counts['c'] += 1 
                        elif util.is_guarded_rbc(barcode):
                            rodentity_wells += 1
                            total_well_counts['r'] += 1   
                    
                        total_unique_samples.add(barcode)
            plate_set_summary.append([dpid,'','','','',str(custom_wells),str(rodentity_wells)])

        if len(plate_set_summary) > 0:
            plate_set_summary.append(['Total','','','','',total_well_counts['c']+total_well_counts['a'], total_well_counts['r']])
                    # 'Total assays': len(total_unique_assays)})
        #print(f"{plate_set_summary=}")
        return plate_set_summary, plate_set_headers

    
    def inputs_as_dataframe(self, caller_id=None):
        """ return the experiment contents as a pandas dataframe """
        rows, headers = self.summarise_inputs()
        if not rows:
            return None
        return pd.DataFrame(rows, columns=headers)


    def summarise_consumables(self, caller_id=None):
        """
        Return a dictionary of all consumables: primers, indexes, taq+water plates
        The descriptions we provide here will likely vary as we find new things to display
        Separated the taq/water plate IDs and volumes according to PCR stages
        TODO: separate custom and NGS primers
        """
        d = {'taqwater_pids_pcr1':[], 'taqwater_pids_pcr2':[], 'taq_vol_pcr1':0, 'taq_vol_pcr2':0,'water_vol_pcr1':0, 
                'water_vol_pcr2':0, 'primer_pids':[], 'primer_count_ngs':0, 'primer_count_custom':0, 'unique_primers':set(), 
                'primer_well_count':0, 'assay_primer_mappings':0, 'reference_files':[], 'unique_references':set(), 
                'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
        
        consumable_plate_purposes = set(['primer', 'taq_water', 'index'])
        for pid in self.plate_location_sample:
            if 'purpose' not in self.plate_location_sample[pid]:
                continue
            if self.plate_location_sample[pid]['purpose'] not in consumable_plate_purposes:
                continue
            plate = transaction.get_plate(self, pid)  # get the plate contents with all usage modifications applied

            if plate['purpose'] == 'taq_water':
                #Separate taq water PCR 1
                if 'pcr_stage' not in plate:
                    # legacy experiments only
                    continue
                elif plate['pcr_stage'] == 1:
                    d['taqwater_pids_pcr1'].append(pid)
                    for well in util.TAQ_WELLS:
                        if well in plate and 'volume' in plate[well]:
                            d['taq_vol_pcr1'] += plate[well]['volume']
                    for well in util.WATER_WELLS:
                        if well in plate and 'volume' in plate[well]:
                            d['water_vol_pcr1'] += plate[well]['volume']

                #Separate taq water PCR 2
                elif plate['pcr_stage'] == 2:
                    d['taqwater_pids_pcr2'].append(pid)
                    for well in util.TAQ_WELLS:
                        if well in plate and 'volume' in plate[well]:
                            d['taq_vol_pcr2'] += plate[well]['volume']
                    for well in util.WATER_WELLS:
                        if well in plate and 'volume' in plate[well]:
                            d['water_vol_pcr2'] += plate[well]['volume']

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
        d['assay_primer_mappings'] = len(self.assayfam_primers)
        return d


    def add_nimbus_outputs(self, nim_outputs, caller_id=None):
        """
        Copy Nimbus output files into the project folder
        """
        transactions={}
        #print(f'{type(nim_outputs)=} {len(nim_outputs)=}', file=sys.stderr)
        try:
            m('Starting add_nimbus_outputs', level='begin', caller_id=caller_id)
            for nim_output in nim_outputs:
                #print(f'{nim_output.name=}', file=sys.stderr) 
                fp = self.get_exp_fn(nim_output.name, trans=True)
                m(f"copying {fp} to experiment folder", level='info', caller_id=caller_id)
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
                final_fp = transaction.convert_pending_to_final(self,fp)
                if final_fp in self.uploaded_files:
                    m(f'{final_fp} already recorded as uploaded, overwriting', level='warning', caller_id=caller_id)
                self.uploaded_files[final_fp] = {'plates': list(plate_set), 'purpose':'DNA'} 
        except Exception as exc:
            m(f'could not upload Hamilton Nimbus output files, {exc}', level='error', caller_id=caller_id)
            return False
        transaction.add_pending_transactions(self, transactions)
        return True


    def get_dna_pids(self, dna_pids=None, echo_ready=False, plate_type=None, caller_id=None):
        """ 
        Return a list of available DNA plate ids 
        If echo_ready then include only those with echo_coc files uploaded
        """
        dpids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'dna']
        if dna_pids:
            dna_pids = [util.guard_pbc(d, silent=True) for d in dna_pids]
            dpids = [d for d in dpids if d in dna_pids]
        if plate_type and plate_type in util.PLATE_TYPES:
            dpids = [d for d in dpids if self.plate_location_sample[d]['plate_type']==util.PLATE_TYPES[plate_type]]
        if echo_ready:
            nfs, efs, xbcs = self.get_nimbus_filepaths()
            echo_ready_pids = []
            for nim in efs:
                echo_filename =  Path(nim).stem
                echo_ready_pids.append(util.guard_pbc(echo_filename.split('_')[-2], silent=True))
            dpids = [d for d in dpids if d in echo_ready_pids]
        return dpids
                    

    def get_dna_records(self, dna_pids=None, caller_id=None):
        """
        dna_fields=['samplePlate','sampleWell','sampleBarcode','strain','sex','alleleSymbol',
                 'alleleKey','assayKey','assays','assayFamilies','clientName','sampleName',
                 'dnaPlate','dnaWell','primer']
        """
        dna_pids = self.get_dna_pids(dna_pids=dna_pids)
        records = []
        for dna_bc in sorted(dna_pids):
            for well in self.plate_location_sample[dna_bc]['wells']:
                if 'barcode' not in self.plate_location_sample[dna_bc][well]:
                    continue
                if 'ngs_assays' not in self.plate_location_sample[dna_bc][well]:
                    continue
                for assay in self.plate_location_sample[dna_bc][well]['ngs_assays']:
                    if assay not in self.assay_assayfam:
                        continue
                    assayfam = self.assay_assayfam[assay]
                    if assayfam not in self.assayfam_primers:
                        continue
                    for primer in self.assayfam_primers[assayfam]:
                        record = {'samplePlate':self.plate_location_sample[dna_bc][well]['samplePlate'], 
                                'sampleWell':self.plate_location_sample[dna_bc][well]['sampleWell'], 
                                'sampleBarcode':self.plate_location_sample[dna_bc][well]['sampleBarcode'],
                                'assays':assay, 'assayFamilies':assayfam, 'primer':primer,
                                'dnaPlate':dna_bc, 'dnaWell':well}
                        record['strain'] = self.plate_location_sample[dna_bc][well].get('strain','')
                        record['sex'] = self.plate_location_sample[dna_bc][well].get('sex','')
                        if 'ngs_assay_records' in self.plate_location_sample[dna_bc][well]:
                            if assay in self.plate_location_sample[dna_bc][well]['ngs_assay_records']:
                                record['alleleSymbol'] = self.plate_location_sample[dna_bc][well]\
                                        ['ngs_assay_records'][assay].get('alleleSymbol','')
                                record['alleleKey'] = self.plate_location_sample[dna_bc][well]\
                                        ['ngs_assay_records'][assay].get('alleleKey','')
                                record['assayKey'] = self.plate_location_sample[dna_bc][well]\
                                        ['ngs_assay_records'][assay].get('assayKey','')
                        else:
                            record['alleleSymbol'] = ''
                            record['alleleKey'] = ''
                            record['assayKey'] = ''
                        record['clientName'] = self.plate_location_sample[dna_bc][well].get('clientName','')
                        record['sampleName'] = self.plate_location_sample[dna_bc][well].get('sampleName','')
                        records.append(record)
        return records    

        
    def add_pcr_plates(self, pcr_plate_list=[], caller_id=None):
        """
        Add one or more empty 384-well PCR plates to self.plate_location_sample
        """
        for pid in pcr_plate_list:
            p = util.guard_pbc(pid, silent=True)
            if p in self.plate_location_sample:
                m(f'barcode {p} already in use, skipping', level='failure', caller_id=caller_id)
                continue
            self.plate_location_sample[p] = {'purpose':'pcr', 
                                             'wells':set(), 
                                             'source':'user', 
                                            'plate_type':util.PLATE_TYPES['PCR384'], 
                                            'barcode':p}
        return True

    def get_pcr_pids(self, caller_id=None):
        """ return a list of available PCR plate ids """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'pcr']


    def get_primer_pids(self, pmr_pids=None, caller_id=None):
        """ return all primer pids, unless restricted by pmr_pids """
        if pmr_pids:
            primer_pids = [util.guard_pbc(ppid, silent=True) for ppid in pmr_pids \
                    if ppid in self.plate_location_sample and self.plate_location_sample[ppid]['purpose'] == 'primer']
        else:
            primer_pids = [util.guard_pbc(ppid, silent=True) for ppid in self.plate_location_sample \
                    if ppid in self.plate_location_sample and self.plate_location_sample[ppid]['purpose'] == 'primer']
        return primer_pids


    def get_available_primer_wells(self, pmr_pids=None, caller_id=None):
        """
        returns dictionary {pmr:[[pid,well,vol,doses],...]
        doses are the number of times a well can be aspirated from
        Used by get_available_primer_vols_doses() and for plating
        """
        primer_pids = self.get_primer_pids(pmr_pids)
        primer_wells_vols_doses = {}
        for ppid in primer_pids:
            for well in self.plate_location_sample[ppid]['wells']:
                if 'primer' in self.plate_location_sample[ppid][well]:
                    pmr = self.plate_location_sample[ppid][well]['primer']
                    if 'volume' in self.plate_location_sample[ppid][well]:
                        raw_vol = self.plate_location_sample[ppid][well]['volume']  # nanolitres
                        usable_vol = util.usable_volume(raw_vol, 'Echo384')
                        doses = util.num_doses(raw_vol, self.transfer_volumes['PRIMER_VOL'], 'Echo384')
                        if pmr not in primer_wells_vols_doses:
                            primer_wells_vols_doses[pmr] = []
                        primer_wells_vols_doses[pmr].append([ppid, well, usable_vol, doses])
        return primer_wells_vols_doses


    def get_available_primer_vols_doses(self, pmr_pids=None, caller_id=None):
        """
        return a dictionary of primer names and total volumes and uses across all wells
        """
        primer_wells_vols_doses = self.get_available_primer_wells(pmr_pids)
        primer_vols = defaultdict(int)
        primer_doses = defaultdict(int)
        for pmr in primer_wells_vols_doses:
            for record in primer_wells_vols_doses[pmr]:
                pid, well, vol, doses = record
                primer_vols[pmr] += vol  # volume is last
                primer_doses[pmr] += doses
        primer_vols_doses = {}
        for pmr in primer_vols:
            primer_vols_doses[pmr] = tuple([primer_vols[pmr], primer_doses[pmr]])
        return primer_vols_doses


    def get_primers_for_assays(self, assays, caller_id=None):
        """
        For each assay in assays, find the matching assay family and the corresponding primers
        Return a dictionary[assay] = [primers] and missing_assays = set()
        """
        assay_primers = {}
        missing_assays = set()
        for assay in assays:           
            if assay in assay_primers:
                continue  # we've already seen this assay
            if assay not in self.assay_assayfam:
                missing_assays.add(assay)
                continue
            assayfam = self.assay_assayfam[assay]
            if assayfam not in self.assayfam_primers:
                missing_assays.add(assay)
                continue
            primers = self.assayfam_primers[assayfam]
            assay_primers[assay] = primers
        return assay_primers, missing_assays


    def get_assay_primer_usage(self, dna_pids, caller_id=None):
        """
        For a given set of DNA plates, return the number of times each assay and primer is needed
        """
        assay_usage = defaultdict(int)
        primer_usage = defaultdict(int)
        missing_assays = set()
        if not dna_pids:
            return assay_usage, primer_usage
        
        dpids = [pid for pid in dna_pids if util.guard_pbc(pid, silent=True) in self.plate_location_sample \
                and self.plate_location_sample[pid]['purpose'] == 'dna']
        
        if not dpids:
            m(f'no DNA plates found matching DNA plate IDs: {dna_pids}', level='error', caller_id=caller_id)
            return assay_usage, primer_usage

        for dpid in dpids:  
            wells = self.plate_location_sample[dpid].get('wells', None)
            if not wells:  # get from source plates
                sources = self.plate_location_sample[dpid].get('source', None)
                #print(f"{type(sources)=} {sources=} ")
                if not sources:
                    continue
                spids = [src for src in sources.strip().split(',') if util.is_guarded_pbc(src)]
                for spid in spids:
                    wells = self.plate_location_sample[spid].get('wells', None)
                    for well in self.plate_location_sample[dpid]['wells']:       
                        assays = self.plate_location_sample[dpid][well].get('ngs_assays', [])
                        for assay in assays:  
                            assay_usage[assay] += 1
                        assays_primers, missing = self.get_primers_for_assays(assays)
                        for missing_assay in missing:
                            missing_assays.add(missing_assay)
                        for a in assays_primers:
                            for primer in assays_primers[a]:
                                primer_usage[primer] += 1
            else:        
                for well in self.plate_location_sample[dpid]['wells']:       
                    assays = self.plate_location_sample[dpid][well].get('ngs_assays', [])
                    for assay in assays:  
                        assay_usage[assay] += 1
                    assays_primers, missing = self.get_primers_for_assays(assays)
                    for missing_assay in missing:
                            missing_assays.add(missing_assay)
                    for a in assays_primers:
                        for primer in assays_primers[a]:
                            primer_usage[primer] += 1
        # only report this on the Primer page
        if missing_assays and caller_id not in [None, 'display_feedback']:
            m(f'The following assays are expected but not found in the assay-primer map (assay list file): {missing_assays}', level='warning', no_log=True, caller_id=caller_id)
        return assay_usage, primer_usage


    def get_index_avail(self, index_pids=None, caller_id=None):
        """
        Provided a list of index PIDs and return dictionaries of fwd and rev indexes
        Returns:
            [fwd_index]={'well_count':int, 'avail_transfers':int, 'avail_vol':nl}
            [rev_index]={'well_count':int, 'avail_transfers':int, 'avail_vol':nl}
        """
        if index_pids:
            idx_pids = [util.guard_pbc(ip, silent=True) for ip in index_pids]
            index_pids = [pid for pid in idx_pids if self.plate_location_sample[pid]['purpose'] == 'index']
        else:
            index_pids = [pid for pid in self.plate_location_sample if self.plate_location_sample[pid]['purpose'] == 'index']
                
        fwd_idx = {}
        rev_idx = {}

        for idx_pid in index_pids:
            idx_plate = transaction.get_plate(self, idx_pid)
          
            for well in idx_plate['wells']:
                if 'idt_name' not in idx_plate[well]:
                    continue
                name = idx_plate[well]['idt_name']
                if 'i7F' in name:
                    if name not in fwd_idx:
                        fwd_idx[name] = {'well_count':0, 'avail_transfers':0, 'avail_vol':0}
                    fwd_idx[name]['well_count'] += 1
                    if 'volume' in idx_plate[well]:
                        vol = util.usable_volume(idx_plate[well]['volume'], 'Echo384')
                        doses = util.num_doses(idx_plate[well]['volume'], self.transfer_volumes['INDEX_VOL'], 'Echo384')
                        fwd_idx[name]['avail_vol'] += vol
                        fwd_idx[name]['avail_transfers'] += doses
                elif 'i5R' in name:
                    if name not in rev_idx:
                        rev_idx[name] = {'well_count':0, 'avail_transfers':0,'avail_vol':0,}
                    rev_idx[name]['well_count'] += 1
                    if 'volume' in idx_plate[well]:
                        vol = util.usable_volume(idx_plate[well]['volume'], 'Echo384')
                        doses = util.num_doses(idx_plate[well]['volume'], self.transfer_volumes['INDEX_VOL'], 'Echo384')
                        rev_idx[name]['avail_vol'] += vol
                        rev_idx[name]['avail_transfers'] += doses
                else:
                    m(f'unexpected index name: {name}', level='warning', caller_id=caller_id)

        return fwd_idx, rev_idx


    def get_index_pairs_avail(self, index_pids, caller_id=None):
        """ returns all available pairings of barcode ends """
        if not index_pids:
            return []
        fwd_idx, rev_idx = self.get_index_avail(index_pids)
        fwd_idx_names = list(fwd_idx.keys())
        rev_idx_names = list(rev_idx.keys())
        pairs_allocated = []
        idx_uses = {}
        total_fwd = len(fwd_idx_names)
        total_rev = len(rev_idx_names)
        # cycle through all combinations in an evenly
        for i in range(total_fwd*total_rev):
            fi = fwd_idx_names[i%total_fwd]
            ri = rev_idx_names[i%total_rev]
            if fi not in idx_uses:
                idx_uses[fi] = 0
            if ri not in idx_uses:
                idx_uses[ri] = 0
            if idx_uses[fi] == fwd_idx[fi]['avail_transfers']:
                continue
            if idx_uses[ri] == rev_idx[ri]['avail_transfers']:
                continue
            pairs_allocated.append((fi,ri))
            idx_uses[fi] += 1
            idx_uses[ri] += 1
        return pairs_allocated
    

    def get_taqwater_avail(self, taqwater_bcs=None, transactions=None, pcr_stage=None, caller_id=None):
        """ 
        Returns (int) taq and (int) water volumes (in nanolitres) loaded as available 
        Will work from a list of barcodes if provided, or will calculate for all available plates
        If pcr_stage is provided, returns the volumes specific to that stage.
        transactions (Optional): a dictionary of plates and their unprocessed changes for accurate calculation
        pcr_stage: 1 (primers) or 2 (index). If None, returns volumes from all taq/water plates
        """
        taq_avail = 0
        water_avail = 0
        pids = []
        if taqwater_bcs is not None:
            bcs = taqwater_bcs
        else:
            bcs = self.plate_location_sample
        
        #Separate based on PCR stage, otherwise get plate IDs from all plates.
        if pcr_stage:
            pids = [p for p in bcs if self.plate_location_sample[p]['purpose'] == 'taq_water'\
                     and self.plate_location_sample[p].get('pcr_stage') == pcr_stage]
        else:
            pids = [p for p in bcs if self.plate_location_sample[p]['purpose'] == 'taq_water']


        for pid in pids:
            tw_plate = transaction.get_plate(self, pid) # get the plate records with adjusted usage
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


    def get_taqwater_volumes_primer(self, num_reactions, caller_id=None):
        """
        Returns the taq and water volumes required for the primer stage in nl
        Args:
            num_reactions(int):
        """
        taq_vol = num_reactions * self.transfer_volumes['PRIMER_TAQ_VOL']
        water_vol = num_reactions * self.transfer_volumes['PRIMER_WATER_VOL']
        return taq_vol, water_vol
        

    def get_taqwater_req_vols_index(self, num_reactions, caller_id=None):
        """
        Returns the taq and water volumes required for the index stages in nl
        Args:
            num_reactions(int): number of reactions for index stage, including amplicon
        """
        water_vol = num_reactions * self.transfer_volumes['INDEX_WATER_VOL']
        taq_vol = num_reactions * self.transfer_volumes['INDEX_TAQ_VOL'] 
        return taq_vol, water_vol
    

    def get_taqwater_pids(self, pcr_stage=None, caller_id=None):
        """
        Return plate IDs for taq/water plates.
        If PCR stage is provided, return the taq/water plate for that stage. 
        pcr_stage: 1 (primer) or 2 (index). If None, return all plate IDs for taq/water plates
        """
        pids = []
        #Plate ID specific to PCR stage (numeric 1 or 2). Otherwise get all taq/water plates. 
        if pcr_stage:
            for p in self.plate_location_sample:
                if self.plate_location_sample[p]['purpose'] == 'taq_water':
                    if 'pcr_stage' not in self.plate_location_sample[p]:
                        pids.append(p)  # legacy experiments only
                    elif self.plate_location_sample[p]['pcr_stage'] == pcr_stage:
                        pids.append(p)       
        else:
            pids = [p for p in self.plate_location_sample
                    if self.plate_location_sample[p]['purpose'] == 'taq_water']
        return pids

    
    def generate_nimbus_inputs(self, caller_id=None):
        success = generate.nimbus_gen(self)
        return success

    
    def get_nimbus_filepaths(self, caller_id=None):
        """ Return the lists of nimbus input files, echo input file (nimbus outputs), 
            and barcodes that are only seen in nimbus """
        #print("In get_nimbus_filepaths")
        nimbus_input_filepaths, echo_input_paths, xbc = generate.match_nimbus_to_echo_files(self)
        return nimbus_input_filepaths, echo_input_paths, xbc

    
    def get_amplicon_pids(self, caller_id=None):
        """ return a list of user supplied amplicon plate ids """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'amplicon']


    def get_num_amplicon_wells(self, amplicon_pids=None, caller_id=None)->int:
        """ 
        Return the number of wells in the chosen amplicon plates. If amplicon_pids is None, get number of wells for
        all amplicons in experiment
        """
        if amplicon_pids:
            num_wells = sum(len(self.plate_location_sample[pid]['wells']) \
                    for pid in self.plate_location_sample if pid in amplicon_pids)
        else:
            num_wells = sum(len(plate['wells']) for plate in self.plate_location_sample.values() if plate['purpose'] == 'amplicon')
        
        return num_wells
    

    def get_num_reactions(self, pcr_pids=None, amplicon_pids=None, caller_id=None):
        num_reactions = 0
        if amplicon_pids:
            num_reactions += self.get_num_amplicon_wells(amplicon_pids)
        if pcr_pids:
            num_reactions += self.get_num_pcr_wells(pcr_pids)
        return num_reactions


    # def generate_echo_primer_survey(self, primer_survey_filename='primer-svy.csv', caller_id=None):
    #     """ 
    #     Generate a primer survey file for use by the Echo. Replaces echovolume.py
    #     Not strictly necessary, but the old code reads this file in making the picklists.
    #     DEPRECATED
    #     """
    #     if self.locked:
    #         m('cannot generate primer survey file while lock is active.', level='failure', caller_id=caller_id)
    #         return False
    #     primer_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'primer']
    #     header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'plate position on Echo 384 PP',
    #             'primer names pooled', 'volume']
    #     transactions = {}                                       
    #     primer_survey_fn = self.get_exp_fn(primer_survey_filename, trans=True)
    #     transactions[primer_survey_fn] = {} # add plates and modifications to this
    #     try:
    #         with open(primer_survey_fn, 'wt') as fout:
    #             print(','.join(header), file=fout)
    #             for i,pid in enumerate(primer_pids):
    #                 transactions[primer_survey_fn][pid] = {}
    #                 plate=self.plate_location_sample[pid]
    #                 for well in util.row_ordered_384:
    #                     if well not in plate['wells']:
    #                         continue
    #                     if 'primer' not in plate[well]:
    #                         continue
    #                     if 'volume' not in plate[well]:
    #                         continue
    #                     outline = ','.join([f"Source[{i+1}]",util.unguard_pbc(pid,silent=True),plate['plate_type'],
    #                             well,plate[well]['primer'],f"{int(plate[well]['volume'])/1000}"])
    #                     print(outline, file=fout)
    #     except Exception as exc:
    #         m(f'could not write primer survey {exc}', level='error', caller_id=caller_id)
    #         return False
    #     transaction.add_pending_transactions(self,transactions)
    #     m(f"written Echo primer survey to {primer_survey_fn}", level='success', caller_id=caller_id)
    #     return True


    def check_plate_presence(self, pids, purpose, caller_id=None):
        """ 
        Check given plate for presence in self.plate_location_sample and compare purpose
        return success (and messages by reference)
        Required by exp.check_ready_pcr1() and exp.check_ready_pcr2()
        """
        if not pids:
            m(f'no plate IDs given for {purpose}', level='warning', 
                    caller_id=caller_id)
            return False
        
        success = True
        for pid in pids:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != purpose:
                    m(f"plate already exists with PID {pid} with purpose "+\
                            f"{self.plate_location_sample[pid]['purpose']}, expected {purpose}", 
                            level='error', caller_id=caller_id)
                    success = False
            else:
                m(f"No plate exists with barcode {pid}", level='critical', caller_id=caller_id)
                success = False
        return success
    

    def check_ready_pcr1(self, selected_pids, caller_id=None):
        """
        Check that everything required to successfully generate PCR1 picklists is available
        TODO: Check that volumes and wells are sufficient!
        Returns success
        """
        success = True
        dna_pids = selected_pids['dna']
        pcr_pids = selected_pids['pcr']
        taqwater1_pids = selected_pids['taqwater1']
        primer_pids = selected_pids['primer']
        if not dna_pids or not pcr_pids or not taqwater1_pids or not primer_pids:
            m('DNA plates, PCR plates, and taq+water plates are all required for primer PCR', level='failure', caller_id=caller_id)

        dna_success = self.check_plate_presence(dna_pids, 'dna', caller_id)
        if not dna_success:
            success = False

        pcr_success = self.check_plate_presence(pcr_pids, 'pcr', caller_id)
        if not pcr_success:
            success = False

        taq_success = self.check_plate_presence(taqwater1_pids, 'taq_water', caller_id)
        if not taq_success:
            success = False
      
        return success


    def check_ready_pcr2(self, selected_pids, amplicon_only=False, caller_id=None):
        """
        Check the everything required to successfully generate PCR2 picklists is available
        TODO: Check that volumes and wells are sufficient!
        Amplicon plates are optional
        Return success
        """
        success = True

        pcr_pids = selected_pids['pcr']
        taqwater2_pids = selected_pids['taqwater2']
        index_pids = selected_pids['index']
        amplicon_pids = selected_pids['amplicon']
        
        if amplicon_only:
            if not amplicon_pids or not taqwater2_pids or not index_pids:
                m('Amplicon plates, taq+water plates, and index plates are all required for indexing',
                        level='error', caller_id=caller_id)
                return False
        else:
            if not pcr_pids or not taqwater2_pids or not index_pids:
                m('PCR plates, taq+water plates, and index plates all required for indexing', 
                        level='error', caller_id=caller_id)
                return False

        if not amplicon_only:
            pcr_success = self.check_plate_presence(pcr_pids, 'pcr', caller_id)
            if not pcr_success:
                success = False
        else:
            amplicon_success = self.check_plate_presence(amplicon_pids, 'amplicon', caller_id)
            if not amplicon_success:
                success = False
                
        taq_success = self.check_plate_presence(taqwater2_pids, 'taq_water', caller_id)
        if not taq_success:
            success = False

        index_success = self.check_plate_presence(index_pids, 'index', caller_id)
        if not index_success:
            success = False

        return success
    

    def generate_echo_PCR1_picklists(self, selected_pids, force=False, caller_id=None):
        """
        Calls generate.generate_echo_PCR1_picklist() to do the work, needs a set of accepted DNA_plates,
        the final assay list, and a set of destination PCR plate barcodes, taq+water plates, primer plates, primer volumes.

        args:
            selected_pids: dictionary of selected plate IDs
            force: if True, will not stop if primers allocations run out
            caller_id: the calling function

        Returns True on success
        """
        dna_pids = selected_pids['dna']
        pcr_pids = selected_pids['pcr']
        taqwater1_pids = selected_pids['taqwater1']
        success = generate.generate_echo_PCR1_picklist(self, dna_pids, pcr_pids, 
                taqwater1_pids,force=force,caller_id=caller_id)
        if not success:
            m('could not generate PCR1 picklists correctly', level='error', caller_id=caller_id)
            return False
        m('generated PCR1 picklists', level='success', caller_id=caller_id)
        return True


    def get_echo_PCR1_picklist_filepaths(self, caller_id=None):
        """
        Return file paths for PCR1_dna-picklist_XXX.csv, 
        PCR1_primer-picklist_XXX.csv, PCR1_taqwater-picklist_XXX.csv
        """
        all_files = os.listdir(self.get_exp_dn())
        dna_picklist_paths = []
        primer_picklist_paths = []
        taqwater_picklist_paths = []
        for f in all_files:
            if f.startswith('PCR1_dna-picklist_') and f.endswith('.csv'):
                dna_picklist_paths.append(self.get_exp_fn(f))
            elif f.startswith('PCR1_primer-picklist_') and f.endswith('.csv'):
                primer_picklist_paths.append(self.get_exp_fn(f))
            elif f.startswith('PCR1_taqwater-picklist_') and f.endswith('.csv'):
                taqwater_picklist_paths.append(self.get_exp_fn(f))
        return dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths


    def generate_echo_PCR2_picklists(self, selected_pids, caller_id=None):
        """
        Calls echo_index.generate_echo_PCR2_picklist() to do the work, 
        needs a set of destination PCR plate barcodes, 
        taq+water plates, index plates, index volumes, 
        and optionally any amplicon plates (post-PCR).
        Returns True on success
        """
        pcr_pids = selected_pids['pcr']
        index_pids = selected_pids['index']
        taqwater2_pids = selected_pids['taqwater2']
        amplicon_pids = selected_pids['amplicon']
            
        success = self.generate_echo_index_survey(index_pids, caller_id=caller_id)
        if not success:
            m('failed to generate index survey file', level='error', caller_id=caller_id)
            return False
        #m('generated index survey file', level='success', caller_id=caller_id)
        success = generate.generate_echo_PCR2_picklist(self, pcr_pids, index_pids, 
                taqwater2_pids, amplicon_pids, caller_id=caller_id)
        if not success:
            m('could not generate PCR2 (index) picklists correctly', level='error', caller_id=caller_id)
            return False
        #m('generated PCR2 (index) picklists', level='success', caller_id=caller_id)
        return success


    def get_echo_PCR2_picklist_filepaths(self, caller_id=None):
        """
        Return file paths for PCR2_index-picklist_XXX.csv, PCR2_taqwater-picklist_XXX.csv
        These picklists should include amplicon plate destinations if provided
        """
        all_files = os.listdir(self.get_exp_dn())
        index_picklist_paths = [] 
        taqwater_picklist_paths = []
        for f in all_files:
            if f.startswith('PCR2_index-picklist_') and f.endswith('.csv'):
                index_picklist_paths.append(self.get_exp_fn(f))
            elif f.startswith('PCR2_taqwater-picklist_') and f.endswith('.csv'):
                taqwater_picklist_paths.append(self.get_exp_fn(f))
        return index_picklist_paths, taqwater_picklist_paths


    def delete_plate(self, pid, caller_id=None):
        """
        Soft-delete plate with the given pid. Note: should we check for matching files and delte those?
        """
        success = True
        if pid in self.plate_location_sample:
            if pid not in self.deleted_plates:
                self.deleted_plates[pid] = []
            self.deleted_plates[pid].append(deepcopy(self.plate_location_sample[pid]))
            # need to find this in reproducible steps and delete
            del self.plate_location_sample[pid]
            m(f'moved plate {pid} to deleted plate bin', level='warning', dest=('noGUI',))
            if pid in self.dest_sample_plates:
                del self.dest_sample_plates[pid]
                m(f'removed 384-well DNA plate entry {pid}', level='warning', dest=('noGUI',))
        elif pid in self.dest_sample_plates:
            del self.dest_sample_plates[pid]
            m(f'removed 384-well DNA plate entry {pid}', level='success', dest=('noGUI',))
        else:
            m(f'{pid} has no definition loaded', level='error', dest=('noGUI',))
            success = False
        return success
    

    def get_stages(self, caller_id=None):
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
                if type(self.pending_steps[file_name]) == 'dict':
                    pids = [util.unguard_pbc(pid, silent=True) for pid in self.pending_steps[file_name].keys()]
                if type(self.pending_steps[file_name]) == 'list':
                    pids = []
                    for thing in self.pending_steps[file_name]:
                        for pid in thing:
                            pids.append(util.unguard_pbc(pid, silent=True))
                else:
                    #pids = []
                    #for pid in self.pending_steps[file_name]:
                    #    if type(pid) is str:
                    #        pids.append(pid)
                        
                    pids = [util.unguard_pbc(pid, silent=True) for pid in \
                            self.pending_steps[file_name] if type(pid) is str]
                stages.append([str(counter), file_name, ', '.join(pids), 'pending'])
        return stages, header


    def get_stage2_pcr_plates(self, caller_id=None):
        """
        Used by Indexing stage to find all the used PCR plate IDs
        """
        stage2_pcr_plates = set()
        stage2_fn = self.get_exp_fn('Stage2.csv')
        if not os.path.exists(stage2_fn):
            return stage2_pcr_plates

        with open(stage2_fn, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue  # header
                cols = line.split(',')
                stage2_pcr_plates.add(cols[8].strip())
        return stage2_pcr_plates
    

    def add_pcr_wells(self, exp, pcr_pids, dna_pids, caller_id=None):
        """
        get the PCR wells based on the Stage2.csv
        """
        dna_records = self.get_dna_records(dna_pids)
        success = generate.allocate_primers_to_dna(exp, dna_records=dna_records)
        dna_records = sorted([d for d in dna_records if d['dnaPlate'] and d['primerPlate']], 
                key=lambda x: (x['samplePlate'],util.padwell(x['sampleWell'])))
        
        # allocate PCR plate wells
        wells = [r+str(c+1) for c in range(24) for r in"ABCDEFGHIJKLMNOP"]
        pcrWells = [(p, w) for p in sorted(pcr_pids) for w in wells]

        for i, (dr, pcrw) in enumerate(zip(dna_records, pcrWells)):
            self.plate_location_sample[pcrw[0]]['wells'].add(pcrw[1])
            if pcrw[1] not in self.plate_location_sample[pcrw[0]]:
                self.plate_location_sample[pcrw[0]][pcrw[1]] = dr   


    def get_num_pcr_wells(self, pcr_pids=None, caller_id=None)->int:
        """ 
        Return the number of wells in the PCR plates. If pcr_pids is None, get number of wells for
        all PCR plates in experiment
        """
        if pcr_pids:
            num_wells = sum(len(self.plate_location_sample[pid]['wells']) \
                            for pid in self.plate_location_sample if pid in pcr_pids)
        else:
            num_wells = sum(len(plate['wells']) for plate in self.plate_location_sample.values() \
                             if plate['purpose'] == 'pcr')
        return num_wells


    def get_file_usage(self, caller_id=None):
        """
        Gather file records for display
        """
        file_usage = {}
        all_files = [fn for fn in self.uploaded_files.keys() if not fn.startswith('_')]
        del_files = []
        for fn in all_files:
            if not Path.exists(Path(fn)):
                del_files.append(fn)
        for fn in del_files:
            m(f'{fn} does not exist! Removing file record', level='warning', caller_id=caller_id)
            self.del_file_record(fn, caller_id=caller_id)
                
        for filename in [fn for fn in self.uploaded_files.keys() if not fn.startswith('_')]:
            file_usage[filename] = {'date modified':None, 'purpose':None}
            if 'purpose' in self.uploaded_files[filename]:
                file_usage[filename]['purpose'] = self.uploaded_files[filename]['purpose']
            if 'date modified' in self.uploaded_files[filename]:
                file_usage[filename]['date modified'] = self.uploaded_files[filename]['date modified']
            else:
                ft = os.path.getmtime(filename)
                ft = datetime.datetime.fromtimestamp(ft).strftime('%Y/%m/%d %H:%M:%S')    
                file_usage[filename]['date modified'] = ft
                self.uploaded_files[filename]['date modified'] = ft
                
            #if 'plates' in self.uploaded_files[filename]:
            #    for pid in self.uploaded_files[filename]['plates']:
            #        if util.is_guarded(pid):
            #            pid = util.unguard_pbc(pid)
            #        file_usage[filename]['plates'].append(pid)
            #file_usage[filename]['plates'] = ', '.join(file_usage[filename]['plates'])
        return file_usage
    

    def get_plate_usage(self, caller_id=None):
        """return list of plate information to display"""
        plate_usage = []
        for p in self.plate_location_sample:
            if 'purpose' not in self.plate_location_sample[p]:
                continue            
            pid = util.unguard_pbc(p, silent=True)
            plate_usage.append([pid, len(self.plate_location_sample[p]['wells']), self.plate_location_sample[p]['purpose']])
        return plate_usage


    def get_index_pids(self, caller_id=None):
        """
        Used by Indexing stage
        """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'index']


    def generate_echo_index_survey(self, user_index_pids, index_survey_filename='index-svy.csv', caller_id=None):
        """ 
        Generate an index survey file for use by the Echo. Replaces echovolume.py 
        Not strictly necessary, but the old code reads this file in making the picklists.
        Called internally by generate_echo_PCR2_picklist_interface()
        Requires the set of user chosen user_index_pids
        """
        if self.locked:
            m('cannot generate index survey while lock is active.', level='failure', caller_id=caller_id)
            return False
        index_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'index' and p in user_index_pids]
        header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'Source well', 
                'Name for IDT ordering', 'index', 'Name', 'Oligo * -phosphothioate modif. against exonuclease', 'volume']
        fn = self.get_exp_fn(index_survey_filename)
        if os.path.exists(fn):
            m(f'overwriting index survey file {fn}', level='warning', caller_id=caller_id)
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
        m(f'index survey file written to {fn}', level='success', caller_id=caller_id)
        return True


    def add_standard_taqwater_plates(self, plate_barcodes, pcr_stage, caller_id=None):  # <- they're the same for primer and barcode stages
        """ These plates have a fixed layout with water in row A and taq in row B of a 6-well reservoir plate.
        See echo_primer.py or echo_barcode.py mytaq2()
        We may need to also store volumes at some stage, but for now that isn't necessary
        pcr_stage(int): which PCR experiment the plate comes from (1 or 2)
        """
        if self.locked:
            m('cannot add taq+water plates while lock is active.', level='failure', caller_id=caller_id)
            return False
        try:
            for plate_barcode in plate_barcodes:
                pid = util.guard_pbc(plate_barcode, silent=True)
                if pid in self.plate_location_sample:
                    if self.plate_location_sample[pid]['purpose'] != 'taq_water':
                        m(f"cannot add taq+water plate {util.unguard(pid)}. "+\
                                f"A plate with this barcode already exists with purpose "+\
                                f"{self.plate_location_sample[pid]['purpose']}", 
                                level='error', caller_id=caller_id)
                        return False
                    elif self.plate_location_sample['pcr_stage'] == pcr_stage:
                        m(f"cannot add taq+water plate {util.unguard(pid)}."+\
                                f"A taq/waterplate with this barcode already exists with for PCR {pcr_stage}",
                                level='error', caller_id=caller_id)
                        return False
                    else:
                        m(f"overwriting existing taq/water plate entry with barcode {util.unguard(pid)}",
                                level='warning', caller_id=caller_id)
                else:
                    self.plate_location_sample[pid] = {}
                self.plate_location_sample[pid]['purpose'] = 'taq_water'
                self.plate_location_sample[pid]['plate_type'] = util.PLATE_TYPES['Echo6']
                self.plate_location_sample[pid]['pcr_stage'] = pcr_stage
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
            m(f"adding taq+water plate barcode failed {plate_barcodes=} {exc}",
                    level='error', caller_id=caller_id)
            return False
        return True


    def get_miseq_samplesheets(self, caller_id=None):
        """ return the MiSeq-XXX.csv samplesheet, if it exists """
        miseq_fps = list(Path(self.get_exp_dn()).glob('MiSeq_*.csv'))
        return miseq_fps  


    def check_sequence_upload_ready(self, caller_id=None):
        """
        Prevent users from uploading sequence files without references being loaded, or without generating Stage3 or MiSeq files
        Requires a list of messages for the GUI which can be appended to (pass by reference)
        """
        success = True
        # check for MiSeq.csv and Stage3 files
        fn1 = self.get_exp_fn(f'MiSeq_{self.name}.csv')
        fn2 = self.get_exp_fn(f'Stage3.csv')
        fns = [fn1, fn2]
        for fn in fns:
            if not Path(fn).exists():
                success = False
                m(f"{fn} not present, cannot proceeed with sequence upload", level='error', caller_id=caller_id)

        # check that at least one reference sequence has been uploaded
        if len(self.reference_sequences) == 0:
            m('No reference sequences have been uploaded yet, please add these before running allele calling', level='warning', 
                    caller_id=caller_id)

        return success

 
    def check_allele_calling_ready(self, caller_id=None):
        """ 
        Return True if everything needed for allele calling is present
            - Stage3.csv is present - done in self.check_sequence_upload_ready()
            - Miseq file is present - done in self.check_sequence_upload_ready()
            - check that reference sequences are uploaded - done in self.check_sequence_upload_ready()
            - the raw directory is present (can be empty)

        Requires a caller_id to be given for returning feedback to the GUI
        """
        success = self.check_sequence_upload_ready(caller_id)

        # check whether the raw directory for FASTQs exists yet - implies at least one FASTQ has been uploaded        
        dn = self.get_exp_fn(f'raw')
        if not Path(dn).exists() or not Path(dn).is_dir():
            success = False
            m(f"{dn} does not exist", level='error', caller_id=caller_id)

        return success


    def add_custom_volumes(self, custom_volumes:dict, caller_id=None) -> bool:
        """
        Modify transfer volume info
        """
        try:
            for value in list(custom_volumes.values()):
                value = int(value)
            self.transfer_volumes = custom_volumes.copy()
            m(f'Custom volumes set', level='success', caller_id=caller_id)
            return True
        except (TypeError, ValueError) as exc:
            m(f'during volume setting {exc}', level='error', caller_id=caller_id)
            return False
    
     
    def save(self):
        """ save experiment details to self.name/experiment.json, returns True on success and False on fail. Assumes the correct working directory """
        try:
            self._finalizer.detach()
            del self._finalizer
            exp = jsonpickle.encode(self, indent=4, keys=True, warn=True, handle_readonly=True)
            if not self.name:
                return False
            with open(os.path.join('run_'+self.name, EXP_FN), 'wt') as f:
                print(exp, file=f)
            print(f'Successfully saved experiment {self.name}', flush=True)
        except Exception as exc:
            print(f"Error saving {self.name=} {exc}", flush=True, file=sys.stderr)
            return False
        return True


    def log(self, message, level, func=None, func_line=None, call_func=None, call_line=None):
        """ 
        Always add the date/time, function where the log was run and the caller function to the log
        level is in {debug, info, warning, error, critical, begin, end, success, failure}
        """
        now = datetime.datetime.now()
        t = f"{now:%Y-%m-%d %H:%M}"
        if func:
            func=func
        else:
            func = sys._getframe(1).f_code.co_name
            
        if func_line:
            func_line = func_line
        else:
            func_line = inspect.getframeinfo(sys._getframe(1)).lineno
            
        if call_func:
            caller = call_func
        else:
            caller = sys._getframe(2).f_code.co_name
            
        if call_line:
            caller_line = call_line
        else:
            caller_line = inspect.getframeinfo(sys._getframe(2)).lineno
        
        self.log_entries.append([t, func, func_line, caller, caller_line, level.capitalize(), message])
        #if (now - self.log_time).seconds > 10 or level in ['Error','Critical','Failure','End']:
        #    self.save()


    def get_log_header(self):
        return ['Time', 'Function name', 'Func line', 'Calling function', 'Call line', 'Level', 'Message']


    def get_log(self, num_entries=100):
        """ return a chunk of the log. -1 gives everything """
        if num_entries == -1:
            return self.log_entries[::-1]
        return self.log_entries[:-(num_entries+1):-1]


    def clear_log(self, message_type=None):
        """ DEPRECATED: needs replacing with better functionality
        delete all the log entries with level='Debug' 
        """
        if message_type is not None:
            if message_type not in ['Debug','Info','Warning','Error','Critical','Begin','End','Success','Failure']:
                self.log_entries.append(f"Error: didn't recognise message type {message_type}")
                return
        clean_log = []
        if message_type is not None:
            for entry in self.log_entries:
                if entry[-2] != message_type:
                    clean_log.append(entry)
        self.log_entries = clean_log


def load_experiment(exp_path):
    """ load experiment details from experiment.json, or return None. exp_path is the folder path only """
    
    exp_file_path = os.path.join(exp_path, EXP_FN)
    print('Attempting to load from path:', exp_file_path)
    if os.path.exists(exp_file_path):
        with open(exp_file_path, 'rt') as f:
            expt = jsonpickle.decode(f.read(), keys=True, handle_readonly=True) # this seems totally flaky, but it must come down to what the contents are?
        if isinstance(expt, Experiment):
            return expt
        elif isinstance(expt, str):
            expt = dict(expt)
        exp = Experiment(expt['name'])
        exp.__setstate__(expt)
        return exp
    else:
        return None