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
from shutil import copyfileobj
import json
from Bio import SeqIO
from copy import deepcopy
from math import ceil, floor
from pathlib import Path
import inspect  
from collections import defaultdict

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
        self.denied_assays = []
        self.denied_primers = []
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
        self.uploaded_files = {}
        self.reproducible_steps = []  # a strict subset of the log that only includes working instructions
        self.pending_steps = set() # reproducible steps (files) that are awaiting user approval to replace existing steps
        self.pending_uploads = set() # not a reproducible step, but needs to be cleared by user
        self.log_entries = []  # use self.log(message, level='') to add to this
        self.log_time = datetime.datetime.now()  # save the log after a certain time has elapsed, to avoid too much IO
        self.transfer_volumes = {'DNA_VOL':util.DNA_VOL, 'PRIMER_VOL':util.PRIMER_VOL, 'PRIMER_TAQ_VOL':util.PRIMER_TAQ_VOL,
                'PRIMER_WATER_VOL':util.PRIMER_WATER_VOL, 'INDEX_VOL':util.INDEX_VOL, 'INDEX_TAQ_VOL':util.INDEX_TAQ_VOL,
                'INDEX_WATER_VOL':util.INDEX_WATER_VOL}  # load the defaults from util.py
        self.dead_volumes = util.DEAD_VOLS
        self.cap_volumes = util.CAP_VOLS
        

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
        #self.save()

    def unlock(self):
        self.log('Info: Locking/unlocking of experiment is currently disabled')
        #self.log('Info: Unlocking experiment. Modification is now possible. There should be good reason for this!')
        #self.locked = False
        #self.save()


    ### functions for returning locally held file paths

    def get_exp_dn(self, subdir=None):
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
                    self.log(f'Critical: failed to make directory {str(dp2)}, which already exists as a file')
            dirname = str(dp2)
            return str(dirname)
        else:                                                       
            self.log(f'Critical: {subdir=} not in {allowed_subdirs=}')

    def get_exp_fn(self, filename, subdir=None, trans=False):
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
      

    def get_raw_fastq_pairs(self):
        """ return a sorted list of tuple(R1_path, R2_path) to raw FASTQ files """
        valid_pairs = []
        rdp = self.get_raw_dirpath()
        r1s = [rdp/Path(f) for f in os.listdir(rdp) if f.endswith('.fastq.gz') and '_R1_' in f]
        for r1 in r1s:
            r2 = Path(str(r1).replace('_R1_001','_R2_001'))
            if r1.is_file() and r2.is_file():
                valid_pairs.append((r1,r2))
            else:
                if not r1.is_file():
                    self.log(f'Warning: {r1} expected raw FASTQ file does not exist')
                elif not r2.is_file():
                    self.log(f'Warning: {r2} expected raw FASTQ file does not exist')
        return sorted(valid_pairs)
    
    ### functions for managing data file records

    def add_file_record(self, file_name, PIDs=None, purpose=None):
        """
        Add an uploaded file to the experiment with associated plates and purpose
        """
        if file_name in self.uploaded_files:
            self.log(f'Warning: file name {file_name} has already been uploaded, overwriting')

        if not PIDs:
            PIDs = []

        file_md5 = util.get_md5(file_name)
        self.uploaded_files[file_name] = {'plates': [util.guard_pbc(PID, silent=True) for PID in PIDs], 
                'purpose': purpose, 'md5':file_md5}
        return True


    def mod_file_record(self, existing_name, new_name=None, new_purpose=None, extra_PIDs=None, PID_name_updates=None):
        """
        Modify and existing file record (self.uploaded_files)
        Mostly used to change a pending file to normal file path, or to modify the list of plates
        new_name (str)
        new_purpose (str)
        extra_pids ([str])
        PID_name_updates ([(str,str)]) - list of (from,to) tuples.
        """
        if existing_name not in self.uploaded_files:
            self.log(f'Error: {existing_name} does not exist in uploaded file records')
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
                    self.log(f'Error: {existing_pid} not present in uploaded_files')
                    return False
                self.uploaded_files[existing_name]['plates'].remove(existing_pid)
                self.uploaded_files[existing_name]['plates'].append(new_pid)
                                                                                    
        if new_name:
            self.uploaded_files[new_name] = self.uploaded_files[existing_name].copy()
            del self.uploaded_files[existing_name]
        
        return True
    
    
    def del_file_record(self, file_name):
        """
        Remove a file record, and the actual file permanently
        """
        if file_name in self.uploaded_files:
            if Path(file_name).exists():
                try:
                    Path(file_name).unlink()
                except Exception as exc:
                    self.log(f'Error: could not delete file {file_name} {exc}')
                    return False
            del self.uploaded_files[file_name]
            self.log(f'Success: removed file record of {file_name}')
            return True
        else:
            self.log(f'Error: {file_name} not present in file records')
            return False

    ### Plate related operations

    def build_dna_plate_entry(self, sample_plate_ids, dna_plate_id, source=None):
        """
        Replaces the rodentity- and custom-specific code for combining 96-well sample plates into
        384-well DNA plates.
        """
        if self.locked:
            self.log('Error: Cannot add DNA plate set while lock is turned on')
            return False
        
        dna_plate_id = util.guard_pbc(dna_plate_id, silent=True)

        if dna_plate_id in self.plate_location_sample:
            self.log(f'Error: plate {dna_plate_id} already exists! Please delete this plate before trying again')
            return False

        sample_plate_ids = sorted([util.guard_pbc(spid, silent=True) for spid in sample_plate_ids if spid])

        if source == 'rodentity':
            self.log(f'Begin: combining Rodentity plate set into 384-well DNA plate {dna_plate_id}')
            for spid in sample_plate_ids:
                purpose = self.plate_location_sample[spid]['purpose']
                source = self.plate_location_sample[spid]['source']
                if purpose != 'sample' or source != 'rodentity':
                    self.log(f'Error: cannot combine {spid} with {purpose=} and {source=}')
                    return False
        elif source == 'custom':
            self.log(f'Begin: combining custom plate set into 384-well DNA plate {dna_plate_id}')
            for spid in sample_plate_ids:
                purpose = self.plate_location_sample[spid]['purpose']
                source = self.plate_location_sample[spid]['source']
                if purpose != 'sample' or source != 'custom':
                    self.log(f'Error: cannot combine {spid} with {purpose=} and {source=}')
                    return False
        else:
            self.log('Error: must choose "rodentity" or "custom" as input source')
            return False

        self.dest_sample_plates[dna_plate_id] = sample_plate_ids
        self.plate_location_sample[dna_plate_id] = {'purpose':'dna', 'source':','.join(sample_plate_ids), 'wells':set()}
        if source == 'rodentity':
            self.unassigned_plates = {1:'', 2:'', 3:'', 4:''}
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
            pid_count = 0
            for i,sample_pid in enumerate(sorted(self.dest_sample_plates[dna_pid])):
                plate_set_details.append(util.unguard_pbc(sample_pid))
                for well in self.plate_location_sample[sample_pid]['wells']:  
                    if util.is_guarded_cbc(self.plate_location_sample[sample_pid][well]['barcode']):
                        custom_wells += 1
                        total_well_counts['c'] += 1 
                    elif util.is_guarded_rbc(self.plate_location_sample[sample_pid][well]['barcode']):
                        rodentity_wells += 1
                        total_well_counts['r'] += 1   
                    #for assay in info['assays']:
                    #    d['Primers required'].add(assay)
                    #    total_unique_assays.add(assay)
                    total_unique_samples.add(self.plate_location_sample[sample_pid][well]['barcode'])
                #print(f'{dna_pid=} {sample_pid=} {custom_wells=} {rodentity_wells=}')
                pid_count += 1
            for j in range(3-pid_count):
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
                if plate['pcr_stage'] == 1:
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
                fp = self.get_exp_fn(nim_output.name, trans=True)
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
                final_fp = transaction.convert_pending_to_final(self,fp)
                if final_fp in self.uploaded_files:
                    self.log(f'Warning: {final_fp} already recorded as uploaded, overwriting')
                self.uploaded_files[final_fp] = {'plates': list(plate_set), 'purpose':'DNA'} 
        except Exception as exc:
            self.log(f'Error: could not upload Hamilton Nimbus output files, {exc}')
            return False
        transaction.add_pending_transactions(self, transactions)
        self.save()
        return True

                                                             
    def get_dna_pids(self, dna_pids=None):
        """ return a list of available DNA plate ids """
        dpids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'dna']
        if dna_pids:
            dna_pids = [util.guard_pbc(d, silent=True) for d in dna_pids]
            dpids = [d for d in dpids if d in dna_pids]
        return dpids
                    

    def get_dna_records(self, dna_pids=None):
        """
        dna_fields=['samplePlate','sampleWell','sampleBarcode','strain','sex','alleleSymbol',
                 'alleleKey','assayKey','assays','assayFamilies','clientName','sampleName',
                 'dnaPlate','dnaWell','primer']
        """
        #print(f'{dna_pids=}')
        dna_pids = self.get_dna_pids(dna_pids=dna_pids)
        #print(f'{dna_pids=}')
        records = []
        for dna_bc in sorted(dna_pids):
            for well in self.plate_location_sample[dna_bc]['wells']:
                #print(f'{dna_bc=} {well=} {self.plate_location_sample[dna_bc][well]=}')
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
        """ return a list of available PCR plate ids """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'pcr']


    def get_primer_pids(self, pmr_pids=None):
        """ return all primer pids, unless restricted by pmr_pids """
        if pmr_pids:
            primer_pids = [util.guard_pbc(ppid, silent=True) for ppid in pmr_pids \
                    if ppid in self.plate_location_sample and self.plate_location_sample[ppid]['purpose'] == 'primer']
        else:
            primer_pids = [util.guard_pbc(ppid, silent=True) for ppid in self.plate_location_sample \
                    if ppid in self.plate_location_sample and self.plate_location_sample[ppid]['purpose'] == 'primer']
        return primer_pids


    def get_available_primer_wells(self, pmr_pids=None):
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


    def get_available_primer_vols_doses(self, pmr_pids=None):
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


    def get_primers_for_assays(self, assays):
        """
        For each assay in assays, find the matching assay family and the corresponding primers
        Return a dictionary[assay] = [primers]
        """
        assay_primers = {}
        for assay in assays:           
            if assay in assay_primers:
                continue  # we've already seen this assay
            if assay not in self.assay_assayfam:
                continue
            assayfam = self.assay_assayfam[assay]
            if assayfam not in self.assayfam_primers:
                continue
            primers = self.assayfam_primers[assayfam]
            assay_primers[assay] = primers
        return assay_primers

    
    def get_assay_primer_usage(self, dna_pids=None):
        """
        For a given set of DNA plates, return the number of times each assay and primer is needed
        If an dna_pids is None, use all available dna plates
        """
        print(f'{dna_pids=}')
        assay_usage = defaultdict(int)
        primer_usage = defaultdict(int)
        if dna_pids:
            dpids = [pid for pid in dna_pids if pid in self.plate_location_sample \
                    and self.plate_location_sample[pid]['purpose'] == 'dna']
        else:
            dpids = [pid for pid in self.plate_location_sample if pid in self.plate_location_sample \
                    and self.plate_location_sample[pid]['purpose'] == 'dna']
        if not dpids:
            self.log(f'Error: no DNA plates found matching DNA plate IDs: {dna_pids}')
            return assay_usage, primer_usage
        print(f'{dpids=}')

        for dpid in dpids:        
            for well in self.plate_location_sample[dpid]['wells']:       
                assays = self.plate_location_sample[dpid][well]['ngs_assays']
                for assay in assays:  
                    assay_usage[assay] += 1
                assays_primers = self.get_primers_for_assays(assays)
                for a in assays_primers:
                    for primer in assays_primers[a]:
                        primer_usage[primer] += 1
        return assay_usage, primer_usage


    def get_index_avail(self):
        """
        Use all index plates available and return dictionaries of fwd and rev indexes
        Returns:
            [fwd_index]={'well_count':int, 'avail_transfers':int, 'avail_vol':nl}
            [rev_index]={'well_count':int, 'avail_transfers':int, 'avail_vol':nl}
        """
        index_pids = []
        for pid in self.plate_location_sample:
            if self.plate_location_sample[pid]['purpose'] == 'index':
                index_pids.append(pid)
                
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
                    self.log(f'Warning: unexpected index name: {name}')

        return fwd_idx, rev_idx
     
    
    def get_index_reactions(self, primer_usage, fwd_idx=None, rev_idx=None):
        """
        Returns the barcode pairs remaining, max available barcode pairs
        Can take the output of get_index_avail to save recalculating them
        """
        reactions_required = sum([v for v in primer_usage.values()])

        if not fwd_idx or not rev_idx:
            fwd_idx, rev_idx = self.get_index_avail()

        if not fwd_idx or not rev_idx:
            return 0, 0

        reactions_possible = 0
        max_i7F = len(fwd_idx)
        max_i5R = len(rev_idx)
        for fwd in fwd_idx:
            if fwd_idx[fwd]['avail_transfers'] < 1:
                continue
            for rev in rev_idx:
                if fwd_idx[fwd]['avail_transfers'] < 1:
                    break
                if rev_idx[rev]['avail_transfers'] < 1:
                    continue
                reactions_possible += 1
                fwd_idx[fwd]['avail_transfers'] -= 1
                rev_idx[rev]['avail_transfers'] -= 1

        return reactions_possible - reactions_required, reactions_possible


    def get_taqwater_avail(self, taqwater_bcs=None, transactions=None, pcr_stage=None):
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
                water_avail += tw_plate[well]['volume'] - self.dead_volumes[util.PLATE_TYPES['Echo6']]
                if transactions is not None:
                    if pid in transactions:
                        if well in transactions[pid]:
                            water_avail += transactions[pid][well] # the transactions are recorded as change, so will be negative values

            for well in ['B1','B2','B3']:
                taq_avail += tw_plate[well]['volume'] - self.dead_volumes[util.PLATE_TYPES['Echo6']]
                if transactions is not None:
                    if pid in transactions:
                        if well in transactions[pid]:
                            taq_avail += transactions[pid][well]
        return taq_avail, water_avail, pids


    def get_taq_water_volumes_required(self, num_reactions):
        """
        Returns the taq and water volumes required for the primer and index stages in nl
        """
        primer_taq_vol = num_reactions * self.transfer_volumes['PRIMER_TAQ_VOL']
        primer_water_vol = num_reactions * self.transfer_volumes['PRIMER_WATER_VOL']
        index_water_vol = num_reactions * self.transfer_volumes['INDEX_WATER_VOL']
        index_taq_vol = num_reactions * self.transfer_volumes['INDEX_TAQ_VOL'] 
        return primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol


    def get_taqwater_pids(self, pcr_stage=None):
        """
        Return plate IDs for taq/water plates.
        If PCR stage is provided, return the taq/water plate for that stage. 
        pcr_stage: 1 (primer) or 2 (index). If None, return all plate IDs for taq/water plates
        """
        pids = []
        #Plate ID specific to PCR stage (1 or 2). Otherwise get all taq/water plates. 
        if pcr_stage:
            for p in self.plate_location_sample:
                if self.plate_location_sample[p]['purpose'] == 'taq_water':
                    if self.plate_location_sample[p]['pcr_stage'] == pcr_stage:
                        pids.append(p)       
        else:
            pids = [p for p in self.plate_location_sample
                    if self.plate_location_sample[p]['purpose'] == 'taq_water']
        return pids


    def generate_nimbus_inputs(self):
        success = generate.nimbus_gen(self)
        return success


    def get_nimbus_filepaths(self):
        """ Return the lists of nimbus input files, echo input file (nimbus outputs), 
            and barcodes that are only seen in nimbus """
        #print("In get_nimbus_filepaths")
        nimbus_input_filepaths, echo_input_paths, xbc = generate.match_nimbus_to_echo_files(self)
        return nimbus_input_filepaths, echo_input_paths, xbc


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
        primer_survey_fn = self.get_exp_fn(primer_survey_filename, trans=True)
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
        transaction.add_pending_transactions(self,transactions)
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
        success = generate.generate_echo_PCR1_picklist(self, dna_plates, pcr_plates, taq_water_plates)
        if not success:
            self.log('Failure: could not generate PCR1 picklists correctly')
            return False
        self.log('Success: generated PCR1 picklists')
        return True

    def get_echo_PCR1_picklist_filepaths(self):
        """
        Return file paths for PCR1_dna-picklist_XXX.csv, PCR1_primer-picklist_XXX.csv, PCR1_taqwater-picklist_XXX.csv
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
        success = generate.generate_echo_PCR2_picklist(self, pcr_plates, index_plates, taq_water_plates, amplicon_plates)
        if not success:
            self.log('Failure: could not generate PCR2 (index) picklists correctly')
            return False
        self.log('Success: generated PCR2 (index) picklists')
        return success

    def get_echo_PCR2_picklist_filepaths(self):
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

    def delete_plates(self, pids):
        """
        Soft-delete plates with the selected pids
        """
        print(f'In delete_plates with {pids=}')
        success = True
        for pid in pids:
            if pid in self.plate_location_sample:
                if pid not in self.deleted_plates:
                    self.deleted_plates[pid] = []
                self.deleted_plates[pid].append(deepcopy(self.plate_location_sample[pid]))
                # need to find this in reproducible steps and delete
                del self.plate_location_sample[pid]
                self.log(f'Warning: moved {pid} to deleted plate bin')
                if pid in self.dest_sample_plates:
                    del self.dest_sample_plates[pid]
                    self.log(f'Warning: removed 384-well DNA plate entry {pid}')
            else:
                self.log(f'Warning: {pid} has no definition loaded')
                success = False
        self.save()
        return success

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
                if type(self.pending_steps[file_name]) == 'dict':
                    pids = [util.unguard_pbc(pid, silent=True) for pid in self.pending_steps[file_name].keys()]
                if type(self.pending_steps[file_name]) == 'list':
                    pids = []
                    for thing in self.pending_steps[file_name]:
                        for pid in thing:
                            pids.append(util.unguard_pbc(pid, silent=True))
                else:
                    pids = [util.unguard_pbc(pid, silent=True) for pid in self.pending_steps[file_name]]
                stages.append([str(counter), file_name, ', '.join(pids), 'pending'])
        return stages, header

    def get_stage2_pcr_plates(self):
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

    def get_file_usage(self):
        """
        Get files in uploaded files that have a purpose or plates. If a file is a reference file, but has been overwritten
        and not being currently used by the experiment, do not include it. 
        """
        file_usage = {}
        for filename in self.uploaded_files:
            #only inclue current reference file using for reference_sequences
            if self.uploaded_files[filename].get('purpose') == 'reference_sequences' \
                        and filename not in self.reference_sequences:
                continue

            file_usage[filename] = {'plates' :[], 'purpose':None}
            if 'purpose' in self.uploaded_files[filename]:
                file_usage[filename]['purpose'] = self.uploaded_files[filename]['purpose']
            
            if 'plates' in self.uploaded_files[filename]:
                for pid in self.uploaded_files[filename]['plates']:
                    if util.is_guarded(pid):
                        pid = util.unguard_pbc(pid)
                    file_usage[filename]['plates'].append(pid)

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
        fn = self.get_exp_fn(index_survey_filename)
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

    def add_standard_taqwater_plates(self, plate_barcodes, pcr_stage):  # <- they're the same for primer and barcode stages
        """ These plates have a fixed layout with water in row A and taq in row B of a 6-well reservoir plate.
        See echo_primer.py or echo_barcode.py mytaq2()
        We may need to also store volumes at some stage, but for now that isn't necessary
        pcr_stage(int): which PCR experiment the plate comes from (1 or 2)
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
                    elif self.plate_location_sample['pcr_stage'] == pcr_stage:
                        self.log(f"Error: cannot add taq+water plate {util.unguard(pid)}."+\
                            f"A taq/waterplate with this barcode already exists with for PCR {pcr_stage}")
                        return False
                    else:
                        self.log(f"Warning: overwriting existing taq/water plate entry with barcode {util.unguard(pid)}")
                else:
                    self.plate_location_sample[pid] = {}
                self.plate_location_sample[pid]['purpose'] = 'taq_water'
                self.plate_location_sample[pid]['plate_type'] = util.PLATE_TYPES['Echo6']
                self.plate_location_sample[pid]['pcr_stage'] = pcr_stage
                self.plate_location_sample[pid]['capacity'] = self.cap_volumes[util.PLATE_TYPES['Echo6']]
                self.plate_location_sample[pid]['dead_vol'] = self.dead_volumes[util.PLATE_TYPES['Echo6']]
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
        miseq_fps = list(Path(self.get_exp_dn()).glob('MiSeq_*.csv'))
        return miseq_fps  


    def check_sequence_upload_ready(self, messages):
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
        dn = self.get_exp_fn(f'raw')
        if not Path(dn).exists() or not Path(dn).is_dir():
            success = False
            msg = f"Error: {dn} does not exist"
            self.log(msg)
            messages.append(msg)

        self.save()
        return success


    #Gabi's code for custom volumes
    def add_custom_volumes(self, custom_volumes:dict) -> bool:
        """
        Edit volumes for DNA, primers, indexes and taq/water.
        Args:
            custom_volumes(dict): Volumes for transfer_volumes
        Return:
            True if successfully changed the transfer volumes. False is there is an error with the type (not an int)
        """
        try:
            for value in list(custom_volumes.values()):
                value = float(value)
            self.transfer_volumes = custom_volumes.copy()
            self.log(f'Success: Custom volumes added')
            self.save()
            return True
        except (TypeError, ValueError) as e:
            self.log(f'Error: {e}')
            return False
    #End of Gabi's code for custom volumes

    def add_custom_plate_volumes(self, custom_plate_vols:dict):
        """
        Update the experiment's dead and cap volumes with the custom volumes. Also update any plates that have been added
        with a dead and cap volume. 
        Args:
            custom_plate_vols(dict): New volumes from user to replace the current dead and cap volumes. Dictionary in 
                                    the form: {'plate_type': [dead_volume, cap_volume]}

        Return:
            True once edits have been added. Returns False if the volumes in the custom_plate_vols aren't numbers.
        """
        try:
            #check the entries are floats
            custom_plate_vols = {plate: [float(vol) for vol in vols] for plate, vols in custom_plate_vols.items()}
        except (TypeError, ValueError) as e:
            self.log(f'Error: {e}')
            return False

        #update dead and cap volumes with the custom volumes
        self.dead_volumes.update({plate: vols[0] for plate, vols in custom_plate_vols.items() \
                                    if plate in self.dead_volumes and vols})
        
        self.cap_volumes.update({plate: vols[1] for plate, vols in custom_plate_vols.items() \
                                    if plate in self.cap_volumes and len(vols) > 1})

        # Update plate_location_sample if any plates have a dead and cap volume
        for plate, info in self.plate_location_sample.items():
            if info.get('dead_vol'):
                info['dead_vol'] = self.dead_volumes[info['plate_type']]
            if info.get('capacity'):
                print(plate,info)
                info['capacity'] = self.cap_volumes[info['plate_type']]

        self.log(f'Success: Custom volumes added')
        self.save()
        return True

    def read_allele_results(self):
        """
        Read in allele_results
        [SampleNo, EPplate, EPwell, sampleBarcode, assays, assayFamilies, dnaPlate, dnaWell, 
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
        if (now - self.log_time).seconds > 10 or level in ['Error','Critical','Failure','End']:
            self.save()

    def get_log_header(self):
        return ['Time', 'Function name', 'Func line', 'Calling function', 'Call line', 'Level', 'Message']

    def get_log(self, num_entries=100):
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