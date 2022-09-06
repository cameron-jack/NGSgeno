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

import pandas as pd
from streamlit import _update_logger

## Ugly hack to find things in ../bin
#from pathlib import Path
#file = Path(__file__).resolve()
#sys.path.append(os.path.join(file.parents[1],'bin'))
import bin.db_io as db_io
import bin.file_io as file_io
from bin.util import CAP_VOLS, DEAD_VOLS, output_error, padwell, unpadwell, row_ordered_96, row_ordered_384, calc_plate_assay_usage
from bin.util import DNA_VOL, PRIMER_VOL, PRIMER_TAQ_VOL, PRIMER_WATER_VOL, BARCODE_VOL, BARCODE_TAQ_VOL, BARCODE_WATER_VOL
from bin.nimbus import read_custom_csv_to_json, nimbus_gen
from bin.echo_primer import generate_echo_PCR1_picklist

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
    Info must have 'purpose' in {'sample', 'primer', 'dna', 'pcr', 'taq_water'}
    We need functions to link results back to sample
    samples are their own dict structure, separate from plates, which have the form plate_locations[pid] = {well:sample_id}

    wells are always unpadded. Echo needs unpadded wells, possibly Miseq?

    Do we read raw inputs and then offer the ability to manually correct them? Maybe only the manifests...
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
        self.reproducible_steps = []  # a strict subset of the log that only includes working instructions
        self.log_entries = []  # use self.log(message, level='') to add to this
        self.transfer_volumes = {'DNA_VOL':DNA_VOL, 'PRIMER_VOL':PRIMER_VOL, 'PRIMER_TAQ_VOL':PRIMER_TAQ_VOL,
                'PRIMER_WATER_VOL':PRIMER_WATER_VOL, 'BARCODE_VOL':BARCODE_VOL, 'BARCODE_TAQ_VOL':BARCODE_TAQ_VOL,
                'BARCODE_WATER_VOL':BARCODE_WATER_VOL}  # load the defaults from util.py
        
        #self.primer_plate_names = {}  # {PID:{well:name}}
        #self.primer_plate_volumes = {}  # {PID:{well:volume}}
        #self.primer_taqwater_plates_names = {}  # {PID:{well:name}}
        #self.primer_taqwater_plate_volumes = {}  # {PID:{well:volume}}
        #self.barcode_plate_names = {}  # {PID:{well:name}}
        #self.barcode_plate_volumes = {}  # {PID:{well:volume}}
        #self.barcode_taqwater_plate_names = {}  # {PID:{well:name}}
        #self.barcode_taqwater_plate_volumes = {}  # {PID:{well:volume}}
        

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

    def get_exp_dir(self):
        """ return the experiment directory name """
        return 'run_' + self.name

    def get_exp_fp(self, filename):
        """ return the expected experiment path to filename """
        dirname = self.get_exp_dir()
        if dirname in filename:
            return fp
        fp = os.path.join(dirname, filename)
        return fp

   
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
            
            sample_plate_ids = sorted([file_io.guard_pbc(spid, silent=True) for spid in sample_plate_ids if spid])
        
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
                    pos = unpadwell(record['wellLocation'])
                    self.plate_location_sample[spid]['wells'].add(pos)
                    if pos not in self.plate_location_sample[spid]:
                        self.plate_location_sample[spid][pos] = {}
                    self.plate_location_sample[spid][pos]['barcode'] = file_io.guard_rbc(record['mouse']['barcode'], silent=True)
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
                    pos = unpadwell(record['wellLocation'])
                    self.plate_location_sample[spid]['wells'].add(pos)
                    self.plate_location_sample[spid][pos] = record['mouse'].copy()  # dict
                    self.plate_location_sample[spid][pos]['mouseId'] = record['mouseId']
                    self.plate_location_sample[spid][pos]['barcode'] = file_io.guard_mbc(record['mouseBarcode'])
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
                        manifest_info[i][key] = file_io.guard_pbc(entry[key], silent=True)
                        dpids_spids[manifest_info[i][key]] = set()
                    if key == 'Plate barcode':
                        manifest_info[i][key] = file_io.guard_pbc(entry[key], silent=True)
                    if key == 'Sample barcode':
                        if not file_io.is_guarded(entry[key]):
                            if default_manifest_type == 'c':
                                manifest_info[i][key] = file_io.guard_cbc(entry[key])
                            elif default_manifest_type == 'r':
                                manifest_info[i][key] = file_io.guard_rbc(entry[key])
                            elif default_manifest_type == 'm':
                                manifest_info[i][key] = file_io.guard_mbc(entry[key])
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
                well = unpadwell(entry['Well'])
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
                if file_io.is_guarded_mbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    musterer_pids.add(pid)
        return musterer_pids

    def get_rodentity_pids(self):
        """ return [pids] with samples that are sourced from Rodentity """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        rodentity_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if file_io.is_guarded_rbc(self.plate_location_sample[pid][well]['sample_barcode']):
                    rodentity_pids.add(pid)
        return rodentity_pids

    def get_custom_pids(self):
        """ return [pids] with samples that are sourced as custom """
        # Plate contents[barcode] : { well_location(row_col):{sample_barcode:str, strain:str, assays: [], possible_gts: [], 
        #   other_id: str, parents: [{barcode: str, sex: str, strain: str, assays = [], gts = []}] } }
        custom_pids = set()
        for pid in self.plate_location_sample:
            for well in self.plate_location_sample[pid]:
                if file_io.is_guarded_cbc(self.plate_location_sample[pid][well]['sample_barcode']):
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
            d = {'DNA PID': file_io.unguard_pbc(dna_pid)}
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
                d['Sample PID'+str(i+1)] = file_io.unguard_pbc(sample_pid)
                for info in (self.plate_location_sample[sample_pid][well] for well in self.plate_location_sample[sample_pid]['wells']):
                    samp_barcode = info['barcode']     

                    if file_io.is_guarded_cbc(samp_barcode):
                        d['Custom samples'] += 1
                        total_well_counts['c'] += 1
                    #elif file_io.is_guarded_mbc(samp_barcode):
                    #    d[' counts']['m'] += 1
                    #    total_well_counts['m'] += 1  
                    elif file_io.is_guarded_rbc(samp_barcode):
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

    def get_valid_assay_counts(self):
        """ valid assays file must not be empty """
        pass

    def get_invalid_assay_counts(self):
        """ valid assays file must not be empty """
        pass

    def get_free_barcode_count(self):
        """ calculate the remaining free i7i5 barcodes available """
        pass

    def add_pcr_plates(self, pcr_plate_list=[]):
        """
        Add one or more empty 384-well plates to self.plate_location_sample
        """
        for pid in pcr_plate_list:
            p = file_io.guard_pbc(pid, silent=True)
            if p in self.plate_location_sample:
                self.log(f'Error: {p} already in use, skipping')
                continue
            self.plate_location_sample[p] = {'purpose':'pcr', 'wells':set(), 'source':'user', 'plate_type':'384PP_AQ_BP'}
        self.save()
        return True
                                                                                            
    def get_pcr_plates(self):
        """ return a list of user supplied PCR plate ids """
        return [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'pcr']

    def get_assay_usage(self, dna_plate_list=[]):
        """ We will want ways to get assay usage from various subsets, but for now do it for everything that 
            has a Nimbus destination plate.
        """
        assay_usage = Counter()
        for dest in self.dest_sample_plates:
            if dna_plate_list and dest not in dna_plate_list:
                continue
            for samp in self.dest_sample_plates[dest]:
                assay_usage.update(calc_plate_assay_usage(self.plate_location_sample[samp]))
        return assay_usage

    def get_volumes_required(self, assay_usage=None, dna_plate_list=[]):
        """
            Standard usage is to call get_assay_usage() first and pass it to this.
        """
        if not assay_usage:
            assay_usage = self.get_assay_usage(dna_plate_list=dna_plate_list)
        reactions = sum([v for v in assay_usage.values()])

        primer_taq_vol = reactions * self.transfer_volumes['PRIMER_TAQ_VOL']
        primer_water_vol = reactions * self.transfer_volumes['PRIMER_WATER_VOL']
        barcode_water_vol = reactions * self.transfer_volumes['BARCODE_WATER_VOL']
        barcode_taq_vol = reactions * self.transfer_volumes['BARCODE_TAQ_VOL']

        primer_vols = {a:assay_usage[a]*self.transfer_volumes['PRIMER_VOL'] for a in assay_usage}
        return primer_vols, primer_taq_vol, primer_water_vol, barcode_taq_vol, barcode_water_vol
        
        
    def get_barcode_remaining_available_volume(self, assay_usage=None):
        """
        Returns the barcode pairs remaining, max available barcode pairs, max barcode pairs allowed by volume
        """
        barcode_pids = []
        for pid in self.plate_location_sample:
            if self.plate_location_sample[pid]['purpose'] == 'i7i5barcodes':
                barcode_pids.append(pid)
        if not barcode_pids:
            return 0, 0, False

        if not assay_usage:
            assay_usage = self.get_assay_usage()
        reactions = sum([v for v in assay_usage.values()])
        
        fwd_barcode_vols = {}  # name=[vol, vol, ...]
        rev_barcode_vols = {}
        
        for bpid in barcode_pids:
            for well in self.plate_location_sample[bpid]['wells']:
                if 'idt_name' not in self.plate_location_sample[bpid][well]:
                    continue
                name = self.plate_location_sample[bpid][well]['idt_name']
                if 'i7F' in name:
                    if name not in fwd_barcode_vols:
                        fwd_barcode_vols[name] = []
                    if 'volume' in self.plate_location_sample[bpid][well]:
                        fwd_barcode_vols[name].append(self.plate_location_sample[bpid][well]['volume'])
                elif 'i5R' in name:
                    if name not in rev_barcode_vols:
                        rev_barcode_vols[name] = []
                    if 'volume' in self.plate_location_sample[bpid][well]:
                        rev_barcode_vols[name].append(self.plate_location_sample[bpid][well]['volume'])
                else:
                    self.log('Unexpected i7i5 barcode name:' + name, level='Warning')
        max_i7F = len(fwd_barcode_vols)
        max_i5R = len(rev_barcode_vols)
        reaction_vol_capacity = 0
        for name in fwd_barcode_vols:
            # get the number of possible reactions and keep the lower number from possible reactions or possible reaction partners
            max_reactions = sum([floor((vol-DEAD_VOLS['384PP_AQ_BP'])/self.transfer_volumes['BARCODE_VOL']) for vol in fwd_barcode_vols[name]])
            reaction_vol_capacity += min(max_reactions, max_i5R)
        max_barcode_pairs = max_i7F * max_i5R

        return max_barcode_pairs-reactions, max_barcode_pairs, reaction_vol_capacity


    def get_primers_avail(self):
        """ returns {primer:count}, {primer:vol} from what's been loaded """
        primer_counts = {}
        primer_vols = {}
        for pid in self.plate_location_sample:
            if self.plate_location_sample[pid]['purpose'] == 'primers':
                for well in self.plate_location_sample[pid]['wells']:
                    if 'primer' in self.plate_location_sample[pid][well]:
                        primer = self.plate_location_sample[pid][well]['primer']
                        if primer not in primer_counts:
                            primer_counts['primer'] = 0
                        if primer not in primer_vols:
                            primer_vols['primer'] = 0
                        primer_counts['primer'] += 1
                        if 'volume' in self.plate_location_sample[pid][well]:
                            primer_vols['primer'] += self.plate_location_sample[pid][well]['volume']
        return primer_counts, primer_vols

    def get_taqwater_avail(self):
        """ returns (int) taq and (int) water volumes loaded as available for a given Echo stage (PCR) in nanolitres """
        taq_avail = 0
        water_avail = 0
        pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'taq_water']
        for pid in pids:
            taq_avail += 3* (CAP_VOLS['6RES_AQ_BP2'] - DEAD_VOLS['6RES_AQ_BP2'])
            water_avail += 3* (CAP_VOLS['6RES_AQ_BP2'] - DEAD_VOLS['6RES_AQ_BP2'])
        return taq_avail, water_avail, pids

    def get_taqwater_plates(self):
        pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'taq_water']
        return pids

    def generate_nimbus_inputs(self):
        success = nimbus_gen(self)
        #print("In generate_nimbus_inputs", success)
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
                print(f"Whoops, {row['DNA PID']=} doesn't actually exist in the experiment!")
                continue
            self.dest_sample_plates.pop(row['DNA PID'])
            #for field in row:
                #if field.startswith('Sample PID') and row[field] in self.plate_location_sample:
                    # Need to check whether this sample plate is used elsewhere, because otherwise we can't get rid of it.
                    # Can we just leave it in here?
                #    self.plate_location_sample.pop(row[field])
            if row['DNA PID'] in self.plate_location_sample:
                self.plate_location_sample.pop(row['DNA PID'])
            #self.plate_source.pop(row['DNA PID'])
            # This also has implications for later stages. Ignore this for now
        self.save()
        return True

    def add_references(self, uploaded_references):
        if self.locked:
            self.log('Error: Cannot add reference sequences while lock is active')
            return False
        for uploaded_reference in uploaded_references:
            if uploaded_reference.name in self.reference_sequences:
                # do we warn? Is this an error?
                self.log(f"Error: Duplicate reference file name: {uploaded_reference.name=} in . Overwriting")
            seq_itr = SeqIO.parse(uploaded_reference.getvalue().decode("utf-8"), "fasta")
            self.reference_sequences[uploaded_reference.name] = {s.id:str(s.seq) for s in seq_itr}
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
            PID = uploaded_primer_plate.name.split('_')[0]  # assumes first field is PID
            if PID in self.plate_location_sample:
                if self.plate_location_sample[PID]['purpose'] == 'primers':
                    self.log(f"Info: Primer file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Primer file PID: {uploaded_primer_plate=} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return
            else:  # create new plate entry and set purpose
                self.plate_location_sample[PID] = {'purpose':'primers', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
                self.log(f"Info: Creating new primer plate record for {PID}")
            # load data into plate
            data = StringIO(uploaded_primer_plate.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0:
                    continue  # header
                if row == '' or row[0] == '' or row[1] == '':
                    continue
                well = unpadwell(row[0])
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
            PID = uploaded_primer_volume.name.split('_')[0]  # assumes first field is PID
            if PID in self.plate_location_sample:
                if self.plate_location_sample[PID]['purpose'] == 'primers':
                    self.log(f"Info: Primer file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Primer file PID: {uploaded_primer_volume=} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return
            else:  # create new plate entry and set purpose
                self.plate_location_sample[PID] = {'purpose':'primers', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
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
                    well = unpadwell(row[0])
                    if well not in self.plate_location_sample[PID]:
                        self.plate_location_sample[PID]['wells'].add(well)
                        self.plate_location_sample[PID][well] = {}
                    self.plate_location_sample[PID][well]['volume'] = float(row[1])*1000
        self.log(f"Success: added primer volumes from {', '.join([upv.name for upv in uploaded_primer_volumes])}")
        self.save()
        return True

    def generate_echo_primer_survey(self, primer_survey_filename='primer-svy.csv'):
        """ 
        Generate a primer survey file for use by the Echo. Replaces echovolume.py
        Not strictly necessary, but the old code reads this file in making the picklists.
        """
        if self.locked:
            self.log('Error: cannot generate primer survey file while lock is active.')
            return False
        primer_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'primers']
        header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'plate position on Echo 384 PP',
                'primer names pooled', 'volume']
        fn = self.get_exp_fp(primer_survey_filename)
        if os.path.exists(fn):
            self.log(f'Warning: overwriting primer survey file {fn}')
        with open(fn, 'wt') as fout:
            print(','.join(header), file=fout)
            for i,pid in enumerate(primer_pids):
                plate=self.plate_location_sample[pid]
                for well in row_ordered_384:
                    # TODO: correct the columns below
                    if well not in plate['wells']:
                        continue
                    if 'primer' not in plate[well]:
                        continue
                    if 'volume' not in plate[well]:
                        continue
                    outline = ','.join([f"Source[{i+1}]",file_io.unguard_pbc(pid,silent=True),plate['plate_type'],
                            well,plate[well]['primer'],f"{int(plate[well]['volume'])/1000}"])
                    print(outline, file=fout)
        self.log(f"Success: written Echo primer survey to {fn}")
        return True

    def check_ready_echo1(self, dna_plates, pcr_plates, taq_water_plates):
        """
        TODO: Should check that everything required to successfully generate PCR1 picklists is available
        """
        return True

    def generate_echo_PCR1_picklist_interface(self, dna_plates, pcr_plates, taq_water_plates):
        """
        Calls echo_primer.generate_echo_PCR1_picklist() to do the work, needs a set of accepted DNA_plates,
        the final assay list, and a set of destination PCR plate barcodes, taq+water plates, primer plates, primer volumes.
        Returns True on success
        """
        print(f"Experiment.generate_echo_PCR1 {dna_plates=} {pcr_plates=} {taq_water_plates=}")
        for pid in pcr_plates:
            if pid in self.plate_location_sample:
                if self.plate_location_sample[pid]['purpose'] != 'pcr':
                    self.log(f"Error: plate already exists with PID {pid} with purpose {self.plate_location_sample[pid]['purpose']}")
                    return False
            self.log(f'Warning: existing entry for PCR plate {pid}. This will be overwritten!')
        success = self.generate_echo_primer_survey()
        if not success:
            self.log('Failure: failed to generate primer survey file')
            return False
        success = generate_echo_PCR1_picklist(self, dna_plates, pcr_plates, taq_water_plates)

    def get_echo_PCR1_picklist_filepaths(self):
        """
        Return file paths for PCR1_dna-picklist_XXX.csv, PCR1_primer-picklist_XXX.csv, PCR1_taqwater-picklist_XXX.csv
        """
        all_files = os.listdir(self.get_exp_dir())
        dna_picklist_path = None
        primer_picklist_path = None
        taqwater_picklist_path = None
        for f in all_files:
            if f.startswith('PCR1_dna-picklist_') and f.endswith('.csv'):
                dna_picklist_path = self.get_exp_fp(f)
            elif f.startswith('PCR1_primer-picklist_') and f.endswith('.csv'):
                primer_picklist_path = self.get_exp_fp(f)
            elif f.startswith('PCR1_taqwater-picklist_') and f.endswith('.csv'):
                taqwater_picklist_path = self.get_exp_fp(f)
        return dna_picklist_path, primer_picklist_path, taqwater_picklist_path


    def delete_plates(self, pids):
        """
        Soft-delete plates with the selected pids
        """
        for pid in pids:
            if pid in self.plate_location_sample:
                if pid not in self.deleted_plates:
                    self.deleted_plates[pid] = []
                self.deleted_plates[pid].append(deepcopy(self.plate_location_sample[pid]))
                del self.plate_location_sample[pid]
            self.log(f'Warning: moved {pid} to deleted plate bin')

    def add_barcode_layouts(self, uploaded_barcode_plates):
        """ add barcode plates with well and barcode columns """
        if self.locked:
            self.log('Error: cannot add barcode plates while lock is active.')
            return False
        for uploaded_barcode_plate in uploaded_barcode_plates:
            PID = uploaded_barcode_plate.name.split('_')[0]  # assumes first field is PID
            if PID in self.plate_location_sample:
                if self.plate_location_sample[PID]['purpose'] == 'i7i5barcodes':
                    self.log(f"Info: Barcode file {PID} already exists, adding data")
                else:
                    self.log(f"Error: Barcode file PID: {uploaded_barcode_plate=} matches "+\
                        f"existing plate entry of different purpose {self.plate_location_sample[PID]}")
                    return
            else:  # create new plate entry and set purpose
                self.plate_location_sample[PID] = {'purpose': 'i7i5barcodes', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
                self.log(f"Info: Creating new barcode plate record for {PID}")
            # load data into plate
            data = StringIO(uploaded_barcode_plate.getvalue().decode("utf-8"), newline='')
            for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
                if i == 0 or row == '':
                    continue  # header or blank
                well = unpadwell(row[0])
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
                    
        self.log(f"Success: added barcode plate layouts from {', '.join([ubp.name for ubp in uploaded_barcode_plates])}")        
        self.save()
        return True

    def add_barcode_volumes(self, uploaded_barcode_volumes):
        """ add volume information to a single barcode plate, with a matched, existing plate name.
        To make life fun there are two possible formats:
        Format 1: plate layout (generated through Echo test software)
        Format 2: column layout (generated through Echo main interface)
        Returns True on success
        """
        if self.locked:
            self.log('Error: cannot add barcode plate volumes while lock is active.')
            return False
        for uploaded_barcode_volume in uploaded_barcode_volumes:
            PID = uploaded_barcode_volume.name.split('_')[0]  # assumes first field is PID
            if PID not in self.plate_location_sample:
                self.log(f"Info: Adding barcode plate {PID}")
                self.plate_location_sample[PID] = {'purpose':'i7i5barcodes', 'source':'user', 'wells':set(), 'plate_type':'384PP_AQ_BP'}
            else:
                if self.plate_location_sample[PID]['purpose'] != 'i7i5barcodes':
                    self.log(f"Error: {PID} plate purpose is "+\
                            f"{self.plate_location_sample[PID]['purpose']}, expected 'i7i5barcodes'")
                    return False
                self.log(f"Info: {PID} exists, appending barcode volumes")
        
            plate_format = False
            data = StringIO(uploaded_barcode_volume.getvalue().decode("utf-8"), newline='')
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
                        well = unpadwell(row[3])
                        self.plate_location_sample[PID]['wells'].add(well)
                        if well not in self.plate_location_sample[PID]:
                            self.plate_location_sample[PID][well] = {}
                        self.plate_location_sample[PID][well]['volume'] = int(float[5])*1000
                self.plate_location_sample[PID]['barcode'] = PID
        self.log(f"Success: added barcode plate volumes from {', '.join([ubv.name for ubv in uploaded_barcode_volumes])}")
        self.save()
        return True

    def generate_echo_i7i5_survey(self, i7i5_survey_filename='i7i5-svy.csv'):
        """ 
        Generate an i7i5 survey file for use by the Echo. Replaces echovolume.py 
        Not strictly necessary, but the old code reads this file in making the picklists.
        """
        if self.locked:
            self.log('Error: cannot generate i7i5 barcode survey while lock is active.')
            return False
        i7i5_pids = [p for p in self.plate_location_sample if self.plate_location_sample[p]['purpose'] == 'i7i5barcodes']
        header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'Source well', 
                'Name for IDT ordering', 'index', 'Name', 'Oligo * -phosphothioate modif. against exonuclease', 'volume']
        fn = self.get_exp_fp(i7i5_survey_filename)
        if os.path.exists(fn):
            self.log(f'Warning: overwriting i7i5 barcode survey file {fn}')
        with open(fn, 'wt') as fout:
            print('\t'.join(header)+'\n', file=fout)
            for i,pid in enumerate(i7i5_pids):
                plate=self.plate_location_sample[pid]
                for well in row_ordered_384:
                    if well not in plate['wells']:
                        continue
                    outline = '\t'.join([f'Source[{i}]',file_io.unguard_pbc(pid,silent=True),plate['plate_type'],
                            well,plate[well]['idt_name'],plate[well]['index'],plate[well]['bcnode'],
                            plate[well]['oligo'],plate[well]['volume']])
                    print(outline, file=fout)
        self.log(f'Success: i7i5 survey file written to {fn}')
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
                pid = file_io.guard_pbc(plate_barcode, silent=True)
                if pid in self.plate_location_sample:
                    if self.plate_location_sample[pid]['purpose'] != 'taq_water':
                        self.log(f"Error: cannot add taq+water plate {file_io.unguard(pid)}. "+\
                            f"A plate with this barcode already exists with purpose {self.plate_location_sample[pid]['purpose']}")
                        return False
                    else:
                        self.log(f"Warning: overwriting existing taq/water plate entry with barcode {file_io.unguard(pid)}")
                else:
                    self.plate_location_sample[pid] = {}
                self.plate_location_sample[pid]['purpose'] = 'taq_water'
                self.plate_location_sample[pid]['wells'] = ['A1','A2','A3','B1','B2','B3']
                self.plate_location_sample[pid]['A1'] = {'name': 'water', 'volume': CAP_VOLS['6RES_AQ_BP2']}
                self.plate_location_sample[pid]['A2'] = {'name': 'water', 'volume': CAP_VOLS['6RES_AQ_BP2']}
                self.plate_location_sample[pid]['A3'] = {'name': 'water', 'volume': CAP_VOLS['6RES_AQ_BP2']}
                self.plate_location_sample[pid]['B1'] = {'name': 'taq', 'volume': CAP_VOLS['6RES_AQ_BP2']}
                self.plate_location_sample[pid]['B2'] = {'name': 'taq', 'volume': CAP_VOLS['6RES_AQ_BP2']}
                self.plate_location_sample[pid]['B3'] = {'name': 'taq', 'volume': CAP_VOLS['6RES_AQ_BP2']}
        except Exception as exc:
            self.log(f"Error: adding taq+water plate barcode failed {plate_barcodes=}")
            return False
        self.save()
        return True

     
    def save(self):
        """ save experiment details to self.name/experiment.json, returns True on success and False on fail. Assumes the correct working directory """
        try:
            exp = jsonpickle.encode(self, indent=4, keys=True, warn=True)
            if not self.name:
                return False
            with open(os.path.join('run_'+self.name, EXP_FN), 'wt') as f:
                print(exp, file=f)
        except Exception as exc:
            output_error(exc, msg=f"Error saving {self.name=}")
            return False
        return True

    def log(self, message, level=''):
        """ Always add the date/time, function where the log was run and the caller function to the log """
        now = datetime.datetime.now()
        t = f"{now:%Y-%m-%d %H:%M}"
        func = sys._getframe(1).f_code.co_name
        caller = sys._getframe(2).f_code.co_name
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
        self.log_entries.append([t, func, caller, level, message])

    def get_log(self, num_entries=10):
        """ return a chunk of the log. -1 gives everything """
        if num_entries == -1:
            return self.log[::-1]
        return self.log_entries[:-(num_entries+1):-1]



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
