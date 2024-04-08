#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import gzip
import csv
from pathlib import Path
import traceback
from collections import Counter, OrderedDict
import itertools
import json
from copy import deepcopy
from io import StringIO, BytesIO
from shutil import copyfileobj
from string import ascii_uppercase
import functools

import openpyxl
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    import bin.util as util
except ModuleNotFoundError:
    import util

try:
    import bin.transaction as transaction
except ModuleNotFoundError:
    import transaction

transact = transaction.transact
untransact = transaction.untransact


"""
@created: Nov 2021
@author: Cameron Jack, Bob Buckley, ANU Bioinformatics Consultancy, 2019-2021

Supports:
 * Hamilton Nimbus generate input and read output files
 * Labcyte Echo generate input files (picklists) for primer PCR reactions and index PCR reactions
 * Illumina Miseq generate samplesheet input files

Contains all file IO methods that involved barcodes, that are not purely for user reporting or robot inputs. 
Trivial file reads which don't involve barcodes should stay with the rest of their functionality for better readability

The following rules must be followed to protect users from typing mistakes and MS Excel weirdness:
* Barcodes should be guarded when first read from a new source
* Barcodes should be written with guards if the file is to be read again by the pipeline
* The pipeline should look for guards in any file that was written by the pipeline - for the pipeline
* Barcodes should be unguarded when written as machine inputs, or for final reporting
* All rows and fields should be stripped of white space immediately
* All file reads should be done with ignore=True to avoid non-ASCII characters
* All file reads should protect against empty rows
* Command-line interfaces and web interfaces expect unguarded barcodes but must be explicit for which type of barcode is needed
* The interface codes are responsible for providing guarded barcodes to all internal functions
* Internal functions only accept guarded barcodes

All uploaded files go through upload(exp, files, purpose) and then process_upload(exp, file, purpose) once
the file has been cleared as a safe transaction. After this the various parsers are called. Finally
we rename any transacted plates that are accepted by the user.

Important note: files are state and are immutable in the eyes of the pipeline. Plates are mutable and are
never held in Experiment.reproducible_steps at the top level, and they are never saved as files by the pipeline.

Behaves like a C++ "friend" of the Experiment class - very tightly coupled.
"""

### Universal upload interface

def upload(exp, streams, purpose, overwrite=False):
    """
    All file stream uploads go through this function. 
    - Each file is uploaded with a pending_ prefix and
    - saved to /exp/uploads/
    - A new exp.uploaded_files record is created and
    - it is added to exp.pending_uploads if required,
    - or is sent to process_upload(exp, file, purpose)
    Return True if all uploaded streams are correctly saved and parsed, else False
    """
    if exp.locked and purpose not in {'primer_assay_map','reference_sequences'}:
        exp.log('Error: Cannot add manifest while lock is active')
        return False
    success = True
    for file_stream in streams:
        fp = exp.get_exp_fn(file_stream.name, subdir='uploads', trans=True)
        if Path(fp).exists():
            try:
                Path(fp).unlink()  # delete the existing pending file
            except Exception as exc:
                exp.log(f'Error: upload to {fp} blocked by existing file, {exc}')
                success = False
                continue

        exp.log(f"Info: uploading {file_stream.name} to {fp}")
        with open(fp, 'wb') as outf:
            try:
                copyfileobj(BytesIO(file_stream.getvalue()), outf)
            except Exception as exc:
                exp.log(f'Could not copy stream {file_stream.name} to file {fp}, {exc}')
                success = False
                continue

        # this is essential for retaining the purpose of this file
        exp.add_file_record(fp, purpose=purpose)
        #print(f'upload: {fp=}')
        final_fp = transaction.untransact(fp)
        if Path(final_fp).exists():
            if overwrite:
                try:
                    Path(final_fp).unlink()
                except Exception as exc:
                    exp.log(f'Error: could not overwrite original file {final_fp}')
                    exp.pending_uploads.add(fp)
                    success = False
                    continue
                exp.log(f'Info: successfully deleted original file {final_fp}')
            else:
                # the user will decide later whether to overwrite what they have
                exp.log(f'Info: adding {fp} to pending upload queue')
                exp.pending_uploads.add(fp)
                continue
        else:
            try:
                Path(fp).rename(final_fp)
            except Exception as exc:
                exp.log(f'Error: failed to rename {fp} to {final_fp}, exc')
                success = False
                continue
            exp.log(f'Success: {final_fp} uploaded. Parsing...')
            exp.mod_file_record(fp, new_name=final_fp) # update file name in records
            s = process_upload(exp, final_fp, purpose)
            if s is False:
                success = False
    print(f'Upload succeeded? {success=}', flush=True)
    return success


def accept_pending_upload(exp, fn):
    """
    Take a pending upload file and rename it, overwriting the original.
    An exp.file_upload record must already exist with the file's purpose recorded
    """
    if exp.locked and purpose not in {'primer_assay_map','reference_sequences'}:
        exp.log('Error: Cannot add manifest while lock is active')
        return False

    if fn not in exp.uploaded_files:
        print('1')
        return False
    if 'purpose' not in exp.uploaded_files[fn]:
        print('2')
        return False
    purpose = exp.uploaded_files[fn]['purpose']
    final_fn = transaction.untransact(fn)
    if Path(final_fn).exists():
        try:
            Path(final_fn).unlink()  # delete the existing pending file
        except Exception as exc:
            exp.log(f'Error: upload to {final_fn} blocked by existing file, {exc}')
            return False
        
    if final_fn in exp.uploaded_files:
        success = exp.del_file_record(final_fn)
        if success:
            exp.log(f'Success: deleted original file record for {final_fn}')
        else:
            exp.log(f'Failure: to remove original file record for {final_fn}')

    try:
        exp.pending_uploads.remove(fn)
    except Exception as exc:
        exp.log(f'Failed: removing {fn} {exc}')
        return False

    exp.log(f'Info: replacing existing file {final_fn}, meant for {purpose}')
    try:
        fn_renamed = Path(fn).rename(final_fn)
    except Exception as exc:
        exp.log(f'Error: failed to rename {fn} to {final_fn}, {exc}')
        return False
    exp.mod_file_record(fn, new_name=final_fn) # update file name in records
    exp.log(f'Success: {final_fn} uploaded. Parsing...')
    
    success = process_upload(exp, final_fn, purpose, overwrite=True)
    return success


def process_upload(exp, filepath, purpose, overwrite=False):
    """
    Called by upload() or accept_pending_upload()
    - Clears the file of its pending status (rename)
    - updates the exp.uploaded_files entry
    - invokes the appropriate parser for the provided purpose and file extension
    purposes: ['amplicon','DNA','pcr','rodentity_sample','custom_sample','primer_layout','primer_volume',
            'index_layout','index_volume','primer_assay_map','reference_sequences','taq_water']
    """
    pids = None  # list of plateIDs associated with a file
    if purpose == 'amplicon':
        success, pids = parse_amplicon_manifest(exp, filepath, overwrite=overwrite)
    elif purpose == 'DNA' or purpose == 'dna':
        success, pids = load_dna_plate(exp, filepath, overwrite=True)
    elif purpose == 'rodentity_sample':  # becomes purpose=sample source=rodentity
        success, pids = parse_rodentity_json(exp, filepath, overwrite=overwrite)
    elif purpose == 'custom_sample':  # becomes purpose=sample source=custom
        success, pids = parse_custom_manifest(exp, filepath, overwrite=overwrite)
    elif purpose == 'pcr':
        exp.log('Critical: no parser implemented for PCR plate records')
        # parse_pcr_record(exp, filepath)
    elif purpose == 'primer_layout':
        success, pids = parse_primer_layout(exp, filepath)
    elif purpose == 'primer_volume':
        success, pids = parse_primer_volume(exp, filepath)
    elif purpose == 'index_layout':
        success, pids = parse_index_layout(exp, filepath)
    elif purpose == 'index_volume':
        success, pids = parse_index_volume(exp, filepath)
    elif purpose == 'assay_primer_map':
        success = parse_assay_primer_map(exp, filepath)
    elif purpose == 'reference_sequences':
        success = parse_reference_sequences(exp, filepath)
    elif purpose == 'taq_water':
        exp.log('Critical: no parser implemented for taq_water plates')
    else:
        exp.log(f'Critical: no parser implemented for {purpose=}')

    if success:
        exp.log(f'Success: parsed {filepath} for {purpose}')
        exp.mod_file_record(filepath, extra_PIDs=pids)
    else:
        exp.log(f'Failure: could not parse {filepath} for {purpose}')
        exp.del_file_record(filepath)
    return success


### parsers that build (mostly) json-like dictionary structures for Experiment

# parse_rodentity_json
# parse_custom_manifest
# parse_primer_assay_mapping
# parse_reference_sequences
# parse_taqwater - not necessary, generated
# parse_targets - use parse_reference_sequences
# parse_primer_layout
# parse_primer_volume
# parse_index_layout
# parse_index_volume
# parse_amplicons

@functools.lru_cache
def load_json(fp):
    """ JSON data loader that is cached for fast reading """
    with open(fp, 'rt', errors='ignore') as src:
        info = json.load(src) 

def parse_rodentity_json(exp, fp, overwrite=False):
    """
    Read a JSON file containing one or more Rodentity plate definitions 
    Write a standardised datastructure into exp.plate_location_sample
    Update entry in exp.uploaded_files to include plate ids
    By default, overwriting plate IDs is not allowed for samples - users must delete the plate ID first
    """
    exp.log(f'Begin: parsing {fp}')
    well_records = {}  # [gpid][pos]
    with open(fp, 'rt', errors='ignore') as src:
        info = json.load(src)  

    for record in info['wells']:
        pid = record['plate']
        gpid = util.guard_pbc(pid, silent=True)
        if gpid not in well_records:
            well_records[gpid] = {'purpose':'sample','source':'rodentity', 'wells':set(), 'plate_type':'96'}
            
        pos = util.unpadwell(record['wellLocation'])
        if pos in well_records[gpid] and well_records[gpid][pos]['wells'] != {}:
            # duplicate entry for a given plate and well position
            exp.log(f'Warning: skipping duplicate sample entry {pid} {pos} in {fp}')
            continue
        well_records[gpid]['wells'].add(pos)
        well_records[gpid][pos] = {}
        well_records[gpid][pos]['barcode'] = util.guard_rbc(record['mouse']['barcode'], silent=True)
        well_records[gpid][pos]['sampleBarcode'] = well_records[gpid][pos]['barcode']  # as the line above
        well_records[gpid][pos]['samplePlate'] = gpid
        well_records[gpid][pos]['sampleWell'] = pos
        well_records[gpid][pos]['strain'] = record['mouse']['mouselineName']
        well_records[gpid][pos]['sex'] = record['mouse']['sex']
        well_records[gpid][pos]['mouse'] = deepcopy(record['mouse'])
            
        # build all the various assay/allele structures required
        well_records[gpid][pos]['ngs_assay_records'] = {}
        well_records[gpid][pos]['other_assay_records'] = {} 
        ngs_assays = []
        other_assays = []
                    
        if 'alleles' in record['mouse']:
            for allele in record['mouse']['alleles']:  # allele is a dict from a list
                for assay in allele['assays']:  # assay is a dict from a list
                    assay_name = str(assay['name'])
                    assay_method = str(assay['method'])
                    assay_key = str(assay['assay_key'])
                    allele_key = str(allele['alleleKey'])
                    allele_symbol = '"'+str(allele['symbol'])+'"'
                    if assay_method == 'NGS' or assay_name.startswith('NGS'):
                        ngs_assays.append(assay_name)
                        ngs_or_other = 'ngs_assay_records'
                    else:
                        other_assays.append(assay_name)
                        ngs_or_other = 'other_assay_records'
                    if assay_name not in well_records[gpid][pos][ngs_or_other]:
                        well_records[gpid][pos][ngs_or_other][assay_name] = {}
                    well_records[gpid][pos][ngs_or_other][assay_name]['alleleKey'] = allele_key
                    well_records[gpid][pos][ngs_or_other][assay_name]['alleleSymbol'] = allele_symbol
                    well_records[gpid][pos][ngs_or_other][assay_name]['assayKey'] = assay_key
                    well_records[gpid][pos][ngs_or_other][assay_name]['assayName'] = assay_name # redundant
                    well_records[gpid][pos][ngs_or_other][assay_name]['assayMethod'] = assay_method # not entirely redundant
                            
        well_records[gpid][pos]['ngs_assays'] = ngs_assays.copy()
        well_records[gpid][pos]['other_assays'] = other_assays.copy()
    
    for gpid in well_records:
        if gpid in exp.plate_location_sample:
            if not overwrite:
                exp.log(f'Error: cannot overwrite existing plate entry {gpid}')
                continue
            elif exp.plate_location_sample[gpid]['purpose'] != 'sample':
                exp.log(f'Error: plate {gpid} already exists with purpose {exp.plate_location_sample[gpid]["purpose"]}')
                continue
        exp.plate_location_sample[gpid] = well_records[gpid].copy()

    exp.log(f'End: parsing {fp} with plates: {list(well_records.keys())}')
    return True, list(well_records.keys())
            

def parse_custom_manifest(exp, fp, overwrite=False):
    """
    Parse a custom manifest files and store them in self.unassigned_plates['custom'] =\
            {plateBarcode={Assay=[],...}}
    Also add an entry for each file into self.uploaded_files
    Example headers:
    sampleNum	plateBarcode	well	sampleBarcode	assay	assay	clientName	sampleName 	alleleSymbol
    Required columns: [plateBarcode, well, sampleBarcode, assay*, clientName] *Duplicates allowed
    Optional columns: [sampleNo, sampleName, alleleSymbol] and anything else
    By default existing plate entries cannot be overwritten
    """
    exp.log(f'Begin: parsing {fp}')
    if fp.lower().endswith('xlsx'):
        with open(fp, 'rb') as bytes:
            workbook = openpyxl.load_workbook(bytes)
            sheet = workbook.active
            file_table = [','.join(map(str,cells)) for cells in sheet.iter_rows(values_only=True)]
    else:
        with open(fp, 'rt', errors='ignore') as src:
            file_table = [line.strip() for line in src]

    well_records = {}  # [plate] = {well:{record}}
    plate_barcode_col = None
    assay_cols = []
    required_cols = ['platebarcode', 'well', 'samplebarcode', 'assay', 'clientname'] 
    for i, line in enumerate(file_table):
        cols = [c.strip() if c is not None else '' for c in line.split(',')]  # lower case column names
        #print(i, cols, file=sys.stderr)
        if i==0:  # process header
            cols_lower = [c.lower() for c in cols]
            header_dict = {k:c for k,c in enumerate(cols)}
            #print(header_dict, file=sys.stderr)
            matching_cols = [col_name in cols_lower for col_name in required_cols]
            if not all(matching_cols):
                exp.log(f'Error: manifest {fp} requires at least columns {required_cols}')
                return False, []
            
            for k in header_dict:
                if header_dict[k].lower() == 'platebarcode':
                    plate_barcode_col = k
                elif header_dict[k].lower() == 'assay':
                    assay_cols.append(k)
            continue
        if len(cols) < 5:  # skip empty rows
            continue
        pid = cols[plate_barcode_col]
        gpid = util.guard_pbc(pid, silent=True)
        
        if gpid not in well_records:
            # create plate_location_sample entries here
            well_records[gpid] = {'purpose':'sample','source':'custom', 'wells':set(), 'plate_type':'96'}  
        assays = [] 
        # do a first pass to set up recording a sample in a well
        for k,c in enumerate(cols):
            if c.lower() == 'none':
                continue
            if k in assay_cols:
                if c != '':  # ignore empty assay entries
                    assays.append(c)
            if header_dict[k].lower() == 'well':
                well = util.unpadwell(c.upper())
                #print(f'well: {well=} {k=} {c=}')
                if well in well_records[gpid]['wells'] and well_records[gpid][well] != {}:
                    exp.log(f"Warning: duplicate {c} in {pid}. Skipping row {i+2} {cols=}")
                    break
                well_records[gpid]['wells'].add(well)
                well_records[gpid][well] = {'sampleWell':well}
        # now do a second pass to collect everything together
        for k,c in enumerate(cols):
            if c.lower() == 'none':
                c = ''
            if header_dict[k].lower() == 'well':
                continue
            if header_dict[k].lower() == 'samplebarcode':
                #if c.startswith('C'):
                #    sid = util.guard_cbc(c, silent=True)
                #elif c.startswith('M'):
                #    sid = util.guard_rbc(c, silent=True)
                #else:
                sid = util.guard_cbc(c, silent=True) # fall back to custom?
                well_records[gpid][well]['barcode'] = sid
                well_records[gpid][well]['sampleBarcode'] = sid
            elif header_dict[k].lower() == 'platebarcode':
                well_records[gpid][well]['platebarcode'] = gpid
                well_records[gpid][well]['samplePlate'] = gpid
            elif header_dict[k].lower() == 'sampleno':
                well_records[gpid][well]['sampleNumber'] = str(c)
            elif header_dict[k].lower() == 'assay':
                continue
            elif header_dict[k].lower() == 'samplename':
                well_records[gpid][well]['sampleName'] = str(c)
            elif header_dict[k].lower() == 'clientname':
                well_records[gpid][well]['clientName'] = str(c)
            else:
                try:
                    well_records[gpid][well][header_dict[k]] = str(c)
                except:
                    print(f'Failed to record well_record {gpid=} {well=} {header_dict[k]=} {c=}', file=sys.stderr)
        well_records[gpid][well]['ngs_assays'] = assays

    for gpid in well_records:
        if gpid in exp.plate_location_sample:
            if not overwrite:
                exp.log(f'Error: cannot overwrite existing plate entry {gpid}')
                continue
            elif exp.plate_location_sample[gpid]['purpose'] != 'sample':
                exp.log(f'Error: plate {gpid} already exists with purpose {exp.plate_location_sample[gpid]["purpose"]}')
                continue
        exp.plate_location_sample[gpid] = well_records[gpid].copy()

    exp.log(f'End: parsing {fp} with plates: {list(well_records.keys())}')
    return True, list(well_records.keys())


def parse_amplicon_manifest(exp, fp, overwrite=False):
    """
    Amplicon plate layouts (containing pre-amplified sequences) are parsed here from manifest files.
    Only column layout (comma separate) is supported. col1: plate barcode; col2: well position; col3: sample barcode; col4 (optional): volume.
    The first row must be a header row: plate, well, sample, (volume in uL, if provided).
    plateBarcode, well, sampleBarcode, amplicon* (*multiple columns with this name are allowed)
    Adding this will result in changing the Miseq and Stage3 output files, so needs to be wrapped as a transaction
    Needs to check whether there are sufficient indexes

    By default existing plate ID records cannot be overwritten and should be first deleted by users
    """
    exp.log(f'Begin: parsing {fp}')
    if fp.lower().endswith('xlsx'):
        with open(fp, 'rb') as bytes:
            workbook = openpyxl.load_workbook(bytes)
            sheet = workbook.active
            file_table = [','.join(map(str,cells)) for cells in sheet.iter_rows(values_only=True)]
    else:
        with open(fp, 'rt', errors='ignore') as src:
            file_table = [line.strip() for line in src]

    well_records = {}     
    plate_barcode_col = None
    well_col = None
    amplicon_cols = []  # combine these
    required_cols = ['platebarcode', 'well', 'samplebarcode']
    for i, line in enumerate(file_table): # csv_reader is super limiting
        cols = [c.strip() if c is not None else '' for c in line.split(',')]  # lower case column names
        if i==0:  # process header
            cols_lower = [c.lower() for c in cols]
            header_dict = {k:c for k,c in enumerate(cols_lower)}
            matching_cols = [col_name in cols_lower for col_name in required_cols] 
            if not all(matching_cols):
                exp.log(f'Error: amplicon manifest {fp} requires at least {required_cols} columns')
                return False, []
            
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
        pid = cols[plate_barcode_col]
        gpid = util.guard_pbc(pid, silent=True)
        if gpid not in well_records:
            well_records[gpid] = {'purpose':'amplicon','source':'custom', 'wells':set(), 
                    'plate_type':'384PP_AQ_BP'}

        # set up well
        pos = util.unpadwell(cols[well_col].upper())
        if pos in well_records[gpid]['wells'] and well_records[gpid][pos] != {}:
            exp.log(f"Warning: duplicate {c} in {pid}. Skipping row {i+2} {cols=}")
            continue
        well_records[gpid]['wells'].add(pos)
        well_records[gpid][pos] = {'sampleWell':pos}
        amplicons = []
        # now collect everything together
        for k,c in enumerate(cols):
            if c.lower() == 'none':
                c = ''
            if header_dict[k] == 'well':
                continue  # we've already got this!
            if header_dict[k] == 'samplebarcode':
                # amplicon guards
                sid = well_records[gpid][pos][header_dict[k]] = util.guard_abc(c, silent=True)
                well_records[gpid][pos]['barcode'] = sid
                well_records[gpid][pos]['sampleBarcode'] = sid
            elif header_dict[k] == 'platebarcode':
                # don't really need it but whatever - record the guarded PID
                well_records[gpid][pos]['platebarcode'] = gpid
                well_records[gpid][pos]['samplePlate'] = gpid
            elif header_dict[k] == 'sampleno':
                well_records[gpid][pos]['sampleNumber'] = str(c)
            elif header_dict[k] == 'amplicon':
                if c != '':  # ignore empty amplicon entries
                    amplicons.append('"'+c+'"')
            elif header_dict[k] == 'volume':
                well_records[gpid][pos]['volume'] = float(c)*1000  # save as nL
            else:
                try:
                    well_records[gpid][pos][header_dict[k]] = str(c)
                except:
                    print(type(header_dict[k]), header_dict[k], str(c), file=sys.stderr)
        well_records[gpid][pos]['amplicons'] = amplicons

    for gpid in well_records:
        if gpid in exp.plate_location_sample:
            if not overwrite:
                exp.log(f'Error: cannot overwrite existing plate entry {gpid}')
                continue
            elif exp.plate_location_sample[gpid]['purpose'] != 'amplicon':
                exp.log(f'Error: plate {gpid} already exists with purpose {exp.plate_location_sample[gpid]["purpose"]}')
                continue
        exp.plate_location_sample[gpid] = well_records[gpid].copy()

    exp.log(f'End: parsing {fp} with plates: {list(well_records.keys())}')
    return True, list(well_records.keys())


def parse_primer_layout(exp, fp, user_gpid=None):
    """
    add primer plate definition with well and name columns
    """
    if user_gpid:
        PID = util.unguard_pbc(user_gpid, silent=True)
        gPID = util.guard_pbc(user_gpid, silent=True)
    else:
        PID = Path(fp).name.split('_')[0]  # assumes first field is PID
        gPID = util.guard_pbc(PID, silent=True)

    # parse the file first to make sure it works
    well_records = {gPID:{}}
    with open(fp, 'rt', errors='ignore') as data:
        for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
            if i == 0:
                continue  # header
            if row == '' or row[0] == '' or row[1] == '':
                continue
            try:
                well = util.unpadwell(row[0])
            except UnboundLocalError as exc:
                exp.log(f'Error: You likely tried to upload volumes instead of layout! Using file {fp}')
                return False, []
            if well in well_records[gPID]:
                exp.log(f'Warning: skipping duplicate well entry {well} in {fp}')
                continue 
            well_records[gPID][well] = {'primer':row[1]}
    
    # check for a legitimate plate match
    if gPID in exp.plate_location_sample:
        if exp.plate_location_sample[gPID]['purpose'] == 'primer':
            exp.log(f"Info: Primer plate {PID} already exists, adding data")
        else:
            exp.log(f"Error: Primer plate PID: {PID} matches "+\
                    f"existing plate entry of different purpose "+\
                    f"{exp.plate_location_sample[gPID]['purpose']}")
            return False, []
    else:  # create new plate entry and set purpose
        exp.plate_location_sample[gPID] = {'purpose': 'primer', 'source':'user', 'wells':set(), 
                'plate_type':util.PLATE_TYPES['Echo384']}
        exp.log(f"Info: Creating new primer plate record for {PID}")

    # copy across data from the temporary plate records                           
    for well in well_records[gPID]:
        if well not in exp.plate_location_sample[gPID]:
            exp.plate_location_sample[gPID]['wells'].add(well)
            exp.plate_location_sample[gPID][well] = {}
        exp.plate_location_sample[gPID][well]['primer'] = well_records[gPID][well]['primer']
                    
    exp.log(f"Success: added primer layout from {fp}")
    return True, [gPID]


def parse_primer_volume(exp, fp, user_gpid=None):
    """
    add primer plate volumes from either plate layout or two-column layout
    """
    if user_gpid:
        PID = util.unguard_pbc(user_gpid, silent=True)
        gPID = util.guard_pbc(user_gpid, silent=True)
    else:
        PID = Path(fp).name.split('_')[0]  # assumes first field is PID
        gPID = util.guard_pbc(PID, silent=True)

    # parse the file first to make sure it works
    well_records = {gPID:{}}
    with open(fp, 'rt', errors='ignore') as data:
        for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
            if i == 0:
                if row[0].lower().startswith('date'):
                    # plate style layout
                    layout = 'plate'
                else:
                    layout = 'columns'
                continue  # header
                
            if layout == 'plate':
                if not row or row[0] == '' or row[0][0] not in ascii_uppercase[:16]:  # A-P
                    continue
                for j,col in enumerate(row):
                    if j==0:
                        continue  # row name cell
                    well = util.unpadwell(row[0] + str(j))
                    if well in well_records[gPID]:
                        exp.log(f'Warning: skipping duplicate well entry {well} in {fp}')
                        continue
                    well_records[gPID][well] = {'volume':float(col)*1000}
            else:  # columns
                well = util.unpadwell(row[0])
                if well in well_records[gPID]:
                    exp.log(f'Warning: skipping duplicate well entry {well} in {fp}')
                    continue
                well_records[gPID][well] = {'volume':float(col)*1000}

    if gPID in exp.plate_location_sample:
        if exp.plate_location_sample[gPID]['purpose'] == 'primer':
            exp.log(f"Info: Primer plate {PID} already exists, adding data")
        else:
            exp.log(f"Error: Primer plate PID: {PID} matches "+\
                    f"existing plate entry of different purpose "+\
                    f"{exp.plate_location_sample[gPID]['purpose']}")
            return False, []
    else:  # create new plate entry and set purpose
        exp.plate_location_sample[gPID] = {'purpose': 'primer', 'source':'user', 'wells':set(), 
                'plate_type':util.PLATE_TYPES['Echo384']}
        exp.log(f"Info: Creating new primer plate record for {PID}")

    for well in well_records[gPID]:
        if well not in exp.plate_location_sample[gPID]:
            exp.plate_location_sample[gPID]['wells'].add(well)
            exp.plate_location_sample[gPID][well] = {}
        exp.plate_location_sample[gPID][well]['volume'] = well_records[gPID][well]['volume']
                    
    exp.log(f"Success: added primer volumes from {fp}")
    return True, [gPID]

 
def parse_index_layout(exp, fp, user_gpid=None):
    """ 
    Add index plates with well and index columns 
    Inputs: experiment, fp - filepath usually containing the PID, user_gpid is an override for the filepath
    Non-sample plates can be loaded straight into experiment.plate_location_sample - save exp when done
    Returns True on success
    """
    if user_gpid:
        PID = util.unguard_pbc(user_gpid, silent=True)
        gPID = util.guard_pbc(user_gpid, silent=True)
    else:
        PID = Path(fp).name.split('_')[0]  # assumes first field is PID
        gPID = util.guard_pbc(PID, silent=True)   
    
    # parse the file first to make sure it works
    well_records = {gPID:{}}
    with open(fp, 'rt', errors='ignore') as data:
        for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_MINIMAL)):
            if i == 0 or row == '':
                continue  # header or blank
            try:
                well = util.unpadwell(row[0])
            except UnboundLocalError as exc:
                exp.log(f'Did you upload a volume file by mistake? {exc}')
                return False, []
            if well in well_records[gPID]:
                exp.log(f'Warning: skipping duplicate well entry {well} in {fp}')
                continue
            well_records[gPID][well] = {}
            idt_name = row[1]
            index = row[2]
            bc_name = row[3]
            oligo = row[4]
            well_records[gPID][well]['idt_name'] = idt_name
            well_records[gPID][well]['index'] = index
            well_records[gPID][well]['bc_name'] = bc_name
            well_records[gPID][well]['oligo'] = oligo

    if gPID in exp.plate_location_sample:
        if exp.plate_location_sample[gPID]['purpose'] == 'index':
            exp.log(f"Info: Index plate {PID} already exists, adding data")
        else:
            exp.log(f"Error: Index plate PID: {PID} matches "+\
                    f"existing plate entry of different purpose "+\
                    f"{exp.plate_location_sample[gPID]['purpose']}")
            return False, []
    else:  # create new plate entry and set purpose
        exp.plate_location_sample[gPID] = {'purpose': 'index', 'source':'user', 'wells':set(), 
                'plate_type':util.PLATE_TYPES['Echo384']}
        exp.log(f"Info: Creating new index plate record for {PID}")

    for well in well_records[gPID]:
        if well not in exp.plate_location_sample[gPID]:
            exp.plate_location_sample[gPID]['wells'].add(well)
            exp.plate_location_sample[gPID][well] = {}
        exp.plate_location_sample[gPID][well]['idt_name'] = well_records[gPID][well]['idt_name']
        exp.plate_location_sample[gPID][well]['index'] = well_records[gPID][well]['index']
        exp.plate_location_sample[gPID][well]['bc_name'] = well_records[gPID][well]['bc_name']
        exp.plate_location_sample[gPID][well]['oligo'] = well_records[gPID][well]['oligo']

    exp.log(f"Success: added index plate layouts from {fp}")        
    return True, [gPID]
    

def parse_index_volume(exp, fp, user_gpid=None):
    """ add volume information to a single barcode plate, with a matched, possibly-existing plate name.
    To make life fun there are two possible formats:
    Format 1: plate layout (generated through Echo test software)
    Format 2: column layout (generated through Echo main interface)
        
    Inputs: experiment, fp - filepath usually containing the PID, user_gpid is an override for the filepath
    Non-sample plates can be loaded straight into experiment.plate_location_sample - save exp when done
    Returns True on success
    """
    if user_gpid:
        PID = util.unguard_pbc(user_gpid, silent=True)
        gPID = user_gpid
    else:
        PID = Path(fp).name.split('_')[0]  # assumes first field is PID
        gPID = util.guard_pbc(PID, silent=True)   
    
    # parse the file first to make sure it works
    well_records = {gPID:{}}

    plate_format = False
    hdr_seen = False
    with open(fp, 'rt', errors='ignore') as data:
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
                        if well in well_records[gPID]:
                            exp.log(f'Warning: duplicate well entry {well} seen in plate {PID}, overwriting')
                        well_records[gPID][well] = {'volume': float(v)*1000}  
            else:                
                #format 2 - one row per well - well ID in r[3], volume in r[5] or r[6]
                #if i == 0 and line.startswith('Run ID'):
                if row.startswith('[DETAILS]'):
                    continue
                if ''.join(row).strip() == '':
                    continue
                if not hdr_seen:
                    hdr_seen = True
                    continue
                
                # finally, the data    
                well = util.unpadwell(row[3])
                if well in well_records[gPID]:
                    exp.log(f'Warning: duplicate well entry {well} seen in plate {PID}, overwriting')
                well_records[gPID][well] = {'volume': int(float(row[5]))*1000}

    if gPID not in exp.plate_location_sample:
        exp.log(f"Info: Adding index plate {PID}")
        exp.plate_location_sample[gPID] = {'purpose':'index', 'source':'user', 'wells':set(), 
                'plate_type':util.PLATE_TYPES['Echo384']}
    else:
        if exp.plate_location_sample[gPID]['purpose'] != 'index':
            exp.log(f"Error: {PID} plate purpose is "+\
                    f"{exp.plate_location_sample[gPID]['purpose']}, expected 'index'")
            return False, []
        exp.log(f"Info: {PID} exists, appending index volumes")

    for well in well_records[gPID]:
        if well not in exp.plate_location_sample[gPID]:
            exp.plate_location_sample[gPID]['wells'].add(well)
            exp.plate_location_sample[gPID][well] = {}
        exp.plate_location_sample[gPID][well]['volume'] = well_records[gPID][well]['volume']
    
    return True, [gPID]


def parse_assay_primer_map(exp, fp):
    """ mapping of assay name, assay family name, and primers
    The assaylist file has two fixed columns: primer and assay
        - every subsequent column is for an extra primer in that assay family
    In Experiment we maintain:
        self.assayfam_primer = {}  # mapping of assay families to list of primers
        self.assay_assayfam = {}  # mapping of assay to assay family
        self.assayfam_assay = {}  # mapping of assay family to list of assays, possibly redundant
        self.primer_assayfam = {}  # reverse mapping
        
    Case is strictly enforced - we don't bother even checking case mistake, just reporting them
    """
    # read the whole file before adding entries
    local_assayfam_primers = {}
    local_assay_assayfam = {}
    local_assayfam_assays = {}
    local_primer_assayfam = {} 
    with open(fp, 'rt', errors='ignore') as data:
        for i, row in enumerate(csv.reader(data, delimiter=',', quoting=csv.QUOTE_ALL)):
            if i == 0:
                continue  # header
            if row == '' or row[0] == '' or row[1] == '':
                continue  # empty row

            assay = row[0].strip()
            if not assay:
                continue  # just skip it
            assay_fam = row[1].strip()
            if not assay_fam:
                exp.log(f'Warning: no assay family associated with {assay} in row {i+1}')
                continue
            primers = [p.strip() for p in row[2:] if p.strip() != '']
            if not primers:
                exp.log(f'Warning: no primers associated with assay {assay} in row {i+1}')
                continue

            # assay to assay family
            if assay in local_assay_assayfam:
                exp.log(f'Warning: {assay} has duplicate entry {local_assay_assayfam[assay]} seen in row {i+1}')
            local_assay_assayfam[assay] = assay_fam

            # assay family to assay
            if assay_fam not in local_assayfam_assays:
                local_assayfam_assays[assay_fam] = []
            local_assayfam_assays[assay_fam].append(assay)

            # assay family to primers
            if assay_fam not in local_assayfam_primers:
                local_assayfam_primers[assay_fam] = primers
            else:
                for primer in primers:
                    if primer not in local_assayfam_primers[assay_fam]:
                        local_assayfam_primers[assay_fam].append(primer)

            # reverse mapping, primers to assayfams
            for primer in primers:
                if primer in local_primer_assayfam:
                    if local_primer_assayfam[primer] != assay_fam:
                        msg = f'Warning: primer {primer} has multiple assayfams: '+\
                                f'{local_primer_assayfam[primer]} and {assay_fam}'
                        exp.log(msg)
                        print(msg)
                        continue
                local_primer_assayfam[primer] = assay_fam

    for a in local_assay_assayfam:
        exp.assay_assayfam[a] = local_assay_assayfam[a]
    for af in local_assayfam_assays:
        exp.assayfam_assays[af] = local_assayfam_assays[af]
    for af in local_assayfam_primers:
        exp.assayfam_primers[af] = local_assayfam_primers[af]
    for primer in local_primer_assayfam:
        exp.primer_assayfam[primer] = local_primer_assayfam[primer]
    
    exp.log(f'Success: added assay/assayfam/primer mapping from {fp}')
    return True


def parse_reference_sequences(exp, fp):
    """
    read in reference (target) IDs and sequences.
    May be added when pipeline is locked.
    Client allows only one reference file to be loaded at a time
    """
    exp.log(f'Begin: parsing reference sequences {fp}')
    ref_seq = {}
    with open(fp, 'rt', errors='ignore') as data:
        for ref_seq_pair in SimpleFastaParser(data):
            ref_seq[ref_seq_pair[0]] = ref_seq_pair[1]

    if fp in exp.reference_sequences:
        msg = f'Warning: {fp} already uploaded. Overwriting records'
        exp.log(msg)
    exp.reference_sequences = {}  # Only one reference file allowed by client
    exp.reference_sequences[fp] = {}
    for ref in ref_seq:
        exp.reference_sequences[fp][ref] = ref_seq[ref]

    exp.log(f'Success: {len(ref_seq)} reference sequences imported from {fp}')
    return True
  

def load_dna_plate(exp, filepath, overwrite=True):
    """
    Filepath should point to an Echo_384_COC file/Nimbus output/Echo input
    """
    # "Echo_384_COC_0001_" + ug_dnaBC + "_0.csv"
    with open(filepath, 'rt') as f:
        #RecordId	TRackBC	TLabwareId	TPositionId	SRackBC	SLabwareId	SPositionId
        #1	p2021120604p	Echo_384_COC_0001	A1	p2111267p	ABg_96_PCR_NoSkirt_0001	A1
        dbc = util.guard_pbc(filepath.split('_')[-2], silent=True)
        exp.plate_location_sample[dbc] = {'purpose':'dna','wells':set(),'source':'','plate_type':'384PP_AQ_BP'}
        source_plate_set = set()
        for i, line in enumerate(f):
            if i == 0:  # header
                # get named columns
                cols = [c.strip() for c in line.split(',')]
                source_plate_bc_col = [j for j,c in enumerate(cols) if c=='"SRackBC"'][0]
                source_well_col = [j for j,c in enumerate(cols) if c=='"SPositionId"'][0]
                source_bc_col = [j for j,c in enumerate(cols) if c=='"SPositionBC"'][0]
                dest_plate_bc_col = [j for j,c in enumerate(cols) if c=='"TRackBC"'][0]
                dest_well_col = [j for j,c in enumerate(cols) if c=='"TPositionId"'][0]
                continue
            if line.strip() == '':
                continue  # skip empty lines
            cols = [c.strip().strip('\"') for c in line.split(',')]
            source_bc = cols[source_bc_col]
            if source_bc == '0M':
                continue  # empty well
            source_pos = util.unpadwell(cols[source_well_col])
            source_plate = util.guard_pbc(cols[source_plate_bc_col], silent=True)
            source_plate_set.add(source_plate)
            dest_plate = util.guard_pbc(cols[dest_plate_bc_col], silent=True)
            dest_pos = util.unpadwell(cols[dest_well_col])
            if dest_plate != dbc:
                exp.log(f"Error: {dbc} doesn't match {dest_plate} as declared in Echo_384_COC file: {filepath}")
                return False, [dbc]         
            try:
                exp.plate_location_sample[dbc][dest_pos] = exp.plate_location_sample[source_plate][source_pos]
            except:
                exp.log(f"Critical: cannot locate {dbc=} {dest_pos=} {source_plate=} {source_pos=}")
                return False, [dbc]
            exp.plate_location_sample[dbc]['wells'].add(dest_pos)
            exp.plate_location_sample[dbc]['samplePlate'] = source_plate
            exp.plate_location_sample[dbc]['sampleWell'] = source_pos
        exp.plate_location_sample[dbc]['source'] = ','.join(source_plate_set)
    return True, [dbc]


def myopen(fn):
    """ Bob's function to handle gzip transparently """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt", errors='ignore')
    return open(fn, errors="ignore")


def read_html_template(html_fn):
    """ 
    Read an html file from the library folder, ready for additional fields 
    Used by cgi-dashboard.py and makehtml.py
    template uses {!field!} instead of python format fields - easier to edit.
    Files should use the .tpl (template) extension as they are not strictly html
    """
    tfn = os.path.join('..', 'library', html_fn)
    with open(tfn, errors='ignore') as src:
        htmlcode = src.read()
    # template uses {!field!} instead of python format fields - easier to edit.
    htmlfmt = htmlcode.replace('{', '{{').replace('}', '}}').replace('{{!', '{').replace('!}}', '}')
    return htmlfmt


def read_csv_or_excel_from_stream(multi_stream):
    """
    Take a stream of one or more table files from CSV or Excel and convert to comma delimited rows
    return a list of names and a list of tables

    DEPRECATED
    """
    manifest_names = []
    manifest_tables = []
    for file_stream in multi_stream:
        manifest_name = file_stream.name
        #print(f'{manifest_name=}', file=sys.stderr)
        if manifest_name.lower().endswith('xlsx'):
            workbook = openpyxl.load_workbook(BytesIO(file_stream.getvalue()))
            sheet = workbook.active
            rows = [','.join(map(str,cells)) for cells in sheet.iter_rows(values_only=True)]
            #print(rows, file=sys.stderr)
        else:
            rows = StringIO(file_stream.getvalue().decode("utf-8"))
        manifest_names.append(manifest_name)
        manifest_tables.append(rows)
    return manifest_names, manifest_tables


if __name__ == '__main__':
    """ library only """
    pass