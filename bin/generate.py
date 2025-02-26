#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#from multiprocessing import popen_spawn_win32
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
    
from stutil import m

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

Important note: files are state and are immutable in the eyes of the pipeline. Plates are mutable and are
never held in Experiment.reproducible_steps at the top level, and they are never saved as files by the pipeline.

Behaves like a C++ "friend" of the Experiment class - very tightly coupled.
"""

### Generators
# nimbus_gen
# mk_picklist
# mk_mytaq_picklist
# generate_echo_PCR1_picklist
# generate_echo_PCR2_picklist
# generate_miseq_samplesheet
# generate_targets

def pcr1_picklists_exist(exp):
    """
    *Stage 3: PCR1*
    Check if the files that are generated for PCR 1 exist
    """
    picklist_files = ['Stage2.csv', f'PCR1_dna-picklist_{exp.name}.csv', 
            f'PCR1_primer-picklist_{exp.name}.csv', f'PCR1_taqwater-picklist_{exp.name}.csv']
    return all(os.path.exists(exp.get_exp_fn(file)) for file in picklist_files)


def pcr2_picklists_exist(exp):
    """
    *Stage 4: PCR2*
    Check if the files that are generated for PCR 2 exist
    """
    picklist_files = ['Stage3.csv', f'PCR2_index-picklist_{exp.name}.csv', 
            f'PCR2_taqwater-picklist_{exp.name}.csv', f'Miseq_{exp.name}.csv']
    
    return all(os.path.exists(exp.get_exp_fn(file)) for file in picklist_files)


@functools.lru_cache
def get_single_pos(resource_dict, name_offset, res_name, rotate=True):
    """
    Generic plate and well chooser for a single resource that manages well capacity 
            with either rotating or serial behaviour
    resource_dict is provided by exp.get_available_primer_wells
    resource_dict[res_name] = [[pid, well, vol, capacity], ...]
    name_offset[res_name] = int is for indexing into the resource_dict[res_name] list
    if rotate, we increment the offset regardless of whether the well we just looked
            in is empty
    """
    if name_offset[res_name] == -1:
        # res_name has run dry
        return None, None
    starting_offset = name_offset[res_name]
    capacity = 0
    offset = name_offset[res_name]
    while capacity == 0:
        pid = resource_dict[res_name][offset][0]
        well = resource_dict[res_name][offset][1]
        capacity = resource_dict[res_name][offset][2]       
        if capacity == 0:  # empty well
            offset += 1
            if offset >= len(resource_dict[res_name]):  # wrap around
                offset = 0
            if offset == starting_offset:  # all wells are dry
                name_offset[res_name] = -1
                return None, None
        else:
            capacity -= 1
            resource_dict[res_name][offset] = [pid, well, capacity]
            break
    if rotate:  # use each well once before coming back
        offset += 1
        if offset >= len(resource_dict[res_name]):
            offset = 0
    name_offset[res_name] = offset
    return pid, well
    

def primer_plating(exp, table, pmr_pids, rotate=True):
    """ 
    table is a Table object
    If rotate is True: use each well of a primer once before moving to the next well
    If rotate is False: use the same well of a primer until it is empty before moving on
    """
    pmr_pos_list = []
    primer_resources = exp.get_available_primer_wells(pmr_pids)
    primer_offsets = {pmr:0 for pmr in primer_resources}
    for row in table:
        pid, well = get_single_pos(primer_resources, primer_offsets, row.primer, rotate=rotate)
        pmr_pos_list.append(tuple([pid, well]))
    return pmr_pos_list


def check_assays_known(exp, sample_pid):
    """
    For a given sample plate ID (96-well ear punch plate or custom plate) find which assays are known
    and which are not. Return Counters for known and unknown.
    """
    known_assays = Counter()
    unknown_assays = Counter()
    
    pid = util.unguard_pbc(sample_pid, silent=True)
    for well in exp.plate_location_sample[pid]['wells']:
        for assay in exp.plate_location_sample[pid][well]['ngs_assays']:
            if assay in exp.assay_assayfam:
                known_assays[assay] += 1
            else:
                unknown_assays[assay] += 1
    return known_assays, unknown_assays


def nimbus_gen(exp, caller_id=None):
    """
    Standardised Nimbus input generator, produces a workfile for the BRF's Nimbus robot.
    
    It needs the name of an output file and sample data for each plate.
    dna_plate_id = barcode or name of target output files, which should be guarded already

    A minimum of 8 well entries must exist in the same column or the Nimbus complains

    As a file generator() it needs to create a transaction (exp.reproducible_steps) entry
    """
    if len(exp.assay_assayfam) == 0:
        m('assay/primer "assay list" must be loaded first', level='error', caller_id=caller_id)
        return False
    transactions = {}
    #try:
    if True:
        
        for dna_BC in exp.dest_sample_plates:
            ug_dnaBC = util.unguard_pbc(dna_BC, silent=True)

            dna_fn = exp.get_exp_fn('Nimbus-'+ug_dnaBC+'.csv', trans=True)
            fnstg = exp.get_exp_fn('Stage1-P'+ug_dnaBC+'.csv', trans=True)

            transactions[dna_fn] = {} # add plates and modifications to this
            transactions[fnstg] = {}
            # rn - row number in Nimbus picklist file - Nimbus needs this
            # wn - count the number of wells needed in the DNA plate
            # plist - list of primer families needed for the samples in the plate
            # unk - unknown assay family names; names that don't map to a primer family name
            row_number, wells_used = 0, 0
            with open(dna_fn, "wt", newline='') as dna_fd, open(fnstg, "wt", newline='') as stage1_fd:
                # Note: the Nimbus is extremely particular about column names.
                nimbus_f  = csv.writer(dna_fd, dialect='unix')
                nimbus_hdr = ['Sample no', 'Plate barcode', 'Well', 'Sample barcode']
                nimbus_f.writerow(nimbus_hdr)
                stage1_f = csv.writer(stage1_fd, dialect='unix')
                stage1_hdr = tuple(['samplePlate','sampleWell','sampleBarcode','strain','sex','alleleSymbol','alleleKey',
                        'assayKey','assays','assayFamilies','clientName','sampleName'])
                stage1_f.writerow(stage1_hdr)
                for pbc in sorted(exp.dest_sample_plates[dna_BC]):
                    shx = exp.plate_location_sample.get(pbc,None)
                    if not shx:
                        m(f'No plate data found for plate {pbc}', level='error', caller_id=caller_id)
                        continue
                    transactions[dna_fn][pbc] = {}
                    transactions[fnstg][pbc] = {}
                    # Get eight wells in the same column at a time
                    #for pos1,pos2,pos3,pos4 in zip(util.nimbus_ordered_96[0::4], util.nimbus_ordered_96[1::4], 
                    #        util.nimbus_ordered_96[2::4], util.nimbus_ordered_96[3::4]): 
                    warnings = set()
                    skipped_assays_pos = {}
                    for p in range(len(util.nimbus_ordered_96)//8):
                        pos_col = [util.nimbus_ordered_96[p*8+offset] for offset in range(8)]
                        # skip empty column
                        if all([posN not in shx for posN in pos_col]):
                            continue
                        #if pos1 not in shx and pos2 not in shx and pos3 not in shx and pos4 not in shx:
                        #    continue  # skip missing blocks of 4 contiguous rows in the same column
                        for pos in pos_col:
                            row_number += 1
                            nim_row_data = {field:'' for field in stage1_hdr}
                            nim_row_data['Sample no'] = row_number
                            if pos in shx:
                                nim_row_data['Sample barcode'] = util.unguard(shx[pos]['barcode'], silent=True)
                            else:
                                nim_row_data['Sample barcode'] = '0'
                            nim_row_data['Well'] = pos
                            nim_row_data['Plate barcode'] = util.unguard_pbc(pbc, silent=True)
                            nimbus_f.writerow([nim_row_data[field] for field in nimbus_hdr])

                            if nim_row_data['Sample barcode'] == '0':
                                continue
                    
                            # stage files are still used, but as interlocks on stages (transactions)
                            stage1_row_data = {field:'' for field in stage1_hdr}
                            if 'sampleNumber' in shx[pos]:
                                stage1_row_data['sampleNumber'] = shx[pos]['sampleNumber']
                            else:
                                stage1_row_data['sampleNumber'] = 0
                            stage1_row_data['samplePlate'] = pbc
                            stage1_row_data['sampleWell'] = pos
                            stage1_row_data['sampleBarcode'] = shx[pos]['barcode']
                            stage1_row_data['strain'] = shx[pos].get('strain','')
                            stage1_row_data['sex'] = shx[pos].get('sex','')
                            alleleSymbols = []
                            alleleKeys = []
                            assayKeys = []
                            assayNames = []
                            assayFamilies = []
                            for assay in shx[pos]['ngs_assays']:
                                if assay not in exp.assay_assayfam:
                                    if assay not in skipped_assays_pos:
                                        skipped_assays_pos[assay] = []
                                    skipped_assays_pos[assay].append(pos)
                                    continue
                                assayNames.append(assay)
                                assayFamilies.append(exp.assay_assayfam[assay])
                                if 'ngs_assay_records' in shx[pos]:
                                    if assay in shx[pos]['ngs_assay_records']:
                                        alleleSymbols.append(shx[pos]['ngs_assay_records'][assay].get('alleleSymbol',''))
                                        alleleKeys.append(shx[pos]['ngs_assay_records'][assay].get('alleleKey',''))
                                        assayKeys.append(shx[pos]['ngs_assay_records'][assay].get('assayKey',''))
                                    
                            stage1_row_data['alleleSymbol'] = ';'.join(alleleSymbols)
                            stage1_row_data['alleleKey'] = ';'.join(alleleKeys)
                            stage1_row_data['assayKey'] = ';'.join(assayKeys)
                            stage1_row_data['assays'] = ';'.join(assayNames)
                            stage1_row_data['assayFamilies'] = ';'.join(assayFamilies)
                            
                            if len(stage1_row_data['assays']) == 0:
                                continue 

                            cn = shx[pos].get('clientName','')
                            if cn:
                                cn = cn.replace(' ','_').replace('"','').replace("'",'').replace(';','').replace('`','')
                            stage1_row_data['clientName'] = cn
                            
                            sn = shx[pos].get('sampleName','')
                            if sn:
                                sn = sn.replace(' ','_').replace('"','').replace("'",'').replace(';','').replace('`','')
                            stage1_row_data['sampleName'] = sn
                            
                            stage1_f.writerow([stage1_row_data[field] for field in stage1_hdr])  # assays are written here
                            wells_used += 1
                            transactions[dna_fn][pbc][pos] = -1000 # 1000 nl of sample is transferred
                            transactions[fnstg][pbc][pos] = -1000
                    ug_dna_BC = util.unguard_pbc(dna_BC, silent=True)
                    if skipped_assays_pos:
                        print(f'{skipped_assays_pos=}', flush=True)
                        assay_line = '; '.join([f'{assay} ({",".join(skipped_assays_pos[assay])})' for assay in skipped_assays_pos])
                        warnings.add(f'skipping assays: {assay_line} in sample plate {ug_dna_BC}')
                    for w in warnings:
                        m(w, level='warning', caller_id=caller_id)
                
    #except Exception as exc:
    #    m(f'Nimbus input file creation {exc}', level='error')
    #    exp.save()
    #    return False
    #print("Transactions in nimbus_gen()", transactions, file=sys.stderr)
    transaction.add_pending_transactions(exp, transactions)
    m(f'Hamilton Nimbus plate definition files have been generated', level='success', caller_id=caller_id)
    return True


class PicklistSrc:
    """ read a survey & contents file (echovolume.py output) """ 

    def __init__(self, fn, idx=0):
        with open(fn) as srcfd:
            src = csv.reader(srcfd, dialect="unix")
            self.hdr = next(src)
            def voldata(xs):
                "last element of each list is an integer"
                # volumes in nanolitres
                v = [[s.strip() for s in map(str,x[:-1])]+[int(float(x[-1].strip())*1000)] for x in xs] # contents in nanolitres
                return v
            self.data = dict((k.strip(), voldata(gs)) for k, gs in \
                    itertools.groupby(sorted(src, key=lambda x:x[idx]), key=lambda x:x[idx]) if k.strip()!='')
       
    
    def xfersrc(self, l, vol, depleted):
        """ transfer data for liquid l - reduces volume """
        dx = self.data
        while dx[l] and dx[l][0][-1]-vol<util.DEAD_VOLS[dx[l][0][2]]:
            dx[l].pop(0) # discard depleted wells
        # assert dx[l], "primer depleted for "+l
        if not dx[l]:
            depleted.update([l])
            return ['']*4
        dx[l][0][-1] -= vol
        return dx[l][0][:4]
    

class PicklistMemorySrc:
    """ read a survey & contents filehandle (StringIO) """ 

    def __init__(self, srcfd, idx=0):
        src = csv.reader(srcfd, dialect="unix")
        self.hdr = next(src)
        def voldata(xs):
            "last element of each list is an integer"
            # volumes in nanolitres
            v = [[s.strip() for s in map(str,x[:-1])]+[int(float(x[-1].strip())*1000)] for x in xs] # contents in nanolitres
            return v
        self.data = dict((k.strip(), voldata(gs)) for k, gs in \
                itertools.groupby(sorted(src, key=lambda x:x[idx]), key=lambda x:x[idx]) if k.strip()!='')
       
    
    def xfersrc(self, l, vol, depleted):
        """ transfer data for liquid l - reduces volume """
        dx = self.data
        while dx[l] and dx[l][0][-1]-vol<util.DEAD_VOLS[dx[l][0][2]]:
            dx[l].pop(0) # discard depleted wells
        # assert dx[l], "primer depleted for "+l
        if not dx[l]:
            depleted.update([l])
            return ['']*4
        dx[l][0][-1] -= vol
        return dx[l][0][:4]
    

def grouper(xs, kf=lambda x:x[0]):
    "group pairs: (x,a), (x,b), ... => (x,[a,b, ...]), ..."
    return itertools.groupby(sorted(xs, key=kf), key=kf)


def mk_picklist(exp, fn, rows, transactions, output_plate_guards=False, caller_id=None):
    """ 
    Output an Echo picklist given the header and rows (transfer spec) 
    Inputs:
        exp - an Experiment instance
        fn - output filename
        rows - data rows to write
        transaction - transaction dictionary of changes from the input plates
        output_plate_guards - whether or not to write guard characters on the plate barcodes
    """
    #try:
    #print("mk_picklist", file=sys.stderr)
    
    if True:
        plhdr_str = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"+\
                "Destination Plate Name,Destination Plate Barcode,Destination Plate Type,"+\
                "Destination Well,Volume"
        plhdr = plhdr_str.split(',')
        if fn not in transactions:
            transactions[fn] = {}
        with open(fn, "wt", newline='') as dstfd:
            dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            dst.writerow(plhdr)
            data = []
            for row in rows:
                if len(row) != len(plhdr):
                    m(f"entry does not match expected header format: {row=} {plhdr=}", level='warning', caller_id=caller_id)
                    continue
                src_plate_barcode = row[1]
                src_plate_well = row[3]
                vol = row[-1]
                #print(f"{src_plate_barcode=} {fn=} {transactions=}", file=sys.stderr)
                if src_plate_barcode not in transactions[fn]:
                    transactions[fn][src_plate_barcode] = {}
                if src_plate_well not in transactions[fn][src_plate_barcode]:
                    transactions[fn][src_plate_barcode][src_plate_well] = 0
                transactions[fn][src_plate_barcode][src_plate_well] -= vol

                if output_plate_guards:  # guard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() \
                            else util.guard_pbc(str(d), silent=True) for (h,d) in zip(plhdr, row)])
                else:  # unguard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() \
                            else util.unguard_pbc(str(d), silent=True) for (h,d) in zip(plhdr, row)])
            dst.writerows(data)
    #except Exception as exc:
    #    m(f"Picklist generation failed for {fn} {exc}", level='error')
    #    if Path(fn).exists():
    #        os.remove(fn)
    #    transactions = None
    #    exp.save()
    #    return False
    m(f'Wrote picklist to {fn}', level='success', caller_id=caller_id)
    return True


def mk_mytaq_picklist(exp, fn, task_wells, pcrPlate_col, pcrWell_col, taqwater_bcs, taq_vol, water_vol, transactions, caller_id=None):
    """ 
    Create a picklist for Mytaq & H2O transfers
    transactions are passed through to mk_picklist and modified there
    Because it is called in both PCR1 and PCR2, we should pass through the column indices we need
    """     
    ppcol = pcrPlate_col
    pwcol = pcrWell_col
    # build a list of all possible taq/water pids and well pair transfers
    pid_ww_tw_list = []
    #print(f"{task_wells=} {taqwater_bcs=}", file=sys.stderr)
    for pid in taqwater_bcs:
        twp = transaction.get_plate(exp, pid, transactions)
        #print(f"{twp=}", file=sys.stderr)
        # iterate over wells forever
        #ww_gen = (x for xs in itertools.repeat(twp['water_wells']) for x in xs)
        #tw_gen = (x for xs in itertools.repeat(twp['taq_wells']) for x in xs)
        empty_wells = set()
        water_well_list = []
        for ww in (x for xs in itertools.repeat(twp['water_wells']) for x in xs):
            if len(empty_wells) == 3:
                break
            if 'volume' not in twp[ww]:
                empty_wells.add(ww)
            if ww in empty_wells:
                continue
            if twp[ww]['volume'] - water_vol - twp['dead_vol'] < 0:
                empty_wells.add(ww)
                continue
            twp[ww]['volume'] -= water_vol
            water_well_list.append((pid,ww))

        empty_wells = set()
        taq_well_list = []
        for tw in (x for xs in itertools.repeat(twp['taq_wells']) for x in xs):
            if len(empty_wells) == 3:
                break
            if 'volume' not in twp[tw]:
                empty_wells.add(tw)
            if tw in empty_wells:
                continue
            if twp[tw]['volume'] - taq_vol - twp['dead_vol'] < 0:
                empty_wells.add(tw)
                continue
            twp[tw]['volume'] -= water_vol
            taq_well_list.append((pid,tw))

        for ww,tw in zip(water_well_list, taq_well_list):
            water_pid, water_well = ww
            taq_pid, taq_well = tw
            if water_pid != taq_pid:
                continue
            pid_ww_tw_list.append((water_pid,water_well,taq_well))

    if len(pid_ww_tw_list) < len(task_wells):
        m(f'Not enough taq or water available from plates {taqwater_bcs}', level='failure', caller_id=caller_id)
        return False

    # now allocate to used wells
    #try:
    if True:
        pcr_pids = sorted(set([w[ppcol] for w in task_wells]))  # w['pcrPlate']
        tw_pids = sorted(set([pid for pid,ww,tw in pid_ww_tw_list]))
        dst_dict = dict((pid, "Destination[{}]".format(i)) for i, pid in enumerate(pcr_pids, start=1))
        src_dict = dict((pid, "Source[{}]".format(i))for i, pid in enumerate(tw_pids, start=1))

        output_rows = []
        for task_well,pid_ww_tw in zip(task_wells,pid_ww_tw_list):
            # Need to generate this
            # plhdr_str = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"+\
            #        "Destination Plate Name,Destination Plate Barcode,Destination Plate Type,"+\
            #        "Destination Well,Volume"
            tw_pid, ww, tw = pid_ww_tw
            water_row = [src_dict[tw_pid], util.unguard(tw_pid, silent=True), util.PLATE_TYPES['Echo6'], ww,
                    dst_dict[task_well[ppcol]], util.unguard(task_well[ppcol], silent=True), 
                    util.PLATE_TYPES['PCR384'], task_well[pwcol], water_vol]
            output_rows.append(water_row)
            taq_row = [src_dict[tw_pid], util.unguard(tw_pid, silent=True), util.PLATE_TYPES['Echo6'], tw,
                    dst_dict[task_well[ppcol]], util.unguard(task_well[ppcol], silent=True), 
                    util.PLATE_TYPES['PCR384'], task_well[pwcol], taq_vol]
            output_rows.append(taq_row)
        mk_picklist(exp, fn, output_rows, transactions)         
    #except Exception as exc:
    #    m(f"failed to generate taq/water picklist {exc}", level='error')
    #    exp.save()
    #    return False
    return True


def allocate_primers_to_dna(exp, dna_records, pmr_pids=None, rotate=True, force=False, caller_id=None):
    """
    dna_records is a list of dictionaries

    args:
        exp - Experiment
        dna_records - a list of dictionaries with the following keys (below)
        pmr_pids - list of primer plate IDs
        rotate - if True, move to next available primer well after each dose
        caller_id - for logging
        force - if True, set the primerPlate and primerWell to None if the primer is not available

    returns:
        success - True if all primers were allocated, False if any primer was
                not available (unless force is True)
    dna_fields=['samplePlate','sampleWell','sampleBarcode','strain','sex','alleleSymbol',
                 'alleleKey','assayKey','assays','assayFamilies','clientName','sampleName',
                 'dnaPlate','dnaWell','primer']
        
    We want to add 'primerPlate' and 'primerWell' to the this dictionary
    if rotate is True, move to next available primer well after each dose
    Note: we carry on if the primer is not available, but we log a warning
    We only fail if the primer was available but runs out (unless force is True)
    """
    success = True
    primer_plate_well_vol_doses = exp.get_available_primer_wells(pmr_pids=pmr_pids)
    primer_uses = {pmr:0 for pmr in primer_plate_well_vol_doses}

    for i,rec in enumerate(dna_records):
        pmr = rec['primer']
        if pmr not in primer_plate_well_vol_doses:
            #print(f'{i=} {pmr=}') 
            m(f'primer {pmr} not available on primer plate(s), ', level='warning', caller_id=caller_id)
            dna_records[i]['primerPlate'] = None
            dna_records[i]['primerWell'] = None
            continue
        if rotate:
            offset = primer_uses[pmr] % len(primer_plate_well_vol_doses[pmr]) 
            while True:
                if sum([doses for pid, well, vol, doses in primer_plate_well_vol_doses[pmr]]) == 0:
                    # all wells for this primer are empty
                    m(f'primer {pmr} has run out of available doses', level='warning', caller_id=caller_id)
                    if force:
                        dna_records[i]['primerPlate'] = None
                        dna_records[i]['primerWell'] = None
                        m(f'Ignoring missing primer {pmr} for {i=}', level='warning', caller_id=caller_id)
                        break
                    else:
                        m(f'stopping allocation of primers to DNA records, use Force to complete or top up well', 
                                level='warning', caller_id=caller_id)
                        success = False
                        break
                pid, well, vol, doses = primer_plate_well_vol_doses[pmr][offset]
                if doses > 0:
                    primer_uses[pmr] += 1
                    primer_plate_well_vol_doses[pmr][offset][-1] -= 1
                    primer_plate_well_vol_doses[pmr][offset][-2] -= exp.transfer_volumes['PRIMER_VOL']
                    dna_records[i]['primerPlate'] = pid
                    dna_records[i]['primerWell'] = well
                    break
                else:
                    offset += 1
                    if offset == len(primer_plate_well_vol_doses[pmr]):
                        offset == 0
        else:  # use up each well before moving on
            offset = 0
            while True:
                if sum([doses for pid, well, vol, doses in primer_plate_well_vol_doses[pmr]]) == 0:
                    # all wells for this primer are empty
                    m(f'primer {pmr} has run out of available doses', level='warning', caller_id=caller_id)
                    if force:
                        m(f'Ignoring missing primer {pmr} for {i=}', level='warning', caller_id=caller_id)
                        dna_records[i]['primerPlate'] = None
                        dna_records[i]['primerWell'] = None
                        break
                    else:
                        m(f'stopping allocation of primers to DNA records, use Force to complete or top up well', 
                                level='warning', caller_id=caller_id)
                        success = False
                        break
                pid, well, vol, doses = primer_plate_well_vol_doses[pmr][offset]
                if doses > 0:
                    primer_uses[pmr] += 1  # not actually used in this mode 
                    primer_plate_well_vol_doses[pmr][offset][-1] -= 1
                    primer_plate_well_vol_doses[pmr][offset][-2] -= exp.transfer_volumes['PRIMER_VOL']
                    dna_records[i]['primerPlate'] = pid
                    dna_records[i]['primerWell'] = well
                    break
                else:
                    offset += 1
    return success


def stage_write(exp, fn, header, data, output_plate_guards=False, ignore_missing_plate_info=False, caller_id=None):
    """ 
    Output stage files
    Check which columns contain plate barcodes, and enforce plate guards or their removal based on output_plate_guards
    """
    try:
        with open(fn, "wt", newline='') as dstfd:
            wx = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_ALL)
            out_lines = []
            for row in data:
                if not ignore_missing_plate_info:
                    failed_row = False
                    for i, (h,d) in enumerate(zip(header, row)):
                        if h in ['primerPlate','dnaPlate','sourcePlate','pcrPlate'] and not d:
                            m(f'Failed to write row: col {i=} for {row=}', level='warning', caller_id=caller_id)
                            failed_row = True
                    if failed_row:
                        continue
                if output_plate_guards:  # guard plate barcodes
                    out_lines.append([util.guard_pbc(str(d), silent=True) \
                        if h in ['primerPlate','dnaPlate','sourcePlate','pcrPlate'] \
                        and d != '' else str(d) for (h,d) in zip(header, row)])
                else:  # unguard plate barcodes
                    out_lines.append([str(d) if h not in ['primerPlate','dnaPlate','sourcePlate','pcrPlate'] 
                            else util.unguard_pbc(str(d), silent=True) for (h,d) in zip(header, row)])
            wx.writerow(header)  
            wx.writerows(out_lines)
    except Exception as exc:
        m(f'could not write to {fn}. {exc}', level='error', caller_id=caller_id)
        return False
    return True


def generate_echo_PCR1_picklist(exp, dna_plate_bcs, pcr_plate_bcs, taq_water_bcs, force=False, caller_id=None):
    """
    Entry point. Takes an experiment instance plus plate barcodes for dna plates, PCR plates, primer plates, taq/water plates.
    """
    #try:
    if not exp.primer_assayfam:
        m('assay list must be loaded before generating primer picklist', level='failure', caller_id=caller_id)
        return False
    if True:
        dna_bcs, dna_success = util.guard_pids(dna_plate_bcs, caller_id=caller_id)
        pcr_bcs, pcr_success = util.guard_pids(pcr_plate_bcs, caller_id=caller_id)
        taq_bcs, taq_success = util.guard_pids(taq_water_bcs, caller_id=caller_id)
        if not dna_success or not pcr_success or not taq_success:
            return False

        primer_pids = exp.get_primer_pids()
        primer_survey_lines = []
        header = ['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'plate position on Echo 384 PP',
                'primer names pooled', 'volume']
        primer_survey_lines.append(','.join(header))
        for i,pid in enumerate(primer_pids):
            plate=exp.plate_location_sample[pid]
            for well in util.row_ordered_384:
                if well not in plate['wells']:
                    continue
                if 'primer' not in plate[well]:
                    continue
                if 'volume' not in plate[well]:
                    continue
                outline = ','.join([f"Source[{i+1}]",util.unguard_pbc(pid,silent=True),plate['plate_type'],
                        well,plate[well]['primer'],f"{int(plate[well]['volume'])/1000}"])
                primer_survey_lines.append(outline)

        with open(exp.get_exp_fn('primer-svy.csv'), 'wt') as fout:
            for psl in primer_survey_lines:
                print(psl, file=fout)
        
        transactions = {}
        outfmt = f"PCR-picklist_{exp.name}.csv"
        # do we also add taq/water plate, taq well, water well?
        pcr1_fields=['samplePlate','sampleWell','sampleBarcode','strain','sex','alleleSymbol',
                 'alleleKey','assayKey','assays','assayFamilies','clientName','sampleName',
                 'dnaPlate','dnaWell','primer','primerPlate','primerWell','pcrPlate','pcrWell']

        dna_records = exp.get_dna_records(dna_bcs)
        success = allocate_primers_to_dna(exp, dna_records,force=force, caller_id=caller_id)
        if not success:
            return False

        dna_records = sorted([d for d in dna_records if d['dnaPlate'] and d['primerPlate']], 
                key=lambda x: (x['samplePlate'],util.padwell(x['sampleWell'])))
                
        # allocate PCR plate wells
        wells = [r+str(c+1) for c in range(24) for r in"ABCDEFGHIJKLMNOP"]
        pcrWells = [(p, w) for p in sorted(pcr_bcs) for w in wells]
        pcr_records = []
        for i, (dr, pcrw) in enumerate(zip(dna_records, pcrWells)):
            record = [dr.get(pf, '') for pf in pcr1_fields[:-2]]  # don't add columns yet for PCR plate and well
            record.append(pcrw[0])
            record.append(pcrw[1])
            pcr_records.append(record)

        # output Stage 2 CSV file - used in Stage 3 below - keeping the plate guards for the next stage (to be safe)
        stage2_fn = exp.get_exp_fn("Stage2.csv", trans=True)
        transactions[stage2_fn] = {}
        success = stage_write(exp, stage2_fn, pcr1_fields, pcr_records, output_plate_guards=True)
        if not success:
            transactions[stage2_fn] = None
            m(f'failed to write {stage2_fn}', level='error', caller_id=caller_id)
            return False
        
        pcr_plate_col = pcr1_fields.index('pcrPlate')
        pcr_well_col = pcr1_fields.index('pcrWell')
        dna_plate_col = pcr1_fields.index('dnaPlate')
        dna_well_col = pcr1_fields.index('dnaWell')
        primer_plate_col = pcr1_fields.index('primerPlate')
        primer_well_col = pcr1_fields.index('primerWell')

        # output DNA picklist
        dna_fn = exp.get_exp_fn(outfmt.replace('PCR','PCR1_dna'), trans=True)
        #plhdr_str = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"+\
        #        "Destination Plate Name,Destination Plate Barcode,Destination Plate Type,"+\
        #        "Destination Well,Volume"
        dst_plate_type = util.PLATE_TYPES['PCR384']
        src_plate_type = util.PLATE_TYPES['Echo384']

        pcr_pid_list = list((set([r[pcr_plate_col] for r in pcr_records if r[pcr_plate_col]])))
        dna_bcs = list((set([r[dna_plate_col] for r in pcr_records if r[dna_plate_col]])))
        primer_pid_list = list((set([r[primer_plate_col] for r in pcr_records if r[primer_plate_col]])))

        ddict = dict((pbc, f"Destination[{i}]") for i, pbc in enumerate(pcr_pid_list, start=1))
        sdict = dict((pbc, f"Source[{i}]")for i, pbc in enumerate(dna_bcs, start=1))
        volume = exp.transfer_volumes['DNA_VOL']
        dna_gen = ((sdict[util.guard_pbc(r[dna_plate_col], silent=True)],
                util.guard_pbc(r[dna_plate_col], silent=True), 
                src_plate_type, r[dna_well_col], 
                ddict[util.guard_pbc(r[pcr_plate_col], silent=True)], 
                util.guard_pbc(r[pcr_plate_col], silent=True), dst_plate_type, 
                r[pcr_well_col], volume) for r in pcr_records if r[dna_plate_col] and r[pcr_plate_col])
        success = mk_picklist(exp, dna_fn, dna_gen, transactions)
        if not success:
            transactions[dna_fn] = None 

        # primer picklist
        primer_fn = exp.get_exp_fn(outfmt.replace('PCR','PCR1_primer'), trans=True)
        transactions[primer_fn] = {}
        
        ddict = dict((pbc, f"Destination[{i}]") for i, pbc in enumerate(pcr_pid_list, start=1))
        sdict = dict((pbc, f"Source[{i}]")for i, pbc in enumerate(primer_pid_list, start=1))
        volume = exp.transfer_volumes['PRIMER_VOL']
        primer_gen = ((sdict[util.guard_pbc(r[primer_plate_col], silent=True)],
                util.guard_pbc(r[primer_plate_col], silent=True),
                src_plate_type, r[primer_well_col], 
                ddict[util.guard_pbc(r[pcr_plate_col], silent=True)],
                util.guard_pbc(r[pcr_plate_col], silent=True), dst_plate_type, 
                r[pcr_well_col], volume) for r in pcr_records if r[primer_plate_col] and r[pcr_plate_col])
        success = mk_picklist(exp, primer_fn, primer_gen, transactions)
        if not success:
            transactions[primer_fn] = None      

        # also PCR1 water and Taq
        taq_fn = exp.get_exp_fn(outfmt.replace('PCR','PCR1_taqwater'), trans=True)
        transactions[taq_fn] = {}
        # exp, fn, task_wells, taqwater_bcs, taq_vol, water_vol, transactions=None
        success = mk_mytaq_picklist(exp, taq_fn, pcr_records, pcr_plate_col, pcr_well_col, 
                taq_bcs, exp.transfer_volumes['PRIMER_TAQ_VOL'], 
                exp.transfer_volumes['PRIMER_WATER_VOL'], transactions)
        if not success:
            transactions[taq_fn] = None

    #except Exception as exc:
    #    m(f"PCR1 Echo picklists could not be created {exc}", level='error')
    #    return False
    transaction.add_pending_transactions(exp, transactions)
    m(f'PCR1 Echo picklists created', level='success', caller_id=caller_id)
    return True


def i7i5alloc_rot(exp, vol, wellcount, index_survey_fn='index-svy.csv', caller_id=None):
    """ allocate vol to wellcount wells with unique index pairs - rotates barcode wells to avoid wells being drained """
    #try:
    if True:
        typebc = util.Table.newtype('BCRecord', "name platebc type well set barcode orderpart oligo volume")
        # typebc = Table.newtype('BCRecord', "well name index indexName oligo volume")
        tab = util.CSVTable(typebc, exp.get_exp_fn(index_survey_fn)) # volumes in file are in uL

        dv = util.DEAD_VOLS[util.PLATE_TYPES['Echo384']]
        #def getvol(x): return x[-1]
        i7s = {(r.barcode, r.set, r.well, r.platebc): float(r.volume)*1000-dv for r in tab.data if '_i7F_' in r.set}
        i5s = {(r.barcode, r.set, r.well, r.platebc): float(r.volume)*1000-dv for r in tab.data if '_i5R_' in r.set}
        i7_idxs = set([i7[0] for i7 in i7s.keys()])
        i5_idxs = set([i5[0] for i5 in i5s.keys()])

        i7gen = (x for xs in itertools.repeat(i7s.keys()) for x in xs)
        i5gen = (x for xs in itertools.repeat(i5s.keys()) for x in xs)
        i7s_empty = set()
        i5s_empty = set()
        index_info_pairs = []
        used_index_pairs = set()
        combos = len(i7_idxs) * len(i5_idxs)
        if combos < wellcount:
            m(f"Not enough barcode combos {combos} for wells requested {wellcount}", level='failure', caller_id=caller_id)
            return
        while len(index_info_pairs) < wellcount:     
            if len(i7s_empty) == len(i7s):
                m(f"i7 wells contain insufficient volume for the number of experiment wells", level='failure', caller_id=caller_id)
                m(f"{i7s_empty=}, {i5s_empty=}, {len(used_index_pairs)=}", level='debug', caller_id=caller_id)
                return
            elif len(i5s_empty) == len(i5s):
                m(f"i5 wells contain insufficient volume for the number of experimental wells", level='failure', caller_id=caller_id)
                m(f"{i7s_empty=}, {i5s_empty=}, {len(used_index_pairs)=}", level='debug', caller_id=caller_id)
                return
            i7 = next(i7gen)
            i5 = next(i5gen)
            if (i7[0],i5[0]) in used_index_pairs:
                if len(i7s) > len(i5s):
                    i7 = next(i7gen)
                else:
                    i5 = next(i5gen)
                if (i7[0],i5[0]) in used_index_pairs:
                    continue
            i7s[i7] -= vol
            if i7s[i7] < 0:
                i7s_empty.add(i7)
                i7s[i7] += vol  # since we didn't actually use the last of it
                continue
            i5s[i5] -= vol
            if i5s[i5] < 0:
                i5s_empty.add(i5)
                i5s[i5] += vol  # since we didn't actually use the last of it
                continue

            index_info_pairs.append((i7,i5))
            used_index_pairs.add((i7[0],i5[0])) 
        return index_info_pairs
      
    #except Exception as exc:
    #    m(f"{exc}", level='error')


def generate_miseq_samplesheet(exp, miseq_fn, s3tab,transactions, caller_id=None):
    """ Now that we have all the required info, generate the Illumina Miseq samplesheet, which drives the sequencing """
    if miseq_fn not in transactions:
        transactions[miseq_fn] = {}
    #try:  
    if True:
        with open(miseq_fn, "wt", newline='') as dstfd:
            dstfd.write(f"""[Header],,,,,,,,,,
IEMFileVersion,4,,,,,,,,,
Investigator Name,,,,,,,,,,
Experiment Name,{'NGSG:'+exp.name},,,,,,,,,
Date,,,,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,,
Application,FASTQ Only,,,,,,,,,
Assay,Nextera XT,,,,,,,,,
Description,Mouse genotyping,,,,,,,,,
Chemistry,Amplicon,,,,,,,,,
,,,,,,,,,,
[Reads],,,,,,,,,,
151,,,,,,,,,,
151,,,,,,,,,,
,,,,,,,,,,
[Settings],,,,,,,,,,
ReverseComplement,0,,,,,,,,,
Adapter,,,,,,,,,,
,,,,,,,,,,
[Data],,,,,,,,,,
""")
            hdr = "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description".split(',')
            dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            dst.writerow(hdr)
            gen = [(util.unguard_pbc(r.pcrPlate, silent=True)+'_'+util.padwell(r.pcrWell), 
                    util.unguard(r.sampleBarcode, silent=True), '', '',
                    r.i7name, r.i7bc,r.i5name, r.i5bc,'NGSgeno') for r in s3tab.data]
            
            dst.writerows(gen)
        
    #except Exception as exc:
    #    m(f"Miseq samplesheet generation {fnmiseq} will be deleted {exc}", level='error')
    #    if Path(fnmiseq).exist():
    #        os.remove(fnmiseq)
    #    exp.save()
    #    return False
    #exp.add_pending_transactions(transactions)
    m(f'Miseq samplesheet written to {miseq_fn}', level='success', caller_id=caller_id)
    return True

    
def generate_echo_PCR2_picklist(exp, pcr_plate_bcs, index_plate_bcs, 
        taq_water_bcs, amplicon_plate_bcs=None, caller_id=None):
    """
    Replacement for main(). Called as a library only.
    User included PCR plates, index plates, taq/water plates, 
    amplicon plates (optional) are user provided
    Also calls generate_miseq_samplesheet()
    Called from Experiment.generate_echo_PCR2_picklists()
    """
    transactions = {}
    #try:
    if True:
        pcr_bcs, pcr_success = util.guard_pids(pcr_plate_bcs, caller_id=caller_id)
        index_bcs, index_success = util.guard_pids(index_plate_bcs, caller_id=caller_id)
        taq_bcs, taq_success = util.guard_pids(taq_water_bcs, caller_id=caller_id)
        amplicon_bcs, amplicon_success = util.guard_pids(amplicon_plate_bcs, caller_id=caller_id)
        
        outfmt = f"PCR-picklist_{exp.name}.csv"
    
        fnstage2 = "Stage2.csv"
        fnstage3 = exp.get_exp_fn("Stage3.csv", trans=True)
        total_wells = 0
        # we need to decompose S2tab so that we can add the amplicon rows to it
        if amplicon_bcs and not pcr_bcs:
            s2_header = ['samplePlate', 'sampleWell', 'sampleBarcode', 'strain', 'sex', 'alleleSymbol',
                                  'alleleKey', 'assayKey', 'assays', 'assayFamilies', 'clientName', 'sampleName',
                                  'dnaPlate', 'dnaWell', 'primer', 'primerPlate', 'primerWell', 'pcrPlate', 'pcrWell']
            s2_data = []
        else:
            try:
            
                s2tab = util.CSVTable('S2Rec', exp.get_exp_fn(fnstage2))
            except FileNotFoundError as fne:
                m(f"Stage 2 file has not been generated from Primer stage. {fne}", level='error', caller_id=caller_id)
                return False
            
 
            s2_header = s2tab.header
            s2_data = s2tab.data # list of S3 records, each is a namedtuple
        total_wells += len(s2_data)
        #print(f'{total_wells=}')      
        
        def check_quotes(string:str) -> str:
            if string.startswith(("'", '"')) and not string.endswith(("'", '"')):
                return string + string[0]
            elif not string.startswith(("'", '"')) and string.endswith(("'", '"')):
                return string[0] + string
            else:
                return string        

        index_vol = exp.transfer_volumes['INDEX_VOL']  # 175 nanolitres
        
        for amp_pid in amplicon_bcs:
            total_wells += len(exp.plate_location_sample[amp_pid]['wells'])
        
        index_alloc = i7i5alloc_rot(exp, index_vol, total_wells)
        #print(f'{len(index_alloc)=}')
        s2_data_rows = []
        for s2_record in s2_data:
            s2_data_rows.append([s2_record[i] for i,h in enumerate(s2_header)])
        #print(f'{s2_data_rows=}')
        #print(f'{s2_data_rows=}', file=sys.stderr)
        for amp_pid in amplicon_bcs:
            for well_str in exp.plate_location_sample[amp_pid]['wells']:
                amp_well = exp.plate_location_sample[amp_pid][well_str]
                amp_row = [# amp_well['sampleNumber'],      # sampleNumber - no longer used
                        amp_pid,                        # samplePlate
                        well_str,                       # sampleWell
                        amp_well['barcode'],            # sampleBarcode
                        '',                             # strain
                        '',                             # sex
                        amp_well['amplicons'][0],       # alleleSymbol
                        '',                             # alleleKey
                        '',                             # assayKey
                        check_quotes(amp_well['amplicons'][0].split('_')[0]), # assay
                        check_quotes(amp_well['amplicons'][0].split('_')[0]), # assayFamily
                        amp_well.get('clientName', ''),     # clientName
                        amp_well.get('sampleName', ''),     # sampleName
                        amp_pid,                        # dnaPlate
                        well_str,                       # dnaWell
                        amp_well['amplicons'][0],       # primer
                        '',                             # primerPlate
                        '',                             # primerWell
                        amp_pid,                        # pcrPlate
                        well_str]                       # pcrWell
                s2_data_rows.append(amp_row)
        #print(s2_data_rows)
        # convert to stream of characters from row*column lists
        s2hdr_str = ','.join(s2_header)
        s2amp_stream = StringIO(s2hdr_str + '\n' + '\n'.join([','.join(row) for row in s2_data_rows]))
        s2amp_tab = util.CSVMemoryTable('S2Rec', s2amp_stream) 
        #print(f'{len(s2amp_tab.data)=}')
        s3flds = ['sampleNo'] + list(s2amp_tab.tt._fields+('i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well', 'index_plate'))
        #print(f'{s3flds=}', file=sys.stderr)
        S3Rec = util.Table.newtype('S3Rec', s3flds)
        # This adds all the index info to the Stage2.csv file and will save it as Stage3.csv
        s3rows = []
        counter = 1
        index_alloc = i7i5alloc_rot(exp, index_vol, total_wells)
        #print(f'{len(s2amp_tab.data)=} {len(index_alloc)=}')
        for p1, (p2, p3) in zip(s2amp_tab.data, index_alloc):
            row = [counter]  # new sampleNo
            for p in p1:
                row.append(p)
            for p in p2[:3]:  # i7bc, i7name, i7well
                row.append(p)
            for p in p3[:4]:  # i5bc, i5name, i5well, index_plate
                row.append(p)
            s3rows.append(row)
            counter += 1

        #print(f'{len(s3rows)=}', file=sys.stderr)
        s3tab = util.Table(S3Rec, s3rows, headers=s3flds)
        #s3tab = util.Table(S3Rec, ([x for xs in (p1, p2[:3], p3[:4]) for x in xs] for p1, (p2, p3) in zip(s2amp_tab.data, index_alloc)), headers=s3flds)
        
        #s3tab = util.Table(S3Rec, ([x for xs in (p1, p2[:3], p3[:4]) for x in xs] for p1, (p2, p3) in zip(s2tab.data, index_alloc)), headers=s3flds) 
        # output Stage 3 CSV file - used for custom and mouse samples. Amplicons are handled separately
        #print(f'{len(s3tab.data)=}')
        transactions[fnstage3] = {}
        try:
            s3tab.csvwrite(fnstage3, output_plate_guards=True)
        except Exception as exc:
            m(f'could not write {fnstage3}. {exc}', level='error', caller_id=caller_id)
            transactions[fnstage3] = None
            return False

        # write index picklists
        # make source dictionary for index plates
        src_dict = {}
        for r in s3tab.data:
            if r.index_plate not in src_dict:
                src_dict[r.index_plate] = f"Source[{len(src_dict)+1}]"
        # make the destination dictionary for pcrPlates
        dest_dict = {}
        for r in s3tab.data:
            if r.pcrPlate not in dest_dict:
                dest_dict[r.pcrPlate] = f"Destination[{len(dest_dict)+1}]"

        pt_src = util.PLATE_TYPES['Echo384']
        pt_dst = util.PLATE_TYPES['PCR384']
        rowi7 = ((src_dict[r.index_plate], r.index_plate, pt_src, r.i7well, dest_dict[r.pcrPlate], r.pcrPlate, pt_dst, r.pcrWell, index_vol) for r in s3tab.data)
        rowi5 = ((src_dict[r.index_plate], r.index_plate, pt_src, r.i5well, dest_dict[r.pcrPlate], r.pcrPlate, pt_dst, r.pcrWell, index_vol) for r in s3tab.data)
        outfmt_index = outfmt
        fn_index = exp.get_exp_fn(outfmt_index.replace('PCR','PCR2_index'), trans=True)
        mk_picklist(exp, fn_index, (r for rs in (rowi7, rowi5) for r in rs), transactions)    
        
        # also PCR2 water and Taq
        outfmt_taq = outfmt
        fn_taqwater = exp.get_exp_fn(outfmt_taq.replace('PCR','PCR2_taqwater'), trans=True)
        transactions[fn_taqwater] = {}
        pcr_plate_col = s3tab.header.index('pcrPlate') 
        pcr_well_col = s3tab.header.index('pcrWell')
        mk_mytaq_picklist(exp, fn_taqwater, s3tab.data, pcr_plate_col, pcr_well_col,
                taq_bcs, exp.transfer_volumes['INDEX_TAQ_VOL'], 
                exp.transfer_volumes['INDEX_WATER_VOL'], transactions)
        
        # MiSeq file - experiment name, name format, ngid(dir), s3 data, verbosity
        fmtmiseq = "MiSeq_{}.csv".format(exp.name)
        miseq_fn = exp.get_exp_fn(fmtmiseq, trans=True)
        transactions[miseq_fn] = {}
        success = generate_miseq_samplesheet(exp, miseq_fn, s3tab, transactions)
    transaction.add_pending_transactions(exp, transactions)
    
    #except Exception as exc:
    #    m(f"{exc}", level='error')
    #    return False
    return True


def myopen(fn):
    """ Bob's function to handle gzip transparently """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt")
    return open(fn, errors="ignore")


def read_html_template(html_fn):
    """ 
    Read an html file from the library folder, ready for additional fields 
    Used by cgi-dashboard.py and makehtml.py
    template uses {!field!} instead of python format fields - easier to edit.
    Files should use the .tpl (template) extension as they are not strictly html
    """
    tfn = os.path.join('..', 'library', html_fn)
    with open(tfn) as src:
        htmlcode = src.read()
    # template uses {!field!} instead of python format fields - easier to edit.
    htmlfmt = htmlcode.replace('{', '{{').replace('}', '}}').replace('{{!', '{').replace('!}}', '}')
    return htmlfmt


def match_nimbus_to_echo_files(exp, caller_id=None):
    """ Make sure there is a 1-1 match of Nimbus to Echo-COC file.
        Assumes plate barcodes are guarded.
        Inputs: an Experiment instance
        Outputs: three lists (nimbus files, echo files, barcodes not matches by echo files)
   """
    try:
        nimbus_files = list(Path(exp.get_exp_dn()).glob('Nimbus-*.csv'))
        echo_files = list(Path(exp.get_exp_dn('uploads')).glob('Echo_384_COC_00??_*_0.csv'))
        #print(f'{nimbus_files=}')
        #print(f'{echo_files=}')
        # cleave off the path and file name junk to just get the barcodes
        nbc = frozenset(fp.name.replace('Nimbus-','').replace('.csv','') for fp in nimbus_files)
        ebc = frozenset(fp.name.replace('.csv','').split('_',5)[4] for fp in echo_files)
        #print(f'{nbc=}')
        #print(f'{ebc=}')
        xbc  = nbc-ebc # any Nimbus plate files that are missing Nimbus run (output) files
        return [str(fp) for fp in nimbus_files], [str(fp) for fp in echo_files], xbc
    except Exception as exc:
        m(f'Problem matching Nimbus to Echo files {exc}', level='error', caller_id=caller_id)
        m(f'{traceback.print_exc(limit=2)=}', level='debug', caller_id=caller_id)
        return [], [], []


def generate_targets(exp, caller_id=None):
    """ create target file based on loaded references """
    #transactions = {}
    #target_fn = exp.get_exp_fn('targets.fa', trans=True)
    target_fn = exp.get_exp_fn('targets.fa')
    #transactions[target_fn] = {} # add plates and modifications to this
    counter = 0
    try:
        with open(target_fn, 'wt') as targetf:
            for group in exp.reference_sequences:
                for id in exp.reference_sequences[group]:
                    print(f'>{id}', file=targetf)
                    print(f'{exp.reference_sequences[group][id]}', file=targetf)
                    counter += 1    
    except Exception as exc:
        m(f'could not write reference sequences to {target_fn} {exc}', level='error', caller_id=caller_id)
        return False
    #transaction.add_pending_transactions(exp, transactions)
    #transaction.accept_pending_transactions(exp)
    m(f'Screated reference sequences file {target_fn} containing {counter} sequences', level='success', caller_id=caller_id)
    return True


def generate_primer_assayfams(exp, caller_id=None):
    """ create a two column  csv file of primer and assayfamily """
    primer_fn = exp.get_exp_fn('primers.csv')
    try:
        with open(primer_fn, 'wt') as primerf:
            for af in exp.assayfam_primers:
                for pmr in exp.assayfam_primers[af]:
                    print(f'{pmr},{af}', file=primerf)
    except Exception as exc:
        m(f'count not write primer and assay families to {primer_fn} {exc}', level='error', caller_id=caller_id)
        return False
    m(f'created primer assay families file {primer_fn}', level='success', caller_id=caller_id)
    return True


if __name__ == '__main__':
    """ library only """
    pass