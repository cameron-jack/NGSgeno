#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Oct 2019
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: adjusted paths relative to app directory
@last_edit: 2022-02-29
@edit_comment: 

Nimbus input file generators

Input is 1-4 96-well sample (ear-punch) plate bar codes ... and a destination 384-well plate bar code.

The program uses either a DB webservice for each ear-punch plate, or has 
the required info passed in from CSV format and converted to JSON.
The webservice returns JSON containing a list of wells and mouse information 
(IDs, barcodes, expected "observables", ...).

The output format is CSV
The output asks the Nimbus Robot to tranfer a sample from each well of 
the sample plates to a well in a 384-well DNA plate. 
This is the first stage of the NGS Genotyping pipeline.

Sample barcodes should remain guarded
Plate barcodes need to be written out as unguarded so that they can match result of barcode scanners
"""

# can we validate barcode entry widgets as GUI focus leaves them, and colour their background accordingly?
# activate go button when ear-punch fields are valid

import os
import sys
import csv

import bin.file_io as file_io
from bin.util import output_error, col_ordered_96, nimbus_ordered_96
from bin.db_io import get_musterer_mouse_info, get_plate_musterer, get_plate_rodentity


app = None # global!

def read_custom_csv_as_table(input_fn, custom_file_contents='custom'):
    """ read a custom sample manifest, Then return as JSON 
        Creates a copy of the original manifest with guarded barcodes.
        v0.15 attempts to download mouse data from Musterer for valid mouse barcodes in manifest
        v0.16 adds DNA plate barcode (dest barcode) to allow for more than 4 sample plates per manifest file
    """
    try:
        data = {}
        errs = []  # error messages
        guarded_fn = input_fn.replace('.csv','_guarded.csv')
        ofh = None
        musterer_barcodes = []
        rodentity_barcodes = []
        dna_plate_samples = {}  # [dna_plate_id] = [sample_plate_barcodes]
        if not os.path.exists(guarded_fn):
            ofh = open(guarded_fn, 'wt')
            out_csv = csv.writer(ofh, dialect='unix')
        with open(input_fn, 'rt', newline=None) as f:
            in_csv = csv.reader(f, delimiter=',')
            for i, row in enumerate(in_csv):
                #print(row, file=sys.stdout)
                #print(file=sys.stdout)
                if i == 0:
                    if ofh:
                        out_csv.writerow(row)
                    continue  # header
                if all(col.strip() == '' for col in row):
                    continue  # blank rows
                if any(col.strip() == '' and j < 5 for j,col in enumerate(row)):
                    msg = f"ERROR. blank manifest entry in entry {i}: {row}"
                    print(msg, file=sys.stdout)
                    errs.append(msg)
                    continue
                    # return
                dnap = str(row[1].strip())
                if file_io.is_guarded_pbc(dnap):
                    gdnap = dnap
                else:
                    gdnap = guard_pbc(dnap)
                plt = str(row[2].strip())
                if file_io.is_guarded_pbc(plt):
                    gsampPID = plt
                else:
                    gsampPID = guard_pbc(plt)
                dest = str(row[3]).strip()
                if file_io.is_guarded(dest):
                    gdestPID = dest
                else:
                    gdestPID = file_io.guard_pbc(dest)
                well = str(row[4]).strip()
                smpl = str(row[5].strip())
                if file_io.is_guarded(smpl):
                    gsampleBarcode = smpl
                elif custom_file_contents == 'custom' or custom_file_contents == 'none':
                    gsampleBarcode = file_io.guard_cbc(smpl)
                elif custom_file_contents == 'musterer':
                    gsampleBarcode = file_io.guard_mbc(smpl)
                    musterer_barcodes.append(gsampleBarcode)
                elif custom_file_contents == 'rodentity':
                    gsampleBarcode = file_io.guard_rbc(smpl)
                    rodentity_barcodes.append(gsampleBarcode)
                else:
                    msg = f"ERROR. unknown custom file contents field: {custom_file_contents}"
                    print(msg, file=sys.stdout)
                    errs.append(msg)

                if ofh:
                    out_row = [r for r in row]
                    out_row[1] = gsampPID
                    out_row[2] = gdestPID
                    out_row[3] = gsampleBarcode
                    out_csv.writerow(out_row)
                assays = [a.strip() for a in row[5:] if a.strip() != '']
                if gsampPID not in data:
                    data[gsampPID] = {'plateId':gsampPID, 'destPID':gdestPID, 'custom':True, 'sampleBarcode':gsampleBarcode, 
                           'wells':[{'wellLocation':well, 'organism':{'sampleId':gsampleBarcode, 
                                    'sampleBarcode':gsampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a.strip()} \
                                    for a in assays]}}]}
                    dna_plate_samples[gdnap].append(gsampPID)
                else:
                    if well in data[gsampPID]:
                        print('duplicate well number', well, 'in plate', gsampPID, file=sys.stdout)
                        return
                    data[gsampPID]['wells'].append({'wellLocation':well, 'organism':{'sampleId':gsampleBarcode, 
                            'sampleBarcode':gsampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a.strip()} \
                            for a in assays]}})
    
        if ofh:
            ofh.close()

        if musterer_barcodes:
             musterer_mouse_info = get_musterer_mouse_info(musterer_barcodes)
             for gsampPID in data:
                 for i,well_info in enumerate(data[gsampPID]['wells']):
                     if well_info['organism']['sampleBarcode'] in musterer_mouse_info:
                         gmbc = well_info['organism']['sampleBarcode']
                         data[gsampPID]['wells'][i]['organism']['sampleId'] = musterer_mouse_info[gmbc]['mouseId']
                         data[gsampPID]['wells'][i]['organism']['strain'] = musterer_mouse_info[gmbc]['strain']
                         data[gsampPID]['wells'][i]['organism']['sex'] = musterer_mouse_info[gmbc]['sex']
        if rodentity_barcodes:
            pass
            # rodentity_mouse_info = get_rodentity_mouse_info(rodentity_barcodes)

        return data, dna_plate_samples, errs

        # below is what the old code needed returned
        #sample_pids = data.keys() # why this structure?
        #data = [(data[p],p) for p in sorted(data)]  # why?
        #return data, sample_pids, dna_plate_samples, errs
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.readCustomCSVtoJSON')

def read_custom_csv_to_json(input_fn, custom_file_contents='custom'):
    """ read a custom sample manifest, Then return as JSON 
        Creates a copy of the original manifest with guarded barcodes.
        v0.15 attempts to download mouse data from Musterer for valid mouse barcodes in manifest
        v0.16 adds DNA plate barcode (dest barcode) to allow for more than 4 sample plates per manifest file
    """
    try:
        data = {}
        errs = []  # error messages
        guarded_fn = input_fn.replace('.csv','_guarded.csv')
        ofh = None
        musterer_barcodes = []
        rodentity_barcodes = []
        dna_plate_samples = {}  # [dna_plate_id] = [sample_plate_barcodes]
        if not os.path.exists(guarded_fn):
            ofh = open(guarded_fn, 'wt')
            out_csv = csv.writer(ofh, dialect='unix')
        with open(input_fn, 'rt', newline=None) as f:
            in_csv = csv.reader(f, delimiter=',')
            for i, row in enumerate(in_csv):
                #print(row, file=sys.stdout)
                #print(file=sys.stdout)
                if i == 0:
                    if ofh:
                        out_csv.writerow(row)
                    continue  # header
                if all(col.strip() == '' for col in row):
                    continue  # blank rows
                if any(col.strip() == '' and j < 5 for j,col in enumerate(row)):
                    msg = f"ERROR. blank manifest entry in entry {i}: {row}"
                    print(msg, file=sys.stdout)
                    errs.append(msg)
                    continue
                    # return
                dnap = str(row[1].strip())
                if file_io.is_guarded_pbc(dnap):
                    gdnap = dnap
                else:
                    gdnap = guard_pbc(dnap)
                plt = str(row[2].strip())
                if file_io.is_guarded_pbc(plt):
                    gsampPID = plt
                else:
                    gsampPID = guard_pbc(plt)
                dest = str(row[3]).strip()
                if file_io.is_guarded(dest):
                    gdestPID = dest
                else:
                    gdestPID = file_io.guard_pbc(dest)
                well = str(row[4]).strip()
                smpl = str(row[5].strip())
                if file_io.is_guarded(smpl):
                    gsampleBarcode = smpl
                elif custom_file_contents == 'custom' or custom_file_contents == 'none':
                    gsampleBarcode = file_io.guard_cbc(smpl)
                elif custom_file_contents == 'musterer':
                    gsampleBarcode = file_io.guard_mbc(smpl)
                    musterer_barcodes.append(gsampleBarcode)
                elif custom_file_contents == 'rodentity':
                    gsampleBarcode = file_io.guard_rbc(smpl)
                    rodentity_barcodes.append(gsampleBarcode)
                else:
                    msg = f"ERROR. unknown custom file contents field: {custom_file_contents}"
                    print(msg, file=sys.stdout)
                    errs.append(msg)

                if ofh:
                    out_row = [r for r in row]
                    out_row[1] = gsampPID
                    out_row[2] = gdestPID
                    out_row[3] = gsampleBarcode
                    out_csv.writerow(out_row)
                assays = [a.strip() for a in row[5:] if a.strip() != '']
                if gsampPID not in data:
                    data[gsampPID] = {'plateId':gsampPID, 'destPID':gdestPID, 'custom':True, 'sampleBarcode':gsampleBarcode, 
                           'wells':[{'wellLocation':well, 'organism':{'sampleId':gsampleBarcode, 
                                    'sampleBarcode':gsampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a.strip()} \
                                    for a in assays]}}]}
                    dna_plate_samples[gdnap].append(gsampPID)
                else:
                    if well in data[gsampPID]:
                        print('duplicate well number', well, 'in plate', gsampPID, file=sys.stdout)
                        return
                    data[gsampPID]['wells'].append({'wellLocation':well, 'organism':{'sampleId':gsampleBarcode, 
                            'sampleBarcode':gsampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a.strip()} \
                            for a in assays]}})
    
        if ofh:
            ofh.close()

        if musterer_barcodes:
             musterer_mouse_info = get_musterer_mouse_info(musterer_barcodes)
             for gsampPID in data:
                 for i,well_info in enumerate(data[gsampPID]['wells']):
                     if well_info['organism']['sampleBarcode'] in musterer_mouse_info:
                         gmbc = well_info['organism']['sampleBarcode']
                         data[gsampPID]['wells'][i]['organism']['sampleId'] = musterer_mouse_info[gmbc]['mouseId']
                         data[gsampPID]['wells'][i]['organism']['strain'] = musterer_mouse_info[gmbc]['strain']
                         data[gsampPID]['wells'][i]['organism']['sex'] = musterer_mouse_info[gmbc]['sex']
        if rodentity_barcodes:
            pass
            # rodentity_mouse_info = get_rodentity_mouse_info(rodentity_barcodes)

        return data, dna_plate_samples, errs

        # below is what the old code needed returned
        #sample_pids = data.keys() # why this structure?
        #data = [(data[p],p) for p in sorted(data)]  # why?
        #return data, sample_pids, dna_plate_samples, errs
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.readCustomCSVtoJSON')


def get_assays(well_info, assay_info):
    """" 
    Combine new (first) and old (second) assay names for well. Also return the set of family names for each well.
    Largely the same as Assay.validate_assay

    Deprecated under standardised pipeline
    """
    try:
        all_assays = frozenset([n['assayName'] for n in well_info.get('assays', []) if assay_info.is_ngs(n['assayName'])] +\
                [n['assayName'] for n in well_info.get('assays', []) if assay_info.is_musterer(n['assayName'])])
        all_families = frozenset([a.split('_',1)[0] for a in all_assays if a.split('_')[0] in assay_info.all_families])
        unknown_assays = frozenset([n['assayName'] for n in well_info.get('assays', []) if n['assayName'] not in assay_info.all_assays])
        return all_assays, all_families, unknown_assays
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.get_assays')


def nimbus(dna_plate_id, platesdata, assay_info):
    """
    Process Musterer punch-plate data to produce a workfile for the BRF's Nimbus robot.
    See nimbus_custom() for custom data version.
    dna_plate_id - the barcode of the output plate
    platesdata - from JSON file
    assay_info - an Assays object from validate_assays.py
    
    It needs the name of an output file and ear-punch data for each plate.
    Output - creates Nimbus_bc.csv and Stage1-P_bc.csv files
    """
    try:
        if not file_io.is_guarded_pbc(dna_plate_id):
            dna_plate_id = file_io.guard_pbc(dna_plate_id)

        # basic validation
        dna_fn = 'Nimbus-'+dna_plate_id+'.csv'
        fnstg = 'Stage1-P'+dna_plate_id+'.csv'
        # rn - row number in Nimbus picklist file - Nimbus needs this
        # wn - count the number of wells needed in the DNA plate
        # plist - list of primer families needed for the mice (samples) in the plate
        # unk - unknown assay family names; names that don't map to a primer family name
        rn, wn, plist, unk = 0, 0, [], []
        with open(dna_fn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
            nimbus_f  = csv.writer(dstfd1, dialect='unix')
            # Note: apparently, the Nimbus is extremely particular about column names.
            nimbus_f.writerow(['Sample no', 'Plate barcode', 'Well', 'Mouse barcode'])
            stage1_f = csv.writer(dstfd2, dialect='unix')
        
            # hdr = tuple(['sampleNo', 'samplePlate','sampleWell', 'strainName', 'barcode', 'assays', 'assayFamilies', 
            #      'mouseID', 'strainName', 'sex', 'parentLitter', 'parentSequence', 'parentStrain', 'parentDate'])

            stage1_f.writerow(file_io.stage1_hdr)
            for pbc in sorted(platesdata.keys()):
                shx = platesdata.get(pbc, None)
                if not shx:
                    continue
                assert 'wells' in shx
                assert all(k in w for w in shx['wells'] for k in ('mouse', 'wellLocation'))
                wells = dict((w['wellLocation'], w['mouse']) for w in shx['wells'])
                # print("Sheet", shx.name) # use status bar
                for c in sorted(frozenset(int(w[1:]) for w in wells.keys())):
                    # rows are two blocks of 4 rows in this order - doesn't actually matter
                    for r in 'ACEGBDFH':
                        row_data = {field:'' for field in file_io.stage1_hdr}
                        nwix = r+str(c) # well index
                        nwixm = r+str(c).zfill(2) # Musterer well ID
                        # nimbus needs non-empty Mouse Barcode so use '0' for 'empty'
                        mbc = file_io.guard_mbc(str(wells[nwixm]["mouseBarcode"])) if nwixm in wells else '0'
                        rn += 1
                        rd = [rn, pbc, nwix, mbc]
                        nimbus_f.writerow(rd)
                        row_data['sampleNo'] = rn
                        row_data['samplePlate'] = pbc
                        row_data['sampleWell'] = nwixm
                        row_data['barcode'] = mbc
                        if nwixm in wells:
                            well_info = wells[nwixm]
                            mouse_assays, mouse_assayFamilies, unknown_assays = get_assays(well_info, assay_info)
                            unk += unknown_assays
                            for field in file_io.stage1_hdr:
                                if 'assay' not in field and field in well_info:
                                    row_data[field] = well_info[field]
                            if mouse_assayFamilies:
                                row_data['assays'] = ';'.join(mouse_assays)
                                row_data['assayFamilies'] = ';'.join(mouse_assayFamilies)
                                stage1_f.writerow([row_data[field] for field in file_io.stage1_hdr])
                                wn += 1
                                plist += mouse_assayFamilies
                            
        # return values for use in cgi-nimbus2.py code
        return wn, sorted(frozenset(plist)), sorted(frozenset(unk))
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.nimbus')


def nimbus_custom(dna_plate_id, platesdata, assay_info):
    """
    Process custom sample-plate data, to produce a workfile for the BRF's Nimbus robot.
    
    It needs the name of an output file and sample data for each plate.
    dna_plate_id = barcode or name of target output files, which should be guarded already        
    """
    try:
        if not file_io.is_guarded_pbc(dna_plate_id):
            dna_plate_id = file_io.guard_pbc(dna_plate_id)

        # basic validation
        dna_fn = 'Nimbus-'+dna_plate_id+'.csv'
        fnstg = 'Stage1-P'+dna_plate_id+'.csv'
        # rn - row number in Nimbus picklist file - Nimbus needs this
        # wn - count the number of wells needed in the DNA plate
        # plist - list of primer families needed for the samples in the plate
        # unk - unknown assay family names; names that don't map to a primer family name
        rn, wn, plist, unk = 0, 0, [], []
        with open(dna_fn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
            nimbus_f  = csv.writer(dstfd1, dialect='unix')
            # Note: apparently, the Nimbus is extremely particular about column names.
            nimbus_f.writerow(['Sample no', 'Plate barcode', 'Well', 'Sample barcode'])
            stage1_f = csv.writer(dstfd2, dialect='unix')
            #stghdr = ['SampleNo', 'EPplate', 'EPwell', 'sampleBarcode', 'assays', 'assayFamilies']
            stage1_f.writerow(file_io.stage1_hdr)
            for pbc in sorted(platesdata):
                shx = platesdata.get(pbc,None)
                if not shx:
                    print('No plate data found for plate', pbc, file=sys.stdout)
                    continue
                #assert 'wells' in shx
                #assert all(k in w for w in shx['wells'] for k in ('organism', 'wellLocation'))
                wells = dict((w['wellLocation'], w['organism']) for w in shx['wells'])
                
                for c in sorted(frozenset(int(w[1:]) for w in wells)):
                    # rows are two blocks of 4 rows in this order - doesn't actually matter
                    for r in 'ACEGBDFH':
                        row_data = {field:'' for field in file_io.stage1_hdr}
                        nwix = r+str(c) # well index
                        #TODO remove this restriction of wells looking like A01, etc
                        nwixm = r+str(c).zfill(2) # Sample plate well ID 
                        # nimbus needs non-empty sample barcode so use '0' for 'empty'
                        mbc = str(wells[nwixm]["sampleBarcode"]) if nwixm in wells else '0'
                        rn += 1
                        rd = [rn, pbc, nwix, mbc]
                        nimbus_f.writerow(rd)
                        row_data['sampleNo'] = rn
                        row_data['samplePlate'] = pbc
                        row_data['sampleWell'] = nwixm
                        row_data['barcode'] = mbc
                        if nwixm in wells:
                            well_info = wells[nwixm]
                            mouse_assays, mouse_assayFamilies, unknown_assays = get_assays(well_info, assay_info)
                            unk += unknown_assays
                            if 'sampleId' in well_info:
                                row_data['sampleId'] = well_info['sampleId']
                            if 'strain' in nwixm:
                                row_data['strainName'] = well_info['strain']
                            if 'sex' in nwixm:
                                row_data['sex'] = well_info['sex']
                            # in a custom run, the assays must be prescribed in the custom manifest
                            #if mouse_assayFamilies:
                                #row_data['assays'] = ';'.join(mouse_assays)
                                #row_data['assayFamilies'] = ';'.join(mouse_assayFamilies)
                            stage1_f.writerow([row_data[field] for field in file_io.stage1_hdr])  # assays are written here
                            wn += 1
                            plist += mouse_assayFamilies
        # return values for use in cgi-nimbus2.py code
        return wn, sorted(frozenset(plist)), sorted(frozenset(unk))
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.nimbus_custom')

def nimbus_gen(exp):
    """
    Standardised Nimbus input generator, produces a workfile for the BRF's Nimbus robot.
    
    It needs the name of an output file and sample data for each plate.
    dna_plate_id = barcode or name of target output files, which should be guarded already        
    """
    try:
    #if True:
        for dna_BC in exp.dest_sample_plates:
            if not file_io.is_guarded_pbc(dna_BC):
                dna_BC = file_io.guard_pbc(dna_BC)

            # basic validation
            dna_fn = os.path.join(exp.get_exp_dir(), 'Nimbus-'+dna_BC+'.csv')
            fnstg = os.path.join(exp.get_exp_dir(), 'Stage1-P'+dna_BC+'.csv')
            # rn - row number in Nimbus picklist file - Nimbus needs this
            # wn - count the number of wells needed in the DNA plate
            # plist - list of primer families needed for the samples in the plate
            # unk - unknown assay family names; names that don't map to a primer family name
            row_number, wells_used = 0, 0
            with open(dna_fn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
                nimbus_f  = csv.writer(dstfd1, dialect='unix')
                # Note: apparently, the Nimbus is extremely particular about column names.
                nimbus_hdr = ['Sample no', 'Plate barcode', 'Well', 'Sample barcode']
                nimbus_f.writerow(nimbus_hdr)
                stage1_f = csv.writer(dstfd2, dialect='unix')
                stage1_hdr = tuple(['sampleNumber', 'samplePlate','sampleWell', 'sampleBarcode', 'assays', 'assayFamilies', 
                        'strain', 'sex'])
            
                stage1_f.writerow(stage1_hdr)
                for pbc in sorted(exp.dest_sample_plates[dna_BC]):
                    shx = exp.plate_location_sample.get(pbc,None)
                    if not shx:
                        print('No plate data found for plate', pbc, file=sys.stdout)
                        continue
                    # Get four wells in the same column at a time
                    for pos1,pos2,pos3,pos4 in zip(nimbus_ordered_96[0::4], nimbus_ordered_96[1::4], nimbus_ordered_96[2::4], nimbus_ordered_96[3::4]):
                        if pos1 not in shx and pos2 not in shx and pos3 not in shx and pos4 not in shx:
                            continue  # skip missing blocks of 4 contiguous rows in the same column
                        for pos in (pos1, pos2, pos3, pos4):
                            row_number += 1
                            row_data = {field:'' for field in file_io.stage1_hdr}
                            row_data['Sample no'] = row_number
                            if pos in shx:
                                row_data['Sample barcode'] = shx[pos]['barcode']
                            else:
                                row_data['Sample barcode'] = '0'
                            row_data['Well'] = pos
                            row_data['Plate barcode'] = pbc
                            nimbus_f.writerow([row_data[field] for field in nimbus_hdr])

                            if row_data['Sample barcode'] == '0':
                                continue
                    
                            # stage files are no longer necessary, but provide useful debugging for now
                            row_data = {field:'' for field in stage1_hdr}
                            if 'sampleNumber' in shx[pos]:
                                row_data['sampleNumber'] = shx[pos]['sampleNumber']
                            else:
                                row_data['sampleNumber'] = 0
                            row_data['samplePlate'] = pbc
                            row_data['sampleWell'] = pos
                            row_data['sampleBarcode'] = shx[pos]['barcode']
                            row_data['assays'] = shx[pos]['assays']
                            row_data['assayFamilies'] = shx[pos]['assayFamilies']
                            if 'strain' in shx[pos]:
                                row_data['strain'] = shx[pos]['strain']
                            else:
                                row_data['strain'] = ''
                            if 'sex' in shx[pos]:
                                row_data['sex'] = shx[pos]['sex']
                            else:
                                row_data['sex'] = ''
                            stage1_f.writerow([row_data[field] for field in stage1_hdr])  # assays are written here
                            wells_used += 1
                
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.nimbus_custom')
        return False
    return True

if __name__ == '__main__':
    pass