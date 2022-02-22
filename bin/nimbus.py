#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Oct 2019
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.15
@version_comment: Uses standardised custom/mouse stage files. Error handling updated.
@last_edit: 2022-02-16
@edit_comment: 

GUI for Nimbus files for the NGSgenotyping project.
Note - The GUI code is no longer used - the nimbus() function is called by cgi-nimbus2.py
in the new HTML interface.

Input is 1-4 96-well sample (ear-punch) plate bar codes ... and a destination 384-well plate bar code.

The program uses either a Musterer webservice for each ear-punch plate, or has 
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
import json
import tkinter as tk
#import tkinter.ttk as ttk
import tkinter.filedialog as fd
import tkinter.messagebox as mb  
import requests
from musterer import get_plate
import file_io
from util import output_error

#import xlrd
#import xlwt

app = None # global!

def readCustomCSVtoJSON(input_fn, custom_file_contents='custom'):
    """ read a custom sample manifest with no more than 4x96 well plates! Then return as JSON 
        Creates a copy of the original manifest with guarded barcodes.
        v0.15 attempts to download mouse data for valid mouse barcodes in manifest
    """
    try:
        data = {}
        errs = []  # error messages
        guarded_fn = input_fn.replace('.csv','_guarded.csv')
        ofh = None
        musterer_barcodes = []
        rodentity_barcodes = []
        if not os.path.exists(guarded_fn):
            ofh = open(guarded_fn, 'wt')
            out_csv = csv.writer(ofh, dialect='unix')
        with open(input_fn, 'rt', newline=None) as f:
            in_csv = csv.reader(f, delimiter=',')
            for i, row in enumerate(in_csv):
                #print(row, file=sys.stderr)
                #print(file=sys.stderr)
                if i == 0:
                    if ofh:
                        out_csv.writerow(row)
                    continue  # header
                if all(col.strip() == '' for col in row):
                    continue  # blank rows
                if any(col.strip() == '' and j < 5 for j,col in enumerate(row)):
                    msg = f"ERROR. blank manifest entry in entry {i}: {row}"
                    print(msg, file=sys.stderr)
                    errs.append(msg)
                    continue
                    # return
                plt = str(row[1].strip())
                if is_guarded_pbc(plt):
                    gplate = plt
                else:
                    gplate = guard_pbc(plt)
                well = str(row[2]).strip()
                smpl = str(row[3].strip())
                if is_guarded(smpl):
                    gsampleBarcode = smpl
                elif custom_file_contents == 'custom':
                    gsampleBarcode = guard_cbc(smpl)
                elif custom_file_contents == 'musterer':
                    gsampleBarcode = guard_mbc(smpl)
                    musterer_barcodes.append(gsampleBarcode)
                elif custom_file_contents == 'rodentity':
                    gsampleBarcode = guard_rbc(smpl)
                    rodentity_barcodes.append(gsampleBarcode)
                else:
                    msg = f"ERROR. unknown custom file contents field: {custom_file_contents}"
                    print(mgs, file=sys.stderr)
                    errs.append(msg)

                if ofh:
                    out_row = [r for r in row]
                    out_row[1] = gplate
                    out_row[3] = gsampleBarcode
                    out_csv.writerow(out_row)
                assays = [a.strip() for a in row[4:] if a.strip() != '']
                if gplate not in data:
                    if len(data) == 4:
                        errs.append('Too many input plates specified in experiment file. Max of 4 plates are allowed per file!')
                        break
                    data[gplate] = {'plateId':gplate, 'custom':True, 'sampleBarcode':gsampleBarcode, 
                           'wells':[{'wellLocation':well, 'organism':{'sampleId':gsampleBarcode, 
                                    'sampleBarcode':gsampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a.strip()} \
                                    for a in assays]}}]}                  
                else:
                    if well in data[gplate]:
                        print('duplicate well number', well, 'in plate', gplate, file=sys.stderr)
                        return
                    data[gplate]['wells'].append({'wellLocation':well, 'organism':{'sampleId':gsampleBarcode, 
                            'sampleBarcode':gsampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a.strip()} \
                            for a in assays]}})
    
        if ofh:
            ofh.close()

        if musterer_barcodes:
             musterer_mouse_info = get_musterer_mouse_info(musterer_barcodes)
             for gplate in data:
                 for i,well_info in enumerate(data[gplate]['wells']):
                     if well_info['organism']['sampleBarcode'] in musterer_mouse_info:
                         gmbc = well_info['organism']['sampleBarcode']
                         data[gplate]['wells'][i]['organism']['sampleId'] = musterer_mouse_info[gmbc]['mouseId']
                         data[gplate]['wells'][i]['organism']['strain'] = musterer_mouse_info[gmbc]['strain']
                         data[gplate]['wells'][i]['organism']['sex'] = musterer_mouse_info[gmbc]['sex']
        if rodentity_barcodes:
            pass
            # rodentity_mouse_info = get_rodentity_mouse_info(rodentity_barcodes)

        pids = data.keys()
        data = [(data[p],p) for p in sorted(data)]
        return data, pids, errs
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.readCustomCSVtoJSON')


def get_assays(well_info, assay_info):
    """" 
    Combine new (first) and old (second) assay names for well. Also return the set of family names for each well.
    Largely the same as Assay.validate_assay()
    """
    try:
        all_assays = frozenset([n['assayName'] for n in well_info.get('assays', []) if assay_info.is_ngs(n['assayName'])] +\
                [n['assayName'] for n in well_info.get('assays', []) if assay_info.is_musterer(n['assayName'])])
        all_families = frozenset([a.split('_',1)[0] for a in all_assays if a.split('_')[0] in assay_info.all_families])
        unknown_assays = frozenset([n['assayName'] for n in well_info.get('assays', []) if n['assayName'] not in assay_info.all_assays])
        return all_assays, all_families, unknown_assays
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.get_assays')


def nimbus(tgtid, platesdata, assay_info):
    """
    Process Musterer punch-plate data to produce a workfile for the BRF's Nimbus robot.
    See nimbus_custom() for custom data version.
    platesdata - from JSON file
    assay_info - an Assays object from validate_assays.py
    
    It needs the name of an output file and ear-punch data for each plate.
    Output - creates Nimbus_bc.csv and Stage1-P_bc.csv files
    """
    try:
        if not file_io.is_guarded_pbc(tgtid):
            tgtid = file_io.guard_pbc(tgtid)

        # basic validation
        tgtfn = 'Nimbus-'+tgtid+'.csv'
        fnstg = 'Stage1-P'+tgtid+'.csv'
        # rn - row number in Nimbus picklist file - Nimbus needs this
        # wn - count the number of wells needed in the DNA plate
        # plist - list of primer families needed for the mice (samples) in the plate
        # unk - unknown assay family names; names that don't map to a primer family name
        rn, wn, plist, unk = 0, 0, [], []
        with open(tgtfn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
            nimbus_f  = csv.writer(dstfd1, dialect='unix')
            # Note: apparently, the Nimbus is extremely particular about column names.
            nimbus_f.writerow(['Sample no', 'Plate barcode', 'Well', 'Mouse barcode'])
            stage1_f = csv.writer(dstfd2, dialect='unix')
        
            # hdr = tuple(['sampleNo', 'samplePlate','sampleWell', 'strainName', 'barcode', 'assays', 'assayFamilies', 
            #      'mouseID', 'strainName', 'sex', 'parentLitter', 'parentSequence', 'parentStrain', 'parentDate'])

            stage1_f.writerow(file_io.stage1_hdr)
            for no, pbc, shx in platesdata: # plate no, barcode, platedata (from JSON format)
                # print('Plate', pbc, 'data:', ', '.join(k+': '+str(v) for k,v in shx.items() if k!=' ))
                # Nope! assert str(shx['barcode']).endswith(pbc)
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


def nimbus_custom(tgtid, platesdata, assay_info):
    """
    Process custom sample-plate data, to produce a workfile for the BRF's Nimbus robot.
    
    It needs the name of an output file and sample data for each plate.
    tgtid = barcode or name of target output files, which should be guarded already        
    """
    try:
        if not file_io.is_guarded_pbc(tgtid):
            tgtid = file_io.guard_pbc(tgtid)

        # basic validation
        tgtfn = 'Nimbus-'+tgtid+'.csv'
        fnstg = 'Stage1-P'+tgtid+'.csv'
        # rn - row number in Nimbus picklist file - Nimbus needs this
        # wn - count the number of wells needed in the DNA plate
        # plist - list of primer families needed for the samples in the plate
        # unk - unknown assay family names; names that don't map to a primer family name
        rn, wn, plist, unk = 0, 0, [], []
        with open(tgtfn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
            nimbus_f  = csv.writer(dstfd1, dialect='unix')
            # Note: apparently, the Nimbus is extremely particular about column names.
            nimbus_f.writerow(['Sample no', 'Plate barcode', 'Well', 'Sample barcode'])
            stage1_f = csv.writer(dstfd2, dialect='unix')
            #stghdr = ['SampleNo', 'EPplate', 'EPwell', 'sampleBarcode', 'assays', 'assayFamilies']
            stage1_f.writerow(file_io.stage1_hdr)
            for no, pbc, shx in platesdata: # plate no, barcode, platedata (from JSON format)
                #print(f"Plate {no}, {pbc}, data: {shx}", file=sys.stderr)
                # Nope! assert str(shx['barcode']).endswith(pbc)
                if not shx:
                    print('No plate data found for plate', no, pbc, file=sys.stderr)
                    continue
                assert 'wells' in shx
                assert all(k in w for w in shx['wells'] for k in ('organism', 'wellLocation'))
                wells = dict((w['wellLocation'], w['organism']) for w in shx['wells'])
                #print(wells.keys(), file=sys.stderr)
                # print("Sheet", shx.name) # use status bar
                for c in sorted(frozenset(int(w[1:]) for w in wells.keys())):
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
                            stage1_f.writerow([row_data[field] for field in file_io.stage1_hdr])
                            wn += 1
                            plist += mouse_assayFamilies
        # return values for use in cgi-nimbus2.py code
        return wn, sorted(frozenset(plist)), sorted(frozenset(unk))
    except Exception as exc:
        output_error(exc, msg='Error in nimbus.nimbus_custom')


if __name__ == '__main__':
    pass