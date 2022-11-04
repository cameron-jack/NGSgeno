#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from experiment import Experiment

try:
    import bin.util as util
except ModuleNotFoundError:
    import util

import sys
import os
import gzip
import csv
from pathlib import Path
import traceback
from collections import Counter, OrderedDict
import itertools
from copy import deepcopy

"""
@created: Nov 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: All machine control file reading/writing functions are now here
@last_edit: 2022-05-05
@edit_comment: relocated guards to util.py

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

The functions below should be organised in the order they are first seen in the operation of the pipeline
"""

### Machine IO functions

def nimbus_gen(exp):
    """
    Standardised Nimbus input generator, produces a workfile for the BRF's Nimbus robot.
    
    It needs the name of an output file and sample data for each plate.
    dna_plate_id = barcode or name of target output files, which should be guarded already
    As a file generator() it needs to create a transaction (exp.reproducible_steps) entry
    """
    transactions = {}
    try:
    #if True:
        
        for dna_BC in exp.dest_sample_plates:
            if not util.is_guarded_pbc(dna_BC):
                dna_BC = util.guard_pbc(dna_BC)

            dna_fn = exp.get_exp_fp('Nimbus-'+dna_BC+'.csv', transaction=True)
            fnstg = exp.get_exp_fp('Stage1-P'+dna_BC+'.csv', transaction=True)

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
                stage1_hdr = tuple(['sampleNumber', 'samplePlate','sampleWell', 'sampleBarcode', 'assays', 'assayFamilies', 
                        'strain', 'sex'])
            
                stage1_f.writerow(stage1_hdr)
                for pbc in sorted(exp.dest_sample_plates[dna_BC]):
                    shx = exp.plate_location_sample.get(pbc,None)
                    if not shx:
                        print('No plate data found for plate', pbc, file=sys.stdout)
                        continue
                    transactions[dna_fn][pbc] = []
                    # Get four wells in the same column at a time
                    for pos1,pos2,pos3,pos4 in zip(util.nimbus_ordered_96[0::4], util.nimbus_ordered_96[1::4], 
                            util.nimbus_ordered_96[2::4], util.nimbus_ordered_96[3::4]):
                        if pos1 not in shx and pos2 not in shx and pos3 not in shx and pos4 not in shx:
                            continue  # skip missing blocks of 4 contiguous rows in the same column
                        for pos in (pos1, pos2, pos3, pos4):
                            row_number += 1
                            row_data = {field:'' for field in stage1_hdr}
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
                    
                            # stage files are still used, but only for historical reasons
                            row_data = {field:'' for field in stage1_hdr}
                            if 'sampleNumber' in shx[pos]:
                                row_data['sampleNumber'] = shx[pos]['sampleNumber']
                            else:
                                row_data['sampleNumber'] = 0
                            row_data['samplePlate'] = pbc
                            row_data['sampleWell'] = pos
                            row_data['sampleBarcode'] = shx[pos]['barcode']
                            row_data['assays'] = ';'.join(shx[pos]['assays'])
                            row_data['assayFamilies'] = ';'.join(shx[pos]['assayFamilies'])
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
                            transactions[dna_fn][pbc].append((pos,-1000)) # 1000 nl of sample is transferred
        
                
    except Exception as exc:
        print("Transactions in nimbus_gen on fail: ", transactions, file=sys.stderr)
        exp.log(f'Error: Nimbus input file creation {exc}')
        exp.save()
        return False
    print("Transactions in nimbus_gen()", transactions, file=sys.stderr)
    exp.add_pending_transactions(transactions)
    exp.log(f'Success: Hamilton Nimbus plate definition files have been generated')
    exp.save()
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
    

def grouper(xs, kf=lambda x:x[0]):
    "group pairs: (x,a), (x,b), ... => (x,[a,b, ...]), ..."
    return itertools.groupby(sorted(xs, key=kf), key=kf)

#def file_get_check(exp, fids, fmt):
#    """ collect Echo DNA plate file names, match plate ids, remove duplicates """
#    try:
#        dups = [fids[i] for i in range(1, len(fids)) if fids[i] in fids[:i]]
#        fx = [fids[i] for i in range(len(fids)) if fids[i] not in fids[:i]]
#        if dups:
#            exp.log(f"Warning: Duplicate DNA plate IDs ignored: {dups}")    
    
#        globs = dict((pid, sorted(glob.glob(fmt.format(pid)))) for pid in fx)
#        nofile = [fid for fid, fns in globs.items() if not fns]
#        if nofile:
#            exp.log(f"Error: files not found {', '.join([fmt.format(fid) for fid in fids])}")
#            return {}
#        return collections.OrderedDict((fid, sorted(ps)[-1]) for fid, ps in globs.items())
#    except Exception as exc:
#        exp.log(f"Error: {exc}")
#        return {}


def mk_picklist(exp, fn, rows, transactions, output_plate_guards=False):
    """ 
    Output an Echo picklist given the header and rows (transfer spec) 
    Inputs:
        exp - an Experiment instance
        fn - output filename
        rows - data rows to write
        transactions - a dictionary of changes from the input plates
        output_plate_guards - whether or not to write guard characters on the plate barcodes
    """
    try:
        plhdr_str = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"+\
                "Destination Plate Name,Destination Plate Barcode,Destination Plate Type,"+\
                "Destination Well,Volume"
        plhdr = plhdr_str.split(',')
       
        with open(fn, "wt", newline='') as dstfd:
            dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            dst.writerow(plhdr)
            data = []
            for row in rows:
                if len(row) != len(plhdr):
                    exp.log(f"Warning: entry does not match expected header format: {row=} {plhdr=}")
                    continue
                src_plate_barcode = row[1]
                src_plate_well = row[3]
                vol = row[-1]
                if src_plate_barcode not in transactions:
                    transactions[src_plate_barcode] = {}
                if src_plate_well not in transactions[src_plate_barcode]:
                    transactions[src_plate_barcode][src_plate_well] = 0
                transactions[src_plate_barcode][src_plate_well] -= vol

                if output_plate_guards:  # guard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() else util.guard_pbc(str(d), silent=True) for (h,d) in zip(plhdr, row)])
                else:  # unguard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() else util.unguard_pbc(str(d), silent=True) for (h,d) in zip(plhdr, row)])
            dst.writerows(data)
    except Exception as exc:
        exp.log(f"Failure: Picklist generation failed for {fn} {exc}")
        print(f"Failure: Picklist generation failed for {fn} {exc}", file=sys.stderr)
        if Path(fn).exists():
            os.remove(fn)
        transactions = None
        exp.save()
        return False
    exp.add_pending_transactions(transactions)
    exp.log(f'Success: Wrote picklist to {fn}')
    exp.save()
    return True


def mk_mytaq_picklist(exp, fn, task_wells, taqwater_bcs, taq_vol, water_vol, transactions=None):
    """ 
    Create a picklist for Mytaq & H2O transfers
    transactions are passed through to mk_picklist and modified there
    """        
    # build a list of all possible taq/water pids and well pair transfers
    pid_ww_tw_list = []
    for pid in taqwater_bcs:
        twp = exp.get_plate(pid, transactions)
        # iterate over wells forever
        ww_gen = (x for xs in itertools.repeat(twp['water_wells']) for x in xs)
        tw_gen = (x for xs in itertools.repeat(twp['taq_wells']) for x in xs)
        empty_well_pairs = 0
        for ww,tw in zip(next(ww_gen),next(tw_gen)):
            if empty_well_pairs >= 3:
                break  # time for a new plate
            if twp[ww] - water_vol - twp['dead_vol'] < 0:
                empty_well_pairs += 1
                continue
            if twp[tw] - taq_vol - twp['dead_vol'] < 0:
                empty_well_pairs += 1
                continue
            twp[ww] -= water_vol
            twp[tw] -= taq_vol
            pid_ww_tw_list.append((pid,ww,tw))

    if len(pid_ww_tw_list) < len(task_wells):
        exp.log(f'Failure: Not enough taq or water available from plates {taqwater_bcs}')
        exp.save()
        return False

    # now allocate to used wells
    try:
        pcr_pids = sorted(set([w.pcrplate for w in task_wells]))
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
                    dst_dict[task_well.pcrplate], util.unguard(task_well.pcrplate, silent=True), 
                    util.PLATE_TYPES['PCR384'], task_well.pcrwell, water_vol]
            output_rows.append(water_row)
            taq_row = [src_dict[tw_pid], util.unguard(tw_pid, silent=True), util.PLATE_TYPES['Echo6'], tw,
                    dst_dict[task_well.pcrplate], util.unguard(task_well.pcrplate, silent=True), 
                    util.PLATE_TYPES['PCR384'], task_well.pcrwell, taq_vol]
            output_rows.append(taq_row)
        mk_picklist(exp, fn, output_rows, transactions)         
    except Exception as exc:
        exp.log(f"Error: failed to generate taq/water picklist {exc}")
        exp.save()
        return False
    return True


def generate_echo_PCR1_picklist(exp, dna_plate_bcs, pcr_plate_bcs, taq_water_bcs):
    """
    Entry point. Takes an experiment instance plus plate barcodes for dna plates, PCR plates, primer plates, taq/water plates.
    """
    try:
        pcr_bcs = []
        dna_bcs = []
        taq_bcs = []
        try:
            pcr_bcs = [util.guard_pbc(p,silent=True) for p in pcr_plate_bcs]
        except Exception as e:
            exp.log(f"{e}")
        try:
            dna_bcs = [util.guard_pbc(d,silent=True) for d in dna_plate_bcs] 
        except Exception as e:
            exp.log(f"{e}")
        try:                          
            taq_bcs = [util.guard_pbc(t, silent=True) for t in taq_water_bcs]
        except Exception as e:
            exp.log(f"{e}")

        primer_survey = exp.get_exp_fp('primer-svy.csv')

        # read Nimbus created DNA plates into experiment and copy source well information
        transactions = {}
        nimbus_outfiles = []
        for d in sorted(dna_bcs):
            fp = exp.get_exp_fp("Echo_384_COC_0001_" + d + "_0.csv")
            if not Path(fp).is_file():
                exp.log(f"{d} has no matching Echo_384_COC file from the Nimbus", level='error')
                continue
            
            nimbus_outfiles.append(fp)
            with open(fp, 'rt') as f:
                #RecordId	TRackBC	TLabwareId	TPositionId	SRackBC	SLabwareId	SPositionId
                #1	p2021120604p	Echo_384_COC_0001	A1	p2111267p	ABg_96_PCR_NoSkirt_0001	A1
                exp.plate_location_sample[d] = {'purpose':'dna','wells':set(),'source':'Nimbus','plate_type':'384PP_AQ_BP'}
                for i, line in enumerate(f):
                    if i == 0:  # header
                        continue
                    cols = [c.strip().strip('\"') for c in line.split(',')]
                    source_pos = util.unpadwell(cols[-1])
                    source_plate = cols[-3]
                    dest_plate = cols[1]
                    dest_pos = util.unpadwell(cols[3])
                    if dest_plate != d:
                        exp.log(f"{util.unguard_pbc(d, silent=True)} doesn't match {util.unguard_pbc(dest_plate, silent=True)} "+\
                            f"as declared in Echo_384_COC file: {fp}", level="error")
                    exp.plate_location_sample[d]['wells'].add(dest_pos)
                    try:
                        exp.plate_location_sample[d][dest_pos] = deepcopy(exp.plate_location_sample[source_plate][source_pos])
                    except:
                        exp.log(f"Critical: cannot locate {d=} {dest_pos=} {source_plate=} {source_pos=}")

        outfmt = f"PCR-picklist_{exp.name}.csv"

        #dna_plate_paths = {bc: exp.get_exp_fp(f"Echo_384_COC_0001_{bc}_0.csv") for bc in sorted(dna_bcs) if \
        #        Path(exp.get_exp_fp(f"Echo_384_COC_0001_{bc}_0.csv")).exists()}
        #dnafns = [fn for fn in dna_plate_paths if Path(fn).exists()]
        #dnafns = file_get_check(exp, dna_bcs, exp.get_exp_fp("Echo_384_COC_0001_{0}_?.csv"))
        nimcolmap = (("TRackBC", "dstplate"), ("TPositionId", 'dstwell'), ('SRackBC', 'srcplate'), ('SPositionId', 'srcwell'))
        typenim = util.Table.csvtype(nimbus_outfiles[0], 'NimRec', hdrmap=nimcolmap)
        nimbusTables = [util.CSVTable(typenim, fn) for fn in nimbus_outfiles]
        
        # read Stage1 files
        dnadict = dict(((x.srcplate, x.srcwell), (x.dstplate, x.dstwell)) for nt in nimbusTables for x in nt.data)    
        dnas = [util.CSVTable('S1Rec', exp.get_exp_fp('Stage1-P{}.csv'.format(dnabc))) for dnabc in dna_bcs]
        
        primerTable = util.CSVTable("PPRec", primer_survey, fields="spn spbc spt well primer volume".split(' '))
        primset = sorted(frozenset(x.primer for x in primerTable.data if x.primer != ''))
        pfdict = dict((k.strip(), list([p.strip() for p in g])) for k, g in itertools.groupby(primset, key=lambda x:x.split('_',1)[0]))

        # Add record for PCR wells into Stage2 files
        wgenflds =  dnas[0].tt._fields + ('dnaplate', 'dnawell') + ('primer',)
        # Is this the plating of DNA samples for each assay? Yes. It is.
        wgen = [xs+dnadict[(xs.samplePlate, xs.sampleWell)]+(x,) for xss in dnas for xs in xss.data \
                for f in xs.assayFamilies.split(';') if f in pfdict for x in pfdict[f]]
        
        # allocate PCR plate wells
        wells = [r+str(c+1) for c in range(24) for r in"ABCDEFGHIJKLMNOP"]
        pcrwellgen = ((p, w) for p in pcr_bcs for w in wells)
        
        # allocate PCR well for each sample
        s2flds = wgenflds+('pcrplate', 'pcrwell')
        S2Rec = util.Table.newtype('S2Rec', s2flds)
        s2tab = util.Table(S2Rec, ([x for xs in rx for x in xs] for rx in zip(wgen, pcrwellgen)), headers=s2flds) 
        # output Stage 2 CSV file - used in Stage 3 below - keeping the plate guards for the next stage (to be safe)
        stage2_fn = exp.get_exp_fp("Stage2.csv", transaction=True)
        try:
            s2tab.csvwrite(stage2_fn, output_plate_guards=True)
        except Exception as exc:
            exp.log(f'Critical: failed to write {stage2_fn} {exc}')
            if Path(stage2_fn).exists():
                os.remove(stage2_fn)
            transactions[stage2_fn] = None
            return False

        # PCR1 picklists: DNA/sample, Primer and Taq/H20
    
        dst_plate_type = util.PLATE_TYPES['PCR384']
        src_plate_type = util.PLATE_TYPES['Echo384']
        # sample DNA picklist
        ddict = dict((pbc, f"Destination[{i}]") for i, pbc in enumerate(pcr_bcs, start=1))
        sdict = dict((pbc, f"Source[{i}]")for i, pbc in enumerate(dna_bcs, start=1))
        volume = exp.transfer_volumes['DNA_VOL']
        dna_fn = exp.get_exp_fp(outfmt.replace('PCR','PCR1_dna'), transaction=True)
        transactions[dna_fn] = {}
        gen = ((sdict[r.dnaplate], r.dnaplate, src_plate_type, r.dnawell,
                ddict[r.pcrplate], r.pcrplate, dst_plate_type, r.pcrwell, volume)
                for r in s2tab.data)
        success = mk_picklist(exp, dna_fn, gen, transactions[dna_fn])
        if not success:
            transactions[dna_fn] = None
        # Primer Picklist
        # [Source[1]', '', '384PP_AQ_BP', 'H6', 'Destination[1]', '3121', 'Hard Shell 384 well PCR Biorad', 'A1', 500]
        primsrc = PicklistSrc(exp.get_exp_fp("primer-svy.csv"), idx=4) # same name as in cgi-nimbus2.py
        volume = exp.transfer_volumes['PRIMER_VOL']        
        primer_fn = exp.get_exp_fp(outfmt.replace('PCR','PCR1_primer'), transaction=True)
        transactions[primer_fn] = {}
        primer_uses = Counter()
        primer_output_rows = []
        for r in s2tab.data:
            if r.primer not in primsrc.data:
                exp.log(f'Warning: Cannot find {r.primer} in known primers')
                continue
            # cycle through available primer wells
            primer_index = primer_uses[r.primer] % len(primsrc.data[r.primer])
            primer_output_rows.append(primsrc.data[r.primer][primer_index][:4] +
                    [ddict[r.pcrplate], r.pcrplate, util.PLATE_TYPES['PCR384'], r.pcrwell, volume])
            primer_uses[r.primer] += 1
            
        success = mk_picklist(exp, primer_fn, primer_output_rows, transactions[primer_fn])
        if not success:
            transactions[primer_fn] = None

        # also PCR1 water and Taq
        taq_fn = exp.get_exp_fp(outfmt.replace('PCR','PCR1_taqwater'), transaction=True)
        transactions[taq_fn] = {}
        print(f"end of generate_pcr1 {taq_bcs=}")
        # exp, fn, task_wells, taqwater_bcs, taq_vol, water_vol, transactions=None
        success = mk_mytaq_picklist(exp, taq_fn, s2tab.data, taq_bcs, exp.transfer_volumes['PRIMER_TAQ_VOL'], 
                exp.transfer_volumes['PRIMER_WATER_VOL'], transactions[taq_fn])
        if not success:
            transactions[taq_fn] = None

    except Exception as exc:
        exp.log(f"Failure: PCR1 Echo picklists could not be created {exc}")
        return False
    exp.add_pending_transactions(transactions)
    exp.log(f'Success: PCR1 Echo picklists created')
    exp.save()
    return True


def i7i5alloc_rot(exp, vol, wellcount, index_survey_fn='index-svy.csv'):
    """ allocate vol to wellcount wells with unique index pairs - rotates barcode wells to avoid wells being drained """
    try:
        typebc = util.Table.newtype('BCRecord', "name platebc type well set barcode orderpart oligo volume")
        # typebc = Table.newtype('BCRecord', "well name index indexName oligo volume")
        tab = util.CSVTable(typebc, exp.get_exp_fp(index_survey_fn)) # volumes in file are in uL

        dv = exp.DEAD_VOLS[util.PLATE_TYPES['Echo384']]
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
            exp.log(f"Error: Not enough barcode combos {combos} for wells requested {wellcount}")
            return
        while len(index_info_pairs) < wellcount:     
            if len(i7s_empty) == len(i7s):
                exp.log(f"Error: i7 wells contain insufficient volume for the number of experiment wells")
                exp.log(f"Info: {i7s_empty=}, {i5s_empty=}, {len(used_index_pairs)=}")
                return
            elif len(i5s_empty) == len(i5s):
                exp.log(f"Error: i5 wells contain insufficient volume for the number of experimental wells")
                exp.log(f"Info: {i7s_empty=}, {i5s_empty=}, {len(used_index_pairs)=}")
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
      
    except Exception as exc:
        exp.log(f"Error: {exc}")


def generate_miseq_samplesheet(exp, fmtmiseq, s3tab,transactions=None):
    """ Now that we have all the required info, generate the Illumina Miseq samplesheet, which drives the sequencing """
    if transactions is None:
        transactions = {}
    try:
        fnmiseq = exp.get_exp_fp(fmtmiseq.format(exp.name), transaction=True)
        with open(fnmiseq, "wt", newline='') as dstfd:
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
            gen = [(r.pcrplate+'_'+util.padwell(r.pcrwell), r.sampleBarcode, '', '',
                    r.i7name, r.i7bc,r.i5name, r.i5bc,'NGSgeno') for r in s3tab.data]
            
            dst.writerows(gen)
        transactions[fnmiseq] = []
        
    except Exception as exc:
        exp.log(f"Failure: Miseq samplesheet generation {fnmiseq} will be deleted {exc}")
        if Path(fnmiseq).exist():
            os.remove(fnmiseq)
        exp.save()
        return False
    exp.add_pending_transactions(transactions)
    exp.log(f'Success: Miseq samplesheet written to {fnmiseq}')
    exp.save()
    return True

    
def generate_echo_PCR2_picklist(exp, pcr_plate_bcs, index_plate_bcs, taq_water_bcs, amplicon_plate_bcs=None):
    """
    Replacement for main(). Called as a library only.
    User included PCR plates, index plates, taq/water plates, amplicon plates (optional) are user provided
    Also calls generate_miseq_samplesheet()
    Called from Experiment.generate_echo_PCR2_picklists()
    """
    transactions = {}
    try:
        pcr_bcs = []
        amplicon_bcs = []
        taq_bcs = []
        try:
            pcr_bcs = [util.guard_pbc(p,silent=True) for p in pcr_plate_bcs]
        except Exception as exc:
            exp.log(f"Error: PCR plate barcode in error {pcr_plate_bcs=} {exc}")
        try:
            index_bcs = [util.guard_pbc(i, silent=True) for i in index_plate_bcs]
        except Exception as exc:
            exp.log(f"Error: Index plate parcode in error {index_plate_bcs=} {exc}")
        try:
            amplicon_bcs = [util.guard_pbc(d,silent=True) for d in amplicon_plate_bcs] 
        except Exception as exc:
            exp.log(f"Error: Amplicon plate barcodes in error {amplicon_plate_bcs=} {exc}")
        try:                          
            taq_bcs = [util.guard_pbc(t, silent=True) for t in taq_water_bcs]
        except Exception as exc:
            exp.log(f"Error: Taq-water plate barcodes in error {taq_water_bcs=} {exc}")
        outfmt = f"PCR-picklist_{exp.name}.csv"
    
        fnstage2, fnstage3 = "Stage2.csv", "Stage3.csv"
    
        ### i7i5 barcodes
        # read Stage2 csv file
        #sampleNumber	samplePlate	sampleWell	sampleBarcode	assays	assayFamilies	strain	sex	dnaplate	dnawell	primer	pcrplate	pcrwell
        #typebc = Table.newtype('S2Record',
        s2tab = util.CSVTable('S2Rec', exp.get_exp_fp(fnstage2))           
            
        index_vol = exp.transfer_volumes['INDEX_VOL']  # 175 nanolitres
        # count up all amplicon sample wells
        total_wells = len(s2tab.data)
        for amp_pid in amplicon_bcs:
            total_wells += len(exp.plate_location_sample[amp_pid]['wells'])

        index_alloc = i7i5alloc_rot(exp, index_vol, total_wells)
        #print(index_alloc, file=sys.stdout)
            
        s3flds = s2tab.tt._fields+('i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well', 'index_plate')
        S3Rec = util.Table.newtype('S3Rec', s3flds)
        
        s3tab = util.Table(S3Rec, ([x for xs in (p1, p2[:3], p3[:4]) for x in xs] for p1, (p2, p3) in zip(s2tab.data, index_alloc)), headers=s3flds) 
        # output Stage 3 CSV file - used for custom and mouse samples. Amplicons are handled separately
        transactions[fnstage3] = {}
        try:
            s3tab.csvwrite(exp.get_exp_fp(fnstage3, transaction=True))
        except Exception as exc:
            exp.log(f'Critical: could not write {fnstage3}')
            transactions[fnstage3] = None
            exp.save()
            return False

        # write index picklists
        # make source dictionary for index plates
        src_dict = {}
        for r in s3tab.data:
            if r.index_plate not in src_dict:
                src_dict[r.index_plate] = f"Source[{len(src_dict)+1}]"
        # make the destination dictionary for pcrplates
        dest_dict = {}
        for r in s3tab.data:
            if r.pcrplate not in dest_dict:
                dest_dict[r.pcrplate] = f"Destination[{len(dest_dict)+1}]"

        pt_src = util.PLATE_TYPES['Echo384']
        pt_dst = util.PLATE_TYPES['PCR384']
        rowi7 = ((src_dict[r.index_plate], r.index_plate, pt_src, r.i7well, dest_dict[r.pcrplate], r.pcrplate, pt_dst, r.pcrwell, index_vol) for r in s3tab.data)
        rowi5 = ((src_dict[r.index_plate], r.index_plate, pt_src, r.i5well, dest_dict[r.pcrplate], r.pcrplate, pt_dst, r.pcrwell, index_vol) for r in s3tab.data)
        fn_index = exp.get_exp_fp(outfmt.replace('PCR','PCR2_index'), transaction=True)
        mk_picklist(exp, fn_index, (r for rs in (rowi7, rowi5) for r in rs), transactions[fn_index])    
        
        # also PCR2 water and Taq
        fn_taqwater = exp.get_exp_fp(outfmt.replace('PCR','PCR2_taqwater'), transaction=True)
        mk_mytaq_picklist(exp, fn_taqwater, s3tab.data, taq_bcs, exp.transfer_volumes['INDEX_TAQ_VOL'], 
                exp.transfer_volumes['INDEX_WATER_VOL'], transactions)
        
        # MiSeq file - experiment name, name format, ngid(dir), s3 data, verbosity
        fmtmiseq = "MiSeq-{}.csv"
        success = generate_miseq_samplesheet(exp, fmtmiseq, s3tab, transactions)
        
        fns3 = "Stage3.html"
        #s2r.build_report(exp, fns3, s3tab.tt._fields, s3tab.data)    
        
    except Exception as exc:
        exp.log(f"Error: {exc}")
        return False
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


### Readers and writers for stage interlock CSV files (Nimbus -> Echo -> Analysis)

# Define headers for all stages. These can handle both mouse and custom samples
stage1_hdr = tuple(['sampleNo', 'samplePlate','sampleWell', 'strainName', 'sampleBarcode', 'assays', 'assayFamilies', 
              'mouseID', 'strainName', 'sex', 'parentLitter', 'parentSequence', 'parentStrain', 'parentDate'])
stage2_hdr = tuple(list(stage1_hdr) + ['dnaPlate', 'dnaWell', 'primer', 'pcrPlate', 'pcrWell'])
stage3_hdr = tuple(list(stage2_hdr) + ['i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well'])


def match_nimbus_to_echo_files(exp: Experiment) -> tuple(list,list,list):
    """ Make sure there is a 1-1 match of Nimbus to Echo-COC file.
        Assumes plate barcodes are guarded.
        Inputs: an Experiment instance
        Outputs: three lists (nimbus files, echo files, barcodes not matches by echo files)
   """
    try:
        nimbus_files = list(Path(exp.get_exp_dir()).glob('Nimbus-*.csv'))
        echo_files = list(Path(exp.get_exp_dir()).glob('Echo_384_COC_00??_*_0.csv'))
        # cleave off the path and file name junk to just get the barcodes
        nbc = frozenset(fp.name.replace('Nimbus-','').replace('.csv','') for fp in nimbus_files)
        ebc = frozenset(fp.name.replace('.csv','').split('_',5)[4] for fp in echo_files)
        xbc  = nbc-ebc # any Nimbus plate files that are missing Nimbus run (output) files
        return [str(fp) for fp in nimbus_files], [str(fp) for fp in echo_files], xbc
    except Exception as exc:
        exp.log('Problem matching Nimbus to Echo files', level='Error')
        exp.log(traceback.print_exc(limit=2), level='Error')    
        return [], [], []

if __name__ == '__main__':
    pass
