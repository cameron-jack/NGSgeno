#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created May 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comments: Functions only - not called as a program
@last_edit: 2022-07-29
@edit_comments: Deprecated main()

Produce picklists for NGS Genotyping pipeline for the Echo robot.
There are two stages: 
    stage 1 takes sample DNA from Nimbus output plates, assay primer plate(s)
    and a Mytaq+H2O plate and produced a stahe 1 PCR plate.
    stage 2 adds barcode primers and more Mytaq+H2O to the plates, and also outputs
    a MiSeq file.

Each of these stages usually follows a couple of runs of the echovolume.py script.
Each such run combines an Echo survey run with a plate description file into
a single file that describes plates and their contents.
    
Read a Nimbus results file and a number of library files, output PCR1, PCR2 picklists and MiSeq file
This is NGS genotyping stage 2 (nimbus.py is stage 1)

It needs primers and barcodes files (see code)
"""

# TO DO:
#   parameterise the various volumes
#   use PCR plates barcodes?
#   tidy up the code - there are too many unused functions.


#import re
import os
import sys
import json
import csv
import glob
import collections
import itertools
import argparse
from pathlib import Path
from copy import deepcopy

import bin.file_io as file_io

from bin.util import padwell, unpadwell
from bin.util import Table, CSVTable
from bin.util import output_error
       
    
class PicklistSrc:
    """ read a survey & contents file (echovolume.py output) """
    deadvol = { '384PP_AQ_BP': 50, '6RES_AQ_BP2': 700 } # Echo dead volume for plate types 

    def __init__(self, fn, idx=0):
        try:
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
        except Exception as exc:
            output_error(exc, msg='echo_primer.PicklistSrc.__init__')
    
    def xfersrc(self, l, vol, depleted):
        """ transfer data for liquid l - reduces volume """
        try:
            assert l in self.data, 'no primer well for '+l
            dx = self.data
            while dx[l] and dx[l][0][-1]-vol<PicklistSrc.deadvol[dx[l][0][2]]:
                dx[l].pop(0) # discard depleted wells
            # assert dx[l], "primer depleted for "+l
            if not dx[l]:
                depleted.update([l])
                return ['']*4
            dx[l][0][-1] -= vol
            return dx[l][0][:4]
        except Exception as exc:
            output_error(exc, msg='Error in echo_primer.PicklistSrc.xfersrc')
    

def grouper(xs, kf=lambda x:x[0]):
    "group pairs: (x,a), (x,b), ... => (x,[a,b, ...]), ..."
    return itertools.groupby(sorted(xs, key=kf), key=kf)
    
global ttno
ttno = 0

def file_get_check(fids, fmt):
    """ collect Echo DNA plate file names, match plate ids, remove duplicates """
    try:
        dups = [fids[i] for i in range(1, len(fids)) if fids[i] in fids[:i]]
        fx = [fids[i] for i in range(len(fids)) if fids[i] not in fids[:i]]
        if dups:
            print("Duplicate DNA plate IDs ignored:", ' '.join(dups), file=sys.stdout)    
    
        globs = dict((pid, sorted(glob.glob(fmt.format(pid)))) for pid in fx)
        nofile = [fid for fid, fns in globs.items() if not fns]
        if nofile:
            print('\n'.join(fmt.format(fid)+': no file found.' for fid in nofile), file=sys.stdout)
            exit(1)
        return collections.OrderedDict((fid, sorted(ps)[-1]) for fid, ps in globs.items())
    except Exception as exc:
        output_error(exc, msg='Error in echo_primer.file_get_check')


def mk_picklist(fndst, rows):
    """ output an Echo picklist given the rows (transfer spec) """
    try:
        plhdr = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"+\
                "Destination Plate Name,Destination Plate Barcode,Destination Plate Type,"+\
                "Destination Well,Volume".split(',')
        def rowchk(r):
            "check the length of a row"
            if len(r)!=len(plhdr):
                print("len(r), len(plhdr) =", (len(r), len(plhdr)))
                print("r =", r)
            assert len(r)==len(plhdr)
            return r
        with open(fndst, "wt", newline='') as dstfd:
            dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            dst.writerow(plhdr)
            dst.writerows(rowchk(r) for r in rows)
        return
    except Exception as exc:
        output_error(exc, msg='Error in echo_primer.mk_picklist')


def mytaq(wellCount, voltaq, volh2o, plateType='6RES_AQ_BP2'):
    """ 
    ***DEPRECATED***
    Create enough source wells for a Mytaq & H2O picklist, uses the same well until it's empty.
    Broken - needs a consistent output with mytaq2
    """
    try:
        wells = [r+c for r in "AB" for c in "123"]
        # volumes in nanolitres
        dv = PicklistSrc.deadvol[plateType]
        # initial well volume is 2800uL per well in a full 6 well plate
        # work in nanolitres
        wc = [(2800-dv)*1000//v for v in (voltaq, volh2o)]
        wct, wcw = [(wellCount+c-1)//c for c in wc]
        plates_required = max([wct//3, wcw//3])
        # Water is A1, A2, ... while Mytaq is B3, B2, ...
        # returns wells: Mytaq list, water list
        return wells[-wct:], wells[:wcw], plates_required
    except Exception as exc:
        output_error(exc, msg='Error in echo_primer.mytaq')


def mytaq2(wellCount, voltaq, volh2o, plateType='6RES_AQ_BP2', plate_barcodes=None):
    """ declare source wells for a Mytaq & H2O picklist, cycles through each well sequentially """
    try:
        water_wells = ['A'+c for c in '123']
        taq_wells = ['B'+c for c in '123']
        #print("water wells", water_wells, file=sys.stdout)
        # calculate whether we have enough wells, volumes in nanolitres
        dv = PicklistSrc.deadvol[plateType]
        # initial well volume is 2800uL per well in a full 6 well plate
        transfers_per_well = [(2800-dv)*1000//v for v in (voltaq, volh2o)]
        max_transfers_per_well = min(transfers_per_well)  # worst case scenario
        taq_wells_required, water_wells_required = [(wellCount+c-1)//c for c in transfers_per_well]
        plates_required = (max([taq_wells_required, water_wells_required])//3)+1
        #print('Plates required:',plates_required, file=sys.stdout)
        # Water is A1, A2, ... while Mytaq is B3, B2, ...
        # returns wells: Mytaq list, water list
        tw_pi = []
        ww_pi = []
        transfer_count = 0
        plate_index = 1  # 1st taq/water plate
        while transfer_count < wellCount:
            if transfer_count != 0 and transfer_count%(max_transfers_per_well*3)==0:
                plate_index += 1
            tw = taq_wells[transfer_count%len(taq_wells)]
            ww = water_wells[transfer_count%len(water_wells)]
            if not plate_barcodes:
                tw_pi.append((tw,plate_index))
                ww_pi.append((ww,plate_index))
            else:
                tw_pi.append((tw,plate_barcodes[plate_index-1]))
                ww_pi.append((ww,plate_barcodes[plate_index-1]))
            transfer_count += 1
        return tw_pi, ww_pi
        #    tw = [taq_wells[i%len(taq_wells)] for i in range(wellCount)]
        #    ww = [water_wells[i%len(water_wells)] for i in range(wellCount)]
            #print(set(tw), set(ww), file=sys.stdout)
        #    return tw, ww
    except Exception as exc:
        output_error(exc, msg='Error in echo_primer.mytaq2')


def mk_mytaq_picklist(fn, wells, taq_plate_barcodes, voltaq, volh2o, plateType='6RES_AQ_BP2'):
    """ create a picklist for Mytaq & H2O transfers """
    try:
        dstPlateType = 'Hard Shell 384 well PCR Biorad'
        well_count = len(wells)
        # mytaq2 return a list of each well and plate for taq wells (tw) and for water wells (ww).
        tw_pbcs, ww_pbcs = mytaq2(well_count, voltaq, volh2o, plate_barcodes=taq_plate_barcodes)
        # now merge the column info together
        tws = (('Source[1]', bc, plateType, tw) for tw,bc in tw_pbcs)
        rowstaq = ([x for xs in xss for x in xs]+[voltaq] for xss in zip(tws, (('', w.pcrplate, dstPlateType, w.pcrwell) for w in wells))) 
        wws = (('Source[1]', bc, plateType, ww) for ww,bc in ww_pbcs)
        rowsh2o = ([x for xs in xss for x in xs]+[volh2o] for xss in zip(wws, (('', w.pcrplate, dstPlateType, w.pcrwell) for w in wells))) 
        mk_picklist(fn, (x for xs in (rowstaq, rowsh2o) for x in xs))
        return
    except Exception as exc:
        output_error(exc, msg='Error in echo_primer.mk_mytaq_picklist')


def generate_echo_PCR1_picklist(exp, dna_plate_bcs, pcr_plate_bcs, taq_water_bcs):
    """
    Entry point. Takes an experiment instance plus plate barcodes for dna plates, PCR plates, primer plates, taq/water plates.
    """
    try:
        pcr_bcs = [file_io.guard_pbc(p,silent=True) for p in pcr_plate_bcs]
    except Exception as e:
        exp.log(e)
    try:
        dna_bcs = [file_io.guard_pbc(d,silent=True) for d in dna_plate_bcs] 
    except Exception as e:
        exp.log(e)
    try:
        primer_bcs = [file_io.guard_pbc(p,silent=True) for p in exp.sample_plate_location if p['purpose'] == 'primers']
    except Exception as e:
        exp.log(e)
    try:
        taq_bcs = [file_io.guard_pbc(t) for t in taq_water_bcs]
    except Exception as e:
        exp.log(e)

    # read Nimbus created DNA plates into experiment
    #RecordId	TRackBC	TLabwareId	TPositionId	SRackBC	SLabwareId	SPositionId
    #1	p2021120604p	Echo_384_COC_0001	A1	p2111267p	ABg_96_PCR_NoSkirt_0001	A1
    for d in dna_bcs:
        fp = exp.get_exp_fp("Echo_384_COC_0001_" + d + "_0.csv")
        if not Path(fp).is_file():
            exp.log(f"{d} has no matching Echo_384_COC file from the Nimbus", level='error')
        else:
            with open(fp, 'rt') as f:
                exp.plate_location_sample[d] = {}
                exp.plate_location_sample[d]['purpose'] = 'dna'
                exp.plate_location_sample[d]['wells'] = set()
                for i, line in enumerate(f):
                    if i == 0:  # header
                        continue
                    cols = [c.strip() for c in line.split(',')]
                    source_pos = unpadwell(cols[-1])
                    source_plate = cols[-3]
                    dest_plate = cols[1]
                    dest_pos = unpadwell(cols[3])
                    if dest_plate != d:
                        exp.log(f"{d} doesn't match {dest_plate} as declared in Echo_384_COC file: {fp}", level="error")
                    exp.plate_location_sample[d]['wells'].add(dest_pos)
                    exp.plate_location_sample[d][dest_pos] = deepcopy(exp.plate_location_sample[source_plate][source_pos])

        pcrfmt = "PCR{{}}-picklist-{}.csv".format(exp.name)

        dnafns = file_get_check(dna_bcs, "Echo_384_COC_0001_{0}_?.csv") 
        nimcolmap = (("TRackBC", "dstplate"), ("TPositionId", 'dstwell'), ('SRackBC', 'srcplate'), ('SPositionId', 'srcwell'))
        typenim = Table.csvtype(dnafns[dna_bcs[0]], 'NimRec', hdrmap=nimcolmap)
        nimbusTables = [CSVTable(typenim, fn) for fid, fn in dnafns.items()]
        def badnimbc(param):
            bc, t = param
            return bc, set(r.dstplate for r in t.data if r.dstplate!=bc)
        badnimdata = [p for p in map(badnimbc, zip(dna_bcs, nimbusTables)) if p[1]]
        if badnimdata:
            message = "Inconsistent barcode(s) in DNA (Nimbus output) file.\n"
            for bc, badset in badnimdata:
                data = ' '.join(sorted(badset))
                message += "expect barcode = {}, file contains {}\n".format(bc, data)
            exp.log(message, level="error")
        dnadict = dict(((x.srcplate, x.srcwell), (x.dstplate, x.dstwell)) for nt in nimbusTables for x in nt.data)    
        dnas = [CSVTable('S1Rec', 'Stage1-P{}.csv'.format(dnabc)) for dnabc in dna_bcs]
        
        primerTable = CSVTable("PPRec", primer_bcs[0], fields="spn spbc spt well primer volume")
        primset = sorted(frozenset(x.primer for x in primerTable.data if x.primer != ''))
        
        # primer family dict to list of primers within family group
        pfdict = dict((k.strip(), list([p.strip() for p in g])) for k, g in itertools.groupby(primset, key=lambda x:x.split('_',1)[0]))
        
        # create record for PCR wells  
        wgenflds =  dnas[0].tt._fields + ('dnaplate', 'dnawell') + ('primer',)
        
        # Is this the plating of DNA samples for each assay? Yes. It is.
        wgen = [xs+dnadict[(xs.EPplate, xs.EPwell)]+(x,) for xss in dnas for xs in xss.data for f in xs.assayFamilies.split(';') if f in pfdict for x in pfdict[f]]

        # allocate PCR plate wells - leave last 3 empty
        wells = [r+str(c+1) for c in range(24) for r in"ABCDEFGHIJKLMNOP"][:-3]
        pcrwellgen = ((p, w) for p in pcr_bcs for w in wells)
        # allocate PCR well for each sample
        
        s2flds = wgenflds+('pcrplate', 'pcrwell')
        S2Rec = Table.newtype('S2Rec', s2flds)
        s2tab = Table(S2Rec, ([x for xs in rx for x in xs] for rx in zip(wgen, pcrwellgen)), headers=s2flds) 
        # output Stage 2 CSV file - used in Stage 3 below
        s2tab.csvwrite("Stage2.csv")

        # PCR1 picklists: DNA/sample, Primer and Taq/H20
    
        dstPlateType = 'Hard Shell 384 well PCR Biorad'
        srcPlateType = "384PP_AQ_BP"
        # sample DNA picklist
        ddict = dict((pbc, "Destination[{}]".format(i)) for i, pbc in enumerate(pcr_bcs, start=1))
        sdict = dict((pbc, "Source[{}]".format(i))for i, pbc in enumerate(dna_bcs, start=1))
        volume = 200
        fn = pcrfmt.format("1dna")
        gen = ((sdict[r.dnaplate], r.dnaplate, srcPlateType, r.dnawell,
                ddict[r.pcrplate], r.pcrplate, dstPlateType, r.pcrwell, volume)
                for r in s2tab.data)
        mk_picklist(fn, gen)
        
        # Primer Picklist
        # [Source[1]', '', '384PP_AQ_BP', 'H6', 'Destination[1]', '3121', 'Hard Shell 384 well PCR Biorad', 'A1', 500]
        primsrc = PicklistSrc("primer-svy.csv", idx=4) # same name as in cgi-nimbus2.py
        volume = 500        
        fn = pcrfmt.format("1prim")
        depleted = collections.Counter()
        primer_uses = collections.Counter()
        primer_output_rows = []
        for r in s2tab.data:
            if r.primer not in primsrc.data:
                print('Cannot find', r.primer, 'in known primers', file=sys.stdout)
                continue
            # cycle through available primer wells
            primer_index = primer_uses[r.primer] % len(primsrc.data[r.primer])
            primer_output_rows.append(primsrc.data[r.primer][primer_index][:4] +
                    [ddict[r.pcrplate], r.pcrplate, dstPlateType, r.pcrwell, volume])
            primer_uses[r.primer] += 1
        # or you can drain each primer well before moving on...
        #gen = ([f for fs in (primsrc.xfersrc(r.primer, volume, depleted), 
        #        (ddict[r.pcrplate], r.pcrplate, dstPlateType, r.pcrwell, volume)) for f in fs]
        #       for r in s2tab.data)
        #mk_picklist(fn, gen)
        mk_picklist(fn, primer_output_rows)

        if len(depleted):
            # Outputting HTML like this seems to display OK - though it shoule be done properly (in cgi-nimbus2.py)
            print("<h1 style='color: red;'>Warning: depleted primer wells ...</h1>\n<pre>"+"\n".join("{:5d} {}".format(n, pr) for pr, n in depleted.items()), "\n</pre>\n")
            advice="""You can fix this by clicking the "back" button, 
            adding the required primer to wells in the primer plate,
            adding relevant information to the primer description file for this sample,
            using the Echo robot to re-survey the primer plate and copying the new
            survey file to the NGSgeno sample folder and choosing the right survey and primer plate
            description files.<br>Good luck!
            
            """
            print("<p>"+advice+"</p>")
        
        
        # also PCR1 water and Taq
        fndst = pcrfmt.format("1TaqWater")
        mk_mytaq_picklist(fndst, s2tab.data, taq_bcs, 1000, 300)
        return True
        

    
def main():
    """
    Main program. Interpret command line.
    Stage 1: Input DNA plate descriptions (Nimbus output)
    Stage 2: Check Assays, output PCR1 Echo picklist

    Returns
    -------
    None.

    """
    global args
    defxname = "NGS experiment"
    fmtmiseq = "MiSeq-{}.csv"
    deflib = os.path.join("..", "library")
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--workdir', default='', help='specify a working directory (default $NGSDIR or $HOME/NGSDIR)')
    parser.add_argument('-l', '--library', default=deflib, help="library directory (default={}), best left alone".format(deflib))
    parser.add_argument('-v', '--verbose', action="store_true", help='lots of reporting')
    parser.add_argument('-t', '--taq', nargs='+', help='Mytaq and water plate survey filename')
    parser.add_argument('--custom', action='store_true', help='Pipeline is running in custom sample mode')
    grp2 = parser.add_argument_group('Stage 2 - dna & assays')
    grp2.add_argument('-p', '--prim', help='primer/assay plate content filename')
    grp2.add_argument('-d', '--pcr', nargs='+',  help='PCR (destination) plate barcodes')
    # grp2.add_argument('-r', '--resume', help='i7i5 plate last well used - for example D15,H19.')
    grp2.add_argument('dna', nargs="+", help='barcode for 1 or more 384-well sample source (DNA - Nimbus output) plates')
    parser.add_argument('-D', '--dnadiff', action="store_true", help='allow different barcode in DNA files')
    args = parser.parse_args()
    try:
        args.pcr = [file_io.guard_pbc(p) for p in args.pcr]
    except file_io.ExistingGuardError as e:
        print(e, file=sys.stdout)
    try:
        args.dna = [file_io.guard_pbc(d) for d in args.dna] 
    except file_io.ExistingGuardError as e:
        print(e, file=sys.stdout)
    try:
        args.taq = [file_io.guard_pbc(t) for t in args.taq]
    except file_io.ExistingGuardError as e:
        print(e, file=sys.stdout)

    
    def getdir(path):
        "find the work directory - looks in $NGSGENO and locality of program if the current directory isn't suitable."
        check = lambda dn: os.path.isdir(dn) and os.path.isdir(os.path.join(dn, args.library))
        if check(path):
            return os.path.normpath(path)
        try:
            # cwd = os.path.join('..', os.path.basename(os.getcwd()))
            parent = os.path.dirname(os.path.dirname(sys.argv[0])) # parent of program's directory
            # directory and a proximate library directory are needed
            import platform
            homestr = "USERPROFILE" if platform.system()=="Windows" else "HOME"
            envopt = os.path.join(os.environ(homestr), os.environ['NGSGENO'] if 'NGSGENO' in os.environ else 'NGSGENO')
            # the following may fail!
            return next(fn for fn in (os.path.normpath(os.path.join(d, args.workdir)) for d in (parent, envopt)) if check(fn))
        except:
            print("target directory {} not found.\n".format(args.workdir), file=sys.stdout)
            return args.workdir
    
    try:
        if args.workdir not in ['', '.']:
            os.chdir(getdir(args.workdir))
    
        ngid = os.path.basename(os.getcwd()) # should be 8 characters -YYYYMMDD
        pcrfmt = "PCR{{}}-picklist-{}.csv".format(ngid)
    
        fnstage2 = "Stage2.csv"
        if args.pcr:
            # read Nimbus output files - record of DNA samples in 384 well plates
            # nimcol = "RecordId","TRackBC","TLabwareId","TPositionId","SRackBC","SLabwareId","SPositionId"
            dnafns = file_get_check(args.dna, "Echo_384_COC_0001_{0}_?.csv")
            if args.verbose:
                print("dnafns =", dnafns)
            nimcolmap = (("TRackBC", "dstplate"), ("TPositionId", 'dstwell'), ('SRackBC', 'srcplate'), ('SPositionId', 'srcwell'))
            typenim = Table.csvtype(dnafns[args.dna[0]], 'NimRec', hdrmap=nimcolmap)
            nimbusTables = [CSVTable(typenim, fn) for fid, fn in dnafns.items()]
            def badnimbc(param):
                bc, t = param
                return bc, set(r.dstplate for r in t.data if r.dstplate!=bc)
            badnimdata = [p for p in map(badnimbc, zip(args.dna, nimbusTables)) if p[1]]
            if badnimdata and not args.dnadiff:
                message = "Inconsistent barcode(s) in DNA (Nimbus output) file.\n"
                for bc, badset in badnimdata:
                    data = ' '.join(sorted(badset))
                    message += "expect barcode = {}, file contains {}\n".format(bc, data)
                print(message, file=sys.stdout)
                raise Exception(message='echo_primer: Inconsistent barcodes in DNA (Nimbus output) file')
        
            dnadict = dict(((x.srcplate, x.srcwell), (x.dstplate, x.dstwell)) for nt in nimbusTables for x in nt.data)
            #print('<<< DNA dictionary >>>', dnadict, file=sys.stdout)
            # find and read Stage1 files - see also primercheck.py and its use in cgi-nimbus2.py
            #s1fns = ['Stage1-P{}.csv'.format(dnabc) for dnabc in args.dna]
            #s1type = Table.csvtype(s1fns[0], 'S1Rec')
            dnas = [CSVTable('S1Rec', 'Stage1-P{}.csv'.format(dnabc)) for dnabc in args.dna]
    
            # read primer plate info - including Echo survey 
            primerTable = CSVTable("PPRec", args.prim, fields="spn spbc spt well primer volume")
            if args.verbose:
                print("Loaded primer plate file:", args.prim)
            primset = sorted(frozenset(x.primer for x in primerTable.data if x.primer != ''))

            # primer family dict to list of primers within family group
            pfdict = dict((k.strip(), list([p.strip() for p in g])) for k, g in itertools.groupby(primset, key=lambda x:x.split('_',1)[0]))

            # create record for PCR wells  
            wgenflds =  dnas[0].tt._fields + ('dnaplate', 'dnawell') + ('primer',)
        
            # Is this the plating of DNA samples for each assay? Yes. It is.
            wgen = [xs+dnadict[(xs.EPplate, xs.EPwell)]+(x,) for xss in dnas for xs in xss.data for f in xs.assayFamilies.split(';') if f in pfdict for x in pfdict[f]]
            #print("wgen",wgen)
            #('9', '80201', 'A2', '105584300122', 'Ahr-fl;Cd79a-cre_MUT;Cd79a-cre_WT;Lef1_MUT;Lef1_WT;Tcf7-fl', 'Ahr-fl;Cd79a-cre;Lef1;Tcf7-fl', '2021090101', 'C1', 'Ahr-fl')
        
            # combine nimbus files if there is more than one
            # nimbusTable = nimbusTables[0] if len(nimbusTables)==1 else Table(typenim, (x for t in nimbusTables for x in t.data)) # fix multi-tables !!!
            # could check that TRackBC matches ID in plate file
            # retrieve mouse data ...
            #    first get the 96-well plate barcodes
               
            # allocate PCR plate wells - leave last 3 empty
            wells = [r+str(c+1) for c in range(24) for r in"ABCDEFGHIJKLMNOP"][:-3]
            pcrwellgen = ((p, w) for p in args.pcr for w in wells)
            # allocate PCR well for each sample
        
            s2flds = wgenflds+('pcrplate', 'pcrwell')
            S2Rec = Table.newtype('S2Rec', s2flds)
            s2tab = Table(S2Rec, ([x for xs in rx for x in xs] for rx in zip(wgen, pcrwellgen)), headers=s2flds) 
            # output Stage 2 CSV file - used in Stage 3 below
            s2tab.csvwrite(fnstage2)                            
            if args.verbose:
                print("created:", fnstage2) 
   
            # PCR1 picklists: DNA/sample, Primer and Taq/H20
    
            dstPlateType = 'Hard Shell 384 well PCR Biorad'
            srcPlateType = "384PP_AQ_BP"
            # sample DNA picklist
            ddict = dict((pbc, "Destination[{}]".format(i)) for i, pbc in enumerate(args.pcr, start=1))
            sdict = dict((pbc, "Source[{}]".format(i))for i, pbc in enumerate(args.dna, start=1))
            volume = 200
            fn = pcrfmt.format("1dna")
            gen = ((sdict[r.dnaplate], r.dnaplate, srcPlateType, r.dnawell,
                    ddict[r.pcrplate], r.pcrplate, dstPlateType, r.pcrwell, volume)
                   for r in s2tab.data)
            mk_picklist(fn, gen)
        
            # Primer Picklist
            # [Source[1]', '', '384PP_AQ_BP', 'H6', 'Destination[1]', '3121', 'Hard Shell 384 well PCR Biorad', 'A1', 500]
            primsrc = PicklistSrc("primer-svy.csv", idx=4) # same name as in cgi-nimbus2.py
            volume = 500        
            fn = pcrfmt.format("1prim")
            depleted = collections.Counter()
            primer_uses = collections.Counter()
            primer_output_rows = []
            for r in s2tab.data:
                if r.primer not in primsrc.data:
                    print('Cannot find', r.primer, 'in known primers', file=sys.stdout)
                    continue
                # cycle through available primer wells
                primer_index = primer_uses[r.primer] % len(primsrc.data[r.primer])
                primer_output_rows.append(primsrc.data[r.primer][primer_index][:4] +
                        [ddict[r.pcrplate], r.pcrplate, dstPlateType, r.pcrwell, volume])
                primer_uses[r.primer] += 1
            # or you can drain each primer well before moving on...
            #gen = ([f for fs in (primsrc.xfersrc(r.primer, volume, depleted), 
            #        (ddict[r.pcrplate], r.pcrplate, dstPlateType, r.pcrwell, volume)) for f in fs]
            #       for r in s2tab.data)
            #mk_picklist(fn, gen)
            mk_picklist(fn, primer_output_rows)

            if len(depleted):
                # Outputting HTML like this seems to display OK - though it shoule be done properly (in cgi-nimbus2.py)
                print("<h1 style='color: red;'>Warning: depleted primer wells ...</h1>\n<pre>"+"\n".join("{:5d} {}".format(n, pr) for pr, n in depleted.items()), "\n</pre>\n")
                advice="""You can fix this by clicking the "back" button, 
                adding the required primer to wells in the primer plate,
                adding relevant information to the primer description file for this sample,
                using the Echo robot to re-survey the primer plate and copying the new
                survey file to the NGSgeno sample folder and choosing the right survey and primer plate
                description files.<br>Good luck!
            
                """
                print("<p>"+advice+"</p>")
            if args.verbose:
                print('created:', fn)
        
            # also PCR1 water and Taq
            fndst = pcrfmt.format("1TaqWater")
            mk_mytaq_picklist(fndst, s2tab.data, args.taq, 1000, 300)
            if args.verbose:
                print('created:', fndst)
        
        return  
    except Exception as exc:
        output_error(exc, msg='Error in echo_primer code body')


if __name__ == '__main__':
    main()
