#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created May 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comments: adjusted paths relative to app directory
@last_edit: 2022-04-29
@edit_comments: 

Produce picklists for NGS Genotyping pipeline for the Echo robot.
Adds i7i5 Illumina barcodes to each sample in an even fashion.

Each of these stages usually follows a couple of runs of the echovolume.py script.
Each such run combines an Echo survey run with a plate description file into
a single file that describes plates and their contents.
    
Read a Nimbus results file and a number of library files, output PCR1, PCR2 picklists and MiSeq file
This is NGS genotyping stage 2 (nimbus.py is stage 1)

"""

import sys
import csv
import itertools

import bin.file_io as file_io
from bin.util import padwell, unpadwell
from bin.util import Table, CSVTable
from bin.util import output_error
from bin.util import DEAD_VOLS, PLATE_TYPES
from bin.echo_primer import mk_picklist, mk_mytaq_picklist
import bin.stage2report as s2r
   
    
class PicklistSrc:
    """ read a survey & contents file (echovolume.py output) """
    deadvol = { '384PP_AQ_BP': 50, '6RES_AQ_BP2': 700 } # Echo dead volume for plate types 

    def __init__(self, exp, fn, idx=0):
        self.exp = exp
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
            return
        except Exception as exc:
            msg = f"{exc}"
            self.exp.log(f"Error: {msg}")
    
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
            output_error(exc, msg='Error in echo_barcode.xfersrc')

def grouper(xs, kf=lambda x:x[0]):
    """ group pairs: (x,a), (x,b), ... => (x,[a,b, ...]), ... """
    return itertools.groupby(sorted(xs, key=kf), key=kf)
    
global ttno
ttno = 0

#def mk_picklist(fndst, rows):
#    """ output an Echo picklist given the rows (transfer spec) """
#    try:
#        plhdr = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Destination Plate Name,Destination Plate Barcode,Destination Plate Type,Destination Well,Volume".split(',')
#        def rowchk(r):
#            "check the length of a row"
#            if len(r)!=len(plhdr):
#                print("len(r), len(plhdr) =", (len(r), len(plhdr)))
#                print("r =", r)
#            assert len(r)==len(plhdr)
#            return r
#        with open(fndst, "wt", newline='') as dstfd:
#            dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
#            dst.writerow(plhdr)
#            dst.writerows(rowchk(r) for r in rows)
#        return  
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode.mk_picklist')


#def mytaq(wellCount, voltaq, volh2o, plateType='6RES_AQ_BP2'):
#    """ 
#    ***DEPRECATED***
#    Create enough source wells for a Mytaq & H2O picklist, uses the same well until it's empty.
#    Broken - needs a consistent output with mytaq2
#    """
#    try:
#        wells = [r+c for r in "AB" for c in "123"]
#        # volumes in nanolitres
#        dv = PicklistSrc.deadvol[plateType]
#        # initial well volume is 2800uL per well in a full 6 well plate
#        # work in nanolitres
#        wc = [(2800-dv)*1000//v for v in (voltaq, volh2o)]
#        wct, wcw = [(wellCount+c-1)//c for c in wc]
#        plates_required = max([wct//3, wcw//3])
#        # Water is A1, A2, ... while Mytaq is B3, B2, ...
#        # returns wells: Mytaq list, water list
#        return wells[-wct:], wells[:wcw], plates_required
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode.mytaq')


#def mytaq2(wellCount, voltaq, volh2o, plateType='6RES_AQ_BP2', plate_barcodes=None):
#    """ declare source wells for a Mytaq & H2O picklist, cycles through each well sequentially """
#    try:
#        water_wells = ['A'+c for c in '123']
#        taq_wells = ['B'+c for c in '123']
#        #print("water wells", water_wells, file=sys.stdout)
#        # calculate whether we have enough wells, volumes in nanolitres
#        dv = PicklistSrc.deadvol[plateType]
#        # initial well volume is 2800uL per well in a full 6 well plate
#        transfers_per_well = [(2800-dv)*1000//v for v in (voltaq, volh2o)]
#        max_transfers_per_well = min(transfers_per_well)  # worst case scenario
#        taq_wells_required, water_wells_required = [(wellCount+c-1)//c for c in transfers_per_well]
#        plates_required = (max([taq_wells_required, water_wells_required])//3)+1
#        #print('Plates required:',plates_required, file=sys.stdout)
#        # Water is A1, A2, ... while Mytaq is B3, B2, ...
#        # returns wells: Mytaq list, water list
#        tw_pi = []
#        ww_pi = []
#        transfer_count = 0
#        plate_index = 1  # 1st taq/water plate
#        while transfer_count < wellCount:
#            if transfer_count != 0 and transfer_count%(max_transfers_per_well*3)==0:
#                plate_index += 1
#            tw = taq_wells[transfer_count%len(taq_wells)]
#            ww = water_wells[transfer_count%len(water_wells)]
#            if not plate_barcodes:
#                tw_pi.append((tw,plate_index))
#                ww_pi.append((ww,plate_index))
#            else:
#                tw_pi.append((tw,plate_barcodes[plate_index-1]))
#                ww_pi.append((ww,plate_barcodes[plate_index-1]))
#            transfer_count += 1
#        return tw_pi, ww_pi
#    #    tw = [taq_wells[i%len(taq_wells)] for i in range(wellCount)]
#    #    ww = [water_wells[i%len(water_wells)] for i in range(wellCount)]
#        #print(set(tw), set(ww), file=sys.stdout)
#    #    return tw, ww
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode.my_taq2')


#def mk_mytaq_picklist(fn, wells, taq_plate_barcodes, voltaq, volh2o, plateType='6RES_AQ_BP2'):
#    """ create a picklist for Mytaq & H2O transfers """
#    try:
#        dstPlateType = 'Hard Shell 384 well PCR Biorad'
#        well_count = len(wells)
#        # mytaq2 return a list of each well and plate for taq wells (tw) and for water wells (ww).
#        tw_pbcs, ww_pbcs = mytaq2(well_count, voltaq, volh2o, plate_barcodes=taq_plate_barcodes)
#        # now merge the column info together
#        tws = (('Source[1]', bc, plateType, tw) for tw,bc in tw_pbcs)
#        rowstaq = ([x for xs in xss for x in xs]+[voltaq] for xss in zip(tws, (('', w.pcrplate, dstPlateType, w.pcrwell) for w in wells))) 
#        wws = (('Source[1]', bc, plateType, ww) for ww,bc in ww_pbcs)
#        rowsh2o = ([x for xs in xss for x in xs]+[volh2o] for xss in zip(wws, (('', w.pcrplate, dstPlateType, w.pcrwell) for w in wells))) 
#        mk_picklist(fn, (x for xs in (rowstaq, rowsh2o) for x in xs))
#        return
#    except Exception as exc:
#        output_error(exc, msg='echo_barcode.mk_mytaq_picklist')


#def i7i5alloc(exp, vol, wellcount, index_survey_fn='index-svy.csv'):
#    """ allocate vol to wellcount wells with unique barcode pairs - Bob's code, exhausts well before moving on """
#    try:
#        typebc = Table.newtype('BCRecord', "name platebc type well set barcode xxx oligo volume")
#        # typebc = Table.newtype('BCRecord', "well name index indexName oligo volume")
#        tab = CSVTable(typebc, index_survey_fn)
                                                      
#        dv = DEAD_VOLS[PLATE_TYPES['Echo384']]
#        def getvol(x): return x[-1]
#        i7s, i5s = (sorted(((r.barcode, r.set, r.well, float(r.volume)) for r in tab.data if x in r.set), key=getvol, reverse=True) for x in ['_i7F_', '_i5R_'])
#        i7len, i5len = (len(x) for x in (i7s, i5s))
        
#        # create lists of transfer counts from each well
#        i7cnt, i5cnt = ([int((r[-1])*1000-dv)//vol for r in ws] for ws in (i7s, i5s))
#        assert wellcount<=sum(i5cnt)
#        # number of imbalanced wells
#        imb = [x-i5cnt[-1] for x in i5cnt]
#        extra = sum(imb)
#        excess = max(extra+i5len-wellcount, 0)
#        while excess:
#            for i in range(len(imb)):
#                if imb[i]:
#                    imb[i] -= 1
#                    excess -= 1
#                    if not excess:
#                        break
#        extra = sum(imb)
#        base = (wellcount-extra)//i5len
#        alloc = [base+i for i in imb]
    
#        # biggest i5 allocation is fewer than the number of i7 wells
#        # if so, there should be no duplicate i7F, i5R pairs - we check below anyway
#        assert alloc[0]<=i7len 
#        assert sum(alloc)==wellcount
#        i7gen = (x[:3] for xs in itertools.repeat(i7s) for x in xs)
#        i5gen = (x[:3] for cnt, xs in zip(alloc, i5s) for x in [xs]*cnt)
#        res = list(zip(i7gen, i5gen))
#        resset = frozenset(tuple(x[0] for x in xs) for xs in res)
#        assert len(resset)==wellcount # all i7-i5 pairs are unique
#        return res
#    except Exception as exc:
#        exp.log(f"Error: {exc}")
#    return []


def i7i5alloc_rot(exp, vol, wellcount, index_survey_fn='index-svy.csv'):
    """ allocate vol to wellcount wells with unique index pairs - rotates barcode wells to avoid wells being drained """
    try:
        typebc = Table.newtype('BCRecord', "name platebc type well set barcode orderpart oligo volume")
        # typebc = Table.newtype('BCRecord', "well name index indexName oligo volume")
        tab = CSVTable(typebc, exp.get_exp_fp(index_survey_fn))

        dv = DEAD_VOLS[PLATE_TYPES['Echo384']]
        def getvol(x): return x[-1]
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
                exp.log(f"Info: i7s_empty:{i7s_empty},i5s_empty{i5s_empty},used index list:{len(used_index_pairs)}")
                return
            elif len(i5s_empty) == len(i5s):
                exp.log(f"Error: i5 wells contain insufficient volume for the number of experimental wells")
                exp.log(f"Info: i7s_empty:{i7s_empty},i5s_empty{i5s_empty},used index list:{len(used_index_pairs)}")
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
                i7s[i7] += vol
                continue
            i5s[i5] -= vol
            if i5s[i5] < 0:
                i5s_emtpy.add(i5)
                i5s[i5] += vol
                continue

            index_info_pairs.append((i7,i5))
            used_index_pairs.add((i7[0],i5[0]))
        print(index_info_pairs, file=sys.stdout)
        return index_info_pairs
      
    except Exception as exc:
        exp.log(f"Error: {exc}")


def generate_miseq_samplesheet(exp, fmtmiseq, s3tab):
    """ Now that we have all the required info, generate the Illumina Miseq samplesheet, which drives the sequencing """
    try:
        fnmiseq = exp.get_exp_fp(fmtmiseq.format(exp.name))
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
            try:
                gen = [(r.pcrplate+'_'+padwell(r.pcrwell), r.mouseBarcode, '', '',
                       r.i7name, r.i7bc,r.i5name, r.i5bc,'NGSgeno') for r in s3tab.data]
            except AttributeError:
                gen = [(r.pcrplate+'_'+padwell(r.pcrwell), r.sampleBarcode, '', '',
                       r.i7name, r.i7bc,r.i5name, r.i5bc,'NGSgeno') for r in s3tab.data]

            dst.writerows(gen)
        
    except Exception as exc:
        exp.log(f"Error: {exc}")
        return False
    return True

    
def generate_echo_PCR2_picklist(exp, pcr_plate_bcs, index_plate_bcs, taq_water_bcs, amplicon_plate_bcs=None):
    """
    Replacement for main(). Called as a library only.
    User included PCR plates, index plates, taq/water plates, amplicon plates (optional) are user provided
    Called from Experiment.generate_echo_PCR2_picklist_interface()
    """
    try:
        pcr_bcs = []
        amplicon_bcs = []
        taq_bcs = []
        try:
            pcr_bcs = [file_io.guard_pbc(p,silent=True) for p in pcr_plate_bcs]
        except Exception as e:
            exp.log(f"{e}")
        try:
            index_bcs = [file_io.guard_pbc(i, silent=True) for i in index_plate_bcs]
        except Exception as e:
            exp.log(f"{e}")
        try:
            amplicon_bcs = [file_io.guard_pbc(d,silent=True) for d in amplicon_plate_bcs] 
        except Exception as e:
            exp.log(f"{e}")
        try:                          
            taq_bcs = [file_io.guard_pbc(t, silent=True) for t in taq_water_bcs]
        except Exception as e:
            exp.log(f"{e}")
        pcrfmt = "PCR{{}}-picklist_{}.csv".format(exp.name)
    
        fnstage2, fnstage3 = "Stage2.csv", "Stage3.csv"
    
        ### i7i5 barcodes
        # read Stage2 csv file
        #sampleNumber	samplePlate	sampleWell	sampleBarcode	assays	assayFamilies	strain	sex	dnaplate	dnawell	primer	pcrplate	pcrwell
        #typebc = Table.newtype('S2Record',
        s2tab = CSVTable('S2Rec', exp.get_exp_fp(fnstage2))           
            
        index_vol = exp.transfer_volumes['INDEX_VOL']  # 175 nanolitres
        index_alloc = i7i5alloc_rot(exp, index_vol, len(s2tab.data))
        #print(index_alloc, file=sys.stdout)
            
        s3flds = s2tab.tt._fields+('i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well', 'index_plate')
        S3Rec = Table.newtype('S3Rec', s3flds)
        
        s3tab = Table(S3Rec, ([x for xs in (p1, p2[:3], p3[:4]) for x in xs] for p1, (p2, p3) in zip(s2tab.data, index_alloc)), headers=s3flds) 
        # output Stage 3 CSV file - used for custom and mouse samples. Amplicons are handled separately
        s3tab.csvwrite(exp.get_exp_fp(fnstage3))
            
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

        # plate barcode for i7i5 plate is unknown
        rowi7 = ((src_dict[r.index_plate], r.index_plate, '384PP_AQ_BP', r.i7well, dest_dict[r.pcrplate], r.pcrplate, '384PP_AQ_BP', r.pcrwell, index_vol) for r in s3tab.data)
        rowi5 = ((src_dict[r.index_plate], r.index_plate, '384PP_AQ_BP', r.i5well, dest_dict[r.pcrplate], r.pcrplate, '384PP_AQ_BP', r.pcrwell, index_vol) for r in s3tab.data)
        fn_index = exp.get_exp_fp(pcrfmt.format("2_index"))
        mk_picklist(exp, fn_index, (r for rs in (rowi7, rowi5) for r in rs))    
        
        # also PCR2 water and Taq
        fn_taqwater = exp.get_exp_fp(pcrfmt.format("2_taqwater"))
        mk_mytaq_picklist(exp, fn_taqwater, s3tab.data, taq_bcs, 2000, 650)
        
        # MiSeq file - experiment name, name format, ngid(dir), s3 data, verbosity
        fmtmiseq = "MiSeq-{}.csv"
        success = generate_miseq_samplesheet(exp, fmtmiseq, s3tab)
        
        fns3 = "Stage3.html"
        #s2r.build_report(exp, fns3, s3tab.tt._fields, s3tab.data)    
        
    except Exception as exc:
        exp.log(f"Error: {exc}")
        return False
    return True

#def main():
#    """
#    Main program. Interpret command line.
#    Stage 1: Input DNA plate descriptions (Nimbus output)
#    Stage 2: Check Assays, output PCR1 Echo picklist
#    Stage 3: Attach sample barcodes to all wells for PCR2 stage, output picklist.

#    Returns
#    -------
#    None.

#    """
#    global args
#    # bcdef = "i7i5_plate_layout_20*.csv"
#    defxname = "NGS experiment"
#    fmtmiseq = "MiSeq-{}.csv"
#    deflib = os.path.join("..", "library")
#    parser = argparse.ArgumentParser()
#    parser.add_argument('-w', '--workdir', default='', help='specify a working directory (default $NGSDIR or $HOME/NGSDIR)')
#    parser.add_argument('-l', '--library', default=deflib, help="library directory (default={}), best left alone".format(deflib))
#    parser.add_argument('-v', '--verbose', action="store_true", help='lots of reporting')
#    parser.add_argument('-t', '--taq', nargs='+', help='Mytaq and water plate survey filename')
#    parser.add_argument('--custom', action='store_true', help='Pipeline is running in custom sample mode')
#    grp2 = parser.add_argument_group('Stage 2 - dna & assays')
#    grp3 = parser.add_argument_group('Stage 3 - Miseq barcodes')
#    grp2.add_argument('-p', '--prim', help='primer/assay plate content filename')
#    grp2.add_argument('-d', '--pcr', nargs='+',  help='PCR (destination) plate barcodes')
#    # grp2.add_argument('-r', '--resume', help='i7i5 plate last well used - for example D15,H19.')
#    grp2.add_argument('dna', nargs="+", help='barcode for 1 or more 384-well sample source (DNA - Nimbus output) plates')
#    parser.add_argument('-D', '--dnadiff', action="store_true", help='allow different barcode in DNA files')
#    grp3.add_argument('-i', '--i7i5', help='i7i5 plate survey filename.')
#    grp3.add_argument('-x', '--xname', default=defxname, help="experiment name - used in MiSeq file (default={})".format(defxname))
#    args = parser.parse_args()
#    # protect plate barcode inputs
#    try:
#        args.taq = [file_io.guard_pbc(t) for t in args.taq]
#    except file_io.ExistingGuardError as e:
#        print(e, file=sys.stdout)
    
    
#    def getdir(path):
#        """ find the work directory - looks in $NGSGENO and locality of program if the current directory isn't suitable. 
#        Almost certainly not necessary. """
#        check = lambda dn: os.path.isdir(dn) and os.path.isdir(os.path.join(dn, args.library))
#        if check(path):
#            return os.path.normpath(path)
#        try:
#            # cwd = os.path.join('..', os.path.basename(os.getcwd()))
#            parent = os.path.dirname(os.path.dirname(sys.argv[0])) # parent of program's directory
#            # directory and a proximate library directory are needed
#            import platform
#            homestr = "USERPROFILE" if platform.system()=="Windows" else "HOME"
#            envopt = os.path.join(os.environ(homestr), os.environ['NGSGENO'] if 'NGSGENO' in os.environ else 'NGSGENO')
#            # the following may fail!
#            return next(fn for fn in (os.path.normpath(os.path.join(d, args.workdir)) for d in (parent, envopt)) if check(fn))
#        except:
#            alert("target directory {} not found.\n".format(args.workdir))
#            return args.workdir
    
#    try:  # wrap whole program
#        if args.workdir not in ['', '.']:
#            os.chdir(getdir(args.workdir))
    
#        ngid = os.path.basename(os.getcwd()) # should be 8 characters -YYYYMMDD
#        pcrfmt = "PCR{{}}-picklist_{}.csv".format(ngid)
    
#        fnstage2, fnstage3 = "Stage2.csv", "Stage3.csv"
    
#        ### i7i5 barcodes
#        # read Stage2 csv file
#        s2tab = CSVTable('Stage2', fnstage2)           
            
#        bcvol = 175 # nanolitres
#        bcalloc = i7i5alloc_rot(bcvol, len(s2tab.data))
#        print(bcalloc, file=sys.stdout)
        
#        def mks3rec(s2rec, bcrec):
#            return 
            
#        s3flds = s2tab.tt._fields+('i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well')
#        S3Rec = Table.newtype('S3Rec', s3flds)
#        s3tab = Table(S3Rec, ([x for xs in (p1, p2[:3], p3[:3]) for x in xs] for p1, (p2, p3) in zip(s2tab.data, bcalloc)), headers=s3flds) 
#        # output Stage 3 CSV file - used in Stage 3 below
#        s3tab.csvwrite(fnstage3)                            
#        if args.verbose:
#            print("created:", fnstage3) 
            
#        # write Barcode Picklists
#        # make the destination dictionary for pcrplates
#        ddict = {}
#        for r in s3tab.data:
#            if r.pcrplate not in ddict:
#                ddict[r.pcrplate] = "Destination[{}]".format(len(ddict)+1)
#        # plate barcode for i7i5 plate is unknown
#        rowi7 = (("Source[1]", "", '384PP_AQ_BP', r.i7well, ddict[r.pcrplate], r.pcrplate, '384PP_AQ_BP', r.pcrwell, bcvol) for r in s3tab.data)
#        rowi5 = (("Source[1]", "", '384PP_AQ_BP', r.i5well, ddict[r.pcrplate], r.pcrplate, '384PP_AQ_BP', r.pcrwell, bcvol) for r in s3tab.data)
#        fnbc = pcrfmt.format("2_barcode")
#        mk_picklist(fnbc, (r for rs in (rowi7, rowi5) for r in rs))
#        if args.verbose:
#            print("created:", fnbc) 
    
        
#        # also PCR2 water and Taq
#        fndst = pcrfmt.format("2_taqwater")
#        mk_mytaq_picklist(fndst, s3tab.data, args.taq, 2000, 650)
#        if args.verbose:
#            print('created:', fndst)
            
#        # MiSeq file - experiment name, name format, ngid(dir), s3 data, verbosity 
#        generate_miseq_samplesheet(args.xname, fmtmiseq, ngid, s3tab, verbose=args.verbose)
        
#        fns3 = "Stage3.html"
#        s2r.build_report(fns3, s3tab.tt._fields, s3tab.data,is_custom=args.custom)    
#        if args.verbose:
#            print('created:', fns3)

#        return 
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode main body')

def main():
    """ library only """
    return


if __name__ == '__main__':
    main()
