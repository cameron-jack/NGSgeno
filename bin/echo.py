#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created May 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.8
@version_comments: extensive changes to seamlessly support custom and mouse pipelines, plate dead volumes increased
@last_edit: 2021-07-08
@edit_comments: better handling of paths, multiple taq/water plates, mytaq2() now handles well and plate counting and allocation tasks, mytaq() is now deprecated

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

def fail(x):
    "display message and stop"
    print(x, file=sys.stderr)
    sys.exit(1)
    
alert = fail # default for command line - rerun after fixing

def multipicker(idxs=None, restype=list):
    """
    pick out results - database select operation (selects columns)
    ts is indexable of indexables (list, tuple, ...) - table rows
    idxs (i,j) pairs - results is all the ts[i][j] values
        note: for namedtuples j can be either a string (field name) or an integer
    restype - result type - usually list or tuple - or a namedtuple
    """
    # should also use getval - see below
    if idxs:
        return lambda ts: restype(ts[i][j] for i,j in idxs)
    return lambda ts: restype(x for xs in ts for x in xs)

def getval(x, i):
    return getattr(x,i) if isinstance(i, str) else x[i]
    
def picker(sel, restype=list):
    return (lambda x: restype(getval(x, i) for i in sel)) if sel else restype

def picker2(sel, restype=list):
    return (lambda x: restype(x[i] for i in sel)) if sel else restype
    
    
class TMap(dict):
    temptableno = 0
    def __init__(self, args):
        dict.__init__(self, *args)
        return
        
    def join(self, t, kc, sel=None):
        " join with Table on a single field"
        kf = lambda x: x[kc]
        mpf = multipicker(sel)
        gen = (mpf((self[kv], jv)) for kv, g in itertools.groupby(sorted(t.data, kf), kf) for jv in g)
        data = itertools.chain([mpf((self.header, t.header))], gen)
        TMap.temptableno += 1
        mytype = Table.newtype('_T'+str(TMap.temptableno).zfill(6), self.tt._fields+t.tt._fields)
        return Table(mytype, data)
        
class Table:
    "Table data type for data collections"
    tt = {}
    
    def __init__(self, clsname, data, headers=None, prefix=[], selector=None):
        """
        read data into a Table object
        data is a sequence/iterable of field iterables, e.g. list or tuple
        they must be the right length - they are put into named tuples type from Table.tt
        
        """
        assert clsname in Table.tt
        self.tt = tt = Table.tt[clsname]
        self.prefix = prefix
        pf = picker(selector)
        def mkrow(x):
            "make a table row based on a picker(selector)"
            return tt(*pf(x))
        self.header = mkrow(headers) if headers else None
        self.data = [mkrow(x) for x in data] # may need to pad rows ... and filter?
        return
    
    @staticmethod
    def newtype(clsname, fields):
        assert clsname not in Table.tt
        Table.tt[clsname] = collections.namedtuple(clsname, fields)
        return clsname
    
    @staticmethod    
    def csvtype(fn, clsname, hdridx=1, hdrmap=None):
        "workout the field names for the table by reading the header in the CSV file"
        # hdrmap lets us check and name selected field
        assert os.path.isfile(fn), "needs file {} in folder {}".format(fn, os.getcwd())
        assert hdridx>0, 'hdridx ({}) must be a positive number.'.format(hdridx)
        with open(fn) as srcfd:
            for i in range(hdridx-1):
                next(srcfd) # skip header lines
            hdr = next(csv.reader(srcfd))
        fields = hdr
        if hdrmap:
            # allow position or names - may work when column names vary but positions are fixed?
            d = dict((hdr[x[0]] if isinstance(x[0], int) else x[0], x[-1]) for x in hdrmap)
            targets = d.values()
            fields = [d.get(k, k) for k in hdr]
            missing = [k for k in targets if k not in fields]
            assert not missing, "Fields missing from Table type "+clsname+': '+' '.join(missing)
        return Table.newtype(clsname, fields)
    
    def csvwrite(self, fn):
        with open(fn, "wt", newline='') as dstfd:
            for r in self.prefix:
                dstfd.write(r)
            wx = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            wx.writerow(self.header)
            wx.writerows(self.data)
        return
    
    def keyidx(self, key):
        "get int index for a field - check for table"
        fs = self.tt._fields
        v = key
        if isinstance(key, str):
            # replace string with int
            v = fs.index(key)
            if v<0 and self.headers:
                v = self.headers.index(key)
        assert isinstance(v, int), " ".join(["key type", str(type(key)), "not supported."])
        assert 0<=v<=len(fs), " ".join(["key value", str(v), "outside range for keys (<"+str(len(fs))+")"])
        return v 
    
    def picker(self, keys=None, restype=list):
        return picker2([self.keyidx(x) for x in keys] if keys else list(range(len(self.data[0]))), restype=restype)
    
    # def getidx(self, k):
    #     "return a function that retrieves required fields"
    #     if k in tuple:
    #         idxs = tuple(map(self.getidx, k))
    #         return lambda x: tuple(x[i] for i in idxs)
    #     if k is str:
    #         kv = self.tt._fields.index(k)
    #         if kv<0:
    #             kv = self.header.index(k)
    #         assert not kv<0, 'unknown table column: '+k
            
    #     assert kv in int, 'improper index type for getidx(): '+type(k)
    #     return lambda x:x[kv]
            
    def joiner(self, kidxs):
        fk =self.getidx(kidxs)
        gen = itertools.groupby(sorted(self.data, key=fk), key=fk)
        return dict((k, tuple(g)) for k, g in gen)
    
    def makemapper(self, k):
        "build a map (unique index) from column value to the data rows"
        d = dict((t[k], t) for t in self.data)
        if len(self.data)!=len(d):
            # this is not a unique map - report error
            # finds none if rows a duplicated
            dups = [t[k] for t in self.data if d[t[k]]!=t]
            assert not dups, "multiple values for keys: "+' '.join(dups)
        return d
    
    def makeindex(self, k):
        "build an index for a column value - key uniqueness isn't expected"
        kf = lambda t:t[k]
        return dict((k, tuple(g)) for k, g in itertools.groupby(sorted(self.data, kf), kf))
            
# def join():
#     """"relational join on Tables - td1 & 2 are table join descriptors
#     a table join descriptor is a Table, a list of join indices, and a list of project indices
#     """
#     js = [t.joiner(kidx) for t, kidx, vidx in tds[1:]]
#     t1, ks1, vs1 = tds
#     kp1, vp1 = map(t1.getidx(fidx) for fidx in (ks1, vs1))
#     lambda k, t, js: tuple(v for )
#     gen = ((tuple(v for k, tx1 in ((kp1(x), x) for x in t1.data) for tx in ())))
 
class CSVTable(Table):
    "Table with a filename - constructor reads the file"
    def __init__(self, clsname, filename, hdridx=1, fields=None, selector=None):
        assert hdridx>=0, "hdridx must be non-negative (>=0)"
        self.filename = filename
        with open(filename, errors="ignore") as src:
            prefix = [x for i, x in zip(range(hdridx-1), src)] # read initial lines
            csvrdr = csv.reader(src)
            hdrrow = next(csvrdr) if hdridx else None
            if clsname not in Table.tt:
                assert fields or hdrrow
                Table.newtype(clsname, fields if fields else hdrrow)
            filtfunc = lambda x: len(x)==len(hdrrow) # could do padding?
            Table.__init__(self, clsname, data=filter(filtfunc, csvrdr), headers=hdrrow, prefix=prefix, selector=selector)
        return

class EchoSurvey(Table):
    "Echo Survey file"
    def __init__(self, filename):
        name = 'EchoSurvey'
        if name not in Table.tt:
            # Echo Plate Survey file headers
            # Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Survey Fluid Height,Survey Fluid Volume,Current Fluid Volume,Fluid Composition,Fluid Units,Fluid Type,Survey Status
            clsname = self.newtype(name, 'srcname srcbc srctype srcwell height volume currvolume composition units ftype status'.split())
        with open(filename, errors="ignore") as src:
            lx = list(src)
        hdrlen = next(n for n, line in enumerate(lx, start=1) if line.startswith('[DETAILS]'))
        self.tail = lx[-5:]
        csvrdr = csv.reader(lx[hdrlen:-5])
        Table.__init__(self, clsname, data=csvrdr, header=next(csvrdr), prefix=lx[:hdrlen])
        return
    
class SourcePlates(dict):
    deadvol = { '384PP_AQ_BP': 50, '6RES_AQ_BP2': 700 } # Echo dead volume for plate types 
    def __init__(self, pairs):
        "add contents and build contents dictionary"
        for table, contents in pairs:
            wdict = dict((r.srcwell, [r.srcbc, r.srctype, r.srcwell, r.volume]) for r in table.data)
            def kf(x):
                return x[0]
            for c, g in itertools.groupby(sorted(contents), key=kf):
                self[c]= [wdict[w] for w in g]
        return
    def source(self, cx, vol):
        ""
        while self[cx][-1] < self.deadvol[self[cx][1]]+vol:
            self[cx].shift() # remove depleted wells
            if not self[cx]:
                fail("run out of "+cx)
        self[cx][-1] -= vol
        return self[cx][:-1]+(vol,)        
    
class PicklistSrc:
    "read a survey & contents file (echovolume.py output)"
    deadvol = { '384PP_AQ_BP': 50, '6RES_AQ_BP2': 700 } # Echo dead volume for plate types 

    def __init__(self, fn, idx=0):
        with open(fn) as srcfd:
            src = csv.reader(srcfd, dialect="unix")
            self.hdr = next(src)
            def voldata(xs):
                "last element of each list is an integer"
                # volumes in nanolitres
                v = [x[:-1]+[int(float(x[-1])*1000)] for x in xs] # contents in nanolitres
                return v
            self.data = dict((k, voldata(gs)) for k, gs in itertools.groupby(sorted(src, key=lambda x:x[idx]), key=lambda x:x[idx]))
        return
    
    def xfersrc(self, l, vol, depleted):
        "transfer data for liquid l - reduces volume"
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

def padwell(w):
    "convert a Nimbus well ID to an EP plate (proper/sortable) well ID"
    return '0'.join([w[0], w[1]]) if len(w)==2 else w

def unpadwell(w):
    "Echo software doesn't like well IDs like A01, it wants A1"
    return w[0]+w[2] if w[1]=='0' else w
        
def getMouseAssays(barcodes):
    """
    Read cached 96-well ear-punch plate descriptors

    Parameters
    ----------
    args : argparse namespace
        commandline arguments.

    Returns
    -------
    list of DNA wells, with mouse and assays list
        DESCRIPTION.

    """
       
    # read the ear-punch plate descriptions for all the relevant plates
    import nimbus
    def loadepbc(bc):
        "load the JSON file for an EP plate with barcode=bc"
        fn = nimbus.eppfn(bc) # consistent filename format
        with open(fn) as src:
            res = json.loads(src.read())
        assert len(res)==1
        rx = res[0]['wells'] # res is a list with one item!
        return ((r['wellLocation'], r['mouse']) for r in rx)
    
    epps = ((bc, loadepbc(bc))for bc in barcodes)    
    return epps # return a generator
    

def grouper(xs, kf=lambda x:x[0]):
    "group pairs: (x,a), (x,b), ... => (x,[a,b, ...]), ..."
    return itertools.groupby(sorted(xs, key=kf), key=kf)

def join2gen(xs, ys):
    xg, yg = grouper(*xs), grouper(*ys)
    xk, xvg = next(xg)
    xvs = list(xvg)
    print("first xk, xv =", xk, xvs)
    yk, yvg = next(yg)
    yvs = list(yvg)
    print("first yk, yv =", yk, yvs)
    xp, yp = (xk, xvs), (yk, yvs)
    
    while True:
        (xk, xvs), (yk, yvs) = xp, yp
        try:
            if xk<yk:
                xp = next(xg)
            elif xk==yk:
                # print("nextpr returns:", (xp, yp))
                yl = list(yvs)
                for xv in xvs:
                    for yv in yl:
                        yield (xv, yv)
                xp, yp = next(xg), next(yg)
            else:
                yp = next(yg)
        except StopIteration:
            break
    return
    
global ttno
ttno = 0

def joiner(tspec1, tspec2):
    "join two table"
    t1, kp1, kv1 = tspec1
    t2, kp2, kv2 = tspec2
    ts = t1, t2
    tspecs = tspec1, tspec2
    global ttno
    ttno += 1
    fx = [ x for t, fp, fv in tspecs for x in fp(t.tt._fields) ]
    newtype = Table.newtype('_TMP'+str(ttno).zfill(6), fx)
    fkv = lambda rs: (z for r, kv in zip(rs, (kv1, kv2)) for z in kv(r))
    for t, fkp, fkv in tspecs:
        print("in joiner for", type(t.data[0]))
        print("t.dict", t.__dict__.keys())
        if t.header:
            hx = fkv(t.header) 
            print("  hx =", hx)
    hx = [ fkv(t.header) for tx in tspecs for t, fkp, fkv in tx ] if all(bool(t.header) for t in ts) else None
    data = list(join2gen((t1.data, kp1), (t2.data, kp2)))
    return Table(newtype, map(fkv, data, header=hx))

def filler(rs):
    "fill blank fields with the value from the previous field - assumes fixed length records"
    rp = itertools.repeat(None)
    for r in rs:
        rp = [v if v else vp for v, vp in zip(r, rp)]
        yield rp
    return

def groupfile(fn):
    "The group file contains two mappings - family to primer, primer to reference"
    with open(fn, errors="ignore") as srcfd:
        rx = csv.reader(srcfd)
        data = [r for r in rx]
    gen1 = filler(r[:2] for r in data if r[1])
    gen2 = filler(r[1:3] for r in data)
    return gen1, gen2

def fileGetCheck(fids, fmt):
    
    dups = [fids[i] for i in range(1, len(fids)) if fids[i] in fids[:i]]
    fx = [fids[i] for i in range(len(fids)) if fids[i] not in fids[:i]]
    if dups:
        print("Duplicate DNA plate IDs ignored:", ' '.join(dups), file=sys.stderr)    
    
    globs = dict((pid, sorted(glob.glob(fmt.format(pid)))) for pid in fx)
    nofile = [fid for fid, fns in globs.items() if not fns]
    if nofile:
        alert('\n'.join(fmt.format(fid)+': no file found.' for fid in nofile))
    return collections.OrderedDict((fid, sorted(ps)[-1]) for fid, ps in globs.items())

def mk_picklist(fndst, rows):
    "output an Echo picklist given the rows (transfer spec)"
    plhdr = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Destination Plate Name,Destination Plate Barcode,Destination Plate Type,Destination Well,Volume".split(',')
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

def mytaq(wellCount, voltaq, volh2o, plateType='6RES_AQ_BP2'):
    """ 
    ***DEPRECATED***
    Create enough source wells for a Mytaq & H2O picklist, uses the same well until it's empty.
    Broken - needs a consistent output with mytaq2
    """
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


def mytaq2(wellCount, voltaq, volh2o, plateType='6RES_AQ_BP2', plate_barcodes=None):
    """ declare source wells for a Mytaq & H2O picklist, cycles through each well sequentially """
    water_wells = ['A'+c for c in '123']
    taq_wells = ['B'+c for c in '123']
    #print("water wells", water_wells, file=sys.stderr)
    # calculate whether we have enough wells, volumes in nanolitres
    dv = PicklistSrc.deadvol[plateType]
    # initial well volume is 2800uL per well in a full 6 well plate
    transfers_per_well = [(2800-dv)*1000//v for v in (voltaq, volh2o)]
    max_transfers_per_well = min(transfers_per_well)  # worst case scenario
    taq_wells_required, water_wells_required = [(wellCount+c-1)//c for c in transfers_per_well]
    plates_required = (max([taq_wells_required, water_wells_required])//3)+1
    #print('Plates required:',plates_required, file=sys.stderr)
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
    #print(set(tw), set(ww), file=sys.stderr)
#    return tw, ww


def mk_mytaq_picklist(fn, wells, taq_plate_barcodes, voltaq, volh2o, plateType='6RES_AQ_BP2'):
    "create a picklist for Mytaq & H2O transfers"
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


def i7i5alloc(vol, wellcount):
    "allocate vol to wellcount wells with unique barcode pairs"
    global args
    
    typebc = Table.newtype('BCRecord', "name platebc type well set barcode xxx oligo volume")
    # typebc = Table.newtype('BCRecord', "well name index indexName oligo volume")
    tab = CSVTable(typebc, args.i7i5)

    dv=PicklistSrc.deadvol['384PP_AQ_BP']
    def getvol(x): return x[-1]
    i7s, i5s = (sorted(((r.barcode, r.set, r.well, float(r.volume)) for r in tab.data if x in r.set), key=getvol, reverse=True) for x in ['_i7F_', '_i5R_'])
    if args.verbose or not (len(i7s) and len(i5s)):
        print("Barcode file:", tab.filename)
        print("   no. of i7Fs =", len(i7s))
        print("   no. of i5Rs =", len(i5s))
        print()
    assert len(i5s)<=len(i7s)
    assert len(i5s)
    # check barcode sets and well sets are unique
    assert all(len(frozenset(x[0] for x in xs))==len(xs) for xs in (i7s, i5s))
    assert all(len(frozenset(x[2] for x in xs))==len(xs) for xs in (i7s, i5s))
    
    i7len, i5len = (len(x) for x in (i7s, i5s))
    assert i5len<=i7len
    # create lists of transfer counts from each well
    i7cnt, i5cnt = ([int((r[-1])*1000-dv)//vol for r in ws] for ws in (i7s, i5s))
    assert wellcount<=sum(i5cnt)
    # number of imbalanced wells
    imb = [x-i5cnt[-1] for x in i5cnt]
    extra = sum(imb)
    excess = max(extra+i5len-wellcount, 0)
    while excess:
        for i in range(len(imb)):
            if imb[i]:
                imb[i] -= 1
                excess -= 1
                if not excess:
                    break
    extra = sum(imb)
    base = (wellcount-extra)//i5len
    alloc = [base+i for i in imb]
    
    # biggest i5 allocation is fewer than the number of i7 wells
    # if so, there should be no duplicate i7F, i5R pairs - we check below anyway
    assert alloc[0]<=i7len 
    assert sum(alloc)==wellcount
    i7gen = (x[:3] for xs in itertools.repeat(i7s) for x in xs)
    i5gen = (x[:3] for cnt, xs in zip(alloc, i5s) for x in [xs]*cnt)
    res = list(zip(i7gen, i5gen))
    resset = frozenset(tuple(x[0] for x in xs) for xs in res)
    assert len(resset)==wellcount # all i7-i5 pairs are unique
    return res
    
def main():
    """
    Main program. Interpret command line.
    Stage 1: Input DNA plate descriptions (Nimbus output)
    Stage 2: Check Assays, output PCR1 Echo picklist
    Stage 3: Attach sample barcodes to all wells for PCR2 stage, output picklist.

    Returns
    -------
    None.

    """
    global args
    # bcdef = "i7i5_plate_layout_20*.csv"
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
    grp3 = parser.add_argument_group('Stage 3 - Miseq barcodes')
    grp2.add_argument('-p', '--prim', help='primer/assay plate content filename')
    grp2.add_argument('-d', '--pcr', nargs='+',  help='PCR (destination) plate barcodes')
    # grp2.add_argument('-r', '--resume', help='i7i5 plate last well used - for example D15,H19.')
    grp2.add_argument('dna', nargs="+", help='barcode for 1 or more 384-well sample source (DNA - Nimbus output) plates')
    parser.add_argument('-D', '--dnadiff', action="store_true", help='allow different barcode in DNA files')
    grp3.add_argument('-i', '--i7i5', help='i7i5 plate survey filename.')
    grp3.add_argument('-x', '--xname', default=defxname, help="experiment name - used in MiSeq file (default={})".format(defxname))
    args = parser.parse_args()
    
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
            alert("target directory {} not found.\n".format(args.workdir))
            return args.workdir
    
    if args.workdir not in ['', '.']:
        os.chdir(getdir(args.workdir))
    
    
    ngid = os.path.basename(os.getcwd()) # should be 8 characters -YYYYMMDD
    pcrfmt = "PCR{{}}-picklist-{}.csv".format(ngid)
    
    fnstage2, fnstage3 = "Stage2.csv", "Stage3.csv"
    if args.pcr:
        # read Nimbus output files - record of DNA samples in 384 well plates
        # nimcol = "RecordId","TRackBC","TLabwareId","TPositionId","SRackBC","SLabwareId","SPositionId"
        dnafns = fileGetCheck(args.dna, "Echo_384_COC_0001_{0}_?.csv")
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
            alert(message)
        
        dnadict = dict(((x.srcplate, x.srcwell), (x.dstplate, x.dstwell)) for nt in nimbusTables for x in nt.data)
        # find and read Stage1 files - see also primercheck.py and its use in cgi-nimbus2.py
        #s1fns = ['Stage1-P{}.csv'.format(dnabc) for dnabc in args.dna]
        #s1type = Table.csvtype(s1fns[0], 'S1Rec')
        dnas = [CSVTable('S1Rec', 'Stage1-P{}.csv'.format(dnabc)) for dnabc in args.dna]
    
        # read primer plate info - including Echo survey 
        primerTable = CSVTable("PPRec", args.prim, fields="spn spbc spt well primer volume")
        if args.verbose:
            print("Loaded primer plate file:", args.prim)
        primset = sorted(frozenset(x.primer for x in primerTable.data))
        pfdict = dict((k, list(g)) for k, g in itertools.groupby(primset, key=lambda x:x.split('_',1)[0]))

        # create recored for PCR wells  
        wgenflds =  dnas[0].tt._fields + ('dnaplate', 'dnawell') + ('primer',)      
        wgen = [xs+dnadict[(xs.EPplate, xs.EPwell)]+(x,) for xss in dnas for xs in xss.data for f in xs.assayFamilies.split(';') if f in pfdict for x in pfdict[f]]
        
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
                print('Cannot find', r.primer, 'in known primers', file=sys.stderr)
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
        
    elif args.i7i5:
        # read Stage2 csv file
        s2tab = CSVTable('Stage2', fnstage2)           
            
        bcvol = 175 # nanolitres
        bcalloc = i7i5alloc(bcvol, len(s2tab.data))
        
        def mks3rec(s2rec, bcrec):
            return 
            
        s3flds = s2tab.tt._fields+('i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well')
        S3Rec = Table.newtype('S3Rec', s3flds)
        s3tab = Table(S3Rec, ([x for xs in (p1, p2[:3], p3[:3]) for x in xs] for p1, (p2, p3) in zip(s2tab.data, bcalloc)), headers=s3flds) 
        # output Stage 3 CSV file - used in Stage 3 below
        s3tab.csvwrite(fnstage3)                            
        if args.verbose:
            print("created:", fnstage3) 
            
        # write Barcode Picklists
        # make the destination dictionary for pcrplates
        ddict = {}
        for r in s3tab.data:
            if r.pcrplate not in ddict:
                ddict[r.pcrplate] = "Destination[{}]".format(len(ddict)+1)
        # plate barcode for i7i5 plate is unknown
        rowi7 = (("Source[1]", "", '384PP_AQ_BP', r.i7well, ddict[r.pcrplate], r.pcrplate, '384PP_AQ_BP', r.pcrwell, bcvol) for r in s3tab.data)
        rowi5 = (("Source[1]", "", '384PP_AQ_BP', r.i5well, ddict[r.pcrplate], r.pcrplate, '384PP_AQ_BP', r.pcrwell, bcvol) for r in s3tab.data)
        fnbc = pcrfmt.format("2bc")
        mk_picklist(fnbc, (r for rs in (rowi7, rowi5) for r in rs))
        if args.verbose:
            print("created:", fnbc) 
    
        
        # also PCR2 water and Taq
        fndst = pcrfmt.format("2TaqWater")
        mk_mytaq_picklist(fndst, s3tab.data, args.taq, 2000, 650)
        if args.verbose:
            print('created:', fndst)
            
        # MiSeq file
        fnmiseq = fmtmiseq.format(ngid)
        with open(fnmiseq, "wt", newline='') as dstfd:
            dstfd.write("""[Header],,,,,,,,,,
IEMFileVersion,4,,,,,,,,,
Investigator Name,,,,,,,,,,
Experiment Name,{xname},,,,,,,,,
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
""".format(xname=args.xname))
            hdr = "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description".split(',')
            dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            dst.writerow(hdr)
            gen =((r.pcrplate+'_'+padwell(r.pcrwell), '', '', '',
                   r.i7name, r.i7bc,
                   r.i5name, r.i5bc,
                   'NGStest') for r in s3tab.data)
            dst.writerows(gen)
        if args.verbose:
            print('created:', fnmiseq)

        import stage2report as s2r
        fns3 = "Stage3.html"
        s2r.output(fns3, s3tab.tt._fields, s3tab.data,custom=args.custom)    
        if args.verbose:
            print('created:', fns3)
    
    return   

if __name__ == '__main__':
    main()
