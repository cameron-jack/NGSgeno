#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Jan 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: Adjusted paths for app folder
@last_edit: 2022-05-05
@edit_comment: 

A collection of code functions and definitions that are used throughout the pipeline
"""

import os
import sys
import glob
import collections
import traceback
from pathlib import PurePath
from collections import defaultdict
import csv

import bin.file_io as file_io 
from bin.file_io import GUARD_TYPES
import itertools

### Set up SSL/TLS certificate
try:
    from OpenSSL import crypto, SSL
except:
    pass

def cert_gen(emailAddress="emailAddress", commonName="commonName", countryName="NT",
        localityName="localityName", stateOrProvinceName="stateOrProvinceName",
        organizationName="organizationName", organizationUnitName="organizationUnitName",
        serialNumber=0, validityStartInSeconds=0, validityEndInSeconds=10*365*24*60*60,
        KEY_FILE = "private.key", CERT_FILE="selfsigned.crt"):
    """
    Generate an SSL/TLS certificate if we don't already have one. Very important as we send passwords around.
    can look at generated file using openssl:
        openssl x509 -inform pem -in selfsigned.crt -noout -text
    """
   
    # create a key pair
    k = crypto.PKey()
    k.generate_key(crypto.TYPE_RSA, 4096)
    # create a self-signed cert
    cert = crypto.X509()
    cert.get_subject().C = countryName
    cert.get_subject().ST = stateOrProvinceName
    cert.get_subject().L = localityName
    cert.get_subject().O = organizationName
    cert.get_subject().OU = organizationUnitName
    cert.get_subject().CN = commonName
    cert.get_subject().emailAddress = emailAddress
    cert.set_serial_number(serialNumber)
    cert.gmtime_adj_notBefore(0)
    cert.gmtime_adj_notAfter(validityEndInSeconds)
    cert.set_issuer(cert.get_subject())
    cert.set_pubkey(k)
    cert.sign(k, 'sha512')
    with open(CERT_FILE, "wt") as f:
        f.write(crypto.dump_certificate(crypto.FILETYPE_PEM, cert).decode("utf-8"))
    with open(KEY_FILE, "wt") as f:
        f.write(crypto.dump_privatekey(crypto.FILETYPE_PEM, k).decode("utf-8"))


# Global
CONFIG_PATH = PurePath('library/config.ini')  # everything is now relative to the app directory
ERROR_FN = 'error_msg.txt'

# default transfer volumes for liquid handling. These are loaded into each experiment upon creation.
DNA_VOL = 200
PRIMER_VOL = 500
PRIMER_TAQ_VOL = 1000
PRIMER_WATER_VOL = 300
INDEX_VOL = 175
INDEX_TAQ_VOL = 2000
INDEX_WATER_VOL = 650
DEAD_VOLS = {'384PP_AQ_BP': 2.5*1000, '6RES_AQ_BP2': 250*1000, 'Hard Shell 384 well PCR Biorad': 2.5*1000}
CAP_VOLS = {'384PP_AQ_BP': 12*1000, '6RES_AQ_BP2': 2800*1000, 'Hard Shell 384 well PCR Biorad': 12*1000}
# convenient name to exact plate model mapping
PLATE_TYPES = {'Echo6': '6RES_AQ_BP2', 'Echo384': '384PP_AQ_BP',
               'PCR384': 'Hard Shell 384 well PCR Biorad'}

    
def output_error(exc, msg=''):
    with open(ERROR_FN, 'wt') as f:
        if msg:
            f.write(f"{exc} {msg}")
        f.write(traceback.format_exc())

row_ordered_96 = []
for r in [chr(65+i) for i in range(8)]:
    for c in [str(j+1) for j in range(12)]:
        row_ordered_96.append(r+c)

row_ordered_384 = []
for r in [chr(65+i) for i in range(16)]:
    for c in [str(j+1) for j in range(24)]:
        row_ordered_384.append(r+c)

col_ordered_96 = []
for c in [str(j+1) for j in range(12)]:
    for r in [chr(65+i) for i in range(8)]:
        col_ordered_96.append(r+c)

col_ordered_384 = []
for c in [str(j+1) for j in range(24)]:
    for r in [chr(65+i) for i in range(16)]:
        col_ordered_384.append(r+c)

nimbus_ordered_96 = []
for c in [str(j+1) for j in range(12)]:
    for r in 'ACEGBDFH':
        nimbus_ordered_96.append(r+c)


def calc_plate_assay_usage(location_sample: dict, denied_assays: list=[], denied_wells: list=[]) -> defaultdict:
    """ 
    Required input: Experiment.plate_location_sample[plateid] dictionary of locations and samples for a single plate
            "wells" is a required field.
    Optional inputs: lists of denied assays and denied wells, which will be ignored in the calculation
    Output: a dictionary of assays and their usage counts
    """
    assay_counts = defaultdict(int)
    for well in location_sample['wells']:
        if well not in denied_wells:
            for assay in location_sample[well]['assays']:
                if assay not in denied_assays:
                    assay_counts[assay] += 1
    return assay_counts
 

# File choosers
def get_mouse_ref():
    return sorted(glob.glob(os.path.join('..','library', "reference_sequences_*.fa")) +        
            glob.glob(os.path.join('..','library',"reference_sequences_*.txt")), reverse=True)[0]

def get_mouse_assaylist():
    return sorted(glob.glob(os.path.join('..','library', 'assay_list_*.csv')), reverse=True)[0]

def get_mouse_conversions():
    return sorted(glob.glob(os.path.join('..','library', 'NGS_assay_conversions_*.xlsx')), reverse=True)[0]

def get_mouse_primer_layout():
    return sorted(glob.glob(os.path.join('..','library', 'primer_layout*_*.csv')), reverse=True)[0]

def get_mouse_i7i5_layout():
    return sorted(glob.glob(os.path.join('..','library', 'i7i5_plate_layout_*.csv')), reverse=True)[0]

def get_mouse_primer_survey():
    return sorted(glob.glob(os.path.join('..','library',"*_platesurvey_fluid_volume_*primerplate.csv")), reverse=True)[0]

def get_mouse_i7i5_survey():
    return sorted(glob.glob(os.path.join('..','library',"*_platesurvey_fluid_volume_I7I5PLATE.csv")), reverse=True)[0]

# I/O for template_files.txt which stores the current pipeline configuration
# DEPRECATED
#class TemplateFiles():
#    """
#    Tracks information on the following file paths (when available)
#    mouse pipeline files:
#    mouse_ref - defaults to library .fa or.txt (.fa takes precedence)
#    assay_list - defaults to library
#    conversions - defaults to library (this file may on borrowed time thanks to Rodentity)

#    custom pipeline files:

#    """
#    def __init__(self, debug=False):
#        """ initialise with values from file, if it exists """
#        self.files = collections.defaultdict(str)
#        if os.path.exists('template_files.txt'):
#            with open('template_files.txt', 'rt') as f:
#                for line in f:
#                    cols = line.strip().split('\t')
#                    self.files[cols[0]] = cols[1]
#        if debug:
#            print('TemplateFiles:', [(k, self.files[k]) for k in self.files], file=sys.stdout)
#        self.cust_format = None # guard type for custom barcodes. Must be one of file_io.GUARD_TYPES

#    def save_config(self):
#        """ save configuration to template_files.txt """
#        with open('template_files.txt', 'wt') as outf:
#            for key in self.files:
#                print(key + '\t' + self.files[key], file=outf)
 
#    def load_config(self, debug=False):
#        """ load existing template file info if available """
#        if os.path.exists('template_files.txt'):
#            with open('template_files.txt', 'rt') as f:
#                for line in f:
#                    cols = line.strip().split('\t')
#                    self.files[cols[0]] = cols[1]
#        else:
#            print("No file found: 'template_files.txt'", file=sys.stdout)
#            return
#        if debug:
#            print('TemplateFiles:', [(k, self.files[k]) for k in self.files], file=sys.stdout)

#    def set_mouse_files(self, ref=None, assay_list=None, conversions=None, save=False):
#        """ set files specific to the mouse pipeline. If arguments are None, use defaults in library folder.
#            variables are prefixed with mouse_ and operate on m or M guarded barcodes """
#        # mouse reference sequence file
#        self.files['mouse_ref'] = get_mouse_ref()
#        if ref:
#            if os.path.exists(ref):
#                self.files['mouse_ref'] = primer_survey
#            else:
#                print('Mouse reference sequences file', ref, 'not found! Using default:', 
#                        self.files['mouse_ref'], file=sys.stdout)    
#        # mouse valid assay list
#        self.files['mouse_assay_list'] = get_mouse_assaylist()
#        if assay_list:
#            if os.path.exists(assay_list):
#                self.files['mouse_assay_list'] = assay_list
#            else:
#                print('Assay list file', assay_list, 'not found! Using default:', self.files['mouse_assay_list'], file=sys.stdout)
#        self.files['mouse_conversions'] = get_mouse_conversions()
#        # mouse assay old/new conversions
#        if conversions:
#            if os.path.exists(conversions):
#                self.files['mouse_conversions'] = conversions
#            else:
#                print('Conversions file', conversions, 'not found! Using default:', self.files['mouse_conversions'], file=sys.stdout)
#        if save:
#            self.save_config()

#    def set_custom_files(self, manifest=None, assays=None, ref=None, save=False):
#        """ set files specific to the custom pipeline, variables have prefix cust_ and work with 'custom' guarded prefices """
#        if manifest:
#            if os.path.exists(manifest):
#                self.files['cust_manifest'] = manifest
#            else:
#                print('Custom manifest file', manifest, 'not found! No custom manifest set', file=sys.stdout)
#        if assays:
#            if os.path.exists(assays):
#                self.files['cust_assay_list'] = assays
#            else:
#                print('Custom assay list file', assays, 'not found! No custom assay list set', file=sys.stdout)
#        if ref:
#            if os.path.exists(ref):
#                self.files['cust_ref'] = cust_ref
#            else:
#                print('Custom reference sequence file', ref, 'not found! No custom reference set', file=sys.stdout)   
#        if save:
#            self.save_config()

#    def set_custom_file_format(self, format):
#        """ sets the default interpretation of unguarded barcodes """
#        if format not in GUARD_TYPES:
#            print('Guard type not recognised for custom manifest format', format,
#                    'defaulting to',self.cust_format, file=sys.stdout)
#        else:
#            self.cust_format = format
    
#    def set_layout_files(self, primer_plate=None, i7i5_plate=None, save=False):
#        """ set layout files for Echo robot processing - primers and barcodes which are shared across runs """
#        self.files['primer_plate'] = get_mouse_primer_layout()
#        if primer_plate:
#            if os.path.exists(primer_plate):
#                self.files['primer_plate'] = primer_plate
#            else:
#                print('Primer plate layout', primer_plate,'not found! Using default:', self.files['primer_plate'], file=sys.stdout)
#        self.files['i7i5_plate'] = get_mouse_i7i5_layout()
#        if i7i5_plate:
#            if os.path.exists(i7i5_plate):
#                self.files['i7i5_plate'] = i7i5_plate
#            else:
#                print('i7i5 barcode plate layout', i7i5_plate,'not found! Using default:', self.files['i7i5_plate'], file=sys.stdout)
#        if save:
#            self.save_config()

#    def set_survey_files(self, primer_survey=None, i7i5_survey=None, save=False):
#        """ set plate survey files for Echo robot processing - primers and barcodes which are shared across runs """
#        self.files['primer_survey'] = get_mouse_primer_survey()
#        if primer_survey:
#            if os.path.exists(primer_survey):
#                self.files['primer_survey'] = primer_survey
#            else:
#                print('Primer survey file', primer_survey, 'not found! Using default:', self.files['primer_survey'], file=sys.stdout)
#        self.files['i7i5_survey'] = get_mouse_i7i5_survey()
#        if i7i5_survey:
#            if os.path.exists(i7i5_survey):
#                self.files['i7i5_survey'] = i7i5_survey
#            else:
#                print('i7i5 barcode survey file', i7i5_survey, 'not found! Using default:', self.files['i7i5_survey'], file=sys.stdout)
#        if save:
#            self.save_config()
         

# Item/attribute selection macros - used for Table selections
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
    
def picker_attr(sel, restype=list):
    return (lambda x: restype(getval(x, i) for i in sel)) if sel else restype

def picker_index(sel, restype=list):
    return (lambda x: restype(x[i] for i in sel)) if sel else restype

def getval(x, i):
    return getattr(x,i) if isinstance(i, str) else x[i]

# Add/remove leading 0's from well ids
def padwell(w):
    "convert a Nimbus well ID to an EP plate (proper/sortable) well ID"
    return '0'.join([w[0], w[1]]) if len(w)==2 else w

def unpadwell(w):
    "Echo software doesn't like well IDs like A01, it wants A1"
    if len(w) == 2:
        return w
    try:
        well = w[0]+w[2] if w[1]=='0' else w
    except:
        print(f"unpadwell failed with {w=}")
    return well 

### Table classes to act like a relational DB

# Never used?
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
        pf = picker_attr(selector)
        def mkrow(x):
            "make a table row based on a picker_attr(selector)"
            return tt(*pf(x))
        self.header = mkrow(headers) if headers else None
        self.data = [mkrow(x) for x in data] # may need to pad rows ... and filter?
        return
    
    @staticmethod
    def newtype(clsname, fields):
        #assert clsname not in Table.tt
        if clsname in Table.tt:
            return clsname
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
    
    def csvwrite(self, fn, output_plate_guards=False):
        """ 
        Check which columns contain plate barcodes, and enforce plate guards or their removal based on output_plate_guards
        """
        with open(fn, "wt", newline='') as dstfd:
            for r in self.prefix:
                dstfd.write(r)
            wx = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
            data = []
            for row in self.data:
                if output_plate_guards:  # guard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() else file_io.guard_pbc(str(d), silent=True) for (h,d) in zip(self.header, row)])
                else:  # unguard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() else file_io.unguard_pbc(str(d), silent=True) for (h,d) in zip(self.header, row)])
            wx.writerow(self.header)
            wx.writerows(data)
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
        return picker_index([self.keyidx(x) for x in keys] if keys else list(range(len(self.data[0]))), restype=restype)
    
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
            
 
class CSVTable(Table):
    "Table with a filename - constructor reads the file"
    def __init__(self, clsname, filename, hdridx=1, fields=None, selector=None):
        assert hdridx>=0, "hdridx must be non-negative (>=0)"
        self.filename = filename
        with open(filename, errors="ignore", newline='') as src:
            prefix = [x for i, x in zip(range(hdridx-1), src)] # read initial lines
            #print(f"{prefix=}")
            csvrdr = csv.reader(src)
            hdrrow = next(csvrdr) if hdridx else None
            #print(f"{hdrrow=}")
            if clsname not in Table.tt:
                assert fields or hdrrow
                Table.newtype(clsname, fields if fields else hdrrow)
            filtfunc = lambda x: len(x)==len(hdrrow) # could do padding?
            # TODO: need a filter-type lambda to apply strip() to everything... covered for now but it's not a safe patch
            Table.__init__(self, clsname, data=filter(filtfunc, csvrdr), headers=hdrrow, prefix=prefix, selector=selector)
        return


## Apparently never used
#class EchoSurvey(Table):
#    "Echo Survey file"
#    def __init__(self, filename):
#        name = 'EchoSurvey'
#        if name not in Table.tt:
#            # Echo Plate Survey file headers
#            # Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Survey Fluid Height,Survey Fluid Volume,Current Fluid Volume,Fluid Composition,Fluid Units,Fluid Type,Survey Status
#            clsname = self.newtype(name, 'srcname srcbc srctype srcwell height volume currvolume composition units ftype status'.split())
#        with open(filename, errors="ignore") as src:
#            lx = list(src)
#        hdrlen = next(n for n, line in enumerate(lx, start=1) if line.startswith('[DETAILS]'))
#        self.tail = lx[-5:]
#        csvrdr = csv.reader(lx[hdrlen:-5])
#        Table.__init__(self, clsname, data=csvrdr, header=next(csvrdr), prefix=lx[:hdrlen])
#        return

## Apparently never used
#class SourcePlates(dict):
#    deadvol = { '384PP_AQ_BP': 50, '6RES_AQ_BP2': 700 } # Echo dead volume for plate types 
#    def __init__(self, pairs):
#        "add contents and build contents dictionary"
#        for table, contents in pairs:
#            wdict = dict((r.srcwell.strip(), [r.srcbc.strip(), r.srctype.strip(), r.srcwell.strip(), r.volume.strip()]) for r in table.data)
#            def kf(x):
#                return x[0]
#            for c, g in itertools.groupby(sorted(contents), key=kf):
#                self[c]= [wdict[w] for w in g]
#        return
#    def source(self, cx, vol):
#        ""
#        while self[cx][-1] < self.deadvol[self[cx][1]]+vol:
#            self[cx].shift() # remove depleted wells
#            if not self[cx]:
#                print("run out of "+cx, file=sys.stdout)
#                exit(1)
#        self[cx][-1] -= vol
#        return self[cx][:-1]+(vol,)

## Apparently never used
#def joiner(tspec1, tspec2):
#    """ join two tables """
#    try:
#        t1, kp1, kv1 = tspec1
#        t2, kp2, kv2 = tspec2
#        ts = t1, t2
#        tspecs = tspec1, tspec2
#        global ttno
#        ttno += 1
#        fx = [ x for t, fp, fv in tspecs for x in fp(t.tt._fields) ]
#        newtype = Table.newtype('_TMP'+str(ttno).zfill(6), fx)
#        fkv = lambda rs: (z for r, kv in zip(rs, (kv1, kv2)) for z in kv(r))
#        for t, fkp, fkv in tspecs:
#            print("in joiner for", type(t.data[0]))
#            print("t.dict", t.__dict__.keys())
#            if t.header:
#                hx = fkv(t.header) 
#                print("  hx =", hx)
#        hx = [ fkv(t.header) for tx in tspecs for t, fkp, fkv in tx ] if all(bool(t.header) for t in ts) else None
#        data = list(join2gen((t1.data, kp1), (t2.data, kp2)))
#        return Table(newtype, map(fkv, data, header=hx))
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode.joiner')

## Apparently never used
#def join2gen(xs, ys):
#    """ a generator that joins two iterables if possible - why? """
#    try:
#        xg, yg = grouper(*xs), grouper(*ys)
#        xk, xvg = next(xg)
#        xvs = list(xvg)
#        print("first xk, xv =", xk, xvs)
#        yk, yvg = next(yg)
#        yvs = list(yvg)
#        print("first yk, yv =", yk, yvs)
#        xp, yp = (xk, xvs), (yk, yvs)
    
#        while True:
#            (xk, xvs), (yk, yvs) = xp, yp
#            try:
#                if xk<yk:
#                    xp = next(xg)
#                elif xk==yk:
#                    # print("nextpr returns:", (xp, yp))
#                    yl = list(yvs)
#                    for xv in xvs:
#                        for yv in yl:
#                            yield (xv, yv)
#                    xp, yp = next(xg), next(yg)
#                else:
#                    yp = next(yg)
#            except StopIteration:
#                break
#        return
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode.join2gen')

## Apparently never used
#def filler(rs):
#    """ fill blank fields with the value from the previous field - assumes fixed length records """
#    try:
#        rp = itertools.repeat(None)
#        for r in rs:
#            rp = [v if v else vp for v, vp in zip(r, rp)]
#            yield rp
#        return
#    except Exception as exc:
#        output_error(exc, msg='Error in echo_barcode.filler')