#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Jan 2022
@author: Cameron Jack, Bob Buckley, ANU Bioinformatics Consultancy, JCSMR, Australian National University

A collection of code functions and definitions that are used throughout the pipeline
"""

import sys
import collections
from pathlib import PurePath, Path, WindowsPath
from collections import defaultdict
import csv
import hashlib
from math import ceil
import threading
import queue
import datetime

import functools

import itertools

from stutil import m

### Create a non-blocking log writer?


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
DNA_VOL = 200  # adds up to 2000
PRIMER_VOL = 250  # adds up to 2000
PRIMER_TAQ_VOL = 1000  # adds up to 2000
PRIMER_WATER_VOL = 550  # adds up to 2000
INDEX_VOL = 175
INDEX_TAQ_VOL = 2000
INDEX_WATER_VOL = 650
DEAD_VOLS = {'384PP_AQ_BP': 30*1000, '6RES_AQ_BP2': 250*1000, 'Hard Shell 384 well PCR Biorad': 0.001*1000}
CAP_VOLS = {'384PP_AQ_BP': 60*1000, '6RES_AQ_BP2': 2800*1000, 'Hard Shell 384 well PCR Biorad': 6*1000}

# convenient name to exact plate model mapping
PLATE_TYPES = {'Echo6': '6RES_AQ_BP2', 'Echo384': '384PP_AQ_BP',
               'PCR384': 'Hard Shell 384 well PCR Biorad','96':'96'}
WATER_WELLS = ['A'+c for c in '123']
TAQ_WELLS = ['B'+c for c in '123']
# global - set of all possible guard types
GUARD_TYPES = set(['m','r','c','p', 'a'])


def usable_volume(vol, plate_type):
    return max(vol - DEAD_VOLS[PLATE_TYPES[plate_type]],0)

def num_doses(raw_vol, per_use_vol, plate_type):
    """ The number of available doses (dispenses) in a well """
    number_of_doses = usable_volume(raw_vol, plate_type) // per_use_vol
    return number_of_doses
                    
def num_req_wells(req_vol, plate_type='Echo384'):
    """ Save a bunch of boilerplate by doing the calculation for us """
    req_wells = ceil(req_vol / (CAP_VOLS[PLATE_TYPES[plate_type]] - DEAD_VOLS[PLATE_TYPES[plate_type]]))
    return req_wells

def num_req_taq_water_plates(req_taq_vol, req_water_vol):
    """ Save some boilerplate """
    per_plate_avail = 3*(CAP_VOLS[PLATE_TYPES['Echo6']] - DEAD_VOLS[PLATE_TYPES['Echo6']])
    taq_plates = ceil(req_taq_vol / per_plate_avail)
    #print(f'{per_plate_avail=}')
    #print(f'{req_taq_vol=}, {req_water_vol=}')
    water_plates = ceil(req_water_vol / per_plate_avail)
    return max([taq_plates, water_plates])

 
def get_md5(fn):
    """ return the md5sum of a file to ensure data safety """
    with open(fn, "rb") as f:
        file_bytes = f.read()  # read file as bytes
        readable_hash = hashlib.md5(file_bytes).hexdigest()
    return readable_hash

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

col_ordered_6 = []
for c in [str(j+1) for j in range(3)]:
    for r in [chr(65+i) for i in range(2)]:
        col_ordered_6.append(r+c)
    
row_ordered_6 = []
for c in [str(j+1) for j in range(2)]:
    for r in [chr(65+i) for i in range(3)]:
        row_ordered_6.append(r+c)



#def calc_plate_assay_primer_usage(exp, dna_pids=None, denied_assays: list=None, denied_wells: list=None, filtered=True, included_guards=GUARD_TYPES) -> defaultdict:
#    """  
#    Optional inputs: 
#        - dna plate pids, if not provided then all will be used
#        - lists of denied assays and denied wells, which will be ignored in the calculation
#        - included_guards is an interable of guard types (r,c,m,etc) which will be included. It defaults to all guard types.
#    Output: a dictionary of assays and their usage counts
#    """                          
#    assay_counts = defaultdict(lambda: 0)
#    primer_counts = defaultdict(lambda: 0)
#    if not dna_pids:
#        dna_pids = exp.get_dna_pids()
   
#    for pid in dna_pids:
#        for well in exp.sample_plate_location[pid]['wells']:
#            if well in denied_wells:
#                continue
#            guard = get_guard_type(exp.plate_location_sample[pid][well]['barcode'])
#            if guard in included_guards:
#                if 'assays' not in exp.plate_location_sample[pid][well]:
#                    continue
#                for assay in exp.platelocation_sample[well]['ngs_assays']:
#                    if assay not in denied_assays:
#                        assay_counts[assay] += 1
#                if not filtered:
#                    for assay in exp.plate_location_sample[pid][well]['unknown_assays']:
#                        if assay not in denied_assays:
#                            assay_counts[assay] += 1
#    for assay in assay_counts:
#        primerfam = choose_primerfam(exp, get_assay_family(assay))
#        if primerfam:
#            primer_counts[primerfam] += assay_counts[assay]

#    return assay_counts, primer_counts
 
### Helper functions for IO

# See top of file for GUARD_TYPES constant

class UnguardedBarcodeError(Exception):
    """ Exception raised for barcodes with degraded or missing guards

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class EmptyBarcodeError(Exception):
    """ Exception raised when a barcode is an empty string

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message


class ExistingGuardError(Exception):
    """ Exception raised when a barcode appears to already have a valid set of guards

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message



# Confirm guards
def is_guarded_mbc(mbc):
    """ Return True if Musterer barcode is correctly guarded, else False """
    mbc = str(mbc)
    if mbc.startswith('m') and mbc.endswith('m'):
        return True
    return False

def is_guarded_rbc(rbc):
    """ Return True if Rodentity barcode is correctly guarded (rNNNNNr), else False """
    rbc = str(rbc)
    if rbc.startswith('r') and rbc.endswith('r'):
        return True
    return False

def is_guarded_cbc(cbc):
    """ Return True if custom barcode is correctly guarded, else False """
    cbc = str(cbc)
    if cbc.startswith('c') and cbc.endswith('c'):
        return True
    return False

def is_guarded_pbc(pbc):
    """ Return True if plate barcode is correctly guarded, else False """
    pbc = str(pbc)
    if pbc.startswith('p') and pbc.endswith('p'):
        return True
    return False

def is_guarded_abc(abc):
    """ Return True if amplicon sample barcode is correctly guarded, else False """
    abc = str(abc)
    if abc.startswith('a') and abc.endswith('a'):
        return True
    return False

def is_guarded(bc):
    """ generic check of guard status """
    if is_guarded_mbc(bc):
        return True
    elif is_guarded_rbc(bc):
        return True
    elif is_guarded_cbc(bc):
        return True
    elif is_guarded_pbc(bc):
        return True
    elif is_guarded_abc(bc):
        return True
    return False

def get_guard_type(bc):
    """ return the character of the guard in a barcode """
    if is_guarded(bc):
        return bc[0]
    return None


# Add guards to barcodes
def guard_mbc(mbc, silent=False):
    """ Add guards to Musterer barcode """
    mbc = str(mbc).strip()
    if mbc == '':
        msg = f"Musterer barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(mbc):
        if silent:
            return mbc
        msg = f"Musterer barcode appears to already have Musterer guards: {mbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(mbc):
        msg = f"Musterer barcode appears to already have Rodentity guards: {mbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(mbc):
        msg = f"Musterer barcode appears to already have custom guards: {mbc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(mbc):
        msg = f"Musterer barcode appears to already have plate guards: {mbc}"
        raise ExistingGuardError(msg)
    return 'm' + mbc + 'm'

def guard_rbc(rbc, silent=False):
    """ Add guards to Rodentity barcode """
    rbc = str(rbc).strip()
    if rbc == '':
        msg = f"Rodentity barcode is zero length string!" 
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(rbc):
        msg = f"Rodentity barcode appears to already have Musterer guards: {rbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(rbc):
        if silent:
            return rbc
        msg = f"Rodentity barcode appears to already have Rodentity guards: {rbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(rbc):
        msg = f"Rodentity barcode appears to already have custom guards: {rbc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(rbc):
        msg = f"Rodentity barcode appears to already have plate guards: {rbc}"
        raise ExistingGuardError(msg)
    return 'r' + rbc + 'r'  # Rodentity barcodes already start with 'M'

def guard_cbc(cbc, silent=False):
    """ Add guards to custom barcode """
    cbc = str(cbc).strip()
    if cbc == '':
        msg = f"Custom barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(cbc):
        msg = f"Custom barcode appears to already have Musterer guards: {cbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(cbc):
        msg = f"Custom barcode appears to already have Rodentity guards: {cbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(cbc):
        if silent:
            return cbc
        msg = f"Custom barcode appears to already have custom guards: {cbc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(cbc):
        msg = f"Custom barcode appears to already have plate guards: {cbc}"
        raise ExistingGuardError(msg)
    return 'c' + cbc + 'c'

def guard_pbc(pbc, silent=False):
    """ Add guards to plate barcode, silent ignores already correctly guarded case """
    pbc = str(pbc).strip()
    if pbc == '':
        msg = f"Plate barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(pbc):
        msg = f"Plate barcode appears to already have Musterer guards: {pbc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(pbc):
        msg = f"Plate barcode appears to already have Rodentity guards: {pbc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(pbc):
        msg = f"Plate barcode appears to already have custom guards: {pbc}"
        raise ExistingGuardError(msg)
    if is_guarded_abc(pbc):
        msg = f"Plate barcode appears to already have amplicon guards: {pbc}"
    if is_guarded_pbc(pbc):
        if silent:
            return pbc
        msg = f"Plate barcode appears to already have plate guards: {pbc}"
        raise ExistingGuardError(msg)
    return 'p' + pbc + 'p'

def guard_abc(abc, silent=False):
    """ Add guards to amplicon sample barcodes, silent ignores already correctly guarded case """
    abc = str(abc).strip()
    if abc == '':
        msg = f"Amplicon sample barcode is zero length string!"
        raise EmptyBarcodeError(msg)
    if is_guarded_mbc(abc):
        msg = f"Amplicon sample barcode appears to already have Musterer guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_rbc(abc):
        msg = f"Amplicon sample barcode appears to already have Rodentity guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_cbc(abc):
        msg = f"Amplicon sample barcode appears to already have custom guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(abc):
        msg = f"Amplicon sample barcode appears to already have plate guards: {abc}"
        raise ExistingGuardError(msg)
    if is_guarded_pbc(abc):
        if silent:
            return abc
        msg = f"Amplicon sample barcode appears to already have amplicon guards: {abc}"
        raise ExistingGuardError(msg)
    return 'a' + abc + 'a'


# Unguard barcodes
def unguard_mbc(mbc, silent=False):
    """ remove guards from a Musterer barcode """
    if type(mbc) != str:
        msg = f"Musterer barcode is not a string! {mbc} is type {type(mbc)}"
        raise AttributeError(msg)
    if not mbc.startswith('m') and not mbc.endswith('m') and silent:  # just return unguarded barcodes as themselves
        return mbc
    if not mbc.startswith('m') or not mbc.endswith('m'):
        msg = f"Musterer barcode guards degraded or missing in: {mbc}"
        raise UnguardedBarcodeError(msg)
    return mbc[1:-1]

def unguard_rbc(rbc, silent=False):
    """ remove guards from a Rodentity barcode """
    if type(rbc) != str:
        msg = f"Rodentity barcode is not a string! {rbc} is type {type(rbc)}"
        raise AttributeError(msg)
    if not rbc.startswith('r') and not rbc.endswith('r') and silent:  # just return unguarded barcodes as themselves
        return rbc
    if not rbc.startswith('r') or not rbc.endswith('r'):
        msg = f"Musterer barcode guards degraded or missing in: {rbc}"
        raise UnguardedBarcodeError(msg)
    return rbc[1:-1]  # Rodentity barcodes keep their 'M' prefix

def unguard_cbc(cbc, silent=False):
    """ remove guards from a custom barcode """
    if type(cbc) != str:
        msg = f"Custom barcode is not a string! {cbc} is type {type(cbc)}"
        raise AttributeError(msg)
    if not cbc.startswith('c') and not cbc.endswith('c') and silent:  # just return unguarded barcodes as themselves
        return cbc
    if not cbc.startswith('c') or not cbc.endswith('c'):
        msg = f"Custom barcode guards degraded or missing in: {cbc}"
        raise UnguardedBarcodeError(msg)
    return cbc[1:-1]

def unguard_pbc(pbc, silent=False):
    """ remove guards from a plate barcode """
    if type(pbc) != str:
        msg = f"Plate barcode is not a string! {pbc} is type {type(pbc)}"
        raise AttributeError(msg)
    if not pbc.startswith('p') and not pbc.endswith('p') and silent:  # just return unguarded barcodes as themselves
        #print(f'{pbc=}', file=sys.stderr)
        return pbc
    if not pbc.startswith('p') or not pbc.endswith('p'):
        msg = f"Plate barcode guards degraded or missing in: {pbc}"
        raise UnguardedBarcodeError(msg)
    return pbc[1:-1]

def unguard_abc(abc, silent=False):
    """ remove guards from an amplicon sample barcode """
    if type(abc) != str:
        msg = f"Amplicon sample barcode is not a string! {abc} is type {type(abc)}"
        raise AttributeError(msg)
    if not abc.startswith('a') and not abc.endswith('a') and silent:  # just return unguarded barcodes as themselves
        return abc
    if not abc.startswith('a') or not abc.endswith('a'):
        msg = f"Amplicon sample barcode guards degraded or missing in: {abc}"
        raise UnguardedBarcodeError(msg)
    return abc[1:-1]

def unguard(bc, silent=False):
    """ shortcut, for when we don't know what kind of guarding we have but need to remove it for output """
    if is_guarded_mbc(bc):
        return unguard_mbc(bc)
    elif is_guarded_rbc(bc):
        return unguard_rbc(bc)
    elif is_guarded_cbc(bc):
        return unguard_cbc(bc)
    elif is_guarded_pbc(bc):
        return unguard_pbc(bc)
    elif is_guarded_abc(bc):
        return unguard_abc(bc)
    else:
        if not silent:
            print('Possibly unguarded already:',bc, file=sys.stdout)
        return bc


# guard/unguard a collection
def guard_pids(pids, silent=True, caller_id=None):
    success = True
    gpids = []
    try:
        gpids = [guard_pbc(p,silent=silent) for p in pids]
    except Exception as exc:
        m(f"{exc}", level='error', caller_id=caller_id)
        success = False
    return gpids, success


def unguard_pids(gpids, silent=True, caller_id=None):
    success = True
    pids = []
    try:
        pids = [unguard_pbc(p,silent=silent) for p in gpids]
    except Exception as exc:
        m(f"{exc}", level='error', caller_id=caller_id)
        success = False
    return pids, success


# File choosers - obsolete
#def get_mouse_ref():
#    return sorted(glob.glob(os.path.join('..','library', "reference_sequences_*.fa")) +        
#            glob.glob(os.path.join('..','library',"reference_sequences_*.txt")), reverse=True)[0]

#def get_mouse_assaylist():
#    return sorted(glob.glob(os.path.join('..','library', 'assay_list_*.csv')), reverse=True)[0]

#def get_mouse_conversions():
#    return sorted(glob.glob(os.path.join('..','library', 'NGS_assay_conversions_*.xlsx')), reverse=True)[0]

#def get_mouse_primer_layout():
#    return sorted(glob.glob(os.path.join('..','library', 'primer_layout*_*.csv')), reverse=True)[0]

#def get_mouse_i7i5_layout():
#    return sorted(glob.glob(os.path.join('..','library', 'i7i5_plate_layout_*.csv')), reverse=True)[0]

#def get_mouse_primer_survey():
#    return sorted(glob.glob(os.path.join('..','library',"*_platesurvey_fluid_volume_*primerplate.csv")), reverse=True)[0]

#def get_mouse_i7i5_survey():
#    return sorted(glob.glob(os.path.join('..','library',"*_platesurvey_fluid_volume_I7I5PLATE.csv")), reverse=True)[0]


# File stuff

def get_soft_delete_path(file_path):
    """
    When soft deleting a file rename it to /path/DEL_yyyymmdd_HHMMSS_filename.ext
    """
    t = datetime.datetime.now()
    t = t.strftime("%y%m%d_%H%M%S_")
    fn_orig = WindowsPath(file_path).name
    fn_root = WindowsPath(file_path).parent
    #print(f'**** {file_path} {fn_orig} {fn_root} ****')
    new_path = Path(fn_root)/('DEL_' + t + str(fn_orig))
    return str(new_path)


def delete_file(file_path, soft=True, caller_id=None):
    """
    By default use soft delete (renamimg)
    """
    file_path = WindowsPath(file_path)
    if str(file_path.name).startswith('pending_'):
        soft = False
    if Path(file_path).exists():
        do_delete = False
        if soft:
            del_path = get_soft_delete_path(file_path)
            if Path(del_path).exists():
                m(f'{del_path} already exists, cannot rename', level='warning', caller_id=caller_id)
                do_delete = True
            else:    
                try:
                    Path(file_path).rename(del_path)
                except Exception as exc:
                    m(f'failed to rename {file_path} to {del_path}, {exc}', level='error', caller_id=caller_id)
                    do_delete = True    
                if not do_delete:
                    m(f'renamed {file_path} to {del_path}', level='success', dest=('noGUI',))
                    return True
                
        if (not soft) or do_delete:
            try:
                Path(file_path).unlink()
            except Exception as exc:
                m(f'could not delete file {file_path} {exc}', level='error', caller_id=caller_id)
                return False
            m(f'deleted file {file_path}', level='info', dest=('noGUI',))
            return True
    else:
        m(f'file {file_path} not present', level='warning', dest=('noGUI',))
        return True


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
    #print(f'{w=}')
    try:
        well = w[0]+w[2] if w[1]=='0' else w
    except:
        print(f"unpadwell failed with {w=}")
    return well 

# Assays/Primers

def choose_primerfam(exp, assayfam, caller_id=None):
    """
    Return the primer family entry that best matches the provided assay family
    """
    if assayfam in exp.assay_primer:
        return exp.assay_primer[assayfam]
    if assayfam.lower() in exp.assay_primer:
        return exp.assay_primer[assayfam.lower()]
    
    m(f'Primer family not found to match assay {assayfam}', level='warning', caller_id=caller_id)

    return [assayfam]
   

def match_assays_to_primers(exp, assays):
    """
    Multiple primers may be needed per assay, and we need to check for manual overrides
    from exp.primer_assay too.
    """
    assay_primers = {}  # each assay will result in a set of primers
    loaded_primers = exp.get_primer_names()
    mapped_primers = exp.primer_assay
    for assay in assays:
        assay_primers[assay] = set()
        for p in mapped_primers:  # first try the manual overrides
            if exp.primer_assay[p] == assay:
                assay_primers[assay].add(p)
                break
            elif '_' in p:  # multiple primers per assay
                pfam = p.split('_')[0]
                if pfam == assay:  # we found the family
                    assay_primers[assay].add(p)
        for p in loaded_primers:  # compare to all the primers in plates
            if p == assay:
                if len(assay_primers[assay]) > 0:
                    break # already found
                assay_primers[assay].add(p)
                break  # got it in one!
            elif '_' in p:  # multiple primers per assay
                pfam = p.split('_')[0]
                if pfam == assay:  # we found the family
                    assay_primers[assay].add(p)
    return assay_primers

    @functools.lru_cache
    def get_taqwater_volumes_primer(self, num_reactions):
        """
        Returns the taq and water volumes required for the primer stage in nl
        Args:
            num_reactions(int):
        """
        taq_vol = num_reactions * self.transfer_volumes['PRIMER_TAQ_VOL']
        water_vol = num_reactions * self.transfer_volumes['PRIMER_WATER_VOL']
        return taq_vol, water_vol
        
    @functools.lru_cache
    def get_taqwater_req_vols_index(self, num_reactions):
        """
        Returns the taq and water volumes required for the index stages in nl
        Args:
            num_reactions(int): number of reactions for index stage, including amplicon
        """
        water_vol = num_reactions * self.transfer_volumes['INDEX_WATER_VOL']
        taq_vol = num_reactions * self.transfer_volumes['INDEX_TAQ_VOL'] 
        return taq_vol, water_vol


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
        #assert os.path.isfile(fn), "needs file {} in folder {}".format(fn, os.getcwd())
        #assert hdridx>0, 'hdridx ({}) must be a positive number.'.format(hdridx)
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
            wx = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_ALL)
            data = []
            for row in self.data:
                if output_plate_guards:  # guard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() else guard_pbc(str(d), silent=True) for (h,d) in zip(self.header, row)])
                else:  # unguard plate barcodes
                    data.append([str(d) if 'plate barcode' not in h.lower() else unguard_pbc(str(d), silent=True) for (h,d) in zip(self.header, row)])
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
        if hdridx < 0:
            print(f'header index must be at least zero {hdridx=}', file=sys.stderr)
            return
        self.filename = filename
        with open(filename, errors="ignore", newline='') as src:
            prefix = [x for i, x in zip(range(hdridx-1), src)] # read initial lines
            csvrdr = csv.reader(src)
            hdrrow = next(csvrdr) if hdridx else None 
            if clsname not in Table.tt:
                if not fields and not hdrrow:
                    return
                Table.newtype(clsname, fields if fields else hdrrow)
            filtfunc = lambda x: len(x)==len(hdrrow) # could do padding?
            # TODO: need a filter-type lambda to apply strip() to everything... covered for now but it's not a safe patch
            Table.__init__(self, clsname, data=filter(filtfunc, csvrdr), headers=hdrrow, prefix=prefix, selector=selector)
        

class CSVMemoryTable(Table):
    "Table in memory from StringIO"
    def __init__(self, clsname, io_string, hdridx=1, fields=None, selector=None):
        if hdridx < 0:
            print(f'header index must be at least zero {hdridx=}', file=sys.stderr)
            return
        if hdridx == 0:
            prefix = []
        else:
            prefix = [x for i, x in zip(range(hdridx-1), io_string)] # read initial lines - why?
        #print(f'{prefix=}', file=sys.stderr)
        csvrdr = csv.reader(io_string)
        hdrrow = next(csvrdr) if hdridx else None 
        #print(f'{hdrrow=}', file=sys.stderr)
        if clsname not in Table.tt:
            if not fields and not hdrrow:
                return
            Table.newtype(clsname, fields if fields else hdrrow)
        filtfunc = lambda x: len(x)==len(hdrrow) # could do padding?
        # TODO: need a filter-type lambda to apply strip() to everything... covered for now but it's not a safe patch
        Table.__init__(self, clsname, data=filter(filtfunc, csvrdr), headers=hdrrow, prefix=prefix, selector=selector)
        
