#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
import os

"""
@created: Nov 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.15
@version_comment: No changes
@last_edit: 
@edit_comment: 

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

### Helper functions for IO

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

# global - set of all possible guard types
GUARD_TYPES = set(['m','M','c','p'])

# Confirm guards
def is_guarded_mbc(mbc):
    """ Return True if Musterer barcode is correctly guarded, else False """
    if mbc.startswith('m') and mbc.endswith('m'):
        return True
    return False

def is_guarded_rbc(rbc):
    """ Return True if Rodentity barcode is correctly guarded (MNNNNNm), else False """
    if rbc.startswith('M') and rbc.endswith('m'):
        return True
    return False

def is_guarded_cbc(cbc):
    """ Return True if custom barcode is correctly guarded, else False """
    if cbc.startswith('c') and cbc.endswith('c'):
        return True
    return False

def is_guarded_pbc(pbc):
    """ Return True if plate barcode is correctly guarded, else False """
    if pbc.startswith('p') and pbc.endswith('p'):
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
    return False


# Add guards to barcodes
def guard_mbc(mbc):
    """ Add guards to Musterer barcode """
    if type(mbc) != str:
        raise TypeException(f"Musterer barcode is not a string! {mbc} is type {type(mbc)}")
    if mbc == '':
        raise EmptyBarcodeError(f"Musterer barcode is zero length string!")
    if is_guarded_mbc(mbc):
        raise ExistingGuardError(f"Musterer barcode appears to already have Musterer guards: {mbc}")
    if is_guarded_rbc(mbc):
        raise ExistingGuardError(f"Musterer barcode appears to already have Rodentity guards: {mbc}")
    if is_guarded_cbc(mbc):
        raise ExistingGuardError(f"Musterer barcode appears to already have custom guards: {mbc}")
    if is_guarded_pbc(mbc):
        raise ExistingGuardError(f"Musterer barcode appears to already have plate guards: {mbc}")
    return 'm' + mbc + 'm'

def guard_rbc(rbc):
    """ Add guards to Rodentity barcode """
    if type(rbc) != str:
        raise TypeError(f"Rodentity barcode is not a string! {rbc} is type {type(rbc)}")
    if rbc == '':
        raise EmptyBarcodeError(f"Rodentity barcode is zero length string!")
    if is_guarded_mbc(rbc):
        raise ExistingGuardError(f"Rodentity barcode appears to already have Musterer guards: {rbc}")
    if is_guarded_rbc(rbc):
        raise ExistingGuardError(f"Rodentity barcode appears to already have Rodentity guards: {rbc}")
    if is_guarded_cbc(rbc):
        raise ExistingGuardError(f"Rodentity barcode appears to already have custom guards: {rbc}")
    if is_guarded_pbc(rbc):
        raise ExistingGuardError(f"Rodentity barcode appears to already have plate guards: {rbc}")
    return rbc + 'M'  # Rodentity barcodes already start with 'M'

def guard_cbc(cbc):
    """ Add guards to custom barcode """
    if type(cbc) != str:
        raise TypeError(f"Custom barcode is not a string! {cbc} is type {type(cbc)}")
    if cbc == '':
        raise EmptyBarcodeError(f"Custom barcode is zero length string!")
    if is_guarded_mbc(cbc):
        raise ExistingGuardError(f"Custom barcode appears to already have Musterer guards: {cbc}")
    if is_guarded_rbc(cbc):
        raise ExistingGuardError(f"Custom barcode appears to already have Rodentity guards: {cbc}")
    if is_guarded_cbc(cbc):
        raise ExistingGuardError(f"Custom barcode appears to already have custom guards: {cbc}")
    if is_guarded_pbc(cbc):
        raise ExistingGuardError(f"Custom barcode appears to already have plate guards: {cbc}")
    return 'c' + cbc + 'c'

def guard_pbc(pbc):
    """ Add guards to plate barcode """
    if type(pbc) != str:
       raise TypeError(f"Plate barcode is not a string! {pbc} is type {type(pbc)}")
    if pbc == '':
        raise EmptyBarcodeError(f"Plate barcode is zero length string!")
    if is_guarded_mbc(pbc):
        raise ExistingGuardError(f"Plate barcode appears to already have Musterer guards: {pbc}")
    if is_guarded_rbc(pbc):
        raise ExistingGuardError(f"Plate barcode appears to already have Rodentity guards: {pbc}")
    if is_guarded_cbc(pbc):
        raise ExistingGuardError(f"Plate barcode appears to already have custom guards: {pbc}")
    if is_guarded_pbc(pbc):
        raise ExistingGuardError(f"Plate barcode appears to already have plate guards: {pbc}")
    return 'p' + pbc + 'p'


# Unguard barcodes
def unguard_mbc(mbc):
    """ remove guards from a Musterer barcode """
    if type(mbc) != str:
        raise TypeError(f"Musterer barcode is not a string! {mbc} is type {type(mbc)}")
    if not mbc.startswith('m') or not mbc.endswith('m'):
        raise UnguardedBarcodeError(f"Musterer barcode guards degraded or missing in: {mbc}")
    return mbc[1:-1]

def unguard_rbc(rbc):
    """ remove guards from a Rodentity barcode """
    if type(rbc) != str:
        raise TypeError(f"Rodentity barcode is not a string! {rbc} is type {type(rbc)}")
    if not rbc.startswith('M') or not rbc.endswith('m'):
        raise UnguardedBarcodeError(f"Musterer barcode guards degraded or missing in: {rbc}")
    return rbc[0:-1]  # Rodentity barcodes keep their 'M' front guard

def unguard_cbc(cbc):
    """ remove guards from a custom barcode """
    if type(cbc) != str:
        raise TypeError(f"Custom barcode is not a string! {cbc} is type {type(cbc)}")
    if not cbc.startswith('c') or not cbc.endswith('c'):
        raise UnguardedBarcodeError(f"Custom barcode guards degraded or missing in: {cbc}")
    return cbc[1:-1]

def unguard_pbc(pbc):
    """ remove guards from a plate barcode """
    if type(pbc) != str:
        raise TypeError(f"Plate barcode is not a string! {pbc} is type {type(pbc)}")
    if not pbc.startswith('p') or not pbc.endswith('p'):
        raise UnguardedBarcodeError(f"Plate barcode guards degraded or missing in: {pbc}")
    return pbc[1:-1]

def unguard(bc, silent=False):
    """ shortcut, for when we don't know what kind of guarding we have but need to remove it for output """
    if is_guarded_mbc(bc):
        return unguard_mbc(bc)
    elif is_guarded_rbc(bc):
        return unguarded_rbc(bc)
    elif is_guarded_cbc(bc):
        return unguard_cbc(bc)
    elif is_guarded_pbc(bc):
        return unguard_pbc(bc)
    else:
        if not silent:
            print('Possibly unguarded already:',bc, file=sys.stderr)
        return bc

def myopen(fn):
    """ Bob's function to handle gzip transparently """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt")
    return open(fn, errors="ignore")

### cgi-echo.py ###




### Readers and writers for stage interlock CSV files (Nimbus -> Echo -> Analysis)

# Define headers for all stages. These can handle both mouse and custom samples
stage1_hdr = tuple(['sampleNo', 'samplePlate','sampleWell', 'strainName', 'barcode', 'assays', 'assayFamilies', 
              'mouseID', 'strainName', 'sex', 'parentLitter', 'parentSequence', 'parentStrain', 'parentDate'])
stage2_hdr = tuple(list(stage1_hdr) + ['dnaPlate', 'dnaWell', 'primer', 'pcrPlate', 'pcrWell'])
stage3_hdr = tuple(list(stage2_hdr) + ['i7bc', 'i7name', 'i7well', 'i5bc', 'i5name', 'i5well'])

#def write_stage1_csv() Is the contents of nimbus.py, so not needed here


def read_stage1_csv():
    """ Called by Stage1report.py, cgi-echo.py and possibly echo_primer.py (currently uses Table() constructor) """
    pass

def write_stage2_csv():
    """ Called by echo_primer.py """
    pass

def read_stage2_csv():
    """ Called by stage2report.py and echo_barcode.py """
    pass

def write_stage3_csv():
    """ Called by echo_barcode.py """
    pass

def read_stage3_csv():
    """ Called by ngsmatch.py and stage2_report.py """
    pass


if __name__ == '__main__':
    pass
