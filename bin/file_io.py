#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

"""
@created: Nov 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.13
@version_comment: Module created
@last_edit: 
@edit_comment: 

Contains all file IO methods that involved barcodes, and are not purely for user reporting or robot inputs. 
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
    return rbc + 'm'  # Rodentity barcodes already start with 'M'

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

### cgi-echo.py ###

def readCustomCSVtoJSON(input_fn):
    """ read a custom sample manifest with no more than 4x96 well plates! Then return as JSON """
    data = {}
    errs = []  # error messages
    with open(input_fn, 'rt', newline=None) as f:
        in_csv = csv.reader(f, delimiter=',')
        for i, row in enumerate(in_csv):
            #print(row, file=sys.stderr)
            #print(file=sys.stderr)
            if i == 0:
                continue  # header
            if all(row == ''):
                continue  # blank rows
            gplate = guard_pbc(row[1].strip())
            well = row[2].strip()
            gsampleBarcode = guard_cbc(row[3].strip())
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
           
    pids = data.keys()
    data = [(data[p],p) for p in sorted(data)]
    return data, pids, errs

