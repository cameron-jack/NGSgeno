#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from experiment import Experiment
import sys
import os
from pathlib import Path
import traceback

"""
@created: Nov 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: Custom manifests now require DNA plate and may support any number of Nimbus runs
@last_edit: 2022-05-05
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
GUARD_TYPES = set(['m','r','c','p', 'a'])

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


### cgi-echo.py ###


### Readers and writers for stage interlock CSV files (Nimbus -> Echo -> Analysis)

# Define headers for all stages. These can handle both mouse and custom samples
stage1_hdr = tuple(['sampleNo', 'samplePlate','sampleWell', 'strainName', 'sampleBarcode', 'assays', 'assayFamilies', 
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


#def match_nimbus_to_echo_files(exp):
#    """ make sure there is a 1-1 match of Nimbus to Echo-COC file. Needs a rewrite. It seems fundamentally borked """
#    try:
#        nfiles = [exp.get_exp_fp(fn) for fn in glob.glob('Nimbus-*.csv')]
#        # barcodes for Nimbus 384-well plate outputs
#        nbc = frozenset(PurePath(fn).name.replace('Nimbus-','').replace('.csv','') for fn in nfiles)
#        gnbc = set([fn for fn in nbc if file_io.is_guarded_pbc(fn)])
#        # identify missing Stage1 files.
#        #missing = [bc for bc in nbc if not os.path.isfile(f"Stage1-P{bc}.html")]
#        #if missing:
#        #    # print("<pre>")
#        #    # print("cwd =", os.getcwd())
#        #    # print("nbc =", ' '.join(nbc))
#        #    # print("missing =", ' '.join(missing))
#        #    # should catch and report output
#        #    if 'cust_manifest' in templates.files:
#        #            subprocess.run(["python", os.path.join("..", "bin", "stage1report.py"), '--custom']+missing)
#        #    else:
#        #        subprocess.run(["python", os.path.join("..", "bin", "stage1report.py")]+missing)
#        #    # print("</pre>")
    
#        # find the Nimbus output files for each DNA plate.
#        efiles = [exp.get_exp_fp(fn) for fn in glob.glob('Echo_384_COC_00??_*.csv')]

#        # check for Echo files without guards, matching Nimbus files with guards, then fix these files and rename original   
#        for e in efiles:
#            e_plate_parts = e.split('_')
#            e_plate = e.replace('.csv','').split('_')[4]  # plate barcode only
#            e_plate_prefix = '_'.join(e.split('_')[:4])
#            if len(e_plate_parts) > 5:
#                e_plate_suffix = '_'.join(e.split('_')[5:])
#            else:
#                e_plate_suffix = '.csv'
#            if file_io.is_guarded_pbc(e_plate):
#                continue
#            ge_plate = file_io.guard_pbc(e_plate)
#            if ge_plate in gnbc:  # unguarded echo_COC matches guarded Nimbus file
#                e_fn = e_plate_prefix + '_' + e_plate + '_' + e_plate_suffix
#                ge_fn = e_plate_prefix + '_' + ge_plate + '_' + e_plate_suffix
#                with open(e_fn, 'rt') as ef, open(ge_fn, 'wt') as gef:
#                    for i, line in enumerate(ef):
#                        if i == 0:
#                            gef.write(line)
#                            continue
#                        cols = line.split(',')
#                        plate_col = cols[1].strip('"')
#                        if not file_io.is_guarded_pbc(plate_col):
#                            cols[1] = '"' + file_io.guard_pbc(plate_col) + '"'
#                        outline = ','.join(cols)
#                        gef.write(outline)
#                os.rename(e, 'original_'+e)        
#        efiles = [exp.get_exp_fp(fn) for fn in glob.glob('Echo_384_COC_00??_*.csv')]
#        # pick out the plate barcode from the filename
#        ebc = frozenset(PurePath(fn).name.replace('.csv','').split('_',5)[4] for fn in efiles)

#        xbc  = nbc-ebc # any Nimbus plate files that are missing Nimbus run (output) files
#        #print(xbc, nbc, ebc, file=sys.stdout)
#        return nfiles, efiles, xbc
#    except Exception as exc:
#        output_error(exc, msg='Problem matching Nimbus to Echo files. match_nimbus_to_echo_files()')
#        return [], [], []

#def match_nimbus_to_echo_files(exp: Experiment) -> tuple(list,list,list):
#    """ Make sure there is a 1-1 match of Nimbus to Echo-COC file.
#        Assumes plate barcodes are guarded.
#        Inputs: an Experiment instance
#        Outputs: three lists (nimbus files, echo files, barcodes not matches by echo files)
#   """
#    try:
#        # Python 3.10
#        #nfiles = [exp.get_exp_fp(fn) for fn in glob.glob('Nimbus-*.csv', root_dir='run_'+exp.name)]
#        #efiles = [exp.get_exp_fp(fn) for fn in glob.glob('Echo_384_COC_00??_*.csv', root_dir='run_'+exp.name)]
#        # Python 3.9
#        #print(glob.glob('*'))
#        nim_str = os.path.join('run_'+exp.name, 'Nimbus-*.csv')
#        echo_str = os.path.join('run_'+exp.name, 'Echo_384_COC_00??_*_0.csv')
#        nfiles = glob.glob(nim_str) #[exp.get_exp_fp(fn) for fn in glob.glob(nim_str)]
#        efiles = glob.glob(echo_str) #[exp.get_exp_fp(fn) for fn in glob.glob(echo_str)]
#        # cleave off the path and file name junk to just get the barcodes
#        nbc = frozenset(PurePath(fn).name.replace('Nimbus-','').replace('.csv','') for fn in nfiles)
#        ebc = frozenset(PurePath(fn).name.replace('.csv','').split('_',5)[4] for fn in efiles)
#        xbc  = nbc-ebc # any Nimbus plate files that are missing Nimbus run (output) files
#        #print(f"{nbc=} {ebc=} {xbc=}")
#        return nfiles, efiles, xbc
#    except Exception as exc:
#        exp.log('Problem matching Nimbus to Echo files', level='Error')
#        exp.log(traceback.print_exc(limit=2), level='Error')    
#        return [], [], []

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
