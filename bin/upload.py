#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
@created: July 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.12
@version_comment: Still not included in production pipeline
@last_edit: 
@edit_comment:

upload.py, a script to:
1) read the results_gt.csv, sanity-checked genotype results file
2) create a 3-column CSV file for upload to Musterer (later Rodentity) via the available webservice
3) Perform the upload and handle any error messages
"""

import sys
import argparse
from collections import defaultdict
import requests
try:
    from bin.db_io import upload_plate
except ModuleNotFoundError:
    from db_io import upload_plate


class SilentFail(Exception):
    """
    Bob Buckley's magic...
    Allows the program to exit an abnormal situation without appearing to have failed.
    Useful if you need a bash script with -e to keep running if a called Python script fails.
    Or maybe if executed from an interactive app?
    ... otherwise a highly questionable design pattern.
    """
    pass


def parse_results(results_fn, required_cols):
    """
    Return for  entries that completely pass genotyping (all tested assays pass and all sanity checks pass):
    ear-punch plate & well position, strain, mouseID, assay, reportable genotype.

    Results_gt.CSV header structure: 
    SampleNo	epplate	epwell	mouseBarcode	mouseID strainName	mouseAssays	assayFamilies	\
    sex	parentLitter	parentSequence	parentStrain	parentDate	dnaplate	dnawell	primer	pcrplate	\
    pcrwell	i7bc	i7name	i7well	i5bc	i5name	i5well	call	gt	alleleRatios	efficiency	pass_heredity	\
    sanity_comment	sire_barcode	sire_strain	sire_gt	dam_barcode	dam_strain	dam_gt	\
    readCount	cleanCount	mergeCount	seqCount	seqName										

    We need mouseID, strainName, EPplate, EPwell, GT, pass_heredity
    """
    records = []
    with open(result_fn, 'rt') as src_fh:
        src = csv.reader(src_fh)
        hdr = next(src)
        names_cols = [(col,i) for i,col in enumerate(hdr) if col in required_cols]
        names = [name for name, col in names_cols]
        cols = [col for name, col in names_cols]
        GTrec = namedtuple('GTrec', names)
        for row in src:
            records.append(GTrec([row[i] for i in cols]))
    return records


def make_dummy_records(required_cols):
    """
    Create some dummy records for an attempted upload
    """
    GTrec = namedtuple('GTrec', required_cols)
    records = []
    records.append(GTrec(['3','ENU24:017:Tmem131:B6:G6','25708','6C','mut/mut','True']))
    records.append(GTrec(['4','ENU24:017:Tmem131:B6:G6','25708','6D','mut/wt','True']))
    return records





def main():
    """
    1) parse results_gt.csv for mouse ID, assay, genotype and ultimately a pass value
    2) output CSV upload file
    3) perform upload
    4) report status
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('resultfile', help='Path to results_gt.csv file from genotyping stage')
    parser.add_argument('-o', '--outputfile', default='musterer_upload.csv', help='Path to output musterer upload CSV')
    parser.add_argument('-t', '--test', action='store_true', help='Create a dummy test and perform it against the Musterer Test DB')
    args = parser.parse_args()

    out_fn = args.outputfile
    required_cols = set(['mouseID', 'strainName', 'EPplate', 'EPwell', 'GT', 'passHeredity'])
    if args.test:
        records = make_dummy_records(required_cols)
        out_fn = 'test_musterer_upload.csv'
    else:
        records = parse_results(args.resultfile, required_cols)
    passing_records = {r for r in records if records.passHeredity.lower() == 'true'}
    write_uploadCSV(args.outputfile, passing_records)
    upload_plate(args.outputfile)


if __name__ == '__main__':
    try:
        main()
    except SilentFail:
        pass
