#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Jan 14 2022
@author: Cameron Jack - ANU Bioinformatics Consultancy JCSMR, ANU
@version: 0.15
@version_comment: New! Shared code for CGI interface.
@last_edit:
@edit_comment: 
"""

# Ugly hack to find things in ../bin
import sys
import os
import glob
import subprocess
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(os.path.join(file.parents[1],'bin'))
from util import output_error, ERROR_FN
import file_io

port=9123
stage1 = "/cgi-bin/cgi-nimbus.py"
stage2 = "/cgi-bin/cgi-echo.py"
stage3 = "/cgi-bin/cgi-analysis.py"

nl = '\n'

def get_error_msg():
    with open(ERROR_FN, 'rt') as ef:
        lines = [l for l in ef]
    return "<p>" + "<br>".join(lines) + "</p>"


def getinfo(ngid):
    inf = glob.glob('*.html')+glob.glob('*.log')
    if inf:
        return "<ul>\n  "+ \
                "\n    ".join("<li><a href='../{ngid}/{fn}'>{fn}</a></li>".format(ngid=ngid, fn=x) for x in sorted(inf))+ \
                "\n  </ul>"
    return '' 


def match_nimbus_to_echo_files(templates):
    try:
        nfiles = glob.glob('Nimbus-*.csv')
        # barcodes for Nimbus 384-well plate outputs
        nbc = frozenset(fn.replace('Nimbus-','').replace('.csv','') for fn in nfiles)
        gnbc = set([fn for fn in nbc if file_io.is_guarded_pbc(fn)])
        # identify missing Stage1 files.
        missing = [bc for bc in nbc if not os.path.isfile(f"Stage1-P{bc}.html")]
        if missing:
            # print("<pre>")
            # print("cwd =", os.getcwd())
            # print("nbc =", ' '.join(nbc))
            # print("missing =", ' '.join(missing))
            # should catch and report output
            if 'cust_manifest' in templates.files:
                    subprocess.run(["python", os.path.join("..", "bin", "stage1report.py"), '--custom']+missing)
            else:
                subprocess.run(["python", os.path.join("..", "bin", "stage1report.py")]+missing)
            # print("</pre>")
    
        # find the Nimbus output files for each DNA plate.
        efiles = glob.glob('Echo_384_COC_00??_*.csv')

        # check for Echo files without guards, matching Nimbus files with guards, then fix these files and rename original   
        for e in efiles:
            e_plate_parts = e.split('_')
            e_plate = e.replace('.csv','').split('_')[4]  # plate barcode only
            e_plate_prefix = '_'.join(e.split('_')[:4])
            if len(e_plate_parts) > 5:
                e_plate_suffix = '_'.join(e.split('_')[5:])
            else:
                e_plate_suffix = '.csv'
            if file_io.is_guarded_pbc(e_plate):
                continue
            ge_plate = file_io.guard_pbc(e_plate)
            if ge_plate in gnbc:  # unguarded echo_COC matches guarded Nimbus file
                e_fn = e_plate_prefix + '_' + e_plate + '_' + e_plate_suffix
                ge_fn = e_plate_prefix + '_' + ge_plate + '_' + e_plate_suffix
                with open(e_fn, 'rt') as ef, open(ge_fn, 'wt') as gef:
                    for i, line in enumerate(ef):
                        if i == 0:
                            gef.write(line)
                            continue
                        cols = line.split(',')
                        plate_col = cols[1].strip('"')
                        if not file_io.is_guarded_pbc(plate_col):
                            cols[1] = '"' + file_io.guard_pbc(plate_col) + '"'
                        outline = ','.join(cols)
                        gef.write(outline)
                os.rename(e, 'original_'+e)        
        efiles = glob.glob('Echo_384_COC_00??_*.csv')
        # pick out the plate barcode from the filename
        ebc = frozenset(fn.replace('.csv','').split('_',5)[4] for fn in efiles)

        xbc  = nbc-ebc # any Nimbus plate files that are missing Nimbus run (output) files
        #print(xbc, nbc, ebc, file=sys.stderr)
        return xbc
    except Exception as exc:
        output_error(exc, msg='Problem matching Nimbus to Echo files. match_nimbus_to_echo_files()')
