#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Jan 14 2022
@author: Cameron Jack - ANU Bioinformatics Consultancy JCSMR, ANU
@version: 0.16
@version_comment: 
@last_edit: 2022-02-23
@edit_comment: Moved match_nimbus_to_echo() to util.py
@description: common components for the CGI interface. This file may disappear under the new layout scheme
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


