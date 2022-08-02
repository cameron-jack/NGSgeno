#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Feb 22 2022
@author: Cameron Jack - ANU Bioinformatics Consultancy JCSMR, ANU
@version: 0.16
@version_comment: New application-style interface for the "pipeline"
@last_edit: 2022-02-22
@edit_comment: New file. Introducing CSS and JS for the first time here.

Define an Experiment and fill the inventory.
"""

import os
import glob
import sys
import cgi
import cgitb
import subprocess
from cgi_util import port, nl, stage1, stage2, stage3, getinfo

# Ugly hack to find things in ../bin
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(os.path.join(file.parents[1],'bin'))

from musterer import get_plate
import validate_assays
from file_io import read_html_template
from util import TemplateFiles, output_error, match_nimbus_to_echo_files
import nimbus


### Compose the interface
# Import fixed components from a template file in the library folder
# The interface is a standard form that is made up of a number of collapseable sections
# We want to start with the main form shape then add in sections:
# 1) Errors
# 2) Run folder
# 3) Inventory - sets of 4 custom/mouse plates + 1 DNA plate; barcode plate + barcode volumes; assay file; import DB info
# - Now run Nimbus
# 4) Primers - plate(s) + matching volume(s)
# - Run Echo PCR 1
# 5) Add extra samples here before PCR 2 barcoding
# 6) Miseq
# 7) Analysis
###

htmlerr = """
<html>
<head> 
<title>NGS Application!</title>
<meta charset="utf-8">
</head>
<body>
  <h1>NGS Genotyping pipeline</h1>
  <div class="errs" style="width: 800px; border: 4px solid red; border-radius: 10px; padding: 10px; background-color: lightpink;">
    <h2>Errors ...</h2>
    {errs}
  </div>
  <p>
    <button onclick="window.history.back();">Back to previous page</button>
  </p>
  <p>
    <button onclick="location.href='http://localhost:{port}/cgi-bin/cgi-nimbus.py';" value="Start again from first pipeline page">Back to start of pipeline</button>
  </p>
</body>
</html>
""".strip()


