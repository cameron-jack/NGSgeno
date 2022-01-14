#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Nov 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.14
@version_comment: Refined interface 
@last_edit: 2022-01-07
@edit_comment: Various refinements to options

Application Webpage for post-sequencing analysis
Need to support match counting of new reference sequences for trialing probes
Need to support various UI options for matching

A single page that calls analysis stages as executed programs
"""

import os
import re
import csv
import glob
import json
import cgi
import cgitb
import requests
import subprocess
import collections
import sys
# Ugly hack to find things in ../bin
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(os.path.join(file.parents[1],'bin'))

# add any /bin imports here
from util import get_mouse_ref, TemplateFiles

port=9123
stage1 = "/cgi-bin/cgi-nimbus.py"
stage2 = "/cgi-bin/cgi-echo.py"
stage3 = "/cgi-bin/cgi-analysis.py"


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


html_analysis = """
<html>
<head> 
<title>NGS Genotyping - Sequencing analysis</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline - Sequencing analysis</h1>
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
  <h2>Sequence matching</h2>
  <h3>Leave these blank to use pipeline defaults</h3>
  <form action="{stage}" method="get">
    <p>
    {rerun}
    <label style="margin-left:25px;" for="outfile">Output filename</label>
      <input type="text" id="outfile" name="outfile" size="50" value="Results.csv"/><br>
    <br>
    <input style="margin-left:25px;" type="checkbox" id="custom" name="custom"></input>
      <label for="custom">Custom. Musterer data columns not present</label><br>
    <input style="margin-left:25px;" type="checkbox" id="exhaustive" name="exhaustive"></input>
      <label for="exhaustive">Exhaustive mode. Match every sequence not matter how low the coverage. Slow!</label>
    <br>
    <br>
    <label style="margin-left:25px;" for="ncpus">Number of task processes to run simultaneously</label>
      <input type="number" min="1" value="4" max="32"></input><br>
    <input style="margin-left:25px" type="checkbox" id="clear_caches" name="clear_caches" title="Clear cache files"></input>
      <label for="clear_caches" title="Clear the caches if you have changed target references">Clear cache files</label><br>
    <input style="margin-left:25px" type="checkbox" id="debug" name="debug"/>
      <label for="debug">Debug - give more detailed debugging info</label><br>
    <input style="margin-left:25px" type="checkbox" id="quiet" name="quiet"/>
      <label for="quiet">Quiet - run with minimal screen output. Normal output is sent to log file instead</label>
    <label for="customref">Add a custom reference file for matching custom samples only (OPTIONAL)</label>
      <input style="margin-left:25px" type="file" id="customref" name="customref"/>
    </p>
    <h2> Genotyping and Sanity Checking</h2>
    <p>                                               
    <label style="margin-left:25px;" for="gtoutfile" title="Defaults to results_gt.csv">Per-row output filename</label>
      <input type="text" id="gtoutfile" name="gtoutfile" size="50" value="results_gt.csv" title="Defaults to results_gt.csv"/><br>
    <label style="margin-left:25px;" for="gtmice" title="Defaults to mice_gts.csv">Per-mouse GTs filename</label>
      <input type="text" id="gtmice" name="gtmice" size="50" value="mice_gts.csv" title="Defaults to mice_gts.csv"/><br>
    <label style="margin-left:25px;" for="gtupload" title="Defaults to mice_uploads.csv">Mouse DB upload filename</label>
      <input type="text" id="gtupload" name="gtupload" size="50" value="mice_uploads.csv" title="Defaults to mice_uploads.csv"/><br>
    </p>
    <p>Click the "Continue" button below when done.</p>
    <input type="hidden" id="existing" name="existing" value="{existingDir}"/>
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}"></input>
  </form>
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



html_no_raw = """
<html>
<head> 
<title>NGS Application!</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline</h1> 
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
  <h2>Awaiting MiSeq results</h2>
  <p>When the MiSeq run completes and you have copied the MiSeq FASTQ files into 
  the "raw" folder for run {ngid} click the Continue button.</p>
  <p>Note: this step involves a lot of processing so expect to wait a while
  for it to complete.</p>
  <form action="{stage}" method="get">
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}" />
    <input type="hidden" id="existing" name="existing" value="{existingDir}"/>
  </form>
  </div>
  {info}
  <p>
    <button onclick="window.history.back();">Back to previous page</button>
  </p>
  <p>
    <button onclick="location.href='http://localhost:{port}/cgi-bin/cgi-nimbus.py';" value="Start again from first pipeline page">Back to start of pipeline</button>
  </p>
</body>
</html>
""".strip()

html_results = """
<html>
<head> 
<title>NGS Application!</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline</h1> 
  <h2>{ngid} results</h2>
  {info}
  <p>
    <button onclick="location.href='http://localhost:{port}/cgi-bin/cgi-analysis.py?ngid={ngid}&analysis=true'">
        Show analysis options</button>
  </p>
  <p>
    <button onclick="window.history.back();">Back to previous page</button>
  </p>
  <p>
    <button onclick="location.href='http://localhost:{port}/cgi-bin/cgi-nimbus.py';" value="Start again from first pipeline page">
        Back to start of pipeline</button>
  </p>
</body>
</html>
""".strip()

class ProgramFail(Exception):
    pass

def eppfn(bc):
    return "Plate{}-EP.json".format(bc) # ear-punch plate (filename) format

mustfix = re.compile(r'(,\s*"(mouseBarcode|plateId)":\s*)(\d+)')



#  
def main():
    def getinfo():
        inf = glob.glob('*.html')+glob.glob('*.log')
        if inf:
            return "<ul>\n  "+ \
                    "\n    ".join("<li><a href='../{ngid}/{fn}'>{fn}</a></li>".format(ngid=ngid, fn=x) for x in sorted(inf))+ \
                    "\n  </ul>"
        return ''  
    fields = cgi.FieldStorage()
    errs, info  = "", []
    cgitb.enable(display=1, logdir=".")
    print("Content-Type: text/html")    # HTML is following
    print()                             # blank line, end of headers
    
    global stage1
    global stage2
    global stage3

    # check for existing experiment folder
    from datetime import date
    now = date.today() # get today's date
    nowstr = str(now.year)+str(now.month).zfill(2)+str(now.day).zfill(2)
    if 'ngid' in fields:
        ngid = fields.getfirst('ngid')
    elif 'existing' in fields:
        ngid = fields.getfirst('existing')
    elif 'projectDir' in fields:
        ngid = 'custom_'+fields.getfirst('projectDir').replace('custom_','')  # in case the user inputs "custom_"
    else:
        ngid = 'mouse_'+nowstr

    #if not os.path.isdir(ngid): # only if dnap is also defined?
    #    os.mkdir(ngid)
    os.chdir(ngid) # change to working directory - probably unnecessary

    ### Everything below this line happens in the run folder ###
    template_files = TemplateFiles()

    default_reference = get_mouse_ref()
    
    global port
    global htmlerr

    if not os.path.isdir("raw") and not os.path.isdir("RAW"):
        print('cgi-analysis:', 'No raw FASTQ directory', file=sys.stderr)
        global html3 # data analysis & reporting step
        print(html3.format(ngid=ngid, info=getinfo(), stage=stage3, existingDir=ngid, port=port))
        return
    if 'analysis' in fields:
        rerun = """
        <input style="margin-left:25px;" type="checkbox" id="rerun" name="rerun"></input>
        <label for="rerun"><b>Rerun analysis</b> - change output names to avoid overwriting</label><br>
        <br>
        """
        if not os.path.exists("Results.csv"):
            rerun = ""
        print(html_analysis.format(ngid=ngid, info=getinfo(), stage=stage3, existingDir=ngid, rerun=rerun, port=port))
        return
    
    if not os.path.isfile("Results.csv") or 'rerun' in fields:
        # parse html_analysis form results
        params = []

        if "outfile" in fields:
            params.append('-o')
            params.append(fields.getfirst('outfile'))
        else:
            params.append('-o')
            params.append('Results.csv')
        #if "targets" in fields:
        #    params.append('-t')
        #    params.append(fields.getfirst('targets'))
        if "custom" in fields:
            params.append('--custom')
        if "exhaustive" in fields:
            params.append('-x')
        if "ncpus" in fields:
            params.append('-n')
            params.append(fields.getfirst('ncpus'))
        if "clear_caches" in fields:
            params.append('-c')
        if "debug" in fields:
            params.append('-d')
        if "quiet" in fields:
            params.append('-q')
        #if "stage3file" in fields:
        #    params.append(fields.getfirst("stage3file"))
        #else:
        #params.append('Stage3.csv')

        cmd = ["python", os.path.join("..", "bin", "ngsmatch.py")] + params
        print('cgi-analysis:', 'Running stage3 FASTQ analysis:', cmd, file=sys.stderr)
        res = subprocess.run(cmd, capture_output=True)
        # capture output to log file
        with open("match.log", "wt") as dst:
            dst.write(res.stdout.decode('utf-8').replace('\r',''))
            if bool(res.stderr):
                dst.write("============ STDERR =============\n")
                dst.write(res.stderr.decode('utf-8').replace('\r',''))
    else:
        params = []
        params.append('-m')
        if "gtmice" in fields:
            params.append(fields.getfirst('gtmice'))
        else:
            params.append('mice_gts.csv')
        params.append('-u')
        if "gtupload" in fields:
            params.append(fields.getfirst('gtupload'))
        else:
            params.append('mice_uploads.csv')
        if "gtconfig" in fields:
            params.append('-k')
            params.append(fields.getfirst('gtconfig'))
        params.append('-o')
        if "gtoutfile" in fields:
            params.append(fields.getfirst('gtoutfile'))
        else:
            params.append('results_gt.csv')
        # positional params
        if "outfile" in fields:
            params.append(fields.getfirst('outfile'))
        else:
            params.append('Results.csv')
     
        cmd = ["python", os.path.join("..", "bin", "gt_mice.py")]+params
        print('Running Genotyping:', cmd, file=sys.stderr)
        res = subprocess.run(cmd, capture_output=True)
        # capture output to log file
        with open("gt_mice.log", "wt") as dst:
            dst.write(res.stdout.decode('utf-8').replace('\r',''))
            if bool(res.stderr):
                dst.write("============ STDERR =============\n")
                dst.write(res.stderr.decode('utf-8').replace('\r',''))
        
                
    print(html_results.format(ngid=ngid, info=getinfo(), stage=stage3, projectDir=ngid, port=port))
    return    
    

if __name__=="__main__":
    try:
        main()
    except ProgramFail:
        pass