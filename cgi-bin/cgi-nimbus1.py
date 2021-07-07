#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Aug 10 2020
@author: Bob Buckley & Cameron Jack - ANU Bioinformatics Consultancy JCSMR, ANU
@version: 0.7
@version_comment: No changes
@last_edit: 2021-07-07
@edit_comment:  Improvements to user interface and error reporting. Allows users to enter custom pipeline folder name

Application Stage 1 webpage.
Display a form for creating Stage 1 Nimbus picklists files.
This should validate ear-punch plate barcodes before doing submit action?
Should check valid DNA plate barcode if there are any EP barcodes?
Allows the optional choosing of a custom input CSV, side-stepping Musterer 
    look-up the mouse specific parts of the pipeline.
"""

import os
import glob
#import cgi
import cgitb

nimbus2 = "/cgi-bin/cgi-nimbus2.py"

html = """
<html>
<head>
<title>NGS Application!</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipelines (Mouse &amp; Custom)</h1>
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
    <h2>Mouse pipeline - Stage 1</h2>
    <p>Stage 1 writes picklist files for the Nimbus robot to tranfer material from 
    96-well ear-punch plates to a 384-well DNA sample plate for use in the Echo robot.
    </p>
    <form action="{stage1}" method="get">
      <p>Please provide the following information:</p>
        {eps}
        <br>
        <label for="dnap">DNA (Echo) plate barcode:</label>
          <input type="text" id="dnap" name="dnap" size="10" required/><br>
        <br>
        <input type="submit" value="Create Nimbus Picklist" style="border-color:red"/>
    </form>
    <hr>
    <form action="{stage1}" method="get">
        <h2>Custom sample pipeline - Stage 1</h2>
        <p>To run the custom sample pipeline, you must first create a "run" folder within the NGSgeno folder for your project. 
        You then need to copy the required files there: Sample file - max 4 plates per file!(CSV), Assay list file (CSV), Custom primer-plate layout (CSV)</p>
        <p>Select project folder</p>
        <label style="margin-left:25px;" for="projectDir">Set project folder:</label>
          <input type="text" id="projectDir" name="projectDir" size="40"/><br><br>

        <h3>Select custom sample, assay and primer plate files (CSV format)</h3>
        <p>These will be used instead of the corresponding Musterer information. These files MUST be in the run folder.</p>
         <label style="margin-left:25px;" for="customSamples">Samples/assays file:</label>
          <input type="file" id="customSamples" name="customSamples" size="80" accept=".csv"/><br><br>
         <label style="margin-left:25px;" for="customAssays">Assay list file:</label>
          <input type="file" id="customAssays" name="customAssays" size="80" accept=".csv"/><br><br>
         <label style="margin-left:25px;" for="customPrimers">Custom primer-plate layout:</label>
          <input type="file" id="customPrimers" name="customPrimers" size="80" accept=".csv" /><br><br>
        <label for="dnap">DNA (Echo) plate barcode:</label>
          <input type="text" id="dnap_cust" name="dnap_cust" size="10" required/><br>
        <br>
        <input type="submit" value="Create Nimbus Picklist" style="border-color:red"/>
    </form>
    <hr>
    <h2>Resume existing pipeline project</h2>
    {optany}
  </div>
  
  <p>
    <button onclick="location.reload();">Refresh form</button>
  </p>
</body>
</html>
""".strip()

#  <p>
#    <button onclick="window.history.back();">Back</button>
#  </p>

class ProgramFail(Exception):
    pass

def main():
    # form = cgi.FieldStorage()
    cgitb.enable(display=1, logdir=".")
    print("Content-Type: text/html")    # HTML is following
    print()                             # blank line, end of headers

    optx = ''
    # dir = form.getfirst('dir')
    # Find optional sample directories - they have P*-EP.json files but not Echo*.csv files
    # If the user doesn't choose a directory, then a new sample directory will be created
    jfiles = frozenset(map(os.path.dirname, glob.glob(os.path.join('*', 'P*-EP.json'))))
    # probably should look at Nimbus*.csv files - not Echo files
    efiles = frozenset(map(os.path.dirname, glob.glob(os.path.join('*', 'Echo_384_COC*.csv'))))
    dirs = sorted(jfiles-efiles, reverse=True)
    if dirs:
        dx = '\n\t\t    '.join("<option value='{0}'>{0}</option>".format(d) for d in dirs)
        optx = """<label for="ngid">NGS Geno run Id:</label>
        <select name="ngid" id="ngid">
            <option value="">New MiSeq run</option>
            {}
        </select><br>
        <p>Leave this (choose 'New MiSeq run') if you are preparing a new MiSeq run. 
        The system will assign a new date-based run Id.
        <br></p>""".format(dx)  
    epfmt = '<label for="ep{0}">Ear-punch plate {0} barcode:</label>\n\t  <input type="text" id="ep{0}" name="ep{0}" size="6" /><br>'
    epx = '\n\t'.join(epfmt.format(i) for i in range(1,5))

    optform = ''
    global nimbus2
    wdirs = sorted(jfiles, reverse=True)
    if wdirs:
        dx = '\n\t\t    '.join("<option value='{0}'>{0}</option>".format(d) for d in wdirs)
        optmiseq = """<label for="ngid2">NGS Geno run Id:</label>
            <select name="existing" id="existing">
                <option value="">None</option>
                {}
            </select><br>""".format(dx)
        optform = """  <p> </p>
          <form action="{stage1}" method="get">
              {optall}
              <input type="submit" value="Go to MiSeq run" />
          </form>""".format(stage1=nimbus2, optall=optmiseq)
    print(html.format(options=optx, optany=optform, eps=epx, stage1=nimbus2 ))
    return

if __name__=="__main__":
    try:
        main()
    except ProgramFail:
        pass