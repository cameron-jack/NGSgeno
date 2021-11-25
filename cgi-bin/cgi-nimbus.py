#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Aug 10 2020
@author: Bob Buckley & Cameron Jack - ANU Bioinformatics Consultancy JCSMR, ANU
@version: 0.12
@version_comment: 
@last_edit:
@edit_comment: 

Application Stage 1 (Nimbus) webpage.
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

stage1 = "/cgi-bin/cgi-nimbus.py"
stage2 = "/cgi-bin/cgi-echo.py"

html = """
<html>
<head>
<title>NGS Application!</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipelines (Mouse &amp; Custom)</h1>
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
    <h2>Mouse pipeline - Stage 1: Nimbus</h2>
    <p>Stage 1 writes picklist files for the Nimbus robot to tranfer material from 
    96-well ear-punch plates to a 384-well DNA sample plate for use in the Echo robot.
    New mouse projects will get a "run" folder created for them which will be today's date prefixed by "mouse_".
    </p>
    <form action="{stage}" method="get">
      <p>Please provide the following information:</p>
        {eps}
        <br>
        <label for="dnap">DNA (Echo) plate barcode:</label>
          <input type="text" id="dnap" name="dnap" size="10" required><br>
        <br>
        <input type="submit" value="Create Nimbus Picklist" style="border-color:red">
    </form>
    <hr>
    <form action="{stage}" method="get">
        <h2>Custom sample pipeline - Stage 1</h2>
        <p>To run the custom sample pipeline, you must first create a custom "run" folder (these MUST be prefixed with "custom_") within the NGSgeno 
        folder for your project. You then need to copy the required files there: Sample file - max 4 plates per file!(CSV), 
        Assay list file (CSV), Custom primer-plate layout (CSV)</p>
        <p>Set project folder (don't type the "custom_" prefix):</p>
        <label style="margin-left:25px;" for="projectDir"><em>custom_</em></label>
          <input type="text" id="projectDir" name="projectDir" size="40"><br>
        <p>Custom samples without guarded sample-name` or barcode are of what type:<br>
        <input type="radio" id="cust_custom" name="cust_type" value="custom" title="Unguarded samples/barcodes should be given custom guards cNNNc">
          <label for="cust_custom" title="Unguarded samples/barcodes should be given custom guards cNNNc">Custom names/barcodes</label><br>
        <input type="radio" id="cust_musterer" name="cust_type" value="musterer" 
            title="Unguarded samples/barcodes should be given Musterer guards mNNNm">
          <label for="cust_musterer" title="Unguarded samples/barcodes should be given Musterer guards mNNNm">Musterer mouse barcodes</label><br>
        <input type="radio" id="cust_rodentity" name="cust_type" value="rodentity" 
            title="Unguarded samples/barcodes should be given Rodentity guards MNNNM">
          <label for="cust_rodentity" title="Unguarded samples/barcodes should be given Rodentity guards MNNNM">Rodentity mouse barcodes</label><br>
        </p>

        <h3>Select custom sample, assay and primer plate files (CSV format)</h3>
        <p>These will be used instead of the corresponding Musterer information. These files MUST be in the run folder.</p>
         <label style="margin-left:25px;" for="customSamples">Manifest (samples/assays) file:</label>
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
    # Find optional sample directories - they have P*-EP.json files but not Nimbus*.csv files
    # If the user doesn't choose a directory, then a new sample directory will be created
    #jsonfiles = frozenset(map(os.path.dirname, glob.glob(os.path.join('*', 'P*-EP.json'))))
    #nimfiles = frozenset(map(os.path.dirname, glob.glob(os.path.join('*', 'Nimbus*.csv'))))
    #echofiles = frozenset(map(os.path.dirname, glob.glob(os.path.join('*', 'Echo*.csv')))) 
    #dirs = sorted(set([jsonfiles-nimfiles,jsonfiles-echofiles]), reverse=True)
    #if dirs:
    #    dx = '\n\t\t    '.join("<option value='{0}'>{0}</option>".format(d) for d in dirs)
    #    optx = """<label for="ngid">NGS Geno run Id:</label>
    #    <select name="ngid" id="ngid">
    #        <option value="">New MiSeq run</option>
    #        {}
    #    </select><br>
    #    <p>Leave this (choose 'New MiSeq run') if you are preparing a new MiSeq run. 
    #    The system will assign a new date-based run Id.
    #    <br></p>""".format(dx)  
    epfmt = '<label for="ep{0}">Ear-punch plate {0} barcode:</label>\n\t  <input type="text" id="ep{0}" name="ep{0}" size="6" /><br>'
    epx = '\n\t'.join(epfmt.format(i) for i in range(1,5))

    optform = ''
    global stage2
    wdirs = [d for d in os.listdir() if (d.startswith('mouse_') or d.startswith('custom_')) and os.path.isdir(d)]
    #wdirs = sorted(jsonfiles, reverse=True)
    if wdirs:
        dx = '\n\t\t    '.join("<option value='{0}'>{0}</option>".format(d) for d in wdirs)
        optmiseq = """<label for="ngid2">NGS Geno run Id:</label>
            <select name="existing" id="existing">
                <option value="">None</option>
                {}
            </select><br>""".format(dx)
        optform = """  <p> </p>
          <form action="{stage}" method="get">
              {optall}
              <input type="submit" value="Go to MiSeq run" />
          </form>""".format(stage=stage2, optall=optmiseq)
    print(html.format(options=optx, optany=optform, eps=epx, stage=stage2 ))
    return

if __name__=="__main__":
    try:
        main()
    except ProgramFail:
        pass