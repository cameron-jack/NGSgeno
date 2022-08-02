#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Aug 10 2020
@author: Bob Buckley & Cameron Jack - ANU Bioinformatics Consultancy JCSMR, ANU
@version: 0.15
@version_comment: Moved Nimbus related code and interface to cgi-nimbus
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
import file_io
from util import TemplateFiles, output_error, match_nimbus_to_echo_files
import nimbus


global htmlerr
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

global html_start 
html_start = """
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
        <input type="hidden" id="nimbus2" name="nimbus2" value="nimbus2">
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
        <p>Custom samples without guarded sample-name or barcode are of what type:<br>
        <input type="radio" id="cust_custom" name="cust_type" value="c" title="Unguarded samples/barcodes should be given custom guards cNNNc">
          <label for="cust_custom" title="Unguarded samples/barcodes should be given custom guards cNNNc">Custom names/barcodes</label><br>
        <input type="radio" id="cust_musterer" name="cust_type" value="m" 
            title="Unguarded samples/barcodes should be given Musterer guards mNNNm">
          <label for="cust_musterer" title="Unguarded samples/barcodes should be given Musterer guards mNNNm">Musterer mouse barcodes</label><br>
        <input type="radio" id="cust_rodentity" name="cust_type" value="M" 
            title="Unguarded samples/barcodes should be given Rodentity guards MNNNM">
          <label for="cust_rodentity" title="Unguarded samples/barcodes should be given Rodentity guards MNNNM">Rodentity mouse barcodes</label><br>
        </p>

        <h3>Select custom sample, assay and primer plate files (CSV format)</h3>
        <p>These will be used instead of the corresponding Musterer information. These files MUST be in the run folder.</p>
         <label style="margin-left:25px;" for="customSamples">Manifest (samples/assays) file:</label>
          <input type="file" id="customSamples" name="customSamples" size="80" accept=".csv"/><br><br>
         <label style="margin-left:25px;" for="customAssays">Assay list file:</label>
          <input type="file" id="customAssays" name="customAssays" size="80" accept=".csv"/><br><br>
        <label for="dnap">DNA (Echo) plate barcode:</label>
          <input type="text" id="dnap_cust" name="dnap_cust" size="10" required/><br>
        <br>
        <input type="submit" value="Create Nimbus Picklist" style="border-color:red"/>
        <input type="hidden" id="nimbus2" name="nimbus2" value="nimbus2"/>
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

global html_start2
html_start2 = """
<html>
<head>
<title>NGS Application!</title>
<meta charset="utf-8">
  <style>
    {{
        box-sizing: border-box;
    }}
    /* Set additional styling options for the columns*/
    .column {{
    float: left;
    width: 50%;
    }}

    .row:after {{
    content: "";
    display: table;
    clear: both;
    }}
  </style>
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipelines (Mouse &amp; Custom)</h1>
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 650px;">
    <h2>Current run folder: </h2> 
    <p>Create new run folder - don't type the "run_" prefix:
        <label style="margin-left:25px;" for="runDir"><em>run_</em></label>
          <input type="text" id="runDir" name="runDir" size="40"><br>
    </p>
    <p>Load existing run folder (either separate page, or Javascript drop-down choice list)</p>
    <hr>
    <div class="row">
   
      <div class="column border-right">
        <form action="{stage}" method="get">
          <h2>Mouse pipeline</h2>
          {eps}
          <br>
          <label for="dnap">DNA (Echo) plate barcode:</label>
            <input type="text" id="dnap" name="dnap" size="10" required><br>
          <br>
          <input type="submit" value="Add" style="border-color:red">
          <input type="hidden" id="nimbus2" name="nimbus2" value="nimbus2">
        </form>
      </div>
      
      <div class="column">
        <h2>Custom pipeline</h2>
        <p>Custom samples without guarded sample-name or barcode are of what type:<br>
        <input type="radio" id="cust_custom" name="cust_type" value="c" title="Unguarded samples/barcodes should be given custom guards cNNNc">
          <label for="cust_custom" title="Unguarded samples/barcodes should be given custom guards cNNNc">Custom names/barcodes</label><br>
        <input type="radio" id="cust_musterer" name="cust_type" value="m" 
            title="Unguarded samples/barcodes should be given Musterer guards mNNNm">
          <label for="cust_musterer" title="Unguarded samples/barcodes should be given Musterer guards mNNNm">Musterer mouse barcodes</label><br>
        <input type="radio" id="cust_rodentity" name="cust_type" value="M" 
            title="Unguarded samples/barcodes should be given Rodentity guards MNNNM">
          <label for="cust_rodentity" title="Unguarded samples/barcodes should be given Rodentity guards MNNNM">Rodentity mouse barcodes</label><br>
        </p>
        <p>
        <label style="margin-left:25px;" for="customSamples">Manifest (samples/assays) file:</label>
          <input type="file" id="customSamples" name="customSamples" size="80" accept=".csv"/><br><br>
         <label style="margin-left:25px;" for="customAssays">Assay list file:</label>
          <input type="file" id="customAssays" name="customAssays" size="80" accept=".csv"/><br><br>
        <label for="dnap">DNA (Echo) plate barcode:</label>
          <input type="text" id="dnap_cust" name="dnap_cust" size="10" required/><br>
        </p>
      </div>
    </div>
    
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


html_nimbus = """
<html>
<head> 
<title>NGS Genotyping</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline</h1>
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
  <h2>Nimbus transfer file(s)</h2>
  <form action="{stage}" method="get">
    <p>Please copy Nimbus output files for the following DNA plate(s):<br>
    <pre>    {bcs}</pre>
    into the work folder: <code>NGSGeno\{ngid}</code><br>
    Note: Nimbus output files are CSV files whose name
    starts with <code>Echo_384_COC_</code><br><br>
    Click the "Continue" button below when done.</p>
    <input type="hidden" id="existing" name="existing" value="{existingDir}"/>
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}"></input>
  </form>
  </div>
  <p>Valid <b>primer usage</b> has been written to <b>"primerlist_valid.csv"</b> and unavailable primers have been written to "primerlist_invalid.csv"</p>
  {info2}
  {info1}
  <p>
    <button onclick="window.history.back();">Back to previous page</button>
  </p>
  <p>
    <button onclick="location.href='http://localhost:{port}/cgi-bin/cgi-nimbus.py';" value="Start again from first pipeline page">Back to start of pipeline</button>
  </p>
</body>
</html>
""".strip()


html_primers = """
<html>
<head> 
<title>NGS Genotyping</title>
<meta charset="utf-8">
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline</h1>
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
  <h2>Nimbus transfer file(s)</h2>
  <form action="{stage}" method="get">
    <p>The following are the number of uses required for each valid primer requested by the given samples:</p>
    {primer_frequency_list}
    <p>Please prepare {numplates} primer {plate_plural}. Place them in the project run folder and then select them below.
      {primer_plate_selection}
      <label style="margin-left:25px;" for="customPrimers1">Custom primer-plate1 layout:</label>
        <input type="file" id="customPrimers1" name="customPrimers1" size="80" accept=".csv" /><br><br>
      
    </p>
</body>
</html>
""".strip()


def main():
    fields = cgi.FieldStorage()
    errs, info  = "", []
    cgitb.enable(display=1, logdir=".")
    print("Content-Type: text/html")    # HTML is following
    print()                             # blank line, end of headers
    
    if "nimbus2" not in fields and "getPrimers" not in fields:
        optx = ''
        epfmt = '<label for="ep{0}">Ear-punch plate {0} barcode:</label>\n\t  <input type="text" id="ep{0}" name="ep{0}" size="6" /><br>'
        epx = '\n\t'.join(epfmt.format(i) for i in range(1,5))

        optform = ''
        
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
        print(html_start.format(options=optx, optany=optform, eps=epx, stage=stage1))
    
    elif "nimbus2" in fields:      
        # create new sample if the sample folder isn't present
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

        if not os.path.isdir(ngid): # only if dnap is also defined?
            os.mkdir(ngid)
        os.chdir(ngid) # change to working directory

        ### Everything below this line happens in the run folder ###
        templates = TemplateFiles()  # load template_file.txt if it exists
        templates.set_mouse_files()
        templates.set_layout_files()
        templates.set_survey_files(save=True)

        global htmlerr
        dnaBC = ''
        if "dnap" in fields:
            dnaBC = fields.getfirst("dnap")
        elif "dnap_cust" in fields:
            dnaBC = fields.getfirst("dnap_cust")
        elif "existing" in fields:
            pass
        else:
            print(htmlerr.format(errs="    <p>\n"+"No samples found, did you enter any plate barcodes?"+"\n    </p>", port=port))
            return
        if dnaBC:
            dnaBC = dnaBC.strip()
            if not file_io.is_guarded_pbc(dnaBC):
                dnaBC = file_io.guard_pbc(dnaBC)

            cust_manifest = None
            cust_assays = None
            cust_ref = None
            # check if we have custom fields from interface
            if 'customSamples' in fields:
                cust_manifest = fields.getfirst('customSamples')
            if 'customAssays' in fields:
                cust_assays = fields.getfirst('customAssays')
            if 'customReference' in fields:
                cust_ref = fields.getfirst('customReference')
            templates.set_custom_files(manifest=cust_manifest, assays=cust_assays, ref=cust_ref, save=True)
            if 'customPrimers1' in fields:
                templates.set_layout_files(primer1=fields.getfirst('customPrimers1'), save=True)
            if 'customPrimers2' in fields:
                templates.set_layout_files(primer2=fields.getfirst('customPrimers2'), save=True)
                
            # Either call getPlate() to get plate data from Musterer, or if custom, get it from a file
            if 'cust_manifest' in templates.files:
                # get plate info from file
                cust_file_contents = 'custom'
                if 'cust_type' in fields:
                    cust_file_contents = fields.getvalue('cust_type')
                pxs, pids, errx = nimbus.readCustomCSVtoJSON(templates.files['customSamples'], cust_file_contents)
                pids = [p if file_io.is_guarded_pbc(p) else file_io.guard_pbc(p) for p in pids]
                if len(pxs) == 0:
                    print(htmlerr.format(errs="    <p>\n"+"No samples found, did you enter any plate barcodes?"+"\n    </p>", port=port))
                    return
                #print('pxs:', pxs, file=sys.stdout)
                #print('pids', pids, file=sys.stdout)
                plates_data = [(n, px[0]['plateId'], px[0]) for n,px in enumerate(pxs)]

            elif 'ep1' in fields:
                # note - nimbus.getPlate_app calls mb.showerror() and app.setstatus and returns None for errors
                # get plate info from input form and lookup in Musterer or cache
                pids = []
                for i in range(1,5):
                    if 'ep'+str(i) not in fields:
                        continue
                    p = fields.getfirst('ep'+str(i))
                    if not file_io.is_guarded_pbc(p):
                        p = file_io.guard_pbc(p)
                    pids.append((i,p))
                pxs = [get_plate(pid) for n, pid in pids]  # guarded plateIDs
                errx = [e for px, es in pxs for e in es]
                if len(pxs) == 0:
                    print(htmlerr.format(errs="    <p>\n"+"No samples found, did you enter any plate barcodes?"+"\n    </p>", port=port))
                    return
                plates_data = [(n, pid, px) for (n, pid), (px, err) in zip(pids, pxs)]

            if errx: # any errors - report them
                errs = "    <p>\n"+"<br>\n    ".join(errx)+"\n    </p>"
                print(htmlerr.format(errs=errs, port=port))
                return

            # load assay list
            if 'cust_assay_list' in templates.files:
                assay_info = validate_assays.Assays(templates.files['cust_assay_list'])
            else:
                assay_info = validate_assays.Assays(templates.files['mouse_assay_list'])
            valid_assays, invalid_assays = assay_info.validate_assays(plates_data)
            info.append(f"<b>Warning</b>: DNA plate {dnaBC} - {len(invalid_assays)} Musterer assays do not map to <b>recognised</b> "+\
                f"target primer families:{nl}<pre>{nl}{nl.join(invalid_assays)}{nl}</pre>")
            info.append(f"<b>The following assays will be used:</b>{nl}<pre>{nl.join(valid_assays)}{nl}</pre>")

            # prepare Nimbus files
            if 'cust_manifest' in templates.files:
                dnacnt, plist, unk = nimbus.nimbus_custom(dnaBC, plates_data, assay_info)
            else:
                dnacnt, plist, unk = nimbus.nimbus(dnaBC, plates_data, assay_info)
            #print("Plates_data", plates_data, file=sys.stdout)
            
        # now wait for the user to identify the Nimbus output file(s)
        xbc = match_nimbus_to_echo_files(templates)
        if xbc:
            # Nimbus output files are not present
            global html_nimbus
            xinfo = ''.join("<p>{}</p>\n".format(x) for x in info)
            print(html_nimbus.format(bcs='\n    '.join(sorted(xbc)), stage=stage2, ngid=ngid, info1=xinfo, info2=getinfo(ngid), existingDir=ngid, port=port))
            return
    elif "getPrimers" in fields:
        # Run validate_primers and show the user the primer lists
        # determine the number of PCR plates needed
        # maybe add reports - no. occupied DNA wells, no. assays per plate, ...
        
        # load assay list - might be already loaded, but should be fine
        if 'customAssays' in templates.files:
            assay_info = validate_assays.Assays(templates.files['cust_assay_list'])
        else:
            assay_info = validate_assays.Assays(templates.files['mouse_assay_list'])
        
        def platePrimers(bc, assay_info):
            """ reads the Stage1-P*.csv file - return tuple: barcode, set of primers per sample """
            fn = f"Stage1-P{bc}.csv"
            with open(fn) as srcfd:
                src = csv.reader(srcfd)  # open each Stage1 file
                hdr = next(src)
                prmidx = hdr.index("assayFamilies")  # get the index for the assayFamilies column
                p_list = [p.strip() for xs in src for p in xs[prmidx].split(';')]  # return the assays that match the families
                res = [frozenset(p for p in p_list if p in assay_info.all_assays)]  # return the matching primer list
                #res = [frozenset(p for p in xs in src for f in map(str.strip, xs[prmidx].split(';')) if p in assay_info.all_assays)]
            return bc, res
      
        # prepare the primer summary at this point - all Stage1-P*.csv files should be present the last time we get here.
        tabdata = [platePrimers(bc, assay_info) for bc in sorted(ebc)]
        cx = collections.Counter(p for bc, pss in tabdata for ps in pss for p in ps)  # counter the primer usage
        with open("PrimerUse.csv", "wt", newline='') as dstfd:
            dst = csv.writer(dstfd)
            dst.writerow(["Primer", "sample count"])
            dst.writerows(cx.most_common())
        info.append("  <table class='tab2'><tr><th>DNA Plate</th><th>Samples</th><th>Primers</th><th>PCR wells</th><tr>\n")
        sampcnt, wellcnt, prmall = 0, 0, set([])
        for bc, pss in tabdata:
            rc = sum(len(ps) for ps in pss)
            prmset = frozenset(p for ps in pss for p in ps)
            info.append(f"    <tr><th>{bc}</th><td style='text-align: right;'>{len(pss)}</td><td style='text-align: right;'>{len(prmset)}</td><td style='text-align: right;'>{rc}</td></tr>\n")
            sampcnt += len(pss)
            wellcnt += rc
            prmall |= prmset
        if len(tabdata)>1:
            info.append(f"    <tr><th>Total</th><td style='text-align: right;'><i>{sampcnt}</i></td><td style='text-align: right;'><i>{len(prmall)}</i></td><td style='text-align: right;'><i>{wellcnt}</i></td></tr>\n")
        info.append("  </table>\n")
        pc = (wellcnt+(384-3-1))//(384-3) # leave 3 empty wells in each plate
        # should also report unknown Assays
        pcrfmt = '<label for="pcr{0}">PCR plate {0} barcode:</label>\n\t  <input type="text" id="pcr{0}" name="pcr" size="8" /><br>'
        pcrx = '\n\t'.join(pcrfmt.format(i) for i in range(1,pc+1))
    else:
        # should never get here
        print('Error in cgi-nimbus.py - no selection made', file=sys.stdout)

if __name__=="__main__":
    try:
        main()
    except Exception as exc:
        output_error(exc, msg='cgi-nimbus.py failed')