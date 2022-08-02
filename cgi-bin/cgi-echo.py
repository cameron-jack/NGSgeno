#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Aug 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.15
@version_comment: Nimbus related code moved to cgi-nimbus.py. 
    Definitions and utility functions moved to cgi_util.py. File related code moved to file_io.py
@last_edit: 2022-01-24
@edit_comment: Rebalanced the stages to improve readability.

Application Webpage.
Nimbus picklists are already created - now gather Nimbus files
or proceed to creating Echo picklists (when all the Nimbus output is present).
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

from cgi_util import port, nl, stage1, stage2, stage3, getinfo, match_nimbus_to_echo_files

# Ugly hack to find things in ../bin
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(os.path.join(file.parents[1],'bin'))
import nimbus
from musterer import getPlate
import validate_assays
import file_io
from util import TemplateFiles, ProgramFail


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

html1 = """
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
  <p>Valid <b>primer usage</b> has been written to <b>"primerlist_valid.csv"</b> and unavailable 
  primers have been written to "primerlist_invalid.csv"</p>
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

htmlpl1 = """
<html>
<head> 
<title>NGS Application!</title>
<meta charset="utf-8">
  <style>
    table {{ border-collapse: collapse; }}
    table, th, td {{ border: 2px solid gray; }}
    table {{ margin: 25px; }}
    th, td {{ padding: 2px 5px; }}
  </style>
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline</h1> 
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
  <h2>Echo Picklists - part 1</h2>
  <form action="{stage}" method="get">
    <p>This step creates an Echo picklist to transfer DNA from DNA plates to PCR reaction plates: <code>{ebcs}</code></p>
    <p>The run needs {wellcount} wells in {noplates} PCR plates for separate assay samples.
    Please provide barcodes for the PCR plates.</p>
    {pcrs}
    <br>
    <hr>
    <p>This step also creates Echo picklists to transfer target primers, Mytaq and H<sub>2</sub>O
    into the PCR wells.</p>
    {primer_survey_html}
    <p>Please provide {tc} Mytaq & H<sub>2</sub>O plates barcodes with:
        <ol>
        <li>wells {taq} filled with Mytaq, and</li>
        <li>wells {h2o} filled with water.</li>
        </ol>
        {taqs}
        <br>
    <input type="hidden" id="existing" name="existing" value="{existingDir}"/>
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}" />
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


#<p>OPTIONAL (Leave this blank except in special circumstances): Use the Echo robot to survey your barcode and Mytaq plates. Ensure you upload the
#    files from the Echo robot into the NGS Gentotyping sample folder for run {ngid}.</p>
#    
#    <label style="margin-left:25px;" for="i7svy">i7i5 plate survey:</label>
#      <input type="file" id="i7svy" name="i7svy" size="40" /><br>

htmlpl2 = """
<html>
<head> 
<title>NGS Application!</title>
<meta charset="utf-8">
  <style>
    table {{ border-collapse: collapse; }}
    table, th, td {{ border: 2px solid gray; }}
    table {{ margin: 25px; }}
    th, td {{ padding: 2px 5px; }}
  </style>
</head>
<body style="padding: 20px;">
  <h1>NGS Genotyping pipeline</h1> 
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 600px;">
  <h2>Echo picklists - part 2</h2>
  <form action="{stage}" method="get">
    <p>Create Echo picklists for sample barcodes and Mytab.</p>
    <p>Create MiSeq file to describe samples for MiSeq run.</p>
    
    <input type="hidden" id="i7stage" name="i7stage" value="i7stage" />
          
    <p>Please provide {tc} Mytaq & H<sub>2</sub>O plates with:
        <ol>
        <li>wells {taq} filled with Mytaq, and</li>
        <li>wells {h2o} filled with water.</li>
        </ol>         
        {taqs}
        <br>
     <input type="hidden" id="existing" name="existing" value="{existingDir}"/>
     <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}" />
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

html3 = """
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


def eppfn(bc):
    return "Plate-{}-EP.json".format(bc) # ear-punch plate (filename) format

def eppfn_old(bc):
    return "Plate{}-EP.json".format(bc) # ear-punch plate (filename) format, unguarded name (older)

mustfix = re.compile(r'(,\s*"(mouseBarcode|plateId)":\s*)(\d+)')


def main(): 
    fields = cgi.FieldStorage()
    errs, info  = "", []
    cgitb.enable(display=1, logdir=".")
    print("Content-Type: text/html")    # HTML is following
    print()                             # blank line, end of headers
    
    def allfiles(fns):
        return all(os.path.isfile(fn) for fn in fns)

    xbc = match_nimbus_to_echo_files()
    if xbc:
        # Nimbus output files are not present, output primer lists if available
        global html_nimbus
        if os.path.exists('primerlist_valid'):
            with open('primerlist_valid','rt') as f:
                info.append('Valid primers requested by samples:')
                info = info + [l.strip() for l in f]
        if os.path.exists('primerlist_invalid'):
            with open('primerlist_invalid','rt') as f:
                if len(info) > 0:
                    info.append('')
                info.append('Invalid primers requested by samples:')
                info = info = [l.strip() for l in f]
        xinfo = ''.join("<p>{}</p>\n".format(x) for x in info)
        print(html_nimbus.format(bcs='\n    '.join(sorted(xbc)), stage=stage2, ngid=ngid, info1=xinfo, info2=getinfo(ngid), existingDir=ngid, port=port))
        return

    prsvy_fn = "primer-svy.csv"
    if 'pcr' in fields:
        # from htmlpl1 form - call the echo.py code
        #prs = [x for xs in zip(templates.files['customPrimers'], fields.getlist('psvy')) for x in xs]
        if 'psvy' in fields:
            templates.files['primer_survey'] = fields.getlist('psvy')
        prs = [fn for fn in templates.files['primer_survey']]
        ### Combine primer layouts and surveys here ###

        if templates.files['cust_primer_layout'] != '':
            prs = [templates.files['cust_primer_layout']] + psvy_fn
        else:
            prs = [templates.files['mouse_primer_layout']] + psvy_fn
        cmd = ["python", os.path.join("..", "bin", "echovolume.py")]+prs+[prsvy_fn]
        print('cgi-echo:', cmd, file=sys.stdout)
        subprocess.run(cmd)
        # expect one only TAQ plate file
        cmd = ["python", os.path.join("..", "bin", "echo_primer.py"), "--prim", prsvy_fn, \
                        "--taq", fields.getfirst('taq'), "--pcr"] + \
                       fields.getlist('pcr')
        if templates.files['customPrimers'] != '':
            cmd += ['--custom']
        cmd += ["--"] + sorted(ebc)
        subprocess.run(cmd)

    #TODO: i5's seem to be uneven in their use
    i7svy = "i7i5-svy.csv"
    if "i7stage" in fields:
        # below handles multiple i7i5 plate surveys, but this is overkill
        #prs = [x for xs in zip(fields.getlist('i7i5'), fields.getlist('i7svy')) for x in xs]
        #cmd = ["python", os.path.join("..", "bin", "echovolume.py")]+prs+[i7svy]
        if "i7svy" in fields:
            i7i5_survey_fn = fields.getfirst('i7sv')
        else:
            i7i5_survey_fn = templates.files['i7i5_survey']
        cmd = ['python', os.path.join('..', 'bin', 'echovolume.py')] + [i7i5_fn, i7i5_survey_fn, i7svy]
        #print(cmd, file=sys.stdout)
        subprocess.run(cmd)
        cmd = ["python", os.path.join("..", "bin", "echo_barcode.py"), "--i7i5",
                        i7svy, "--taq"] + fields.getlist('taq')
        if templates.files['customPrimers'] != '':
            cmd += ['--custom']
        cmd += ["--"]+sorted(ebc)
        subprocess.run(cmd) 

    # TODO: shift this later - we don't know the primer layout at the start of the pipeline, but we do need the PrimerUse.csv file
    # collect data for Stage 2 run
    # need to work out the number of PCR plates required.
    # Create input entries for the required number of plates.
    # also, we need names of the various primer, TAQ+water plates ...
    if not os.path.isfile("MiSeq-{}.csv".format(ngid)):
        global htmlpl1, htmlpl2
        # determine the number of PCR plates needed
        # maybe add reports - no. occupied DNA wells, no. assays per plate, ...
        
        # load assay list - might be already loaded, but should be fine
        if 'customAssays' in templates.files:
            assay_info = validate_assays.Assays(templates.files['cust_assay_list'])
        else:
            assay_info = validate_assays.Assays(templates.files['mouse_assay_list'])
        
        def platePrimers(bc, assay_info):
            "reads the Stage1-P*.csv file - return tuple: barcode, set of primers per sample"
            fn = f"Stage1-P{bc}.csv"
            with open(fn) as srcfd:
                src = csv.reader(srcfd)
                hdr = next(src)
                prmidx = hdr.index("assayFamilies")
                p_list = [p.strip() for xs in src for p in xs[prmidx].split(';')]
                res = [frozenset(p for p in p_list if p in assay_info.all_assays)]
                #res = [frozenset(p for p in xs in src for f in map(str.strip, xs[prmidx].split(';')) if p in assay_info.all_assays)]
            return bc, res
      
        # prepare the primer summary at this point - all Stage1-P*.csv files should be present the last time we get here.
        tabdata = [platePrimers(bc, assay_info) for bc in sorted(ebc)]
        cx = collections.Counter(p for bc, pss in tabdata for ps in pss for p in ps)
        with open("PrimerUse.csv", "wt", newline='') as dstfd:
            dst = csv.writer(dstfd)
            dst.writerow(["Primer", "sample count"])
            dst.writerows(cx.most_common())
        info.append("  <table class='tab2'><tr><th>DNA Plate</th><th>Samples</th><th>Primers</th><th>PCR wells</th><tr>\n")
        sampcnt, wellcnt, prmall = 0, 0, set([])
        for bc, pss in tabdata:
            rc = sum(len(ps) for ps in pss)
            prmset = frozenset(p for ps in pss for p in ps)
            info.append("    <tr><th>{}</th><td style='text-align: right;'>{}</td><td style='text-align: right;'>{}</td><td style='text-align: right;'>{}</td></tr>\n".format(bc, len(pss), len(prmset), rc))
            sampcnt += len(pss)
            wellcnt += rc
            prmall |= prmset
        if len(tabdata)>1:
            info.append("    <tr><th>Total</th><td style='text-align: right;'><i>{}</i></td><td style='text-align: right;'><i>{}</i></td><td style='text-align: right;'><i>{}</i></td></tr>\n".format(sampcnt, len(prmall), wellcnt))
        info.append("  </table>\n")
        pc = (wellcnt+(384-3-1))//(384-3) # leave 3 empty wells in each plate
        # should also report unknown Assays
        
        pcrfmt = '<label for="pcr{0}">PCR plate {0} barcode:</label>\n\t  <input type="text" id="pcr{0}" name="pcr" size="8" /><br>'
        pcrx = '\n\t'.join(pcrfmt.format(i) for i in range(1,pc+1))
    
        #TODO: Need to be able detect and use multiple Taq/water plates
        import echo_primer
        if not allfiles([prsvy_fn]):
            taq, h2o = echo_primer.mytaq2(wellcnt, 1000, 300)
            tc = max(max([p for w,p in h2o]),max([p for t,p in taq]))
            #print(len(taq), len(h2o), file=sys.stdout)
            taq_set = sorted(set([w for w,p in taq]))
            h2o_set = sorted(set([w for w,p in h2o]))
            taqfmt = '<label for="taq">Mytaq/Water plate {0} barcode:</label>\n\t  <input type="text" id="taq" name="taq" size="8" /><br>'
            taqx ='\n\t'.join(taqfmt.format(i) for i in range(1,tc+1))
            primer_survey_html=''
            if template_files['customPrimers']:
                primer_survey_html="""
                <p>Custom workflows: Use the Echo robot to survey your primer plates. Upload the primer plate survey
    files from the Echo robot into the NGS Gentotyping sample folder and select them here.</p>
                <label>Primer plate -</label><br>
                <label style="margin-left:25px;" for="psvy">survey file(s):</label>
                 <input type="file" multiple="true" id="psvy" name="psvy" size="80" /><br>
                """
            print(htmlpl1.format(ngid=ngid, taq=taq_set, h2o=h2o_set, ebcs=' '.join(sorted(ebc)), 
                    wellcount=wellcnt, noplates=pc, primer_survey_html=primer_survey_html, pcrs=pcrx, tc=tc, 
                    taqs=taqx, info=getinfo(), stage=stage2, existingDir=ngid,port=port))
            return
        if not allfiles([i7svy]):
            taq, h2o = echo_primer.mytaq2(wellcnt, 2000, 650)
            tc = max(max([p for w,p in h2o]),max([p for t,p in taq]))
            taq_set = sorted(set([w for w,p in taq]))
            h2o_set = sorted(set([w for w,p in h2o]))
            taqfmt = '<label for="taq">Mytaq/Water plate {0} barcode:</label>\n\t  <input type="text" id="taq" name="taq" size="8" /><br>'
            taqx ='\n\t'.join(taqfmt.format(i) for i in range(1,tc+1))
            print(htmlpl2.format(ngid=ngid, taq=taq_set, h2o=h2o_set, tc=tc, taqs=taqx, info=getinfo(), stage=stage2, existingDir=ngid,port=port))
            return
        
    global stage3
    if not os.path.isdir("raw") and not os.path.isdir("RAW"):
        print('cgi-echo:', 'No raw FASTQ directory', file=sys.stdout)
        global html3 # data analysis & reporting step
        print(html3.format(ngid=ngid, info=getinfo(), stage=stage3, existingDir=ngid, port=port))
        return
    else:
        if not os.path.isfile("Results.csv"): # we've yet to run the analysis stage
            print(html3.format(ngid=ngid, info=getinfo(), stage=stage3, existingDir=ngid, port=port))
            return
        else:
            print(html_results.format(ngid=ngid, info=getinfo(), stage=stage2, existingDir=ngid, port=port))
            return
    return           
       

if __name__=="__main__":
    try:
        main()
    except ProgramFail:
        pass