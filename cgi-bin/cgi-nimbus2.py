#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Aug 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.7
@version_comment: all html reporting now working for custom pipeline
@last_edit: 2021-07-04
@edit_comment: Calls stage1report.py and stage2report/echo with --custom option if required

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
# Ugly hack to find things in ../bin
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(os.path.join(file.parents[1],'bin'))
import nimbus
from musterer import getPlate
import primercheck


nimbus2 = "/cgi-bin/cgi-nimbus2.py"
stage2  = "/cgi-bin/cgi-stage2.py"
i7i5_fn = os.path.join('..', 'library', 'i7i5_plate_layout_20200226.csv')

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
    <button onclick="window.history.back();">Back</button>
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
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 400px;">
  <h2>Nimbus transfer file(s)</h2>
  <form action="{stage1}" method="get">
    <p>Please copy Nimbus output files for the following DNA plate(s):<br>
    <pre>    {bcs}</pre>
    into the work folder: <code>NGSGeno\{ngid}</code><br>
    Note: Nimbus output files are CSV files whose name
    starts with <code>Echo_384_COC_</code><br><br>
    Click the "Continue" button below when done.</p>
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}"></input>
  </form>
  </div>
  {info2}
  {info1}
  <p>
    <button onclick="window.history.back();">Back</button>
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
  <form action="{stage2}" method="get">
    <p>This step creates an Echo picklist to transferring DNA from DNA plates: <code>{ebcs}</code></p>
    <p>The run needs {wellcount} wells in {noplates} PCR plates for separate assay samples.
    Please provide barcodes for the PCR plates.</p>
    {pcrs}
    <br>
    <hr>
    <p>This step also creates Echo picklists to transfer target primers, Mytaq and H<sub>2</sub>O
    into the PCR wells.</p>
    <p>Use the Echo robot to survey your primer plates. Upload the primer plate survey
    files from the Echo robot into the NGS Gentotyping sample folder for run <code>{ngid}</code>.</p>
          <label>Primer plate -</label><br>
          <label style="margin-left:25px;" for="psvy">survey file:</label>
          <input type="file" id="psvy" name="psvy" size="80" /><br>
    <p>Please provide a Mytaq & H<sub>2</sub>O plate with:
        <ol>
        <li>wells {taq} filled with Mytaq, and</li>
        <li>wells {h2o} filled with water.</li>
        </ol>
          
          <label for="taq">Mytaq #1 - plate barcode:</label>
          <input type="text" id="taq" name="taq" size="10" /><br><br>
          
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}" />
  </form>
  </div>
  {info}
  <p>
    <button onclick="window.history.back();">Back</button>
  </p>
</body>
</html>
""".strip()

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
  <form action="{stage2}" method="get">
    <p>Create Echo picklists for sample barcodes and Mytab.</p>
    <p>Create MiSeq file to describe samples for MiSeq run.</p>
    <p>Use the Echo robot to survey your barcode and Mytaq plates. Ensure you upload the
    files from the Echo robot into the NGS Gentotyping sample folder for run {ngid}.</p>
    <p>Please provide filenames, etc. for the following plates ...</p>
        <label style="margin-left:25px;" for="i7svy">i7i5 plate survey:</label>
        <input type="file" id="i7svy" name="i7svy" size="40" /><br>
          
    <p>Please provide a Mytaq & H<sub>2</sub>O plate with:
        <ol>
        <li>wells {taq} filled with Mytaq, and</li>
        <li>wells {h2o} filled with water.</li>
        </ol>         
          <label for="taq">Mytaq #2 - plate barcode:</label>
          <input type="text" id="taq" name="taq" size="10" /><br><br>
     <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}" />
  </form>
  </div>
  {info}
  <p>
    <button onclick="window.history.back();">Back</button>
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
  <div style="border: 2px solid black; border-radius: 10px; padding:10px; width: 400px;">
  <h2>Awaiting MiSeq results</h2>
  <p>When the MiSeq run completes and you have copied the MiSeq FASTQ files into 
  the "raw" folder for run {ngid} click the Continue button.</p>
  <p>Note: this step involves a lot of processing so expect to wait a while
  for it to complete.</p>
  <form action="{stage2}" method="get">
    <input type="submit" value="Continue" />
    <input style="display: none;" type="text" id="ngid" name="ngid" value="{ngid}" />
  </form>
  </div>
  {info}
  <p>
    <button onclick="window.history.back();">Back</button>
  </p>
</body>
</html>
""".strip()

html4 = """
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
    <button onclick="window.history.back();">Back</button>
  </p>
</body>
</html>
""".strip()

class ProgramFail(Exception):
    pass

def eppfn(bc):
    return "Plate{}-EP.json".format(bc) # ear-punch plate (filename) format

mustfix = re.compile(r'(,\s*"(mouseBarcode|plateId)":\s*)(\d+)')

def readCustomCSVtoJSON(input_fn, pids):
    """ read a custom plate and make sure the plate IDs match, then return as JSON """
    data = {}
    with open(input_fn, 'rt', newline=None) as f:
        in_csv = csv.reader(f, delimiter=',')
        for i, row in enumerate(in_csv):
            #print(row, file=sys.stderr)
            #print(file=sys.stderr)
            if i == 0:
                continue  # header
            plate = row[1]
            well = row[2]
            sampleBarcode = row[3]
            assays = [a.strip() for a in row[4:]]
            if plate not in pids:
                continue
            if plate not in data:
                data[plate] = {'plateId':plate, 'custom':True, 'sampleBarcode':sampleBarcode, 
                       'wells':[{'wellLocation':well, 'organism':{'sampleId':sampleBarcode, 
                                'sampleBarcode':sampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a} for a in assays]}}]}                  
            else:
                if well in data[plate]:
                    print('duplicate well number', well, 'in plate', plate, file=sys.stderr)
                    return
                data[plate]['wells'].append({'wellLocation':well, 'organism':{'sampleId':sampleBarcode, 
                        'sampleBarcode':sampleBarcode, 'assays':[{'assayMethod':'NGS','assayName':a} for a in assays]}})
           
    used_plates = set([p for p in data])
    for u in used_plates:
        assert u in pids, 'Plate '+str(u)+ ' not in requested plate ids' + '; '.join(map(str,pids))
    err = []
    for p in pids:
        if p not in data:
            err.append('Plate {p1} not found in {ifile}'.format(p1=p, ifile=input_fn))
    data = [(data[p],p) for p in sorted(data)]
    #print(data, file=sys.stderr)
    #with open('test_horizo_json.txt', 'wt') as out:
    #    print(data, file=out)
    return data, err

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
    
    global nimbus2
    # create new sample if the sample folder isn't present
    from datetime import date
    now = date.today() # get today's date
    nowstr = str(now.year)+str(now.month).zfill(2)+str(now.day).zfill(2)
    ngid = fields.getfirst("ngid") if "ngid" in fields else nowstr
    if not os.path.isdir(ngid): # only if dnap is also defined?
        os.mkdir(ngid)
    os.chdir(ngid) # change to working directory

    
    global htmlerr
    if "dnap" in fields:
        dnaBC = fields.getfirst("dnap")
        # get plate info from cache or from Musterer                               
        pids = [(i, fields.getfirst('ep'+str(i))) for i in range(1,5) if 'ep'+str(i) in fields]
        # Either call getPlate() to get plate data from Musterer, or if custom, get it from a file
        if 'customSamples' in fields:
            pxs, errx = readCustomCSVtoJSON(fields.getfirst('customSamples'),[pid for n,pid in pids])
            if fields.getfirst('customAssays') in (None, ''):
                errx.append('Custom assay list file required')
        else:
            # note - nimbus.getPlate_app calls mb.showerror() and app.setstatus and returns None for errors
            pxs = [getPlate(pid) for n, pid in pids]
            errx = [e for px, es in pxs for e in es]
        if len(pxs) == 0:
            print(htmlerr.format(errs="    <p>\n"+"No samples found, did you enter any plate barcodes?"+"\n    </p>"))
            return
        if errx: # any errors - report them
            errs = "    <p>\n"+"<br>\n    ".join(errx)+"\n    </p>"
            print(htmlerr.format(errs=errs))
            return
        if 'customSamples' in fields:
            plates_data = [(n, px[0]['plateId'], px[0]) for n,px in enumerate(pxs)]
        else:
            plates_data = [(n, pid, px) for (n, pid), (px, err) in zip(pids, pxs)]
        #print("Plates_data", plates_data, file=sys.stderr)
        
        # error checks - len(plates_data)>0
        custom_assay_fn = None
        if 'customAssays' in fields:
            custom_assay_fn = fields.getfirst('customAssays')
        custom_primer_fn = None
        if 'customPrimers' in fields:
            custom_primer_fn = fields.getfirst('customPrimers')
        
        #tgtfn = 'Nimbus'+dnaBC+'.csv'
        if 'customSamples' in fields:
            dnacnt, plist, unk = nimbus.nimbus_custom(dnaBC, plates_data, custom_assay_fn, custom_primer_fn)
        else:
            dnacnt, plist, unk = nimbus.nimbus(dnaBC, plates_data)

        # save "state" of custom resources to file
        if 'customSamples' in fields:
            with open('template_files.txt', 'wt') as outf:
                print('customSamples\t' + fields.getfirst('customSamples'), file=outf)
                if 'customAssays' in fields:
                    print('customAssays\t' + fields.getfirst('customAssays'), file=outf)
                if 'customPrimers' in fields:
                    print('customPrimers\t' + fields.getfirst('customPrimers'), file=outf)

        # acnt = len(plist)
        # infostr = ("Nimbus file for DNA plate ID {} uses {} wells.<br>".format(dnaBC, dnacnt) +  
        #           " {} wells will be used for assays in PCR plates.<br><br>".format(acnt))
        # info.append(infostr)
        if unk:
            info.append("Please ensure that all required primer description files are present in the work folder.")                    
            info.append("<b>Warning</b>: DNA plate {} - {} Musterer assays do not map to <b>recognised</b> target primer families:\n<pre>\n{}\n</pre>".format(dnaBC, len(unk), "\n".join(unk)))

    nfiles = glob.glob('Nimbus*.csv')
    # barcodes for Nimbus 384-well plate outputs
    nbc = frozenset(fn[6:-4] for fn in nfiles)
    # identify missing Stage1 files.
    missing = [bc for bc in nbc if not os.path.isfile("Stage1-P{}.html".format(bc))]
    if missing:
        # print("<pre>")
        # print("cwd =", os.getcwd())
        # print("nbc =", ' '.join(nbc))
        # print("missing =", ' '.join(missing))
        # should catch and report output
        if 'customSamples' in fields:
             subprocess.run(["python", os.path.join("..", "bin", "stage1report.py"), '--custom']+missing)
        else:
            subprocess.run(["python", os.path.join("..", "bin", "stage1report.py")]+missing)
        # print("</pre>")
    
    # find the Nimbus output files for each DNA plate.
    efiles = glob.glob('Echo_384_COC_00??_*.csv')
    # pick out the plate barcode from the filename
    ebc = frozenset(fn[:-4].split('_',5)[4] for fn in efiles)
    xbc  = nbc-ebc # any Nimbus plate files that are missing Nimbus run (output) files
    #print(xbc, nbc, ebc, file=sys.stderr)
    
    # now wait for the user to identify the Nimbus output file(s)
    if xbc:
        # Nimbus output files are not present
        global html1
        xinfo = ''.join("<p>{}</p>\n".format(x) for x in info)
        print(html1.format(bcs='\n    '.join(sorted(xbc)), stage1=nimbus2, ngid=ngid, info1=xinfo, info2=getinfo()))
        return
    
    def allfiles(fns):
        return all(os.path.isfile(fn) for fn in fns)
    
    #library = os.path.join("..", "library")
    
    template_files = collections.defaultdict(str)
    if os.path.exists('template_files.txt'):
        with open('template_files.txt', 'rt') as f:
            for line in f:
                cols = line.strip().split('\t')
                template_files[cols[0]] = cols[1]
        print('cgi-nimbus2:', [(k, template_files[k]) for k in template_files], file=sys.stderr)
    else:
        # set the most recent standard library
        template_files['assays'] = sorted(glob.glob(os.path.join('..', 'library', 'assay_list_20*.csv')), reverse=True)[0]
        template_files['primers'] = sorted(glob.glob(os.path.join('..', 'library', 'primer_layout*_20*.csv')), reverse=True)[0]
        template_files['ref'] = sorted(glob.glob(os.path.join('..', 'library', "reference_sequences_*.txt")), reverse=True)[0]

    prsvy = "primer-svy.csv"
    if 'pcr' in fields:
        # from htmlpl1 form - call the echo.py code
        #prs = [x for xs in zip(template_files['customPrimers'], fields.getlist('psvy')) for x in xs]
        if template_files['customPrimers'] != '':
            prs = [template_files['customPrimers']] + fields.getlist('psvy')
        else:
            prs = [template_files['primers']] + fields.getlist('psvy')
        cmd = ["python", os.path.join("..", "bin", "echovolume.py")]+prs+[prsvy]
        print('cgi-nimbus2:', cmd, file=sys.stderr)
        subprocess.run(cmd)
        # expect one only TAQ plate file
        cmd = ["python", os.path.join("..", "bin", "echo.py"), "--prim", prsvy, \
                        "--taq", fields.getfirst('taq'), "--pcr"] + \
                       fields.getlist('pcr')
        if template_files['customPrimers'] != '':
            cmd += ['--custom']
        cmd += ["--"] + sorted(ebc)
        subprocess.run(cmd)

    #TODO: i5's seem to be uneven in their use
    i7svy = "i7i5-svy.csv"
    if "i7svy" in fields:
        global i7i5_fn
        # from htmlpl2 form - call the echo.py code
        # should handle all the other options!!!
        # below handles multiple i7i5 plate surveys, but this is overkill
        #prs = [x for xs in zip(fields.getlist('i7i5'), fields.getlist('i7svy')) for x in xs]
        #cmd = ["python", os.path.join("..", "bin", "echovolume.py")]+prs+[i7svy]
        cmd = ['python', os.path.join('..', 'bin', 'echovolume.py')] + [i7i5_fn, fields.getfirst('i7svy'), i7svy]
        #print(cmd, file=sys.stderr)
        subprocess.run(cmd)
        cmd = ["python", os.path.join("..", "bin", "echo.py"), "--i7i5",
                        i7svy, "--taq", fields.getfirst('taq')]
        if template_files['customPrimers'] != '':
            cmd += ['--custom']
        cmd += ["--"]+sorted(ebc)
        subprocess.run(cmd) 

   
    # collect data for Stage 2 run
    # need to work out the number of PCR plates required.
    # Create input entries for the required number of plates.
    # also, we need names of the various primer, TAQ+water plates ...
    if not os.path.isfile("MiSeq-{}.csv".format(ngid)):
        global htmlpl1, htmlpl2
        # determine the number of PCR plates needed
        # maybe add reports - no. occupied DNA wells, no. assays per plate, ...
        
        prm = primercheck.PrimerLookup(fnmm=template_files['customAssays'], fnpp=template_files['customPrimers'])
        
        def platePrimers(bc):
            "reads the Stage1-P*.csv file - return tuple: barcode, set of primers per sample"
            fn = "Stage1-P{}.csv".format(bc)
            with open(fn) as srcfd:
                src = csv.reader(srcfd)
                hdr = next(src)
                prmidx = hdr.index("assayFamilies")
                res = [frozenset(p for f in map(str.strip, xs[prmidx].split(';')) if f for p in prm.pfdict.get(f, [])) for xs in src ]
            return bc, res
      
        # prepare the primer summary at this point - all Stage1-P*.csv files should be present the last time we get here.
        tabdata = [platePrimers(bc) for bc in sorted(ebc)]
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
        import echo
        if not allfiles([prsvy]):
            taq, h2o = echo.mytaq2(wellcnt, 1000, 300)
            print(htmlpl1.format(ngid=ngid, taq=sorted(set(taq)), h2o=sorted(set(h2o)), ebcs=' '.join(sorted(ebc)), wellcount=wellcnt, noplates=pc, pcrs=pcrx, info=getinfo(), stage2=nimbus2))
            return
        if not allfiles([i7svy]):
            taq, h2o = echo.mytaq2(wellcnt, 2000, 650)
            print(htmlpl2.format(ngid=ngid, taq=sorted(set(taq)), h2o=sorted(set(h2o)), info=getinfo(), stage2=nimbus2))
            return
        
        # the form calls echo.py (indirectly?) then displays the form for collecting MiSeq results
        # print(html2.format(wellcount=wellcnt, noplates=pc, pcrs=pcrx, stage2=nimbus2, ngid=ngid, info=''.join(info)+getinfo()))
        # return

    if not os.path.isdir("raw") and not os.path.isdir("RAW"):
        print('cgi-nimbus2:', 'No raw FASTQ directory', file=sys.stderr)
        global html3 # data analysis & reporting step
        print(html3.format(ngid=ngid, info=getinfo(), stage2=nimbus2))
        return
    
    if not os.path.isfile("Results.csv"):
        print('cgi-nimbus2:', 'Running stage3 FASTQ analysis', file=sys.stderr)
        res = subprocess.run([os.path.join("..", "bin", "stage3.bat"), "Stage3.csv"], capture_output=True)
        # capture output to log file
        with open("match.log", "wt") as dst:
            dst.write(res.stdout.decode('utf-8').replace('\r',''))
            if bool(res.stderr):
                dst.write("============ STDERR =============\n")
                dst.write(res.stderr.decode('utf-8').replace('\r',''))
                
    print(html4.format(ngid=ngid, info=getinfo(), stage2=nimbus2))
    return    
    

if __name__=="__main__":
    try:
        main()
    except ProgramFail:
        pass