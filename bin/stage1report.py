#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: June 2020
@author: Bob Buckley
@version: 0.15
@version_comment: New standardised stage file for custom/mouse samples. Error handling updated.
@last_edit: 2022-02-16
@edit_comment:

Draw 96-well plates using Chart.js

This is run after Nimbus.py completes (and the Nimbusxxx.csv files were written), and generates HTML files to report where everything has been placed.
Because it isn't part of the interface per-se it isn't in the cgi-bin folder. It is called by cgi-echo.py and operates on Stage1-PXXXX.csv files
        written by Nimbus.py.
"""

import os
import glob
import sys
import csv
import argparse
import collections
import file_io
from util import output_error

# Template for the HTML file
fmtReport = """
<!DOCTYPE html>
<html lang="en" dir="ltr" >
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type">
<title>NGS Genotyping - sample input report</title>
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/chart.js/dist/Chart.min.js"></script>
</head>
<body>
<h1>Mouse Ear-punch Plates Report</h1>
<p>Following are plate layouts for the Nimbus ear-punch plates.</p>
{!figspace!}
<p>You can "mouse over" a well to see which mouse ear-punch is in a well.</p>
<script>
{!charts!}
</script>
</body>
</html>
"""

fmtReportCustom = """
<!DOCTYPE html>
<html lang="en" dir="ltr" >
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type">
<title>NGS Genotyping - sample input report</title>
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/chart.js/dist/Chart.min.js"></script>
</head>
<body>
<h1>Custom Sample Plates Report</h1>
<p>Following are plate layouts for the Nimbus ear-punch plates.</p>
{!figspace!}
<p>You can "mouse over" a well to see which sample ID/barcode is in a well.</p>
<script>
{!charts!}
</script>
</body>
</html>
"""

# the following might go into a template file
fmtChart = """
var myChart{!figno!} = new Chart("fig{!figno!}", {
    type: 'scatter',
    data: {
        datasets: [{!data!},]
    },
    options: {
        title: {
                 display: true,
                 text: '{!title!}',
                 fontFamily: ['Georgia', 'Times', 'serif' ],
                 fontSize: 16, // padding: 35,
               },
        legend: { display: false },
        responsive: true,
        aspectRatio: 1.25, 
        scales: {
            xAxes: [{
                type: 'linear', offset: false, 
                ticks: { min: 0.5, max: 12.5, 
                         stepSize: 1,
                        callback: function(v,i,vs) { return (''+v).substr(-2)==".5"?'':v;}
                },
                // ticks: { beginAtZero: true },
                display: true, autoskip: false, 
                // scaleLabel: { labelString: 'raw read count', display: true, },
                gridLines: { display: false },
            }],
            yAxes: [{
                type: 'linear',
                ticks: {
                    min:0.5, max:8.5,
                    stepSize: 1,
                    callback: function(v, i, vs) { return 1<=i && i<=8?"ABCDEFGH".charAt(i-1):''; }, 
                },
                scaleLabel: { display: false },
                gridLines: { display: false },
            }],
        },
        tooltips: {
            callbacks: {
                label: function(item, data) {
                    //var values = item.xLabel+' '+item.yLabel+'%';
                    //return data.datasets[item.datasetIndex].labels[item.index]+' '+values;
                    return data.datasets[item.datasetIndex].labels[item.index];
                }
            }
        }
    }
});

"""
def makeplate(fig, datagens, title='???'):
    try:
        data = [(c, l, tuple(d)) for c, l, d in datagens]
        global fmtChart
        dsfmt = """{{
                backgroundColor: '{color}',
                pointRadius: 25, pointHoverRadius: 30,
                label: "{datalab}",
                data: [{pairs}],
                labels: [{labels}],
            }}"""
        dss = (dsfmt.format(color=c, datalab=l,
                            pairs=', '.join(("{x:%s, y:%s}" % (x[2]+1, x[1]+1) for x in d)), 
                            labels=', '.join('"%s"'%x[3] for x in d))
               for c, l, d in data)
        datasets = ",\n            ".join(dss)
        return fmtChart.format(figno=fig, title=title, data=datasets)
    except Exception as exc:
        output_error(exc, msg='Error in stage1report.makeplate')

fmtFigSpace = """
<figure><canvas id="fig{!figno!}" class="chartjs"></canvas>
<figcaption>Figure {!figno!}. {!caption!}</figcaption>
</figure>
"""

def makefigure(fig, caption=''):
    global fmtFigSpace
    return fmtFigSpace.format(figno=fig, caption=caption)

# use {!tag!} formatting - this makes it easier to edit javascript code
# in the templates. 
for s1, s2 in (('{', '{{'), ('{{!', '{'), ('}', '}}'), ('!}}', '}')):
    fmtReport, fmtReportCustom, fmtChart, fmtFigSpace = (x.replace(s1, s2) \
            for x in (fmtReport, fmtReportCustom, fmtChart, fmtFigSpace))

def getvalues(fn):
    try:
        fnx = fn.rsplit('/',1)[-1][:-4]
        d = {}
        fldnames = 'raw:', 'Pairs:', 'Joined:', 'Ambiguous:'
        with open(fn, "rt") as src:
            # print("fn=", fn)
            for rx in src:
                fld = rx.rstrip().split(None,2)
                if fld and fld[0] in fldnames:
                    d[fld[0]] = fld[1]
        return (fnx,)+tuple(d[k] for k in fldnames)
    except Exception as exc:
        output_error(exc, msg='Error in stage1report.getvalues')


def makelabel(x, n):
    # print(x)
    pct = "%.2f"%(int(x[n])/int(x[1])*100)
    return '\\n'.join((x[0], ' '.join((x[1], x[n], pct+'%')))), x[1], pct


def build_nimbus_reports(plate_bcs, is_custom=False):
    """ Programmatic interface - expects plate_bcs var to be guarded already
        Read from Stage1-PXXX.csv and produce an HTML report 
    """
    try:
        # sample number, guarded plate barcode, well, sampleid, guarded musterer barcode
        if is_custom:
            Rec = collections.namedtuple('Rec', ['sno', 'gpbc', 'well', 'gsid'])
        else:
            Rec = collections.namedtuple('Rec', ['sno', 'gpbc', 'well', 'gmbc'])

        for tgtp in plate_bcs:
            is_old_format = False
            nimbus_fn_old = f"Nimbus{file_io.unguard_pbc(tgtp)}.csv"
            if os.path.exists(nimbus_fn_old):
                nimbus_fn = nimbus_fn_old
                is_old_format = True
            else:
                nimbus_fn = f"Nimbus-{tgtp}.csv"
            records = []
            with open(nimbus_fn) as srcfd:
                src = csv.reader(srcfd, dialect="unix")
                next(src) # hdr = next(src) # drop header
                for row in src:
                    fields = [field.strip() for field in row]
                    if any(field == '' for field in fields):
                        print(f"Incomplete sample entry in {nimbus_fn}: {fields}", file=sys.stderr)
                        continue
                    #if not file_io.is_guarded(fields[3]) and not is_old_format:
                    #    print('WARNING: Possibly unguarded or incorrect sample barcode.', fields[3], file=sys.stderr)
                    try:
                        records.append(Rec(*fields))
                    except TypeError:
                        print('Bogus data line:',fields, file=sys.stderr)
                        exit(1)
                
            plate_set = sorted(set([x.gpbc for x in records]))
        
            if is_custom:
                pdata = dict((gpbc, [x for x in records if x.gpbc==gpbc and x.gsid!='0']) for gpbc in plate_set)
            else:
                pdata = dict((gpbc, [x for x in records if x.gpbc==gpbc and x.gmbc!='0']) for gpbc in plate_set)
            fs = []
            ch = []
            for figno, gpbc in enumerate(plate_set, start=1):
                if is_custom:
                    fs.append(makefigure(figno, 'Custom sample plate '+file_io.unguard_pbc(gpbc)))
                else:
                    fs.append(makefigure(figno, "Ear-punch Plate "+file_io.unguard_pbc(gpbc)))
                welz = set((ord('H')-ord(x.well[0]), int(x.well[1:])-1) for x in pdata[gpbc])
                if is_custom:
                    argx = [("red", "Sample", [(None, ord('H')-ord(x.well[0]), int(x.well[1:])-1, 
                                                "Well: "+x.well+'; Sample barcode: '+file_io.unguard(x.gsid)) for x in pdata[gpbc]]),
                             ("white", "empty", [(None, r, c, chr(ord('H')-r)+str(c+1)+'; empty')
                                                 for c in range(12) for r in range(8) if (r,c) not in welz])]
                else:
                    argx = [("red", "Sample", [(None, ord('H')-ord(x.well[0]), int(x.well[1:])-1,
                                                "Well: "+x.well+'; Mouse barcode: '+file_io.unguard(x.gmbc, silent=True)) for x in pdata[gpbc]]),
                             ("white", "empty", [(None, r, c, chr(ord('H')-r)+str(c+1)+'; empty')
                                                 for c in range(12) for r in range(8) if (r,c) not in welz])]

                ch.append(makeplate(figno, argx, title="Plate "+file_io.unguard_pbc(gpbc)))
        
            html_fn = f"Stage1-P{file_io.unguard_pbc(tgtp)}.html"
            with open(html_fn, "wt") as dst:
                if is_custom:
                    print(fmtReportCustom.format(figspace='\n'.join(fs), charts='\n'.join(ch)), file=dst)
                else:
                    print(fmtReport.format(figspace='\n'.join(fs), charts='\n'.join(ch)), file=dst)
    except Exception as exc:
        output_error(exc, msg='Error in stage1report.build_nimbus_reports')


def main():
    """ provide a command line interface. We'd prefer guarded plate barcodes, but we'll try to work with unguarded pbcs too """
    try:
        parse = argparse.ArgumentParser(description="NGS Genotyping - Summary Report Program")
        parse.add_argument("barcodes", nargs='+', help="barcodes of Nimbus target plates")
        parse.add_argument('--custom', action='store_true', help='True if running the custom samples pipeline')
        args = parse.parse_args()

        guarded_pbcs = []
        for pbc in args.barcodes:
            if file_io.is_guarded_pbc(pbc):
                guarded_pbcs.append(pbc)
            else:
                try:
                    guarded_pbcs.append(file_io.guard_pbc(pbc))
                except TypeException as t_ex:
                    print(t_ex, file=sys.stderr)
                    exit(1)
                except file_io.EmptyBarcodeError as ebe_ex:
                    print(ebe_ex, file=sys.stderr)
                    exit(1)
                except file_io.ExistingGuardError as ege_ex:
                    # We might want to check if they already have plate guards and let them through
                    print(ege_ex, file=sys.stderr)
                    exit(1)
        build_nimbus_reports(guarded_pbcs, is_custom=args.custom)
    except Exception as exc:
        output_error(exc, msg='Error in stage1report.main')


if __name__=="__main__":
     main()

