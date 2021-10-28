#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: June 2020
@author: Bob Buckley
@version: 0.12
@version_comment:
@last_edit:
@edit_comment:

Draw 96-well plates using Chart.js

This is run after Nimbus.py completes (and the Nimbusxxx.csv files were written), and generates files to report where everything has been placed.
"""

import os
import glob
import sys
import csv
import argparse
import collections

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

def makelabel(x, n):
    # print(x)
    pct = "%.2f"%(int(x[n])/int(x[1])*100)
    return '\\n'.join((x[0], ' '.join((x[1], x[n], pct+'%')))), x[1], pct

def main():
    parse = argparse.ArgumentParser(description="NGS Genotyping - Summary Report Program")
    parse.add_argument("barcodes", nargs='+', help="barcodes of Nimbus target plates")
    parse.add_argument('--custom', action='store_true', help='True if running the custom samples pipeline')
    args = parse.parse_args()
    
    if args.custom:
        Rec = collections.namedtuple('Rec', 'sno pbc well sid')
    else:
        Rec = collections.namedtuple('Rec', 'sno pbc well mbc')
    
    for pbc in args.barcodes:
        nimbusfn = "Nimbus{}.csv".format(pbc)
        with open(nimbusfn) as srcfd:
            src = csv.reader(srcfd, dialect="unix")
            next(src) # hdr = next(src) # drop header 
            data = [Rec(*x) for x in src]
        
        plates = sorted(set(x.pbc for x in data))
        if args.custom:
            pdata = dict((bc, [x for x in data if x.pbc==bc and x.sid!='0']) for bc in plates)
        else:
            pdata = dict((bc, [x for x in data if x.pbc==bc and x.mbc!='0']) for bc in plates)
        fs = []
        ch = []
        for figno, pn in enumerate(plates, start=1):
            if args.custom:
                fs.append(makefigure(figno, 'Custom sample plate '+pn))
            else:
                fs.append(makefigure(figno, "Ear-punch Plate "+pn))
            welz = set((ord('H')-ord(x.well[0]), int(x.well[1:])-1) for x in pdata[pn])
            if args.custom:
                argx = [("red", "Sample", [(None, ord('H')-ord(x.well[0]), int(x.well[1:])-1, 
                                            "Well: "+x.well+'; Sample barcode: '+x.sid) for x in pdata[pn]]),
                         ("white", "empty", [(None, r, c, chr(ord('H')-r)+str(c+1)+'; empty') 
                                             for c in range(12) for r in range(8) if (r,c) not in welz])]
            else:
                argx = [("red", "Sample", [(None, ord('H')-ord(x.well[0]), int(x.well[1:])-1, 
                                            "Well: "+x.well+'; Mouse barcode: '+x.mbc) for x in pdata[pn]]),
                         ("white", "empty", [(None, r, c, chr(ord('H')-r)+str(c+1)+'; empty') 
                                             for c in range(12) for r in range(8) if (r,c) not in welz])]
            ch.append(makeplate(figno, argx, title="Plate "+pn))
        
        htmlfn = "Stage1-P{}.html".format(pbc)
        with open(htmlfn, "wt") as dst:
            if args.custom:
                print(fmtReportCustom.format(figspace='\n'.join(fs), charts='\n'.join(ch)), file=dst)
            else:
                print(fmtReport.format(figspace='\n'.join(fs), charts='\n'.join(ch)), file=dst)

    return

if __name__=="__main__":
     main()

