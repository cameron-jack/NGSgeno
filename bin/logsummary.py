#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created Apr 2020
@author: Bob Buckley
@version: 0.10
@version_comment: 
@last_edit:
@edit_comment:

Program to summarise merge log file info.
Uses the log files in the mclean directory created by the merge.sh script. 
It assumes that paired end read were "cleaned".
"""
defdir = "."
defrpt = "Summary.htm"
defcsv = "Summary.csv"

import os
import glob
import sys
import csv
import argparse
#import collections

# Template for the HTML file
fmtReport = """
<!DOCTYPE html>
<html lang="en" dir="ltr" >
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type">
<title>NGS genotyping - MiSeq base report</title>
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/chart.js@2.9.3/dist/Chart.min.js"></script>
</head>
<body>
<h1>Base Report</h1>
<p>The following figure shows trimmed, merged and ambiguous read rates relative to the number of
raw reads for each sample reported.</p>
{!figspace!}
<p>You can click on items in the legend to hide or reveal them.</p>
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
                 fontSize: 18, // padding: 24,
               },
        legend: { display: true },
        scales: {
            xAxes: [{
                type: 'logarithmic',
                // display: true,
                scaleLabel: { labelString: 'raw read count', display: true, },
                gridLines: { display: true },
            }],
            yAxes: [{
                type: 'linear',
                ticks: {
                    beginAtZero: true,
                    callback: function(v, i, vs) { return v+'%'; }, 
                },
                scaleLabel: { display: false, labelString: 'percent' },
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
def makemchart(fig, datagens, title='???'):
    data = [(c, l, tuple(d)) for c, l, d in datagens]
    global fmtChart
    dsfmt = """{{
            backgroundColor: '{color}',
            label: "{datalab}",
            data: [{pairs}],
            labels: [{labels}],
        }}"""
    dss = (dsfmt.format(color=c, datalab=l, 
                        pairs=', '.join(("{x:%s, y:%s}" % (x[1], x[2]) for x in d)), 
                        labels=', '.join('"%s"'%x[0] for x in d))
           for c, l, d in data)
    datasets = ",\n            ".join(dss)
    return fmtChart.format(figno=fig, title=title, data=datasets)

fmtFigSpace = """
<figure><canvas id="fig{!figno!}"></canvas>
<figcaption>Figure {!figno!}. {!caption!}</figcaption>
</figure>
"""

def makefigure(fig, caption=''):
    global fmtFigSpace
    return fmtFigSpace.format(figno=fig, caption=caption)

# use {!tag!} formatting - this makes it easier to edit javascript code
# in the templates. 
for s1, s2 in (('{', '{{'), ('{{!', '{'), ('}', '}}'), ('!}}', '}')):
    fmtReport, fmtChart, fmtFigSpace = (x.replace(s1, s2) for x in (fmtReport, fmtChart, fmtFigSpace))

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

def mkchart(fn, data, separate=False):
    
    # construct the HTML report using Chart.js for charts
    fs = [makefigure(1, "combined")]
    args0 = (("blue", "trimmed", [makelabel(x, 2) for x in data]),
             ("green", "merged", [makelabel(x, 3) for x in data]),
             #("brown", "ambiguous",[makelabel(x, 4) for x in datavec])
            )
    ch = [makemchart(1, args0, title="Combined")]
    if separate:
        for i, x in enumerate(args0, start=2):
            c, t, d = x
            fs.append(makefigure(i, t))
            ch.append(makemchart(i, [x], title=t[0].upper()+t[1:]))
        
    with open(fn, "wt") as dst:
        print(fmtReport.format(figspace='\n'.join(fs), charts='\n'.join(ch)), file=dst)

    return   

def main():
    global defdir, defrpt, defcsv
    parse = argparse.ArgumentParser(description="NGS Genotyping - Summary Report Program")
    parse.add_argument("-m", "--mergedir", default=defdir, help="data directory (contains *.log files) default="+defdir)
    parse.add_argument("-r", "--report", default=defrpt, help="name of summary report (HTML file) default="+defrpt)
    parse.add_argument("-c", "--csv", default=defcsv, help="name of output CSV file, default="+defcsv)
    parse.add_argument("-s", "--separate", action='store_true', help="include separate charts")
    # parse.add_argument("files", nargs="*", help="list of R1 files (default is clean/*_R1_*.fq.gz)")
    args = parse.parse_args()

    fns = glob.glob(os.path.join(args.mergedir, "*.log"))
    
    print(len(fns), "log files", file=sys.stderr)
    datavec = tuple(getvalues(fn) for fn in fns)
    
    # output the CSV file
    with open(args.csv, "wt") as dstfd:
        w = csv.writer(dstfd,  quoting=csv.QUOTE_MINIMAL, dialect='unix')
        w.writerow(['Well/Sample', 'Read Count', 'Cleaned', 'Joined', 'Ambiguous'])
        for r in datavec:
            w.writerow(r)
    mkchart(args.report, datavec, args.separate)
    return

if __name__=="__main__":
     main()