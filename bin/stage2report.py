#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created June 2020
@author: Bob Buckley and Cameron Jack
@version: 0.12
@version_comment:
@last_edit:
@edit_comment:

Draw 384-well plates using Chart.js

This runs after all the Echo stages and Miseq setup stages have completed
"""
defdir = "."
defrpt = "Stage2.htm"
defcsv = "Stage2.csv"

#import os
#import glob
#import sys
import csv
import argparse
import itertools
import collections

# Template for the HTML file
fmtReport = """
<!DOCTYPE html>
<html lang="en" dir="ltr" >
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type">
<title>NGS genotyping - MiSeq base report</title>
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/chart.js/dist/Chart.min.js"></script>
</head>
<body>
<h1>Combined DNA/Sample and PCR Plates Report</h1>
<p>Following are plate layouts for the PCR plates for a NGS Genotyping MiSeq run.</p>
{!figspace!}
<p>You can "mouse over" a well to see which mouse and assay is in a well, which source plate and combined-DNA plate well it came through.</p>
<p>The following table shows the detail for all the samples in the MiSeq NGS Genotyping run.
It shows the source plates, the various assays names, etc.
</p>
{!table!}
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
        maintainAspectRatio: true,
        scales: {
            xAxes: [{
                type: 'linear', offset: false, 
                ticks: { min: 0.5, max: 24.5, 
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
                    // beginAtZero: true,
                    min: 0.5, max: 16.5,
                    stepSize: 1,
                    callback: function(v, i, vs) { return 1<=v && v<=16?"ABCDEFGHIJKLMNOP".charAt(16-v):''; }, 
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
            pointRadius: 12, pointHoverRadius: 15,
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
<div class="chartjs">
<figure><canvas id="fig{!figno!}" class="chartjs"></canvas>
<figcaption>Figure {!figno!}. {!caption!}</figcaption>
</figure>
</div>
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


def output(filename, hdr, data, custom=False):
    """ output the HTML reports, custom==True for custom sample pipeline runs """
    dnaplates = sorted(set(x.dnaplate for x in data))
    dnadata = dict((dnabc, dict((w, list(g)) for w, g in 
            itertools.groupby(filter(lambda x:x.dnaplate==dnabc, data), key=lambda x: x.dnawell))) for dnabc in dnaplates)
    pcrplates = sorted(set(x.pcrplate for x in data))
    pcrdata = dict((bc, [ x for x in data if x.pcrplate==bc]) for bc in pcrplates)
    
    fs = []
    ch = []
    figno = 0
    
    def mklabdna(rs, custom=False):
        segs = ["Well: "+rs[0].dnawell+" from EP well: "+rs[0].EPwell]
        if custom:
            segs.append('sampleBarcode: '+rs[0].sampleBarcode)
        else:
            segs.append("mouseBarcode: "+rs[0].mouseBarcode)
            segs.append("strain: "+rs[0].strainName)
        segs.append("assay primers: "+', '.join(r.pcrplate+'-'+r.pcrwell+"="+r.primer for r in rs))
        return '; '.join(segs)
    
    for pn in dnaplates:
        # print("DNA/Sample Plate", pn, "has", len(dnadata[pn]), "wells.")
        figno += 1
        fs.append(makefigure(figno, "DNA/Sqample Plate "+pn))
        wdict = dnadata[pn]
        argx = [
                 ("lightblue", "Sample", [(None, ord('P')-ord(w[0]), int(w[1:])-1, mklabdna(wdict[w],custom)) for w in wdict.keys()]),
                 ("white", "empty", [(None, r, c, chr(ord('P')-r)+str(c+1)+'; empty') for c in range(24) for r in range(16) if chr(ord('P')-r)+str(c+1) not in wdict]),
               ] 
        ch.append(makeplate(figno, argx, title="DNA/Sample Plate "+pn))
        
    def mklab(r, custom=False):
        """ Make column labels """
        segs = ["Well: "+r.pcrwell]
        if custom:
            segs.append('sampleBarcode: '+r.sampleBarcode)
        else:
            segs.append('mouseBarcode: '+r.mouseBarcode)
            segs.append('strain: '+r.strainName)
            if r.mouseAssays!=r.primer:
                segs.append("Musterer assays: "+r.mouseAssays)            
        segs.append("assay: "+r.primer)
        segs.append("DNA plate, well: "+r.dnaplate+', '+r.dnawell)
        segs.append("EP plate, well: "+r.EPplate+", "+r.EPwell)
        return '; '.join(segs)
              
    for pn in pcrplates:
        # print("PCR Plate", pn, "has", len(pcrdata[pn]), "wells.")
        figno += 1
        fs.append(makefigure(figno, "PCR Plate "+pn))
        welz = set((ord('P')-ord(x.pcrwell[0]), int(x.pcrwell[1:])-1) for x in pcrdata[pn])
        argx = [
                 ("lightgreen", "Sample", [(None, ord('P')-ord(x.pcrwell[0]), int(x.pcrwell[1:])-1, mklab(x,custom)) for x in pcrdata[pn]]),
                 ("white", "empty", [(None, r, c, chr(ord('P')-r)+str(c+1)+'; empty') for c in range(24) for r in range(16) if (r,c) not in welz]),
                ]
        ch.append(makeplate(figno, argx, title="PCR Plate "+pn))
            
    tablefmt = """<table>
    {hdr}
    {rows}
    </table>"""
    def mkrow(r, element='td'):
        elfmt = "<{0}>{{}}</{0}>".format(element)
        return "<tr>"+''.join(elfmt.format(x) for x in r)+'</tr>'
    hdrrow = mkrow(hdr, element="th")
    tabdata = "\n".join(mkrow(r) for r in data)
    table = tablefmt.format(hdr=hdrrow, rows=tabdata)
        
    with open(filename, "wt") as dst:
        print(fmtReport.format(
                figspace = '\n'.join(fs), 
                charts = '\n'.join(ch),
                table = table),
              file=dst)
    return

def main():
    global defrpt, defcsv
    parse = argparse.ArgumentParser(description="NGS Genotyping - Summary Report Program")
    parse.add_argument("-r", "--report", default=defrpt, help="name of summary report (HTML file) default="+defrpt)
    parse.add_argument('--custom', action='store_true', help='Running custom pipeline rather than mouse pipeline')
    parse.add_argument("csv", default=defcsv, help="name of Stage2 CSV file")
    args = parse.parse_args()
    
    # read the CSV file
    with open(args.csv) as srcfd:
        src = csv.reader(srcfd, dialect="unix")
        hdr = next(src)
        Rec = collections.namedtuple('Rec', hdr)
        data = [Rec(*x) for x in src]
        
    output(args.report, hdr, data, custom=args.custom)
    return

if __name__=="__main__":
     main()

