#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Jul 2020
@author: Bob Buckley, Cameron Jack & Gabrielle Ryan, ANU Bioinformatics Consultancy

Create SVG plate images, Summary and Result Table HTML files for sequence counts info.
Swap r.assays/r.strainName on custom/mouse data
"""

import os
import sys         
import csv
import json
import collections
import functools
import jsonpickle

try:
    import bin.transaction as transaction
except ModuleNotFoundError:
    import transaction

try:
    import bin.util as util 
except ModuleNotFoundError:
    import util 
     
templatefn = os.path.join('library','ResultPlate.tpl') # name of template file
rtfn = os.path.join('library','ResultTable.tpl')


@functools.lru_cache
def generate_heatmap_html(jsonpickle_plate, pid, scaling=0.5):
    """ Given a given plate taken from transaction.get_plate(), generate a heatmap style plate image
    In future, add options for which things to highlight or scale, etc.
    Return a string containing html for the heatmap
    """
    plate = jsonpickle.decode(jsonpickle_plate, keys=True)
    purpose = 'Unknown purpose'
    if 'purpose' in plate:
        purpose = plate['purpose']

    chart_data_str = "chart.data = ["
    well_entries = []

    well_fields = []
    #default plate size: 384-wells
    plate_order = util.col_ordered_384
    width_num = 980
    height_num = 550
    left_padding = str(24)
    tool_text = "";
    if purpose == 'dna':
        #384-well plates
        well_fields = ['barcode', 'sex', 'strain', 'assays', 'gts']
        tool_text = '{y}{x}: ID: {barcode}\\nStrain: {strain}\\nAssays: {assays}\\nSex: {sex}\\nObserved GT: {gts}';
        
    elif purpose == 'sample':
        #sample plate: 96-wells
        plate_order = util.col_ordered_96
        well_fields = ['barcode', 'sex', 'strain', 'assays', 'gts']
        tool_text = '{y}{x}: ID: {barcode}\\nStrain: {strain}\\nAssays: {assays}\\nSex: {sex}\\nObserved GT: {gts}';
        width_num = 600
        height_num = 400
        left_padding = str(33)
    
    elif purpose == 'pcr':
        #384-well plates
        #purpose = 'PCR'
        well_fields = ['sampleBarcode', 'sex', 'strain', 'assays', 'primer']
        tool_text = '{y}{x}:\\nsample ID: {sampleBarcode}\\nStrain: {strain}\\nAssays: {assays}\\nPrimer: {primer}\\nSex: {sex}';

    elif purpose == 'taq_water':
        #taq/water plate: 6 wells
        plate_order = util.col_ordered_6
        well_fields = ['name', 'volume']
        tool_text = '{y}{x}: Type: {name}\\nVolume: {volume} μL';
        width_num = 400
        height_num = 300
        left_padding = str(37)

    elif purpose == 'primer':
        #384-well plates
        well_fields = ['primer', 'volume']
        tool_text = '{y}{x}: Primer: {primer}\\nVolume: {volume} μL';

    elif purpose == 'index':
        #384-well plates
        well_fields = ['idt_name', 'index', 'bc_name', 'oligo', 'volume']
        tool_text = '{y}{x}: Index Name: {idt_name}\\nIndex: {index}\\nBarcode Name: {bc_name}\\nVolume: {volume} μL';
    
    elif purpose == 'amplicon':
        #384-well plates
        well_fields = ['sampleNumber', 'sampleBarcode']
        tool_text = '{y}{x}:\\nSample Barcode: {sampleBarcode}\\nNumber: {sampleNumber}';

    for well in plate_order:
        entry_str = f'"y" : "{well[0].upper()}", "x" : "{well[1:]}"'
        if well not in plate:
            entry_str += ', "color" : colors.empty'
        elif purpose == 'index' and 'idt_name' not in plate[well]:
            entry_str += ', "color" : colors.empty'
        elif purpose == 'primer' and 'primer' not in plate[well]:
            entry_str += ', "color" : colors.empty'
        else:
            entry_str += ', "color" : colors.passed'
            sample_data = plate[well]
            for field in well_fields:
                if field in sample_data:
                    if field in ['gts']:
                        field_data = '; '.join(sample_data[field])
                    elif field in ['volume']:
                        field_data = sample_data[field]/1000
                    else:
                        field_data = sample_data[field]
                    if field_data:
                        entry_str += f', "{field}" : "{str(field_data)}"'
        well_entries.append(entry_str)
    chart_data_str += '\n{' + '},\n{'.join(well_entries) + '}\n' + '];\n'

    width = str(int(width_num*scaling)) + 'px'
    height = str(int(height_num*scaling)) + 'px'
    font_axis = str(int(20*scaling))
    font_label = str(int(12*scaling))
    font_popup = str(int(14*scaling))
    header_str = """
    <!-- Styles -->
    <style>
    #chartdiv {
      width: WWWWWW;
      height: HHHHHH;
      margin: 0 auto;
      display: flex;
      align-items: center;
      background-color: white;
      background-size: 95%;
      border: 5px;
      padding: 1px;
      border-style: solid;
      border-color: #bfdbf2;
      border-radius: 25px;
      font-family: helvetica;
    }
    #plate_container {
        height: auto !important
    }
    h1 {
      font-family: helvetica;
    }
    h2 {
      font-family: helvetica;
    }
    p {
      font-family: helvetica;
    }

    .lbox{
      display:inline-block;
      height: 90vh;
       width: 59%;
        background-color: white;
        overflow: hidden;
    }

    .rbox{
      display:inline-block;
      height: 80vh;
      width: 39%;
      transform: translate(0%, -12%);
      background-color: lightgrey;
      overflow: scroll;
      font-size: FPFPFP;
      font-family: helvetica;
    }


    </style>

    <!-- Resources -->
    <script src="https://www.amcharts.com/lib/4/core.js"></script>
    <script src="https://www.amcharts.com/lib/4/charts.js"></script>
    <script src="https://www.amcharts.com/lib/4/themes/animated.js"></script>

    <!-- Chart code -->
    <script>
    am4core.ready(function() {

    // Themes begin
    am4core.useTheme(am4themes_animated);
    // Themes end

    var chart = am4core.create("chartdiv", am4charts.XYChart);
    chart.hiddenState.properties.opacity = 0; // this creates initial fade-in

    chart.maskBullets = false;

    var xAxis = chart.xAxes.push(new am4charts.CategoryAxis());
    var yAxis = chart.yAxes.push(new am4charts.CategoryAxis());

    xAxis.dataFields.category = "x";
    yAxis.dataFields.category = "y";

    xAxis.renderer.grid.template.disabled = true;
    xAxis.renderer.minGridDistance = 50;
    xAxis.renderer.fontSize=FAFAFA;
    xAxis.renderer.opposite=true;


    yAxis.renderer.grid.template.disabled = true;
    yAxis.renderer.inversed = true;
    yAxis.renderer.minGridDistance = 50;
    yAxis.renderer.fontSize=FAFAFA;

    var series = chart.series.push(new am4charts.ColumnSeries());
    series.dataFields.categoryX = "x";
    series.dataFields.categoryY = "y";
    series.dataFields.value = "value";
    series.sequencedInterpolation = true;
    series.defaultState.transitionDuration = 1000;

    // Set up column appearance
    var column = series.columns.template;
    column.strokeWidth = 3;
    column.strokeOpacity = 1;
    column.stroke = am4core.color("lightgrey");
    //column.tooltipText = "{y}{x}: {value.workingValue.formatNumber('#.')}";
    column.tooltipText = "TTTTTT";
    column.width = am4core.percent(85);
    column.height = am4core.percent(85);
    column.column.cornerRadius(60, 60, 60, 60);
    column.propertyFields.fill = "color";

    // Set up bullet appearance
    var bullet1 = series.bullets.push(new am4charts.CircleBullet());
    bullet1.circle.propertyFields.radius = "value";
    bullet1.circle.fill = am4core.color("#FFF");
    bullet1.circle.strokeWidth = 0;
    bullet1.circle.fillOpacity = 0.4;
    bullet1.interactionsEnabled = false;

    var bullet2 = series.bullets.push(new am4charts.LabelBullet());
    bullet2.label.text = "{y}{x}";
    bullet2.label.fill = am4core.color("#000");
    bullet2.zIndex = 1;
    bullet2.fontSize = FLFLFL;
    bullet2.interactionsEnabled = false;

    // define colors
    var colors = {
        "empty": "lightgrey",
        "passed": "#94f043",
        "problem": "#ffff33",
        "failed": "#ff1111",
        "critical": "#ca0101",
        "bad": "#e17a2d",
        "medium": "#e1d92d",
        "good": "#5dbe24",
        "verygood": "#0b7d03"
    };
                                       +
    // Set the container height dynamically
    function setContainerHeight() {
        var container = document.getElementById("plate_container");
        var chartdiv = document.getElementById("chartdiv");
        var containerHeight = chartdiv.scrollHeight; // Get the height of chartdiv element
        var windowHeight = window.innerHeight;

        // Adjust the container height as needed
        var desiredContainerHeight = containerHeight + 100; // Add an offset if necessary
        var newContainerHeight = Math.min(desiredContainerHeight, windowHeight * 0.8); // Limit the container height to a maximum of 80% of window height
        container.style.height = newContainerHeight + "px";
    }

    //setContainerHeight(); // Initial height calculation

    // Recalculate container height on window resize
    window.addEventListener("resize", function() {
        setContainerHeight();
    });

    """.replace('WWWWWW',width).replace('HHHHHH',height).replace('FAFAFA',font_axis).\
                replace('FLFLFL', font_label).replace('FPFPFP',font_popup).replace('LLLLLL', left_padding).\
                    replace('TTTTTT', tool_text)

    footer_str = """
    var baseWidth = Math.min(chart.plotContainer.maxWidth, chart.plotContainer.maxHeight);
    var maxRadius = baseWidth / Math.sqrt(chart.data.length) / 2 - 2; // 2 is jast a margin
    series.heatRules.push({ min: 10, max: maxRadius, property: "radius", target: bullet1.circle });

    chart.plotContainer.events.on("maxsizechanged", function() {
        var side = Math.min(chart.plotContainer.maxWidth, chart.plotContainer.maxHeight);
        bullet1.circle.clones.each(function(clone) {
            clone.scale = side / baseWidth;
        })
    })

    }); // end am4core.ready()
    </script>

    <!-- HTML -->
    <body>
        <!-- <h1 style="text-align:center;color:#458cde">PUPUPU Plate PPPPPP</h1> Save the height-->
        <div id="plate_container">
            <div id="chartdiv"></div>
        </div>
    </body>
    """.replace('PUPUPU', purpose.capitalize()).replace('PPPPPP', str(util.unguard_pbc(pid, silent=True)))

    return header_str + chart_data_str + footer_str

def rowdata(r):
    # this has to match the headers in the table template file
    try:
        rcol = [5,6,7] # right align these columns
        rx = [r.pcrplate+'&nbsp;'+r.pcrwell, r.EPplate+'&nbsp;'+r.EPwell,
              r.dnaplate+'&nbsp;'+r.dnawell, r.sampleBarcode, r.assays]
        #rx.append(r.assay if r.assay==r.MustAssay else r.assay+'<br/>('+r.MustAssay+')')
        rx.append(r.primer)
        # rx += [r.readCount, r.cleanCount, r.mergeCount]
        rx.append("m: "+r.mergeCount+"<br/>c: "+r.cleanCount+"<br/>r: "+r.readCount if r.readCount else None)
        rz = r.matches
        rx.append('<br/>'.join(rz[::2]))
        rx.append('<br/>'.join("%.1f"%(float(x)*100.0/int(r.mergeCount)) for x in rz[::2]))
        rx.append('<br/>'.join(rz[1::2]))
        # trim missing fields from the end of the data
        while rx and not rx[-1]:
            rx.pop()
        return ''.join(('<td class="r">' if c in rcol else '<td>')+d+'</td>' for c, d in enumerate(rx))
    except Exception as exc:
        print(exc, msg='Error in makehtml.rowdata')


def main():
    """
    create an SVG plate image file 
    based on a CSV list of well descriptors
    """
    try:
        assert sys.argv[1].endswith('.csv')
    
        # read the Results CSV file into a table
        with open(sys.argv[1]) as srcfd:
            src = csv.reader(srcfd)
            hdr = next(src)
            rowlen = len(hdr)-2 # last two are first of repeated pairs
            Rec = collections.namedtuple('Rec', hdr[:rowlen]+['matches'])
            # print('fields =', Rec._fields)
            # print(len(Rec._fields), 'fields')
            def mkrow(rfs):
                p1, p2 = rfs[:rowlen], [rfs[rowlen:]]
                if len(p1)<rowlen:
                    p1 += [None]*(rowlen-len(p1))
                return Rec(*(p1+p2))
            data = [mkrow(x) for x in src]
    
        # create a dictionary with 16x24=384 entries for each PCR plate, and an empty array for the data
        pd = dict((pk, [[None for j in range(16)] for i in range(24)]) for pk in frozenset(r.pcrplate for r in data))
        for r in data:
            wx, pos = pd[r.pcrplate], r.pcrwell
            j, i = ord(pos[0])-ord('A'), int(pos[1:])-1
            wx[i][j] = {'pos':pos, 'reads': r.readCount, 'clean': r.cleanCount, 'merge': r.mergeCount}
            wx[i][j]['source']  = r.EPplate+'-'+r.EPwell+' => '+r.dnaplate+'-'+r.dnawell
            wx[i][j]['sampleBarcode'] = r.sampleBarcode
            wx[i][j]['strain']  = r.assays
            wx[i][j]['primer']   = r.primer
            wx[i][j]['results'] = r.matches        
    
        # output the first results file based on its template file
        global templatefn
        tfn = os.path.join('..', 'library', templatefn)
        with open(tfn) as src:
            htmlfmt = src.read()
        # template uses {!field!} instead of python format fields - easier to edit.
        htmlfmt = htmlfmt.replace('{', '{{').replace('}', '}}').replace('{{!', '{').replace('!}}', '}')
    
        x = "[\n    "+',\n    '.join(json.dumps({'pid':pk, 'wells':pd[pk]}) for pk in sorted(pd.keys()))+'  ]'
        html = htmlfmt.format(plates=x) # just one field at present
    
        dfn = sys.argv[1][:-4]+".html"
        with open(dfn, "wt") as dst:
            dst.write(html)
        
        # outputs the second results file (Table file) using its template file
        rows = data
        bclist = dict((k, str(i)) for i, k in enumerate(frozenset(r.sampleBarcode for r in rows), start=1))
        slist  = dict((k, str(i)) for i, k in enumerate(frozenset(r.assays for r in rows), start=1))
        alist  = dict((k, str(i)) for i, k in enumerate(frozenset(r.primer for r in rows), start=1))
    
        def mkclasses(r):
            return ' '.join((r.pcrplate, 'bc'+bclist[r.sampleBarcode], 
                       's'+slist[r.assays], 'a'+alist[r.primer]))
   
        tabrows = "".join("  <tr class='{classes}'>\n    ".format(classes=mkclasses(r))+rowdata(r)+"\n  </tr>\n" for r in rows)
    
        global rtfn
        tfn = os.path.join('..', 'library', rtfn)
        with open(tfn) as src:
            htmlfmt = src.read()
        # template uses {!field!} instead of python format fields - easier to edit.
        htmlfmt = htmlfmt.replace('{', '{{').replace('}', '}}').replace('{{!', '{').replace('!}}', '}')
    
        sxstr = '\n'.join("    <option value='{}'>{}</option>".format(n, s) for n, s in sorted(enumerate(slist, start=1), key=lambda x:x[1]))
        axstr = '\n'.join("    <option value='{}'>{}</option>".format(n, s) for n, s in sorted(enumerate(alist, start=1), key=lambda x:x[1]))
        bcstr = '\n'.join("    <option value='{}'>{}</option>".format(n, s) for n, s in sorted(enumerate(bclist, start=1), key=lambda x:x[1]))
        html = htmlfmt.format(tabrows=tabrows, sx=sxstr, ax=axstr, bcx=bcstr)
    
        dfn = sys.argv[1][:-4]+"-Tab.html"
        with open(dfn, "wt") as dst:
            dst.write(html)
        
        # create the Summary.html file - code in another module
        import logsummary as ls
        #print('row =', rows[0])
        chdata = [(' '.join((r.pcrplate, r.pcrwell, r.sampleBarcode, r.assays, r.primer)), r.readCount, r.cleanCount, r.mergeCount) for r in rows if r.readCount]
        ls.mkchart("Summary.html", chdata)
    
        return
    except Exception as exc:
        print(exc, msg='Error in makehtml.main')


if __name__ == '__main__':
    main()   