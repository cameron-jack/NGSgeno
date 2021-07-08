#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Jul 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy
@version: 0.8
@version_comment:
@last_edit:
@edit_comment:

Create SVG plate images, Summary and Result Table HTML files.
"""
templatefn = 'ResultPlate.templ' # name of template file!
rtfn = 'ResultTable.templ'

import os
import sys
import csv
import json
import collections

def rowdata(r):
    # this has to match the headers in the table template file
    rcol = [5,6,7] # right align these columns
    rx = [r.pcrplate+'&nbsp;'+r.pcrwell, r.EPplate+'&nbsp;'+r.EPwell,
          r.dnaplate+'&nbsp;'+r.dnawell, r.mouseBarcode, r.strainName]
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

def main():
    """
    create an SVG plate image file 
    based on a CSV list of well descriptors
    """
   
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
        wx[i][j]['mouseBarcode'] = r.mouseBarcode
        wx[i][j]['strain']  = r.strainName
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
    bclist = dict((k, str(i)) for i, k in enumerate(frozenset(r.mouseBarcode for r in rows), start=1))
    slist  = dict((k, str(i)) for i, k in enumerate(frozenset(r.strainName for r in rows), start=1))
    alist  = dict((k, str(i)) for i, k in enumerate(frozenset(r.primer for r in rows), start=1))
    
    def mkclasses(r):
        return ' '.join((r.pcrplate, 'bc'+bclist[r.mouseBarcode], 
                   's'+slist[r.strainName], 'a'+alist[r.primer]))
   
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
    chdata = [(' '.join((r.pcrplate, r.pcrwell, r.mouseBarcode, r.strainName, r.primer)), r.readCount, r.cleanCount, r.mergeCount) for r in rows if r.readCount]
    ls.mkchart("Summary.html", chdata)
    
    return

if __name__ == '__main__':
    main()   
