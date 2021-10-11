# -*- coding: utf-8 -*-
"""
@created: Nov 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.10
@version_comment:
@last_edit: 
@edit_comment: 

Program to combine Echo Survey files with coreresponding plate description files 
to produce a single file that describes a set of surveyed Echo source plates
for use in creating an Echo picklist.

The command line arguments are pairs of filenames:
    a description CSV file whose rows list a well and its contents
    an Echo survey file in Echo Survey format (including the well volume)
The final argument is the name of the destination file for the rows of combined data. 
"""

import sys
import csv
import itertools

def takeuntil(predicate, iterable):
    "like takewhile() but returns last value"
    for x in iterable:
        yield x;
        if predicate(x):
            break
    return

def widfix(wid):
    return wid[0]+wid[2] if wid[1]=='0' else wid
    
class EchoSurvey(dict):
    """
    dictionary subclass from Echo Survey file.
    Map well ID to volume in the well
    and retain some of the header information.
    """
    def __init__(self, filename):
        "Echo survey files have two formats"
        self.pid = None
        with open(filename, errors="ignore") as src:
            line1 = next(src)
            if line1.startswith('Date'):
                # format 1 - volume in r,c matrix
                # matrix format - 2 initial lines, a blank line then a matrix
                next(src), next(src) # skip two more lines
                csvrdr = csv.reader(src)
                hdr = next(csvrdr) # first row of matrix - cells contain columns
                assert hdr[0]==''
                assert all(c==int(idx) for c, idx in enumerate(hdr) if idx)
                for r in csvrdr:
                    assert r[0] # row has row ID (a letter)
                    for c, v in zip(hdr[1:], r[1:]):
                        if v:
                            self[r[0]+c] = v
                self.pid = ["Source[1]", "", "384PP_AQ_BP"]
                self.hdrrow = ["Source Plate Name",	"Source Plate Barcode",	"Source Plate Type"]
            else:                
                #format 2 - one row per well - well ID in r[3], volume in r[5] or r[6]
                assert line1.startswith('Run ID')
                list(takeuntil(lambda x:x.startswith('[DETAILS]'), src))
                # stop at blank line
                csvrdr = itertools.takewhile(lambda x:bool(len(x)), csv.reader(src))
                self.hdrrow = next(csvrdr)
                for r in csvrdr:
                    if not self.pid:
                        # self.pid = dict(zip('srcname srcbc srctype'.split(), r))
                        self.pid = r[:3]
                    if len(r)<6:
                        continue
                    if not r[:3]==self.pid:
                        print("bad description: all rows should have the same Plate ID info")
                        print("  pid  =", self.pid)
                        print("  r[:3]=", r[:3])
                    assert r[:3]==self.pid # all rows should have the same Plate ID info
                    # extract src well and volume - field 5 or field 6 (they seem to be the same)
                    self[widfix(r[3])] = r[5]
        return
    
def main():
    """
    Combine Echo Plate survey files with content list files
    to create an input file for the Echo picklist program.
    """
    if len(sys.argv)%2 or len(sys.argv)<4: # we need an odd number of filenames
        print("usage:", sys.argv[0], "[platedesc survey] ... destination", file=sys.stderr)
        print("Provided: ", sys.argv, file=sys.stderr)
        sys.exit(1)
    
        # open the destination file (last filename)
    with open(sys.argv[-1], "wt", newline='') as dstfd:
        dst = csv.writer(dstfd, dialect='unix', quoting=csv.QUOTE_MINIMAL)
        
        # step through filename pairs
        g = iter(sys.argv[1:-1])
        for fnum, dscfn, svyfn in zip(range(len(sys.argv)//2), g, g):
            dictsvy = EchoSurvey(svyfn)
            missing, cnt = [], 0
            with open(dscfn) as srcfd:
                src = csv.reader(srcfd)
                hdr = next(src)
                # drop trailing blank fields
                while hdr and not hdr[-1]:
                    del hdr[-1]
                rlen = len(hdr)
                assert bool(rlen)
                # construct the header row from the header of the first file
                if not fnum:
                    myhdr = dictsvy.hdrrow[:3]+hdr+["volume"]
                    dst.writerow(myhdr)
                else:
                    assert myhdr==dictsvy.hdrrow[:3]+hdr+["volume"]
                for r in src:
                    w = widfix(r[0])
                    if w in dictsvy:
                        dst.writerow(dictsvy.pid+[w]+r[1:rlen]+[dictsvy[w]])
                        cnt += 1
                    else:
                        missing.append(w)
                # probably should stop if any missing volumes
                # print("Echo survey file =", svyfn)
                # print("Plate description file =", dscfn)
                # print("matched wells:", cnt)
                if missing:
                     print(len(missing), "missing primer names for wells:", missing, file=sys.stderr)
                #     print("wells =", sorted(dictsvy.keys()))
                # print()
    return

if __name__ == '__main__':
    main()