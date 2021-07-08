# -*- coding: utf-8 -*-
"""
@created Aug 2020
@author: Bob Buckley and Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.8
@version_comment: supports custom and mouse pipelines
@last_edit: 2021-07-08
@edit_comment: removed all local library paths

This code checks assay data (from the "library") for the NGS Geno application
The challenges are:
    1. Assays names from the Musterer data need to map to assay names in the assay primer plate description.
    2. Assay primer plate names need to match sequence names in the reference file.
"""

import os
import glob
import csv
import itertools
import sys
import argparse

def grouper(xs, kf=lambda x:x[0]):
    "group pairs: a, b, c ... => (kf(a),[a,b, ...]), (kf(c), [c, ...]), ..."
    return itertools.groupby(sorted(xs, key=kf), key=kf)

class PrimerLookup:
    """
    Class to represent the mappings and relationships between "observable names"
    in Musterer and the assay names used in primer plate descriptions ... and the
    names of the sequences in the reference file.
    For custom runs, the observable/assay names will be the same in most instances, but by keeping
    the dual label structure we have the flexibility to have different names for these.
    
    holds useful dictionaries:
        self.mmdict - Musterer mapping dictionary - maps Musterer observable names to primer names
        self.pfdict - maps primer family names to a list of primers and their wells
        self.rfdict - dictionary of references
        
    self.isna is a set of N/A Musterer observables
    """
    
    def __init__(self, fnmm=None, fnpp=None, fnref=None):
        """
        Read the Musterer mapping file - map Musterer names to sets of assay family names 
        fnmm - assay_list_<date>.csv file, maps Musterer names to assay/sequence names (custom inputs should
        have duplicated columns)
        fnpp - file name primer layout plate 'primer_layout_<date>.csv'
        fnref - file name reference sequences 'reference_sequences_<date?>.txt
        """
        
        fns = sorted(glob.glob(fnmm), reverse=True)
        self.fnmm = fns[0]
        with open(fns[0]) as srcfd:
            src = csv.reader(srcfd)
            hdr = next(src)
            data = sorted([a.strip() for a in x] for x in src if len(x)==2)
        # what to do if first column is empty?
        self.mmdict = dict((k, frozenset(x[1].split('_',1)[0] for x in g if x[1]!='N/A')) for k, g in grouper(data))
        #self.mmdict = dict((m, f.split('_',1)[0]) for m, f in data if m and f!='N/A')
        self.isna = frozenset(x[0] for x in data if x[1]=='N/A')
        del data

        # read the set of known primers, create map from primer to the PP wells where it can be found
        # known assays are names in the "primer names pooled" column of the primer plate
        # the following should allow multiple files - but how that will be done has not been decided
        # look for local primer descriptions first
        fns = sorted(glob.glob(fnpp), reverse=True)
        #if not fns: # no local file - use library
        #    fns = sorted(glob.glob(os.path.join('..', 'library', fnpp)), reverse=True)
        self.fnpp = fns[0]
        with open(fns[0]) as srcfd:
            src = csv.reader(srcfd)
            hdr = next(src)
            # the following is probably (1,0)
            anidx, widx = hdr.index("primer names pooled"), hdr.index("plate position on Echo 384 PP")
            def echofix(w):
                # fix wells for Echo format A01 => A1
                return w[0]+w[2] if w[1]=='0' else w
            pws = sorted((x[anidx].strip(), echofix(x[widx].strip())) for x in src)
        # dictionary of primers => list of wells
        self.ppdict = dict((k, list(v[1] for v in vs)) for k, vs in grouper(pws))
        # primer families derived from primer names (keys from ppdict)
        self.pfdict = dict((k, list(g)) for k, g in grouper(sorted(self.ppdict.keys()), kf=lambda x:x.split('_',1)[0]))
        del pws
        
        # read names of the reference sequences from the reference file
        # build a dictionary indexed by primer names from reference names
        if fnref:
            fns = sorted(glob.glob(fnref), reverse=True)
            self.fnref = fns[0]
            with open(fns[0]) as src:
                hdr = next(src)
                refgen = (x[1:].strip() for x in src if x.startswith('>'))
                # gen = itertools.groupby(sorted(refgen), key=lambda s:keyfix(s.split('_',1)[0]))
                gen = itertools.groupby(sorted(refgen), key=lambda s:s.split('_',1)[0])
                self.rfdict = dict((k, list(g)) for k, g in gen)
        
        return
    
    def checkmust(self):
        "check consistency of assay_list (Musterer) names with primer plate names"
        mustunk = frozenset(x for xs in self.mmdict.values() for x in xs if x not in self.pfdict)
        return mustunk
    
    def checkrefs(self):
        "check consistency of reference names to primer names"
        refnames = frozenset(self.rfdict.keys())
        primnames = frozenset(self.pfdict.keys())
        return len(refnames&primnames), refnames-primnames, primnames-refnames
        
def main():
    "check library files for NGS Geno pipeline"
    parser = argparse.ArgumentParser()
    parser.add_argument('fnref', help='Path to reference file (FASTA format)')
    parser.add_argument('fnpp', help='Path to primer plate layout file')
    parser.add_argument('fnmm', help='Path to Musterer to primer family names')
    args = parser.parse_args()
    
    pl = PrimerLookup(fnref=args.fnref,fnpp=args.fnpp,fnmm=args.fnmm)
        
    print(len(pl.ppdict), 'primers in the primer plate:', pl.fnpp)
    print(len(pl.pfdict), 'primer families.')
    mpfs = [f for f, ps in pl.pfdict.items() if len(ps)>1]
    if mpfs:
        print('   ', len(mpfs), "multi-primer families:", ' '.join(mpfs))
    print()
    print(sum(len(x) for x in pl.rfdict.values()), 'sequences in file:', pl.fnref)
    print(len(pl.rfdict), 'family names found in sequence names.')
    print()
    
    mxs = pl.checkmust()
    print(len(pl.mmdict), "mappings from Musterer to primer family names:", pl.fnmm)
    if pl.isna:
        print("N/A Musterer assay names:", ' '.join(pl.isna))
    if mxs:
        print(len(mxs), "Musterer map targets with no matching primer target family in primer plate:")
        print(" ", "\n  ".join(sorted(mxs)))
    print()
   
    n, pxr, pxg = pl.checkrefs()
    if n:
        print(n, 'sequence family names match primer family names (from primer plates).')
        if not pxr and not pxg: print()
    if pxr:
        print(len(pxr), "sequence family names are not names in a primer plate:")

        pnlen = max(len(pn) for pn in pxr)+1
        for pn in sorted(pxr):
            rns = pl.rfdict[pn]
            print(' ', pn.ljust(pnlen), rns[0])
            for x in rns[1:]:
                print(' ', ' '.ljust(pnlen), x)
        print()
    if pxg:
        xstr = "primer name from the primer plate does not match any sequence family name:" \
            if len(pxg)==1 else "primer names from the primer plate do not match sequence family names:"
        print(len(pxg), xstr)
        print(" ", "\n  ".join(sorted(pxg)))
        print()
    
    return    
        
if __name__ == '__main__':
    main()
