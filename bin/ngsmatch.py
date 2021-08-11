#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Mar 2020
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.9
@version_comment: Match all possible target sequences. Drop per-well 'calling'.
@last_edit: 2021-07-22
@edit_comment: Inexact match was broken - minimising scoring rather than maximising

Read a NGS Genotyping Stage3.csv file
Read an NGS Genotyping reference file containing assay target sequences.
Process raw FASTQ files and check the sequences they contain against all known target sequences.
Raw MiSeq reads shold be in the 'raw' directory.
It cleans reads into files in the 'clean' directory, then merges them into an 'mclean'
directory ... along with a log of the cleaning and merging that tells us how many
reads there are at each stage.
This code uses BBtools to clean and merge paired-end read files.
"""

import os
import sys
import argparse
import re
import csv
import collections
from collections import Counter
import itertools
import subprocess
import gzip
import glob

import Bio.pairwise2 as pw2
import Bio.SeqIO

import echo

bclen = 8 # length of pseudo barcode
matchcutoff = 1.5 # should control via command line - must be >1.x

basetrans = str.maketrans("ACGT-", "TGCA-") # include '-' to allow for translating matched sequences

class ProgramFail(Exception):
    pass


def printv(*x, **kw):
    global args
    if args.verbose:
        print(*x, **kw)
    return


class Target(collections.namedtuple("TgtSeqTuple", "hdr seq")):
    """ Target sequence: name/hdr and R1 and R2 versions """
    def __new__(cls, idx, seq):
        """ construct from FASTA record - tuple of header and sequence """
        return super(Target, cls).__new__(cls, idx, seq)

    
    def offmatch(self, s, off):
        """ exact match with target - allow for an initial offset """
        ts = self.seq if not off else self.seq[off:]
        return ts.startswith(s) if len(s)<=len(ts) else s.startswith(ts)


def myopen(fn):
    """ handles gz files seamlessly """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt")
    return open(fn, errors="ignore")


def domatch(mcnt, s, n, cs):
    """
    exact matching for a (merged) read with all known target sequences

    Parameters
    ----------
    mcnt : TYPE
        DESCRIPTION. list counting number of matches to match candidates
    s : TYPE
        DESCRIPTION. sequence to match (merged read)
    n : TYPE
        DESCRIPTION. number of reads that have this (merged) sequence - if there 
        is a match add this to the appropriate value in mcnt
    cs : TYPE
        DESCRIPTION. Match candidates - list of sequences for the assay

    Returns
    -------
    True on success, else False

    """
    # it would be good to count the offsets below to see what they look like
    for j, cx in enumerate(cs):
        c =cx.seq
        # match with small offset - len(a) <= len(b)
        a, b = s, c
        if len(c)<len(s):
            a, b = c, s
            # note: 'in' operator does exact substring matching
        if len(b)<=len(a)+5 and a in b:
            mcnt[j] += n
            return True
    return False


def exact_match(match_cnt, s, n, targets):
    """
    exact matching for a (merged) read with all known target sequences.
    Length of seq to target must be within 5 bases.

    Parameters
    ----------
    match_cnt : Counter
        DESCRIPTION. list counting number of matches to match candidates
    s : sequence (string)
        DESCRIPTION. sequence to match (merged read)
    n : int
        DESCRIPTION. number of reads that have this (merged) sequence - if there 
        is a match add this to the appropriate value in mcnt
    targets : list(Target)
        DESCRIPTION. Match candidates - list of sequences for the assay

    Returns
    -------
    True on success, else False

    """
    # it would be good to count the offsets below to see what they look like
    for j, cx in enumerate(targets):
        c = cx.seq
        # match with small offset - len(a) <= len(b)
        a, b = s, c
        if len(c)<len(s):
            a, b = c, s
            # note: 'in' operator does exact substring matching
        if len(b)<=len(a)+5 and a in b:
            match_cnt[cx.hdr] += n
            return True
    return False


def best_match(match_cnt, s, n, targets):
    """
    Sequence name - inexact match of sequence to targets - returns the best matching target name
    s is a sequence, ts -s a list (set) or targets to base the name on
    Pairwise match - base name on matching name pluse difference
    
    The name is the most similar sequence with a description of the sequence 
    difference appended ... or just the sequence if there is no sufficiently
    similar sequence.

    Parameters
    ----------
    match_cnt : Counter
    s : string
        DNA sequence.
    n : int 
        counts of this sequence
    ts : list of Targets
        namedtuple with hdr and sequence.

    Returns
    -------
    True if an inexact within 75% was found, False if not.

    """
    def bestalign(naxs):
        score, best_align, bn = max((x[0].score, x[0], n) for n, x in naxs)
        return bn, best_align
    def mkstr(s1, s2, start):
        pos = -10
        for p, (ca, cb) in enumerate(zip(s1, s2), start=start):
            if ca!=cb:
                # only provide a number for first of a sequence
                yield ('' if pos+1==p else str(p))+ca+cb
                pos = p
        return
    varsep = " // "
    #print('before zs. Number of targets:', len(targets), 'Counts:', n, file=sys.stderr)
    #print(targets, file=sys.stderr)
    zs = [(t.hdr, pw2.align.localms(s, t.seq, 1, -1, -1, -1,one_alignment_only=True)) for t in targets if varsep not in t.hdr]
    #print(zs, file=sys.stderr)
    if zs:
        best_name, ax = bestalign(zs)
        #print(best_name, ax, file=sys.stderr)
        pos = ax.start
        sa, sb = [s[pos:ax.end] for s in (ax.seqA, ax.seqB)]
        if ax.score>len(sa)*0.75: # only use a name if we have a close-ish (75%) match
            dx = ''.join(mkstr(sa, sb, pos))
            match_cnt[varsep.join((best_name, dx))] += n
            #print('Got a hit for', n, 'counts', file=sys.stderr)
            return True
    return False


def main():
    """
    Read background data: target reference sequences file
    Then processes merged pairs files producing NGS report.
    
    Raises
    ------
    ProgramFail - fail silently...
        if the program fails for any (programmed) reason.
    """
    global args
    parse = argparse.ArgumentParser(description="NGS Reporting Program")
    parse.add_argument("-v", "--verbose", action="store_true", help="more reporting/output")
    #parse.add_argument("-d", "--datadir", default="mclean", help="data directory")
    parse.add_argument("-t", "--targets", default=None, help="file of targets in FASTA format")
    # parse.add_argument("-a", "--assays", default=None, help="file linking mice/samples to assay names")
    parse.add_argument('-o','--outfn', default='Results.csv', help='Path to output file name (CSV format)')
    parse.add_argument("stagefile", default="Stage3.csv", help="NGS genotyping Stage 3 file (default=Stage3.csv")
    args = parse.parse_args()
    
    # read the wells data
    fn = args.stagefile
    with open(fn) as srcfd:
        src = csv.reader(srcfd, dialect="unix")
        hdr = next(src)
        WRec = collections.namedtuple("WRec", hdr)
        wdata = sorted((WRec(*r) for r in src), key=lambda x:(x.pcrplate, x.pcrwell[0], int(x.pcrwell[1:]), x.primer))
        
    print(len(wdata), "sample wells to process.")
        
    ## get a set of assay family names - family name ends with first underscore char  
    #assays = frozenset(r.primer.split('_',1)[0] for r in wdata)
        
    if not args.targets:
        ts = [t for rpat in ("*ref*seq*.fa", "*ref*seq*.fasta", "*ref*seq*.txt") for t in glob.glob(os.path.join('..', 'library', rpat))]
        if not ts:
            print("no target reference file found.")
            raise ProgramFail
        if len(ts)>1:
            print("Multiple target reference files found:")
            for t in ts:
                print('   ', t)
            print("use the --targets option to avoid this message.")
        args.targets = ts[0]
        print("using targets reference:", args.targets)
    # read the target sequence file into a dictionary of target sequences
    #global fadict

    # make sure there are no lowercase chars in target sequences            
    with myopen(args.targets) as src:
        fadict = dict((x.id, Target(x.id, str(x.seq).upper())) for x in Bio.SeqIO.parse(src, "fasta"))
    # fagen = ((h, s.upper()) for h,s in fasta.fasta_reader(args.targets))
    # fadict = dict(sorted((x.hdr, x) for x in map(Target, fagen)))
    print("read", len(fadict), "target sequences.")
    # print("Seq IDs:", sorted(fadict.keys()))
    print()
    
    # associate sequences with the known assays
    #adict = {}
    #for a in sorted(assays):
    #    al = []
    #    for k in fadict.keys():
    #        if k.startswith(a):
    #            al.append(k)
    #    adict[a] = al
    
    # we may need more testing - are there enough sequences for the type of test?
    #noseq = [a for a in sorted(assays) if a not in adict]
    #if noseq:
    #    print('Assays with no associated sequences ...')
    #    for a in noseq:
    #        print(a)
    #    print()
            
    #sdict = dict((a, [fadict[k] for k in ks]) for a, ks in adict.items())
    # should check that sequences in each group are different
    
    if not os.path.isdir("raw"):
        print("raw directory is absent - please transfer MiSeq data.")
        raise ProgramFail
    
    for d in ("clean", "mclean"):
        if not os.path.isdir(d):
            os.mkdir(d)
                
    # process data for each sample well
    with open(args.outfn, "wt", buffering=1) as dstfd:
        dst = csv.writer(dstfd, dialect="unix")
        hdrres1 = ("readCount", "cleanCount", "mergeCount")
        hdrres2 = ("seqCount", "seqName")
        dst.writerow(x for xs in (hdr, hdrres1, hdrres2) for x in xs)
        # process samples wells - code handles bad (manually prepared) plate picklists with >1 assay in a well
        # wdata is already sorted.
        for k, gwrs in itertools.groupby(wdata, key=lambda x:(x.pcrplate, x.pcrwell)):
            wrs = list(gwrs)
            if len(wrs)>1:
                print('-'.join(k).ljust(len(k[0])+4), 'has', len(wrs), "assays:", ', '.join(w.primer for w in wrs))
            wr = wrs[0]
            # find source files
            print(wr.pcrplate, wr.pcrwell, file=sys.stderr)
            rfn= "*{}-{}_*_R1_*.fastq.gz".format(wr.pcrplate, echo.padwell(wr.pcrwell))
            print(rfn, file=sys.stderr)
            fn1s = glob.glob(os.path.join("raw", rfn))
            if not fn1s:
                print("no data for", '-'.join((wr.pcrplate, wr.pcrwell)))
                dst.writerow(wr)
                continue
            
            if len(fn1s)>1:
                print("too many reads files for", '-'.join((wr.pcrplate, wr.pcrwell)))
                print('   R1 files=', fn1s)
                # continue
            
            fnr1 = fn1s[0]
            fn1 = os.path.basename(fnr1) # just the one file
            fn2 = fn1.replace('_R1_', '_R2_')
            fnr2 = os.path.join('raw', fn2) 
            
            # find the data file
            if not os.path.isfile(fnr2):
                print("missing file:", fnr2)
                raise ProgramFail
                continue
            # use Windows file separators
            fncs = tuple(os.path.join("clean", fn) for fn in (fn1, fn2))
            bbmapd = os.path.join('..', 'bbmap', 'current')
            # unpick the file name a bit
            fnparts = fn1.split('_', 2)
            fnmfmt = "{}-{}_{}{{}}.{{}}".format(wr.pcrplate, echo.padwell(wr.pcrwell), fnparts[1])
            fnms = tuple(os.path.join("mclean", fnmfmt.format('_'+tx, "fastq.gz")) for tx in ('M', 'U1', 'U2'))
            fnlog = os.path.join("mclean", fnmfmt.format('', 'log'))
            if not all(os.path.isfile(fn) for fn in (fnms[0], fnlog)):
                # java and bbmap need to be properly installed
                # run BBDuk to clean reads - note: bash isn't available on Windows
                cmd = r"java -ea -Xmx1g -cp {} jgi.BBDuk in1={} in2={} out1={} out2={} t=2 qtrim=rl trimq=20 minlen=50 k=23 ktrim=r mink=11 hdist=1 overwrite=true ref=..\bbmap\resources\adapters.fa".format(bbmapd, fnr1, fnr2, fncs[0], fncs[1])
                pres1 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
            
                # run BBMerge to join paired-end reads
                # cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} verystrict=t".format(*((bbmapd,)+fncs+fnms)).split()
                cmd = r"java -ea -Xmx1g -cp {} jgi.BBMerge in1={} in2={} out={} outu1={} outu2={} overwrite=true pfilter=1".format(*((bbmapd,)+fncs+fnms)).split()
                pres2 = subprocess.run(cmd, check=True, text=True, timeout=30, capture_output=True)
                # could delete cleaned data once it's been merged.
                
                # keep the log output as record counts get used
                with open(fnlog, "wt") as log:
                    for s in (pres1.stdout, pres1.stderr):
                        if s:
                            print(s, file=log)
                    print('======', file=log)
                    for s in (pres2.stdout, pres2.stderr):
                        if s:
                            print(s, file=log)
                log1, log2 = ['\n'.join((p.stdout, p.stderr)) for p in (pres1, pres2)]

            else:
                print('Skipping cleaning and merging of reads for', fnmfmt.replace('{}.{}',''), file=sys.stderr)
                with open(fnlog) as log:
                    logdata = log.read()
                log1, log2 = logdata.split("======\n", 1)
            
            # get counts from BBDuk & BBmerge output
            readCount, cleanCount, cleanCount2, joinCount = -1, -1, -1, -1
            m = re.search(r'Input:\s+(\d+) reads', log1)
            if m:
                readCount = int(m.group(1))//2
            m = re.search(r'Result:\s+(\d+) reads', log1) 
            if m:
                cleanCount = int(m.group(1))//2
            m = re.search(r'Pairs:\s+(\d+)\s', log2)
            if m:
                cleanCount2 = int(m.group(1))
            m = re.search(r'Joined:\s+(\d+)\s', log2)
            if m:
                joinCount = int(m.group(1))
                
            if cleanCount!=cleanCount2:
                print("Pair count mismatch.", cleanCount, '!=', cleanCount2)
                raise ProgramFail
                
            fn = fnms[0] # merged reads file
            with myopen(fn) as src:
                # read a FASTQ file - sequence in in second line of each 4 line
                gs = (r.strip() for i, r in enumerate(src) if i%4==1)
                seqcnt = collections.Counter(gs)
            mergeCount = sum(seqcnt.values())
            assert mergeCount==joinCount # check data is "good"
            
            printv((wr.pcrplate+'-'+echo.padwell(wr.pcrwell)).ljust(6), str(mergeCount).rjust(6), wr.primer)
            
            # we are looking for sequences that match any target in fadict
            # what's the best matching algorithm?
            #assayfam = wr.primer.split('_',1)[0] # assay family name for this well/sample
            #assayfams = [r.primer.split('_',1)[0] for r in wrs] # multiple assay families in same well!!!
            targets = [fadict[k] for k in sorted(fadict)]
            match_cnt = Counter()
            low_cov_cutoff = max(mergeCount/20, 50)
            other_count = 0
            for seq, num in seqcnt.most_common():
                if num >= low_cov_cutoff:
                    if exact_match(match_cnt, seq, num, targets):
                        print('Exact match of', num, 'counts', file=sys.stderr)
                        continue
                    if best_match(match_cnt, seq, num, targets):
                        print('Inexact match of', num, 'counts', file=sys.stderr)
                        continue
                    else:
                        match_cnt[seq] += num
                        print('Missed matching', num, 'counts', file=sys.stderr)
                else:
                    match_cnt['other'] += num
                    other_count += 1
            # rename 'other' to include number of separate sequence variations
            match_cnt['other ('+str(other_count)+')'] = match_cnt['other']
            match_cnt['other'] = 0
            res1 = readCount, cleanCount, mergeCount
            resx = sorted(((match_cnt[key], key) for key in match_cnt if match_cnt[key] >= low_cov_cutoff), reverse=True)
            print(resx, file=sys.stderr)
            #resx = sorted(((cnt, target.hdr) for cnt, target in zip(mcnt, mseq) if cnt), reverse=True)
            for n, a in resx:
                printv('      ', str(n).rjust(6), a)
            res2 = tuple(x for p in resx for x in p)
            dst.writerow(wr+res1+res2)
            
            # should process extras and output to another file?  
    return
   
if __name__=="__main__":
    try:
        main()
    except ProgramFail:
        pass
