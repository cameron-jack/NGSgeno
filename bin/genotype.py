#!/usr/bin/env python3

import csv
import collections
import sys
from argparse import ArgumentParser as AP
from itertools import chain

"""
    Take read counts for observed alleles in the report so far and create a new report with called genotypes.
    Types of assays, and thus the possible genotypes, are defined in the groups_XXX.csv file
    these may be simplex (wt/wt, mut/wt, mut/mut), complex (present/absent), or new (previously unobserved)
"""

class Assay(object):
    """
        Refers to an entry from groups_xxx.csv
    """
    def __init__(self, family, primer, reference, assay_type):
        self.family = family
        self.primer = primer  # primer name
        self.references = {reference}  # set of assay names, also use primer
        self.assay_type = assay_type.lower()  # simplex, complex, or new
        self.references.add(primer)


    def __str__(self):
        outstr = ','.join([self.family, self.primer, ';'.join(self.references), self.assay_type])
        return outstr


    def add_reference(self, reference):
        self.references.add(reference)


def read_assay_groups(fn):
    """
        Read the groups.csv file to for determining assays and assay type
    """
    groups = dict()
    with open(fn, 'rt', encoding='windows-1252') as f:
        csvr = csv.reader(f, delimiter=',')
        previous_family = ''
        for file_row in csvr:
            row = [c.strip() for c in file_row]  # remove emtpy spaces and such
            assert len(row) == 4, '4 fields per row required in groups.csv'
            #family, primer, reference, assay_type = row
            family = row[0]
            if family is None or family == '':
                if previous_family is None or previous_family == '':
                    print('cannot have empty family entry when no previous family '+\
                            'has been seen', file=sys.stderr)
                    sys.exit(1)
                reference = row[2]
                if reference is None or reference == '':
                    continue  # dud entry?
                assay = groups[previous_family]
                assay.add_reference(reference)
            else:
                previous_family = family
                # New assay
                if family in groups:
                    reference = row[2]
                    if reference is None or reference == '':
                        continue  # dud entry?
                    assay.add_reference(reference)
                else:
                    assay = Assay(*row)
                    groups[assay.family] = assay
    return groups


def genotype(assays, matches, minreads=100, assay_type='new', must_assay=None, family=None):
    """
        Report observed genotype from counts data
        inputs: dictionary of Assays, list of tuples (counts, allele), min number of reads,
            assay_type is {'simplex', 'complex', 'new'},
            must_assay is the name of the assay that Musterer expects,
            family is the assay family that must_assay belongs to
        outputs: genotype (str), list of allele ratios for each retained allle,
            list of tuples(count,allele)
        For simplex types, return wt/wt, wt/mut, mut/wt, or mut/mut
        For complex types, return the alleles present
        For new types, return 'NEW'
    """
    if must_assay is None:
        return '', [], []

    # group together counts and their assay reference name
    if not matches or len(matches) == 0:
        return '', [], []
    #counts_alleles = list(zip(map(int, matches[::2]), matches[1::2]))
    counts_alleles = [(int(c),a) for c,a in matches]
    max_allele_count = max([c for c,a in counts_alleles])
    # only keep assays in meaningful proportions
    counts_alleles = sorted([(c,a) for c,a in counts_alleles if c > max_allele_count/10], reverse=True)
    #print (counts_alleles, file=sys.stderr)

    # fail if we have less than some sensible minimum reads across all alleles
    all_counts = [c for c,a in counts_alleles]
    if sum(all_counts) < minreads:
        return '', [], []

    total_count = sum(all_counts)
    allele_ratios = [c/total_count for c in all_counts]

    # fail on new assays
    if assay_type == 'new':
        return 'NEW', allele_ratios, counts_alleles

    # decide genotype based on observed alleles
    family_counts_alleles = [(c,a) for c,a in counts_alleles if a.split('_')[0] == family]
    num_family_alleles = len(family_counts_alleles)
    if num_family_alleles == 0:
        return 'Assay mismatch! Expected: ' + must_assay, [], counts_alleles
    elif num_family_alleles == 1:  # hom - should be wt/wt or mut/mut
        cnt1, al1 = family_counts_alleles[0]
        al1 = al1.lower()
        if assay_type == 'simplex':
            if 'wt' in al1.lower():
                return 'wt/wt', [1.0], counts_alleles
            else:
                return 'mut/mut', [1.0], counts_alleles
        else:  # complex
            if '_pos_' in al1:
                return 'pos', [1.0], counts_alleles
            else:
                return 'neg', [1.0], counts_alleles
    elif num_family_alleles == 2:  # hom or het
        cnt1, al1 = family_counts_alleles[0]
        cnt2, al2 = family_counts_alleles[1]
        on_target_cnt = cnt1 + cnt2
        AR1 = cnt1/on_target_cnt
        al1 = al1.lower()
        al2 = al2.lower()
        # if less than 15% of contributing reads, it's a hom
        if cnt2 < 0.15*(cnt1 + cnt2):
            if assay_type == 'simplex':
                if 'wt' in al1:
                    return 'wt/wt', [AR1], counts_alleles
                else:
                    return 'mut/mut', [AR1], counts_alleles
            else:  # complex
                if '_neg_' in al1:
                    return 'neg', [AR1], counts_alleles
                else:
                    return 'pos', [AR1], counts_alleles
        # if allele2 has more than 35% of contributing reads, it's a het
        elif cnt2 > 0.35*(on_target_cnt):
            if assay_type == 'simplex':
                if 'wt' in al1:
                    return 'wt/mut', [AR1, 1-AR1], counts_alleles
                else:
                    return 'mut/wt', [AR1, 1-AR1], counts_alleles
            else:  # complex - should NOT occur as each should be in their own well
                if '_neg_' in al1 and '_neg_' in al2:
                    return al1 + ': neg +' + al2 + ': neg', [AR1, 1-AR1], counts_alleles
                elif '_neg_' in al1:
                    return al1 + ': neg +' + al2 + ': pos', [AR1, 1-AR1], counts_alleles
                elif '_neg_' in al2:
                    return al1 + ': pos +' + al2 + ': neg', [AR1, 1-AR1], counts_alleles
                else:
                    return al1 + ': pos +' + al2 + ': pos', [AR1, 1-AR1], counts_alleles
        # ambiguous in between
        else:
            return 'ambig', [AR1, 1-AR1], counts_alleles
    # more than 2 alleles
    else:
        return 'multi', allele_ratios, counts_alleles


def main(result_fn, assay_group_fn, out_fn=None):
    """
        For each row in the results table, extract alleles and counts, then call the genotypes
    """
    assay_groups = read_assay_groups(assay_group_fn)
    print('Read', len(assay_groups), 'assay families from ', assay_group_fn, file=sys.stderr)
    # print('\n'.join([':'.join([k,v.__str__()]) for k,v in assay_groups.items()]), file=sys.stderr)

    with open('Results.csv', 'rt') as srcfd:
        src = csv.reader(srcfd)
        hdr = next(src)
        # build new header for output now
        hdrstr = ','.join(hdr[:23] + ['assay_type','genotype', 'alleleRatios', 'efficiency'] + hdr[23:])
        rowlen = len(hdr)-2 # last two are first of repeated pairs
        Rec = collections.namedtuple('Rec', hdr[:rowlen]+['matches'])
        def mkrow(r):
            g = (x for x in r[rowlen:])
            matches = list(zip(g,g))
            p1, p2 = r[:rowlen]+['']*max(0, rowlen-len(r)), [matches]
            return Rec(*(p1+p2))
        records = list(map(mkrow, src))

    # set up output
    if out_fn:
        out = open(out_fn, 'wt')
    else:
        out = sys.stdout
    print(hdrstr, file=out)

    for i,r in enumerate(records):
        #print(r.MustAssay, r.assay,r.readCount, r.cleanCount, r.mergeCount, r.matches, file=sys.stderr)
#        if i<5:
#            print(r.pcrplate, r.pcrwell, r.assay, r.matches, file=sys.stderr)
        family = r.mouseAssays.split('_')[0]
        if r.matches is not None and r.matches != []:
            assay_group = assay_groups.get(family, None)
            if not assay_group:
                gt, allele_ratios, counts_alleles = genotype(assay_groups, r.matches)
            else:
                gt, allele_ratios, counts_alleles = genotype(assay_groups, r.matches,
                        assay_type=assay_group.assay_type, must_assay=r.MustAssay, family=family)
            eff = "{:.2f}".format(sum([c for c,a in counts_alleles])/int(r.readCount))
            allele_ratios = '/'.join(["{:.2f}".format(a) for a in allele_ratios])
            #print (gt, allele_ratios, counts_alleles, file=sys.stderr)
        else:
            gt = ''
            allele_ratios = ''
            counts_alleles = ''
            eff = ''
#        if i<5:
#            print(gt, allele_ratios, counts_alleles, eff, file=sys.stderr)
        assay_type = ''
        if assay_group and assay_group.assay_type:
            assay_type = assay_group.assay_type
        original_fields = '","'.join(map(str, list(r[:23])))
        new_fields = '","'.join(map(str, [assay_type, gt, allele_ratios, eff]))
        later_fields = '","'.join(map(str, r[23:25]))
        match_fields = '","'.join(map(str, list(chain(*r[26]))))
        outstr = '"' + original_fields + '","' + new_fields + '","' + later_fields + '","' + match_fields + '"'
        outstr = outstr.replace('None', '').replace('none','')
        print(outstr, file=out)
    out.close()
    print('\nGenotyping complete.', file=sys.stderr)


if __name__ == '__main__':
    parser = AP()
    parser.add_argument('resultfile', help='Path to Results.csv file')
    parser.add_argument('-g', '--groups', required=True, help='Path to groups_XXX.csv file describing assay groups')
    parser.add_argument('-o', '--outfile', help='Path to output file. If not present will output to stdout')
    args = parser.parse_args()
    main(args.resultfile, args.groups, out_fn=args.outfile)
