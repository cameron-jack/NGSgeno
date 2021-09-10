#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: Dec 2020
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.10
@version_comment:
@last_edit:
@edit_comment:

Take read counts for observed alleles in the report so far and create a new report with
called genotypes.

Rather than treating wells as independent, this new code collects records for each mouse.
This allows for complex assays - with a mutant and wt presence/absence - observations to be
combined and sanity checked in a meaningful way. We will reconstruct the final table to resemble
the allele counting table.

New NGS-specific assays are coming which are a single observable per assay, rather than an
assay family. These also have only mut/wt/pos/neg as reportable options. We still need to
cope with old mappings though, so there is an Excel file for this. We also need to convert from
parental assays to NGS assays done for the first time when it comes to 'sanity checking'.

NGS genotyping sanity checking:
    Genotyping checklist is set up with basic Boolean operators: AND, NOT & OR
In some cases the breeding configuration will be a pair or trio, this will need to be taken into
account. Breeding trios are set up as 2 females to 1 male ratio, conversely breeding pairs are a
1:1 ratio.
Keywords:
Wildtype = WT OR wt/wt OR wt/Y
Heterozygous = HET OR wt/mut OR mut/wt
Mutant = MUT OR mut/mut or mut/Y
Positive = POS is the presence of transgene
Negative= NEG is the absence of transgene

1.	Assays for mutations on autosomal chromosomes including transgenic or knock-in assays that
are location specific (results are WT, HET or MUT):
    In some cases, especially breeding trios, the gender of the mouse is an important
    factor in a pass/fail result.

"""

import csv
from collections import defaultdict, namedtuple
import sys
import gzip
import requests
from pathlib import PurePath
from configparser import ConfigParser as CP
from argparse import ArgumentParser as AP
from itertools import chain
from openpyxl import load_workbook
import inheritance_rules as IR
import warnings
from musterer import get_musterer_mouse_info


def myopen(fn):
    """ Bob's function to handle gzip transparently """
    if fn.endswith('.gz') :
        return gzip.open(fn, "rt")
    return open(fn, errors="ignore")


class Observable:
    """ plain old data, supports conversion from old assays types to NGS """
    def __init__(self, old_family, old_assay, old_obs, ngs_family, ngs_assay, ngs_obs):
        self.old_family = old_family
        self.old_assay = old_assay
        self.old_observables = tuple(old_obs)  # the possible genotypes
        self.ngs_family = ngs_family
        self.ngs_assay = ngs_assay
        self.ngs_observables = tuple(ngs_obs)  # the possible genotypes


    def __str__(self):
        parts = ['Musterer family:', self.old_family, 'Musterer assay:', self.old_assay,
               'Musterer observables:', ';'.join(self.old_observables),
               'NGS family:', self.ngs_family, 'NGS assay:', self.ngs_assay,
               'NGS observables:', ';'.join(self.ngs_observables)]
        return ' '.join(parts)


    def __repr__(self):
        return 'Observable()'


    def get_old_obs(self):
        """ Needed because the old assays often have weird observables.
            Instead, look up the new observables and use the matchd indices to return the old ones
        """
        wt_ob = self.old_observables[[i for i,o in enumerate(self.ngs_observables) if 'wt' in o.lower() and not 'mut' in o.lower() and not 'y' in o.lower()][0]]
        mut_ob = self.old_observables[[i for i,o in enumerate(self.ngs_observables) if 'mut' in o.lower() and not 'wt' in o.lower() and not 'y' in o.lower()][0]]
        het_ob = self.old_observables[[i for i,o in enumerate(self.ngs_observables) if 'wt'  in o.lower() and 'mut' in o.lower() and not 'y' in o.lower()][0]]
        wtY_indices = [i for i,o in enumerate(self.ngs_observables) if 'wt' in o.lower() and not 'mut' in o.lower() and 'y' in o.lower()]
        mutY_indices = [i for i,o in enumerate(self.ngs_observables) if 'mut' in o.lower() and not 'wt' in o.lower() and 'y' in o.lower()]
        wtY_ob = None
        if wtY_indices:
            wtY_ob = self.old_observables[wtY_indices[0]]
        mutY_ob = None
        if mutY_indices:
            mutY_ob = self.old_observables[mutY_indices[0]]
        return wt_ob, het_ob, mut_ob, wtY_ob, mutY_ob


def _flush_conversion_buffer(buffer):
    """ yield one or more new Observable objects from a line buffer """
    for i, row in enumerate(buffer):
        if i == 0:
            old_assay = row[1]
            old_obs = [row[2]]
            old_family = row[3]
            ngs_assay = row[4]
            ngs_obs = [row[5]]
            ngs_family = row[6]
        elif row[1] == '' and row[3] == '' and row[4] == '' and row[6] == '':
            # add observables
            if row[2] != '':
                old_obs.append(row[2])
            if row[5] != '':
                ngs_obs.append(row[5])
        elif row[1] != '' or row[3] != '' or row[4] != '' or row[6] != '':
            # record info seen so far
            yield Observable(old_family, old_assay, old_obs, ngs_family, ngs_assay, ngs_obs)
            # new assay - potentially multiple NGS assays for each old assay
            if row[4] != '':
                if row[2] != '':
                    old_obs = [row[2]]
                if row[5] != '':
                    ngs_obs = [row[5]]
                ngs_assay = row[4]
            else:
                print('Multiple assays detected without new entry',i,row, file=sys.stderr)
    yield Observable(old_family, old_assay, old_obs, ngs_family, ngs_assay, ngs_obs)


def parse_conversions(fn):
    """ ugly Excel file for converting NGS assays to Amplifluor assay names and vice versa """
    old_families = defaultdict(list)  # dict[family] = list(assays)
    ngs_families = defaultdict(list)  # dict[family] = list(assays)
    old_assays = defaultdict(list)  # dict[assay] = [Observables]
    ngs_assays = defaultdict(list)  # dict[assay] = [Observables]
    buffer = []  # row buffer
    wb = load_workbook(fn, read_only=True)
    for sheet in wb:
        for i,rowobjs in enumerate(sheet.rows):
            row = [str(cell.value).replace('None','').strip() for cell in rowobjs]
            #print(row, file=sys.stderr)
            if i == 0:  # header
                continue
            if row[4] != '' and buffer != []:
                # start a new entry and populate dicts with buffer info
                for obs in _flush_conversion_buffer(buffer):
                    old_families[obs.old_family].append(obs)
                    ngs_families[obs.ngs_family].append(obs)
                    old_assays[obs.old_assay].append(obs)
                    ngs_assays[obs.ngs_assay].append(obs)
                buffer = [row]
            else:
                if ''.join(row).strip() == '':
                    continue
                buffer.append(row)
        for obs in _flush_conversion_buffer(buffer):
            old_families[obs.old_family].append(obs)
            ngs_families[obs.ngs_family].append(obs)
            old_assays[obs.old_assay].append(obs)
            ngs_assays[obs.ngs_assay].append(obs)
        buffer = []
        break
    print('Musterer to NGS conversions successfully read.\n', file=sys.stderr)
    print('Musterer family names:', ','.join(old_families), '\n', file=sys.stderr)
    print('Musterer assays:', ','.join(old_assays), '\n', file=sys.stderr)
    print('NGS family names:', ','.join(ngs_families), '\n', file=sys.stderr)
    print('NGS assays:', ','.join(ngs_assays), '\n', file=sys.stderr)
    return old_families, old_assays, ngs_families, ngs_assays


def get_old_new_assay(query, ngs_families, ngs_assays, old_families, old_assays):
    """ 
    There is a chance of getting a new NGS assay and needing to figure out the matching old assay name.
    Alternatively, get the new ngs_assay from the old assay or family name 
    """
    #print('Assay or family to look up', assay_or_fam, file=sys.stderr)
    #print('ngs_families', ngs_families, file=sys.stderr)
    #print('ngs_assays', ngs_assays, file=sys.stderr)
    #print('old_families', old_families, file=sys.stderr)
    #print('old_assays', old_assays, file=sys.stderr)
    ngs_query = 'NGS::' + query
    amplifluor_query = 'Amplifluor' + query
    if ngs_query in ngs_assays:
        a = ngs_assays[ngs_query][0].old_assay
        print('Found in ngs assays!', a, file=sys.stderr)
        return a
    elif ngs_query + query in ngs_families:
        a = ngs_families[ngs_query][0].old_assay
        print('Found in ngs families!', a, file=sys.stderr)
        return a
    elif amplifluor_query in old_assays:
        a = old_assays[amplifluor.query][0].ngs_assay
        print('Found in old assays!', a, file=sys.stderr)
        return a
    elif amplifluor_query in old_families:
        a = old_families[amplifluor_query][0].ngs_assay
        print('Found in old families!', a, file=sys.stderr)
        return a
    else:
        print('Cannot find', query, 'in known assay or family conversions', file=sys.stderr)
        return ''


def read_results(result_fn):
    """ return header and records from Bob's matching of sequences """
    with open(result_fn, 'rt') as srcfd:
        src = csv.reader(srcfd)
        hdr = next(src)
        rowlen = len(hdr)-2 # last two are first of repeated pairs
        Rec = namedtuple('Rec', hdr[:rowlen]+['matches'])
        def mkrow(r):
            g = (x for x in r[rowlen:])
            matches = list(zip(g,g))
            p1, p2 = r[:rowlen]+['']*max(0, rowlen-len(r)), [matches]
            return Rec(*(p1+p2))
        records = list(map(mkrow, src))
    return hdr, records


def standardise_allele(config, allele):
    """ take an observed or Musterer allele (half an observable) and try to convert it to standard terms """
    allele = allele.lower()
    if allele == 'other':
        return ''
    if allele == 'y':
        allele = '/y'
    mut_set = set(config['genotypes']['mut'].lower().split(','))
    wt_set = set(config['genotypes']['wt'].lower().split(','))
    y_set = set(config['genotypes']['y'].lower().split(','))
    pos_set = set(config['genotypes']['pos'].lower().split(','))
    neg_set = set(config['genotypes']['neg'].lower().split(','))
    multi_set = set(config['genotypes']['multi'].lower().split(','))

    if allele in y_set:
        return '/y'
    if allele in mut_set:
        return 'mut'
    elif allele in wt_set:
        return 'wt'
    elif allele in pos_set:
        return 'pos'
    elif allele in neg_set:
        return 'neg'
    elif allele in multi_set:
        return allele
    else:
        print('Unknown allele', allele, file=sys.stderr)
        return ''


def match_assay(my_assay, target_assays, ngs_families, ngs_assays, old_families, old_assays):
    """ try to find matching assays in a target set. Try converting old/new if required """
    my_assay = my_assay.lower()
    musterer_assays = [ta for ta in target_assays if ta.lower() in my_assay]
    if musterer_assays:
        #print('(1)Primers matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    musterer_assays = [ta for ta in target_assays if ta.lower().replace('-','') in my_assay.replace('-','')]
    if musterer_assays:
        #print('(2)Primers matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    musterer_assays = [ta for ta in target_assays if ta.lower().split('-')[0] in my_assay.split('-')[0]]
    if musterer_assays:
        #print('(3)Primers matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    
    # try at family level
    f = my_assay.split('_')[0].lower()
    musterer_assays = [ta for ta in target_assays if f in ta.lower()]
    if musterer_assays:
        #print('(1)Family matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    musterer_assays = [ta for ta in target_assays if f.replace('-','') in ta.lower().replace('-','')]
    if musterer_assays:  # sometimes just the first part of a name is all that is in common
        #print('(2)Family matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    musterer_assays = [ta for ta in target_assays if f.split('-')[0] in ta.split('-')[0].lower()]
    if musterer_assays:
        #print('(3)Family matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays

    # try converting to old/new assay
    musterer_assays = [get_old_new_assay(my_assay, ngs_families, ngs_assays, old_families, old_assays)]
    if musterer_assays:
        #print('(1)Assay conversion matching Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    musterer_assays = [get_old_new_assay(f, ngs_families, ngs_assays, old_families, old_assays)]
    if musterer_assays:
        #print('(1)Family conversion Musterer assays:', musterer_assays, file=sys.stderr)
        return musterer_assays
    print('No match found for', my_assay, file=sys.stderr)
    return []             


class MObs():
    """ Mouse Observations, encapsulates per-mouse results """
    def __init__(self, idx_rec, old_assays, ngs_assays): #assay_conversion=None):
        self.barcode = ''
        self.indexed_records = defaultdict(list)  # [idx] = record
        self.family_genotype = defaultdict(str)  # should be a dictionary of assay family:genotype
        self.family_to_musterer_assays = defaultdict(list)  # for each recorded family, these are the associated assay names
        #self.family_heredity_check = defaultdict(lambda: False)  # does the assay pass heredity sanity checks
        self.musterer_assay_genotype = defaultdict(str) # musterer assay name, associated with an observable as its genotype
        self.musterer_assay_heredity_check = defaultdict(lambda: False)  # does each assay pass its heredity checks
        self.family_indices = defaultdict(list)
        self.mouseAssays = []  # directly from Rec.mouseAssays
        self.indexed_primer = defaultdict(str)   # [idx] = primer
        self.family_sanity_pass = defaultdict(str)  # dictionary assay:true/false
        self.family_sanity_comment = defaultdict(str)  # dictionary assay:''
        self.family_allele_ratio = defaultdict(str)  # propagated to self.indexed_allele_ratio
        self.indexed_matches = defaultdict(list)  # all observed matches within a well
        #self.indexed_sum_counts = defaultdict(lambda: 0)  # for calculating efficiency, etc
        #self.indexed_primer_conversion = defaultdict(None)  # one assay conversion object per index
        self.indexed_allele_ratio = defaultdict(lambda: ['NA'])
        self.indexed_efficiency = defaultdict(lambda: 'NA')
        self.indexed_contaminated = defaultdict(lambda: False)  # [i] = True/False
        #self.indexed_seqs = defaultdict(list)  # the sequence names that match the expected primer
        #self.indexed_empty = defaultdict(lambda: False)  # True if no observed DNA fragments
        self.indexed_new = defaultdict(lambda: False)
        self.indexed_family_match = defaultdict(lambda: False)
        self.indexed_genotype = defaultdict(str)
        self.indexed_musterer_assays = defaultdict(list)
        self.indexed_chosen_matches = defaultdict(list)  # contains list of (count,assay) which match the primer
        self.upload_records = []  # created on calling genotypes
        self.checked_upload_records = []  # filtered on heredity checking
        if idx_rec:
            self.barcode = idx_rec[1].mouseBarcode
            self.upload_identity = idx_rec[1].EPplate + '_' + idx_rec[1].EPwell +\
                    '_' + idx_rec[1].strainName + '#' + idx_rec[1].mouseID
            self.add_indexed_record(idx_rec, old_assays, ngs_assays) #assay_conversion)
            #if not success:
            #    print('Failed to add indexed record', idx_rec, file=sys.stderr)
        #else:
        #    print('Successfully added indexed record', idx_rec, file=sys.stderr)
        # these three are for data we pull from Musterer
        self.musterer_info = defaultdict(None)
        self.sire_barcode = None
        self.sire_info = defaultdict(None)  # [sire_bc] = {dict:sire_info}
        self.sire_gt = None
        self.dams_barcodes = None
        self.dams_info = defaultdict(None)  # [dam_bc] = {dict:dam_info}
        self.dam_bc_gt = None  # barcode of dam that passed sanity check

    def __str__(self):
        return '\t'.join([self.barcode, 'records:', str(len(self.indexed_records)), 'matches:', str(len(self.indexed_matches))])


    # def get_parents_strains(self):
    #     """ return a list of tuples, sire_bc/strain, dam_bc/strain, ... """
    #     parents_strains = []
    #     for sire_bc in self.sire_info:
    #         parents_strains.append((sire_bc, self.sire_info[sire_bc][strain]))
    #     for dam_bc in self.dams_info:
    #         parents_strains.append((dam_bc, self.dams_info[dam_bc][strain]))
    #     return parents_strains

    #def _assays_found(self, idx, our_assays, old_assays, ngs_assays):
    #    """ common code to check if we have correctly identified the expected assay """

    #    # attempt to collapse any assays in both ngs and old formats to ngs only
    #    #print('Original assays', our_assays, file=sys.stderr)
    #    our_assays2 = set()
    #    for oa in our_assays:
    #        if oa.replace('NGS::','') in our_assays2:
    #            continue
    #        old_name = 'Amplifluor::' + oa
    #        if old_name.lower() in old_assays:
    #            ngs_a = old_assays[old_name.lower()][0].ngs_assay
    #            #print(ngs_a, file=sys.stderr)
    #            if ngs_a.replace('NGS::','').lower() in our_assays2:
    #                continue
    #            our_assays2.add(ngs_a.replace('NGS::',''))
    #            continue
    #        our_assays2.add(oa)
    #    selected_assays = list(our_assays2)

    #    #if len(selected_assays) == 0:
    #    #    print('No assays found from:', our_assays, file=sys.stderr)
    #    #    return False

    #    #print('Collapsed assays', selected_assays, file=sys.stderr)

    #    rec = self.indexed_records[idx]
    #    #if idx == 7:
    #    #    print('Our assays:', our_assays, rec.primer, file=sys.stderr)
    #    if len(selected_assays) == 1:
    #        # name = 'Amplifluor::' + our_assays[0].split('_')[0]
    #        # if name.lower() in old_assays:
    #        #     for oa in old_assays[name.lower()]: # list of conversion objects... probably shouldn't be
    #        #         self.indexed_primer[idx] = our_assays[0]
    #        # else:
    #        self.indexed_primer[idx] = selected_assays[0]
    #        return True
    #    elif len(selected_assays) > 1: # 2 ore more possible assays
    #        # is there more primer info we can use to distinguish the correct assay
    #        if '_' in rec.primer:
    #            primer_parts = rec.primer.split('_')
    #            for pp in primer_parts[1:]:
    #                for assay in selected_assays:
    #                    if pp.lower() in assay.lower():
    #                        self.indexed_primer[idx] = assay
    #                        return True
    #    else:
    #        for oa in our_assays:
    #            if 'ngs::' in oa.lower():
    #                if rec.primer.lower() in oa.lower():
    #                    self.indexed_primer[idx] = oa
    #                    return True
    #        # for oa in our_assays:
    #        #     if rec.primer.lower() in oa.lower():
    #        #         self.indexed_primer[idx] = oa
    #        #         return True
    #    # Or... is one an NGS assay?
    #    # ngs = []
    #    # for a in our_assays:
    #    #     #if idx==7:
    #    #     #    print(a, ngs_assays, file=sys.stderr)
    #    #     if 'NGS::'+a in ngs_assays:
    #    #         ngs.append(a)
    #    #     name = 'Amplifluor::' + a
    #    #     if name.lower() in old_assays:
    #    #         for oa in old_assays[name.lower()]: # list of conversion objects... probably shouldn't be
    #    #             ngs.append(oa.ngs_assay)
    #    # if len(ngs) == 1:
    #    #     self.indexed_primer[idx] = a
    #    #     return True
    #    # if len(ngs) > 1:
    #    #     print('multiple NGS assays... still', file=sys.stderr)
    #    # okay I'm out of ideas
    #    return False


    def add_indexed_record(self, idx_rec, old_assays, ngs_assays, minreads=50): # assay_conversion=None, minreads=100):
        """
            Add a record for the given index and filter assays/counts
            Record the filtered assay matches by index
            Record the allele ratio by index
        """
        idx = idx_rec[0]
        rec = idx_rec[1]
        self.indexed_records[idx] = rec
        if len(self.mouseAssays) == 0:
            self.moouseAssays = rec.mouseAssays
        #self.indexed_primer_conversion[idx] = assay_conversion
        # group together counts and their assay reference name
        self.indexed_matches[idx] = []
        self.indexed_primer[idx] = rec.primer
        matches = rec.matches
            
        if not matches:
            print('No matchable sequence reads found for record', idx, file=sys.stderr)
            return

        counts_assays = sorted([(int(c),a) for c,a in matches], reverse=True)
        all_counts = [c for (c,a) in counts_assays]
        if not all_counts or sum(all_counts) < minreads:
            print('Insufficient sequence reads found for record:', idx, 'Minimum required reads:', minreads, file=sys.stderr)
        
        self.indexed_matches[idx] = counts_assays
        self.indexed_primer[idx] = rec.primer
        #return True


    def call_genotypes(self, ngs_families, ngs_assays, old_families, old_assays, config):
        """ Establish which assay we are using then compare the sequenced allele tags against the available observables """
        
        # gather rows belonging to the same assay family
        for i in self.indexed_primer:
            fam = self.indexed_primer[i].split('_')[0]
            if not fam:
                continue
            self.family_indices[fam].append(i)
        
        for f in sorted(self.family_indices):
            print('family being genotyped',f, file=sys.stderr)
            all_assays = [p for p in [self.indexed_primer[i] for i in self.family_indices[f]]]
            if len(all_assays) == 0:
                print('No assay primers found for family:', f, file=sys.stderr)
                continue 
            all_musterer_assays = self.musterer_info['assay_value_options'].keys()
            # find all primer assays in musterer assay options
            musterer_assays = list(set(chain(*[match_assay(aa, all_musterer_assays, ngs_families, ngs_assays, old_families, old_assays) for aa in all_assays])))
            print('Musterer assays for this mouse:', all_musterer_assays, file=sys.stderr)

            # apply count filtering here and calculate efficiency, allele ratio, etc
            #all_matches = [m for m in [self.indexed_matches[i] for i in self.family_indices[f]]]
            family_matches = []
            for i in self.family_indices[f]:
                primer_matches = []
                well_total = sum([c for c,a in self.indexed_matches[i]])
                well_matches = []
                for c,a in self.indexed_matches[i]:  # m = (count, assay)
                    # efficiency is calculated per-well
                    # assay ratio is calculated across all wells for all valid alleles
                    # valid alleles are found per-well and combined together
                    # alleles are already sorted by decreasing count
                    if f.lower() in a.lower():
                        if not well_matches:
                            well_matches.append((c,a))
                            self.indexed_chosen_matches[i].append((c,a))
                        elif c > 0.3*(sum([c for c,a in well_matches])/len(well_matches)):
                            # catch underperforming alleles with 1/3rd the counts of others
                            well_matches.append((c,a))
                            self.indexed_chosen_matches[i].append((c,a))
                if well_total:
                    efficiency = 100*sum([c for c,a in well_matches])/well_total
                    self.indexed_efficiency[i] = f"{efficiency:.2f}%"
                family_matches += well_matches

            standard_alleles_present = {standardise_allele(config, a.split('_')[1]):c for c,a in family_matches}
            standard_alleles_present = {sa:standard_alleles_present[sa] for sa in standard_alleles_present if sa != ''}
            total_counts_present = sum([c for c,a in family_matches])

            print('Standardised alleles present', standard_alleles_present, file=sys.stderr)

            for assay in musterer_assays:
                assay = assay_long.replace('Amplifluor::','').replace('NGS::','')
                if assay not in self.musterer_info['assay_value_options']:
                    print("This was meant to come from self.musterer_info['assay_value_options']... the code is wrong!", assay, 
                            self.musterer_info['assay_value_options'], file=sys.stderr)
                    continue
                obs = self.musterer_info['assay_value_options'][assay]
                sex_linked = False
                if any(['/y' in o.lower() for o in obs]):
                    sex_linked = True
                print('Assay:',assay, 'Observables:',obs, file=sys.stderr)
                for ob in obs:
                    std_ob_list = [standardise_allele(config, ob)]
                    print('ob', ob, 'obs', obs, file=sys.stderr)
                    if not '/' in ob:
                        # CRE or other single insert mutations
                        o = standardise_allele(config, ob)
                        # is any match good enough, or does it need to be the ONLY match?
                        if o in standard_alleles_present or (o == 'neg' and not standard_alleles_present):
                            # if CREneg then there won't be any alleles present
                            self.family_to_musterer_assays[f].append(assay)
                            self.musterer_assay_genotype[assay] = ob
                            if self.family_genotype[f]:
                                self.family_genotype[f] = self.family_genotype[f] + '/' + ob
                                allele_ratios = [1/self.family_genotype[f].count('/')+1]
                                self.family_allele_ratio[f] = [f"{ar:.2f}" for ar in allele_ratios]
                            else:
                                self.family_genotype[f] = ob
                                self.family_allele_ratio[f] = ['1.00']
                            break 
                    elif ob.count('/') == 1:
                        # standard wt/mut style assay
                        o1, o2 = ob.split('/')
                        so1 = standardise_allele(config, o1)
                        so2 = standardise_allele(config, o2)
                        print ('Ob1:',so1, 'Ob2:',so2, 'Standard alleles present:',standard_alleles_present, file=sys.stderr)
                        # deal with sex-linked /Y genotypes
                        if (so1 == '/y' or so2 == '/y') and self.musterer_info['sex'] == 'M' and \
                                len(standard_alleles_present) == 1 and list(standard_alleles_present.keys())[0] in [so1,so2]:
                            self.family_to_musterer_assays[f].append(assay)
                            self.musterer_assay_genotype[assay] = ob
                            self.family_genotype[f] = ob
                            self.family_allele_ratio[f] = ['1.00']
                            break
                        # test that we have a bijection 
                        elif all([sa in set([o1,o2]) for sa in standard_alleles_present]) and all([o in standard_alleles_present for o in [o1,o2]]):
                            if self.musterer_info['sex'] == 'M' and any(['/y' in o.lower() for o in obs]) and '/y' not in ob.lower():
                                #print('Skipping', file=sys.stderr)
                                continue
                            self.family_to_musterer_assays[f].append(assay)
                            self.musterer_assay_genotype[assay] = ob
                            self.family_genotype[f] = ob
                            allele_ratios = [standard_alleles_present[sa]/total_counts_present for sa in standard_alleles_present]
                            self.family_allele_ratio[f] = [f"{ar:.2f}" for ar in allele_ratios]
                            break
                        #else:
                        #    print('Missed:', ob, 'with:', standard_alleles_present, file=sys.stderr)
                    else:
                        # multi-allelic fun
                        all_o = ob.split('/')
                        all_standard_o = set([standardise_allele(config, o) for o in all_o])
                        if all([sa in all_standard_o for sa in standard_alleles_present]) and \
                                all(o in standard_alleles_present for o in all_standard_o):
                            self.family_to_musterer_assays[f].append(assay)
                            self.musterer_assay_genotype[assay] = ob
                            self.family_genotype[f] = ob
                            allele_ratios = [standard_alleles_present[sa]/total_counts_present for sa in standard_alleles_present]
                            self.family_allele_ratio[f] = [f"{ar:.2f}" for ar in allele_ratios]
                            break
            if not self.family_genotype[f]:
                self.family_genotype[f] = 'unknown'
                print('No genotype found for assay family',f, 'Family matches:', family_matches, 'Standard alleles present:', 
                        standard_alleles_present, file=sys.stderr)
            
            #print(self.family_to_musterer_assays[f], self.family_genotype[f], alleles_present, file=sys.stderr)
            
            # set all indexed genotypes and indexed_musterer_assay to be the same for each record in the same assay family
            for i in self.family_indices[f]:
                if f in self.family_genotype: 
                    self.indexed_genotype[i] = self.family_genotype[f]
                    self.indexed_musterer_assays[i] = self.family_to_musterer_assays[f]
                    self.indexed_allele_ratio[i] = self.family_allele_ratio[f]
                else:
                    print('Mismatched assay family name:', f, self.family_genotype.keys(), file=sys.stderr)

        ### This should be done AFTER sanity checking
        # Create upload records based on accepted musterer assays
        for assay in self.musterer_assay_genotype:
            if self.musterer_assay_genotype[assay] not in set(['unknown','ambig']):
                self.upload_records.append(self.upload_identity + '\t' + assay + '\t' + self.musterer_assay_genotype[assay])
                   
            


    def check_family_match(self):
        """ aka Sanity Checking
        we need: mMA, mgt, msex, sire_gt, dam_gt

        Make sure mMA, mgt, msex make sense given the sire and dams
        Only one dam is required to pass for this mouse to pass the checks

        Filter self.upload_records into self.checked_upload_records
        """
        msex = self.musterer_info['sex']
        for fam in self.family_genotypes:
            for i in self.family_indices[fam]:
                mgtl = self.indexed_genotype[i].lower()
                sire_assays = {a for a in self.sire_info['assay_names_values'] \
                        if a.split('_')[0] == self.indexed_records.primer.split('_')[0]}
                dam_assays = {a for a in chain(*[self.dams_info[dbc]['assay_names_values'] \
                                 for dbc in self.dams_barcodes]) \
                              if a.split('_') == self.indexed_records.primer.split('_')[0]}
                for sa in sire_assays:
                    pass
                for da in dam_assays:
                    pass

            obs_gt = self.genotypes[expected_assay]
            sgtl = self.sire_info['assay_names_values'][expected_assay]

            mgtl = obs_gt.lower()
            # Do sex linked assays first
            if '/y' in mMA.lower():
                if mgtl == 'wt/wt' or mgtl == 'wt/y':
                    self.family_pass[ea] = 'Pass' if IR.test_true_XY_wt(msex, sgtl, dgtl) else 'Fail'
                elif mgtl == 'mut/wt' or mgtl =='wt/mut':
                    sanity_result = 'Pass' if IR.test_true_XY_het(msex, sgtl, dgtl) else 'Fail'
                elif mgtl == 'mut/mut' or mgtl =='mut/y':
                    sanity_result = 'Pass' if IR.test_true_XY_mut(msex, sgtl, dgtl) else 'Fail'
                if sanity_result == 'Fail':
                    sanity_comment = 'Sex impossible'
            # now handle regular simplex assays
            elif mgtl == 'mut/mut':
                sanity_result = 'Pass' if IR.test_true_mut(sgtl, dgtl) else 'Fail'
            elif mgtl == 'mut/wt' or mgtl == 'wt/mut':
                sanity_result = 'Pass' if IR.test_true_het(sgtl, dgtl) else 'Fail'
            elif mgtl == 'wt/wt':
                sanity_result = 'Pass' if IR.test_true_wt(sgtl, dgtl) else 'Fail'
            elif 'pos' in mgtl:
                sanity_result = 'Pass' if IR.test_true_pos(sgtl, dgtl) else 'Fail'
            elif 'neg' in mgtl:
                sanity_result = 'Pass' if IR.test_true_neg(sgtl, dgtl) else 'Fail'
            return sanity_result, sanity_comment


    def write_upload_records_to_CSV(self, csv_fh):
        """
        Inputs: Open file handle (tab delimited entries)
        Column #1: mouse 
        Column #2: Asssay name
        Column #3: result
 
        25708_6C_ENU24:017:Tmem131:B6:G6#3  Amplifluor::01-36811427-IGL561 (Tmem131) mut/mut
        25708_6D_ENU24:017:Tmem131:B6:G6#4   Amplifluor::01-36811427-IGL561 (Tmem131)  mut/wt
 
        In the first column where we identify the mouse the breakdown is as follows:
        plate_wellPosition_strainName#mouseNumber
        """
        for ur in self.upload_records:
            print(ur, file=csv_fh)


    def upload_records_to_musterer():
        """ How do we interact with this web service? """
        pass


def gather_mouse_results(records, ngs_assays, old_assays, ngs_families, old_families):
    """ combine records on a per-mouse basis and return a dictionary of mouse objects """
    mice = defaultdict(None)  # [barcode] = MObs()
    for i,r in enumerate(records):
        if r.mouseBarcode in mice:
            mice[r.mouseBarcode].add_indexed_record((i,r), old_assays, ngs_assays)
        else:
            mice[r.mouseBarcode] = MObs((i,r), old_assays, ngs_assays)
    return mice


def output_results_csv(out_fn, records, mice, in_hdr, ngs_families, ngs_assays, old_families, old_assays):
    """ output a header and all mouse records to CSV file, matching results.csv input """
    # build new header for output
    out_hdr = ','.join(in_hdr[:23] + ['gt','alleleRatio',
            'efficiency','sanity_result', 'sanity_comment','sire_barcode','sire_strain', 'sire_gt',
            'dam_barcode','dam_strain','dam_gt'] + in_hdr[23:])
    #print(out_hdr, file=sys.stderr)
    with open(out_fn, 'wt') as out:
        out.write(out_hdr + '\n')
        for i,r in enumerate(records):
            mbc = r.mouseBarcode
            m = mice[mbc]
            gt = m.indexed_genotype[i]
            assays = m.indexed_musterer_assays[i]
            allele_ratio = ':'.join(map(str, m.indexed_allele_ratio[i]))
            # TODO: efficiency needs to be based off mergedCount
            efficiency = m.indexed_efficiency[i]
            #sanity_result = m.indexed_family_match[i]
            sire_bc = m.sire_barcode
            
            def set_sire_strain_gt(sire_info, assays, ngs_families, ngs_assays, old_families, old_assays):
                sire_gt = 'NA'
                sire_strain = 'NA'
                if 'strain' in sire_info:
                    sire_strain = m.sire_info['strain']
                if 'assay_names_values' in sire_info:
                    sire_assays = sire_info.get('assay_names_values')
                    # Will need to match Ptprc against Ly5A/B
                    # sire_assays = chain(*[match_assay(sa,assays, ngs_families, ngs_assays, old_families, old_assays) for sa in sire_assays])
                    for a in assays:
                        if a in sire_assays:
                            sire_gt = sire_info['assay_names_values'][a]
                            break
                return sire_strain, sire_gt
            if not m.sire_info:
                sire_strain = 'NA'
                sire_gt = 'NA'
            else:
                sire_strain, sire_gt = set_sire_strain_gt(m.sire_info, assays, ngs_families, ngs_assays, old_families, old_assays)
            
            #dam_bc = m.dam_bc_gt # this needs to be added
            dam_bcs = m.dams_barcodes

            dams_strain = ';'.join(map(str, [m.dams_info[dam_bc]['strain'] for dam_bc in dam_bcs if m.dams_info[dam_bc] != 'NA']))
            try:
                dams_gt = ';'.join(map(str, [m.dams_info[dam_bc]['assay_names_values'][m.indexed_primer[i]] for dam_bc in dam_bcs]))
            except:
                dams_gt = 'NA'

            #print('Record:', i, r, file=sys.stderr)
            outcols = list(r[:23]) + [gt, allele_ratio, efficiency, '', '', sire_bc,
                                      sire_strain, sire_gt, ';'.join(map(str,dam_bcs)), dams_strain, dams_gt] + \
                    list(r[23:26]) + list(chain(*r[26]))
            outline = ','.join(map(str, outcols))
            #print(outcols, file=sys.stderr)
            print(outline, file=out)


def output_mouse_gts(mice_gt_fn, mice):
    """ output a CSV of mouse barcodes with genotypes
    mouseID	EPplate	EPwell	mouseBarcode	strainName	mouseAssays	assayFamilies	dnaplate	dnawell	primer	pcrplate	pcrwell	gt	alleleRatios	efficiency	seqCount	seqName		
    33	71561	A5	105210000027	ASD792:RCE:LoxP::F1rF4N6	RCEloxP	RCE-LoxP	20210706	A3	RCE-LoxP_MUT; RCE-LoxP_WT	3121; 3121	H1; I1	mut/wt			5193; 1715 & 7 	RCE-LoxP_MUT_Neo_100bp; RCE-LoxP_WT_100bp & RCE-LoxP_MUT_Neo_100bp		
    197	71590	B2	105366300045	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	B13	Apc	3121	F5	mut/wt	NA	NA	4128 & 3037	Apc_WT_120bp & Apc_MUT_PM_TtoA_120bp		
    197	71590	B2	105366300045	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	B13	Ku70_MUT; Ku70_WT	3121; 3121	G5; H5 	wt/wt	NA	NA	; 3600	; Ku70_WT_86bp		no reads found for mutant reaction
    193	71590	A2	105366300050	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	A13	Apc	3121	M4	mut/wt	NA	NA	1710 & 1478	Apc_WT_120bp & Apc_MUT_PM_TtoA_120bp		
    193	71590	A2	105366300050	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	A13	Ku70_MUT; Ku70_WT	3121; 3121	N4; O4	Ku70 KO/wt	NA	NA	1176 & 8 ; 1591	Ku70_MUT_Neo_86bp & Ku70_WT_86bp ; Ku70_WT_86bp		
    """
    mice_gt_fn = 'mice_gts.csv'
    hdrs = 'mouseID EPplate EPwell mouseBarcode strainName mouseAssays assayFamilies dnaplate '+\
            'dnawell primer pcrplate pcrwell gt alleleRatios efficiency seqCount seqName'
    col_hdrs = hdrs.split(' ')
    
    #print(out_hdr, file=sys.stderr)
    with open(mice_gt_fn, 'wt') as out:
        print('\t'.join(col_hdrs), file=out)
        for mouse in mice:
            m = mice[mouse]
            for f in m.family_indices:
                info = defaultdict(list)
                for x,i in enumerate(m.family_indices[f]):
                    if x == 0:
                        info['mouseID'] = m.indexed_records[i].mouseID
                        info['EPplate'] = m.indexed_records[i].EPplate
                        info['EPwell'] = m.indexed_records[i].EPwell
                        info['mouseBarcode'] = m.indexed_records[i].mouseBarcode
                        info['strainName'] = m.indexed_records[i].strainName
                        info['mouseAssays'] = m.indexed_records[i].mouseAssays
                        info['assayFamilies'] = m.indexed_records[i].assayFamilies
                        info['dnaplate'] = m.indexed_records[i].dnaplate
                        info['dnawell'] = m.indexed_records[i].dnawell
                    info['primer'].append(m.indexed_records[i].primer)
                    info['pcrplate'].append(m.indexed_records[i].pcrplate)
                    info['pcrwell'].append(m.indexed_records[i].pcrwell)
                    info['gt'].append(m.indexed_genotype[i])
                    if x == 0:
                        info['alleleRatios'] = ':'.join(m.indexed_allele_ratio[i])
                        info['efficiency'] = m.indexed_efficiency[i]
                    info['seqCount'].append('; '.join([str(c) for c,a in m.indexed_chosen_matches[i]]))
                    info['seqName'].append('; '.join([a for c,a in m.indexed_chosen_matches[i]]))
                # add placeholders for missing data
                for field in info:
                    if not info[field]:
                        info[field] = 'NA'
                    
                out_line = '\t'.join([' & '.join(map(str,info[col])) if type(info[col]) is list else str(info[col]) for col in col_hdrs])
                print(out_line, file=out)




def main(result_fn, reference_fn, config_fn, conversions_fn, out_fn=None, mice_fn=None, upload_fn=None):
    """
        Read reference file
        Read results file...
        For each row in the results table:
            1) extract alleles and counts
            2) decide the genotype
            3) sanity check the genotype
            4) output per-well table
            5) output per-mouse table
            6) output upload table
    """

    warnings.filterwarnings("ignore")
    config = CP()
    config.read(PurePath(config_fn))

    # read Amplifluor/NGS assay name conversion document in different indexings
    old_families, old_assays, ngs_families, ngs_assays = parse_conversions(conversions_fn)
    in_hdr, records = read_results(result_fn)

    print('Parsing mouse results...', file=sys.stderr)
    mice = gather_mouse_results(records, ngs_assays, old_assays, ngs_families, old_families)
    
    # look up mouse assays from Musterer by mouse barcode and add this to each mouse object
    print('Looking up Musterer records...', file=sys.stderr)
    mouse_barcodes = {r.mouseBarcode for r in records}
    mouse_musterer_info = get_musterer_mouse_info(list(mouse_barcodes))
    if not mouse_musterer_info:
        print('No Musterer info', file=sys.stderr)
    else:
        # collect parent info now, so we can do all checking in one pass
        print('Collecting mouse sire/dam records from Musterer....', file=sys.stderr)
        sires_barcodes = {mouse_musterer_info[mbc]['sire_barcode'] for mbc in mouse_musterer_info}
        sires_info = get_musterer_mouse_info(list(sires_barcodes))

        # need a little magic to flatten out a list of lists since there are multiple dams
        dams_barcodes = list(chain(*[mouse_musterer_info[mbc]['dams_barcodes'] \
                                 for mbc in mouse_musterer_info]))
        dams_info = get_musterer_mouse_info(list(dams_barcodes))

        # add info to each experiment mouse object
        for mbc in sorted(mice):  # sorted by barcode
            m = mice[mbc]
            m.musterer_info = mouse_musterer_info[mbc]
            m.sire_barcode = m.musterer_info['sire_barcode']
            m.sire_info = sires_info[m.sire_barcode]  
            m.dams_barcodes = m.musterer_info['dams_barcodes']
            m.dams_info = {qdbc:dams_info[qdbc] for qdbc in m.dams_barcodes}

    #  call genotypes
    print('Genotyping mice', file=sys.stderr)
    #loggable = []
    
    if upload_fn:
        ur = open(upload_fn, 'wt')
    for j, mbc in enumerate(mice):
        m = mice[mbc]
        print(m, file=sys.stderr)
        print(j+1,'/',len(mice),mbc, file=sys.stderr)
        #print('matches and assay', [(i, m.indexed_matches[i], m.indexed_primer[i]) \
        #        for i in m.indexed_matches], file=sys.stderr)
        
        #debug = m.call_genotypes(ngs_families, ngs_assays, old_families, old_assays, config)
        #loggable.append(debug)
        m.call_genotypes(ngs_families, ngs_assays, old_families, old_assays, config)
        print('GTs by index:', [(i, m.indexed_genotype[i]) for i in m.indexed_records], file=sys.stderr)
        #if mbc == '105371500058':
        #    print(m.sire_barcode, m.sire_info)
        #    print(m.dams_info)
        #    sys.exit(0)
        #print(file=sys.stderr)
        #m.check_family_match()
        if upload_fn:
            m.write_upload_records_to_CSV(ur)
    ur.close()

    output_results_csv(out_fn, records, mice, in_hdr, ngs_families, ngs_assays, old_families, old_assays)
    output_mouse_gts('mice_gts.csv', mice)
    print('\nGenotyping complete.', file=sys.stderr)


if __name__ == '__main__':
    parser = AP()
    parser.add_argument('resultfile', help='Path to Results.csv file')
    parser.add_argument('-r', '--reference', required=True, help='Path to reference sequence file')
    parser.add_argument('-c', '--conversions', required=True, help='Path to NGS observables '+\
                        'conversionsassay description file')
    parser.add_argument('-o', '--outfile', default='results_gt.csv', help='Path to output file. If not present will '+\
                        'output to stdout')
    parser.add_argument('-m', '--mice', help='Path to output per-mouse genotype CSV')
    parser.add_argument('-k','--config',default='../bin/config.ini',help='Path to config file')
    parser.add_argument('-u','--upload', help='Path to create file for upload of genotypes to Musterer'
                        , default='upload_records.csv')
    args = parser.parse_args()
    main(args.resultfile, args.reference, args.config, args.conversions, out_fn=args.outfile, mice_fn=args.mice, upload_fn=args.upload)
