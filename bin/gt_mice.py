#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: Dec 2020
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.7
@version_comment: finally incorporated
@last_edit: 2021-07-02
@edit_comment: all test cases genotyping correctly

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

We need to output 2 colums: "call" and "GT". Call is a row-based decision based on the
observed counts vs the expected primer, while GT is of the form that matches what is expected 
by Musterer.
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
        self.old_observables = tuple(old_obs)  # the per-reaction calls/gts
        self.ngs_family = ngs_family
        self.ngs_assay = ngs_assay
        self.ngs_observables = tuple(ngs_obs)  # the per-reaction calls/gts


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


class MObs():
    """ Mouse Observations, encapsulates per-mouse results """
    def __init__(self, idx_rec, old_assays, ngs_assays): #assay_conversion=None):
        self.barcode = ''
        self.indexed_records = defaultdict(list)  # [idx] = record
        self.family_genotype = defaultdict(list)  # should be a dictionary of assay family:genotype
        self.family_to_musterer_assays = defaultdict(list)  # for each recorded family, these are the assays names to go with the genotypes
        self.family_heredity_check = defaultdict(lambda: False)  # does the assay pass heredity sanity checks
        self.musterer_assay_genotype = defaultdict(str) # musterer assay name, associated with an observable as its genotype
        self.musterer_assay_heredity_check = defaultdict(lambda: False)
        self.family_indices = defaultdict(list)
        self.mouseAssays = []  # directly from Rec.mouseAssays
        self.indexed_primer = {}   # [idx] = primer
        self.family_pass = defaultdict(str)  # dictionary assay:true/false
        self.family_comment = defaultdict(str)  # dictionary assay:''
        self.indexed_reactions = defaultdict(list)  # all observed reactions within a well
        self.indexed_sum_counts = defaultdict(lambda: 0)  # for calculating efficiency, etc
        #self.indexed_primer_conversion = defaultdict(None)  # one assay conversion object per index
        self.indexed_allele_ratios = defaultdict(lambda: ['NA'])
        self.indexed_efficiencies = defaultdict(lambda: 'NA')
        self.indexed_contaminated = defaultdict(lambda: False)  # [i] = True/False
        self.indexed_calls = defaultdict(str)  # string 'plus'/'minus'/'new'/''
        self.indexed_seqs = defaultdict(list)  # the sequence names that match the expected primer
        self.indexed_empty = defaultdict(lambda: False)  # True if no observed DNA fragments
        self.indexed_new = defaultdict(lambda: False)
        self.indexed_family_match = defaultdict(lambda: False)
        self.indexed_genotype = defaultdict(str)
        self.upload_records = []  # created on calling genotypes
        self.checked_upload_records = []  # filtered on heredity checking
        if idx_rec:
            self.barcode = idx_rec[1].mouseBarcode
            self.upload_identity = idx_rec[1].EPplate + '_' + idx_rec[1].EPwell +\
                    '_' + idx_rec[1].strainName + '#' + idx_rec[1].SampleNo
            self.add_indexed_record(idx_rec, old_assays, ngs_assays) #assay_conversion)
            #if not success:
            #    print('Failed to add indexed record', idx_rec, file=sys.stderr)
        #else:
        #    print('Successfully added indexed record', idx_rec, file=sys.stderr)
        # these three are for data we pull from Musterer
        self.musterer_info = None
        self.sire_barcode = None
        self.sire_info = defaultdict(lambda: 'NA')  # [sire_bc] = {dict:sire_info}
        self.dams_barcodes = None
        self.dams_info = defaultdict(lambda: 'NA')  # [dam_bc] = {dict:dam_info}
        self.dam_bc_gt = None  # barcode of dam that passed sanity check

    def __str__(self):
        return '\t'.join([self.barcode, 'records:', str(len(self.indexed_records)), 'reactions:', str(len(self.indexed_reactions))])


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
            Record the filtered assay reactions by index
            Record the allele ratios by index
        """
        idx = idx_rec[0]
        rec = idx_rec[1]
        self.indexed_records[idx] = rec
        if len(self.mouseAssays) == 0:
            self.moouseAssays = rec.mouseAssays
        #self.indexed_primer_conversion[idx] = assay_conversion
        # group together counts and their assay reference name
        self.indexed_reactions[idx] = []
        self.indexed_primer[idx] = rec.primer
        matches = rec.matches
            
        if not matches:
            print('No matchable sequence reads found for record', idx, file=sys.stderr)
            return

        counts_assays = sorted([(int(c),a) for c,a in matches], reverse=True)
        all_counts = [c for (c,a) in counts_assays]
        if all_counts:
            max_count = max(all_counts)

        if sum(all_counts) < minreads:
            print('Insufficient sequence reads found for record:', idx, 'Minimum required reads:', minreads, file=sys.stderr)
            return

        # filter, only keep assays in meaningful proportions
        counts_assays = [(c,a) for c,a in counts_assays \
                if c > max_count/10 and c > minreads/len(all_counts)]
        if len(counts_assays) == 0:  # can this even happen? Probably not...
            print('No sequence reads found of significant proportion', idx, file=sys.stderr)
            return
        
        self.indexed_reactions[idx] = counts_assays
        self.indexed_primer[idx] = rec.primer
        #return True


    def call_reactions(self, minreads=100):
        """
            'Call' whether we have successfully observed a primer sequence, based on rec.primer
            Record all indexed calls as 'plus' or 'minus'
        """   #
        #if self.musterer_info is None:
        #    print('Musterer assay information unavailable, calling impossible', file=sys.stderr)
        #    sys.exit(-1)
        # If assay is only unknown sequence then it's 'new'
        #dna_chars = {'A','C','G','T'}
        for i in self.indexed_reactions:
            self.indexed_calls[i] = 'minus'
            assay = self.indexed_primer[i]

            observed_assays = [a for (c,a) in self.indexed_reactions[i]]
            # check that both the assay family name and the specific variation e.g. Apc_WT, are present in at least 1 sequence name
            for oa in observed_assays:
                if all(part.lower() in oa.lower() for part in assay.split('_')):
                    self.indexed_calls[i] = 'plus'
                    break

    def call_obs(self):
        """
            'Call' the observed alleles if the reference sequence name is related to the primer.
            Record the sequence names that match
        """
        for i in self.indexed_reactions:
            expected_assay = self.indexed_primer[i]
            observed_assays = [a for (c,a) in self.indexed_reactions[i]]
            # check that both the assay family name and the specific variation e.g. Apc_WT, are present in at least 1 sequence name
            for oa in observed_assays:
                if all(part.lower() in oa.lower() for part in expected_assay.split('_')):
                    self.indexed_seqs[i].append(oa)


    def call_genotypes(self, old_families, ngs_families, old_assays, ngs_assays, config):
        """ Establish which assay we are using then compare the sequenced allele tags against the available observables """

        for i in self.indexed_primer:
            fam = self.indexed_primer[i].split('_')[0]
            if not fam:
                continue
            self.family_indices[fam].append(i)
        #print(self.assay_indices, file=sys.stderr)
        #print(ngs_assays.keys(), file=sys.stderr)
        # check if sequenced assays match expected assay
        for f in sorted(self.family_indices):
            #print('family being genotyped',f, file=sys.stderr)
            #matching_reactions = set()
            passing_assays = []
            passing_seqs = []
            for i in sorted(self.family_indices[f]):
                if self.indexed_calls[i] == 'plus':
                    passing_assays.append(self.indexed_primer[i])
                    passing_seqs.append(self.indexed_seqs[i])
            passing_seqs = list(chain(*passing_seqs))
            all_musterer_assays = self.musterer_info['assay_value_options'].keys()
            musterer_assays = [ma for ma in all_musterer_assays if f.lower().replace('-','') in ma.lower().replace('-','')]
            print('Barcode:',self.barcode,'Family:', f, 'Musterer assays:', musterer_assays, passing_assays, passing_seqs, file=sys.stderr)
            if not musterer_assays:  # sometimes just the first part of a name is all that is in common
                musterer_assays = [ma for ma in all_musterer_assays if f.split('-')[0].lower() in ma.split('-')[0].lower()]
            if len(musterer_assays) == 0:
                print('No recognised Musterer assay for family:', f, 'Possible Musterer assays:', all_musterer_assays, file=sys.stderr)
                #self.family_genotype[f] = '?'
                #self.family_to_musterer_assays[f] = ['NA']
                continue
            if len(passing_assays) == 0:
                print('No successful assays found for family:', f, musterer_assays, file=sys.stderr)
                #self.family_genotype[f] = '?'
                #self.family_to_musterer_assays[f] = ['NA']
                continue

            alleles_present = set([ps.split('_')[1] for ps in passing_seqs])
            self.family_to_musterer_assays[f] = musterer_assays
            if len(musterer_assays) == 2: # chose one assay if possible
                if musterer_assays[0].split('_')[0].lower().replace('-','') in musterer_assays[1].split('_')[0].lower().replace('-','') or \
                        musterer_assays[1].split('_')[0].lower().replace('-','') in musterer_assays[0].split('_')[0].lower().replace('-',''):
                    if musterer_assays[0] in self.musterer_info['assay_value_options']:
                        musterer_assays = [musterer_assays[0]]
                    elif musterer_assays[1] in self.musterer_info['assay_value_options']:
                        musterer_assays = [musterer_assays[1]]
                    else:
                        print('Assays not found in Musterer records', musterer_assays, self.musterer_info['assay_value_options'], file=sys.stderr)

            # we can get the allele from the 2nd field of the reference sequence. We should be able to gather these and largely be done with it
            possible_observables = set()
            for ma in musterer_assays:
                try:
                    obs = self.musterer_info['assay_value_options'][ma]
                except KeyError:
                    continue
                for o in obs:
                    possible_observables.add(o)
                self.family_to_musterer_assays[f] = [ma]
                break

            chosen_obs = [o for o in obs if all(alleles_present in o)]
            if len(chosen_obs) == 0:
                print('No genotype found for assay family',f, 'Passing seqs:', passing_seqs, 'Alleles present:', alleles_present, file=sys.stderr)
                self.family_genotype[f] = 'unknown'
            elif len(chosen_obs) == 1:
                self.family_genotype[f] = chosen_obs[0]
            else:
                print('Multiple genotypes possible for assay family', f,'Passing seqs:', passing_seqs, 'Alleles present:', alleles_present, file=sys.stderr)
                self.family_genotype[f] = 'ambig'
            self.musterer_assay_genotype[musterer_assay] = self.family_genotype[f]
            # set all indexed genotypes to be the same for each record in the same assay family
            for i in self.family_indices[f]:
                if f in self.family_genotype: 
                    self.indexed_genotype[i] = self.family_genotype[f]
                else:
                    print('Mismatched assay family name:', f, self.family_genotype.keys(), file=sys.stderr)

            if self.family_genotype[f] not in set(['unknown','ambig']):
                # create upload record
                self.upload_records.append(self.upload_identity + '\t' + musterer_assay + '\t' + self.family_genotype[f])            
            


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

    def write_upload_records_to_CSV(csv_fh):
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



def identify_assay(record, ngs_assays, old_assays, ngs_families, old_families):
    """
        Based on the expected assay family in the primer column of the result table:
            return conversion object(s) expected, or None
    """
    parts = record.primer.split('_')
    expected_assay_family = parts[0]
    additional_info = []
    if len(parts) > 1:
        additional_info = parts[1:]
    #print('NGS assays:', ngs_assays, file=sys.stderr)
    #print('NGS families:', ngs_families, file=sys.stderr)
    if expected_assay_family in ngs_families:
        #print(expected_assay_family)
        #print(ngs_families)
        assays = ngs_families[expected_assay_family]
        #print(assays)
        #sys.exit(0)
        if len(assays) > 1:
            # try to correctly identify the assay
            for ai in additional_info:
                for assay in assays:
                    if ai in assay:
                        return ngs_assays[assay]
            print('Oops, assay is ambiguous from primer info:', record.primer, file=sys.stderr)
            return None
        else:
            return ngs_assays[assays[0]]
    elif expected_assay_family in old_families:
        assays = old_families[expected_assay_family]
        if len(assays) > 1:
            # try to correctly identify the assay
            for ai in additional_info:
                for assay in assays:
                    if ai in assay:
                        return ngs_assays[assay]
            print('Oops, assay is ambiguous from primer info:', record.primer, file=sys.stderr)
            return None
        else:
            return ngs_assays[assays[0]]
    else:
        print('No recognised conversion for assay:', expected_assay_family,
              ' ... Defaulting to Musterer search info', file=sys.stderr)
        return None


def gather_mouse_results(records, ngs_assays, old_assays, ngs_families, old_families):
    """ combine records on a per-mouse basis and return a dictionary of mouse objects """
    mice = {}  # [barcode] = MObs()
    for i,r in enumerate(records):
        #assay_conversion = identify_assay(r, ngs_assays, old_assays, ngs_families, old_families)
        if r.mouseBarcode in mice:
            mice[r.mouseBarcode].add_indexed_record((i,r), old_assays, ngs_assays) #assay_conversion)
        else:
            mice[r.mouseBarcode] = MObs((i,r), old_assays, ngs_assays) #assay_conversion)
    return mice


def output_results_csv(out_fn, records, mice, in_hdr):
    """ output a header and all mouse records to CSV file, matching results.csv input """
    # build new header for output
    out_hdr = ','.join(in_hdr[:23] + ['call','gt','alleleRatios',
            'efficiency','sanity_result', 'sanity_comment','sire_barcode','sire_strain', 'sire_gt',
            'dam_barcode','dam_strain','dam_gt'] + in_hdr[23:])
    #print(out_hdr, file=sys.stderr)
    with open(out_fn, 'wt') as out:
        out.write(out_hdr + '\n')
        for i,r in enumerate(records):
            mbc = r.mouseBarcode
            m = mice[mbc]
            call = m.indexed_calls[i]
            gt = m.indexed_genotype[i]
            allele_ratios = ';'.join(map(str, m.indexed_allele_ratios[i]))
            # TODO: efficiency needs to be based off mergedCount
            efficiency = m.indexed_efficiencies[i]
            #sanity_result = m.indexed_family_match[i]
            sire_bc = m.sire_barcode
            if m.sire_info == 'NA':
                sire_gt = 'NA'
                sire_strain = 'NA'
            else:
                sire_strain = m.sire_info['strain']
                try:
                    sire_gt = m.sire_info['assay_names_values'][m.indexed_primer[i]]
                except:
                    sire_gt = 'NA'
            #dam_bc = m.dam_bc_gt # this needs to be added
            dam_bcs = m.dams_barcodes

            dams_strain = ';'.join(map(str, [m.dams_info[dam_bc]['strain'] for dam_bc in dam_bcs if m.dams_info[dam_bc] != 'NA']))
            try:
                dams_gt = ';'.join(map(str, [m.dams_info[dam_bc]['assay_names_values'][m.indexed_primer[i]] for dam_bc in dam_bcs]))
            except:
                dams_gt = 'NA'

            #print('Record:', i, r, file=sys.stderr)
            outcols = list(r[:23]) + [call, gt, allele_ratios, efficiency, '', '', sire_bc,
                                      sire_strain, sire_gt, ';'.join(map(str,dam_bcs)), dams_strain, dams_gt] + \
                    list(r[23:26]) + list(chain(*r[26]))
            outline = ','.join(map(str, outcols))
            #print(outcols, file=sys.stderr)
            print(outline, file=out)


def output_mouse_gts(mice_gt_fn, mice):
    """ output a CSV of mouse barcodes with genotypes"""
    mice_gt_fn = 'mice_gts.csv'
    pass


def main(result_fn, reference_fn, config_fn, conversions_fn, out_fn=None, mice_fn=None):
    """
        Read reference file
        Read results file...
        For each row in the results table:
            1) extract alleles and counts
            2) then call the alleles
            3) decide the genotype
            4) sanity check the genotype
            5) output genotype table
            6) output per-mouse table
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

    # "call" each well, then call genotypes
    print('"Calling" per-well reactions, and genotyping mice', file=sys.stderr)
    for j, mbc in enumerate(mice):
        m = mice[mbc]
        m.call_reactions()
        m.call_obs()
        print(m, file=sys.stderr)
        print(j+1,'/',len(mice),mbc, file=sys.stderr)
        print('Reactions, assay, and calls:', [(i, m.indexed_reactions[i], m.indexed_primer[i],
                m.indexed_calls[i]) for i in m.indexed_reactions if not m.indexed_empty[i]], file=sys.stderr)
        m.call_genotypes(old_families, ngs_families, old_assays, ngs_assays, config)
        print('GTs by index:', [(i, m.indexed_genotype[i]) for i in m.indexed_records \
                                if not m.indexed_empty[i]], file=sys.stderr)
        #if mbc == '105371500058':
        #    print(m.sire_barcode, m.sire_info)
        #    print(m.dams_info)
        #    sys.exit(0)
        #print(file=sys.stderr)
        #m.check_family_match()

    output_results_csv(out_fn, records, mice, in_hdr)
    #output_mouse_gts('mice_gts.csv', mice)
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
    args = parser.parse_args()
    main(args.resultfile, args.reference, args.config, args.conversions, out_fn=args.outfile)
