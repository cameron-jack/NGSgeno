#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: Dec 2020
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: Corrected paths for app folder operation
@last_edit: 2022-05-05
@edit_comment: 

TODO: We need a "lite" version or staged operation which is able to assign standardised genotypes based upon observed count ratios.
    These can then be updated if further info is available to let them match expected naming conventions, etc.

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
from pathlib import PurePath
from configparser import ConfigParser as CP
from argparse import ArgumentParser as AP
from itertools import chain
from openpyxl import load_workbook
import inheritance_rules as IR
import warnings

from bin.db_io import get_musterer_mouse_info
import bin.file_io as file_io
from bin.util import CONFIG_PATH, output_error


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
        try:
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
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.Observable.get_old_obs')


def _flush_conversion_buffer(buffer):
    """ yield one or more new Observable objects from a line buffer """
    try:
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
                    print('Multiple assays detected without new entry',i,row, file=sys.stdout)
        yield Observable(old_family, old_assay, old_obs, ngs_family, ngs_assay, ngs_obs)
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice._flush_conversion_buffer')


def parse_conversions(fn):
    """ ugly Excel file for converting NGS assays to Amplifluor assay names and vice versa """
    try:
        old_families = defaultdict(list)  # dict[family] = list(assays)
        ngs_families = defaultdict(list)  # dict[family] = list(assays)
        old_assays = defaultdict(list)  # dict[assay] = [Observables]
        ngs_assays = defaultdict(list)  # dict[assay] = [Observables]
        buffer = []  # row buffer
        wb = load_workbook(fn, read_only=True)
        for sheet in wb:
            for i,rowobjs in enumerate(sheet.rows):
                row = [str(cell.value).replace('None','').strip() for cell in rowobjs]
                #print(row, file=sys.stdout)
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
        print('Musterer to NGS conversions successfully read.\n', file=sys.stdout)
        print('Musterer family names:', ','.join(old_families), '\n', file=sys.stdout)
        print('Musterer assays:', ','.join(old_assays), '\n', file=sys.stdout)
        print('NGS family names:', ','.join(ngs_families), '\n', file=sys.stdout)
        print('NGS assays:', ','.join(ngs_assays), '\n', file=sys.stdout)
        return old_families, old_assays, ngs_families, ngs_assays
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.parse_conversions')


def get_old_new_assay(query, ngs_families, ngs_assays, old_families, old_assays):
    """ 
    There is a chance of getting a new NGS assay and needing to figure out the matching old assay name.
    Alternatively, get the new ngs_assay from the old assay or family name 
    """
    #print('Assay or family to look up', assay_or_fam, file=sys.stdout)
    #print('ngs_families', ngs_families, file=sys.stdout)
    #print('ngs_assays', ngs_assays, file=sys.stdout)
    #print('old_families', old_families, file=sys.stdout)
    #print('old_assays', old_assays, file=sys.stdout)
    try:
        ngs_query = 'NGS::' + query
        amplifluor_query = 'Amplifluor::' + query
        pcr_query = 'PCR::' + query
        for na in ngs_assays:
            if query.lower().replace('-','') in na.lower().replace('-',''):
                a = ngs_assays[na][0].old_assay
                print('Found assay in NGS assays, returning old assay', a, file=sys.stdout)
                return a
        for oa in old_assays:
            if query.lower().replace('-','') in oa.lower().replace('-',''):
                a = old_assays[oa][0].ngs_assay
                print('Found assay in old assays, returning NGS assay', a, file=sys.stdout)
                return a
        # the code below is basically redundant
        if ngs_query in ngs_assays:
            a = ngs_assays[ngs_query][0].old_assay
            print('Found assay in NGS assays', a, file=sys.stdout)
            return a
        elif ngs_query + query in ngs_families:
            a = ngs_families[ngs_query][0].old_assay
            print('Found assay in NGS families', a, file=sys.stdout)
            return a
        elif amplifluor_query in old_assays:
            a = old_assays[amplifluor_query][0].ngs_assay
            print('Found assay in old assays(Amplifluor)', a, file=sys.stdout)
            return a
        elif pcr_query in old_assays:
            a = old_assay[pcr_query][0].ngs_assay
            print('Found assay in old assays (PCR)', a, file=sys.stdout)
        elif amplifluor_query in old_families:
            a = old_families[amplifluor_query][0].ngs_assay
            print('Found assay in old families (Amplifluor)', a, file=sys.stdout)
            return a
        elif pcr_query in old_families:
            a = old_families[pcr_query][0].ngs_assay
            print('Found assay in old families (PCR)', a, file=sys.stdout)
            return a
        else:
            print('Cannot find', query, 'in known assay or family conversions', file=sys.stdout)
    except Exception:
        output_error(exc, msg='Error in gt_mice.get_old_new_assay')


def read_results(result_fn):
    """ return header and records from Bob's matching of sequences """
    try:
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
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.read_results')


def standardise_allele(config, allele):
    """ take an observed or Musterer allele (half an observable) and try to convert it 
        to standard terms, which are always lower case! """
    try:
        allele = allele.lower()
        #print(f"allele={allele}",file=sys.stdout)
        if allele == 'other':
            return ''
        mut_set = set(config['genotypes']['mut'].lower().split(','))
        wt_set = set(config['genotypes']['wt'].lower().split(','))
        y_set = set(config['genotypes']['y'].lower().split(','))
        pos_set = set(config['genotypes']['pos'].lower().split(','))
        neg_set = set(config['genotypes']['neg'].lower().split(','))
        tm1a_set = set(config['genotypes']['tm1a'].lower().split(','))
        tm1b_set = set(config['genotypes']['tm1b'].lower().split(','))
        tm1c_set = set(config['genotypes']['tm1c'].lower().split(','))
        tm1d_set = set(config['genotypes']['tm1d'].lower().split(','))
        #multi_set = set(config['genotypes']['multi'].lower().split(','))
        unk_set = set(['?','unk'])
        nr_set = set(['nr'])

        if allele in y_set:
            return 'y'
        if allele in mut_set:
            return 'mut'
        elif allele in wt_set:
            return 'wt'
        elif allele in pos_set:
            return 'pos'
        elif allele in neg_set:
            return 'neg'
        elif allele in tm1a_set:
            return 'tm1a'
        elif allele in tm1b_set:
            return 'tm1b'
        elif allele in tm1c_set:
            return 'tm1c'
        elif allele in tm1d_set:
            return 'tm1d'
        #elif allele in multi_set:
        #    return allele
        elif allele in unk_set:
            return '?'
        elif allele in nr_set:
            return 'nr'
        else:
            print('Unknown allele', allele, file=sys.stdout)
            return ''
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.standardise_allele')


def match_assay(my_assay, target_assays, ngs_families, ngs_assays, old_families, old_assays):
    """ try to find matching assays in a target set. Try converting old/new if required """
    try:
        my_assay = my_assay.lower()

        musterer_assays = [ta for ta in target_assays if ta.lower() in my_assay]
        if musterer_assays:
            print('(1)Primers matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
        musterer_assays = [ta for ta in target_assays if ta.lower().replace('-','') in my_assay.replace('-','')]
        if musterer_assays:
            print('(2)Primers matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
        musterer_assays = [ta for ta in target_assays if ta.lower().split('-')[0] in my_assay.split('-')[0]]
        if musterer_assays:
            print('(3)Primers matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
    
        # try at family level
        f = my_assay.split('_')[0].lower()
        musterer_assays = [ta for ta in target_assays if f in ta.lower()]
        if musterer_assays:
            print('(4)Family matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
        musterer_assays = [ta for ta in target_assays if f.replace('-','') in ta.lower().replace('-','')]
        if musterer_assays:  # sometimes just the first part of a name is all that is in common
            print('(5)Family matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
        musterer_assays = [ta for ta in target_assays if f.split('-')[0] in ta.split('-')[0].lower()]
        if musterer_assays:
            print('(6)Family matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays

        # try converting to old/new assay
        musterer_assays = [get_old_new_assay(my_assay, ngs_families, ngs_assays, old_families, old_assays)]
        if musterer_assays:
            print('(7)Assay conversion matching Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
        musterer_assays = [get_old_new_assay(f, ngs_families, ngs_assays, old_families, old_assays)]
        if musterer_assays:
            print('(8)Family conversion Musterer assays:', musterer_assays, file=sys.stdout)
            return musterer_assays
        print('No match found for', my_assay, file=sys.stdout)            
        return []
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.match_assay')


def db_obs_from_std_obs(assay, std_ob):
    """ for a given assay and a standardised observable, return the database-given observable """
    try:
        obs = self.musterer_info['assay_value_options'][assay]
        for ob in obs:
            so1, so2 = [standardise_allele(config, o) for o in ob.split('/')]
            if so1 in std_ob and so2 in std_ob and all([sa in [so1,so2] for sa in std_ob]):
                return ob  # the database observable matching our std_ob
        return ''  # we didn't find a match somehow?!
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.db_obs_from_std_obs')


class MObs():
    """ 
    Mouse Observations, encapsulates per-mouse results. 
    Loads of 'state' is stored here. It needs some cleaning up 
    """
    def __init__(self, idx_rec, old_assays, ngs_assays): #assay_conversion=None):
        try:
            self.barcode = ''
            self.indexed_records = defaultdict(list)  # [idx] = record
            self.family_std_genotype = defaultdict(str)  # dictionary of assay family:(standardised alleles)
            self.family_to_musterer_assays = defaultdict(list)  # for each recorded family, these are the associated assay names
            #self.family_heredity_check = defaultdict(lambda: False)  # does the assay pass heredity sanity checks
            self.musterer_assay_genotype = defaultdict(str) # musterer assay name, associated with an observable as its genotype
            self.musterer_assay_indices = defaultdict(list)  # lets us link back to self.indexed_sanity_results
            self.family_indices = defaultdict(list)
            self.family_sex_linked = defaultdict(lambda: False)  # the assays of this family contain '/Y' alleles
            self.mouseAssays = []  # directly from Rec.mouseAssays
            self.indexed_primer = defaultdict(str)   # [idx] = primer
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
            self.indexed_sanity_result = defaultdict(defaultdict)  # we should always check if this exists
            self.indexed_musterer_assays = defaultdict(list)
            self.indexed_chosen_matches = defaultdict(list)  # contains list of (count,assay) which match the primer
            self.upload_records = []  # created on calling genotypes
            self.checked_upload_records = []  # filtered on heredity checking
            if idx_rec:
                self.barcode = idx_rec[1].mouseBarcode
                self.upload_identity = file_io.unguard_pbc(idx_rec[1].EPplate) + '_' + idx_rec[1].EPwell +\
                        '_' + idx_rec[1].strainName + '#' + idx_rec[1].mouseID
                self.add_indexed_record(idx_rec, old_assays, ngs_assays) #assay_conversion)
                #if not success:
                #    print('Failed to add indexed record', idx_rec, file=sys.stdout)
            #else:
            #    print('Successfully added indexed record', idx_rec, file=sys.stdout)
            # these three are for data we pull from Musterer
            self.musterer_info = defaultdict(None)
            self.sire_barcode = None
            self.sire_info = defaultdict(None)  # [sire_bc] = {dict:sire_info}
            self.sire_gt = None
            self.dams_barcodes = None
            self.dams_info = defaultdict(None)  # [dam_bc] = {dict:dam_info}
            self.dam_bc_gt = None  # barcode of dam that passed sanity check
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.__init__')


    def __str__(self):
        return '\t'.join([self.barcode, 'records:', str(len(self.indexed_records)), 'matches:', str(len(self.indexed_matches))])


    def long_to_short_assay_name(self, assay_long):
        """ helper function used by genotype() and sanity_check(). DEPRECATED? """
        try:
            ngs_assay = assay_long.replace('Amplifluor::','NGS::').replace('PCR::','NGS::')
            amp_assay = ngs_assay.replace('NGS::','Amplifluor::')
            pcr_assay = ngs_assay.replace('NGS::','PCR::')
            if ngs_assay in self.musterer_info['assay_value_options'].keys():
                assay = ngs_assay
                print(f"{assay_long} is NGS", file=sys.stdout)
            elif amp_assay in self.musterer_info['assay_value_options'].keys():
                assay = amp_assay
                print(f"{assay_long} is Amplifluor")
            elif pcr_assay in self.musterer_info['assay_value_options'].keys():
                assay = pcr_assay
                print(f"{assay_long} is PCR")
            else:
                print("This was meant to come from self.musterer_info['assay_value_options']... the code is wrong!", 
                        assay_long, self.musterer_info['assay_value_options'].keys(), file=sys.stdout)
                return ''
            return assay
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.long_to_short_assay_name')


    def add_indexed_record(self, idx_rec, old_assays, ngs_assays, minreads=50): # assay_conversion=None, minreads=100):
        """
            Add a record for the given index and filter assays/counts
            Record the filtered assay matches by index
            Record the allele ratio by index
            Never fails to record something - the final sanity check requires that every record attributable to a
            genotype receives a "Pass"
        """
        try:
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
                print('No matchable sequence reads found for record', idx, file=sys.stdout)
                return

            counts_assays = sorted([(int(c),a) for c,a in matches], reverse=True)
            all_counts = [c for (c,a) in counts_assays]
            if not all_counts or sum(all_counts) < minreads:
                print('Insufficient sequence reads found for record:', idx, 'Minimum required reads:', minreads, file=sys.stdout)
        
            self.indexed_matches[idx] = counts_assays
            self.indexed_primer[idx] = rec.primer
            #return True
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.add_indexed_record')


    def call_genotypes(self, ngs_families, ngs_assays, old_families, old_assays, config):
        """ Establish which assay we are using then compare the sequenced allele tags against the available observables """
        try:
            # gather rows belonging to the same assay family
            for i in self.indexed_primer:
                fam = self.indexed_primer[i].split('_')[0]
                if not fam:
                    continue
                self.family_indices[fam].append(i)
        
            for f in sorted(self.family_indices):
                print('family being genotyped',f, file=sys.stdout)
                all_assays = [p for p in [self.indexed_primer[i] for i in self.family_indices[f]]]
                if len(all_assays) == 0:
                    print('No assay primers found for family:', f, file=sys.stdout)
                    continue 
                all_musterer_assays = self.musterer_info['assay_value_options'].keys()
                # find all primer assays in musterer assay options
                musterer_assays = list(set(chain(*[match_assay(aa, all_musterer_assays, ngs_families, ngs_assays, old_families, old_assays) for aa in all_assays])))
                print('Musterer assays for this mouse:', all_musterer_assays, file=sys.stdout)

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
                            elif c > 0.2*(sum([c for c,a in well_matches])/len(well_matches)):
                                # catch underperforming alleles with 1/3rd the counts of others
                                well_matches.append((c,a))
                                self.indexed_chosen_matches[i].append((c,a))
                    if well_total:
                        efficiency = 100*sum([c for c,a in well_matches])/well_total
                        self.indexed_efficiency[i] = f"{efficiency:.2f}%"
                    family_matches += well_matches
                std_alleles_pres = {standardise_allele(config, a.split('_')[1]):c for c,a in family_matches}
                if not std_alleles_pres:
                    # this could be Cre that's missing
                    print(f"No alleles seen from sequencing. {f}", file=sys.stdout)
                std_alleles_pres = {sa:std_alleles_pres[sa] for sa in std_alleles_pres if sa != ''}
                total_counts_present = sum([c for c,a in family_matches])

                print('Standardised alleles present', std_alleles_pres, file=sys.stdout)

                for assay_long in sorted(musterer_assays):
                    short_assay = assay_long.split('::')[1]
                    assay = ''
                    for ma in self.musterer_info['assay_value_options']:
                        if short_assay.lower() in ma.lower():
                            assay = ma
                    if assay == '':
                        continune

                    obs = self.musterer_info['assay_value_options'][assay]
                    ob_to_std_obs = {}
                    #assay_type = 'simplex'  # complex, Cre and EUCOMM (Tm1x)
                    for ob in obs:
                        # classify our assay
                        #print(f"ob={ob}",file=sys.stdout)
                        #for o in ob.split('/'):
                        #    if 'tm1' in o.lower():
                        #        assay_type = 'EUCOMM'
                        #    elif 'cre' in o.lower():
                        #        assay_type = 'Cre'
                        #    # what does 'complex' look like?
                        ob_to_std_obs[ob] = '/'.join([standardise_allele(config, o) for o in ob.split('/')])
                    print(f"Assay: {assay}, Observables/standardised: {ob_to_std_obs}, Alleles: {std_alleles_pres}", file=sys.stdout)
                
                    # check now if we'll need to worry about sex-linked assays
                    sex_linked_assay = False
                    for ob in ob_to_std_obs:
                        sos = ob_to_std_obs[ob].split('/')
                        if any([so == 'y' for so in sos]):
                            sex_linked_assay = True
                    self.family_sex_linked[f] = sex_linked_assay
                
                    matched = False
                    for ob in ob_to_std_obs:
                        std_ob = ob_to_std_obs[ob]
                        sos = std_ob.split('/')
                        if len(std_alleles_pres) == 0 and len(sos) == 1: # Didn't find Cre!
                            if sos[0] == 'neg':
                                matched = True
                                self.family_std_genotype[f] = std_ob
                                self.family_to_musterer_assays[f].append(assay)
                                self.musterer_assay_genotype[assay] = ob
                                self.musterer_assay_indices[assay].append(i)
                                self.family_allele_ratio[f] = ['1.00']
                                break
                        elif len(std_alleles_pres) == 1 and len(sos) == 1:  # Cre and friends
                            print(f"Assay: {assay}, Observables/standardised: {ob_to_std_obs}, Alleles: {std_alleles_pres}", file=sys.stdout)
                            if sos[0] in std_alleles_pres and all([sa == sos[0] for sa in std_alleles_pres]):
                                matched = True
                                self.family_std_genotype[f] = std_ob
                                self.family_to_musterer_assays[f].append(assay)
                                self.musterer_assay_genotype[assay] = ob
                                self.musterer_assay_indices[assay].append(i)
                                self.family_allele_ratio[f] = ['1.00']
                                break
                            else:
                                continue
                                #print('obs:', ob_to_std_obs, 'alleles:', std_alleles_pres, file=sys.stdout)
                                #exit()
                        elif len(std_alleles_pres) == 1 and len(sos) == 2:
                            so1, so2 = sos
                            if self.musterer_info['sex'] == 'M' and sex_linked_assay:
                                if (so1 == 'y' and so2 in std_alleles_pres) or \
                                        (so2 == 'y' and so1 in std_alleles_pres):
                                    # we matched a sex-linked genotype
                                    matched = True
                                    self.family_std_genotype[f] = std_ob
                                    self.family_to_musterer_assays[f].append(assay)
                                    self.musterer_assay_genotype[assay] = ob
                                    self.musterer_assay_indices[assay].append(i)
                                    self.family_allele_ratio[f] = ['1.00']
                                    break
                            else:  # either female or not sex-linked
                                if so1 in std_alleles_pres and so2 in std_alleles_pres \
                                        and all(sa in (so1,so2) for sa in std_alleles_pres):
                                    # we have an exact match with a homozygous genotype
                                    matched = True
                                    self.family_std_genotype[f] = std_ob
                                    self.family_to_musterer_assays[f].append(assay)
                                    self.musterer_assay_genotype[assay] = ob
                                    self.musterer_assay_indices[assay].append(i)
                                    self.family_allele_ratio[f] = ['1.00']
                                    break
                        elif len(std_alleles_pres) == 2 and len(sos) == 2:
                            # presumably a het
                            so1, so2 = sos
                            if so1 != so2 and so1 in std_alleles_pres and so2 in std_alleles_pres:
                                matched = True
                                self.family_std_genotype[f] = std_ob
                                self.family_to_musterer_assays[f].append(assay)
                                self.musterer_assay_genotype[assay] = ob
                                self.musterer_assay_indices[assay].append(i)
                                allele_ratios = [std_alleles_pres[sa]/total_counts_present for sa in std_alleles_pres]
                                self.family_allele_ratio[f] = [f"{ar:.2f}" for ar in allele_ratios]
                        else:  # 3 or more alleles observed
                            # check for an exact match - might be too restrictive?
                            if all([so in std_alleles_pres for so in sos]) and all([sa in sos for sa in std_alleles_pres]):
                                matched = True
                                self.family_std_genotype[f] = std_ob
                                self.family_to_musterer_assays[f].append(assay)
                                self.musterer_assay_genotype[assay] = ob
                                self.musterer_assay_indices[assay].append(i)
                                allele_ratios = [std_alleles_pres[sa]/total_counts_present for sa in std_alleles_pres]
                                self.family_allele_ratio[f] = [f"{ar:.2f}" for ar in allele_ratios]
                                break
                    if not matched:
                        if f in self.family_std_genotype and self.family_std_genotype[f] != 'Unk':
                            continue  # don't replace a meaningful result!
                        self.family_std_genotype[f] = 'Unk'
                        self.family_to_musterer_assays[f].append(assay)
                        self.musterer_assay_genotype[assay] = 'Unk'
                        self.musterer_assay_indices[assay].append(i)
                        self.family_allele_ratio[f] = ['NA']
                        if len(std_alleles_pres) > 1:
                            print("Sequenced alleles don't match expectation", f, std_alleles_pres, ob_to_std_obs, file=sys.stdout)
                        else:
                            print("Skipping, apparent failed sequencing", f, std_alleles_pres, ob_to_std_obs, file=sys.stdout)
                        continue
            
                # set all indexed genotypes and indexed_musterer_assay to be the same for each record in the same assay family
                for i in self.family_indices[f]:
                    if f in self.family_std_genotype: 
                        #self.indexed_genotype[i] = self.family_std_genotype[f]  # we need this for sanity checking?
                        self.indexed_musterer_assays[i] = self.family_to_musterer_assays[f]
                        assays = [a for a in self.family_to_musterer_assays[f] if a != 'Unk']
                        if len(assays) > 0:
                            assay = assays[0]
                            self.indexed_genotype[i] = self.musterer_assay_genotype[assay]
                        else:
                            self.indexed_genotype[i] = self.family_std_genotype[f]
                        self.indexed_allele_ratio[i] = self.family_allele_ratio[f]
                    else:
                        print('Mismatched assay family name:', f, self.family_std_genotype.keys(), file=sys.stdout) 
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.call_genotypes')


    def do_sanity_tests(self, msex, mgt, sgt, dgt, sex_linked=False):
        """
        msex = query mouse sex
        mgt = query mouse standardised G/T
        sgt = sire standardised G/T
        dgt = dam standardised G/T
        sex_linked [T/F] = True if assay is sex linked
        returns a sanity_result string (Pass/Fail) and a sanity_comment string (usually empty)
        """
        try:
            sanity_comment = ''
            sanity_result = 'Fail'
            # Do sex linked assays first
            print(f"msex={msex} mgt={mgt} sgt={sgt}, dgt={dgt}, sex_linked={sex_linked}", file=sys.stdout)
            msex = msex.lower()
            mgt = mgt.lower()
            sgt = sgt.lower()
            dgt = dgt.lower()
            if mgt == '?':
                sanity_comment = 'Unknown genotype'
                return sanity_result, sanity_comment
            if sex_linked:
                if mgt in {'wt/wt', 'wt/y', 'y/wt'}:
                    sanity_result = 'Pass' if IR.test_true_XY_wt(msex, sgt, dgt) else 'Fail'
                elif mgt in {'mut/wt', 'wt/mut'}:
                    sanity_result = 'Pass' if IR.test_true_XY_het(msex, sgt, dgt) else 'Fail'
                elif mgt in {'mut/mut', 'mut/y', 'y/mut'}:
                    sanity_result = 'Pass' if IR.test_true_XY_mut(msex, sgt, dgt) else 'Fail'
                if sanity_result == 'Fail':
                    sanity_comment = 'Sex impossible'
            # now handle regular simplex assays
            elif mgt == 'mut/mut':
                sanity_result = 'Pass' if IR.test_true_mut(sgt, dgt) else 'Fail'
            elif mgt == 'mut/wt' or mgt == 'wt/mut':
                sanity_result = 'Pass' if IR.test_true_het(sgt, dgt) else 'Fail'
            elif mgt == 'wt/wt':
                sanity_result = 'Pass' if IR.test_true_wt(sgt, dgt) else 'Fail'
            elif 'pos' in mgt:
                sanity_result = 'Pass' if IR.test_true_pos(sgt, dgt) else 'Fail'
            elif 'neg' in mgt:
                sanity_result = 'Pass' if IR.test_true_neg(sgt, dgt) else 'Fail'
            elif 'tm1a' in mgt or 'tm1b' in mgt or 'tm1c' in mgt or 'tm1d' in mgt:
                sanity_result = 'Pass' if IR.test_eucomm_assay(mgt, sgt, dgt) else 'Fail'
            print(f"sanity_result:{sanity_result} sanity_comment:{sanity_comment}", file=sys.stdout)
            return sanity_result, sanity_comment
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.do_sanity_checks')


    def check_family_match(self, ngs_families, ngs_assays, old_families, old_assays, config):
        """ aka Sanity Checking
        IMPORTANT: Any failure to get complete info is a sanity checking fail. We must get 'Pass' marks for 
        every indexed_record.
        we need: mMA, mgt, msex, sire_gt, dam_gt

        Make sure mMA, mgt, msex make sense given the sire and dams
        Only one dam is required to pass for this mouse to pass the checks

        Convert genotypes to standardised form and compare

        Save parental info at this point to self.indexed_sanity_result[i], 
            where i is the index from sle.fsamily_indices
        """
        def get_parental_assay(assays, parental_info, ngs_families, ngs_assays, old_families, old_assays):
            """ return the parental assay appropriate for the query mouse assay - it must have observables """
            try:
                parental_assay = ''
                for a_long in assays:
                    if a_long in parental_info['assay_names_values'] and parental_info['assay_names_values'][a_long] is not None:
                        parental_assay = a_long
                        #print('first',  a_long, parental_info['assay_names_values'][a_long], file=sys.stdout)
                        break
                    elif '::' in a_long:
                        a = a_long.split('::')[1]
                        if a in parental_info['assay_names_values'] and parental_info['assay_names_values'][a] is not None:
                            parental_assay = a
                            #print('second', file=sys.stdout)
                            break
                    else:
                        alt_assay = get_old_new_assay(a_long, ngs_families, ngs_assays, old_families, old_assays)
                        #print('converted:',parental_assay, alt_assay, file=sys.stdout)
                        if not alt_assay:
                            continue
                        if alt_assay in parental_info['assay_names_values'] and parental_info['assay_names_values'][alt_assay] is not None:
                            parental_assay = alt_assay
                            #print('third', file=sys.stdout)
                            break
                        elif '::' in alt_assay:
                            a = alt_assay.split('::')[1]
                            if a in parental_info['assay_names_values'] and parental_info['assay_names_values'][a] is not None:
                                parental_assay = a
                                #print('fourth', file=sys.stdout)
                                break
                return parental_assay
            except Exception as exc:
                output_error(exc, msg='Error in gt_mice.MObs.check_family_match.get_parental_assay')

        try:
            msex = self.musterer_info['sex']
            sanity_dict = {'sanity_result': 'Fail', 'sanity_comment': '',
                    'sire_gt':'Unk', 'sire_strain':'Unk', 'dam_bc': 'Unk', 'dam_gt': 'Unk', 'dam_strain':'Unk'}
            if not self.sire_info:
                sanity_dict['sanity_comment'] = 'No sire information'
            else:
                for fam in self.family_std_genotype:
                    m_std_gt = self.family_std_genotype[fam]
                    sex_linked = self.family_sex_linked[fam]
                    sanity_dict = {'sanity_result': 'Fail', 'sanity_comment': '',
                            'sire_gt':'Unk', 'sire_strain':self.sire_info['strain'], 'dam_bc': 'Unk', 'dam_gt': 'Unk', 'dam_strain':'Unk'}
                    if m_std_gt in {'?', 'Unk'}:
                        sanity_dict['sanity_comment'] = 'Unknown GT'
                        continue
                    if 'assay_names_values' not in self.sire_info:
                        sanity_dict['sanity_comment'] = 'Sire missing assays'
                        continue
                    sire_assays = {a for a in self.sire_info['assay_names_values']}
                    print('Attempting to match:', fam, sire_assays)
                    assays = match_assay(fam, sire_assays, ngs_families, ngs_assays, old_families, old_assays)
                    if len(assays) == 0 or assays[0] is None:
                        sanity_dict['sanity_comment'] = 'Sire assay family not matched:' +fam
                        continue

                    sire_assay = get_parental_assay(assays, self.sire_info, ngs_families, ngs_assays, old_families, old_assays)
                    
                    print(f"sire_assay:{sire_assay}", file=sys.stdout)
                    if not sire_assay:
                        sanity_dict['sanity_comment'] = 'Sire assay has no observables'
                        continue
                    sire_gt = self.sire_info['assay_names_values'][sire_assay]
                    sanity_dict['sire_gt'] = sire_gt
                    print(sire_assay, 'GT:',sire_gt)
                    sire_std_gt = '/'.join([standardise_allele(config, o) for o in sire_gt.split('/')])
                     
                    dam_gt = 'Unk'
                    dam_strain = 'Unk'
                    sanity_result = 'Fail'
                    sanity_comment = 'Dam not found'
                    for dam_bc in self.dams_barcodes:
                        if 'assay_names_values' not in self.dams_info[dam_bc]:
                            sanity_dict['sanity_comment'] = 'Dam missing assays'
                            continue
                        dam_assays = {a for a in self.dams_info[dam_bc]['assay_names_values']}
                        # after getting the g/t, check each dam/sire pairing to see if it works with the mouse in question.
                        # We only need one pair to match for success
                        assays = match_assay(fam, dam_assays, ngs_families, ngs_assays, old_families, old_assays)
                        if len(assays) == 0 or assays[0] is None:
                            sanity_dict['sanity_comment'] = 'Dam assay family not matched:'+fam
                            continue
                        dam_assay = get_parental_assay(assays, self.dams_info[dam_bc], ngs_families, ngs_assays, old_families, old_assays)
                
                        if not dam_assay:
                            sanity_dict['sanity_comment'] = 'Dam assay has no observables'
                            continue
                        dam_gt = self.dams_info[dam_bc]['assay_names_values'][dam_assay]
                        dam_std_gt = '/'.join([standardise_allele(config, o) for o in dam_gt.split('/')])
                        dam_strain = self.dams_info[dam_bc]['strain']
                        sanity_dict['dam_bc'] = dam_bc
                        sanity_dict['dam_gt'] = dam_gt
                        sanity_dict['dam_strain'] = dam_strain

                        print(f"Going to sanity tests with msex:{msex} m_std_gt:{m_std_gt} sire_std_gt:{sire_std_gt} " + 
                                f"dam_std_gt:{dam_std_gt} sex_linked:{sex_linked}", file=sys.stdout)
                        sanity_result, sanity_comment = self.do_sanity_tests(msex, m_std_gt, sire_std_gt, dam_std_gt, sex_linked)
                        sanity_dict['sanity_result'] = sanity_result
                        sanity_dict['sanity_comment'] = sanity_comment
                        print(sanity_dict, file=sys.stdout)
                        if sanity_result == 'Pass':
                            break

                    #sanity_dict = {'sanity_result': sanity_result, 'sanity_comment': sanity_comment,
                    #                'sire_gt':sire_gt, 'sire_strain': sire_strain, 'dam_bc': dam_bc, 'dam_gt': dam_gt, 'dam_strain':dam_strain}

                    for i in self.family_indices[fam]:
                        self.indexed_sanity_result[i] = sanity_dict
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.check_family_match')


    def create_upload_records(self):
        """
        Create records after sanity checking and before writing for efficiency
        """
        try:
            records = {}
            sanity_dict = {'sanity_result': 'Fail', 'sanity_comment': '',
                    'sire_gt':'Unk', 'dam_bc': 'Unk', 'dam_gt': 'Unk', 'dam_strain':'Unk'}
            for assay in sorted(self.musterer_assay_genotype):
                if self.musterer_assay_genotype[assay] not in set(['Unk', 'unknown', 'ambig']):
                    for i in self.musterer_assay_indices[assay]:
                        if i not in self.indexed_sanity_result:
                            self.indexed_sanity_result[i] = sanity_dict
                        elif self.indexed_sanity_result[i]['sanity_result'] == 'Pass':
                            line = '\t'.join([self.upload_identity,assay,self.musterer_assay_genotype[assay]])
                            records[line] = 0  # keep unique records only
            for r in records:
                self.upload_records.append(r)
        except Exception as exc:
            output_error(exc, msg='Error in gt_mice.MObs.create_upload_records')


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
    try:
        mice = defaultdict(None)  # [barcode] = MObs()
        for i,r in enumerate(records):
            if r.mouseBarcode in mice:
                mice[r.mouseBarcode].add_indexed_record((i,r), old_assays, ngs_assays)
            else:
                mice[r.mouseBarcode] = MObs((i,r), old_assays, ngs_assays)
        return mice
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.gather_mouse_results')


def output_results_csv(out_fn, records, mice, in_hdr, ngs_families, ngs_assays, old_families, old_assays):
    """ output a header and all mouse records to CSV file, matching results.csv input """
    try:
        # build new header for output
        out_hdr = ','.join(in_hdr[:23] + ['gt','alleleRatio',
                'efficiency','sanity_result', 'sanity_comment','sire_barcode','sire_strain', 'sire_gt',
                'dam_barcode','dam_strain','dam_gt'] + in_hdr[23:])
        #print(out_hdr, file=sys.stdout)
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
                sire_strain = m.sire_info['strain']

                sanity_result = 'Fail'
                sanity_comment = 'Unknown'
                sire_gt = 'Unk'
                dam_bc = 'Unk'
                dam_gt = 'Unk'
                dam_strain = 'Unk'
                if i in m.indexed_sanity_result:
                    sanity_result = m.indexed_sanity_result[i]['sanity_result']
                    sanity_comment = m.indexed_sanity_result[i]['sanity_comment']
                    sire_gt = m.indexed_sanity_result[i]['sire_gt']
                    dam_bc = m.indexed_sanity_result[i]['dam_bc']
                    dam_gt = m.indexed_sanity_result[i]['dam_gt']
                    dam_strain = m.indexed_sanity_result[i]['dam_strain']

                #print('Record:', i, r, file=sys.stdout)
                outcols = list(r[:23]) + [gt, allele_ratio, efficiency, sanity_result, sanity_comment, sire_bc,
                                          sire_strain, sire_gt, ';'.join(map(str,[dbc for dbc in m.dams_barcodes])), 
                                          dam_strain, dam_gt] + \
                        list(r[23:26]) + list(chain(*r[26]))
                outline = ','.join(map(str, outcols))
                #print(outcols, file=sys.stdout)
                print(outline, file=out)
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.output_results_csv')


def output_mouse_gts(mice_gt_fn, mice):
    """ output a CSV of mouse barcodes with genotypes
    mouseID	EPplate	EPwell	mouseBarcode	strainName	mouseAssays	assayFamilies	dnaplate	dnawell	primer	pcrplate	pcrwell	gt	alleleRatios	efficiency	seqCount	seqName		
    33	71561	A5	105210000027	ASD792:RCE:LoxP::F1rF4N6	RCEloxP	RCE-LoxP	20210706	A3	RCE-LoxP_MUT; RCE-LoxP_WT	3121; 3121	H1; I1	mut/wt			5193; 1715 & 7 	RCE-LoxP_MUT_Neo_100bp; RCE-LoxP_WT_100bp & RCE-LoxP_MUT_Neo_100bp		
    197	71590	B2	105366300045	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	B13	Apc	3121	F5	mut/wt	NA	NA	4128 & 3037	Apc_WT_120bp & Apc_MUT_PM_TtoA_120bp		
    197	71590	B2	105366300045	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	B13	Ku70_MUT; Ku70_WT	3121; 3121	G5; H5 	wt/wt	NA	NA	; 3600	; Ku70_WT_86bp		no reads found for mutant reaction
    193	71590	A2	105366300050	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	A13	Apc	3121	M4	mut/wt	NA	NA	1710 & 1478	Apc_WT_120bp & Apc_MUT_PM_TtoA_120bp		
    193	71590	A2	105366300050	Ku70 ApcMin DKO::::F3	Apc-min;Apc_18:34312602-T>A;Ku70	Apc;Ku70	20210706	A13	Ku70_MUT; Ku70_WT	3121; 3121	N4; O4	Ku70 KO/wt	NA	NA	1176 & 8 ; 1591	Ku70_MUT_Neo_86bp & Ku70_WT_86bp ; Ku70_WT_86bp		
    """
    try:
        #mice_gt_fn = 'mice_gts.csv' - must be provided with correct folder in path
        hdrs = 'mouseID EPplate EPwell mouseBarcode strainName mouseAssays assayFamilies dnaplate '+\
                'dnawell primer pcrplate pcrwell gt alleleRatios efficiency seqCount seqName'
        col_hdrs = hdrs.split(' ')
    
        #print(out_hdr, file=sys.stdout)
        with open(mice_gt_fn, 'wt') as out:
            print(','.join(col_hdrs), file=out)
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
                        info['seqCount'].append('; '.join([str(c) if c != '' else str(0) for c,a in m.indexed_chosen_matches[i]]))
                        info['seqName'].append('; '.join([a if a != '' else 'Empty' for c,a in m.indexed_chosen_matches[i]]))
                    # add placeholders for missing data
                    for field in info:
                        if not info[field]:
                            info[field] = 'NA'
                    info['gt'] = list({gt:0 for gt in info['gt']}.keys())
                    
                    out_line = ','.join([' & '.join(map(str,info[col])) if type(info[col]) is list else str(info[col]) for col in col_hdrs])
                    print(out_line, file=out)
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.output_mouse_gts')


def main(result_fn, config_fn, out_fn=None, mice_fn=None, upload_fn=None):
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
    try:
        warnings.filterwarnings("ignore")
        config = CP()
        config.read(PurePath(config_fn))

        template_files = TemplateFiles()

        # read Amplifluor/NGS assay name conversion document in different indexings
        old_families, old_assays, ngs_families, ngs_assays = parse_conversions(template_files['conversions'])
        in_hdr, records = read_results(result_fn)

        print('Parsing mouse results...', file=sys.stdout)
        mice = gather_mouse_results(records, ngs_assays, old_assays, ngs_families, old_families)
    
        # look up mouse assays from Musterer by mouse barcode and add this to each mouse object
        print('Looking up Musterer records...', file=sys.stdout)
        mouse_barcodes = {r.mouseBarcode for r in records}
        mouse_musterer_info = get_musterer_mouse_info(list(mouse_barcodes))
        if not mouse_musterer_info:
            print('No Musterer info', file=sys.stdout)
        else:
            # collect parent info now, so we can do all checking in one pass
            print('Collecting mouse sire/dam records from Musterer....', file=sys.stdout)
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
        print('Genotyping mice', file=sys.stdout)
        #loggable = []
    
        if upload_fn:
            ur = open(upload_fn, 'wt')
        for j, mbc in enumerate(mice):
            m = mice[mbc]
            print(m, file=sys.stdout)
            print(j+1,'/',len(mice),mbc, file=sys.stdout)
            #print('matches and assay', [(i, m.indexed_matches[i], m.indexed_primer[i]) \
            #        for i in m.indexed_matches], file=sys.stdout)
        
            #debug = m.call_genotypes(ngs_families, ngs_assays, old_families, old_assays, config)
            #loggable.append(debug)
            m.call_genotypes(ngs_families, ngs_assays, old_families, old_assays, config)
            print('GTs by index:', [(i, m.indexed_genotype[i]) for i in m.indexed_records], file=sys.stdout)
            #if mbc == '105371500058':
            #    print(m.sire_barcode, m.sire_info)
            #    print(m.dams_info)
            #    sys.exit(0)
            #print(file=sys.stdout)
            m.check_family_match(ngs_families, ngs_assays, old_families, old_assays, config)
            m.create_upload_records()
            if upload_fn:
                m.write_upload_records_to_CSV(ur)
        if upload_fn:
            ur.close()

        # This is the row-by-row report
        output_results_csv(out_fn, records, mice, in_hdr, ngs_families, ngs_assays, old_families, old_assays)
        # easier to read GT report
        output_mouse_gts('mice_gts.csv', mice)
        print('\nGenotyping complete.', file=sys.stdout)
    except Exception as exc:
        output_error(exc, msg='Error in gt_mice.main')


if __name__ == '__main__':
    parser = AP()
    parser.add_argument('resultfile', help='Path to Results.csv file')
    parser.add_argument('-o', '--outfile', default='results_gt.csv', help='Path to output file. If not present will '+\
                        'output to stdout')
    parser.add_argument('-m', '--mice', help='Path to output per-mouse genotype CSV')
    parser.add_argument('-k','--config',default=CONFIG_PATH,help='Path to config file')
    parser.add_argument('-u','--upload', help='Path to create file for upload of genotypes to Musterer'
                        , default='upload_records.csv')
    args = parser.parse_args()
    main(args.resultfile, args.config, out_fn=args.outfile, mice_fn=args.mice, upload_fn=args.upload)
