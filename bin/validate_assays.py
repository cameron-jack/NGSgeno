# -*- coding: utf-8 -*-
"""
@created Dec 2021
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.15
@version_comment:
@last_edit: 2022-01-13
@edit_comment: output in sorted order

This code checks assay data (from the "library") for the NGS Geno application
It then outputs a list of primers used and the number of times each one is needed (to aid 
    in the production of the primer plate)
It also outputs a second list of primers that were not seen in the assay list
Cannot be run as a script. Program should import the Assays class only.
"""

import os
import glob
import csv
import itertools
import sys
import argparse
from collections import Counter

class Assays:
    def __init__(self, assay_fn):
        self.musterer_to_ngs = {}
        self.ngs_to_musterer = {}
        self.all_assays = set()
        self.all_families = set()
        self.valid_assay_counts = Counter()
        self.invalid_assay_counts = Counter()
        with open(assay_fn, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue  # header line
                assay_names = [n.strip() for n in line.strip().split(',')]
                if len(assay_names) != 2:
                    continue  # weird...
                musterer_name, ngs_name = assay_names
                self.musterer_to_ngs[musterer_name] = ngs_name
                self.ngs_to_musterer[ngs_name] = musterer_name
                self.all_assays.add(musterer_name)
                self.all_assays.add(ngs_name)
                self.all_families.add(musterer_name.split('_',1)[0])
                self.all_families.add(ngs_name.split('_',1)[0])

    def convert_assay(self, query_assay):
        """
        Based on the assay list, if the query matches to the Musterer name, return the NGS name, and vice versa
        """
        if query_assay in self.musterer_to_ngs:
            return self.musterer_to_ngs[query_assay]
        elif query_assay in self.ngs_to_musterer:
            return self.ngs_to_musterer[query_assay]
        else:
            return None

    def is_musterer(self, query_assay):
        """ True if a Musterer name """
        if query_assay in self.musterer_to_ngs:
            return True

    def is_ngs(self, query_assay):
        """ True if an NGS name """
        if query_assay in self.ngs_to_musterer:
            return True

    def to_family(self, query_assay):
        """ reduce an assay to it's family name """
        return query_assay.split('_',1)[0]

    def validate_assays(self, plates_data, valid_fn="primerlist_valid.csv", invalid_fn="primerlist_invalid.csv"):
        """
        Match assay families. Use both names for counters and treat them separately, then combine as NGS only.
        Adds values to self.valid_assays and self.invalid_assays
        """
        for n, pid, p_data in plates_data: 
            for well in p_data['wells']:
                for assay in well['mouse']['assays']:
                    if 'assayName' in assay:
                        assay_name = assay['assayName']
                        family_name = self.to_family(assay_name)
                        if family_name in self.all_families:
                            if self.is_musterer(assay_name):
                                self.valid_assay_counts[self.convert_assay(assay_name)] += 1
                            else:
                                self.valid_assay_counts[assay_name] += 1
                        else:
                            if self.is_musterer(assay_name):
                                self.invalid_assay_counts[self.convert_assay(assay_name)] += 1
                            else:
                                self.invalid_assay_counts[assay_name] += 1
        if valid_fn:
            with open(valid_fn, 'wt') as out:
                for vac in sorted(self.valid_assay_counts):
                    print(vac + ',' + str(self.valid_assay_counts[vac]), file=out)
        if invalid_fn:
            with open(invalid_fn, 'wt') as out:
                for invac in sorted(self.invalid_assay_counts):
                    print(invac + ',' + str(self.invalid_assay_counts[invac]), file=out)

        return self.valid_assay_counts, self.invalid_assay_counts
 
def main():
    return

        
if __name__ == '__main__':
    main()
