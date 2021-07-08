# -*- coding: utf-8 -*-
"""
inheritance_rules.py

@created Nov 2020
@author: Cameron Jack, cameron.jack@anu.edu.au, ANU Bioinformatics Consultancy
@version: 0.8
@version_comment: untested, some cases not implemented
@last_edit:
@edit_comment: 
"""

def test_true_wt(sire_gt, dam_gt):
    """
    Sample result is a TRUE WT:
        Pass = if breeding pair are WT AND WT
        Pass = if breeding trio are WT AND WT AND WT
        Pass = if breeding pair are WT AND HET
        Pass = if breeding trio are WT AND HET AND WT
        Pass = if breeding trio are WT AND WT AND MUT female NOT MUT male
        Pass = if breeding pair are HET AND HET
        Pass = if breeding trio are HET AND HET AND HET
        Pass = if breeding trio are WT AND HET AND HET
        Pass = if breeding trio are HET AND HET AND MUT female NOT MUT male
        Fail = if breeding pair are MUT AND MUT
        Fail = if breeding trio are MUT AND MUT AND MUT
        Fail = if breeding trio are MUT AND MUT AND HET
        Fail = if breeding pair are MUT AND HET
        Fail = if breeding trio are WT AND WT AND MUT male NOT MUT female
        Fail = if breeding trio are HET AND HET AND MUT male NOT MUT female

        Litter result within a cage: the percentage of a WT from a breeding pair WT x WT is 100 %,
        WT x HET is 50 % and HET x HET is 25% of the offspring.
    """
    if sire_gt == 'mut/mut' or dam_gt == 'mut/mut':
        return False
    return True


def test_true_het(sire_gt, dam_gt):
    """
    Sample result is TRUE HET:
        Pass = if breeding pair are WT AND MUT
        Pass = if breeding trio are WT AND WT AND MUT
        Pass = if breeding trio are MUT AND WT AND MUT
        Pass = if breeding pair are MUT AND HET
        Pass = if breeding pair are WT AND HET
        Pass = if breeding trio are WT AND HET AND HET
        Pass = if breeding trio are WT AND HET AND WT
        Pass = if breeding pair are HET AND HET
        Pass = if breeding trio are HET AND HET AND HET
        Pass = if breeding trio are MUT AND HET AND MUT
        Pass = if breeding trio are MUT AND HET AND HET
        Pass = if breeding pair are MUT AND HET
        Fail = if breeding pair are MUT AND MUT
        Fail = if breeding pair are WT AND WT
        Fail = if breeding trio are WT AND WT AND WT
        Fail = if breeding trio are MUT AND MUT AND MUT

    Litter result within a cage: the percentage of a HET from a breeding pair WT x HET is 50%,
        HET x HET is 50 % and HET x MUT is 50% of the offspring.
    """
    if sire_gt == 'mut/mut' and dam_gt == 'mut/mut':
        return False
    elif sire_gt == 'wt/wt' and dam_gt == 'wt/wt':
        return False
    return True


def test_true_mut(sire_gt, dam_gt):
    """
    Sample result is TRUE MUT:
        Pass = if breeding pair are MUT AND MUT
        Pass = if breeding trio are MUT AND MUT AND MUT
        Pass = if breeding trio are MUT AND MUT AND HET
        Pass = if breeding trio are MUT AND HET AND HET
        Pass = if breeding pair are MUT AND HET
        Pass = if breeding pair are HET AND HET
        Pass = if breeding trio are HET AND HET AND HET
        Pass = if breeding trio are WT AND HET AND HET
        Pass = if breeding trio are MUT AND MUT AND WT female NOT WT male
        Fail = if breeding pair are WT AND WT
        Fail = if breeding trio are WT AND WT AND WT
        Fail = if breeding pair are WT AND HET
        Fail = if breeding trio are WT AND HET AND WT
        Fail = if breeding pair are WT AND MUT
        Fail = if breeding trio are WT AND WT AND MUT
        Fail = if breeding trio are MUT AND MUT AND WT male NOT WT female

    Litter result within a cage: the percentage of a MUT from a breeding pair HET x HET is 25%,
        MUT x MUT is 100 % and HET x MUT is 50% of the offspring.
    """
    if sire_gt == 'mut/mut':  # mut sire
        if dam_gt == 'mut/mut' or dam_gt == 'wt/mut' or dam_gt == 'mut/wt':
            return True
        else:
            return False
    elif sire_gt == 'mut/wt' or sire_gt == 'wt/mut':  # het sire
        if dam_gt == 'mut/mut' or dam_gt == 'wt/mut' or dam_gt == 'mut/wt':
            return True
        else:
            return False
    return False  # WT sire


def test_true_pos(sire_gt, dam_gt):
    """
    2. Transgenic or knock-in assays that only detect Positive (POS) or Negative (NEG):
    Sample result is a TRUE POS:
        Pass = if breeding pair are POS AND POS
        Pass = if breeding pair are POS AND NEG
        Pass = if breeding trio are POS AND NEG AND POS
        Pass = if breeding trio are POS AND NEG AND NEG
        Fail = if breeding pair are NEG AND NEG
    """
    if 'NEG' in sire_gt and 'NEG' in dam_gt:
        return False
    return True


def test_true_neg(sire_gt, dam_gt):
    """
    Sample result is a TRUE NEG:
        no check based on parent genotypes possible.
        Surely if either of the parents are HOM then this is a problem!
    """
    if 'pos' in sire_gt and 'pos' in dam_gt:
        return False
    return True


###
# Sex linked assays
###

def test_true_XY_wt(msex, sire_gt, dam_gt):
    """
    3.	X-chromosomal assays:
    In Musterer are easily recognised by the notation for the male mouse using “ /Y”.
        In this assay type the gender of the samples as well as the parents need to be checked
        as this will affect the result and whether they passed.

    This will employ two columns Well “Call” and
    Musterer “Genotype” (or similar) in which:
        1) All reactions will be called as a general call, as “wt/wt”, “wt/mut”
    or “mut/mut” in the “Call” column.
        2) A final “genotype” will be made in a secondary column that will take into
    accumulative account of the gender and make a genotype “wt/Y” or “mut/Y” in a male’s case.

    Samples from female mice:
    If sample is a TRUE WT:
        Pass = if Mother is WT AND Father is WT OR wt/Y
        Pass = if Mother is HET AND Father is WT OR wt/Y
        Fail = if Mother is HET AND Father is MUT OR mut/Y
        Fail = if Mother is MUT AND Father is MUT OR mut/Y
        Fail = if Mother is WT AND Father is MUT OR mut/Y
        Fail = if Mother is MUT AND Father is WT OR wt/Y
    Samples from male mice:
        If sample is a TRUE WT:
        Pass = if Mother is WT AND Father is WT OR wt/Y
        Pass = i f Mother is HET AND Father is WT OR wt/Y
        Pass = if Mother is WT AND Father is MUT OR mut/Y
        Pass = if Mother is HET AND Father is MUT OR mut/Y
        Fail = if Mother is MUT AND Father is MUT OR mut/Y
        Fail = if Mother is MUT AND Father is WT OR wt/Y
    """
    if msex == 'f':
        if dam_gt == 'wt/wt' and (sire_gt == 'wt/wt' or sire_gt == 'wt/y'):
            return True
        if (dam_gt == 'mut/wt' or dam_gt == 'wt/mut') and (sire_gt =='wt' or sire_gt == 'wt/y'):
            return True
        return False
    elif msex == 'm':
        if dam_gt == 'wt/wt' or dam_gt == 'mut/wt' or dam_gt == 'wt/mut':
            if sire_gt == 'wt/wt' or sire_gt == 'wt/y' or sire_gt == 'mut':
                return True
        return False
    return False


def test_true_XY_het(msex, sire_gt, dam_gt):
    """
    Samples from female mice
    If sample is a TRUE HET:
        Pass = if Mother is HET AND Father is WT OR wt/Y
        Pass = if Mother is WT AND Father is MUT OR mut/Y
        Pass = if Mother is HET AND Father is MUT OR mut/Y
        Pass = if Mother is MUT AND Father is WT OR wt/Y
        Fail = if Mother is MUT AND Father is MUT OR mut/Y
        Fail = if Mother is WT AND Father is WT OR wt/Y
    Sample from male mice
Alternatively, Pass if Mother is WT or HET, Fail if Mother is MUT. No check for genotype of father is required.

    """
    if msex == 'f':
        if dam_gt == 'mut/mut':
            return False
        return True
    elif msex == 'm':
        return False


def test_true_XY_mut(msex, sire_gt, dam_gt):
    """
    Samples from female mice:
    If sample is a TRUE MUT:
        Pass = if Mother is MUT AND Father is MUT OR mut/Y
        Pass = if Mother is HET AND Father is MUT OR mut/Y
        Fail = if Mother is MUT AND Father is WT OR wt/Y
        Fail = if Mother is WT AND Father is WT OR wt/Y
        Fail = if Mother is HET AND Father is WT OR wt/Y
        Fail = if Mother is WT AND Father is MUT OR mut/Y
    Samples from male mice:
    If sample is a TRUE MUT OR Mut/Y:
        Pass = if Mother is HET AND Father is WT OR wt/Y
        Pass = if Mother is HET AND Father is MUT OR mut/Y
        Pass = if Mother is MUT AND Father is MUT OR mut/Y
        Pass = if Mother is MUT AND Father is WT OR t/Y
        Fail = if Mother is WT AND Father is WT OR wt/Y
        Fail = if Mother is WT AND Father is MUT OR mut/Y

    Alternatively, Pass if Mother is HET or MUT, Fail if Mother is WT. No check for genotype of father required.
    """
    if msex == 'f':
        if sire_gt == 'wt/wt' or sire_gt == 'wt/y':
            return False
        elif (sire_gt == 'mut/mut' or sire_gt == 'mut/y') and dam_gt == 'wt/wt':
            return False
        return True
    elif msex == 'm':
        if dam_gt == 'wt/wt':
            return False
        return True

"""
NGS genotyping checking for complex assays
1.	Assays for mutations on autosomal chromosomes including transgenic or knock-out assays that
have 1 mutation but have a two-part reaction (results are WT, HET or MUT). In some cases,
especially breeding trios, the gender of the mouse is an important factor in a pass/fail result.

A)	This scenario encapsulates a NGS family group that requires a two-part reaction to get
    one genotyping result. This is the primary strategy for knock-out assays.
•	We have discussed a method on the 14.10.2020 to take a two-part approach to analysis.
This will employ two columns Well “Call” and Musterer “Genotype” (or similar) in which:
1) All reactions will be called as “pos” and “neg” in the “Call” column,
    meaning the allele is present or absent.
2) A final “genotype” will be made in a secondary column that will take into accumulative
    account and make a genotype “wt/wt”, “wt/mut” or “mut/mut”.
"""

def test_cmplx_true_wt(msex, sire_gt, dam_gt):
    """ need to handle trios - hence move to gts for dams

    Sample result is a TRUE WT
    Using the two-part approach to analysis a WT will be called IF the MUT reaction is
    NEG and the WT reaction is POS.

    Pass = if breeding pair are WT AND WT
    Pass = if breeding trio are WT AND WT AND WT
    Pass = if breeding pair are WT AND HET
    Pass = if breeding trio are WT AND HET AND WT
    Pass = if breeding trio are WT AND WT AND MUT female NOT MUT male
    Pass = if breeding pair are HET AND HET
    Pass = if breeding trio are HET AND HET AND HET
    Pass = if breeding trio are WT AND HET AND HET
    Pass = if breeding trio are HET AND HET AND MUT female NOT MUT male
    Fail = if breeding pair are MUT AND MUT
    Fail = if breeding trio are MUT AND MUT AND MUT
    Fail = if breeding trio are MUT AND MUT AND HET
    Fail = if breeding pair are MUT AND HET
    Fail = if breeding trio are WT AND WT AND MUT male NOT MUT female
    Fail = if breeding trio are HET AND HET AND MUT male NOT MUT female
    """
    if 'wt' in sire_gt and 'wt' in dam_gt:
        return True
    else:
        return False


def test_cmplx_true_wt_litter_proportion(sire_gt, dam_gt, proportion_wt, tolerance=0.15):
    """
    Litter result within a cage: the percentage of a WT from a breeding pair WT x WT
    is 100 %, WT x HET is 50 % and HET x HET is 25% of the offspring. Need percentages for trios.
    If the proportion is outside the given tolerance then we flag this.
    """
    pass


def test_cmplx_true_het(sire_gt, dam_gt):
    """
    Sample result is TRUE HET
    Using the two-part approach to analysis a HET will be called IF the MUT reaction
    is POS and the WT reaction is POS.

    Pass = if breeding pair are WT AND MUT
    Pass = if breeding trio are WT AND WT AND MUT
    Pass = if breeding trio are MUT AND WT AND MUT
    Pass = if breeding pair are MUT AND HET
    Pass = if breeding pair are WT AND HET
    Pass = if breeding trio are WT AND HET AND HET
    Pass = if breeding trio are WT AND HET AND WT
    Pass = if breeding pair are HET AND HET
    Pass = if breeding trio are HET AND HET AND HET
    Pass = if breeding trio are MUT AND HET AND MUT
    Pass = if breeding trio are MUT AND HET AND HET
    Pass = if breeding pair are MUT AND HET
    Fail = if breeding pair are MUT AND MUT
    Fail = if breeding pair are WT AND WT
    Fail = if breeding trio are WT AND WT AND WT
    Fail = if breeding trio are MUT AND MUT AND MUT
    """
    if sire_gt != 'wt/wt' and 'wt' in sire_gt and dam_gt != 'wt/wt' and 'wt' in dam_gt:
        return True
    return False


def test_cmplx_true_het_litter_proportion(sire_gt, dam_gt, proportion_het, tolerance=0.15):
    """
    Litter result within a cage: the percentage of a HET from a breeding pair WT x HET is 50%,
    HET x HET is 50 % and HET x MUT is 50% of the offspring.
    """
    pass


def test_cmplx_true_mut(sire_gt, dam_gt):
    """
    Sample result is TRUE MUT
    Using the two-part approach to analysis a MUT will be called IF the MUT reaction is
    POS and the WT reaction is NEG.

    Pass = if breeding pair are MUT AND MUT
    Pass = if breeding trio are MUT AND MUT AND MUT
    Pass = if breeding trio are MUT AND MUT AND HET
    Pass = if breeding trio are MUT AND HET AND HET
    Pass = if breeding pair are MUT AND HET
    Pass = if breeding pair are HET AND HET
    Pass = if breeding trio are HET AND HET AND HET
    Pass = if breeding trio are WT AND HET AND HET
    Pass = if breeding trio are MUT AND MUT AND WT female NOT WT male
    Fail = if breeding pair are WT AND WT
    Fail = if breeding trio are WT AND WT AND WT
    Fail = if breeding pair are WT AND HET
    Fail = if breeding trio are WT AND HET AND WT
    Fail = if breeding pair are WT AND MUT
    Fail = if breeding trio are WT AND WT AND MUT
    Fail = if breeding trio are MUT AND MUT AND WT male NOT WT female
    """
    if 'wt' not in sire_gt and 'wt' not in dam_gt:
        return True
    return False


def test_cmplx_true_mut_litter_proportion(sire_gt, dam_gt, proportion_het, tolerance=0.15):
    """
    Litter result within a cage: the percentage of a MUT from a breeding pair
    HET x HET is 25%, MUT x MUT is 100 % and HET x MUT is 50% of the offspring.
    """
    pass


"""
2.	Transgenic or knock-in assays that only detect Positive (POS) or Negative (NEG)

A)	This scenario is easily recognisable in groups_SB_10.07.2020.csv with the amount
    of references (1 only), for example below

In this scenario, there may have different result inputs, especially in the Musterer
case of CRE, using CREpos and CREneg – but this can be possibly changed when moving
to NGS observables in Musterer.
•	Above in section 1, we have discussed a two-part approach to analysis.
This could be used here. This will employ two columns Well “Call” and
Musterer “Genotype” (or similar) in which:
1) All reactions will be called as a general call, as “pos” or “neg” in the “Call” column.
2) A final “genotype” will be made in a secondary column that will take
    into accumulative account of the gender and make a genotype “CREneg” or “CREpos”.
"""

def test_true_transgene(msex, sire_gt, dam_gts):
    """
    Sample result is a TRUE POS
    Pass = if breeding pair are POS AND POS
    Pass = if breeding pair are POS AND NEG
    Pass = if breeding trio are POS AND NEG AND POS
    Pass = if breeding trio are POS AND NEG AND NEG
    Fail = if breeding pair are NEG AND NEG

    Sample result is a TRUE NEG
    No check based on parent genotypes possible.
    """
    if sire_gt == 'neg' and {gt for gt in dam_gts} == {'neg'}:
        return False
    return True

"""
### 4.	Family groups that encapsulate more than 1 mutant allele (results WT, MUT, HET)

A)	This scenario encapsulates a NGS family group that has more than one Musterer
    observable. It is easily recognisable in groups_SB_10.07.2020.csv with the amount of
    references (more than 2), for example below. IRF4 contains 3 Musterer observables and
    a further 2 other possible mutations in the amplified region.


This will require looking at the reference name for the alignment and go from
here to import the correct results. The naming of the reference can be altered to ensure
the correct result is imported into Musterer.
"""
def test_multi_allele(mgt, sire_gt, dam_gt):
    """
        Multiple (overlapping) primer combos across a region (same family)
        Each separate assays is of form wt/wt,wt/mut,mut/mut
        Do these share primers?
    """
    sire_alleles = set(sire_gt.split('/'))
    dam_alleles = set(dam_gt.split('/'))
    m_alleles = set(mgt.split('/'))
    for ma in m_alleles:
        if ma not in sire_alleles and ma not in dam_alleles:
            return False
    for sa in sire_alleles:
        if sa not in m_alleles and sa in dam_alleles:
            return False
    for da in dam_alleles:
        if da not in m_alleles and da in sire_alleles:
            return False
    return True


"""
### 5. EUCOMM assays

In Musterer have different notation compared to most other assays and are easily recognised using
“tm1c/ ”,“tm1a/ ”, “tm1d/ ”, “tm1b/ ” etc. This section of the complex rules has been
separated as it can be considered an exception to the rule and extremely convoluted
– I would suggest that this type of assays remain as manual analysis for now.

These “Family groups” so far include:
-	Lifr
-	Arpc1b
-	Dd
"""
def test_eucomm_assay(mgt, sire_gt, dam_gt):
    sire_alleles = set(sire_gt.split('/'))
    dam_alleles = set(dam_gt.split('/'))
    m_alleles = set(mgt.split('/'))
    for ma in m_alleles:
        if ma not in sire_alleles and ma not in dam_alleles:
            return False
    for sa in sire_alleles:
        if sa not in m_alleles and sa in dam_alleles:
            return False
    for da in dam_alleles:
        if da not in m_alleles and da in sire_alleles:
            return False
    return True