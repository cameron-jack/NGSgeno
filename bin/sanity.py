#! /usr/bin/env python

from collections import namedtuple
import sys
import csv
import json
import requests
from argparse import ArgumentParser as AP
from itertools import chain

"""

NGS genotyping checking:
Genotyping checklist is set up with basic Boolean operators: AND, NOT & OR
    In some cases the breeding configuration will be a pair or trio, this will need to be taken into account.
    Breeding trios are set up as 2 females to 1 male ratio, conversely breeding pairs are a 1:1 ratio.
Keywords:
	Wildtype = WT OR wt/wt OR wt/Y
	Heterozygous = HET OR wt/mut OR mut/wt
	Mutant = MUT OR mut/mut or mut/Y
	Positive = POS is the presence of transgene
	Negative= NEG is the absence of transgene

1.	Assays for mutations on autosomal chromosomes including transgenic or knock-in assays that are location
    specific (results are WT, HET or MUT):
In some cases, especially breeding trios, the gender of the mouse is an important factor in a pass/fail result.
"""

def simplify_mouse_JSON(json_data):
    """
    Musterer mouse lookup web service returns JSON. Pick out the info we need and return:
    mouse_data = {}  # {mouseBarcode: {sex:M/F, strain:strainName, assays_names_values: {assayName:assayValue},
                     # sire: mouseID, dams:[mouseID]}
    from:
    [
       {
          "barcode":"105254500035",
          "strainId":52545,
          "strainName":"ASD552:Cahill:3::N3F4",
          "mouseId":35,
          "sex":"F",
          "bornDate":"09/04/2020",
          "mouseState":"Alive",
          "dams":[
             {
                "mouseId":14,
                "barcode":105139600013,
                "strainName":"ASD552:Cahill:3::N3F3",
                "strainId":51396
        },
         { 
            "mouseId":14,
            "barcode":105139600018,
            "strainName":"ASD552:Cahill:3::N3F3",
            "strainId":51396
    }

    ],
          "sire":{ 
             "mouseId":14,
             "barcode":105139600014,
             "strainName":"ASD552:Cahill:3::N3F3",
             "strainId":51396
    },
          "assays":[ 
             { 
                "assayMethod":"Amplifluor",
                "assayValueOptions":["wt/wt","mut/wt","mut/mut"],
                "assayId":2134,
                "assayName":"Ikbkb-V203I",
                "assayValue":"mut/wt"
    }
    ]
    }
    ]

    """
    mice = {}
    #print(json_data[0], file=sys.stderr)
    for bunch in json_data:
        for mouse in bunch:
            barcode = str(mouse.get('barcode'))
            sire = mouse.get('sire')
            dams = mouse.get('dams')
            assays = mouse.get('assays')
            m = {'sex' : mouse.get('sex'),
                'strain' : mouse.get('strainName')}
            if sire:
                m['sire_barcode'] = str(sire.get('barcode'))
            if dams:
                m['dams_barcodes'] = [str(d.get('barcode')) for d in dams]
            if assays:
                m['assay_names_values'] = {a.get('assayName'):a.get('assayValue') for a in assays}
                m['assay_value_options'] = {a.get('assayName'):a.get('assayValueOptions') for a in assays}
            mice[barcode] = m
    return mice


def get_musterer_mouse_info(mouse_barcodes):
    """
        Contact the Musterer web service for mouse parental barcodes and Assay/Observables
        inputs: list of mouse barcodes
        intermediate: Send JSON blocks to simplify_mouse_JSON(json_data)
        output: dictionary of simplified mouse data
            mouse_data = {}  # {mouseID: {sex:M/F, strain:strainName, assays: [], sire: mouseID, dams:[mouseID]}
        Notes: No caching of parental information. Not sure how many mice we can ask for info in one transaction.
        https://musterer.apf.edu.au/musterer2/webservice/getMouseGenotypes?barcode=105254500035&barcode=105254500037
        # This code is now broken for security reasons but is obsolete anyway, replaced by gt_mice.py
    """
    barcode_limit = 50 # max chars in GET request is 2048
    barcode_groups = [mouse_barcodes[i:i+barcode_limit] for i in range(0, len(mouse_barcodes), barcode_limit)]
    #print(barcode_groups)
    results = []
    for bg in barcode_groups:
        #print(bg)
        #url = 
        #url += '&barcode='.join(bg)
        #print(url)
        #username=
        #passwd=
        r = requests.get(url, auth=(username, passwd), verify=False)

        assert r.headers['content-type']=="application/json;charset=UTF-8"
        if not r.ok:
            #app.setstatus("Return status =", r.status_code)
            #app.setstatus("when accessing URL:", url)
            return None
        assert r.ok
        if not r.status_code==200:
           #app.setstatus("Request for {service} failed. Return code={status}".format(service=service, status=r.status_code))
           return None
        res = r.json()
        if len(res)==0:
            msg = "Musterer returned no results for barcodes:" + ','.join(bg)
            #app.setstatus(msg)
            #mb.showerror("Musterer error", msg)
        results.append(res)
    return simplify_mouse_JSON(results)


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


def test_true_XY_wt(msex, sire_gt, dam_gt):
    """
    3.	X-chromosomal assays:
    In Musterer are easily recognised by the notation for the male mouse using “ /Y”.
        In this assay type the gender of the samples as well as the parents need to be checked
        as this will affect the result and whether they passed.
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
4.	EUCOMM assays:
In Musterer have different notation compared to most other assays and are easily recognized using “tm1c/ ”, “tm1a/ ”, “tm1d/ ”, “tm1b/ ”  etc. 

-	To be completed later.

"""



def do_tests(mMA, mgt, msex, sire_gt, dam_gt):
    """
        Make sure mMA, mgt, msex make sense given the sire and dams
        Note - only one dam is sent to this function at a time. We are assuming that if
        one valid dam allows a positive test then this is good enough.
    """
    sanity_result = ''
    sanity_comment = ''
    mgtl = mgt.lower()
    sgtl = sire_gt.lower()
    dgtl = dam_gt.lower()
    # Do sex linked assays first
    if '/y' in mMA.lower():
        if mgtl == 'wt/wt' or mgtl == 'wt/y':
            sanity_result = 'Pass' if test_true_XY_wt(msex, sgtl, dgtl) else 'Fail'
        elif mgtl == 'mut/wt' or mgtl =='wt/mut':
            sanity_result = 'Pass' if test_true_XY_het(msex, sgtl, dgtl) else 'Fail'
        elif mgtl == 'mut/mut' or mgtl =='mut/y':
            sanity_result = 'Pass' if test_true_XY_mut(msex, sgtl, dgtl) else 'Fail'
        if sanity_result == 'Fail':
            sanity_comment = 'Sex impossible'
            return sanity_result, sanity_comment
    # now handle regular simplex assays
    elif mgtl == 'mut/mut':
        sanity_result = 'Pass' if test_true_mut(sgtl, dgtl) else 'Fail'
    elif mgtl == 'mut/wt' or mgtl == 'wt/mut':
        sanity_result = 'Pass' if test_true_het(sgtl, dgtl) else 'Fail'
    elif mgtl == 'wt/wt':
        sanity_result = 'Pass' if test_true_wt(sgtl, dgtl) else 'Fail'
    elif 'pos' in mgtl:
        sanity_result = 'Pass' if test_true_pos(sgtl, dgtl) else 'Fail'
    elif 'neg' in mgtl:
        sanity_result = 'Pass' if test_true_neg(sgtl, dgtl) else 'Fail'
    return sanity_result, sanity_comment


def correct_gt(mgt, assay_value_options):
    """
        assay_value_options are Musterer acceptable genotypes. We should correct mgt to match these
    """
    mut_options = {'hom', 'mut', 'pos', 'poscre', '+'}
    wt_options = {'wt', 'neg', 'negcre', '-'}
    avo = [o.lower() for o in assay_value_options]
    if '/' not in mgt:
        return mgt
    try:
        a,b = mgt.split('/')
    except:
        print('Reported mouse genotype is unexpected:', mgt, file=sys.stderr)
        return mgt

    new_gt = ''
    if a in mut_options and b in mut_options: # homozygous mutant
        for o in avo:
            c,d = o.split('/')
            if c in mut_options and d in mut_options:
                return o
        print("Couldn't find", o, 'in', mut_options, wt_options)
        return mgt
    elif (a in mut_options and b in wt_options) or (a in wt_options and b in mut_options): # het
        for o in avo:
            c,d = o.split('/')
            if (c in mut_options and d in wt_options) or (c in wt_options and d in mut_options):
                return o
        print("Couldn't find", o, 'in', mut_options, wt_options)
        return mgt
    elif a in wt_options and b in wt_options: # homozygous wildtype
        for o in avo:
            c,d = o.split('/')
            if c in wt_options and d in wt_options:
                return o
        print("Couldn't find", o, 'in', mut_options, wt_options, file=sys.stderr)
        return mgt
    else:
        print("Couldn't find", mgt, 'in', mut_options, wt_options, file=sys.stderr)
        return mgt


def write_result(out, info, result, comment='', sire_bc='', sire_strain='', sire_gt='', 
        dam_bcs='', dam_strain='', dam_gt='', mgt=None):
    """ write result info to out """
    if mgt:
        original_fields = '","'.join(map(str, list(info[:24]) + [mgt]))
    else:
        original_fields = '","'.join(map(str, list(info[:25])))
    new_fields = '","'.join(map(str, [result, comment, sire_bc, sire_strain, sire_gt, dam_bcs, dam_strain, dam_gt]))
    later_fields ='","'.join(map(str, info[25:28]))
    match_fields = '","'.join(map(str, info[28]))
    outstr = '"' + original_fields + '","' + new_fields + '","' + later_fields + '","' + match_fields + '"'
    outstr = outstr.replace('None', '').replace('none','')
    print(outstr, file=out)


def main():
    """ For each mouse assay, check the parental genotypes to make sure the subject mouse genotype makes sense """
    # Need to have an interface set up to read the genotyping stage output
    # Need graphical interface
    # For each listed mouse assay, look up the parental genotypes
    parser = AP()
    parser.add_argument('results', help='Path to CSV table of genotyping results')
    parser.add_argument('-d', '--debug', action='store_true', help='Display debugging printouts')
    parser.add_argument('-o', '--output', help='Path to output CSV')
    args = parser.parse_args()

    # Read in genotyping results file
    records = []
    with open(args.results, 'rt') as srcfd:
        src = csv.reader(srcfd)
        hdr = next(src)
        hdrstr = ','.join(hdr[:25] + ['Sanity_result', 'Sanity_comment','Sire_BC','Sire_strain',
                 'Sire_GT','Dam_BC','Dam_strain','Dam_GT'] + hdr[25:])
        #print(hdrstr, file=sys.stderr)
        Rec = namedtuple('Rec', hdr[:28]+['matches'])
        #print('fields =', Rec._fields, file=sys.stderr)
        #print(len(Rec._fields), 'fields', file=sys.stderr)
        for rfs in src:
            p1, p2 = rfs[:28], [rfs[27:]]
            if len(p1)<27:
                p1 += [None]*(27-len(p1))
            records.append(Rec(*(p1+p2)))

    # hold all in the info we need in a list of tuples (barcode, assay, genotype, sex)
    mice_bc_assay_gt_sex = []
    for i, r in enumerate(records):
#        if i == 20:
#            break
        mice_bc_assay_gt_sex.append((r.mouseBarcode, r.MustAssay, r.genotype, r.sex))

    mouse_barcodes = [mbc for mbc,mMA,mgt,msex in mice_bc_assay_gt_sex]
    mice_info = get_musterer_mouse_info(mouse_barcodes)
    #if args.debug:
    #    for m in mice_info:
    #        print('Query', m, mice_info[m],'\n', file=sys.stderr)

    ### We need parent information to compare to. We'll link them all up again afterwards
    # sire details for assay info
    sires_barcodes = [mice_info[m]['sire_barcode'] for m in mice_info]
    sires_info = get_musterer_mouse_info(sires_barcodes)

    # Find out the range of possible observables
    #obs = set()
    #for s in sires_info:
    #    if 'assay_value_options' in sires_info[s]:
    #        for a in sires_info[s]['assay_value_options']:
    #            for o in sires_info[s]['assay_value_options'][a]:
    #                obs.add(o)
    #print(obs)

    #if args.debug:
    #    for s in sires_info:
    #        print('Sire', s, sires_info[s],'\n', file=sys.stderr)

    # Now look up dams details for assay info
    dams_barcodes = list(chain(*[mice_info[m]['dams_barcodes'] for m in mice_info]))
    dams_info = get_musterer_mouse_info(dams_barcodes)

    mouse_sanity_result = {}  # [barcode] = ('result','comment')
    # mbc = mouse barcode
    # mgt = mouse genotype
    # mMA = mouse Muster Assay
    num_pass = 0
    num_unsure = 0
    num_fail = 0

    out_fn = 'Results_sanity.csv'
    if args.output:
        out_fn = args.output
    with open(out_fn, 'wt') as out:
        print(hdrstr, file=out)
        for i,r in enumerate(records):
            mbc = r.mouseBarcode
            mMA = r.MustAssay
            msex = r.sex
            mgt = r.genotype
            # deal with weird mouse GTs
            if not mgt or mgt == '' or 'mismatch' in mgt.lower():
                write_result(out, r, 'Fail', 'No GT')
                continue
            if '+' in mgt:  # multiplexed assays
                # work with the first assay only - TODO: make work for multiplexed assays
                mgt = mgt.split('+')[0].strip()
            if ':' in mgt:  # clean away assay info - TODO: make work for multiplex assays
                mgt = mgt.split(':')[1].strip()

            # get parent info
            query_sire_bc = mice_info[mbc]['sire_barcode']
            query_sire_info = sires_info[query_sire_bc]
            query_dams_bcs = mice_info[mbc]['dams_barcodes']
            query_dams_info = [dams_info[qdbc] for qdbc in query_dams_bcs]

            # check that all required information is available
            query_sire_assays_genotypes = query_sire_info.get('assay_names_values')
            query_dams_assays_genotypes = [qdi.get('assay_names_values') for qdi in query_dams_info]
            if args.debug:
                print(i, 'mbc:', mbc, 'mMA:', mMA, 'mgt:', mgt, file=sys.stderr)
                print('>  sire_assays_gts:', query_sire_assays_genotypes, file=sys.stderr)
                print('>  dams_assays_gts:', query_dams_assays_genotypes, file=sys.stderr)
            if not query_sire_assays_genotypes:
                write_result(out, r, 'Unsure', 'No sire assays')
                if r.assay_type == 'simplex':
                    num_unsure += 1
                continue
            if mMA not in query_sire_assays_genotypes:
                write_result(out,r, 'Unsure', 'Sire missing assay')
                if r.assay_type == 'simplex':
                    num_unsure += 1
                continue
            sire_gt = query_sire_assays_genotypes.get(mMA)
            if not sire_gt:
                write_result(out,r, 'Unsure', 'No sire GT')
                if r.assay_type == 'simplex':
                    num_unsure += 1
                continue
            # as long as one dam is okay, we should continue
            result = None
            comment = None
            for dam_assay_gt in query_dams_assays_genotypes:
                if dam_assay_gt:
                    if mMA and dam_assay_gt.get(mMA):
                        dam_gt = dam_assay_gt[mMA]
                        if args.debug:
                            print('mMA:',mMA, 'mgt:', mgt, 'sire_gt:',sire_gt, 'dam_gt:',dam_gt, file=sys.stderr)
                        result,comment = do_tests(mMA, mgt, msex, sire_gt, dam_gt)
                        if args.debug:
                            print(result,comment, file=sys.stderr)
                        if result == 'Pass':
                            break
            # correct GTs to match Musterer allowable values where possible
            assay_value_options = query_sire_info.get('assay_value_options')
            valid_options = None
            if assay_value_options:
                valid_options = assay_value_options.get(mMA)
            if valid_options:
                mgt = correct_gt(mgt, valid_options)

            if not result:
                if r.assay_type == 'simplex':
                    num_unsure += 1
                write_result(out, r, 'Unsure', 'No dam assay', sire_bc=query_sire_bc,
                        sire_strain=';'.join(query_sire_assays_genotypes.keys()), sire_gt=sire_gt, mgt=mgt)
            if result == 'Pass':
                write_result(out, r, 'Pass', '', sire_bc=query_sire_bc,
                        sire_strain=';'.join(query_sire_assays_genotypes.keys()),
                        sire_gt=sire_gt, dam_bcs=';'.join(query_dams_bcs),
                        dam_strain=';'.join(dam_assay_gt.keys()), dam_gt=dam_gt, mgt=mgt)
                if r.assay_type == 'simplex':
                    num_pass += 1
            if result == 'Fail':
                if r.assay_type == 'simplex':
                    num_fail += 1
                write_result(out, r, 'Fail', comment, sire_bc=query_sire_bc,
                        sire_strain=';'.join(query_sire_assays_genotypes.keys()),
                        sire_gt=sire_gt, dam_bcs=';'.join(query_dams_bcs),
                        dam_strain=';'.join(dam_assay_gt.keys()), dam_gt=dam_gt, mgt=mgt)


    print('Passes:',num_pass, 'Unsure:', num_unsure, 'Fails:', num_fail)


if __name__ == '__main__':
    main()
