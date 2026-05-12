#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 31 Oct 2024
@author: Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Streamlit application to read a genotyping report and output a custom manifest,
and matching Echo_COC_XXX.csv files for retypening failed assays. 

It would be much nicer if it just created DNA plate entries for existing 
DNA plates, but that would also require changing the pipeline, which for 
now is considered frozen. This should be a goal of NGS Genotyping version 2.

        Run with: streamlit --run ngsgeno.py

"""
#from csv import field_size_limit
#from msilib.schema import File
#from xmlrpc.client import boolean
#import jsonpickle
#import os
#import sys
from pathlib import Path  
#from subprocess import check_output, CalledProcessError, STDOUT
#import subprocess
import pandas as pd
import chardet  # For automatic encoding detection
import streamlit as st

import util as util
from stutil import hline, add_css
#import extra_streamlit_components as stx
#from time import sleep

###
# Use Eslam Ibrahim's code below to check for non-ASCII characters
###

######################################################################################################
# This script is used to validate the format of FASTA files.
# It is designed to be used in NGS genotyping pipelines.
# The script utilizes Streamlit to provide a web interface for the end-user.
# To run the script, use the command: #"python -m streamlit run fasta_checker.py"
# This script was developed and tested by Eslam Ibrahim (Eslam.ibrahim@anu.edu.au) on OCT. 10, 2024.
######################################################################################################

def check_non_ascii(file_content):
    """
    A function to check for non-ASCII characters.
    Inputs: character stream
    Outputs: list of dictionaries
    """
    issue_num = 0
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line_num, line in enumerate(lines):
        for idx, char in enumerate(line):
            if ord(char) > 127:  # ASCII characters have values from 0 to 127
                issue_num += 1
                issues.append({
                    'Issue Number': issue_num,
                    'Line Number': line_num,
                    'Issue': f"Invalid character in sequence header: '{char}' at position {idx + 1}"
                })
        line_num += 1  # Move to the next line number
    return issues  


def check_valid_sequence(file_content):
    """
    Function to check that sequences contain only A, T, C, G, N
    Spaces and tabs will not be reported here, only checked in check_gaps
    Inputs: character stream
    Outputs: list of dictionaries
    """
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line in lines:
        if line.startswith(">") or not line.strip():  # Ignore headers and blank lines
            line_num += 1
            continue
        
        sequence = line.strip().replace(" ", "").replace("\t", "")  # Clean the sequence

        for idx, char in enumerate(sequence):
            if char not in "AaTtCcGgNn":  # Check for valid characters
                global issue_num
                issue_num += 1
                issues.append({
                    'Issue Number': issue_num,
                    'Line Number': line_num,
                    'Issue': f"Invalid character in sequence: '{char}' at position {idx + 1}"
                })
        
        line_num += 1  # Move to the next line number

    return issues
#=====================================
# 
def check_gaps(file_content):
    """
    Function to check for gaps (spaces or tabs) within sequences.
    Inputs: character stream
    Outputs: list of dictionaries
    """
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line in lines:
        if line.startswith(">") or not line.strip():  # Ignore headers and blank lines
            line_num += 1
            continue

        if " " in line or "\t" in line:
            global issue_num
            issue_num += 1
            issues.append({
                'Issue Number': issue_num,
                'Line Number': line_num,
                'Issue': "Gap (space or tab) in sequence"
            })
        
        line_num += 1  # Move to the next line number

    return issues


def check_blank_lines(file_content):
    """
    FASTA files should not contain blank lines
    Inputs: character stream
    Outputs: list of dictionaries
    """
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line in lines:
        if not line.strip():  # Blank line
            global issue_num
            issue_num += 1
            issues.append({
                'Issue Number': issue_num,
                'Line Number': line_num,
                'Issue': "Blank line found"
            })
        line_num += 1  # Move to the next line number
    return issues


def check_fasta_file(file_content):
    """
    Issue feature checks for FASTA specific file features
    Inputs: character stream
    Outputs: list of dictionaries
    """
    issues = []
    issues.extend(check_non_ascii(file_content))
    issues.extend(check_valid_sequence(file_content))
    issues.extend(check_gaps(file_content))
    issues.extend(check_blank_lines(file_content))
    return issues   


def parse_failed_genotyping_results(file_content):
    """
    Genotyping returns a tab-delimited report file with the results of the genotyping assessment for each assay.
    The header looks like this:
        barcode	code_assays	plate	wellLocation	sex	alleleSymbol	alleleKey	assayKey\
        passFail	seqName1	seqName2	seqNameNNN  efficiency	alleleRatio\
        alleleRatioAdjusted	genotype	args	reason
    Presumably we are working with passFail?
    
    We need code_assays (sample barcode plus assay), plate, wellLocation, passFail
    """
    #failed_assays = {}  # (barcode,plate,wellLocation,sex) = [failed_assays]
    try:
        seqName_cols = []
        genotype_col = []
        for line in file_content:
            if line.startswith('barcode'):
                # look for seqName, which is an assay column header
                cols = [c.strip() for c in line.split(',')]
                for i,c in enumerate(cols):
                    if 'seqName' in c:
                        seqName_cols.append(i)
                    if 'genotype' in c.lower():
                        genotype_col.append(i)
                continue  # header
            #print(f'{seqName_cols=} {genotype_col=}')
            cols = [c.strip() for c in line.split(',')]
            #print(f'{cols=}')
            if all([not c for c in cols]):
                continue
            code_assays = cols[1]
            barcode, assay = code_assays.split(';')
            plate = cols[2]
            wellLocation = cols[3]
            sex = cols[4]
            alleleSymbol = cols[5].strip()
            alleleKey = cols[6].strip()
            assayKey = cols[7].strip()
            passFail = cols[8].strip()
            genotype = cols[genotype_col[0]]
            seqNames = [cols[i] for i in seqName_cols]
            reason = cols[-1]
            ident = (barcode, plate, wellLocation, sex)
            #print(f'{ident=} {passFail=} {genotype=}')
            if '?' in genotype:
                if ident not in st.session_state['failed_assays']:
                    st.session_state['failed_assays'][ident] = []
                st.session_state['failed_assays'][ident].append((assay,alleleSymbol))
    except Exception as exc:
        st.write(f'Failed to parse results workbook {exc}')


# def generate_retype_manifest(failed_assays, manifest_name='../Downloads/retype_failures.csv'):
#     """
#     We need to output a CSV (comma delimited) with the following headers:
#         Sample no	plateBarcode	well	sampleBarcode	Assay	Assay	Assay	Assay	\
#         Assay	Assay	Assay	clientName	sampleName 	alleleSymbol
        
#     Sample no. is automatically incremented from 1 for each sample.
    
#     inputs: failed_assays {(barcode,plate,well,sex)=[(assay,alleleSymbol)]}
#     """
#     with open(manifest_name, 'wt') as f:
#         # header
#         print(','.join(['Sample no','plateBarcode','well','sampleBarcode','Assay','Assay',
#                 'Assay','Assay','Assay','Assay','Assay','clientName','sampleName','alleleSymbol']), file=f)
        
#         for i, ident in enumerate(failed_assays):
#             barcode, plate, wellLocation, sex = ident
#             all_assays = []
#             all_alleleSymbols = []
#             for assay, alleleSymbol in failed_assays[ident]:
#                 all_assays.append(assay)
#                 all_alleleSymbols.append(alleleSymbol)
#             row_array = [str(i+1),plate, wellLocation, barcode] + all_assays + ['']*(7-len(all_assays)) + ['retype','',';'.join(all_alleleSymbols)]
#             print(','.join(row_array), file=f)
#     return manifest_name


def generate_echo_coc_files():
    """
    Not wanted by the genotyping team.
    """
    pass


def parse_stage3_csv(file_content):
    """
    Read csv file stream and return a dictionary

    sampleNo	samplePlate	sampleWell	sampleBarcode	strain	sex	alleleSymbol	alleleKey	assayKey	assays	assayFamilies	clientName	sampleName	dnaPlate	dnaWell	primer	primerPlate	primerWell	pcrPlate	pcrWell	i7bc	i7name	i7well	i5bc	i5name	i5well	index_plate
    1	p230123VRP1p	A1	cM014502400050c						Dst-ex75	Dst-ex75	Rodentity_Validation	M014502400050	p2306D04p	A1	Dst-ex75	pP1p	B12	p2306R2PCR1p	A1	AACCTTGG	96_SET_1_i7F_1	A1	AACCGAAC	96_SET_1_i5R_1	A7	Index

    return a dictionary[(barcode,assay)] = [(dnaPlate,dnaWell,alleleSymbol,alleleKey,assayKey,assay,assayFamily,clientName,sampleName,)]
    """
    counter = 0
    try:
        for i, line in enumerate(file_content):
            l = line.replace('"','')
            if i == 0:
                if l.startswith('sampleNo'):
                    continue  # expected header
                else:
                    print(f'Failed to detect header in file')
            cols = l.strip().split(',')
            if len(cols)<5 or cols[0] == '':
                continue
            barcode = util.unguard(cols[3], silent=True)
            assay = cols[9]
            assayFamily = cols[10]
            dnaPlate = cols[13]
            dnaWell = cols[14]
            alleleSymbol = cols[6]
            alleleKey = cols[7]
            assayKey = cols[8]
            primer = cols[15]
            clientName = cols[11]
            sampleName = cols[12]
            if barcode not in st.session_state['stage3_dict']:
                st.session_state['stage3_dict'][barcode] = []
            st.session_state['stage3_dict'][barcode].append((dnaPlate,dnaWell,alleleSymbol,alleleKey,assayKey,assay,assayFamily,clientName,sampleName))
            counter += 1
    except Exception as exc:
        st.write(f'Parsing Stage3.csv content failed {exc}')
    st.write(f'Read {counter} entries')


def collate_manifest_entries(failed_assays, stage3_dict):
    """
    Helper function that matches failed assays to DNA plate sample entries
    Also returns the maximum number of assays across all samples
    """
    # collect entries by plate first
    # {(barcode,samplePlate,sampleWell,sex)=[(assay,alleleSymbol)]} from failed_assays
    plate_set = {}  # [PID] = {well:(sample, assay...., clientName,sampleName,alleleSymbol...)}
    #print(f'{failed_assays=}')
    #print(f'{stage3_dict=}')
    max_assays = 0
    try:
        for ident in failed_assays:
            barcode,_,_,_ = ident
            assays = [a[0] for a in failed_assays[ident]]
            if len(assays) > max_assays:
                max_assays = len(assays)
            alleleSymbols = [a[1] for a in failed_assays[ident]]
            sampleName = ''
            clientName = ''
            dnaPlate = ''
            dnaWell = ''
            if barcode in stage3_dict:
                for result in stage3_dict[barcode]:
                    dnaPlate_,dnaWell_,_,_,_,_,_,clientName_,sampleName_ = result
                    if clientName_ and not clientName:
                        clientName = clientName_
                    if sampleName_ and not sampleName:
                        sampleName = sampleName_
                    if dnaPlate_ and not dnaPlate:
                        dnaPlate = dnaPlate_
                    if dnaWell_ and not dnaWell:
                        dnaWell = dnaWell_
                    if dnaPlate and dnaWell:
                        if dnaPlate not in plate_set:
                            plate_set[dnaPlate] = {}
                        if dnaWell not in plate_set[dnaPlate]:          
                            plate_set[dnaPlate][dnaWell] = {'barcode':barcode,'assays':assays,'alleleSymbols':alleleSymbols,'clientName':'','sampleName':''}
                        if clientName and not plate_set[dnaPlate][dnaWell]['clientName']:
                            plate_set[dnaPlate][dnaWell]['clientName'] = clientName
                        if sampleName and not plate_set[dnaPlate][dnaWell]['sampleName']:
                            plate_set[dnaPlate][dnaWell]['sampleName'] = sampleName
    except Exception as exc:
        st.write(f'Failed to collate records {exc}')
    return plate_set, max_assays


def generate_manifest_384(failed_assays, stage3_dict, manifest_name='../Downloads/retype_manifest_384.csv'):
    """
    Much like a custom manifest but 384 well and potentially more metadata
    Use information from the Stage3.csv file to get the existing DNA plate locations.
    Inputs: Stage3.csv dict, failed assay dict[(barcode,plate,well,sex)] = [assays]

    We need to output a CSV (comma delimited) with the following headers:
        Sample no	plateBarcode	well	sampleBarcode	Assay...	\
                clientName	sampleName 	alleleSymbol
    There should be as many assay columns as needed
        
    Sample no. is automatically incremented from 1 for each sample.
    """
    success = False
    plate_set, max_assays = collate_manifest_entries(failed_assays, stage3_dict)
    print(f'{plate_set=} {max_assays=}')
    if not plate_set:
        st.write('Failed to collate any entries across the two files. Are they from the same experiment?')
        return success
    sample_counter = 0
    try:
        with open(manifest_name, 'wt') as fout:
            # print the header 
            print(','.join(['Sample no','plateBarcode','well','sampleBarcode'] + ['Assay']*max_assays +\
                    ['clientName','sampleName','alleleSymbol']), file=fout)
            for gpid in sorted(plate_set):
                for pos in util.col_ordered_384:
                    if pos in plate_set[gpid]:
                        barcode = plate_set[gpid][pos]['barcode']
                        pid = util.unguard_pbc(gpid, silent=True)
                        assays = plate_set[gpid][pos]['assays']
                        padding_cols = ['']*(max_assays-len(plate_set[gpid][pos]['assays']))
                        if padding_cols:
                            assays.extend(padding_cols)
                        clientName = plate_set[gpid][pos]['clientName']
                        sampleName = util.unguard(plate_set[gpid][pos]['sampleName'], silent=True)
                        alleleSymbol = ';'.join(plate_set[gpid][pos]['alleleSymbols'])
                        sample_counter += 1
                        print(f','.join([str(sample_counter),pid,pos,barcode]+assays+[clientName,sampleName,alleleSymbol]), file=fout)
            st.write(f'Wrote {sample_counter} 384-well entries')
            success = True 
    except Exception as exc:
        st.write(f'Failed to write to {manifest_name}. Do you have an older version of the file open for viewing? {exc}')
    return success


def display_file_character_issues(issues):
    """
    Common GUI widget for displaying any file check problems
    """
    # Create DataFrame and sort by Line Number if necessary
    df = pd.DataFrame(issues)
    # df = df.sort_values(by="Line Number")  # Uncomment if you want to sort by line number
        
    # Display the number of issues found and the table with custom size
    st.write(f"There are {len(issues)} issues found in the file:")
    st.dataframe(df, height=300, width=500)


def read_text_file(file_stream):
    """
    Common text file handling to search for character issues before parsing
    """
    # Read the raw file
    rawfile = file_stream.read()
    
    # Detect encoding and Decode the file content
    result = chardet.detect(rawfile)
    if not result:
        st.write('Could not interpret file format')
        return None
    charenc = result['encoding']
    if not charenc:
        st.write('Could not interpret file format')
        return None

    file_content = rawfile.decode(charenc)
    
    # Check for issues in the file content
    issues = check_non_ascii(file_content)
    if issues:
        display_file_character_issues(issues)
    else:
        return ''.join(file_content).split('\n')
    return None


def get_report_cb():
    """
    Callback to read the genotyping report file
    """ 
    if st.session_state['gt_report_uploader'] is not None:
        file_content = read_text_file(st.session_state['gt_report_uploader'])
        if file_content:
            parse_failed_genotyping_results(file_content)
   

def get_stage3_cb():
    """
    Callback to read the Stage3.csv file
    """
    stage3f = read_text_file(st.session_state['gt_stage3_uploader'])
    if stage3f:
        parse_stage3_csv(stage3f)


#=====================================
# Streamlit UI

def main():
    """
    The NGSgeno "Xplorer" application. Allows full control of all sections of the pipeline,
    and displays all aspects of the experiment state at any time.
    """    
    st.set_page_config(
        page_title="NGSG: Retype failed assays",
        page_icon="ngsg_icon.png",
        layout="wide"
    )
    add_css()
    st.title("NGS Genotyping: Retype failed assays")
    st.subheader("Upload a genotyping report file (results workbook) in the left hand panel, and the matching Stage3.csv for this run in the right hand panel")
    hline()
    if 'failed_assays' not in st.session_state:
        st.session_state['failed_assays'] = {}
    if 'stage3_dict' not in st.session_state:
        st.session_state['stage3_dict'] = {}
    screen_col1, screen_col2 = st.columns(2)
    
    with screen_col1:
        with st.container(border=True):
            #st.write("Upload a genotyping report file and a custom "+\
            #        "manifest will be generated to retype failed assays")
            st.file_uploader("Choose a genotyping report file", on_change=get_report_cb, key='gt_report_uploader')
            if st.session_state['failed_assays']:
                st.write('The following sample/assay combinations appear to have failed from the run:')
                scrollable = st.container(height=400, border=True)
                for ident in st.session_state['failed_assays']:
                    barcode, plate, wellLocation, sex = ident
                    assays = '; '.join([a.strip() for a,ak in st.session_state['failed_assays'][ident] if a.strip() != ''])
                    allele_keys = '; '.join([ak.strip() for a,ak in st.session_state['failed_assays'][ident] if ak.strip() != ''])
                    row = ' '.join(['Barcode:', barcode, '  plate:', plate, '  well:', wellLocation, '  assays:', assays, '  allele_keys:', allele_keys])
                    scrollable.text(row)


    with screen_col2:
        with st.container(border=True):
            if st.session_state['failed_assays']:
                st.write("Upload matching Stage3.csv file here")
                st.file_uploader("Select the matching Stage3.csv file", on_change=get_stage3_cb, key='gt_stage3_uploader')
            else:
                st.write("No failed assays seen")
                
    hline()
    if st.session_state['failed_assays'] and st.session_state['stage3_dict']:
        manifest_name='../Downloads/retype_manifest_384.csv'
        success = generate_manifest_384(st.session_state['failed_assays'], 
                st.session_state['stage3_dict'], manifest_name=manifest_name)
        if success:
            st.write(f'Successfully wrote new 384-well plate manifest to: {manifest_name}')
        else:
            st.write(f'Attempt failed')
        st.session_state['failed_assays'] = None
        st.session_state['stage3_dict'] = None

    refresh = st.button('Clear session')
    if refresh:
        st.session_state['failed_assays'] = {}
        st.session_state['stage3_dict'] = {}
        st.rerun()
    

if __name__ == '__main__':
    main()