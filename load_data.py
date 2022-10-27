#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New streamlit interactive web interface
@last_edit: 2022-08-29
@edit_comment:

A web based interactive GUI with Streamlit. Plate and sample barcodes here are unguarded.
In all other code they must be guarded. We guard them here before we send them to any external function.
"""

from ctypes.wintypes import WIN32_FIND_DATAA
from msilib.schema import File
import os
import sys
from pathlib import PurePath, Path
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO
from shutil import copy2

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components

from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode
#try:
#    from bin.experiment import Experiment, EXP_FN, load_experiment
#except ModuleNotFoundError:
#    from experiment import Experiment, EXP_FN, load_experiment
try:
    import bin.util as util
except ModuleNotFoundError:
    import util
try:
    import bin.file_io as file_io
except ModuleNotFoundError:
    import file_io
try:
    import bin.db_io as db_io
except ModuleNotFoundError:
    import db_io
try:
    from bin.makehtml import generate_heatmap_html
except ModuleNotFoundError:
    from makehtml import generate_heatmap_html
import display_components as dc


credits="""
@created: March 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: New interactive web application interface
@last_edit: 
@edit_comment: 

Web application style interface using Streamlit. Present information as both a dashboard and a workflow in one page.

The GUI interacts with a single Experiment object at one time. Methods are called on this to activate pipeline
functionality. The Experiment then deals directly with the pipeline logic.
"""


def load_rodentity_data():
    """
    Home page 
    """
    rodentity_exp = st.expander('Add data from Rodentity JSON files',expanded=True)
    rod_col1, rod_col2= rodentity_exp.columns(2)

    rod_dp = ''
    plates_to_clear = [False, False, False, False]
    rod_up_disabled = False  # disable uploading more plates until space available
    exp = st.session_state['experiment']
    #print(f"{exp.unassigned_plates=}")
    #if 'prev_rod_upload' in st.session_state:
    #    print(f"{st.session_state['prev_rod_upload']=}")
    #if 'rod_dp_key' in st.session_state:
    #    print(f"{st.session_state['rod_dp_key']=}")
    if all(exp.unassigned_plates[k] for k in exp.unassigned_plates):
        rod_up_disabled = True              

    rodentity_epp = rodentity_exp.file_uploader('Rodentity JSON file', type='json', disabled=rod_up_disabled)
    if rodentity_epp:
        rodentity_plate_name = util.guard_pbc(rodentity_epp.name.rstrip('.json'))
        if 'prev_rod_upload' not in st.session_state or \
            st.session_state['prev_rod_upload'] != rodentity_plate_name:
            if rodentity_epp.name.rstrip('.json') in exp.unassigned_plates:
                print('We see this already')
            else: 
                for k in sorted(exp.unassigned_plates):
                    if not exp.unassigned_plates[k]:
                        exp.unassigned_plates[k] = rodentity_plate_name
                        # Save a copy we can edit and reload
                        fn = os.path.join('run_'+exp.name, db_io._eppfn_r(rodentity_plate_name))
                        with open(fn, 'wt') as outf:
                            outf.write(rodentity_epp.getvalue().decode("utf-8"))
                        print('here', k)
                        st.session_state['prev_rod_upload'] = rodentity_plate_name
                        st.experimental_rerun()     


    rod_col1.write('Checkboxes required for clearing plate IDs')
                    
    for i,p in enumerate(plates_to_clear):
        plates_to_clear[i] = rod_col1.checkbox(f"P{str(i+1)}: {util.unguard_pbc(exp.unassigned_plates[i+1], silent=True)}", 
                help='Click the checkbox to allow a set plate ID to be cleared', key='check'+str(i+1))

    #Accept rodentity button
    accept_disabled = False
    rod_dp = rod_col2.text_input('Destination plate barcode', max_chars=20, key='rod_dp_key') 

    if rod_dp and any(exp.unassigned_plates):
        rod_dp = util.guard_pbc(rod_dp, silent=True)
        if rod_dp in st.session_state['experiment'].dest_sample_plates or \
                rod_dp in st.session_state['experiment'].plate_location_sample or \
                rod_dp in st.session_state['experiment'].unassigned_plates:
            rod_col1.markdown('<p style="color:#FF0000">Destination plate barcode already in use: ' +\
                    util.unguard_pbc(rod_dp) + '</p>', unsafe_allow_html=True)
    else:
        accept_disabled = True

    accept_rod_button = rod_col2.button('Accept', help='Read and confirm plates', 
                    disabled=accept_disabled, key='accept_rod_button_key')

    if accept_rod_button:
        success = exp.add_rodentity_plate_set([exp.unassigned_plates[k] for k in exp.unassigned_plates], rod_dp)
        if not success:
            rod_col2.markdown('<p style="color:#FF0000">Failed to incorporate plate set. Please read the log.</p>', unsafe_allow_html=True)
        else:
            rod_col2.write('Successfully added plate set')
            rod_dp = ''
            exp.unassigned_plates = {1:'',2:'',3:'',4:''}
            plates_to_clear = [False, False, False, False]
            accept_rod_button = False
            exp.save()
            #?
            st.experimental_rerun()

    #Clear ID button
    clear_disabled = True
    for i,plate in enumerate(plates_to_clear):
        if plate and exp.unassigned_plates[i+1]:
            clear_disabled = False
            break  # only need to enable the button once

    clear_plates_button = rod_col2.button('Clear IDs', help='Clear selected Rodentity plate IDs', disabled=clear_disabled)

    cpb_activated = False
    if clear_plates_button:
        for i, plate in enumerate(plates_to_clear):
            if plate and exp.unassigned_plates[i+1]:
                if 'prev_rod_upload' in st.session_state:
                    if st.session_state['prev_rod_upload'] == plate:
                        st.session_state['prev_rod_upload'] = ''
                exp.unassigned_plates[i+1] = ''
                plates_to_clear[i] = False
                cpb_activated = True
        clear_plates_button = False
        if cpb_activated:
            exp.save()
            st.experimental_rerun()
                
             
def load_custom_csv(expanded=False):
    '''
    Upload data from a custom csv file
    '''
    custom_csv_exp = st.expander('Add data from a custom manifest CSV file', expanded=expanded)
    custom_col1,custom_col2 = custom_csv_exp.columns(2)
    
    
    st.session_state['default_manifest_type'] = 'c'
    manifest_choice = custom_col1.radio('The default contents of my manifest file (if not declared) are', 
                        ['Custom', 'Rodentity mouse', 'Musterer mouse'])

    uploaded_file = custom_col2.file_uploader('', type='csv')

    if manifest_choice == 'Rodentity mouse':
        st.session_state['default_manifest_type'] = 'r'
    elif manifest_choice == 'Musterer mouse':
        st.session_state['default_manifest_type'] = 'm'

    if uploaded_file:
        if 'custom_upload' not in st.session_state or st.session_state['custom_upload'] != uploaded_file.name:
            success = st.session_state['experiment'].add_manifest(uploaded_file, st.session_state['default_manifest_type'])
            if not success:
                custom_col2.markdown('<p style="color:#FF0000">Failed to incorporate manifest. Please read the log.</p>', unsafe_allow_html=True)
            st.session_state['custom_upload'] = uploaded_file.name
            #?
            st.experimental_rerun()


def load_database_data():
    database_exp = st.expander('Add data from a database')
    #??
    db_r1cols = database_exp.columns(4)
    db_r2cols = database_exp.columns(4)

    epp1 = epp2 = epp3 = epp4 = dnap = ''
    epp1_str = epp2_str = epp3_str = epp4_str = None
    epp_files = [epp1, epp2, epp3, epp4]
    epp_strings = [epp1_str, epp2_str, epp3_str, epp4_str]
    dnap_msg = 'DNA plate barcode'

    #Create buttons for the four Ear punch plate and file (row 1 columns) - not sure if it works yet
    for i in range(4):
        num=str(i+1)
        epp_msg = f"Ear punch plate {num} barcode"
        key_txt = f"epp{num}input"
        epp_str = epp_strings[i]

        epp_str = db_r1cols[i].text_input(label='', placeholder=epp_msg, key=key_txt).strip()
        if epp_str:
            try:
                epp_files[i] = util.guard_pbc(epp_str)
            except Exception as exc:
                print(f'load_database_data() {epp_msg=} {exc}', file=sys.stderr)

    #Row2 column 1: DNA plate ID input
    dnap_str = db_r2cols[0].text_input(label='', placeholder=dnap_msg, key='dnap_input').strip()
    if dnap_str:
        try:
            dnap = file_io.guard_pbc(dnap_str)
        except Exception as exc:
            print(f'load_database_data() {dnap_msg=} {exc}', file=sys.stderr)

    # nothing in row2col2

    #Row2 column 3: Search button
    db_r2cols[2].markdown('')
    db_r2cols[2].markdown('')
    add_musterer_button = db_r2cols[3].button('Search Musterer')
    if add_musterer_button:
        epps = [epp.strip() for epp in epp_files if epp.strip() != '']
        for epp in epps:
            success = db_io.get_plate_musterer(st.session_state['experiment'].name, epp, True)
        success = st.session_state['experiment'].add_musterer_plate_set(epps,dnap)
        if not success:
            db_r2cols[2].markdown('<p style="color:#FF0000">Failed to incorporate manifest. Please read the log.</p>', unsafe_allow_html=True)
        st.experimental_rerun() #?

    db_r2cols[3].markdown('')
    db_r2cols[3].markdown('')



def upload_pcr_files(PCR_stage = 'all'):
    """
    Make big changes to this
    Home page: load data

    Upload options for extra custom data including reference files, assay lists, primer plates, index plates, taq/water plates, and PCR plates (barcode only)

    PCR stage = [1, 2, 'all']
    As stage 1, primer plate layouts and volumes, taq and water plates and/or barcodes, and PCR plate barcodes are required.. 
    Users can also upload custom reference files, assay lists and extra taq and water plates. 

    At stage 2, user provides amplicon plates, index plates and taq/water plates (?)

    """
    upload_col1, upload_col2 = st.columns(2)

    if PCR_stage == 1:
        
        uploaded_primer_layouts = upload_col1.file_uploader('Upload Primer Plate Layouts', \
            key='primer_uploader', type='csv', accept_multiple_files=True)

        uploaded_primer_volumes = upload_col2.file_uploader('Upload Primer Plate Volumes', \
            key='primer_vol_uploader', type='csv', accept_multiple_files=True)

        pcr_plate_barcode = upload_col2.text_input('PCR Plate Barcode', \
            placeholder='Type in the barcode for PCR plate', key='pcr_plate_barcode')


        if uploaded_primer_layouts:
            if 'primer_layouts_upload' not in st.session_state or \
                st.session_state['primer_layouts_upload'] != [upl.name for upl in uploaded_primer_layouts]:
                success = st.session_state['experiment'].add_primer_layouts(uploaded_primer_layouts)
                st.session_state['primer_layouts_upload'] = [upl.name for upl in uploaded_primer_layouts]
        
        if uploaded_primer_volumes:
            if 'primer_volumes_upload' not in st.session_state or \
                st.session_state['primer_volumes_upload'] != [upv.name for upv in uploaded_primer_volumes]:
                success = st.session_state['experiment'].add_primer_volumes(uploaded_primer_volumes)
                st.session_state['primer_volumes_upload'] = [upv.name for upv in uploaded_primer_volumes]

        if pcr_plate_barcode:
            # compare file name against the previously loaded file. Only runs if the file changes
            if 'pcr_barcode_upload' not in st.session_state or \
                    st.session_state['pcr_barcode_upload'] != pcr_plate_barcode:
                success = st.session_state['experiment'].add_pcr_plates([pcr_plate_barcode])
                # set the uploaded file in cache so we don't accidentally add it again
                st.session_state['pcr_barcode_upload'] = pcr_plate_barcode


    if PCR_stage == 2:
        uploaded_index_layouts = upload_col1.file_uploader('Upload i7i5 Index Plate Layout', \
            key='index_uploader', type='csv', accept_multiple_files=True)

        uploaded_index_volumes = upload_col2.file_uploader('Upload i7i5 Index Plate Volumes', \
            key='index_vol_uploader', type='csv', accept_multiple_files=True)

        uploaded_amplicon_plates = upload_col2.file_uploader('Upload extra amplicon plates', \
            key='amplicon_plate_uploader', type='csv', accept_multiple_files=True)


        if uploaded_index_layouts:
            if 'index_layout_upload' not in st.session_state or \
                st.session_state['index_layout_upload'] != [uil.name for uil in uploaded_index_layouts]:
                success = st.session_state['experiment'].add_index_layouts(uploaded_index_layouts)
                st.session_state['barcode_layout_upload'] = [uil.name for uil in uploaded_index_layouts]
                            
        
        if uploaded_index_volumes:
            if 'index_volume_upload' not in st.session_state or \
                    st.session_state['index_volume_upload'] !=\
                    [uiv.name for uiv in uploaded_index_volumes]:
                success = st.session_state['experiment'].add_index_volumes(uploaded_index_volumes)
                st.session_state['index_volume_upload'] = [uiv.name for uiv in uploaded_index_volumes]
        
        #improv
        if uploaded_amplicon_plates:
            if 'amplicon_plate_upload' not in st.session_state or \
                    st.session_state['amplicon_plate_upload'] != \
                    [uap.name for uap in uploaded_amplicon_plates]:
                success = st.session_state['experiment'].add_amplicon_manifests(uploaded_amplicon_plates)
                st.session_state['amplicon_plate_upload'] = [uap.name for uap in uploaded_amplicon_plates]
     
    taqwater_barcode = upload_col1.text_input('Taq and Water Plate Barcode', \
        placeholder='Type in the barcode for taq+water plate', key='taq_water_barcode')

    upload_col1.write('')
    upload_col1.write('')
    upload_col1.write('')
    upload_col1.write('')


    uploaded_taqwater_plates = upload_col1.file_uploader('Upload Taq and Water Plates', \
        key='taq_water_upload', type='csv', accept_multiple_files=True)

    upload_col2.write('')
    upload_col2.write('')
    upload_col2.write('')
    upload_col2.write('')

    uploaded_references = upload_col1.file_uploader('Upload Custom Reference File', \
        key='ref_uploader', type=['txt','fa','fasta'], accept_multiple_files=True)

    uploaded_assaylists = upload_col2.file_uploader('Upload Assay List', \
        key='assaylist_uploader', type=['txt','csv'], accept_multiple_files=True)


    # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    if uploaded_references:
        if 'reference_upload' not in st.session_state or \
                st.session_state['reference_upload'] != [upr.name for upr in uploaded_references]:
            success = st.session_state['experiment'].add_references(uploaded_references)
            # set the list of uploaded files so we don't accidentally add them again
            st.session_state['reference_upload'] = [upr.name for upr in uploaded_references]
        
    if uploaded_assaylists:
        if 'assaylist_upload' not in st.session_state or \
            st.session_state['assaylist_upload'] != [ual.name for ual in uploaded_assaylists]:
            success = st.session_state['experiment'].add_assaylists(uploaded_assaylists)
            st.session_state['assaylist_upload'] = [ual.name for ual in uploaded_assaylists]
         
    if taqwater_barcode:
        if 'taqwater_barcode_upload' not in st.session_state or \
                st.session_state['taqwater_barcode_upload'] != taqwater_barcode:
            success = st.session_state['experiment'].add_standard_taqwater_plates([taqwater_barcode])
            st.session_state['taqwater_barcode_upload'] = taqwater_barcode

    if uploaded_taqwater_plates:
        if 'taqwater_upload' not in st.session_state or \
                st.session_state['taqwater_upload'] != [utp.name for utp in uploaded_taqwater_plates]:  
            success = st.session_state['experiment'].add_standard_taqwater_plates(uploaded_taqwater_plates)
            st.session_state['taqwater_upload'] = [utp.name for utp in uploaded_taqwater_plates]


    #with st.container():
    #    with st.expander('Upload reference sequences, assay lists, taq and water plates, primer plates, barcode plates, PCR plate barcodes', expanded=False):

    #        upload_option = st.radio('Upload Options:', ('Custom Reference File', 'Assay List', 'Taq and Water Plates', 
    #                'Primer Plates', 'Barcode Plates', 'PCR Plates (barcodes only)'))

    #        if upload_option == 'Custom Reference File':
    #            uploaded_references = st.file_uploader('Upload custom reference file', key='ref_uploader', type=['txt','fa','fasta'], accept_multiple_files=True)
    #            if uploaded_references:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'reference_upload' not in st.session_state or \
    #                        st.session_state['reference_upload'] != [upr.name for upr in uploaded_references]:
    #                    success = st.session_state['experiment'].add_references(uploaded_references)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['reference_upload'] = [upr.name for upr in uploaded_references]
                

    #        if upload_option == 'Assay List':
    #            uploaded_assaylists = st.file_uploader('Upload assay list', key='assaylist_uploader', type=['txt','csv'], accept_multiple_files=True)
    #            if uploaded_assaylists:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'assaylist_upload' not in st.session_state or \
    #                        st.session_state['assaylist_upload'] != [upa.name for upa in uploaded_assaylists]:
    #                    success = st.session_state['experiment'].add_assaylists(uploaded_assaylists)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['assaylist_upload'] = [upa.name for upa in uploaded_assaylists]

    #        if upload_option == 'Primer Plates':
    #            uploaded_primer_layouts = st.file_uploader('Upload primer plate layouts', 
    #                    key='primer_uploader', type='csv', accept_multiple_files=True)
    #            if uploaded_primer_layouts:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'primer_layouts_upload' not in st.session_state or \
    #                        st.session_state['primer_layouts_upload'] != [upl.name for upl in uploaded_primer_layouts]:
    #                    success = st.session_state['experiment'].add_primer_layouts(uploaded_primer_layouts)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['primer_layouts_upload'] = [upl.name for upl in uploaded_primer_layouts]
    #            uploaded_primer_volumes = st.file_uploader('Upload primer plate volumes', key='primer_vol_uploader', type='csv', accept_multiple_files=True)
    #            if uploaded_primer_volumes:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'primer_volumes_upload' not in st.session_state or \
    #                        st.session_state['primer_volumes_upload'] != [upv.name for upv in uploaded_primer_volumes]:
    #                    success = st.session_state['experiment'].add_primer_volumes(uploaded_primer_volumes)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['primer_volumes_upload'] = [upv.name for upv in uploaded_primer_volumes]

    #        if upload_option == 'Index Plates':
    #            uploaded_barcode_layouts = st.file_uploader('Upload i7i5 barcode plate layout', 
    #                    key='barcode_uploader', type='csv', accept_multiple_files=True)
    #            if uploaded_barcode_layouts:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'barcode_layout_upload' not in st.session_state or \
    #                        st.session_state['barcode_layout_upload'] != [ubl.name for ubl in uploaded_barcode_layouts]:
    #                    success = st.session_state['experiment'].add_barcode_layouts(uploaded_barcode_layouts)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['barcode_layout_upload'] = [ubl.name for ubl in uploaded_barcode_layouts]
                            
    #            uploaded_barcode_volumes = st.file_uploader('Upload i7i5 barcode plate volumes', key='barcode_vol_uploader',
    #                    type='csv', accept_multiple_files=True)
    #            if uploaded_barcode_volumes:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'barcode_volume_upload' not in st.session_state or \
    #                        st.session_state['barcode_volume_upload'] != [ubv.name for ubv in uploaded_barcode_volumes]:
    #                    success = st.session_state['experiment'].add_barcode_volumes(uploaded_barcode_volumes)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['barcode_volume_upload'] = [ubv.name for ubv in uploaded_barcode_volumes]

    #        #Option to add taq water
    #        if upload_option == 'Taq and Water Plates':
    #            taqwater_barcode = st.text_input('Taq+Water plate barcode', placeholder='Type in the barcode for taq+water plate', key='taq_water_barcode')
    #            if taqwater_barcode:
    #                # compare file name against the previously loaded file. Only runs if the file changes
    #                if 'taqwater_barcode_upload' not in st.session_state or \
    #                        st.session_state['taqwater_barcode_upload'] != taqwater_barcode:
    #                    success = st.session_state['experiment'].add_standard_taqwater_plates([taqwater_barcode])
    #                    # set the uploaded file in cache so we don't accidentally add it again
    #                    st.session_state['taqwater_barcode_upload'] = taqwater_barcode

    #            uploaded_taqwater_plates = st.file_uploader('Upload taq and water resevoir plates', 
    #                    key='taq_water_upload', type='csv', accept_multiple_files=True)
    #            if uploaded_taqwater_plates:
    #                # compare list of all names against the previously loaded list. Only runs if the whole file list changes
    #                if 'taqwater_upload' not in st.session_state or \
    #                        st.session_state['taqwater_upload'] != [utp.name for utp in uploaded_taqwater_plates]:  
    #                    success = st.session_state['experiment'].add_standard_taqwater_plates(uploaded_taqwater_plates)
    #                    # set the list of uploaded files so we don't accidentally add them again
    #                    st.session_state['taqwater_upload'] = [utp.name for utp in uploaded_taqwater_plates]

    #        if upload_option == 'PCR Plates (barcodes only)':
    #            pcr_plate_barcode = st.text_input('PCR plate barcode', placeholder='Type in the barcode for PCR plate', key='pcr_plate_barcode')
    #            if pcr_plate_barcode:
    #                # compare file name against the previously loaded file. Only runs if the file changes
    #                if 'pcr_barcode_upload' not in st.session_state or \
    #                        st.session_state['pcr_barcode_upload'] != pcr_plate_barcode:
    #                    success = st.session_state['experiment'].add_pcr_plates([pcr_plate_barcode])
    #                    # set the uploaded file in cache so we don't accidentally add it again
    #                    st.session_state['pcr_barcode_upload'] = pcr_plate_barcode

def upload_miseq_fastqs(exp):
    st.markdown('<h4 style="color:#000000">Add Miseq FASTQ files to experiment</h4>', unsafe_allow_html=True)
    fastq_path = dc.st_directory_picker("Select location of Miseq FASTQ files")
    fastq_files = [f for f in fastq_path.glob('*.fastq*')] + [f for f in fastq_path.glob('*.fq*')]
    if len(fastq_files) > 0:
        st.write(f"1. {str(fastq_files[0])}")
        if len(fastq_files) > 1:
            st.write(f"...")
            st.write(f"{len(fastq_files)}. {str(fastq_files[-1])}")
    import_fastqs = st.button('Import FASTQs')
    file_field = st.empty()
    copy_progress = st.progress(0)
    if import_fastqs and len(fastq_files) > 0:
        for i,fp in enumerate(fastq_files):
            file_field.markdown(f'<hr5>Copying {fp} to experiment</h5>', unsafe_allow_html=True)
            copy_progress.progress(i/len(fastq_files))
            copy2(fp, exp.get_raw_dirpath())
        file_field.markdown('<h5>Done</h5>', unsafe_allow_html=True)
        copy_progress.progress(100)