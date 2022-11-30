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

import os
import sys
from pathlib import PurePath, Path
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO
from shutil import copy2
from time import sleep

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components
import extra_streamlit_components as stx

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
    Manage combining up to four rodentity plates (from JSON) with one destination PID 
    """
    exp = st.session_state['experiment']
    plates_to_clear = [False, False, False, False]
    with st.expander('Add data from Rodentity JSON files',expanded=True):
        if all([exp.unassigned_plates[unass_pid] != '' for unass_pid in exp.unassigned_plates]):
            disable_uploads = True
        else:
            disable_uploads = False

        with st.form('rodentity_upload_form'):
            epps_col1, epps_col2 = st.columns([1,2])
            rodentity_epps = epps_col1.file_uploader('Choose up to four Rodentity JSON files', type='json', accept_multiple_files=True, disabled=disable_uploads)
            rod_upload_submit_button = st.form_submit_button('Submit')
            if rod_upload_submit_button and rodentity_epps:
                used_keys = set()
                for rod_epp in rodentity_epps:
                    rod_pid = util.guard_pbc(rod_epp.name.rstrip('.json'), silent=True)
                    if all([exp.unassigned_plates[slot] != '' for slot in [1,2,3,4]]):
                        st.markdown(f'<p style="color:#FF0000">Ran out of free slots for {util.unguard_pbc(rod_pid, silent=True)}</p>', unsafe_allow_html=True)
                        continue
                    if rod_pid in exp.dest_sample_plates or rod_pid in exp.plate_location_sample or \
                            rod_pid in exp.unassigned_plates[1] or rod_pid in exp.unassigned_plates[2] or \
                            rod_pid in exp.unassigned_plates[3] or rod_pid in exp.unassigned_plates[4] or \
                            rod_pid in exp.unassigned_plates['custom']:
                        st.markdown('<p style="color:#FF0000">Rodentity plate barcode already in use: ' +\
                                util.unguard_pbc(rod_pid, silent=True) + '</p>', unsafe_allow_html=True)
                        sleep(0.5)
                        continue                     
                    for key in exp.unassigned_plates:
                        if key not in used_keys and key != 'custom' and key != 'None':
                            exp.unassigned_plates[key] = rod_pid
                            used_keys.add(key)
                            # Save a copy we can edit and reload
                            rod_fn = os.path.join('run_'+exp.name, db_io._eppfn_r(rod_pid))
                            with open(rod_fn, 'wt') as outf:
                                outf.write(rod_epp.getvalue().decode("utf-8"))
                            exp.log(f'Info: uploaded copy Rodentity plate {util.unguard_pbc(rod_pid, silent=True)} to {rod_fn}')
                            break
                exp.save()
    
        rod_col1, rod_col2, _ = st.columns([1,1,1])
        with rod_col1.form('set_rod_plates_form', clear_on_submit=True):
            for i in range(4):
                plates_to_clear[i] = st.checkbox(f"P{str(i+1)}: {util.unguard_pbc(exp.unassigned_plates[i+1], silent=True)}", 
                        help='Click the checkbox to allow a set plate ID to be cleared', key='check'+str(i+1))
            clear_plates_button = st.form_submit_button('Clear IDs', help='Clear selected Rodentity plate IDs')
            if clear_plates_button:
                for i, plate in enumerate(plates_to_clear):
                    if plate and exp.unassigned_plates[i+1]:
                        exp.unassigned_plates[i+1] = ''
                        plates_to_clear[i] = False
                exp.save()
                

        with rod_col2.form('rod_destination_form', clear_on_submit=True):
            rod_dp = st.text_input('Destination plate ID (barcode)', max_chars=30, key='rod_dp_key')
            if rod_dp:
                rod_dp = util.guard_pbc(rod_dp, silent=True)
            accept_rod_dest_button = st.form_submit_button('Accept')

            if accept_rod_dest_button and rod_dp:
                if rod_dp in exp.dest_sample_plates or rod_dp in exp.plate_location_sample or \
                        rod_dp in exp.unassigned_plates[1] or rod_dp in exp.unassigned_plates[2] or \
                        rod_dp in exp.unassigned_plates[3] or rod_dp in exp.unassigned_plates[4] or \
                        rod_dp in exp.unassigned_plates['custom']:
                    st.markdown('<p style="color:#FF0000">Destination plate barcode already in use: ' +\
                            util.unguard_pbc(rod_dp, silent=True) + '</p>', unsafe_allow_html=True)
                    sleep(1.5)
                else:
                    success = exp.add_rodentity_plate_set([exp.unassigned_plates[k] for k in [1,2,3,4] if exp.unassigned_plates[k]], rod_dp)
                    if not success:
                        st.markdown('<p style="color:#FF0000">Failed to incorporate plate set. Please read the log.</p>', unsafe_allow_html=True)
                        sleep(1.5)
                    else:
                        st.write('Successfully added plate set')
                        sleep(1)
                        exp.unassigned_plates = {1:'',2:'',3:'',4:''}
                        exp.save()
                        st.experimental_rerun()
                        
                   
             
#def load_custom_csv(expanded=False):
#    '''
#    Upload data from a custom csv file - DEPRECATED
#    '''
#    custom_csv_exp = st.expander('Add data from a custom manifest CSV file', expanded=expanded)
#    custom_col1,custom_col2 = custom_csv_exp.columns(2)
    
    
#    st.session_state['default_manifest_type'] = 'c'
#    manifest_choice = custom_col1.radio('The default contents of my manifest file (if not declared) are', 
#                        ['Custom', 'Rodentity mouse', 'Musterer mouse'])

#    uploaded_file = custom_col2.file_uploader('', type='csv')

#    if manifest_choice == 'Rodentity mouse':
#        st.session_state['default_manifest_type'] = 'r'
#    elif manifest_choice == 'Musterer mouse':
#        st.session_state['default_manifest_type'] = 'm'

#    if uploaded_file:
#        if 'custom_upload' not in st.session_state or st.session_state['custom_upload'] != uploaded_file.name:
#            success = st.session_state['experiment'].add_manifest(uploaded_file, st.session_state['default_manifest_type'])
#            if not success:
#                custom_col2.markdown('<p style="color:#FF0000">Failed to incorporate manifest. Please read the log.</p>', unsafe_allow_html=True)
#            st.session_state['custom_upload'] = uploaded_file.name
#            #?
#            st.experimental_rerun()

def load_custom_manifests():
    """ Demonstration upload panel for custom manifests """
    exp = st.session_state['experiment']
    if 'custom' not in exp.unassigned_plates or not exp.unassigned_plates['custom']:
        exp.unassigned_plates['custom'] = {'None':{}}
    with st.expander('Upload custom manifests', expanded=True):
        with st.form(key='manifest_upload_form', clear_on_submit=True):
            col1, col3 = st.columns([1,2])
            with col1:
                manifest_uploads = st.file_uploader(label="Custom manifest upload", accept_multiple_files=True)
                #st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                submit_manifest = st.form_submit_button()
                if submit_manifest and manifest_uploads:
                    success = exp.read_custom_manifests(manifest_uploads)
                    if not success:
                        st.markdown('<p style="color:#FF0000">Failed to read custom manifests.'+\
                                'Please see the log</p>', unsafe_allow_html=True)
                    if success:
                        st.write('Successfully read custom manifests')
       
        with st.form(key='selection_form', clear_on_submit=True):
            st.write('Select up to four 96-well sample plate IDs (barcodes) to combine into a 384-well DNA plate')
            col1, col2, _, col3, _ = st.columns([2,2,1,2,2])
            plate_options = [util.unguard_pbc(pid, silent=True) for pid in exp.unassigned_plates['custom'].keys()]
            with col1:
                p1 = util.guard_pbc(st.selectbox(label='Plate 1', options=plate_options, key='s1'), silent=True)
                p2 = util.guard_pbc(st.selectbox(label='Plate 2', options=plate_options, key='s2'), silent=True)
            with col2:
                p3 = util.guard_pbc(st.selectbox(label='Plate 3', options=plate_options, key='s3'), silent=True)
                p4 = util.guard_pbc(st.selectbox(label='Plate 4', options=plate_options, key='s4'), silent=True)
            with col3:
                dest_pid = st.text_input(label='Destination plate ID (barcode)', key='dest1')
                if dest_pid:
                    dest_pid = util.guard_pbc(dest_pid, silent=True)
                    if dest_pid in exp.dest_sample_plates or dest_pid in exp.plate_location_sample or \
                            dest_pid in exp.unassigned_plates[1] or dest_pid in exp.unassigned_plates[2] or \
                            dest_pid in exp.unassigned_plates[3] or dest_pid in exp.unassigned_plates[4] or \
                            dest_pid in exp.unassigned_plates['custom']:
                        st.markdown('<p style="color:#FF0000">WARNING! Destination plate barcode already in use: ' +\
                                util.unguard_pbc(dest_pid, silent=True) + '</p>', unsafe_allow_html=True)
                    else:
                        st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                submit_selection = st.form_submit_button("Accept")
                if submit_selection and dest_pid:
                    sample_plates = [pid for pid in [p1,p2,p3,p4] if pid != util.guard_pbc('None', silent=True) and pid is not None]
                    print(f'{sample_plates=}', file=sys.stderr)
                    success = exp.accept_custom_manifests(dest_pid, sample_plates)
                    if not success:
                        st.markdown('<p style="color:#FF0000">Failed to assign custom plates.'+\
                                'Please see the log</p>', unsafe_allow_html=True)
                    if success:
                        st.write('Successfully assigned custom plates')
    exp.save()


#def load_database_data():
#    """
#    Import Musterer data from the webservice. 
#    Keep as a template for eventual Rodentity web service 
#    """
#    database_exp = st.expander('Add data from a database')
#    db_r1cols = database_exp.columns(4)
#    db_r2cols = database_exp.columns(4)

#    epp1 = epp2 = epp3 = epp4 = dnap = ''
#    epp1_str = epp2_str = epp3_str = epp4_str = None
#    epp_files = [epp1, epp2, epp3, epp4]
#    epp_strings = [epp1_str, epp2_str, epp3_str, epp4_str]
#    dnap_msg = 'DNA plate barcode'

#    #Create buttons for the four Ear punch plate and file (row 1 columns) - not sure if it works yet
#    for i in range(4):
#        num=str(i+1)
#        epp_msg = f"Ear punch plate {num} barcode"
#        key_txt = f"epp{num}input"
#        epp_str = epp_strings[i]

#        epp_str = db_r1cols[i].text_input(label='', placeholder=epp_msg, key=key_txt).strip()
#        if epp_str:
#            try:
#                epp_files[i] = util.guard_pbc(epp_str)
#            except Exception as exc:
#                print(f'load_database_data() {epp_msg=} {exc}', file=sys.stderr)

#    #Row2 column 1: DNA plate ID input
#    dnap_str = db_r2cols[0].text_input(label='', placeholder=dnap_msg, key='dnap_input').strip()
#    if dnap_str:
#        try:
#            dnap = file_io.guard_pbc(dnap_str)
#        except Exception as exc:
#            print(f'load_database_data() {dnap_msg=} {exc}', file=sys.stderr)

#    # nothing in row2col2

#    #Row2 column 3: Search button
#    db_r2cols[2].markdown('')
#    db_r2cols[2].markdown('')
#    add_musterer_button = db_r2cols[3].button('Search Musterer')
#    if add_musterer_button:
#        epps = [epp.strip() for epp in epp_files if epp.strip() != '']
#        for epp in epps:
#            success = db_io.get_plate_musterer(st.session_state['experiment'].name, epp, True)
#        success = st.session_state['experiment'].add_musterer_plate_set(epps,dnap)
#        if not success:
#            db_r2cols[2].markdown('<p style="color:#FF0000">Failed to incorporate manifest. Please read the log.</p>', unsafe_allow_html=True)
#        st.experimental_rerun() #?

#    db_r2cols[3].markdown('')
#    db_r2cols[3].markdown('')

def provide_barcodes(key):
    """
    :param key: unique key for the input widgets
    """
    exp = st.session_state['experiment']
    pcr_col, taqwater_col = st.columns(2)
    with pcr_col.form('Add PCR plate barcode', clear_on_submit=True):
        pcr_plate_barcode = st.text_input('PCR Plate Barcode', \
                placeholder='Type in the barcode for PCR plate', key='pcr_plate_barcode' + key)
        pcr_plate_entered = st.form_submit_button('Add')

        if pcr_plate_entered and pcr_plate_barcode:
            guarded_pcr_plate_barcode = util.guard_pbc(pcr_plate_barcode, silent=True)
            if guarded_pcr_plate_barcode not in exp.plate_location_sample:
                success = exp.add_pcr_plates([guarded_pcr_plate_barcode])
                if success:
                    st.write(f'Successfully added PCR plate barcode {pcr_plate_barcode}')
                    sleep(1.5)
                    st.experimental_rerun()
                else:
                    st.write(f'Could not add PCR plate barcode {pcr_plate_barcode}, please see the log')
            else:
                st.write(f'This plate barcode {pcr_plate_barcode} appears to already be in use')

    with taqwater_col.form('Add taq+water plate barcode', clear_on_submit=True):
        taqwater_plate_barcode = st.text_input('Taq and Water Plate Barcode', \
                placeholder='Type in the barcode for taq+water plate', key='taq_water_barcode' + key)
        taqwater_plate_entered = st.form_submit_button('Add')

        if taqwater_plate_entered and taqwater_plate_barcode:
            guarded_taqwater_plate_barcode = util.guard_pbc(taqwater_plate_barcode, silent=True)
            if guarded_taqwater_plate_barcode not in exp.plate_location_sample:
                success = exp.add_standard_taqwater_plates([taqwater_plate_barcode])
                if success:
                    st.write(f'Successfully added taq+water plate barcode {taqwater_plate_barcode}')
                    sleep(1.5)
                    st.experimental_rerun()
                else:
                    st.write(f'Could not add taq+water plate barcode {taqwater_plate_barcode}, please see the log')
            else:
                st.write(f'This plate barcode {taqwater_plate_barcode} appears to already be in use')




def upload_pcr1_files(key):
    """
    Uploads for primer layout and volumes
    """
    exp = st.session_state['experiment']
    with st.form('primer plate upload'+key, clear_on_submit=True):
        col1, col2 = st.columns(2)
        uploaded_primer_layouts = col1.file_uploader('Upload Primer Plate Layouts - the barcode must be the first part of the filename e.g. 12345_primer_layout.csv', \
            key='primer_layout_uploader'+key, type='csv', accept_multiple_files=True) 
        uploaded_primer_volumes = col2.file_uploader('Upload Primer Plate Volumes - the barcode must be the first part of the filename e.g. 12345_primer_volume.csv', \
            key='primer_vol_uploader'+key, type='csv', accept_multiple_files=True)
        upload_button = st.form_submit_button("Upload Files")

        if upload_button:
            if uploaded_primer_layouts:
                upl_pids = [upl.name for upl in uploaded_primer_layouts]
                success = exp.add_primer_layouts(uploaded_primer_layouts)
                if success:
                    st.write(f'Successfully added primer layouts for plates {upl_pids}')
                else:
                    st.write(f'Failed to write at least one primer layout, please see the log')
                     
            if uploaded_primer_volumes:
                upv_pids = [upv.name for upv in uploaded_primer_volumes]
                success = exp.add_primer_volumes(uploaded_primer_volumes)
                if success:
                    st.write(f'Successfully added primer volumes for plates {upv_pids}')
                else:
                    st.write(f'Failed to write at least one set of primer volumes, please see the log')

 
def upload_pcr2_files(key):
    """
    Upload inputs for indexing layout and volume. Extra option for amplicon plate upload.
    """
    exp = st.session_state['experiment']
    with st.form('index plate upload'+key, clear_on_submit=True):
        col1, col2 = st.columns(2)
        uploaded_index_layouts = col1.file_uploader('Upload i7i5 Index Plate Layout - the barcode must be the first part of the filename e.g. 12345_index_layout.csv',
            key='index_layout_uploader'+key, type='csv', accept_multiple_files=True)
 
        uploaded_index_volumes = col2.file_uploader('Upload i7i5 Index Plate Volumes - the barcode must be the first part of the filename e.g. 12345_index_volume.csv', 
            key='index_vol_uploader'+key, type='csv', accept_multiple_files=True)
 
        uploaded_amplicon_plates = col1.file_uploader('Upload Extra Amplicon Plates', 
                key='amplicon_plate_uploader'+key, type='csv', accept_multiple_files=True)
        upload_button = st.form_submit_button("Upload Files")

        if upload_button:
            if uploaded_index_layouts:
                uil_pids = [uil.name for uil in uploaded_index_layouts]
                success = st.session_state['experiment'].add_index_layouts(uploaded_index_layouts)
                if success:
                    st.write(f'Successfully added index layouts for plates {uil_pids}')
                else:
                    st.write(f'Failed to write at least one index layout, please see the log')    

            if uploaded_index_volumes:
                uiv_pids = [uiv.name for uiv in uploaded_index_volumes]
                success = st.session_state['experiment'].add_index_volumes(uploaded_index_volumes)
                if success:
                    st.write(f'Successfully added index volumes for plates {uiv_pids}')
                else:
                    st.write(f'Failed to write at least one set of index volumes, please see the log')    

            if uploaded_amplicon_plates:
                uap_pids = [uap.name for uap in uploaded_amplicon_plates]
                success = exp.add_amplicon_manifests(uploaded_amplicon_plates)
                if success:
                    st.write(f'Successfully added amplicon manifests from files {uap_pids}')
                else:
                    st.write(f'Failed to upload at least one amplicon manifest, please see the log')


def upload_extra_consumables(key):
    """
    Uploads for extra files such as custom references, assay lists, and taq/water plates
    """
    exp = st.session_state['experiment']
    with st.form('Consumables upload'+key, clear_on_submit=True):
        col1, col2 = st.columns(2)
        uploaded_references = col1.file_uploader('Upload Custom Reference Files', key='ref_uploader'+key, 
                type=['txt','fa','fasta'], accept_multiple_files=True)
                                                                        
        uploaded_assaylists = col2.file_uploader('Upload Assay Lists', key='assaylist_uploader'+key, 
                type=['txt','csv'], accept_multiple_files=True)
                                                                                     
        #uploaded_taqwater_plates = col1.file_uploader('Upload Taq and Water Plates', key='taq_water_upload'+key, 
        #        type='csv', accept_multiple_files=True)
        upload_button = st.form_submit_button("Upload Files")

        if upload_button:
            if uploaded_references:
                success = exp.add_references(uploaded_references)
                ref_names = [ur.name for ur in uploaded_references]
                if success:
                    st.write(f'Successfully added reference sequences from files {ref_names}')
                else:
                    st.write(f'Failed to upload at least one reference sequence file, please see the log')
            if uploaded_assaylists:
                success = exp.add_assaylists(uploaded_assaylists)
                assaylist_names = [ual.name for ual in uploaded_assaylists]
                if success:
                    st.write(f'Successfully added assay/primer lists from files {assaylist_names}')
                else:
                    st.write(f'Failed to upload at least one assay/primer list file, please see the log')
            # TODO: add the ability to add a custom taq+water plate
            #if uploaded_taqwater_plates:
            #    success = exp.add_taqwater_layout_volume


 
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


#def upload_rodentity_demo():
#    """ Demonstration upload panel for Rodentity plates - not in use for now """
#    exp = st.session_state['experiment']
#    with st.expander('Upload Rodentity Plates', expanded=False):
#        with st.form(key='rodentity_upload_form', clear_on_submit=True):
#            col1, col2, col3 = st.columns([4,4,1])
#            with col1:
#                rod_up1 = st.file_uploader(label='Rodentity plate upload', key='rod1', accept_multiple_files=False)
#                rod_up2 = st.file_uploader(label='Rodentity plate upload', key='rod2', accept_multiple_files=False)
#            with col2:
#                rod_up3 = st.file_uploader(label='Rodentity plate upload', key='rod3', accept_multiple_files=False)
#                rod_up4 = st.file_uploader(label='Rodentity plate upload', key='rod4', accept_multiple_files=False)
#            with col1:
#                dest_pid = st.text_input(label='Destination plate ID (barcode)', key='dest1')
#            with col3:
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
#                submit_rodentity = st.form_submit_button()
#            if submit_rodentity:
#                if not dest_pid:
#                    with col2:
#                        st.markdown('<p style="color:#FF0000">Please provide the destination plate barcode</p>', unsafe_allow_html=True)
#                if not rod_up1 and not rod_up2 and not rod_up3 and not rod_up4:
#                    with col2:
#                        st.markdown('<p style="color:#FF0000">Please provide at least one Rodentity plate file</p>', unsafe_allow_html=True)
#                else:
#                    if dest_pid:
#                        rod_ups = [rod_up for rod_up in [rod_up1, rod_up2, rod_up3, rod_up4] if rod_up is not None]
#                        success = exp.add_rodentity_plate_set(rod_ups, dest_pid)
#                    if not success:
#                        with col2:
#                            st.markdown('<p style="color:#FF0000">Failed to incorporate plate set. Please read the log.</p>', unsafe_allow_html=True)
#                    else:
#                        with col2:
#                            st.write('Successfully added plate set')






