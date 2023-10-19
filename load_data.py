#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University

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
from unittest import expectedFailure

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components
import extra_streamlit_components as stx
import stutil as stutil

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
    import bin.parse as parse
except ModuleNotFoundError:
    import parse
try:
    import bin.generate as generate
except ModuleNotFoundError:
    import generate
try:
    import bin.db_io as db_io
except ModuleNotFoundError:
    import db_io
try:
    from bin.makehtml import generate_heatmap_html
except ModuleNotFoundError:
    from makehtml import generate_heatmap_html

try: 
    import bin.transaction as trans
except ModuleNotFoundError:
    import transaction as trans

import display_components as dc


credits="""
@created: March 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy, 2019-2023

Web application style interface using Streamlit. Present information as both a dashboard and a workflow in one page.

The GUI interacts with a single Experiment object at one time. Methods are called on this to activate pipeline
functionality. The Experiment then deals directly with the pipeline logic.
"""

HELP_COLOUR = '#7f8b8f'
FORM_BUTTON_COLOUR = '#4287f5'

def do_pending(combined_pending):
    """
    Action the pending_file_widget
    """
    exp = st.session_state['experiment']
    warning_area = st.container()
    expected_keys = ['pending_file_checkbox_'+str(pending) for pending in combined_pending]
    pending_checked = []
    # test all the checkboxes to see if they were ticked
    for key in expected_keys:
        if key in st.session_state:
            if st.session_state[key]:
                 pending_checked.append(True)
            else:
                pending_checked.append(False)
        else:
             pending_checked.append(False)
    # reset all the checkbox keys in case they are used again
    for key in expected_keys:
        if key in st.session_state:
            st.session_state[key] = False

    # go through each pending item to see if it is kept or discarded
    for pending, pc in zip(combined_pending, pending_checked):
        if pc:
            clashing_filenames, clashing_pids = trans.check_for_clashing_transactions(exp, filenames=pending)
            msg = f'Overwriting the previous version of {str(pending)} with the new file'
            warning_area.info(msg)
            exp.log(f'Info: {msg}')
            if pending in exp.pending_uploads:
                success = parse.accept_pending_upload(exp, pending)
            elif pending in exp.pending_steps:
                success = trans.accept_pending_transaction(exp, pending)
            if success:
                exp.log(f'Success: overwrote existing file {pending}')
                if clashing_filenames:
                    msg = f'{", ".join(clashing_filenames)} potentially affected by change'
                    exp.log(msg)
                    warning_area.warning(msg)
                if clashing_pids:
                    msg = f'{", ".join(clashing_pids)} potentially affected by change'
                    exp.log(msg)
                    warning_area.warning(msg)
            else:
                msg = f'Failure: Could not overwrite existing file {pending}'
                exp.log(msg)
                warning_area.error(msg)
            
        else:
            warning_area.info(f'Keeping the old version of the file {pending}')
            exp.log(f'Info: removing unwanted pending file {pending}')
            if pending in exp.pending_uploads:
                exp.pending_uploads.remove(pending)
                if pending in exp.uploaded_files:
                    success = exp.del_file_record(pending)
            elif pending in exp.pending_steps:
                success = trans.clear_pending_transaction(exp, pending)
            if success:
                exp.log(f'Success: pending file or transaction {pending} cleared')
            else:
                msg = f'Failure: pending file or transaction {pending} could not be removed'
                exp.log(msg)
                warning_area.error(msg)

    exp.save()          


def pending_file_widget(key):
    """
    Allow user input via checkboxes to deal with clashes.
    Requires st.session_state['clashes'] from previously calling set_pending_file_clashes()
    """
    exp = st.session_state['experiment']
    warning_area = st.container()
    combined_pending = []  # both pending_uploads and pending_steps
    
    pending_uploads = {}
    if exp.pending_uploads:
        pending_uploads = exp.pending_uploads.copy()
    pending_steps = {}
    if exp.pending_steps:
        pending_steps = exp.pending_steps.copy()
    if not pending_uploads and not pending_steps:
        return
    
    # clear missing files, and make checkboxes for all pending files
    st.session_state['clashing_filenames'] = {}
    st.session_state['clashing_pids'] = {}
    for pending_upload in pending_uploads:
        if not Path(pending_upload).exists():
            exp.pending_uploads.remove(pending_upload)
            continue
    for pending_step in pending_steps:
        if not Path(pending_step).exists():
            del exp.pending_steps[pending_step]
            continue
    pending_uploads = {}
    if exp.pending_uploads:
        pending_uploads = exp.pending_uploads.copy()
    pending_steps = {}
    if exp.pending_steps:
        pending_steps = exp.pending_steps.copy()
    if not pending_uploads and not pending_steps:
        return
    
    warning_area.warning('The following files have been uploaded previously. '+
                        'Select the ones you want to overwrite and submit')

    with warning_area.form('warning_area_'+key, clear_on_submit=True):
        for pending_upload in pending_uploads:
            combined_pending.append(pending_upload)
            #affected_pids = trans.get_affected_pid_chain(pending_upload)
            #affected_fns = trans.get_affect_fn_chain(pending_upload)
            clashing_filenames, clashing_pids = trans.check_for_clashing_transactions(exp, filenames=pending_upload)
            st.checkbox(pending_upload + ' - ' + ','.join(clashing_filenames) + \
                    ' - ' + ','.join([util.unguard_pbc(pid, silent=True) for pid in clashing_pids]),
                          key=f'pending_file_checkbox_{pending_upload}')    
        for pending_step in pending_steps:
            combined_pending.append(pending_step)
            st.checkbox(pending_step, key=f'pending_file_checkbox_{pending_step}')
        st.form_submit_button('Submit', on_click=do_pending, args=[combined_pending])
    

def upload_echo_inputs(key):
    """
    Echo inputs are Nimbus outputs
    """
    exp = st.session_state['experiment']
    if exp.locked:
        st.warning(f'Experiment {exp.name} locked from further modification')
    else:
        title_area = st.container()
        title_area.write('')
        title = ''
        title_colour = '#f63366'
        info_area = st.container()
        _,upload_col,_ = st.columns([1,3,1])
        with upload_col:                
            with st.form("Echo input upload", clear_on_submit=True):  
                nim_outputs = st.file_uploader('Upload files: e.g. Echo_384_COC_0001_....csv', type='csv', 
                        accept_multiple_files=True, help='You can upload more than one file at once')
                submitted = st.form_submit_button("Upload files")

            if submitted and nim_outputs is not None:
                success = parse.upload(exp, nim_outputs, purpose='dna', overwrite=True)   

        nfs, efs, xbcs = exp.get_nimbus_filepaths()
        if not efs and not xbcs:
            title = "Load data inputs to enable Nimbus input file generation."
        if efs and not xbcs:
            title = 'All expected Echo inputs/Nimbus outputs now uploaded'
            title_colour = '#83b3c9'
            uploaded_nims = '</br>'.join([Path(ef).name for ef in efs])
            with info_area:
                stutil.custom_text('p', 'green', uploaded_nims)         
        if xbcs:
            title = 'Further Echo inputs/Nimbus outputs are expected'
            uploaded_nims = '</br>'.join([Path(ef).name for ef in efs])
            missing_nims = '</br>'.join(['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs])
            with info_area:
                stutil.custom_text('p', 'green', uploaded_nims)
                stutil.custom_text('p', 'black', 'The following Echo input files are expected (from the Nimbus):')
                stutil.custom_text('p', '#cf3276', missing_nims)

        with title_area:
            stutil.custom_text('h5', title_colour, title)
            stutil.add_vertical_space(1)
                        

def upload_pcr1_files(key):
    """
    Upload form for primer layout and volumes. Copies input files into the directory 
    and manages transactions
    """
    stutil.add_vertical_space(1)
    st.write('**Upload Primer Files (PCR round 1)**')
    exp = st.session_state['experiment']
    warning_area = st.container()
    st.session_state['upload_option'] = ''  # do we display pending files here
    primer_form = st.form('primer plate upload'+key, clear_on_submit=True)
    with primer_form:
        col1, col2 = st.columns(2)
        uploaded_primer_layouts = col1.file_uploader(
                'Upload Primer Plate Layouts - the barcode must be the first' +
                'part of the filename e.g. 123_primer_layout.csv', \
                key='primer_layout_uploader'+key,
                type='csv',
                accept_multiple_files=True) 
        uploaded_primer_volumes = col2.file_uploader(
                'Upload Primer Plate Volumes - the barcode must be the first part' +
                'of the filename e.g. 123_primer_volume.csv', \
                key='primer_vol_uploader'+key, 
                type='csv', 
                accept_multiple_files=True)

        upload_button = st.form_submit_button("Upload Files")

        if upload_button:
            st.session_state['upload_option'] = 'pcr1'
            if uploaded_primer_layouts:                          
                upl_pids = ''.join(upl.name.split('_')[0] for upl in uploaded_primer_layouts)
                success = parse.upload(exp, uploaded_primer_layouts, purpose='primer_layout')
                if success and not trans.is_pending(exp):
                    st.success(f'Added primer layouts for plates {",".join(upl_pids)}')
                elif not success:
                    st.error(f'Failed to write at least one primer layout, please see the log')
                     
            if uploaded_primer_volumes:
                upv_pids = ''.join(upv.name.split('_')[0] for upv in uploaded_primer_volumes)
                success = parse.upload(exp, uploaded_primer_volumes, purpose='primer_volume')
                if success and not trans.is_pending(exp):
                    st.success(f'Added primer volumes for plates {",".join(upv_pids)}')
                elif not success:
                    st.error(f'Failed to write at least one set of primer volumes, please see the log')
    
    #manage transactions:
    with warning_area:
        if trans.is_pending(exp) and st.session_state['upload_option'] == 'pcr1':
            pending_file_widget(key)
            st.session_state['upload_option'] = ''


def accept_amplicons(uploaded_amplicon_files, miseq_fn, stage3_fn):
    exp = st.session_state['experiment']
    warning_area = st.session_state['message_area']
    
    if Path(miseq_fn).exists():
        Path(miseq_fn).unlink()
    if Path(stage3_fn).exists():
        Path(stage3_fn).unlink()
    trans.enforce_file_consistency(exp)
    uaf_ids = ''.join(uaf.name for uaf in uploaded_amplicon_files)
    success = parse.upload(exp, uploaded_amplicon_files, purpose='amplicon')
    if success:
        warning_area.success(f'Added amplicon manifests from files {uaf_ids}')
    else:          
        warning_area.error(f'Failed to upload at least one amplicon manifest, please see the log')


def cancel_amplicons():
    stutil.add_vertical_space(1)
    st.write('**Upload Amplicon Files (PCR round 2)**')
    exp = st.session_state['experiment']
    warning_area = st.session_state['message_area']
    st.session_state['upload_option'] = ''  # do we display pending files here

                               
def upload_pcr2_files(key):
    """
    Upload inputs for indexing layout and volume. Extra option for amplicon plate upload.
    Copy uploaded files into the run folder and manage subsequent transactions
    """
    stutil.add_vertical_space(1)
    st.write('**Upload Index Files (PCR round 2)**')
    exp = st.session_state['experiment']
    warning_area = st.session_state['message_area']
    st.session_state['upload_option'] = ''  # do we display pending files here
    with st.form('index plate upload'+key, clear_on_submit=True):
        col1, col2 = st.columns(2)
        uploaded_index_layouts = col1.file_uploader(
                    'Upload i7i5 Index Plate Layout - the barcode must be'+
                    'the first part of the filename e.g. 123_index_layout.csv',
                    key='index_layout_uploader'+key, 
                    type='csv', 
                    accept_multiple_files=True)
 
        uploaded_index_volumes = col2.file_uploader(
                    'Upload i7i5 Index Plate Volumes - the barcode must be'+
                    'the first part of the filename e.g. 123_index_volume.csv', 
                    key='index_vol_uploader'+key, 
                    type='csv', 
                    accept_multiple_files=True)
 
        uploaded_amplicon_plates = col1.file_uploader(
                    'Upload Extra Amplicon Plates - CSV or XLSX.'+
                    'Invalidates MiSeq and Stage3 CSV files if they exist', 
                    key='amplicon_plate_uploader'+key, 
                    type=['csv', 'xlsx'], 
                    accept_multiple_files=True)
        
        upload_button = st.form_submit_button("Upload Files")

        if upload_button:
            st.session_state['upload_option'] = 'pcr2'
            if uploaded_index_layouts:
                uil_pids = ''.join(uil.name for uil in uploaded_index_layouts)
                success = parse.upload(exp, uploaded_index_layouts, purpose='index_layout')
                if success and not trans.is_pending(exp):
                    st.success(f'Added index layouts for plates {uil_pids}')
                elif not success:
                    st.error(f'Failed to write at least one index layout, please see the log')    

            if uploaded_index_volumes:
                uiv_pids = ''.join(uiv.name for uiv in uploaded_index_volumes)
                success = parse.upload(exp, uploaded_index_volumes, purpose='index_volume')
                if success and not trans.is_pending(exp):
                    st.success(f'Added index volumes for plates {uiv_pids}')
                elif not success:
                    st.error(f'Failed to write at least one set of index volumes, please see the log')   

            if uploaded_amplicon_plates:
                # Try to remove obsolete files
                miseq_fn = exp.get_exp_fn('MiSeq_'+exp.name+'.csv')
                if Path(miseq_fn).exists():
                    success = exp.del_file_record(miseq_fn)
                    if success:
                        exp.log(f'Info: removed obsolete Miseq.csv file')
                    else:
                        exp.log(f'Warning: could not remove obsolete Miseq.csv file')
                    
                stage3_fn = exp.get_exp_fn('Stage3.csv')
                if Path(stage3_fn).exists():
                    success = exp.del_file_record(stage3_fn)
                    if success:
                        exp.log(f'Info: removed obsolete Stage3.csv file')
                    else:
                        exp.log(f'Warning: could not remove obsolete Stage3.csv file')

                success = parse.upload(exp, uploaded_amplicon_plates, purpose='amplicon')
                uap_ids = [uap.name for uap in uploaded_amplicon_plates]
                if success:
                    warning_area.success(f'Added amplicon manifests from files {uap_ids}')
                else:
                    warning_area.write(f'Failed to upload at least one amplicon manifest, please see the log')
                
    
    # manage transactions
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'pcr2':
        pending_file_widget(key)
        st.session_state['upload_option'] = ''   
    stutil.add_vertical_space(1)


def upload_extra_consumables(key):
    """
    Uploads for extra files such as custom references, assay lists, and taq/water plates.
    Copy the uploaded files into the run folder and deal with subsequent transactions.
    """
    stutil.add_vertical_space(1)
    st.write('**Upload Custom Reference / Assay Lists**')
    exp = st.session_state['experiment']
    st.session_state['upload_option'] = ''  # do we display pending files here
    warning_area = st.container()
    #with st.form('Consumables upload'+key, clear_on_submit=True):
    upload_form = st.form('Consumables upload'+key, clear_on_submit=True)
    
    col1, col2 = upload_form.columns(2)
    uploaded_assaylists = col1.file_uploader('Upload Assay Lists', key='assaylist_uploader'+key, 
            type=['txt','csv'], accept_multiple_files=True)
    
    uploaded_references = col2.file_uploader('Upload Custom Reference Files', key='ref_uploader'+key, 
            type=['txt','fa','fasta'], accept_multiple_files=True)
                                                                            
    # maybe someday?
    #uploaded_taqwater_plates = col1.file_uploader('Upload Taq and Water Plates', key='taq_water_upload'+key, 
    #        type='csv', accept_multiple_files=True)
    upload_button = upload_form.form_submit_button("Upload Files")
    if upload_button:
        st.session_state['upload_option'] = 'consumables'
        if uploaded_references:
            success = parse.upload(exp, uploaded_references, 'reference_sequences')
            ref_names = [ur.name for ur in uploaded_references]
            if success and not trans.is_pending(exp):
                st.success(f'Added reference sequences from files {ref_names}')
            elif not success:
                st.error(f'Failed to upload at least one reference sequence file, please see the log')
        if uploaded_assaylists:
            #success = parse.upload(exp, uploaded_assaylists, 'primer_assay_map')
            # Genotyping team use the assay->primer direction
            success = parse.upload(exp, uploaded_assaylists, 'assay_primer_map')
            assaylist_names = ''.join(ual.name for ual in uploaded_assaylists)
            if success and not trans.is_pending(exp):
                st.success(f'Added assay/primer lists from files {assaylist_names}')
            elif not success:
                st.error(f'Failed to upload at least one assay/primer list file, please see the log')
        # TODO: add the ability to add a custom taq+water plate, perhaps?
        #if uploaded_taqwater_plates:
        #    success = exp.add_taqwater_layout_volume
   
    #manage transactions:
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'consumables':
        with warning_area:
            pending_file_widget(key)
        st.session_state['upload_option'] = ''
    stutil.add_vertical_space(1)


def custom_volumes(exp):
    volumes_dict = exp.transfer_volumes.copy()

    st.markdown('**Volumes for DNA and Primers**', help='To edit the volume, check the box.')
    col1, col2, col3,_ = st.columns([1,1,1,4])
    col1.write(f'Current volumes (nL):')
    allow_edit = col2.checkbox('Edit current volumes')
    reset_vols = col3.button("Reset to default")
    transfer_vol_df = pd.DataFrame.from_dict(data=exp.transfer_volumes, 
                                                orient='index', 
                                                columns=['Volume'])
    
    custom_vol_editor = st.data_editor(transfer_vol_df.T,
                                                    height=80,
                                                    use_container_width=True, 
                                                    key="volume_df", 
                                                    disabled=not allow_edit)
    tip_col, _ = st.columns(2)
    if allow_edit:
        tip_col.info('To change the volume, click or double click on the cell, type in the new value and hit Enter')
    
    if reset_vols:
        success = exp.add_custom_volumes({'DNA_VOL':util.DNA_VOL, 
                                          'PRIMER_VOL':util.PRIMER_VOL, 
                                          'PRIMER_TAQ_VOL':util.PRIMER_TAQ_VOL,
                                          'PRIMER_WATER_VOL':util.PRIMER_WATER_VOL, 
                                          'INDEX_VOL':util.INDEX_VOL, 
                                          'INDEX_TAQ_VOL':util.INDEX_TAQ_VOL,
                                          'INDEX_WATER_VOL':util.INDEX_WATER_VOL})
        if success:
            st.experimental_rerun()
        else:
            st.error('Volumes were not reset. Check the log for more information')

    #if edit is made to dataframe
    if st.session_state['volume_df']['edited_rows']:
        for vol in list(volumes_dict.keys()):
            volumes_dict[vol] = float(custom_vol_editor[vol]['Volume'])
        success = exp.add_custom_volumes(volumes_dict)
        if success:
            st.success("Added new volume")
            st.experimental_rerun()
        else:
            st.warning("Volume could not be added. Please enter in integers only.")

    stutil.add_vertical_space(1)


def upload_reference_sequences(key):
    """
    Upload form for custom references
    Args:
        key(str): Key for streamlit widget
    """
    exp = st.session_state['experiment']
    #TODO: check if they have already uploaded a reference file. 
    col, _ = st.columns([2,1])
    with col.form('References upload'+key, clear_on_submit=True):
        uploaded_references = st.file_uploader('Upload Custom Reference Files',
                                                key='ref_uploader'+key,
                                                type=['txt','fa','fasta'],    
                                                accept_multiple_files=True)
        
        upload_button = st.form_submit_button("Upload Files")

    if upload_button and uploaded_references:      
        success = parse.upload(exp, uploaded_references, 'reference_sequences')
        ur_names = [ur.name for ur in uploaded_references] 
        if success:
            msg = f'Succcess: added rodentity plate data from files {ur_names}'
            st.success(msg)
            exp.log(msg)
        else:
            msg = f'Failure: to upload at least one rodentity plate file, please see the log'
            st.error(msg)
            exp.log(msg)


def upload_rodentity_data(key):
    """
    Manage combining up to four rodentity plates (from JSON) with one destination PID 
    """
    exp = st.session_state['experiment']
    warning_area = st.session_state['message_area']
    st.session_state['upload_option'] = ''  # do we display pending files here
    plates_to_clear = [False, False, False, False]
    with st.expander('Add data from Rodentity JSON files',expanded=True):
        stutil.custom_text('p', '#7f8b8f', 'Upload up to 4 of your .json files and '+\
                           f'<span style="color: {FORM_BUTTON_COLOUR}">Submit </span>. Enter a desination plate ID '+\
                           f'and <span style="color: {FORM_BUTTON_COLOUR}">Accept </span> to assign them to a DNA plate. '+\
                            '<br> You have the option of clearing a plate ID '+\
                            'before you enter the desintation plate ID.', 'left')
        
        with st.form('rodentity_upload_form', clear_on_submit=True):
            epps_col1, _ = st.columns([1,1])
            rodentity_epps = epps_col1.file_uploader(
                        'Choose up to four Rodentity JSON files',
                        type='json',
                        accept_multiple_files=True)
            rod_upload_submit_button = st.form_submit_button('Submit')

        if rod_upload_submit_button and rodentity_epps:
            st.session_state['upload_option'] = 'rodentity'
            success = parse.upload(exp, rodentity_epps, 'rodentity_sample')
            rod_names = ', '.join(rod.name for rod in rodentity_epps)
            if success:
                st.success(f'Added rodentity plate data from files {rod_names}')
            else:
                st.error(f'Failed to upload at least one rodentity plate file, please see the log')
               
            if trans.is_pending(exp) and st.session_state['upload_option'] == 'rodentity':
                with warning_area:
                    pending_file_widget(key)
                st.session_state['upload_option'] = ''

            for rod_epp in rodentity_epps:
                rod_pid = util.guard_pbc(rod_epp.name.rstrip('.json'), silent=True)
                if all([exp.unassigned_plates[slot] != '' for slot in [1,2,3,4]]):
                    st.markdown(
                        '<p style="color:#FF0000">'+
                        f'Ran out of free slots for {util.unguard_pbc(rod_pid, silent=True)}</p>', 
                        unsafe_allow_html=True)
                    continue
                if rod_pid in exp.dest_sample_plates or \
                        rod_pid in exp.unassigned_plates[1] or rod_pid in exp.unassigned_plates[2] or \
                        rod_pid in exp.unassigned_plates[3] or rod_pid in exp.unassigned_plates[4]:
                    st.markdown('<p style="color:#FF0000">'+
                                'Warning: Rodentity plate barcode already used at least once</p>',
                                unsafe_allow_html=True)

                if rod_pid in exp.unassigned_plates['custom'] or \
                        (rod_pid in exp.plate_location_sample and 
                                exp.plate_location_sample[rod_pid]['purpose'] != 'sample'):
                    st.markdown('<p style="color:#FF0000">'+
                                'Error: Rodentity plate barcode already in use for a different purpose</p>', 
                                unsafe_allow_html=True)
                    continue                     
                for plate_key in [1,2,3,4]:
                    if exp.unassigned_plates[plate_key] not in exp.plate_location_sample:
                        exp.unassigned_plates[plate_key] = ''
                    if exp.unassigned_plates[plate_key] != '':
                        continue
                    exp.unassigned_plates[plate_key] = rod_pid
                    break
        
        # this is so ugly. Never use this form style again!
        for plate_key in [1,2,3,4]:
            if exp.unassigned_plates[plate_key] not in exp.plate_location_sample:
                exp.unassigned_plates[plate_key] = ''
        exp.save()
    
        rod_col1, _, rod_col2, _ = st.columns([3,1,2,2])
        with rod_col1.form('set_rod_plates_form', clear_on_submit=True):
            for i in range(4):
                plates_to_clear[i] = st.checkbox(
                            f"P{str(i+1)}: {util.unguard_pbc(exp.unassigned_plates[i+1], silent=True)}", 
                            help='Click the checkbox to allow a set plate ID to be cleared', 
                            key='check'+str(i+1))
            clear_plates_button = st.form_submit_button('Clear IDs', 
                    help='Clear selected Rodentity plate IDs')
            if clear_plates_button:
                for i, plate in enumerate(plates_to_clear):
                    if plate and exp.unassigned_plates[i+1]:
                        exp.unassigned_plates[i+1] = ''
                        plates_to_clear[i] = False
                exp.save()
                st.experimental_rerun()
                

        with rod_col2.form('rod_destination_form', clear_on_submit=True):
            rod_dp = st.text_input('Destination plate ID (barcode)', 
                            max_chars=30, 
                            key='rod_dp_key')
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
                    sample_plate_ids = [exp.unassigned_plates[k] for k in [1,2,3,4] if exp.unassigned_plates[k]]
                    success = exp.build_dna_plate_entry(sample_plate_ids, rod_dp, source='rodentity')
                    if not success:
                        st.markdown('<p style="color:#FF0000">' +
                                    'Failed to incorporate plate set. Please read the log.</p>', 
                                    unsafe_allow_html=True)
                        sleep(1.5)
                    else:
                        st.success('Added plate set')
                        exp.unassigned_plates = {1:'',2:'',3:'',4:''}
                        exp.save()
                        sleep(2)
                        st.experimental_rerun()
                        
                   
def upload_custom_manifests(key):
    """ Demonstration upload panel for custom manifests """
    exp = st.session_state['experiment']
    warning_area = st.session_state['message_area']
    st.session_state['upload_option'] = ''  # do we display pending files here
    if 'custom' not in exp.unassigned_plates or not exp.unassigned_plates['custom']:
        exp.unassigned_plates['custom'] = {'None':{}}
    with st.expander('Upload custom manifests', expanded=True):
        stutil.custom_text('p', '#7f8b8f', 'Upload the custom manifests and '+\
                          f'<span style="color: {FORM_BUTTON_COLOUR}">Submit </span>. Assign each file to a plate '+\
                           'using the selectboxes, then enter a destination plate ID and '+\
                          f'<span style="color: {FORM_BUTTON_COLOUR}">Accept </span>.', 'left')
        
        with st.form(key='manifest_upload_form', clear_on_submit=True):
            col1, _ = st.columns([1,1])
            with col1:
                manifest_uploads = st.file_uploader(label="Custom manifest upload", accept_multiple_files=True)
                #st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                submit_manifest = st.form_submit_button()
                if submit_manifest and manifest_uploads:
                    success = parse.upload(exp, manifest_uploads, 'custom_sample')
                    st.session_state['upload_option'] = 'custom'

        if trans.is_pending(exp) and st.session_state['upload_option'] == 'custom':
            with st.form(f'clash form {key}', clear_on_submit=True):
                pending_file_widget(key)
            st.session_state['upload_option'] = ''
       
        with st.form(key='selection_form', clear_on_submit=True):
            st.write('Select up to four 96-well sample plate IDs (barcodes) to combine into a 384-well DNA plate')
            col1, col2, _, col3, _ = st.columns([2,2,1,2,2])
            plate_options = [util.unguard_pbc(pid, silent=True) for pid in exp.plate_location_sample \
                    if exp.plate_location_sample[pid]['purpose'] == 'sample' \
                    and exp.plate_location_sample[pid]['source'] == 'custom']
            plate_options = [None] + plate_options
            #plate_options = [util.unguard_pbc(pid, silent=True) for pid in exp.unassigned_plates['custom'].keys()]
            with col1:
                p1 = util.guard_pbc(st.selectbox(label='Plate 1', options=plate_options, key='s1',), silent=True)
                p2 = util.guard_pbc(st.selectbox(label='Plate 2', options=plate_options, key='s2'), silent=True)
            with col2:
                p3 = util.guard_pbc(st.selectbox(label='Plate 3', options=plate_options, key='s3'), silent=True)
                p4 = util.guard_pbc(st.selectbox(label='Plate 4', options=plate_options, key='s4'), silent=True)
            with col3:
                dest_pid = st.text_input(label='Destination plate ID (barcode)', key='dest1')
                if dest_pid:
                    dest_pid = util.guard_pbc(dest_pid, silent=True)
                    
                submit_selection = st.form_submit_button("Accept")
                if submit_selection and dest_pid:
                    if dest_pid in exp.dest_sample_plates or dest_pid in exp.plate_location_sample or \
                            dest_pid in exp.unassigned_plates[1] or dest_pid in exp.unassigned_plates[2] or \
                            dest_pid in exp.unassigned_plates[3] or dest_pid in exp.unassigned_plates[4] or \
                            dest_pid in exp.unassigned_plates['custom']:
                        st.markdown('<p style="color:#FF0000">Error: Destination plate barcode already in use: ' +\
                                util.unguard_pbc(dest_pid, silent=True) + '</p>', unsafe_allow_html=True)
                    else:
                        st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                        sample_pids = [pid for pid in [p1,p2,p3,p4] 
                                if pid != util.guard_pbc('None', silent=True) and pid is not None]

                        success = exp.build_dna_plate_entry(sample_pids, dest_pid, source='custom')
                        if not success:
                            st.markdown('<p style="color:#FF0000">Failed to assign custom plates.'+\
                                    'Please see the log</p>', unsafe_allow_html=True)
                        if success:
                            st.success('Assigned custom plates')
    exp.save()


def add_barcodes(key, pcr_stage):
    """
    :param key: unique key for the input widgets
    """
    stutil.add_vertical_space(1)
    #st.write('**Add plate barcodes**')
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
                    st.success(f'Added PCR plate barcode {pcr_plate_barcode}')
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
                success = exp.add_standard_taqwater_plates([taqwater_plate_barcode], pcr_stage)
                if success:
                    st.success(f'Added taq+water plate barcode {taqwater_plate_barcode}')
                    sleep(1.5)
                    st.experimental_rerun()
                else:
                    st.write(f'Could not add taq+water plate barcode {taqwater_plate_barcode}, please see the log')
            else:
                st.write(f'This plate barcode {taqwater_plate_barcode} appears to already be in use')
    exp.save()


def upload_miseq_fastqs():
    """
    Widget for finding .fastq files in the user's local computer. 
    """
    exp = st.session_state['experiment']
    if not exp.locked:
        st.warning('Uploading sequence files will lock previous stages'+
                   'of the pipeline, preventing changes to plate layouts')
        
    st.subheader('Add Miseq FASTQ files to the experiment', 
                 help='Use the arrows and selectbox to navigate your local directory and '+\
                      'find the folder containing your FASTQ files.')
    
    fastq_path = dc.st_directory_picker("Select location of Miseq FASTQ files")
    fastq_files = [f for f in fastq_path.glob('*.fastq*')] + [f for f in fastq_path.glob('*.fq*')]
    if len(fastq_files) > 0:
        st.write(f"1. {str(fastq_files[0])}")
        if len(fastq_files) > 1:
            st.write(f"...")
            st.write(f"{len(fastq_files)}. {str(fastq_files[-1])}")

    import_fastqs = st.button('Import FASTQs', type='primary')
    file_field = st.empty()
    copy_progress = st.progress(0)
    if import_fastqs and len(fastq_files) > 0:
        for i,fp in enumerate(fastq_files):
            file_field.markdown(f'<hr5>Copying {fp} to experiment</h5>', unsafe_allow_html=True)
            copy_progress.progress(i/len(fastq_files))
            copy2(fp, exp.get_exp_dn('raw'))
        file_field.markdown('<h5>Done</h5>', unsafe_allow_html=True)
        copy_progress.progress(100)
        if not exp.locked:
            exp.lock()
            st.warning(f'Experiment {exp.name} is now locked from changes to plate layouts')

