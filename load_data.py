#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University

A web based interactive GUI with Streamlit. Plate and sample barcodes here are unguarded.
In all other code they must be guarded. We guard them here before we send them to any external function.
"""

import os
from pickle import TRUE
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

from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode

from stutil import custom_text, add_vertical_space, m, init_state, set_state, flip_state, mq
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
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021

Web application style interface using Streamlit. Present information as both a dashboard and a workflow in one page.

The GUI interacts with a single Experiment object at one time. Methods are called on this to activate pipeline
functionality. The Experiment then deals directly with the pipeline logic.

All callbacks must be labelled _cb. Because callbacks are executed first, they must create message
entries to be displayed, rather than writing the messages themselves
"""

HELP_COLOUR = '#7f8b8f'
FORM_BUTTON_COLOUR = '#4287f5'


def do_pending_cb(combined_pending, caller_id):
    """
    Action the pending_file_widget
    """
    exp = st.session_state['experiment']
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
            #m(f'Overwriting the previous version of {str(pending)} with the new file', 
            #        level='info', dest=('log','console','toast'), caller_id=caller_id)
            if pending in exp.pending_uploads:
                success = parse.accept_pending_upload(exp, pending)
            elif pending in exp.pending_steps:
                success = trans.accept_pending_transactions(exp, pending)
            if success:
                m(f'Success: overwrote existing file {pending}', level='success', dest=('log','console'), caller_id=caller_id)
                if clashing_filenames:
                    m(f'{", ".join(clashing_filenames)} potentially affected by change', 
                            level='warning', dest=('console','log'), caller_id=caller_id)
                if clashing_pids:
                    m(f'{", ".join(clashing_pids)} potentially affected by change', dest=('log','console'), 
                            level='warning', caller_id=caller_id)
            else:
                m(f'Failure: Could not overwrite existing file {pending}', level='failure', dest=('log','debug'),
                        caller_id=caller_id)
        else:
            # m(f'Keeping the old version of the file {pending}', level='info', dest=('log','console','toast'),
            #         caller_id=caller_id)
            success = exp.del_file_record(pending)
            if success:
                trans.clear_pending_transaction(exp, pending)
                m(f'Success: pending file or transaction {pending} cleared', level='failure', dest=('log'), caller_id=caller_id)
            else:
                m(f'Failure: pending file or transaction {pending} could not be removed', level='failure', dest=('log','debug'),
                    caller_id=caller_id)


def pending_file_widget(key, caller_id):
    """
    Allow user input via checkboxes to deal with clashes.
    Requires st.session_state['clashes'] from previously calling set_pending_file_clashes()
    args:
        key (str): a unique identifier for this widget
        caller_id (str): a unique identifier for the calling widget, for message passing
    """
    exp = st.session_state['experiment']
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
    
    m('The following files have been uploaded previously. Select the ones you want to overwrite and submit', 
            level = 'warning')
    
    with st.form('message_container_'+key, clear_on_submit=True):
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
        st.form_submit_button('Submit', on_click=do_pending_cb, args=[combined_pending, caller_id])
    

def upload_echo_inputs(key):
    """
    Echo inputs are Nimbus outputs
    """
    exp = st.session_state['experiment']
    caller_id = 'upload_echo_inputs'
    if exp.locked:
        m(f'Experiment {exp.name} locked from further modification', level='warning')
    else:
        title_area = st.container()
        title_area.write('')
        _, expect_col, _ = st.columns([1,3,1])
        _,upload_col,_ = st.columns([1,3,1])

        nfs, efs, xbcs = exp.get_nimbus_filepaths()

        if xbcs:
            with upload_col:                
                with st.form("Echo input upload", clear_on_submit=True):  
                    nim_outputs = st.file_uploader('Upload files: e.g. Echo_384_COC_0001_....csv', type='csv', 
                            accept_multiple_files=True, help='You can upload more than one file at once')
                    submitted = st.form_submit_button("Upload files")

                if submitted and nim_outputs is not None:
                    success = parse.upload(exp, nim_outputs, purpose='dna', overwrite=True)  
                    nim_names = ', '.join(nim.name for nim in nim_outputs)
                    if success:
                        m(f'{nim_names} successfully uploaded', level='info', dest=('log',))

        nfs, efs, xbcs = exp.get_nimbus_filepaths()

        title = ''
        title_colour = '#f63366'
        if not efs and not xbcs:
            title = "Load data inputs to enable Nimbus input file generation."
        if efs and not xbcs:
            title = 'All expected Echo inputs/Nimbus outputs now uploaded'
            title_colour = '#83b3c9'
            uploaded_nims = '</br>'.join([Path(ef).name for ef in efs])
            with expect_col:
                m(uploaded_nims, dest=('css',), size='p', color='green')         
        if xbcs:
            title = 'Echo inputs/Nimbus outputs are expected'
            title_colour = '#83b3c9'
            uploaded_nims = '</br>'.join([Path(ef).name for ef in efs])
            missing_nims = '</br>'.join(['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs])
            with expect_col:
                if len(efs) > 0:
                    m('Uploaded: ', dest=('css',))
                    m(uploaded_nims, dest=('css',), size='p', color='green')
                m('Upload the following Echo input files (from the Nimbus):', dest=('css',))
                m(missing_nims, dest=('css',), color='#cf3276')
        with title_area:
            m(title, dest=('css',), size='h5', color=title_colour)
            add_vertical_space(1)
            
    # manage file transactions
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'pcr1':
        pending_file_widget(key, caller_id)
        st.session_state['upload_option'] = ''
        
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
                       

def do_upload_primer_layout(upl, caller_id):
    """ helper method for clean deferred uploads of primer layouts """
    exp = st.session_state['experiment']
    upl_pid = upl.name.split('_')[0]
    success = parse.upload(exp, [upl], purpose='primer_layout')
    if success and not trans.is_pending(exp):
        m(f'Added primer layouts for plate {upl_pid} from file {upl.name}',
                level='success', dest=('log','console'), caller_id=caller_id)
    elif not success:
        m(f'Failed to read primer plate layout from {upl.name}, please see the log',
                level='error', dest=('debug',), caller_id=caller_id)


def do_upload_primer_volume(upv, caller_id):
    """ helper method for clean deferred uploads of primer volumes """
    exp = st.session_state['experiment']
    upv_pid = upv.name.split('_')[0]
    success = parse.upload(exp, [upv], purpose='primer_volume')
    if success and not trans.is_pending(exp):
        m(f'Added primer volumes for plate {upv_pid} from file {upv.name}',
                level='success', dest=('log','console'), caller_id=caller_id)
    elif not success:
        m(f'Failed to read primer plate volumes from {upv.name}, please see the log',
                level='error', dest=('debug',), caller_id=caller_id)


def upload_pcr1_files(key):
    """
    Upload form for primer layout and volumes. Copies input files into the directory 
    and manages transactions
    """  
    m('**Upload Primer Files (PCR round 1)**', dest=('mkdn',))
    m('Ensure the first prefix of the file name is the same for layouts and volumes (e.g. P1_layout.csv, P1_volume.csv)', 
      level = 'info')
    
    exp = st.session_state['experiment']
    caller_id = 'upload_pcr1_files'
    st.session_state['upload_option'] = ''  # do we display pending files here
    primer_form = st.form('primer plate upload'+key, clear_on_submit=True)
    init_state('deferred_primer_layout_upload_list', [])
    init_state('deferred_primer_volume_upload_list', [])
    with primer_form:
        col1, col2 = st.columns(2)
        uploaded_primer_layouts = col1.file_uploader(
                'Upload Primer Plate Layouts', \
                key='primer_layout_uploader'+key,
                type='csv',
                accept_multiple_files=True) 
        uploaded_primer_volumes = col2.file_uploader(
                'Upload Primer Plate Volumes', \
                key='primer_vol_uploader'+key, 
                type='csv', 
                accept_multiple_files=True)

        upload_button = st.form_submit_button("Upload Files")
        
        # check whether the files are potentially incorrect based on file names
        if upload_button:
            st.session_state['upload_option'] = 'pcr1'
            if uploaded_primer_layouts:
                for upl in uploaded_primer_layouts:
                    if 'index' in upl.name.lower() or 'vol' in upl.name.lower():
                        st.session_state['deferred_primer_layout_upload_list'].append(upl)
                    else:
                        do_upload_primer_layout(upl, caller_id)
                 
            if uploaded_primer_volumes:
                for upv in uploaded_primer_volumes:   
                    if 'index' in upv.name.lower() or 'layout' in upv.name.lower():
                        st.session_state['deferred_primer_volume_upload_list'].append(upv)
                    else:
                        do_upload_primer_volume(upv, caller_id)

    init_state('confirm_primer_upload_button', False)
    primer_confirmation_window = st.container()
    if st.session_state.get('confirm_primer_upload_button', False):
        for dl in st.session_state['deferred_primer_layout_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dl.name, False):
                do_upload_primer_layout(dl, caller_id)
            else:
                m(f'Skipping upload of {dl.name}')
        for dv in st.session_state['deferred_primer_volume_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dv.name, False):
                do_upload_primer_volume(dv, caller_id)
            else:
                m(f'Skipping upload of {dv.name}')
        st.session_state['deferred_primer_layout_upload_list'] = []
        st.session_state['deferred_primer_volume_upload_list'] = []
        st.session_state['confirm_primer_upload_button'] = False
        
    with primer_confirmation_window:
        make_button = False
        if st.session_state['deferred_primer_layout_upload_list']:
            make_button = True
            col1, col2 = st.columns(2)
            for dl in st.session_state['deferred_primer_layout_upload_list']:
                col1.checkbox(f'Upload {dl.name}', key=key+'checkbox_confirm_'+dl.name)
                if 'index' in dl.name.lower(): 
                    col2.write(f'Might be an index file, rather than primer layout')
                elif 'vol' in dl.name.lower():
                    col2.write(f'Might be a volume file, rather than primer layout')
                else:
                    col2.write(f'File flagged as potentially other purpose than primer layout')
        if st.session_state['deferred_primer_volume_upload_list']:
            make_button = True
            col1, col2 = st.columns(2)
            for dv in st.session_state['deferred_primer_volume_upload_list']:
                col1.checkbox(f'Upload {dv.name}', key=key+'checkbox_confirm_'+dv.name)
                if 'index' in dv.name.lower():
                    col2.write(f'Might be an index file, rather than primer volume')
                elif 'layout' in dv.name.lower():
                    col2.write(f'Might be a layout file, rather than primer volume')
                else:
                    col2.write(f'File flagged as potentially other purpose than primer volume')
        if make_button:
            st.button('Confirm', key='deferred_primer_list_confirm'+key, on_click=set_state, args=['confirm_primer_upload_button', True])
            
    #manage transactions:
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'pcr1':
        pending_file_widget(key, caller_id)
        st.session_state['upload_option'] = ''
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
    

def load_amplicons(key):
    """
    Upload only amplicon plates here
    """
    exp = st.session_state['experiment']
    caller_id = 'load_amplicons'
    add_vertical_space(2)
    m('**Upload Amplicon Plate Files**', dest=('mkdn',))
    m('You may add extra 384-well plates of pre-prepared amplicons. These will be incorporated during indexing (stage 4)', level='info')
    with st.form('index plate upload'+key, clear_on_submit=True): 
        uploaded_amplicon_plates = st.file_uploader(
                    'Upload Extra Amplicon Plates - CSV or XLSX.'+
                    'Invalidates MiSeq and Stage3 CSV files if they exist', 
                    key='amplicon_plate_uploader'+key, 
                    type=['csv', 'xlsx'], 
                    accept_multiple_files=True)
        
        upload_button = st.form_submit_button("Upload Files")
    
        if upload_button:
            st.session_state['upload_option'] = 'amplicons'
            if uploaded_amplicon_plates:
                # Try to remove obsolete files
                miseq_fn = exp.get_exp_fn('MiSeq_'+exp.name+'.csv')
                if Path(miseq_fn).exists():
                    success = exp.del_file_record(miseq_fn)
                    if success:
                        m(f'Removed existing file {miseq_fn}', level='info', dest=('log','console'))
                    else:
                        m(f'Could not remove obsolete file {miseq_fn}', level='warning', dest=('log','debug'))
                stage3_fn = exp.get_exp_fn('Stage3.csv')
                if Path(stage3_fn).exists():
                    success = exp.del_file_record(stage3_fn)
                    if success:
                        m(f'Removed obsolete file {stage3_fn}', level='info', dest=('log','console'))
                    else:
                        m(f'Could not remove obsolete file {stage3_fn}', level='warning', 
                                dest=('log','debug'))
                success = parse.upload(exp, uploaded_amplicon_plates, purpose='amplicon')
                uap_ids = ','.join(uap.name for uap in uploaded_amplicon_plates)
                if success:
                    if not trans.is_pending(exp):
                        m(f'Added amplicon manifests from {uap_ids}', level='success', 
                                dest=('log','console'))
                else:
                    m(f'Failed to upload at least one amplicon manifest, please see the log', 
                            level='error', dest=('debug',))
    # manage transactions
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'amplicons':
        pending_file_widget(key, caller_id)
        st.session_state['upload_option'] = ''           
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
        

def do_upload_index_layout(uil, caller_id):
    """ helper method for clean deferred uploads of index layouts """
    exp = st.session_state['experiment']
    uil_pid = uil.name.split('_')[0]
    success = parse.upload(exp, [uil], purpose='index_layout')
    if success and not trans.is_pending(exp):
        m(f'Added index layouts for plate {uil_pid} from file {uil.name}',
                level='success', dest=('log','console'), caller_id=caller_id)
    elif not success:
        m(f'Failed to read index plate layout from {uil.name}, please see the log',
                level='error', dest=('debug',), caller_id=caller_id)


def do_upload_index_volume(uiv, caller_id):
    """ helper method for clean deferred uploads of index volumes """
    exp = st.session_state['experiment']
    uiv_pid = uiv.name.split('_')[0]
    success = parse.upload(exp, [uiv], purpose='index_volume')
    if success and not trans.is_pending(exp):
        m(f'Added index volumes for plate {uiv_pid} from file {uiv.name}',
                level='success', dest=('log','console'), caller_id=caller_id)
    elif not success:
        m(f'Failed to read index plate volumes from {uiv.name}, please see the log',
                level='error', dest=('debug',), caller_id=caller_id)
        

def upload_pcr2_files(key):
    """
    Upload inputs for indexing layout and volume. Extra option for amplicon plate upload.
    Copy uploaded files into the run folder and manage subsequent transactions
    """
    m('**Upload Index Files (PCR round 2)**', dest=('mkdn',))
    exp = st.session_state['experiment']
    caller_id = 'upload_pcr2_files'
    st.session_state['upload_option'] = ''  # do we display pending files here
    init_state('deferred_index_layout_upload_list', [])
    init_state('deferred_index_volume_upload_list', [])
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
                for uil in uploaded_index_layouts:
                    if 'primer' in uil.name.lower() or 'vol' in uil.name.lower():
                        st.session_state['deferred_index_layout_upload_list'].append(uil)
                    else:
                        do_upload_index_layout(uil, caller_id)

            if uploaded_index_volumes:
                for uiv in uploaded_index_volumes:
                    if 'primer' in uiv.name.lower() or 'layout' in uiv.name.lower():
                        st.session_state['deferred_index_volume_upload_list'].append(uiv)
                    else:
                        do_upload_index_volume(uiv, caller_id)

            if uploaded_amplicon_plates:
                # Try to remove obsolete files
                miseq_fn = exp.get_exp_fn('MiSeq_'+exp.name+'.csv')
                if Path(miseq_fn).exists():
                    success = exp.del_file_record(miseq_fn)
                    if success:
                        m(f'Removed obsolete file {miseq_fn}', level='success', dest=('log','console'))
                    else:
                        m(f'Could not remove obsolete file {miseq_fn}', 
                                level='warning', dest=('log','debug'))
                    
                stage3_fn = exp.get_exp_fn('Stage3.csv')
                if Path(stage3_fn).exists():
                    success = exp.del_file_record(stage3_fn)
                    if success:
                        m(f'Removed obsolete file {stage3_fn}', level='success', dest=('log','console'))
                    else:
                        m(f'Could not remove obsolete file {stage3_fn}',
                                level='warning', dest=('log','debug'))

                success = parse.upload(exp, uploaded_amplicon_plates, purpose='amplicon')
                uap_ids = ','.join(uap.name for uap in uploaded_amplicon_plates)
                if success:
                    if not trans.is_pending(exp):
                        m(f'Added amplicon manifests from {uap_ids}',
                                level='success', dest=('log','console'))
                else:
                    m(f'Failed to upload at least one amplicon manifest, please see the log',
                            level='error', dest=('debug',))
      
    # handle deferred index file uploads
    init_state('confirm_index_upload_button', False)
    index_confirmation_window = st.container()
    if st.session_state.get('confirm_index_upload_button', False):
        for dl in st.session_state['deferred_index_layout_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dl.name, False):
                do_upload_index_layout(dl, caller_id)
            else:
                m(f'Skipping upload of {dl.name}')
        for dv in st.session_state['deferred_index_volume_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dv.name, False):
                do_upload_index_volume(dv, caller_id)
            else:
                m(f'Skipping upload of {dv.name}')
        st.session_state['deferred_index_layout_upload_list'] = []
        st.session_state['deferred_index_volume_upload_list'] = []
        st.session_state['confirm_index_upload_button'] = False
        
    with index_confirmation_window:
        make_button = False
        if st.session_state['deferred_index_layout_upload_list']:
            make_button = True
            col1, col2 = st.columns(2)
            for dl in st.session_state['deferred_index_layout_upload_list']:
                col1.checkbox(f'Upload {dl.name}', key=key+'checkbox_confirm_'+dl.name)
                if 'primer' in dl.name.lower(): 
                    col2.write(f'Might be a primer file, rather than index layout')
                elif 'vol' in dl.name.lower():
                    col2.write(f'Might be a volume file, rather than index layout')
                else:
                    col2.write(f'File flagged as potentially other purpose than index layout')
        if st.session_state['deferred_index_volume_upload_list']:
            make_button = True
            col1, col2 = st.columns(2)
            for dv in st.session_state['deferred_index_volume_upload_list']:
                col1.checkbox(f'Upload {dv.name}', key=key+'checkbox_confirm_'+dv.name)
                if 'primer' in dv.name.lower():
                    col2.write(f'Might be a primer file, rather than index volume')
                elif 'layout' in dv.name.lower():
                    col2.write(f'Might be a layout file, rather than index volume')
                else:
                    col2.write(f'File flagged as potentially other purpose than primer volume')
        if make_button:
            st.button('Confirm', key='deferred_index_list_confirm'+key, on_click=set_state, args=['confirm_index_upload_button', True])

    # manage transactions
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'pcr2':
        #with st.session_state['message_container']:
        pending_file_widget(key, caller_id)
        st.session_state['upload_option'] = ''

    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
        

def upload_assaylist(key):
    """
    Upload just the assay list (assay-primer mapping) file. Return True on success, False on failure
    """
    m('**Upload Assay List file (assay-primer mapping)**', dest=('mkdn',))
    exp = st.session_state['experiment']
    caller_id = ('upload_assaylist')
    upload_form = st.form('Consumables upload'+key, clear_on_submit=True)
    
    col1, col2 = upload_form.columns(2)
    uploaded_assaylists = col1.file_uploader('Upload Assay Lists', key='assaylist_uploader'+key, 
            type=['txt','csv'], accept_multiple_files=True)

    upload_button = upload_form.form_submit_button("Upload Files")
    success = False
    if upload_button:
        if uploaded_assaylists:
            success = parse.upload(exp, uploaded_assaylists, 'assay_primer_map')
            assaylist_names = ''.join(ual.name for ual in uploaded_assaylists)
            if success and not trans.is_pending(exp):
                m(f'Added assay/primer lists from files {assaylist_names}',
                        level='success', dest=('log','console'))
            elif not success:
                m(f'Failed to upload at least one assay/primer list file, please see the log',
                        level='error', dest=('debug',))
   
    #manage transactions:
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'consumables':
        #with st.session_state['message_container']:
        pending_file_widget(key, caller_id)
        st.session_state['upload_option'] = ''
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
    return success
    

def upload_extra_consumables(key):
    """
    Uploads for extra files such as custom references, assay lists, and taq/water plates.
    Copy the uploaded files into the run folder and deal with subsequent transactions.
    """
    m('**Upload Custom Reference / Assay Lists**', dest=('mkdn',))
    exp = st.session_state['experiment']
    caller_id = 'upload_extra_consumables'
    st.session_state['upload_option'] = ''  # do we display pending files here
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
            ref_names = ', '.join(ur.name for ur in uploaded_references)
            if success and not trans.is_pending(exp):
                m(f'Added reference sequences from files {ref_names}', 
                        level='success', dest=('log','console'))
            elif not success:
                m(f'Failed to upload at least one reference sequence file, please see the log', 
                        level='error', dest=('debug',))
        if uploaded_assaylists:
            #success = parse.upload(exp, uploaded_assaylists, 'primer_assay_map')
            # Genotyping team use the assay->primer direction
            success = parse.upload(exp, uploaded_assaylists, 'assay_primer_map')
            assaylist_names = ''.join(ual.name for ual in uploaded_assaylists)
            if success and not trans.is_pending(exp):
                m(f'Added assay/primer lists from files {assaylist_names}',
                        level='success', dest=('log','console'))
            elif not success:
                m(f'Failed to upload at least one assay/primer list file, please see the log',
                        level='error', dest=('debug',))
        # TODO: add the ability to add a custom taq+water plate, perhaps?
        #if uploaded_taqwater_plates:
        #    success = exp.add_taqwater_layout_volume
   
    #manage transactions:
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'consumables':
        #with st.session_state['message_container']:
        pending_file_widget(key, caller_id)
        st.session_state['upload_option'] = ''
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def check_assay_file(exp):
    return any(file.get('purpose', '') == 'assay_primer_map' for file in exp.uploaded_files.values())


def custom_volumes(exp):
    volumes_dict = exp.transfer_volumes.copy()
    caller_id = 'custom_volumes'
    m('**Volumes for DNA and Primers**', dest=('mkdn',))
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
    if allow_edit:
        m('To change the volume, click or double click on the cell, type in the new value and hit Enter', 
                level='info')
    
    if reset_vols:
        success = exp.add_custom_volumes({'DNA_VOL':util.DNA_VOL, 
                                          'PRIMER_VOL':util.PRIMER_VOL, 
                                          'PRIMER_TAQ_VOL':util.PRIMER_TAQ_VOL,
                                          'PRIMER_WATER_VOL':util.PRIMER_WATER_VOL, 
                                          'INDEX_VOL':util.INDEX_VOL, 
                                          'INDEX_TAQ_VOL':util.INDEX_TAQ_VOL,
                                          'INDEX_WATER_VOL':util.INDEX_WATER_VOL}, caller_id)
        if success:
            m('Transfer volumes reset', level='info', dest=('log',))
            sleep(1.5)
            st.rerun()
        else:
            m('Transfer volumes were not reset. Check the log for more information', 
                    level='error', dest=('debug',))

    #if edit is made to dataframe
    if st.session_state['volume_df']['edited_rows']:
        for vol in list(volumes_dict.keys()):
            volumes_dict[vol] = float(custom_vol_editor[vol]['Volume'])
        success = exp.add_custom_volumes(volumes_dict, caller_id)
        if success:
            m(f'Modified transfer volumes {volumes_dict}', level='success', dest=('log','console'))
            sleep(1.5)
            st.rerun()
        else:
            m("Volume could not be modified. Please use integers only", level='error', dest=('debug',))

    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def upload_reference_sequences(key):
    """
    Upload custom references
    """
    exp = st.session_state['experiment']
    caller_id = 'upload_reference_sequences'
    with st.form('References upload'+key, clear_on_submit=True):
        uploaded_references = st.file_uploader(
                'Upload Custom Reference Files',
                key='ref_uploader'+key,
                type=['txt','fa','fasta'],    
                accept_multiple_files=True)
        upload_button = st.form_submit_button("Upload Files")

    if upload_button and uploaded_references:      
        success = parse.upload(exp, uploaded_references, 'reference_sequences')
        ur_names = [ur.name for ur in uploaded_references] 
        if success:
            if not trans.is_pending(exp):
                m(f'Added rodentity plate data from files {ur_names}', 
                        level='success', dest=('log','console'))
        else:
            m(f'Could not upload at least one rodentity plate file, please see the log',
                    level='error', dest=('debug',))
    #manage transactions:
    if trans.is_pending(exp) and st.session_state['upload_option'] == 'consumables':
        #with st.session_state['message_container']:
        pending_file_widget(key)
        st.session_state['upload_option'] = ''
        
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
    

def load_rodentity_data(key):
    """
    Manage combining up to four rodentity plates (from JSON) with one destination PID 
    """
    exp = st.session_state['experiment']
    caller_id = 'load_rodentity_data'
    st.session_state['upload_option'] = ''  # do we display pending files here
    plates_to_clear = [False, False, False, False]
    with st.expander('Add data from Rodentity JSON files',expanded=True):
        m('Add up to four Rodentity JSON ear-punch plate (96-well) files at a time and click on **Submit**. '+\
                'This will move the plate identifiers into the four slots below (P1-4). ', level='info')
        
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
                m(f'Added rodentity plate data from files {rod_names}', 
                        level='success', dest=('log','console'))
            else:
                m(f'Failed to upload at least one rodentity plate file, please see the log',
                        level='error', dest=('debug',))
               
            if trans.is_pending(exp) and st.session_state['upload_option'] == 'rodentity':
                #with message_container:
                pending_file_widget(key, caller_id)
                st.session_state['upload_option'] = ''

            for rod_epp in rodentity_epps:
                rod_pid = util.guard_pbc(rod_epp.name.rstrip('.json'), silent=True)
                if all([exp.unassigned_plates[slot] != '' for slot in [1,2,3,4]]):
                    m(f'Ran out of free slots for {util.unguard_pbc(rod_pid, silent=True)}',
                            level='warning', dest=('log',))
                    continue
                if rod_pid in exp.dest_sample_plates or \
                        rod_pid in exp.unassigned_plates[1] or rod_pid in exp.unassigned_plates[2] or \
                        rod_pid in exp.unassigned_plates[3] or rod_pid in exp.unassigned_plates[4]:
                    m('Warning: Rodentity plate barcode already used at least once'
                            , level='warning', dest=('log',))

                if rod_pid in exp.unassigned_plates['custom'] or \
                        (rod_pid in exp.plate_location_sample and  
                        exp.plate_location_sample[rod_pid]['purpose'] != 'sample'):
                    m('Error: Rodentity plate barcode already in use for a different purpose', 
                            level='warning', dest=('log',))
                    continue                     
                for plate_key in [1,2,3,4]:
                    if exp.unassigned_plates[plate_key] not in exp.plate_location_sample:
                        exp.unassigned_plates[plate_key] = ''
                    if exp.unassigned_plates[plate_key] != '':
                        continue
                    exp.unassigned_plates[plate_key] = rod_pid
                    break
        m('Then assign a DNA plate (384-well) barcode identifier for plate combining with the '+\
                'Hamilton Nimbus, and click on **Accept**',level='info')
        # this is so ugly. Never use this form style again!
        for plate_key in [1,2,3,4]:
            if exp.unassigned_plates[plate_key] not in exp.plate_location_sample:
                exp.unassigned_plates[plate_key] = ''
    
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
                st.rerun()
                

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
                    m(f'Destination plate barcode already in use: {util.unguard_pbc(rod_dp, silent=True)}',
                            level='warning', dest=('log',))
                    sleep(1.5)
                else:
                    sample_plate_ids = [exp.unassigned_plates[k] for k in [1,2,3,4] if exp.unassigned_plates[k]]
                    success = exp.build_dna_plate_entry(sample_plate_ids, rod_dp, source='rodentity')
                    if not success:
                        m('Failed to incorporate plate set. Please read the log', 
                                level='error', dest=('debug',))
                        sleep(1.5)
                    else:
                        m(f'Added plate set for {util.unguard_pbc(rod_dp, silent=True)}',
                               level='success', dest=('log','console'))
                        exp.unassigned_plates = {1:'',2:'',3:'',4:''}
                        sleep(2)
                        #st.rerun()
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
                        
                   
def load_custom_manifests(key):
    """ Demonstration upload panel for custom manifests """
    exp = st.session_state['experiment']
    caller_id = 'load_custom_manifest'
    st.session_state['upload_option'] = ''  # do we display pending files here
    if 'custom' not in exp.unassigned_plates or not exp.unassigned_plates['custom']:
        exp.unassigned_plates['custom'] = {'None':{}}
    with st.expander('Upload custom manifests', expanded=True):
        m('For non-Rodentity samples in 96-well plates, you can add these here. '+\
                'Upload them first and then choose the layout you want.', level='info')
        with st.form(key='manifest_upload_form', clear_on_submit=True):
            col1, _ = st.columns([1,1])
            with col1:
                manifest_uploads = st.file_uploader(label="Custom manifest upload", accept_multiple_files=True)
                #st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                submit_manifest = st.form_submit_button()
                if submit_manifest and manifest_uploads:
                    success = parse.upload(exp, manifest_uploads, 'custom_sample')
                    if success:
                        m(f'Custom file uploaded', level='success')
                    else:
                        m('Custom file upload failed, see log', level='error', dest=('debug',))
                    st.session_state['upload_option'] = 'custom'

        if trans.is_pending(exp) and st.session_state['upload_option'] == 'custom':
            #with st.form(f'clash form {key}', clear_on_submit=True):
            pending_file_widget(key, caller_id)
            st.session_state['upload_option'] = ''
        m('Choose the desired layout for the 384-well DNA plate (made with the Hamilton Nimbus). '+\
                'Then provide the DNA plate name and press **Accept**', level='info')
        
        with st.form(key='selection_form', clear_on_submit=True):
            m('Select up to four 96-well sample plate IDs (barcodes) to combine into a 384-well DNA plate', level='normal')
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
                        m(f'Destination plate barcode already in use: {util.unguard_pbc(dest_pid, silent=True)}',
                                level='warning', dest=('log','debug'))
                    else:
                        #st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                        sample_pids = [pid for pid in [p1,p2,p3,p4] 
                                if pid != util.guard_pbc('None', silent=True) and pid is not None]
                        success = exp.build_dna_plate_entry(sample_pids, dest_pid, source='custom')
                        if success:
                            m(f'Assigned custom plates to {util.unguard_pbc(dest_pid)}', level='success', dest=('log','console'))
                        else:
                            m('Failed to assign custom plates', level='failure', dest=('log','debug'))
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
        

def add_pcr_barcodes(key, dna_pids):
    caller_id = 'add_pcr_barcodes'
    exp = st.session_state['experiment']

    with st.form('Add PCR plate barcode', clear_on_submit=True):
        pcr_plate_barcode = st.text_input('PCR Plate Barcode', \
                    placeholder='Type in the barcode for PCR plate', key='pcr_plate_barcode' + key)
        
        pcr_plate_entered = st.form_submit_button('Add')

        if pcr_plate_entered and pcr_plate_barcode: 
            guarded_pcr_plate_barcode = util.guard_pbc(pcr_plate_barcode, silent=True)
            if guarded_pcr_plate_barcode not in exp.plate_location_sample:
                success = exp.add_pcr_plates([guarded_pcr_plate_barcode])
                if success:
                    m(f'Added PCR plate barcode {pcr_plate_barcode}', level='success', 
                            dest=('log','console'), caller_id=caller_id)
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add PCR plate barcode {pcr_plate_barcode}, please see the log',
                            level='error', dest=('debug','log'), caller_id=caller_id)
            else:
                m(f'This plate barcode {pcr_plate_barcode} appears to already be in use', 
                        level='error', dest=('debug','log'), caller_id=caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def add_taqwater_barcodes(key, pcr_stage):
    exp = st.session_state['experiment']
    caller_id = 'add_taqwater_barcodes'
    with st.form('Add taq+water plate barcode', clear_on_submit=True):
        taqwater_plate_barcode = st.text_input('Taq and Water Plate Barcode', \
                placeholder='Type in the barcode for taq+water plate', key='taq_water_barcode' + key)
        taqwater_plate_entered = st.form_submit_button('Add')

        if taqwater_plate_entered and taqwater_plate_barcode:
            guarded_taqwater_plate_barcode = util.guard_pbc(taqwater_plate_barcode, silent=True)
            if guarded_taqwater_plate_barcode not in exp.plate_location_sample:
                success = exp.add_standard_taqwater_plates([taqwater_plate_barcode], pcr_stage)
                if success:
                    m(f'Added taq+water plate {taqwater_plate_barcode}', level='success', 
                            dest=('log','console'), caller_id=caller_id)
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add taq+water plate barcode {taqwater_plate_barcode}, please see the log',
                            level='error', dest=('debug',), caller_id=caller_id)
            else:
                m(f'This plate barcode {taqwater_plate_barcode} appears to already be in use',
                        level='error', dest=('log','debug'), caller_id=caller_id)
                
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def provide_barcodes(key, pcr_stage, dna_pids):
    """
    :param key: unique key for the input widgets
    """
    caller_id = 'provide_barcodes'
    #add_vertical_space(2)
    m('**Add plate barcodes**', dest=('mkdn',))
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
                    m(f'Added PCR plate barcode {pcr_plate_barcode}', level='success', 
                            dest=('log','console'), caller_id=caller_id)
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add PCR plate barcode {pcr_plate_barcode}, please see the log',
                            level='error', dest=('debug','log'), caller_id=caller_id)
            else:
                m(f'This plate barcode {pcr_plate_barcode} appears to already be in use', 
                        level='error', dest=('debug','log'), caller_id=caller_id)

    with taqwater_col.form('Add taq+water plate barcode', clear_on_submit=True):
        taqwater_plate_barcode = st.text_input('Taq and Water Plate Barcode', \
                placeholder='Type in the barcode for taq+water plate', key='taq_water_barcode' + key)
        taqwater_plate_entered = st.form_submit_button('Add')

        if taqwater_plate_entered and taqwater_plate_barcode:
            guarded_taqwater_plate_barcode = util.guard_pbc(taqwater_plate_barcode, silent=True)
            if guarded_taqwater_plate_barcode not in exp.plate_location_sample:
                success = exp.add_standard_taqwater_plates([taqwater_plate_barcode], pcr_stage)
                if success:
                    m(f'Added taq+water plate {taqwater_plate_barcode}', level='success', 
                            dest=('log','console'), caller_id=caller_id)
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add taq+water plate barcode {taqwater_plate_barcode}, please see the log',
                            level='error', dest=('debug',), caller_id=caller_id)
            else:
                m(f'This plate barcode {taqwater_plate_barcode} appears to already be in use',
                        level='error', dest=('log','debug'), caller_id=caller_id)
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def upload_miseq_fastqs():
    exp = st.session_state['experiment']
    caller_id = 'upload_miseq_fastqs'
    if not exp.locked:
        pass
        #st.warning('Uploading sequence files will lock previous stages '+
        #        'of the pipeline, preventing changes to plate layouts')
    st.markdown('<h4 style="color:#000000">Add Miseq FASTQ files to experiment</h4>', 
            unsafe_allow_html=True)
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
            copy2(fp, exp.get_exp_dn('raw'))
        file_field.markdown('<h5>Done</h5>', unsafe_allow_html=True)
        copy_progress.progress(100)
        #if not exp.locked:
        #    exp.lock()
        #    m(f'Experiment {exp.name} is now locked from changes to plate layouts', level='info', dest=('log',))
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []



