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
from io import StringIO, BytesIO
from shutil import copy2, copyfileobj
from time import sleep
from unittest import expectedFailure
import uuid

import pandas as pd
from pandas.core.arrays.period import raise_on_incompatible

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

def queue_upload(streams, purpose, caller_id=None, overwrite_plates=True):
    """
    All file stream uploads go here first and are queued for parsing  
    - Each file is uploaded with a pending_ prefix and
    - saved to /exp/uploads/
    - pending file uploads are saved to exp['uploaded_files']['_upload_queue'] =\
            (pending_fn, purpose, caller_id, overwrite_plates)
    The queue is processed by parse.process_pending_uploads(exp)
    """
    exp = st.session_state['experiment']
    if exp.locked and purpose not in {'primer_assay_map','reference_sequences'}:
        m(f'Cannot add files other than assay file or references while lock is active', level='failure', caller_id=caller_id)
        return False
    
    for stream in streams:
        pfp = exp.get_exp_fn(stream.name, subdir='uploads', trans=True)
        if Path(pfp).exists():
            success = util.delete_file(pfp, caller_id=caller_id)
            if not success:
                continue
            
        m(f"uploading {stream.name} to {pfp}", level='info', dest=('noGUI',))
        with open(pfp, 'wb') as outf:
            try:
                copyfileobj(BytesIO(stream.getvalue()), outf)
            except Exception as exc:
                m(f'Could not copy stream {stream.name} to file {pfp}, {exc}', level='error', caller_id=caller_id)
                continue

        # queue the upload for later processing
        if '_upload_queue' not in exp.uploaded_files:
            exp.uploaded_files['_upload_queue'] = {}
        
        exp.uploaded_files['_upload_queue'][pfp] = (purpose, caller_id, overwrite_plates)
        m(f'{pfp} queued for parsing', level='info', dest=('noGUI',))


def do_pending_cb(pending_cbs, caller_id):
    """
    Action the pending_file_widget
    """
    exp = st.session_state['experiment']
    for pcb in pending_cbs:
        pfp = pending_cbs[pcb]
        if pcb in st.session_state and st.session_state[pcb]:
            clashing_files, clashing_pids = trans.check_for_clashing_transactions(exp, filenames=pfp)
            if pfp in exp.uploaded_files['_upload_pending']:
                purpose, caller_id, overwrite_plates = exp.uploaded_files['_upload_pending'][pfp]
                success = parse.accept_pending_upload(exp, pfp, purpose, caller_id=caller_id, overwrite_plates=overwrite_plates)
            elif pfp in exp.pending_steps:
                success = trans.accept_pending_transactions(exp, pfp, caller_id=caller_id)
            if success:
                if clashing_files:
                    m(f'{", ".join(clashing_files)} potentially affected by change', 
                            level='warning', dest=('noGUI',))
                if clashing_pids:
                    m(f'{", ".join(clashing_pids)} potentially affected by change', 
                            level='warning', dest=('noGUI',))
            else:
                m(f'could not overwrite existing file {pfp}', level='error',
                        dest=('noGUI',))
            if pfp in exp.uploaded_files['_upload_pending']:
                del exp.uploaded_files['_upload_pending'][pfp]


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
    
    pending_uploads = list(exp.uploaded_files['_upload_pending'].keys())
    pending_steps = {}
    if exp.pending_steps:
        pending_steps = exp.pending_steps.copy()
    if not pending_uploads and not pending_steps:
        return
    
    # clear missing files, and make checkboxes for all pending files
    st.session_state['clashing_filenames'] = {}
    st.session_state['clashing_pids'] = {}
    for pfp in pending_uploads:
        if not Path(pfp).exists():
            del exp.uploaded_files['_upload_pending'][pfp]
            continue
    for pending_step in pending_steps:
        if not Path(pending_step).exists():
            del exp.pending_steps[pending_step]
            continue
    pending_uploads = list(exp.uploaded_files['_upload_pending'].keys())
    pending_steps = {}
    if exp.pending_steps:
        pending_steps = exp.pending_steps.copy()
    if not pending_uploads and not pending_steps:
        return
    
    st.warning('The following files have been uploaded previously. Select the ones you want to overwrite and submit')
    
    with st.form('message_container_'+key, clear_on_submit=True):
        pending_cbs = {}
        for pfp in pending_uploads:
            key = f'pending_file_checkbox_{uuid.uuid4()}'
            pending_cbs[key] = pfp
            #clashing_filenames, clashing_pids = trans.check_for_clashing_transactions(exp, filenames=pfp)
            #st.checkbox(pfp + ' - ' + ','.join(clashing_filenames) + \
            #        ' - ' + ','.join([util.unguard_pbc(pid, silent=True) for pid in clashing_pids]),
            #              key=f'pending_file_checkbox_{pfp}')
            st.checkbox(pfp, key=key)    
        for pending_step in pending_steps:
            key = f'pending_file_checkbox_{uuid.uuid4()}'
            pending_cbs[key] = pending_step
            st.checkbox(pending_step, key=key)
        st.form_submit_button('Submit', on_click=do_pending_cb, args=[pending_cbs, caller_id])
    

def upload_echo_inputs(key):
    """
    Echo inputs are Nimbus outputs
    """
    exp = st.session_state['experiment']
    caller_id = 'upload_echo_inputs'
    if exp.locked:
        st.warning(f'Experiment {exp.name} locked from further modification')
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
                    queue_upload(nim_outputs, purpose='dna', caller_id=caller_id, overwrite_plates=True)  

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
                m(uploaded_nims, level='display', dest=('css',), size='p', color='green')         
        if xbcs:
            title = 'Echo inputs/Nimbus outputs are expected'
            title_colour = '#83b3c9'
            uploaded_nims = '</br>'.join([Path(ef).name for ef in efs])
            missing_nims = '</br>'.join(['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs])
            with expect_col:
                if len(efs) > 0:
                    m('Uploaded: ', level='display', dest=('css',))
                    m(uploaded_nims, level='display', dest=('css',), size='p', color='green')
                m('Upload the following Echo input files (from the Nimbus):', level='display', dest=('css',))
                m(missing_nims, level='display', dest=('css',), color='#cf3276')
        with title_area:
            m(title, level='display', dest=('css',), size='h5', color=title_colour)
            add_vertical_space(1)

    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
        
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
                       

def do_upload_primer_layout(upl, caller_id):
    """ helper method for clean deferred uploads of primer layouts """
    queue_upload([upl], 'primer_layout', caller_id=caller_id)


def do_upload_primer_volume(upv, caller_id):
    """ helper method for clean deferred uploads of primer volumes """
    queue_upload([upv], 'primer_volume', caller_id=caller_id)


def upload_pcr1_files(key):
    """
    Upload form for primer layout and volumes. Copies input files into the directory 
    and manages transactions
    """  
    st.markdown('**Upload Primer Files (PCR round 1)**')
    st.info('Ensure the first prefix of the file name is the same for layouts and volumes (e.g. P1_layout.csv, P1_volume.csv)')
    
    exp = st.session_state['experiment']
    caller_id = 'upload_pcr1_files'
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
                m(f'Skipping upload of {dl.name}', level='info')
        for dv in st.session_state['deferred_primer_volume_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dv.name, False):
                do_upload_primer_volume(dv, caller_id)
            else:
                m(f'Skipping upload of {dv.name}', level='info')
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
     
    parse.process_upload_queue(exp)
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
    

def load_amplicons(key):
    """
    Upload only amplicon plates here
    """
    exp = st.session_state['experiment']
    caller_id = 'load_amplicons'
    add_vertical_space(2)
    st.markdown('**Upload Amplicon Plate Files**')
    st.info('You may add extra 384-well plates of pre-prepared amplicons. These will be incorporated during indexing (stage 4)')
    with st.form('index plate upload'+key, clear_on_submit=True): 
        uploaded_amplicon_plates = st.file_uploader(
                    'Upload Extra Amplicon Plates - CSV or XLSX.'+
                    'Invalidates MiSeq and Stage3 CSV files if they exist', 
                    key='amplicon_plate_uploader'+key, 
                    type=['csv', 'xlsx'], 
                    accept_multiple_files=True)
        
        upload_button = st.form_submit_button("Upload Files")
    
        if upload_button:
            if uploaded_amplicon_plates:
                # Try to remove obsolete files
                miseq_fn = exp.get_exp_fn('MiSeq_'+exp.name+'.csv')
                if Path(miseq_fn).exists():
                    success = exp.del_file_record(miseq_fn)
                    if success:
                        m(f'Removed existing file {miseq_fn}', level='info')
                    else:
                        m(f'Could not remove obsolete file {miseq_fn}', level='warning')
                stage3_fn = exp.get_exp_fn('Stage3.csv')
                if Path(stage3_fn).exists():
                    success = exp.del_file_record(stage3_fn)
                    if success:
                        m(f'Removed obsolete file {stage3_fn}', level='info')
                    else:
                        m(f'Could not remove obsolete file {stage3_fn}', level='warning')
                queue_upload(uploaded_amplicon_plates, 'amplicon', caller_id=caller_id)
               
    parse.process_upload_queue(exp)
    # manage transactions
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)          
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
        

def do_upload_index_layout(uil, caller_id):
    """ helper method for clean deferred uploads of index layouts """
    queue_upload([uil], 'index_layout', caller_id=caller_id)
    

def do_upload_index_volume(uiv, caller_id):
    """ helper method for clean deferred uploads of index volumes """
    queue_upload([uiv], 'index_volume', caller_id=caller_id)
            

def upload_pcr2_files(key):
    """
    Upload inputs for indexing layout and volume. Extra option for amplicon plate upload.
    Copy uploaded files into the run folder and manage subsequent transactions
    """
    st.markdown('**Upload Index Files (PCR round 2)**')
    exp = st.session_state['experiment']
    caller_id = 'upload_pcr2_files'
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
                        m(f'Removed obsolete file {miseq_fn}', level='success')
                    else:
                        m(f'Could not remove obsolete file {miseq_fn}', 
                                level='warning')
                    
                stage3_fn = exp.get_exp_fn('Stage3.csv')
                if Path(stage3_fn).exists():
                    success = exp.del_file_record(stage3_fn)
                    if success:
                        m(f'Removed obsolete file {stage3_fn}', level='success')
                    else:
                        m(f'Could not remove obsolete file {stage3_fn}',
                                level='warning')
                queue_upload(uploaded_amplicon_plates, purpose='amplicon', caller_id=caller_id)              
      
    # handle deferred index file uploads
    init_state('confirm_index_upload_button', False)
    index_confirmation_window = st.container()
    if st.session_state.get('confirm_index_upload_button', False):
        for dl in st.session_state['deferred_index_layout_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dl.name, False):
                do_upload_index_layout(dl, caller_id)
            else:
                m(f'Skipping upload of {dl.name}', level='warning')
        for dv in st.session_state['deferred_index_volume_upload_list']:
            if st.session_state.get(key+'checkbox_confirm_'+dv.name, False):
                do_upload_index_volume(dv, caller_id)
            else:
                m(f'Skipping upload of {dv.name}', level='warning')
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

    parse.process_upload_queue(exp)
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)

    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
        

def upload_assaylist(key):
    """
    Upload just the assay list (assay-primer mapping) file. Return True on success, False on failure
    """
    st.markdown('**Upload Assay List file (assay-primer mapping)**')
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
            queue_upload(uploaded_assaylists, 'assay_primer_map', caller_id=caller_id)
    
    parse.process_upload_queue(exp)
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
    return success
    

def upload_extra_consumables(key):
    """
    Uploads for extra files such as custom references, assay lists, and taq/water plates.
    Copy the uploaded files into the run folder and deal with subsequent transactions.
    """
    st.markdown('**Upload Custom Reference / Assay Lists**')
    exp = st.session_state['experiment']
    caller_id = 'upload_extra_consumables'
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
        if uploaded_references:
            queue_upload(uploaded_references, 'reference_sequences', caller_id=caller_id)
            
        if uploaded_assaylists:
            queue_upload(uploaded_assaylists, 'assay_primer_map', caller_id=caller_id)
            
        # TODO: add the ability to add a custom taq+water plate, perhaps?
        #if uploaded_taqwater_plates:
        #    success = exp.add_taqwater_layout_volume

    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []


def check_assay_file():
    exp = st.session_state['experiment']
    return any(file.get('purpose', '') == 'assay_primer_map' for file in exp.uploaded_files.values())


def custom_volumes(key):
    exp = st.session_state['experiment']
    volumes_dict = exp.transfer_volumes.copy()
    caller_id = 'custom_volumes'
    st.markdown('**Volumes for DNA and Primers**')
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
        st.info('To change the volume, click or double click on the cell, type in the new value and hit Enter')
    
    if reset_vols:
        success = exp.add_custom_volumes({'DNA_VOL':util.DNA_VOL, 
                'PRIMER_VOL':util.PRIMER_VOL, 
                'PRIMER_TAQ_VOL':util.PRIMER_TAQ_VOL,
                'PRIMER_WATER_VOL':util.PRIMER_WATER_VOL, 
                'INDEX_VOL':util.INDEX_VOL, 
                'INDEX_TAQ_VOL':util.INDEX_TAQ_VOL,
                'INDEX_WATER_VOL':util.INDEX_WATER_VOL}, 
                caller_id=caller_id)
        if success:
            st.info('Transfer volumes reset')
            sleep(1.5)
            st.rerun()
        else:
            st.error('Transfer volumes were not reset. Check the log for more information')

    #if edit is made to dataframe
    if st.session_state['volume_df']['edited_rows']:
        for vol in list(volumes_dict.keys()):
            volumes_dict[vol] = float(custom_vol_editor[vol]['Volume'])
        success = exp.add_custom_volumes(volumes_dict, caller_id)
        if success:
            m(f'Modified transfer volumes {volumes_dict}', level='success')
            sleep(1.5)
            st.rerun()
        else:
            m("Volume could not be modified. Please use integers only", level='error')

    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
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
        queue_upload(uploaded_references, 'reference_sequences', caller_id=caller_id)
        
    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
    

def load_rodentity_data(key):
    """
    Select up to four JSON Rodentity plate definition files for upload
    """
    exp = st.session_state['experiment']
    caller_id = 'load_rodentity_data'
    
    st.info('Add up to four Rodentity JSON ear-punch plate (96-well) files at a time and click on **Submit**. '+\
            'This will move the plate identifiers into the four slots below (P1-4). ')
    free_slots = len([slot for slot in exp.unassigned_plates.keys() if exp.unassigned_plates[slot] == ''])
    with st.form('rodentity_upload_form', clear_on_submit=True):
        epps_col1, _ = st.columns([1,1])
        rodentity_epps = epps_col1.file_uploader(
                f'Choose up to {free_slots} Rodentity JSON files to fill available DNA plate positions',
                type='json',
                accept_multiple_files=True)
        rod_upload_submit_button = st.form_submit_button('Submit')

    if 'custom' not in exp.unassigned_plates:
        exp.unassigned_plates['custom'] = {}
    if rod_upload_submit_button and rodentity_epps:
        if len(rodentity_epps) > free_slots:
            m(f'Please select {free_slots} or less Rodentity JSON files at a time',
                    level='warning', caller_id=caller_id)
            rodentity_epps = []
        else:
            rod_pids = [util.guard_pbc(epp.name.rstrip('.json'), silent=True) for epp in rodentity_epps]
            allowed_pids = []
            for rp in rod_pids:
                if rp in exp.unassigned_plates['custom'] or \
                        (rp in exp.plate_location_sample and   
                        exp.plate_location_sample[rp]['purpose'] != 'sample'):
                    m(f'Rodentity plate barcode {util.unguard_pbc(rp)} already in use for a different purpose, skipping', level='warning')
                else:
                    allowed_pids.append(rp)
            # put Rodentity plates in DNA plate slots
            for rp in allowed_pids:
                for slot in exp.unassigned_plates:
                    if exp.unassigned_plates[slot] == '':
                        exp.unassigned_plates[slot] = rp
                        break
                            
            epps = [ep for ep in rodentity_epps if util.guard_pbc(ep.name.rstrip('.json'), silent=True) in allowed_pids]
            queue_upload(epps, 'rodentity_sample', caller_id=caller_id)
    parse.process_upload_queue(exp)
    #manage transactions:
    if trans.is_pending(exp):
        pending_file_widget(key+'1a', caller_id)
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []


def assign_rodentity_dna_plate(key):
    """
    Manage combining up to four rodentity plates (from JSON) with one destination DNA PID 
    """
    exp = st.session_state['experiment']
    caller_id = "assign_rodentity_dna_plate"
    st.info('Then assign a DNA plate (384-well) barcode identifier for plate combining with the '+\
            'Hamilton Nimbus, and click on **Accept**')
    # this is so ugly. Never use this form style again!
    # declare the checkbox space but update at end of function to make sure it's the latest view
    plates_to_clear = [False, False, False, False]
    rod_col1, _, rod_col2, _ = st.columns([3,1,2,2])
    slotbox_container = rod_col1.container()
         
    with rod_col2.form('rod_destination_form', clear_on_submit=True):
        rod_dp = st.text_input('Destination plate ID (barcode)', 
                        max_chars=30, 
                        key='rod_dp_key')
        if rod_dp:
            rod_dp = util.guard_pbc(rod_dp, silent=True)
        accept_rod_dest_button = st.form_submit_button('Accept')

        if accept_rod_dest_button and rod_dp:
            if rod_dp in exp.plate_location_sample or \
                    rod_dp in exp.unassigned_plates[1] or rod_dp in exp.unassigned_plates[2] or \
                    rod_dp in exp.unassigned_plates[3] or rod_dp in exp.unassigned_plates[4] or \
                    rod_dp in exp.unassigned_plates['custom']:
                m(f'Destination plate barcode already in use: {util.unguard_pbc(rod_dp, silent=True)}',
                        level='warning')
                sleep(1.5)
            else:
                sample_plate_ids = [exp.unassigned_plates[k] for k in [1,2,3,4] if exp.unassigned_plates[k]]
                success = exp.build_dna_plate_entry(sample_plate_ids, rod_dp, source='rodentity')
                if not success:
                    m('Failed to incorporate plate set. Please read the log', 
                            level='error')
                    sleep(1.5)
                else:
                    m(f'Added plate set for {util.unguard_pbc(rod_dp, silent=True)}',
                            level='success')
                    sleep(1)
                    st.rerun()
                        
        
    with slotbox_container.form('set_rod_plates_form', clear_on_submit=True):
        for i in range(4):
            plates_to_clear[i] = st.checkbox(
                    f"P{str(i+1)}: {util.unguard_pbc(exp.unassigned_plates[i+1], silent=True)}", 
                    help='Click the checkbox to allow a set plate ID to be cleared', 
                    key='check'+str(i+1))
        clear_plates_button = st.form_submit_button('Clear IDs', help='Clear selected Rodentity plate IDs')
        if clear_plates_button:
            for i, plate in enumerate(plates_to_clear):
                if plate and exp.unassigned_plates[i+1]:
                    exp.unassigned_plates[i+1] = ''
                    plates_to_clear[i] = False
            st.rerun()

    #manage transactions:
    if trans.is_pending(exp):
        pending_file_widget(key+'1b', caller_id)

    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []
                        
                   
def load_custom_manifests(key):
    """ Demonstration upload panel for custom manifests """
    exp = st.session_state['experiment']
    caller_id = 'load_custom_manifest'
    if 'custom' not in exp.unassigned_plates or not exp.unassigned_plates['custom']:
        exp.unassigned_plates['custom'] = {}
    with st.expander('Upload custom manifests', expanded=True):
        st.info('For non-Rodentity samples in 96-well plates, you can add these here. '+\
                'Upload them first and then choose the layout you want.')
        with st.form(key='manifest_upload_form', clear_on_submit=True):
            col1, _ = st.columns([1,1])
            with col1:
                manifest_uploads = st.file_uploader(label="Custom manifest upload", accept_multiple_files=True)
                #st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                submit_manifest = st.form_submit_button()
                if submit_manifest and manifest_uploads:
                    queue_upload(manifest_uploads, 'custom_sample', caller_id=caller_id)

        parse.process_upload_queue(exp)

        if trans.is_pending(exp):
            pending_file_widget(key+'1a', caller_id)

        if caller_id in mq:
            for msg, lvl in mq[caller_id]:
                m(msg, level=lvl, no_log=True)
            sleep(0.3)
            mq[caller_id] = []

        st.info('Choose the desired layout for the 384-well DNA plate (made with the Hamilton Nimbus). '+\
                'Then provide the DNA plate name and press **Accept**')
        caller_id = 'assign_custom_dna'
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
                        st.warning(f'Destination plate barcode already in use: {util.unguard_pbc(dest_pid, silent=True)}')
                    else:
                        #st.markdown('<p style="color:#FEFEFE">.</p>', unsafe_allow_html=True)
                        sample_pids = [pid for pid in [p1,p2,p3,p4] 
                                if pid != util.guard_pbc('None', silent=True) and pid is not None]
                        success = exp.build_dna_plate_entry(sample_pids, dest_pid, source='custom')
                        if success:
                            m(f'Assigned custom plates to {util.unguard_pbc(dest_pid)}', level='success')
                        else:
                            m('Failed to assign custom plates', level='error')
    
        parse.process_upload_queue(exp)
    
        if trans.is_pending(exp):
            pending_file_widget(key, caller_id)
    
        if caller_id in mq:
            for msg, lvl in mq[caller_id]:
                m(msg, level=lvl, no_log=True)
            sleep(0.3)
            mq[caller_id] = []
        

def add_pcr_barcodes(key):
    caller_id = 'add_pcr_barcodes'
    exp = st.session_state['experiment']

    with st.form('Add PCR plate barcode', clear_on_submit=True):
        pcr_plate_barcode = st.text_input('PCR Plate Barcode', \
                    placeholder='Type in the barcode for PCR plate', key='pcr_plate_barcode' + key)
        
        pcr_plate_entered = st.form_submit_button('Add')

        if pcr_plate_entered and pcr_plate_barcode: 
            guarded_pcr_plate_barcode = util.guard_pbc(pcr_plate_barcode, silent=True)
            if guarded_pcr_plate_barcode not in exp.plate_location_sample:
                success = exp.add_pcr_plates([guarded_pcr_plate_barcode], caller_id=caller_id)
                if success:
                    m(f'Added PCR plate barcode {pcr_plate_barcode}', level='success')
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add PCR plate barcode {pcr_plate_barcode}, please see the log',
                            level='error')
            else:
                m(f'This plate barcode {pcr_plate_barcode} appears to already be in use', 
                        level='error')
    
    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
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
                    m(f'Added taq+water plate {taqwater_plate_barcode}', level='success')
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add taq+water plate barcode {taqwater_plate_barcode}, please see the log',
                            level='error')
            else:
                m(f'This plate barcode {taqwater_plate_barcode} appears to already be in use',
                        level='error')
     
    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []                


def provide_barcodes(key, pcr_stage, dna_pids):
    """
    :param key: unique key for the input widgets
    """
    caller_id = 'provide_barcodes'
    #add_vertical_space(2)
    st.markdown('**Add plate barcodes**')
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
                            caller_id=caller_id)
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add PCR plate barcode {pcr_plate_barcode}, please see the log',
                            level='error', caller_id=caller_id)
            else:
                m(f'This plate barcode {pcr_plate_barcode} appears to already be in use', 
                        level='error', caller_id=caller_id)

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
                            caller_id=caller_id)
                    sleep(1.5)
                    st.rerun()
                else:
                    m(f'Could not add taq+water plate barcode {taqwater_plate_barcode}, please see the log',
                            level='error', caller_id=caller_id)
            else:
                m(f'This plate barcode {taqwater_plate_barcode} appears to already be in use',
                        level='error', caller_id=caller_id)
    
    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []


def upload_miseq_fastqs(key):
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

    parse.process_upload_queue(exp)
    
    if trans.is_pending(exp):
        pending_file_widget(key, caller_id)
    
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = []

