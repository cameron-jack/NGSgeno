#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Core interface script. Run with: streamlit --run ngsgeno.py

Needs load_data.py for GUI functions that are responsible for incorporating data into an
experiment, and display_components.py for functions dedicated to the presentation of GUI
elements
"""
import select
from telnetlib import theNULL
import jsonpickle
import os
#from ssl import SSLSession  # We may want this for secure logins in future
import sys
from pathlib import Path  
from math import fabs, floor, ceil  # leave these incase they're needed later
from subprocess import check_output, CalledProcessError, STDOUT
import subprocess
import queue
import threading
import atexit

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components

from stutil import add_vertical_space, custom_text, hline, init_state, \
        upper_info, upper_height, lower_info, lower_height, m, mq, set_state

from bin.experiment import Experiment, EXP_FN, load_experiment
try:
    import bin.util as util
except ModuleNotFoundError:
    import util
try:
    import bin.transaction as trans
except ModuleNotFoundError:
    import transaction as trans

try:
    import bin.generate as generate
except ModuleNotFoundError:
    import generate

try:
    import bin.parse as parse
except ModuleNotFoundError:
    import parse

#import bin.file_io as file_io
#import bin.db_io as db_io

import extra_streamlit_components as stx
import display_components as dc
import load_data as ld
import asyncio
from time import sleep

def get_status_output(cmd):
    """ calls a process and returns the run status and any output """
    try:
        data = check_output(cmd, shell=True, universal_newlines=True, stderr=STDOUT)
        status = 0
    except CalledProcessError as ex:
        data = ex.output
        status = ex.returncode
    if data[-1:] == '\n':
        data = data[:-1]
    return status, data

def get_run_folders():
    """ return an alphabetically sorted list of run folders, without their run_ prefix """
        
    # get run folder names without 'run_'
    run_folders = [''] + [d[4:] for d in os.listdir('.') if d.startswith('run_') and os.path.isdir(d)]
    return sorted(run_folders)

def create_run_folder(newpath):
    """ returns Experiment or None, and a message string. Requires a full path and generates a new experiment file """
    m(f'Attempting to create: {newpath}')
    if not os.path.exists(newpath):
        try:
            os.mkdir(newpath)
        except Exception as exc:
            m(f'Could not create new folder: {newpath} {exc}', dest=('debug'))
            return None, 'Failed to create new folder: ' + newpath
    else:
        if os.path.exists(os.path.join(newpath, EXP_FN)):
            return None, f'Experiment already exists with this name: {newpath}'
    print('Generating experiment: ', newpath)
    if 'experiment' in st.session_state and st.session_state['experiment']:
        st.session_state['experiment'].save()
    exp = Experiment(name=newpath.lstrip('run_'))
    return exp, ''


async def report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog):
    """
    Allows the interface to keep running while a background process (ngsmatch.py) progress is tracked
    """
    while True:
        progress_files = list(Path(rundir).glob('match_progress_*'))
        if len(progress_files) == 0:
            launch_progress = 0
            match_progress = 0
        elif len(progress_files) == 1:
            launch_progress = int(str(progress_files[0]).split('_')[-2])
            match_progress = int(str(progress_files[0]).split('_')[-1])
        
        if launch_progress > 100:
            launch_progress = 100
        if match_progress > 100:
            match_progress = 100
        launch_msg.write('Allele calling task launch progress: '+str(launch_progress)+'%')
        launch_prog.progress(launch_progress)
        completion_msg.write('Allele calling task completion progress: '+str(match_progress)+'%')
        match_prog.progress(match_progress)

        if launch_progress == 100 and match_progress == 100:
            st.write('Analysis completed')
            st.session_state['match_running'] = False
            return
        await asyncio.sleep(1)


def run_generate(exp, target_func, *args, **kwargs):
    """
    Run target_func() and ask user to respond to any file/pipeline clashes in the given context, if provided
    """
    trans.clear_pending_transactions(exp)
    trans.enforce_file_consistency(exp)
    success = target_func(*args, **kwargs)
    if not success:
        trans.clear_pending_transactions(exp)
    else:
        if exp.pending_steps is not None and len(exp.pending_steps) > 0:
            clashes = trans.clashing_pending_transactions(exp)
            if len(clashes) > 0:
                st.warning(f'The following output files already exist {clashes}')
                st.warning(f'Click "Accept" to replace older files and clear all output files from subsequent pipeline stages')
                st.button('Accept', on_click=trans.accept_pending_transactions, args=[exp], key='accept_overwrite_button')
                st.button('Cancel', on_click=trans.clear_pending_transactions, args=[exp], key='cancel_overwrite_button')
            else:
                trans.accept_pending_transactions(exp)
    return success

def load_experiment_screen():
    """
    Landing screen
    """
    experiment_title = 'Current Experiment: None'
    logo_col, ver_col,_, new_folder_col, create_button_col, ex_folder_col, _ = st.columns([2,2,2,2,1,2,1])
    current_status, current_ver = get_status_output("git describe")
    ver_col.markdown(f'<p style="color:#83b3c9; font-size: 90%"> {current_ver}</p>', 
            unsafe_allow_html=True)

    logo_col.image('ngsg_explorer.png', caption=f'{experiment_title}')
        
    add_run_folder = new_folder_col.text_input('Create new run folder')
    with create_button_col:
        add_vertical_space(2)

    add_run_folder_button = create_button_col.button('Create')
    #create_run_folder_button = ftab1.button(label='Create', key='create_run_folder_button')

    try:
        existing_run_folders = get_run_folders()
    except Exception as exc:
        print(f'Cannot locate NGSgeno folder {exc}', file=sys.stderr)
        return

    run_folder = ex_folder_col.selectbox("Select a run folder to open", existing_run_folders)
    error_msg=''

    #error message comes up now after entering a new folder name - need to fix
    if add_run_folder and add_run_folder_button:
        add_run_folder_str = 'run_' + add_run_folder
        exp, msg = create_run_folder(add_run_folder_str)
        if exp:
            print(f'Saving experiment {exp.name}', file=sys.stderr, flush=True)
            exp.save()                                                                      
            st.session_state['experiment'] = exp
            st.experimental_rerun()
        else:
            if 'already exists' in msg:
                error_msg = "Folder name already exists"
            else:
                error_msg = "Fatal path error: " + msg
                    
    if run_folder:
        if st.session_state['experiment'] == None or st.session_state['experiment'].name != run_folder:
            
            ch_run_path = 'run_' + run_folder
            if os.path.exists(ch_run_path):
                exp = load_experiment(ch_run_path)
                #st.session_state['folder'] = 'existing'
                if not exp:
                    error_msg = "Could not load experiment from: "+ ch_run_path
                elif ch_run_path.endswith(exp.name):
                    # success!
                    st.session_state['experiment'] = exp
                    st.session_state['pipeline_stage'] = 0
                    st.session_state['load_tab'] = 1
                    st.session_state['nimbus_tab'] = 1
                    st.session_state['primer_tab'] = 1
                    st.session_state['index_tab'] = 1
                    st.session_state['miseq_tab'] = 1
                    st.session_state['allele_tab'] = 1
                    st.experimental_rerun()
                else:
                    error_msg = "Invalid experiment file in: " + ch_run_path
    new_folder_col.markdown(f'<p style="color:#FF0000; text-align:center">{error_msg}</p>',\
            unsafe_allow_html=True)


def display_pipeline_header(exp):
    """
    main pipeline header section (constant across pipeline stages)
    exp - Experiment (st.session_state['experiment'])
    """
    experiment_title = 'Current Experiment: ' + exp.name
    logo_col, info_col, pipe_col = st.columns([1,2,9])
    logo_col.image('ngsg_icon.png')
    current_status, current_ver = get_status_output("git describe")
    info_col.markdown(f'<p style="color:#83b3c9; font-size: 90%"> {current_ver}</p>', 
            unsafe_allow_html=True)
    info_col.markdown(f'<p style="color:#83b3c9; font-size: 90%"> {experiment_title}</p>', 
            unsafe_allow_html=True)
    unload_button = info_col.button('🏠', type='primary', help='Go back and change experiment')
    if unload_button:
        if 'experiment' in st.session_state and st.session_state['experiment']:
            st.session_state['experiment'].save()
        st.session_state['experiment'] = None
        st.experimental_rerun()

    if 'stage' not in st.session_state:
        st.session_state['stage'] = None
            
    exp = st.session_state['experiment']
    pipeline_stages=["Load", "Nimbus", "Primers", "Index", "Miseq", "Alleles", "Reports"]
    pipe_stage = None
    with pipe_col:
        pipe_stage = stx.stepper_bar(steps=pipeline_stages, lock_sequence=False)
    return pipe_stage


def unlocked(exp):
    """
    Checks whether the experiment is locked. An experiment will be locked if it the user has uploaded sequence files
    from Miseq. 
    Args:
        exp (st.session_state['experiment']):
    Returns:
        (boolean) True if the experiment is not locked, False if it is
    """
    if exp.locked:
        st.warning(f'Experiment {exp.name} locked from further modification')
        return False
    return True


def last_rites(exp_list):
    """
    Save experiments before the application closes
    """
    for exp in exp_list:
        try:
            print(f'Saving experiment {exp.name}', file=sys.stderr, flush=True) 
            exp.save()
        except Exception as exc:
            print(f'Saving experiment failed! {exc}', file=sys.stderr, flush=True)
    print('Closing application, goodbye.', file=sys.stderr, flush=True)


def main():
    """
    The NGSgeno "Xplorer" application. Allows full control of all sections of the pipeline,
    and displays all aspects of the experiment state at any time.
    """    
    st.set_page_config(
        page_title="NGS Genotyping",
        page_icon="ngsg_icon.png",
        layout="wide"
    )

    dc.add_css()

    init_state('experiment', None)
    init_state('last_rites', False)  # lock to ensure we only register this function once
    init_state('save_list', [])  # a list of experiments that need to be saved at closing time
    if not st.session_state['last_rites']:
        atexit.register(last_rites, st.session_state['save_list'])
        st.session_state['last_rites'] = True

    if 'experiment' not in st.session_state or st.session_state['experiment'] is None:
        load_experiment_screen()
        
    #================================================ START EXPERIMENT =================================================
    else:  # main program
        exp = st.session_state['experiment']
        pipeline_stage = display_pipeline_header(exp)
        
        if pipeline_stage:
            st.session_state['info_expand'] = False

        if not pipeline_stage and pipeline_stage != 0: # not pipeline_stage evaluates to 0!
            if 'pipeline_stage' not in st.session_state or st.session_state['pipeline_stage'] is None:
                if exp.locked:
                    st.session_state['pipeline_stage'] = 5
                else:
                    st.session_state['pipeline_stage'] = 0
            pipeline_stage = st.session_state['pipeline_stage']
        
        #info_bar = dc.info_bar('central')
        init_state('info_expand', False)
        upper_container = st.container()
        message_container = st.container()
        main_body_container = st.container()
        add_vertical_space(2)
        hline()
        lower_container = st.container()
        
        # required for interactive content
        st.session_state['message_container'] = message_container
        
        # callbacks can't write directly as the callbacks go out of scope
        init_state('messages_temp', [])  # messages are tuples of (message, level:info/warning/error/None)
        init_state('messages_persist', [])
        
        # define four info panels, two upper, two lower
        init_state('info_panel1', 'None')
        init_state('info_panel2', 'None')
        init_state('info_panel3', 'Files')
        init_state('info_panel4', 'Log')
        
        # define default heights for these panels
        init_state('upper_panel_height', 250)
        init_state('lower_panel_height', 350)
        
        # standard upper info viewer code for each tab
        def upper_info_viewer_code(tab_col3, tab_col2, widget_key, default_view1='None', default_view2='None', checked=False):
            with tab_col3:
                add_vertical_space(2)
                dc.show_upper_info_viewer_checkbox(widget_key, value=checked)
            if st.session_state['show_upper_info_viewer']:
                with tab_col2:
                    ignore = dc.info_selection(widget_key+"top_viewer", 'info_panel1', 'info_panel2', 
                            'upper_panel_height', default_view1=default_view1, default_view2=default_view2, 
                            default_height=st.session_state.get('upper_panel_height',250))


        # info panel displays are updated at the bottom of the script, so that they reflect any changes
        with lower_container:
            success = dc.info_selection("bottom_viewer", 'info_panel3', 'info_panel4', 
                    'lower_panel_height', default_view1='Files', default_view2='Log', 
                    default_height=st.session_state.get('lower_panel_height',350))
            
        #============================================== STAGE 1: Load data =============================================
        if pipeline_stage == 0:          
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
            with tab_col1:
                load_data_tab = dc.create_tabs([("Load Samples", ""),("Load Consumables", "")])   
            if not load_data_tab:
                init_state('load_tab', 1)
                load_data_tab = st.session_state['load_tab']
            #------------------------------------ Load ~ TAB 1: Load sample data  --------------------------------------
            if load_data_tab == 1:
                st.session_state['load_tab'] = 1
                if unlocked(exp):
                    init_state('run queue', [])
                    with main_body_container:
                        st.subheader('Upload Sample Files')
                        m('**Upload Rodentity plate files, custom manifests, or amplicon plate definition files '+\
                                'then move to the next tab to upload pipeline consumables and other important files**',
                                dest=('mkdn',))
                        ld.load_rodentity_data('rodentity_load1')
                        hline()
                        ld.load_custom_manifests('custom_load1')
                        hline()
                        ld.load_amplicons('amp_load1')

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_load', default_view1='Samples', 
                        default_view2='Files')

            #------------------------------------ Load ~ TAB 2: Load consumables ---------------------------------------
            if load_data_tab == 2:
                st.session_state['load_tab'] = 2

                if unlocked(exp):
                    init_state('upload stage', None)
                    with main_body_container:
                        st.subheader('Upload Consumables')
                        ld.upload_extra_consumables('consumables_load2')
                        ld.upload_pcr1_files('pcr1_load2')
                        ld.upload_pcr2_files('pcr2_load2')
                        add_vertical_space(1)
                        st.subheader('Custom Volumes')
                        ld.custom_volumes(exp)
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_consumables', default_view1='Consumables', 
                        default_view2='Samples')
        
        #=============================================== STAGE 2: Nimbus ===============================================
        if pipeline_stage == 1:
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])

            with tab_col1:
                nimbus_tab = dc.create_tabs([("Download", "Nimbus input files"),("Upload", "Echo input files")])
            if not nimbus_tab:
                init_state('nimbus_tab', 1)
                nimbus_tab = st.session_state['nimbus_tab']

            nfs, efs, xbcs = exp.get_nimbus_filepaths()

            #------------------------------------ Nimbus ~ TAB 1: Download Nimbus --------------------------------------
            if nimbus_tab == 1:
                st.session_state['nimbus_tab'] = 1
                if unlocked(exp):
                    with main_body_container:
                        st.subheader('Generate Echo Files')
                        dc.set_nimbus_title(exp, nfs, efs)
                        add_vertical_space(2)
                        if not ld.check_assay_file(exp):
                            m('Assay list file (assay to primer mapping) is required to proceed. '+\
                                    'Please upload this to continue', level='Error')
                            success = ld.upload_assaylist('nimbus1_assaylist')
                            if success:
                                st.experimental_rerun()
                        else:
                            #Generate files button
                            _, btn_col,_ = st.columns([2,2,1])
                            if st.session_state['experiment'].dest_sample_plates:
                                with btn_col:
                                    run_gen_nimbus = st.button('Generate Nimbus input files', type="primary")

                                #Generate files
                                if run_gen_nimbus:
                                    success = run_generate(exp, exp.generate_nimbus_inputs)    
                                    if not success:
                                        m('Failed to generate the Nimbus files. See the log for details', 
                                                level='error', dest=('debug',))
                                    else:
                                        add_vertical_space(2)
                                        nfs, efs, xbcs = exp.get_nimbus_filepaths()
                    
                        add_vertical_space(3)
                        dc.get_echo_download_buttons(nfs)
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_nimbus1', default_view1='Status', 
                        default_view2='Files')

            #---------------------------------- Nimbus ~ TAB 2: Upload echo input files --------------------------------
            if nimbus_tab == 2:
                st.session_state['nimbus_tab'] = 2
                if unlocked(exp):
                    with main_body_container:
                        _, header_col, _ = st.columns([2,2,1])
                    
                        header_col.subheader('Upload Echo Input Files')
                        ld.upload_echo_inputs('1')

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_nimbus2', default_view1='Status', 
                        default_view2='Files')

        #=========================================== STAGE 3: PCR 1 Primers ============================================
        if pipeline_stage == 2:
            st.session_state['assay_filter'] = True
            pcr_stage = 1
            
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    
            #Tabs
            with tab_col1:
                primer_tab = dc.create_tabs([("PCR 1", "Components"), ("Generate", "Picklists")])
            if not primer_tab:
                init_state("primer_tab", 1)
                primer_tab = st.session_state['primer_tab']
            
            #nimbus fp, echo fp, barcodesnot in echo
            nfs, efs, xbcs = exp.get_nimbus_filepaths()
            missing_nims = ['Echo_384_COC_0001_'+util.unguard(xbc, silent=True)+'_0.csv' for xbc in xbcs]

            #------------------------------------ Primers ~ TAB 1: PCR 1 Components ------------------------------------
            if primer_tab == 1:
                st.session_state['primer_tab'] = 1
                if unlocked:
                    with main_body_container:
                        st.subheader('Primer (PCR 1) Components')
                        if efs:
                            st.info('Provide the barcodes for PCR plates and Taq/water plates and '+\
                                    'upload primer layouts and volumes here, then move to the *Generate Picklists* tab')
                        else:
                            st.error('Load Echo input files (output from Nimbus) to enable PCR1')
                        checkbox_keys = dc.display_plate_checklist('pmr1_checklist', 
                                ['dna','pcr','taqwater1','primer'])
                        selected_pids = dc.collect_plate_checklist(checkbox_keys)
                        hline()    
                        dc.display_pcr_common_components(selected_pids)
                        dc.display_pcr1_components(selected_pids)
                        hline()
                        
                        #barcodes for PCR and taq and water, upload files for primer, adjust volumes
                        add_vertical_space(1)
                        st.subheader('Add Barcodes', help='Add barcodes for plates')
                        ld.provide_barcodes('barcodes_tab1', pcr_stage=pcr_stage)
                        add_vertical_space(1)

                        st.subheader('Upload Files')
                        ld.upload_pcr1_files(key='pcr1_primer1')
                        add_vertical_space(1)

                        st.subheader('Custom Volumes')
                        ld.custom_volumes(exp)

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_pcr1', default_view1='Primers', 
                        default_view2='Files')   

            #-------------------------------- Primers ~ TAB 2: Generate PCR 1 picklists --------------------------------
            if primer_tab == 2:
                st.session_state['primer_tab'] = 2
                if unlocked:
                    with main_body_container:
                        caller_id = 'pcr1 generate'
                        st.subheader('Generate Primer (PCR 1) Echo Picklists')
                        st.info('Select the resources you wish to include for primer PCR and then click on the '+\
                                '**Generate Echo Picklists** button below, to create Echo picklist files')
                        checkbox_keys = dc.display_plate_checklist('pmr1_checklist', 
                                ['dna','pcr','taqwater1','primer'])
                        selected_pids = dc.collect_plate_checklist(checkbox_keys)
                        hline()    
                        dc.display_pcr_common_components(selected_pids)
                        dc.display_pcr1_components(selected_pids)
                        hline()
                        
                        if selected_pids['dna']:
                            if exp.check_ready_pcr1(selected_pids, caller_id):
                                _,button_col,_ = st.columns([2, 2, 1])
                                echo_picklist_go = button_col.button('Generate Echo Picklists',
                                            key='echo_pcr1_go_button',
                                            type='primary')
                                if echo_picklist_go:
                                    success = run_generate(exp, exp.generate_echo_PCR1_picklists,
                                            selected_pids)
                                    if not success:
                                        st.error('Picklist generation failed. Please see the log')
                        for msg,lvl in mq[caller_id]:
                            m(msg, lvl)
                        mq[caller_id] = []

                        if generate.pcr1_picklists_exist(exp):
                            dc.get_echo1_download_btns()
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_pcr2', default_view1='Primers', 
                        default_view2='Consumables') 

        #============================================ STAGE 4: PCR 2 Index =============================================
        if pipeline_stage == 3:
            pcr_stage = 2
            
            tab_col1, tab_col2,tab_col3 = upper_container.columns([5,5,1])

            #Tab setup
            with tab_col1:
                index_tab = dc.create_tabs([("PCR 2", "Components"), ("Generate", "Picklists")])
            if not index_tab:
                init_state('index_tab', 1)
                index_tab = st.session_state['index_tab']
            
            #------------------------------------- Index ~ TAB 1: PCR 2 Components -------------------------------------
            if index_tab == 1:
                st.session_state['index_tab'] = 1
                if unlocked(exp):
                    with main_body_container:
                        st.subheader('Indexing (PCR 2) Components')
                        st.info('Provide the resources needed to perform a sufficient number of indexing reactions '+\
                                'for your experiment, then move to the *Generate Picklists* tab')
                        checkbox_keys = dc.display_plate_checklist('idx_checklist1',
                                ['pcr','taqwater2','amplicon','index'])
                        hline()
                        selected_pids = dc.collect_plate_checklist(checkbox_keys)
                        dc.display_pcr_common_components(selected_pids)
                        dc.display_pcr2_components(selected_pids)
                        hline()
                        add_vertical_space(1)
                        st.info('Indexing (PCR 2) requires index index plates, taq/water plates, '+\
                                'and either Echo plates (prepared by the Nimbus) or amplicon plates')

                        #provide barcodes for pcr & taq/water, upload index files, adjust volumes
                        st.subheader('Add Barcodes', help='Add barcodes for plates')
                        ld.provide_barcodes('index_barcodes', 2)
                        add_vertical_space(1)

                        st.subheader('Upload Files')
                        ld.upload_pcr2_files(key='pcr2_index1')
                        add_vertical_space(1)

                        st.subheader('Custom Volumes')
                        ld.custom_volumes(exp)
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_index1', default_view1='Indexes', 
                        default_view2='Consumables')

            #--------------------------------- Index ~ TAB 2: Generate PCR 2 picklists ---------------------------------
            if index_tab == 2:                
                st.session_state['index_tab'] = 2 
                init_state('run_index', False)
                if unlocked(exp):
                    with main_body_container:
                        caller_id = 'pcr2 generate'
                        st.subheader('Generate Index (PCR 2) Echo Picklists')
                        st.info('Select the resources you wish to include for indexing PCR and then click on the '+\
                                '**Generate Echo Picklists** button below, to create Echo picklist files')
                        do_generate = False
                        checkbox_keys = dc.display_plate_checklist('idx_checklist2', 
                                ['pcr','taqwater2','amplicon','index'])
                        selected_pids = dc.collect_plate_checklist(checkbox_keys)
                        hline()
                        dc.display_pcr_common_components(selected_pids)
                        dc.display_pcr2_components(selected_pids)
                        hline()
                        add_vertical_space(1)
                                                
                        if selected_pids['pcr']:
                            success = exp.check_ready_pcr2(selected_pids, caller_id)
                            if success:
                                do_generate = True
                        elif selected_pids['amplicon']:
                            success = exp.check_ready_pcr2(selected_pids, caller_id, amplicon_only=True)
                            if success:
                                do_generate = True
                        else:
                            m('Either PCR plates or Amplicon plates must exist for indexing to begin')
                                    
                        for msg,lvl in mq[caller_id]:
                            m(msg, lvl)
                        mq[caller_id] = []

                        if do_generate:
                            _,picklist_button_col,_ = st.columns([2, 2, 1])
                            echo_picklist_go = picklist_button_col.button('Generate Echo Picklists',\
                                    key='echo_pcr2_go_button')
                            picklist_button_col.write('')
                            if echo_picklist_go:
                                if selected_pids['pcr']:
                                    st.session_state['run_index'] = True
                                elif selected_pids['amplicon'] and not selected_pids['pcr']:
                                    _, amp1, amp2, amp3, _ = st.columns([3, 3, 1, 1, 3])
                                    amp1.warning('Create picklists with only amplicons?')
                                    amp2.button('Yes', on_click=set_state, args=('run_index', True))
                                    amp3.button('No', on_click=set_state, args=('run_index', False))
                        
                        if st.session_state['run_index']:
                            success = run_generate(exp, exp.generate_echo_PCR2_picklists,
                                    selected_pids)
                            set_state('run_index', False)
                            if not success:
                                m('Picklist generation failed. Please see the log')
                                    
                        if generate.pcr2_picklists_exist(exp):
                            dc.show_echo2_outputs() 
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_index2', default_view1='Status', 
                        default_view2='Files')

        #=============================================== STAGE 5: Miseq ================================================
        if pipeline_stage == 4:
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
            info_holder = st.container()

            with tab_col1:
                with tab_col1:
                    miseq_tab = dc.create_tabs([("Download", "Miseq Samplesheet"), ("Upload", "Miseq Sequence Files")])
                if not miseq_tab:
                    init_state('miseq_tab', 1)
                    miseq_tab = st.session_state['miseq_tab']

            add_vertical_space(1)
            
            #-------------------------------- Miseq ~ TAB 1: Download Miseq Samplesheet --------------------------------
            if miseq_tab == 1:
                st.session_state['miseq_tab'] = 1
                with main_body_container:
                    _,header_col,_ = st.columns([2,2,1])

                    if exp.locked:
                        st.warning(f'Experiment {exp.name} locked from further modification')
                
                    with header_col:
                        st.subheader('Download MiSeq File')
                        add_vertical_space(1)
                    #hline()

                    #ld.upload_reference_sequences('reference_miseq1')
                    if exp.get_miseq_samplesheets():
                        dc.get_miseq_download_btn(exp)
                        add_vertical_space(4)
                    
                    else:
                        st.warning(f'No MiSeq Samplesheet available for download')
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_miseq1', default_view1='Status', 
                        default_view2='Plates', checked=False)

            #------------------------ Miseq ~ TAB 2: Upload Reference File and Miseq Sequences -------------------------
            if miseq_tab == 2:
                st.session_state['miseq_tab'] = 2
                with main_body_container:
                    st.subheader('Upload Custom Reference Files')
                    ld.upload_reference_sequences('reference_miseq2')
                    add_vertical_space(2)
                
                    ready_messages = []
                    if not exp.check_sequence_upload_ready(ready_messages):
                        for msg in ready_messages:
                            st.error(msg)
                        st.warning('These resources are required for allele calling and must be present before FASTQs can be uploaded')
                    else:
                        ld.upload_miseq_fastqs()

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_miseq2', default_view1='Status', 
                        default_view2='Plates', checked=False)

        #=========================================== STAGE 6: Allele Calling ===========================================
        if pipeline_stage == 5:
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])

            with tab_col1:
                allele_tab = dc.create_tabs([("Allele Calling", "")])
            if not allele_tab:
                init_state('allele_tab', 1)
                allele_tab = st.session_state['allele_tab']

            # Only offer upload in the Miseq pipeline section
            #-------------------------------- Allele ~ TAB 1: Allele calling --------------------------------
            if allele_tab == 1:
                st.session_state['allele_tab'] = 1
                with main_body_container:
                    rundir = exp.get_exp_dn()
                    caller_id = 'pre-execute-analysis'    
                    success = exp.check_sequence_upload_ready(caller_id)
                    for msg,lvl in mq[caller_id]:
                        m(msg, level=lvl)
                    mq[caller_id] = []
                    if not success:
                        st.warning('Resources are required for allele calling and must be present before FASTQs can be uploaded')
                        st.subheader('Upload reference sequences')
                        ld.upload_reference_sequences('reference_allele1')
                        ld.upload_miseq_fastqs()
                    else:
                        ld.upload_miseq_fastqs()
                        # check whether output files already exist and are opened elsewhere
                        target_fn = exp.get_exp_fn('target.fa')
                        results_fn = exp.get_exp_fn('results.csv')
                        matchlog_fn = exp.get_exp_fn('match.log')
                        for fn in [target_fn, results_fn, matchlog_fn]:
                            if Path(fn).exists():
                                try:
                                    os.rename(fn, fn.replace('.','_tmp_swap.'))
                                    os.rename(fn.replace('.','_tmp_swap.'), fn)
                                except PermissionError:
                                    st.error(f'{fn} appears to be in use. Please close this file before continuing')
                            
                        with st.form('allele_calling_form', clear_on_submit=True):
                            #num_unique_seq = st.number_input("Number of unique sequences per work unit", value=1)
                            cpus_avail = os.cpu_count() -1
                            num_cpus = st.number_input(\
                                    label=f"Number of processes to run simultaneously, default: {cpus_avail}",\
                                            value=cpus_avail)
                            mincov = st.number_input(label="Do not match unique sequences with less than this "+\
                                    "many reads coverage, default 5", format='%i',min_value=0, step=1,value=5)
                            minprop = st.number_input(label="Do not match unique sequences with less than this "+\
                                    "proportion of the reads seen for the most observed (expected) allele, default 0.1. Must be between 0.0 and 1.0",
                                    format='%f',min_value=0.0, max_value=1.0, value=0.1, step=0.05)
                            exact_only = st.checkbox("Exact only: disable inexact matching")
                            exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                            debug_mode = st.checkbox('Turn on debugging for allele calling')
                            do_matching = st.form_submit_button("Run allele calling")

                        if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                            st.info('Analysis in progress')
                        if do_matching:
                            success = generate.generate_targets(exp)
                            if not success:
                                m('Critical: failed to save reference sequences to target file', dest=('log',))
                                sleep(0.5)
                            success = generate.generate_primer_assayfams(exp)
                            if not success:
                                m('Critical: failed to save primers and assay families to file', dest=('log'))
                                sleep(0.5)
                            else:
                                matching_prog = os.path.join('bin','ngsmatch.py')
                                cmd_str = f'python {matching_prog} --ncpus {num_cpus} --rundir {rundir} --mincov {mincov} --minprop {minprop}'
                                if exact_only:
                                    cmd_str += ' --exact'
                                if exhaustive_mode:
                                    cmd_str += ' --exhaustive'
                                if debug_mode:
                                    cmd_str += ' --debug'
                                exp.log(f'Info: {cmd_str}')
                                #print(f"{cmd_str=}", file=sys.stderr)
                                subprocess.Popen(cmd_str.split(' '))
                                st.write(f'Calling {cmd_str}')
                                sleep(1.0)

                    if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                        launch_msg = st.empty()
                        launch_prog = st.progress(0)
                        completion_msg = st.empty()
                        match_prog = st.progress(0)
                        asyncio.run(report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog))

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_allele1', default_view1='Status', 
                        default_view2='Plates', checked=False)
        
        #=============================================== STAGE 7: Reports ==============================================
        if pipeline_stage == 6:
            results_fp = exp.get_exp_fn('results.csv')
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
            with tab_col1:
                allele_tab = dc.create_tabs([("Sequencing results", "")])
            if not allele_tab:
                init_state('results_tab', 1)
                allele_tab = st.session_state['results_tab']
                
            with main_body_container:
                if not Path(results_fp).exists():
                    m('**No allele calling results available**', dest=('mkdn'))
                else:
                    rodentity_results = []
                    custom_results = []
                    other_results = []
                    with open(results_fp, 'rt') as rfn:
                        for i, line in enumerate(rfn):
                            l = line.replace('"','')
                            cols = [c.strip() for c in l.split(',')]
                            #if len(cols) != 30:
                            #    print(cols)
                            if i == 0:
                                hdr = cols
                            else:
                                sample = cols[3]
                            
                                if util.is_guarded_cbc(sample):
                                    custom_results.append(cols)
                                elif util.is_guarded_rbc(sample):
                                    rodentity_results.append(cols)
                                else:
                                    other_results.append(cols)

                    #print(hdr)
                    #print(custom_results[0:3])

                    m(f'{len(custom_results)=} {len(rodentity_results)=} {len(other_results)=}', dest=('css',))
                    rodentity_view = st.expander('Rodentity results: '+str(len(rodentity_results)))
                    with rodentity_view:
                        dfr = pd.DataFrame(rodentity_results, columns=hdr)
                        dc.aggrid_interactive_table(dfr, key='rodentity_view_key')

                    custom_view = st.expander('Custom results: '+str(len(custom_results)))
                    with custom_view:
                        dfc = pd.DataFrame(custom_results, columns=hdr)
                        dc.aggrid_interactive_table(dfc, key='custom_view_key')
                    other_view = st.expander('Other results: '+str(len(other_results)))
                    with other_view:
                        dfo = pd.DataFrame(other_results, columns=hdr)
                        dc.aggrid_interactive_table(dfo, key='other_view_key')
            # ** Info viewer **
            upper_info_viewer_code(tab_col3, tab_col2, 'upper_report1', default_view1='Status', 
                    default_view2='Plates', checked=False)

        st.session_state['pipeline_stage'] = pipeline_stage

        #=============================================== UPPER INFO SECTION ============================================
        
        upper_panels = [v for v in (st.session_state.get('info_panel1', 'None'), 
                st.session_state.get('info_panel2', 'None')) if v != 'None']
        if any(upper_panels) and st.session_state.get('show_upper_info_viewer'):
            with upper_container:
                dc.show_info_viewer(upper_panels, st.session_state.get('upper_panel_height',250), 'upper_view_panels')  
        
        #=============================================== LOWER INFO SECTION ============================================
        
        lower_panels = [v for v in (st.session_state.get('info_panel3','None'),
                st.session_state.get('info_panel4','None')) if v != "None"]
        if any(lower_panels):
            with lower_container:
                dc.show_info_viewer(lower_panels, st.session_state.get('lower_panel_height',350), 'lower_view_panels')

        #=============================================== UPDATE MESSAGES ==============================================
        with message_container:
            #dc.display_temporary_messages()
            dc.display_persistent_messages(key='main1')


        ### End of main display ###

if __name__ == '__main__':
    main()


