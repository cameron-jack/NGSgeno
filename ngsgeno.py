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
import jsonpickle
import os
#from ssl import SSLSession  # We may want this for secure logins in future
import sys
from pathlib import Path  
from math import fabs, floor, ceil  # leave these incase they're needed later
from subprocess import check_output, CalledProcessError, STDOUT
import subprocess

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components

from stutil import add_vertical_space, custom_text, hline, init_state, \
        upper_info, upper_height, lower_info, lower_height

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
    print('Attempting to create: ', newpath)
    if not os.path.exists(newpath):
        try:
            os.mkdir(newpath)
        except Exception as exc:
            print(f'Could not create new folder: {newpath} {exc}', file=sys.stderr)
            return None, 'Failed to create new folder: ' + newpath
    else:
        if os.path.exists(os.path.join(newpath, EXP_FN)):
            return None, 'Experiment already exists with this name: ' + newpath
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
    #print(f"launching {target_func=}", file=sys.stderr)
    success = target_func(*args, **kwargs)
    #print(f"completed {target_func=} {success=}", file=sys.stderr)
    if not success:
        trans.clear_pending_transactions(exp)
    else:
        if exp.pending_steps is not None and len(exp.pending_steps) > 0:
            clashes = trans.clashing_pending_transactions(exp)
            if len(clashes) > 0:
                #if context:
                #    context.warning(f'The following output files already exist {clashes}')
                #    context.warning(f'Click "Accept" to replace older files and clear all output files from subsequent pipeline stages')
                #    if context.button('Accept'):
                #        exp.accept_pending_transactions()
                #    if context.button('Cancel'):
                #        exp.clear_pending_transactions()
                #else:
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


def pcr1_picklists_exist(exp):
    """
    *Stage 3: PCR1*
    Check if the files that are generated for PCR 1 exist
    """
    picklist_files = ['Stage2.csv', 
                        'PCR1_dna-picklist_test_picklists.csv', 
                        'PCR1_primer-picklist_test_picklists.csv', 
                        'PCR1_taqwater-picklist_test_picklists.csv']
    
    return all(os.path.exists(exp.get_exp_fn(file)) for file in picklist_files)

def get_echo_picklist_btn_pcr1(exp, DNA_plates, PCR_plates, taqwater_plates):
    """
    *Stage 3: PCR1*
    Display for generate echo pcr 1 picklist
    Args:
        exp (st.session_state['experiment])
        DNA_plates: included DNA plates
        PCR_plates: included PCR plates
        taqwater_plates: included Taq/water plates
    """
    pcr1_messages = []
    if not exp.check_ready_pcr1(DNA_plates,\
                                PCR_plates, \
                                taqwater_plates, \
                                pcr1_messages):
        
        for msg in pcr1_messages:
            st.warning(msg)
    
    else:
        _,button_col,_ = st.columns([2, 2, 1])
        echo_picklist_go = button_col.button('Generate Echo Picklists',
                                              key='echo_pcr1_go_button',
                                              type='primary')
        if echo_picklist_go:
            success = run_generate(exp, 
                                   exp.generate_echo_PCR1_picklists, 
                                   DNA_plates, 
                                   PCR_plates, 
                                   taqwater_plates)
            if not success:
                st.error('Picklist generation failed. Please see the log')
            else:
                st.session_state['pcr1 picklist'] = True



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
        def upper_info_viewer_code(tab_col3, tab_col2, widget_key, default_view1='None', default_view2='None'):
            with tab_col3:
                add_vertical_space(2)
                dc.show_upper_info_viewer_checkbox(widget_key)
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
                    stx.TabBarItemData(id=1, title="Load Samples", description=""),
                init_state('load_tab', 1)
                    #stx.TabBarItemData(id=3, title="View Data", description="")
            #------------------------------------ Load ~ TAB 1: Load sample data  --------------------------------------
            if not load_data_tab:
                st.session_state['load_tab'] = 1
                if unlocked(exp):
                    init_state('run queue', [])
                    with main_body_container:
                        st.subheader('Upload Sample Files')
                        ld.load_rodentity_data('rodentity_load1')
                        ld.load_custom_manifests('custom_load1')
                        ld.load_amplicons('amp_load1')

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_load', default_view1='Samples', 
                        default_view2='Files')

            #------------------------------------ Load ~ TAB 2: Load consumables ---------------------------------------
                    ld.load_amplicons('amp_load1')
                st.session_state['load_tab'] = 2

                if unlocked(exp):
                    init_state('upload stage', None)
                    with main_body_container:
                        st.subheader('Custom Volumes')
                        ld.custom_volumes(exp)
                        add_vertical_space(1)

                        st.subheader('Upload Consumables')
                        ld.upload_extra_consumables('consumables_load2')
                        ld.upload_pcr1_files('pcr1_load2')
                        ld.upload_pcr2_files('pcr2_load2')
                else:
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_consumables', default_view1='Consumables', 
                        default_view2='Samples')
                st.session_state['load_tab'] = 2

            with tab_col2:
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])

            with tab_col1:
                nimbus_tab = dc.create_tabs([("Download", "Nimbus input files"),("Upload", "Echo input files")])
                nimbus_tab = stx.tab_bar(data=[
                init_state('nimbus_tab', 1)
                    stx.TabBarItemData(id=2, title="Upload", description="Echo input files")

            if not ld.check_assay_file(exp):
                with message_container:
                    st.warning("Upload assay list file before generating Echo files")

            info_holder = st.container()
            
            
            if not nimbus_tab:
                st.session_state['nimbus_tab'] = 1
                if unlocked(exp):
                    with main_body_container:
                        _, header_col, _ = st.columns([2,2,1])
                        header_col.subheader('Generate Echo Files')
            nfs, efs, xbcs = exp.get_nimbus_filepaths()
                        #Subtitle
                        nimbus_title = dc.set_nimbus_title(exp, nfs, efs)
                        if nimbus_title:
                            _, title_col,_ = main_body_container.columns([2,2,1])
                            with title_col:
                                custom_text('h5', '#83b3c9', dc.set_nimbus_title(exp, nfs, efs), align="left")
                        add_vertical_space(1)

                        #Generate files button
                        _, btn_col,_ = st.columns([2,2,1])
                        with btn_col:
                            run_gen_nimbus = st.button('Generate Nimbus input files', type="primary")
                        nimbus_title = "Load data inputs to enable Nimbus input file generation."
                        #Generate files
                        if run_gen_nimbus:
                            success = run_generate(exp, exp.generate_nimbus_inputs)    
                            if not success:
                                message_container.error('Failed to generate the Nimbus files. Check the log for information.')
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
                        st.write(f'<p style="color:#FF0000">{nim_tab1_err}</p>', unsafe_allow_html=True)
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_nimbus2', default_view1='Status', 
                        default_view2='Files')

                    if (i+1) % 2 != 0:
                        dl_col1.markdown('<p style="text-align:left;color:#4b778c;padding:5px">'+
                                        unsafe_allow_html=True)

            
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    
            #Tabs
            with tab_col1:
                primer_tab = dc.create_tabs([("PCR 1", "Components"), ("Generate", "Picklists")])
                        dl_col4.download_button("Download ", 
                init_state("primer_tab", 1)
                                                key='nimbus_input'+str(i), 
                                                help=f"Download Nimbus input file {nf}") 
                st.session_state['nimbus_tab'] = 1

            #Upload nimbus
            if nimbus_tab == 2:
                ld.upload_echo_inputs('1')                
                st.session_state['nimbus_tab'] = 2
                st.session_state['primer_tab'] = 1
                if unlocked:
                    with main_body_container:
                        tip_col, _ = st.columns([2,1])
                        tip_col.info('Provide the barcodes for PCR plates and Taq/water plates and '+\
                                'upload primer layouts and volumes to generate the picklists.')
                        primer_checklist = st.container()
                        primer_checklist.subheader('Plate Checklist')
                        hline()

                        if efs:
                            with primer_checklist:
                                included_DNA_pids, included_PCR_pids, included_taqwater_pids =\
                                            dc.plate_checklist_pcr1(exp)
        if pipeline_stage == 2:
                            if included_DNA_pids:
                                st.subheader('PCR 1 Components', help='Required plates and volumes for the PCR reaction')
                                dc.display_pcr_common_components(dna_pids=included_DNA_pids)
                                dc.display_pcr1_components(dna_pids=included_DNA_pids)
                                hline()
                        else:
                            st.error("Load Nimbus output files to enable PCR stages")
            tab_col1, tab_col2 = st.columns([9,1])
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
                    stx.TabBarItemData(id=2, title="Generate", description="Picklists")
                ], return_type=int)
            info_holder = st.container()
            
            if not primer_tab:
                if 'primer_tab' not in st.session_state:
                    st.session_state['primer_tab'] = 1
                st.session_state['primer_tab'] = 2
                if unlocked:
                    with main_body_container:
                        init_state('pcr1 picklist', False)

                        st.subheader('Plate checklist')
                        primer_checklist_exp = st.container()
                        hline()
                        if not efs:
                            st.warning('No DNA plate information available. Have you uploaded Echo input files yet?')
                        else:
                            with primer_checklist_exp:
                                included_DNA_pids, included_PCR_pids,\
                                        included_taqwater_pids = dc.plate_checklist_pcr1(exp)
                                    
                                st.session_state['included_DNA_pids'] = included_DNA_pids
                    st.write('***')
                            if included_DNA_pids:
                                get_echo_picklist_btn_pcr1(exp, included_DNA_pids, included_PCR_pids,\
                                                                    included_taqwater_pids)

                if st.session_state['pcr1 picklist'] or pcr1_picklists_exist(exp):
                    dc.get_echo1_download_btns()
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_pcr2', default_view1='Primers', 
                        default_view2='Consumables') 
                        default_view2='Files')   

            #-------------------------------- Primers ~ TAB 2: Generate PCR 1 picklists --------------------------------
            pcr_stage = 2
            
            tab_col1, tab_col2,tab_col3 = upper_container.columns([5,5,1])
                else:
            #Tab setup
            with tab_col1:
                index_tab = dc.create_tabs([("PCR 2", "Components"), ("Generate", "Picklists")])
                        with primer_checklist_exp:
                init_state('index_tab', 1)
                                    dc.plate_checklist_pcr1(exp)
            
            #------------------------------------- Index ~ TAB 1: PCR 2 Components -------------------------------------
                        if included_DNA_plates:
                st.session_state['index_tab'] = 1
                if unlocked(exp):
                    with main_body_container:
                        index_checklist = st.container()
                        index_checklist.subheader('Plate Checklist')
                        hline()
                    
                        title_holder = st.empty()
                        pcr_comp_holder = st.container()
                                echo_picklist_go = picklist_button_col.button('Generate Echo Picklist',
                        #provide barcodes for pcr & taq/water, upload index files, adjust volumes
                        st.subheader('Add Barcodes', help='Add barcodes for plates')
                        ld.provide_barcodes('index_barcodes', 2)
                        add_vertical_space(1)

                        st.subheader('Upload Files')
                        ld.upload_pcr2_files(key='pcr2_index1')
                        add_vertical_space(1)

                        st.subheader('Custom Volumes')
                        ld.custom_volumes(exp)
                    
                        with index_checklist:
                            checkbox_keys = dc.display_plate_checklist('idx_checklist', inc_pcr=True, 
                                    inc_taqwater=True, inc_amplicon=True, inc_index=True)
                            # included_PCR_pids, included_taqwater_pids, included_index_pids, \
                            #         included_amplicon_pids = dc.plate_checklist_pcr2(exp)
                        
                        with pcr_comp_holder:
                            st.subheader('PCR 2 Components')
                            selected_pids = dc.collect_checklists(checkbox_keys)
                            dc.display_pcr_common_components(pcr_pids=selected_pids['pcr'],
                                    amplicon_pids=selected_pids['amplicon'])
                            dc.display_pcr2_components(pcr_pids=selected_pids['pcr'], 
                                    amplicon_pids=selected_pids['amplicon'], 
                                    taqwater_pids=selected_pids['taqwater'])
                            hline()
                            add_vertical_space(1)
                st.write('')
                        if not selected_pids['pcr'] and not selected_pids['amplicon']:
                            title_holder.error("Load Nimbus output files to enable PCR stages. For amplicon "+\
                                    "only, upload amplicon and index files and provide a taq water plate barcode.")
                
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_index1', default_view1='Indexes', 
                        default_view2='Consumables')
        #============================================ STAGE 4: PCR 2 Index =============================================
            #--------------------------------- Index ~ TAB 2: Generate PCR 2 picklists ---------------------------------
            if index_tab == 2:                
                st.session_state['index_tab'] = 2 
                if unlocked(exp):
                    with main_body_container:
                        index_checklist = st.container()
                        index_checklist.subheader('Plate Checklist')
                        hline()
        
                        with index_checklist:
                            do_generate = False
                            included_PCR_pids, included_taqwater_pids, included_index_pids,\
                                    included_amplicon_pids = dc.plate_checklist_pcr2(exp)
                            if included_PCR_pids:  # standard run
                                if not exp.check_ready_pcr2(included_PCR_pids, included_taqwater_pids, 
                                        included_index_pids, included_amplicon_pids):
                                    st.warning('Cannot generate PCR2 picklists, please see the log for details')
                                else:
                                    do_generate = True
                            else:
                                if included_amplicon_pids:
                                    if not exp.check_ready_pcr2_amplicon_only(included_taqwater_pids,
                                            included_index_pids, included_amplicon_pids):
                                        st.warning('Cannot generate PCR2 picklists, please see the log for details')
                                    else:
                                        do_generate = True
                    title_holder = st.empty()
                        if do_generate:
                            show_generate = False
                            if 'amplicon_only' in st.session_state:
                                show_generate = st.session_state['amplicon_only']
                            if included_PCR_pids:
                                show_generate = True

                            if included_amplicon_pids and not included_PCR_pids:
                                _, amp1, amp2, amp3, _ = st.columns([3, 3, 1, 1, 3])
                                amp1.warning('Is this an amplicon only run?')
                                yes_amplicon = amp2.button('Yes')
                                no_amplicon = amp3.button('No')
                        dc.display_pcr_components(pcr_stage=pcr_stage)
                                if yes_amplicon:
                                    st.session_state['amplicon_only'] = True
                                    show_generate = True
                                elif no_amplicon:
                                    st.session_state['amplicon_only'] = False
                                    show_generate = False
                    
                            if show_generate:
                                _,picklist_button_col,_ = st.columns([2, 2, 1])
            if index_tab == 2:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    st.write('Plate checklist')
                    index_checklist = st.container()
                                if echo_picklist_go:
                                    st.session_state['idx_picklist'] = True
                        
                                    success = run_generate(exp, exp.generate_echo_PCR2_picklists, included_PCR_pids,
                                            included_index_pids, included_taqwater_pids, included_amplicon_pids)
                                    if not success:
                                        st.write('Picklist generation failed. Please see the log')
                        if included_PCR_plates:  # standard run
                        dc.show_echo2_outputs()           
                                for msg in pcr2_messages:
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_index2', default_view1='Status', 
                        default_view2='Files')
                                    for msg in pcr2_messages:
                                        st.warning(msg)
                                else:
            tab_col1, tab_col2, tab_col3 = main_body_container.columns([9,3,1])
                        if do_generate:

            with tab_col1:
                with tab_col1:
                    miseq_tab = dc.create_tabs([("Download", "Miseq Samplesheet"), ("Upload", "Miseq Sequence Files")])
                if not miseq_tab:
                    init_state('miseq_tab', 1)
                    miseq_tab = st.session_state['miseq_tab']
                            if included_amplicon_plates and not included_PCR_plates:
            add_vertical_space(1)
            
            #-------------------------------- Miseq ~ TAB 1: Download Miseq Samplesheet --------------------------------
            if miseq_tab == 1:
                st.session_state['miseq_tab'] = 1
                _,header_col,_ = st.columns([2,2,1])


                                if yes_amplicon:
                                    st.session_state['amplicon_only'] = True
                                    show_generate = True
                                elif no_amplicon:
                                    st.session_state['amplicon_only'] = False
                                    show_generate = False
                        
                            if show_generate:
                                _,picklist_button_col,_ = st.columns([2, 2, 1])

                            echo_picklist_go = picklist_button_col.button('Generate Echo Picklists',\
                                        key='echo_pcr2_go_button')

                            picklist_button_col.write('')
                                if echo_picklist_go:
                                    st.session_state['idx_picklist'] = True
                            
                st.session_state['miseq_tab'] = 2
                st.subheader('Upload Custom Reference Files')
                                    success = run_generate(exp, exp.generate_echo_PCR2_picklists, included_PCR_plates,
                                            included_index_plates, included_taqwater_plates, included_amplicon_plates)
                                    if not success:
                                        st.write('Picklist generation failed. Please see the log')
                        

                dc.show_echo2_outputs()
                st.session_state['index_tab'] = 2            
                
            with tab_col2:
                st.write('')
             # ** Info viewer **
            with tab_col3:
                add_vertical_space(2)
                dc.show_upper_info_viewer_checkbox()
            if st.session_state['show_upper_info_viewer']:
                with tab_col2:
                    selection, height = dc.info_selection(5)
                if selection:
                    with info_holder:
                        dc.show_info_viewer(selection, height)
            tab_col1, tab_col2 = st.columns([9,1])
            with tab_col1:
                miseq_tab = stx.tab_bar(data=[
            tab_col1, tab_col2, tab_col3 = main_body_container.columns([9,3,1])
                ], return_type=int)
            info_holder = st.container()
            
   
            if not miseq_tab:
                init_state('allele_tab', 1)
                    st.session_state['miseq_tab'] = 1
                miseq_tab = st.session_state['miseq_tab']

            st.write('')
            exp = st.session_state['experiment']
            if miseq_tab == 1:
                _, miseq_col1, miseq_col2, _ =  st.columns([2,1,1,2])
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
                st.session_state['miseq_tab'] = 1

            #------------------------ Miseq ~ TAB 2: Upload Reference File and Miseq Sequences -------------------------
            if miseq_tab == 2:
                ld.upload_reference_sequences('reference_miseq2')
                add_vertical_space(2)
                
                ready_messages = []
                if not exp.check_sequence_upload_ready(ready_messages):
                    for msg in ready_messages:
                        st.error(msg)
                    st.warning('These resources are required for allele calling and must be present before FASTQs can be uploaded')
                else:
                    ld.upload_miseq_fastqs()
                st.session_state['miseq_tab'] = 2

            with tab_col2:
                st.write('')
                st.write('')
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #=========================================== STAGE 6: Allele Calling ===========================================
        if pipeline_stage == 5:
            exp = st.session_state['experiment']
            tab_col1, tab_col2 = st.columns([9,1])
            with tab_col1:
                allele_tab = stx.tab_bar(data=[
                    #stx.TabBarItemData(id=1, title="Upload", description="Sequence Files"),
                    stx.TabBarItemData(id=1, title="Allele Calling", description="")
                    #stx.TabBarItemData(id=3, title="View Data", description=""),
                ], key='allele_tab_bar' , return_type=int)
            info_holder = st.container()

            with tab_col1:
                allele_tab = dc.create_tabs([("Allele Calling", "")])
            if not allele_tab:
                if 'allele_tab' not in st.session_state:
                    st.session_state['allele_tab'] = 1
                allele_tab = st.session_state['allele_tab']

            # Only offer upload in the Miseq pipeline section
            #-------------------------------- Allele ~ TAB 1: Allele calling --------------------------------
            if allele_tab == 1:
                st.session_state['allele_tab'] = 1
                rundir = exp.get_exp_dn()
                seq_ready_messages = []
                call_ready_messages = []
                ld.upload_reference_sequences('reference_allele1')
                if not exp.check_sequence_upload_ready(seq_ready_messages):
                    for msg in seq_ready_messages:
                        st.error(msg)
                    st.warning('These resources are required for allele calling and must be present before FASTQs can be uploaded')

                elif not exp.check_allele_calling_ready(call_ready_messages):
                    for msg in call_ready_messages:
                        st.error(msg)
                    st.warning('These resources are required before allele calling can proceed')
                    ld.upload_miseq_fastqs()
                else:
                    ld.upload_miseq_fastqs()
                    # check whether output files already exist and are opened elsewhere

             # ** Info viewer **
            with tab_col3:
                add_vertical_space(2)
                dc.show_upper_info_viewer_checkbox()
            if st.session_state['show_upper_info_viewer']:
                with tab_col2:
                    selection, height = dc.info_selection(6)
                if selection:
                    with info_holder:
                        dc.show_info_viewer(selection, height)
                            
                    with st.form('allele_calling_form', clear_on_submit=True):
                        #num_unique_seq = st.number_input("Number of unique sequences per work unit", value=1)
                        cpus_avail = os.cpu_count() -1
            tab_col1, tab_col2, tab_col3 = main_body_container.columns([9,3,1])
                                label=f"Number of processes to run simultaneously, default: {cpus_avail}",\
                                        value=cpus_avail)
                        mincov = st.number_input(label="Do not match unique sequences with less than this "+\
                                "many reads coverage, default 50", format='%i',min_value=0, step=1,value=50)
                        minprop = st.number_input(label="Do not match unique sequences with less than this "+\
                                "proportion of the reads seen for the most observed (expected) allele, default 0.2. Must be between 0.0 and 1.0",
                                format='%f',min_value=0.0, max_value=1.0, value=0.2)
                        exact_only = st.checkbox("Exact only: disable inexact matching")
                        exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                        debug_mode = st.checkbox('Turn on debugging for allele calling')
                        do_matching = st.form_submit_button("Run allele calling")

                    if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                        st.info('Analysis in progress')
                    if do_matching:
                        success = generate.generate_targets(exp)
                        if not success:
                            msg = 'Critical: failed to save reference sequences to target file'
                            exp.log(msg)
                            st.error(msg)
                            sleep(0.5)
                        success = generate.generate_primer_assayfams(exp)
                        if not success:
                            msg = 'Critical: failed to save primers and assay families to file'
                            exp.log(msg)
                            st.error(msg)
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
            dc.display_temporary_messages()
            dc.display_persistent_messages(key='main1')

            with tab_col2:
                st.write('')
                show_info_viewer_checkbox()
                
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)
        
        #=============================================== STAGE 7: Reports ==============================================
        if pipeline_stage == 6:
            results_fp = exp.get_exp_fn('results.csv')
            tab_col1, tab_col2 = st.columns([9, 1])
            info_holder = st.container()
            
            
            if not Path(results_fp).exists():
                st.markdown('**No allele calling results available**')
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

                print(f'{len(custom_results)=} {len(rodentity_results)=} {len(other_results)=}')
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

            with tab_col2:
                st.write('')
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        ### End of main display ###

if __name__ == '__main__':
    main()


