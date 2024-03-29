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


def show_info_viewer_checkbox():
    """
    Allows the user to turn the info viewer panel on and off
    """
    if 'show_info_viewer' not in st.session_state:
        st.session_state['show_info_viewer'] = False
    
    if st.checkbox('Info Viewer'):
        st.session_state['show_info_viewer'] = True
    else:
        st.session_state['show_info_viewer'] = False


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

    st.markdown('''
    <style>
        .stApp [data-testid="stToolbar"]{
            display:none;
        }
        #root > div:nth-child(1) > div > div > div > div > section > div {padding-top: 1rem;}
     </style>
     ''', unsafe_allow_html=True)

   

    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')

    if 'folder' not in st.session_state:
        st.session_state['folder'] = None
    #Folder inputs
    
    #remove drag and drop labels from upload buttons. Class name 'css-9ycgxx' could change in future streamlit versions
    hide_label = """
    <style>
        .css-9ycgxx {
            display: none;
        }
    </style>
    """
    st.markdown(hide_label, unsafe_allow_html=True)

    logo_col, ver_col,_, new_folder_col, create_button_col, ex_folder_col, _ = st.columns([2,2,2,2,1,2,1])
    # current_experiment_col, _ = st.columns(2)
    # with title_column:
    #     st.markdown('<h3 style="align:center; color:#2040a1">NGS Genotyping</h3>', unsafe_allow_html=True)
    #     
    #     if 'experiment' in st.session_state and st.session_state['experiment'] is not None:
    #         experiment_title += (st.session_state['experiment'].name)
    #     else:
    #         experiment_title += ('None')

    #     st.markdown(f'<p style="color:#83b3c9">{experiment_title}</p>', unsafe_allow_html=True)
    experiment_title = 'Current Experiment: '
    if 'experiment' in st.session_state and st.session_state['experiment'] is not None:
             experiment_title += (st.session_state['experiment'].name)
    else:
             experiment_title += ('None')

    current_status, current_ver = get_status_output("git describe")
    ver_col.markdown(f'<p style="color:#83b3c9; font-size: 90%"> {current_ver}</p>', unsafe_allow_html=True)

    logo_col.image('ngsg_explorer.png', caption=f'{experiment_title}')



    add_run_folder = new_folder_col.text_input('Create new run folder')
    create_button_col.write('')
    create_button_col.write('')
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

    if 'stage' not in st.session_state:
        st.session_state['stage'] = None

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

    if st.session_state['experiment']:
        exp = st.session_state['experiment']
        pipeline_stages=["Load", "Nimbus", "Primers", "Index", "Miseq", "Alleles", "Reports"]
        pipeline_stage = stx.stepper_bar(steps=pipeline_stages, lock_sequence=False)
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

        if 'info_expand' not in st.session_state:
            st.session_state['info_expand'] = False
        
        subsection = st.container()
        message_area = st.container()
        st.session_state['message_area'] = message_area

        #Load data
        if pipeline_stage == 0:
            tab_col1, tab_col2 = subsection.columns([9,1])
            with tab_col1:
                load_data_tab = stx.tab_bar(data=[
                    stx.TabBarItemData(id=1, title="Load Samples", description=""),
                    stx.TabBarItemData(id=2, title="Load Consumables", description="")
                    #stx.TabBarItemData(id=3, title="View Data", description="")
                ], return_type=int)                            
            #dc.info_viewer()  
            if not load_data_tab:
                if 'load_tab' not in st.session_state:
                    st.session_state['load_tab'] = 1
                load_data_tab = st.session_state['load_tab']

            info_holder = st.container()
            
            # load sample data
            if load_data_tab == 1:
                summary_holder = st.container()
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    if 'run queue' not in st.session_state:
                        st.session_state['run queue'] = []
                    ld.load_rodentity_data('rodentity_load1')
                    ld.load_custom_manifests('custom_load1')
                    ld.load_amplicons('amp_load1')
                st.session_state['load_tab'] = 1
                with summary_holder:
                    summary = exp.summarise_inputs()
                    if len(summary) > 1:
                        dc.display_samples('load_data_tab1', height=180)

            #load consumables, references and assays/primer mappings
            if load_data_tab == 2:
                dc.display_consumables('load_data_tab2')
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    if 'upload stage' not in st.session_state:
                        st.session_state['upload stage'] = None 
                    ld.custom_volumes(exp)
                    ld.upload_extra_consumables('consumables_load2')
                    ld.upload_pcr1_files('pcr1_load2')
                    ld.upload_pcr2_files('pcr2_load2')
                st.session_state['load_tab'] = 2

            with tab_col2:
                st.write('')
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)
        
        #Nimbus
        if pipeline_stage == 1:
            nim_tab1_err = ''
            nimbus_title = ''

            tab_col1, tab_col2 = st.columns([9,1])
            with tab_col1:
                nimbus_tab = stx.tab_bar(data=[
                    stx.TabBarItemData(id=1, title="Download", description="Nimbus input files"),
                    stx.TabBarItemData(id=2, title="Upload", description="Echo input files")
                    #stx.TabBarItemData(id=3, title="View Data", description="")
                ], return_type=int)
            info_holder = st.container()
            
            
            if not nimbus_tab:
                if 'nimbus_tab' not in st.session_state:
                    st.session_state['nimbus_tab'] = 1
                nimbus_tab = st.session_state['nimbus_tab']

            exp = st.session_state['experiment']
            nfs, efs, xbcs = exp.get_nimbus_filepaths()

            #download nimbus
            if nimbus_tab == 1:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    if not st.session_state['experiment'].dest_sample_plates:
                        nimbus_title = "Load data inputs to enable Nimbus input file generation."
                    else:
                        # do we have any Nimbus inputs to generate + download
                        yet_to_run = len(exp.dest_sample_plates) - len(nfs)
                        if len(efs) == len(nfs) and len(efs) != 0:  # already exists
                            nimbus_title = 'All Echo inputs received.'
                        if yet_to_run > 0: 
                            nimbus_title += ' ' + str(yet_to_run) + " 96-well plate sets need Nimbus input file generation"

                            plates_to_run = [dest_plate for dest_plate in exp.dest_sample_plates\
                                         if all([dest_plate not in nf for nf in nfs])]
                            plates_to_run_str = '\n'.join(plates_to_run)
                                
                    st.markdown(f'<h5 style="text-align:center;color:#f63366">{nimbus_title}</h5>',\
                             unsafe_allow_html=True)
                    st.write('')

                    run_gen_nimbus = st.button('Generate Nimbus input files')
                    if run_gen_nimbus:
                        success = run_generate(exp, exp.generate_nimbus_inputs)
                        if not success:
                            nim_tab1_err = "Failed to generate Nimbus files. Please read the log."
                        else:
                            nfs, efs, xbcs = exp.get_nimbus_filepaths()

                    if nim_tab1_err:
                        st.write(f'<p style="color:#FF0000">{nim_tab1_err}</p>', unsafe_allow_html=True)

                _,dl_col1,dl_col2,dl_col3,dl_col4,_= st.columns([1,9,6,9,6,1])
                
                #Generate file names to download + download buttons
                #print(f"{nfs=} {efs=} {xbcs=}")
                for i,nf in enumerate(nfs):
                    nimbus_fn=Path(nf).name

                    if (i+1) % 2 != 0:
                        dl_col1.markdown('<p style="text-align:left;color:#4b778c;padding:5px">'+
                                        f'{nimbus_fn}</p>',
                                        unsafe_allow_html=True)

                        dl_col2.download_button("Download ", 
                                                open(nf), 
                                                file_name=nimbus_fn, 
                                                key='nimbus_input'+str(i), 
                                                help=f"Download Nimbus input file {nf}")
                
                    else:
                        dl_col3.markdown('<p style="text-align:left;color:#4b778c;padding:5px">'+
                                        f'{nimbus_fn}</p>',
                                        unsafe_allow_html=True)
                 
                        dl_col4.download_button("Download ", 
                                                open(nf), file_name=nimbus_fn,\
                                                key='nimbus_input'+str(i), 
                                                help=f"Download Nimbus input file {nf}") 
                st.session_state['nimbus_tab'] = 1

            #Upload nimbus
            if nimbus_tab == 2:
                ld.upload_echo_inputs('1')                
                st.session_state['nimbus_tab'] = 2

            with tab_col2:
                st.write('')
                st.write('')
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #Primer PCR
        if pipeline_stage == 2:

            exp = st.session_state['experiment']
            st.session_state['assay_filter'] = True
            pcr_stage = 1

            #set up tabs
            tab_col1, tab_col2 = st.columns([9,1])
            with tab_col1:
                primer_tab = stx.tab_bar(data=[
                    stx.TabBarItemData(id=1, title="PCR 1", description="Components"),
                    #stx.TabBarItemData(id=2, title="Provide", description="Plates"),
                    stx.TabBarItemData(id=2, title="Generate", description="Picklists")
                ], return_type=int)
            info_holder = st.container()
            
            if not primer_tab:
                if 'primer_tab' not in st.session_state:
                    st.session_state['primer_tab'] = 1
                primer_tab = st.session_state['primer_tab']
            
            #nimbus fp, echo fp, barcodesnot in echo
            nfs, efs, xbcs = exp.get_nimbus_filepaths()
            missing_nims = ['Echo_384_COC_0001_'+util.unguard(xbc, silent=True)+'_0.csv' for xbc in xbcs]

            #PCR 1 info & adding in barcodes / files
            if primer_tab == 1:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    st.write('Plate checklist')
                    primer_checklist = st.container()
                    st.write('***')
                    if efs:
                        with primer_checklist:
                            included_DNA_plates, included_PCR_plates, included_taqwater_plates =\
                                        dc.plate_checklist_pcr1(exp)
                            
                        if included_DNA_plates:
                            dc.display_pcr_components(dna_pids=included_DNA_plates, pcr_pids=included_PCR_plates)
                            dc.display_pcr1_components(dna_pids=included_DNA_plates)
                    else:
                        no_nimbus_msg = "Load Nimbus output files to enable PCR stages"
                        st.markdown(f'<h5 style="text-align:center;color:#f63366">{no_nimbus_msg}</h5',\
                                unsafe_allow_html=True)

                    #barcodes for PCR and taq and water, upload files for primer, adjust volumes
                    ld.provide_barcodes('barcodes_tab1', pcr_stage=pcr_stage)
                    ld.upload_pcr1_files('pcr1_primer1')
                    ld.custom_volumes(exp)

                st.session_state['primer_tab'] = 1
                

            #generate PCR 1 picklists
            if primer_tab == 2:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    st.write('Plate checklist')
                    primer_checklist_exp = st.container()
                    st.write('***')
                    if not efs:
                        st.warning('No DNA plate information available. Have you uploaded Echo input files yet?')
                    else:
                        with primer_checklist_exp:
                            included_DNA_plates, included_PCR_plates, included_taqwater_plates =\
                                    dc.plate_checklist_pcr1(exp)
                            st.session_state['included_DNA_pids'] = included_DNA_plates
                        
                        if included_DNA_plates:
                            pcr1_messages = []
                            if not exp.check_ready_pcr1(included_DNA_plates,\
                                        included_PCR_plates, included_taqwater_plates, pcr1_messages):
                                for msg in pcr1_messages:
                                    st.warning(msg)
                            else:
                                _,picklist_button_col,_ = st.columns([2, 2, 1])
                                picklist_button_col.write('')
                                echo_picklist_go = picklist_button_col.button('Generate Echo Picklist',
                                        key='echo_pcr1_go_button')
                                if echo_picklist_go:
                                    success = run_generate(exp, exp.generate_echo_PCR1_picklists, 
                                            included_DNA_plates, included_PCR_plates, included_taqwater_plates)
                                    if not success:
                                        st.write('Picklist generation failed. Please see the log')
                        
                dc.show_echo1_outputs() 
                st.session_state['primer_tab'] = 2
            
            with tab_col2:
                st.write('')
                st.write('')
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #Index PCR
        if pipeline_stage == 3:

            exp = st.session_state['experiment']
            tab_col1, tab_col2 = st.columns([9,1])
            pcr_stage = 2

            #tab setup
            with tab_col1:
                index_tab = stx.tab_bar(data=[
                    stx.TabBarItemData(id=1, title="PCR 2", description="Components"),
                    #stx.TabBarItemData(id=2, title="Provide", description="Plates"),
                    stx.TabBarItemData(id=2, title="Generate", description="Picklists")
                ], return_type=int)
            info_holder = st.container()
            
            if not index_tab:
                if 'index_tab' not in st.session_state:
                    st.session_state['index_tab'] = 1
                index_tab = st.session_state['index_tab']

            #PCR 2 components
            if index_tab == 1:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    st.write('Plate checklist')
                    index_checklist = st.container()
                    st.write('***')
                    title_holder = st.empty()
                    pcr_comp_holder = st.container()

                    #provide barcodes for pcr & taq/water, upload index files, adjust volumes
                    ld.provide_barcodes('index_barcodes', 2)
                    ld.upload_pcr2_files('pcr2_index1')
                    ld.custom_volumes(exp)

                    st.session_state['index_tab'] = 1
                    with index_checklist:
                        included_PCR_plates, included_taqwater_plates, included_index_plates, \
                                included_amplicon_plates = dc.plate_checklist_pcr2(exp)
                    with pcr_comp_holder:
                        dc.display_pcr_components(pcr_pids=included_PCR_plates, amplicon_pids=included_amplicon_plates)
                        dc.display_pcr2_components(pcr_pids=included_PCR_plates,amplicon_pids=included_amplicon_plates)
                    
                    if not included_PCR_plates and not included_amplicon_plates:
                        no_nimbus_msg = "Load Nimbus output files to enable PCR stages"
                        title_holder.markdown(f'<h5 style="text-align:center;color:#f63366">{no_nimbus_msg}</h5',\
                                unsafe_allow_html=True)
                    
                st.session_state['index_tab'] = 1

            #generate PCR 2 picklist
            if index_tab == 2:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                else:
                    st.write('Plate checklist')
                    index_checklist = st.container()
                    st.write('***')
                    with index_checklist:
                        do_generate = False
                        included_PCR_plates, included_taqwater_plates, included_index_plates,\
                                included_amplicon_plates =\
                                dc.plate_checklist_pcr2(exp)
                        pcr2_messages = []  # pass by reference
                        if included_PCR_plates:  # standard run
                            if not exp.check_ready_pcr2(included_PCR_plates, included_taqwater_plates, 
                                    included_index_plates, included_amplicon_plates, pcr2_messages):
                                for msg in pcr2_messages:
                                    st.warning(msg)
                            else:
                                do_generate = True
                        else:
                            if included_amplicon_plates:
                                if not exp.check_ready_pcr2_amplicon_only(included_taqwater_plates,
                                        included_index_plates, included_amplicon_plates, pcr2_messages):
                                    for msg in pcr2_messages:
                                        st.warning(msg)
                                else:
                                    do_generate = True
                                    
                        if do_generate:
                            show_generate = False
                            if 'amplicon_only' in st.session_state:
                                show_generate = st.session_state['amplicon_only']
                            if included_PCR_plates:
                                show_generate = True

                            if included_amplicon_plates and not included_PCR_plates:
                                _, amp1, amp2, amp3, _ = st.columns([3, 3, 1, 1, 3])
                                amp1.warning('Is this an amplicon only run?')
                                yes_amplicon = amp2.button('Yes')
                                no_amplicon = amp3.button('No')

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
                            
                                    success = run_generate(exp, exp.generate_echo_PCR2_picklists, included_PCR_plates,
                                            included_index_plates, included_taqwater_plates, included_amplicon_plates)
                                    if not success:
                                        st.write('Picklist generation failed. Please see the log')
                        

                dc.show_echo2_outputs()
                st.session_state['index_tab'] = 2            
                
            with tab_col2:
                st.write('')
                st.write('')
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #Miseq
        if pipeline_stage == 4:
            tab_col1, tab_col2 = st.columns([9,1])
            with tab_col1:
                miseq_tab = stx.tab_bar(data=[
                    stx.TabBarItemData(id=1, title="Download", description="Miseq Samplesheet"),
                    stx.TabBarItemData(id=2, title="Upload", description="Miseq Sequence Files"),
                ], return_type=int)
            info_holder = st.container()
            
   
            if not miseq_tab:
                if 'miseq_tab' not in st.session_state:
                    st.session_state['miseq_tab'] = 1
                miseq_tab = st.session_state['miseq_tab']

            st.write('')
            exp = st.session_state['experiment']
            if miseq_tab == 1:
                _, miseq_col1, miseq_col2, _ =  st.columns([2,1,1,2])
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                ld.upload_reference_sequences('reference_miseq1')
                if exp.get_miseq_samplesheets():
                    for fp in exp.get_miseq_samplesheets():
                        fp_name = str(Path(fp).name)
                        miseq_col1.markdown(f'<strong style="color:#486e7a">{fp_name}</strong>', unsafe_allow_html=True)
                    
                        download_miseq = miseq_col2.download_button(label='Download', data=open(fp, 'rt'), 
                                file_name=fp_name, mime='text/csv', key='dnld_samplesheet_'+str(fp))
                else:
                    st.warning(f'No MiSeq Samplesheet available for download')
                st.session_state['miseq_tab'] = 1

            if miseq_tab == 2:
                ld.upload_reference_sequences('reference_miseq2')
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

        #Allele Calling
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
            
            
            if not allele_tab:
                if 'allele_tab' not in st.session_state:
                    st.session_state['allele_tab'] = 1
                allele_tab = st.session_state['allele_tab']

            # Only offer upload in the Miseq pipeline section
            if allele_tab == 1:
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
                    completion_msg = st.empty()
                    match_prog = st.progress(0)
                    asyncio.run(report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog))
                    
                st.session_state['allele_tab'] = 1

            with tab_col2:
                st.write('')
                show_info_viewer_checkbox()
                
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)
        
        # Reports
        if pipeline_stage == 6:
            exp = st.session_state['experiment']
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

        st.session_state['pipeline_stage'] = pipeline_stage

if __name__ == '__main__':
    main()


