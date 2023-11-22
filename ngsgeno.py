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
from stutil import add_vertical_space, custom_text, hline
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


def plate_checklist_expander(available_nimbus, pcr_stage=1):
    """
    Allows the selection/deselection all plates involved in a reaction stage
    """
    exp = st.session_state['experiment']
    included_PCR_plates = set()
    included_taqwater_plates = set()
    checklist_col = st.columns(4)
    pcr_plate_title = "PCR Plates"
    taqwater_plate_title = "Taq/Water Plates"

    if pcr_stage == 1:
        included_DNA_plates = set()
        dna_plate_title = "DNA Plates"
        
        checklist_col[0].markdown(f'**{dna_plate_title}**')
        for nim in available_nimbus:
            echo_filename=Path(nim).stem
            inc_dna = checklist_col[0].checkbox(echo_filename, value=True, key='chk_box_dna_'+nim)
            if inc_dna:
                included_DNA_plates.add(util.guard_pbc(echo_filename.split('_')[-2], silent=True))
    
        checklist_col[1].markdown(f'**{pcr_plate_title}**')
        for pcr_pid in exp.get_pcr_pids():
            inc_pcr = checklist_col[1].checkbox(util.unguard_pbc(pcr_pid, silent=True),\
                            value=True, key='chk_box_pcr_'+pcr_pid)
            if inc_pcr:
                included_PCR_plates.add(util.guard_pbc(pcr_pid, silent=True))

        checklist_col[2].markdown(f'**{taqwater_plate_title}**')
        for taqwater_pid in exp.get_taqwater_pids(pcr_stage):
            inc_taqwater = checklist_col[2].checkbox(util.unguard_pbc(taqwater_pid, silent=True), 
                        value=True, key='chk_box_taqwater_'+taqwater_pid)
            if inc_taqwater:
                included_taqwater_plates.add(util.guard_pbc(taqwater_pid, silent=True))
        
        return included_DNA_plates, included_PCR_plates, included_taqwater_plates

    if pcr_stage == 2:
        included_index_plates = set()
        included_amplicon_plates = set()
        #could make a for loop

        index_plate_title = "Index Plates"
        amplicon_plate_title = "Amplicon Plates"

        checklist_col[0].markdown(f'**{pcr_plate_title}**')
        for pcr_pid in exp.get_pcr_pids():
            inc_pcr = checklist_col[0].checkbox(util.unguard_pbc(pcr_pid, silent=True),
                    value=True, key='chk_box_pcr_'+pcr_pid)
            if inc_pcr:
                included_PCR_plates.add(pcr_pid)

        checklist_col[1].markdown(f'**{taqwater_plate_title}**')
        for taqwater_pid in exp.get_taqwater_pids(pcr_stage):
            inc_taqwater = checklist_col[1].checkbox(util.unguard_pbc(taqwater_pid, silent=True), 
                    value=True, key='chk_box_taqwater_'+taqwater_pid)
            if inc_taqwater:
                included_taqwater_plates.add(taqwater_pid)

        checklist_col[2].markdown(f'**{index_plate_title}**')
        for index_pid in exp.get_index_pids():
            inc_index = checklist_col[2].checkbox(util.unguard_pbc(index_pid, silent=True),
                                                value=True, key='chk_box_index_'+index_pid)
            if inc_index:
                included_index_plates.add(index_pid)
        
        checklist_col[3].markdown(f'**{amplicon_plate_title}**')
        for amplicon_pid in exp.get_amplicon_pids():
            amplicon_index = checklist_col[3].checkbox(util.unguard_pbc(amplicon_pid, silent=True), 
                    value=True, key='chk_box_amplicon_'+amplicon_pid)
            if amplicon_index:
                included_amplicon_plates.add(amplicon_pid)

        return included_PCR_plates, included_taqwater_plates, included_index_plates, included_amplicon_plates
 
    
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
    set_session_state('show_info_viewer', False)
    
    if st.toggle('Info Viewer'):
        st.session_state['show_info_viewer'] = True
    else:
        st.session_state['show_info_viewer'] = False

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

def set_session_state(key, value):
    """
    Initialises a session state with the value if the key is not in session state.
    Args:
        key (str): for the st.session_state dictionary
        value (str or None): value for the key in st.session_state
    """
    if key not in st.session_state:
        st.session_state[key] = value

def create_tabs(tab_data):
    """
    Create tabs from streamlit_extra_components. Assigns ID through enumerating given list.
    Args:
        tab_data (list): list of tuples containg the name and description of each tab
    Returns
        Create tab bar
    """
    return stx.tab_bar(data=[
        stx.TabBarItemData(id=i+1, title=title, description=desc)
        for i, (title, desc) in enumerate(tab_data)
    ], return_type=int)

def add_css():
    #CSS
    st.markdown('''
    <style>
        .stApp [data-testid="stToolbar"]{
            display:none;
        }
        #root > div:nth-child(1) > div > div > div > div > section > div {padding-top: 1rem;}
     </style>
     ''', unsafe_allow_html=True)
    
    #remove drag and drop labels from upload buttons. Class name 'css-9ycgxx' could change in future streamlit versions
    hide_label = """
    <style>
        .css-9ycgxx {
            display: none;
        }
    </style>
    """
    st.markdown(hide_label, unsafe_allow_html=True)

    #css for all form_submit_buttons
    form_button_css = """
    <style>
    div[data-testid="stFormSubmitButton"] button {
        background-color: #4287f5;
        color: white;
        padding: 0.25rem 0.75rem;
        margin: 8px 0;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        font-weight: 400;
        width: fit-content;
        height: auto
    }
    div[data-testid="stFormSubmitButton"] button:hover {
        opacity: 0.8;
        background-color: #cf3276;
        color:white;

    }
    </style>
    """
    st.markdown(form_button_css, unsafe_allow_html=True)

    #css for all buttons with type primary
    primary_button_css = """
    <style>
    button[data-testid="baseButton-primary"] {
        background-color: #4287f5;
        color: white;
        padding: 0.25rem 0.75rem;
        margin: 8px 0;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        font-weight: 400;
        width: fit-content;
        height: auto
    }
    button[data-testid="baseButton-primary"]:hover {
        opacity: 0.8;
    }
    </style>
    """
    st.markdown(primary_button_css, unsafe_allow_html=True)
    #css for all buttons with type secondary
    secondary_button_css = """
    <style>
    button[data-testid="baseButton-secondary"] {
        background-color: #83b3c9;
        color: black;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        font-weight: 400;
    }
    button[data-testid="baseButton-secondary"]:hover {
        opacity: 0.8;
    }
    </style>
    """
    st.markdown(secondary_button_css, unsafe_allow_html=True)
    

def add_new_folder(folder_name):
    """
    *Setting Up*
    Args:
        folder_name (str): name of new folder from user
    Return:
        str or None: An error message if folder exists or there's a fatal path error. None if there's no error.
    """
    add_run_folder_str = 'run_' + folder_name
    exp, msg = create_run_folder(add_run_folder_str)
    print(f'{exp=}, {msg=}')
    if exp:
        exp.save()                                                                      
        st.session_state['experiment'] = exp
        st.experimental_rerun()    
    else:
        if 'already exists' in msg:
            return "Folder name already exists"
        else:
            return "Fatal path error: " + msg
    
def load_run_folder(run_folder):
    """
    *Setting Up*
    Load experiment from existing folder and initilialise pipeline stages.
    Args:
        run_folder (str): Name of existing folder selected to load.
    Return:
        str or None: Error message if could not load experiment or invalid file. None if there's no error.
    """
    if st.session_state['experiment'] == None or st.session_state['experiment'].name != run_folder:
        ch_run_path = 'run_' + run_folder
        if os.path.exists(ch_run_path):
            if not os.path.exists(ch_run_path):
                return "Folder not found: " + ch_run_path
            
            exp = load_experiment(ch_run_path)
            if not exp:
                return "Could not load experiment from: " + ch_run_path
            elif ch_run_path.endswith(exp.name):
                #experiment loaded
                st.session_state['experiment'] = exp
                #Set session state for pipeline stages.
                st.session_state['pipeline_stage'] = 0
                st.session_state.update({key: 1 for key in ['load_tab', 'nimbus_tab', 'primer_tab', 
                                                            'index_tab', 'miseq_tab', 'allele_tab']})
 
                st.experimental_rerun()
            else:
                return "Invalid experiment file in: " + ch_run_path
    
    return None

def set_nimbus_title(exp, efs, nfs):
    """
    *Stage 2: Nimbus*
    Title for nimbus stage
    Args:
        exp (st.session_state['experiment'])
        efs (str): file path to echo files
        nfs (str): file path for nimbus files
    Return
        str or None: title
    """
    #first stage sample files haven't been loaded
    if not st.session_state['experiment'].dest_sample_plates:
        return "Load data inputs to enable Nimbus input file generation."
    else:
        # do we have any Nimbus inputs to generate + download
        echo_files_exist = len(efs) == len(nfs) and len(efs) != 0
        yet_to_run = len(exp.dest_sample_plates) - len(nfs)

        if echo_files_exist:
            return 'All Echo inputs received.'
        if yet_to_run > 0:
            return f'{str(yet_to_run)} 96-well plate sets need Nimbus input file generation' 

def generate_download_buttons(nfs):
    """
    *Stage 2: Nimbus*
    Generates the echo file download buttons
    Args:
        nfs (str): nimbus file paths
    """

    _,dl_col1,dl_col2,dl_col3,dl_col4,_= st.columns([1,9,6,9,6,1])
    
    #print(f"{nfs=} {efs=} {xbcs=}")
    for i,nf in enumerate(nfs):
        nimbus_fn=Path(nf).name

        if (i+1) % 2 != 0:
            with dl_col1:
                custom_text("p", "#4b778c", nimbus_fn, "left")

            dl_col2.download_button("Download ", 
                                    open(nf), 
                                    file_name=nimbus_fn, 
                                    key='nimbus_input'+str(i), 
                                    help=f"Download Nimbus input file {nf}")
    
        else:
            with dl_col3:
                custom_text("p", "#4b778c", nimbus_fn, "left")
        
            dl_col4.download_button("Download ", 
                                    open(nf), file_name=nimbus_fn,\
                                    key='nimbus_input'+str(i), 
                                    help=f"Download Nimbus input file {nf}")

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

def get_echo_picklist_btn_pcr2(exp, PCR_plates, taqwater_plates, index_plates, amplicon_plates):
    """
    *Stage 4: PCR2*
    TODO: Combine this with the PCR1 function
    Display for generate echo pcr 1 picklist
    Args:
        exp (st.session_state['experiment])
        PCR_plates:
        taqwater_plates:
        amplicon_plates
    """
    pcr2_messages = []  # pass by reference
    if not exp.check_ready_pcr2(PCR_plates, taqwater_plates, 
                            index_plates, amplicon_plates, pcr2_messages):
        
        for msg in pcr2_messages:
            st.warning(msg)

    _,button_col,_ = st.columns([2, 2, 1])

    echo_picklist_go = button_col.button('Generate Echo Picklists',\
                                         key='echo_pcr2_go_button', 
                                         type='primary')

    if echo_picklist_go:
        success = run_generate(exp, 
                               exp.generate_echo_PCR2_picklists, 
                               PCR_plates,
                               index_plates, 
                               taqwater_plates, 
                               amplicon_plates)
        if not success:
            st.error('Picklist generation failed. Please see the log')
        else:
            st.session_state['pcr2 picklist'] = True
    

def get_miseq_download_btn(exp):
    """
    *Stage 5: Miseq*
    Args:
        exp (st.session_state['experiment])
    """
    _, miseq_col1, miseq_col2, _ =  st.columns([2,1,1,2])
    for fp in exp.get_miseq_samplesheets():
        
        fp_name = str(Path(fp).name)
        with miseq_col1:
            add_vertical_space(1)
            custom_text('strong', 'black', fp_name, align='right')
        
        with miseq_col2:
            download_miseq = st.download_button(label='Download', 
                                                data=open(fp, 'rt'), 
                                                file_name=fp_name, 
                                                mime='text/csv', 
                                                key='dnld_samplesheet_'+str(fp), 
                                                type='primary')
    

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

    add_css()

    #Initialise experiment
    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')

    set_session_state('folder', None)
    set_session_state('stage', None)

    #Title
    experiment_title = 'Current Experiment: '
    if 'experiment' in st.session_state and st.session_state['experiment'] is not None:
             experiment_title += (st.session_state['experiment'].name)
    else:
             experiment_title += ('None')
    
    #Version (git)
    current_status, current_ver = get_status_output("git describe")

    #Experiment folders
    try:
        existing_run_folders = get_run_folders()
    except Exception as exc:
        print(f'Cannot locate NGSgeno folder {exc}', file=sys.stderr)
        return

    #Header set up
    logo_col, ver_col,_, new_folder_col, create_button_col, ex_folder_col, _ = st.columns([2,2,2,2,1,2,1])
    ver_col.markdown(f'<p style="color:#83b3c9; font-size: 90%"> {current_ver}</p>', unsafe_allow_html=True)
    logo_col.image('ngsg_explorer.png', caption=f'{experiment_title}')

    #Run folder buttons
    new_folder = new_folder_col.text_input('Create new run folder')
    with create_button_col:
        add_vertical_space(2)
    add_new_folder_button = create_button_col.button('Create')
    run_folder = ex_folder_col.selectbox("Select a run folder to open", existing_run_folders)
    
    #Load experiment / create new
    error_msg = None
    if new_folder and add_new_folder_button:
        print(f'{new_folder=}')
        error_msg = add_new_folder(new_folder)
    if run_folder:
        error_msg = load_run_folder(run_folder)
    if error_msg:
        with new_folder_col:
            custom_text("p", "#FF0000", error_msg)

#================================================== START EXPERIMENT ===================================================
    if st.session_state['experiment']:
        exp = st.session_state['experiment']

        #Pipeline stepper bar set up
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

        #Initialise info viewer checkbox
        set_session_state('info_expand', False)
        
        subsection = st.container() 
        message_area = st.container()
        st.session_state['message_area'] = message_area

        #---------------------------------------------- STAGE 1: Load data ---------------------------------------------
        if pipeline_stage == 0:
            _,help_col,_ = subsection.columns([1,4,1]) 
            tab_col1, tab_col2 = subsection.columns([9,1])
            tip_col, _ = st.columns(2)
            info_holder = st.container()
            summary_holder = st.container()
        
            with help_col:
                custom_text(size='p', 
                            color=ld.HELP_COLOUR, 
                            text='Step 1: Upload sample data and assign to a DNA plate. '+\
                                 'Then upload the assay lists, and optionally other files.', 
                            align='left',
                            style="italic")    
            
            with tab_col1:
                load_data_tab = create_tabs([("Load Samples", ""),("Load Consumables", ""), ("Volumes", "")])                   
            if not load_data_tab:
                set_session_state('load_tab', 1)
                load_data_tab = st.session_state['load_tab']
            
            # ******** TAB 1: Load sample data  ***********
            if load_data_tab == 1:
                if unlocked(exp):
                    set_session_state('run queue', [])
                    
                    st.subheader("Upload Sample Files")
                    ld.upload_rodentity_data(key='rodentity_load1')
                    ld.upload_custom_manifests(key='custom_load1')

                with summary_holder:
                    tip_col.info('Upload your sample files, either from Rodentity files or as custom manifests.')
                    st.subheader("Sample Summary")
                    summary = exp.summarise_inputs()
                    if len(summary) > 1:
                        dc.display_samples('load_data_tab1', height=180)
                
                st.session_state['load_tab'] = 1

            # ******** TAB 2: Load consumables, references and assays/primer mappings, edit volumes ********
            if load_data_tab == 2:
                st.info('Upload extra files (consumables). **Ensure that you have at least uploaded '+\
                        'your assay list in order to generate Echo files in Step 2. (Nimbus)**. ')
                
                st.subheader("Summary")
                dc.display_consumables('load_data_tab2')
                if unlocked(exp):
                    set_session_state('upload stage', None)
                    add_vertical_space(1)
                    


                    st.subheader("Upload Consumables")
                    ld.upload_extra_consumables('consumables_load2')
                    ld.upload_pcr1_files('pcr1_load2')
                    ld.upload_pcr2_files('pcr2_load2')
                
                st.session_state['load_tab'] = 2
            
            #******** TAB 3: Editable volumes ********
            if load_data_tab == 3:
                st.info("You can edit the volumes for DNA, primers, indexes, and taq and water, "+
                          "as well as dead volumes, and volume capacity, for each plate type")
                st.subheader("Custom Volumes")
                dc.custom_volumes(exp)
                add_vertical_space(1)
                dc.custom_plate_volumes(exp)

                st.session_state['load_tab'] = 3

            with tab_col2:
                add_vertical_space(1)
                show_info_viewer_checkbox()

            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)
        
        #------------------------------------------------ STAGE 2: Nimbus ----------------------------------------------
        if pipeline_stage == 1:
            exp = st.session_state['experiment']

            _,help_col,_ = subsection.columns([1,4,1])
            tab_col1, tab_col2 = subsection.columns([9,1])
            info_holder = st.container()

            with help_col:
                custom_text(size='p', 
                            color=ld.HELP_COLOUR, 
                            text='Step 2: Generate and download files for the Nimbus. After running these, upload '+\
                                  'the outputs (ie Echo input files.)', 
                            align='left',
                            style="italic")
                
            with tab_col1:
                nimbus_tab = create_tabs([("Download", "Nimbus input files"),("Upload", "Echo input files")])
            if not nimbus_tab:
                set_session_state('nimbus_tab', 1)
                nimbus_tab = st.session_state['nimbus_tab']

            #nimbus fp, echo pb, not barcodes
            nfs, efs, xbcs = exp.get_nimbus_filepaths()

            # ******** TAB 1: Download nimbus ********
            if nimbus_tab == 1:
                if unlocked(exp):
                    st.info('At this Nimbus stage, generate the files to input into the Nimbus. '+\
                            '**For generation, make sure you have at least uploaded sample plates and an assay list**')
                    
                    _, header_col, _ = st.columns([2,2,1])
                    header_col.subheader('Generate Echo Files')

                    #Subtitle
                    nimbus_title = set_nimbus_title(exp, nfs, efs)
                    if nimbus_title:
                        _, title_col,_ = st.columns([1, 4, 1])
                        with title_col:
                            custom_text('h5', '#83b3c9', set_nimbus_title(exp, nfs, efs))
                    add_vertical_space(1)
                    
                    #Generate files button
                    _, btn_col,_ = st.columns([2,2,1])
                    with btn_col:
                        run_gen_nimbus = st.button('Generate Nimbus input files', type="primary")
                    
                    #Generate files
                    if run_gen_nimbus:
                        success = run_generate(exp, exp.generate_nimbus_inputs)    
                        if not success:
                            custom_text("p", "#FF0000", 
                                        "Failed to generate Nimbus files. Please read the log.", 
                                        align="left")
                        else:
                            add_vertical_space(2)
                            nfs, efs, xbcs = exp.get_nimbus_filepaths()
                            generate_download_buttons(nfs)

                st.session_state['nimbus_tab'] = 1

            #******** TAB 2: Upload echo input files ********
            if nimbus_tab == 2:
                if unlocked(exp):
                    st.info('After running the Nimbus, upload the csv files that are required to run the Echo. '+\
                                'These should match the DNA plate barcodes.')
                    _, header_col, _ = st.columns([2,2,1])
                    with header_col:
                        st.subheader('Upload Echo Input Files')

                    ld.upload_echo_inputs('1')

                st.session_state['nimbus_tab'] = 2

            #Info viewer
            with tab_col2:
                add_vertical_space(2)
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #-------------------------------------------- STAGE 3: PCR 1 Primers ------------------------------------------
        if pipeline_stage == 2:
            exp = st.session_state['experiment']
            pcr_stage = 1
            st.session_state['assay_filter'] = True
            #nimbus fp, echo fp, barcodes not in echo
            nfs, efs, xbcs = exp.get_nimbus_filepaths()
            #missing_nims = ['Echo_384_COC_0001_'+util.unguard(xbc, silent=True)+'_0.csv' for xbc in xbcs]

            _,help_col, _ = subsection.columns([1,4,1])
            tab_col1, tab_col2 = subsection.columns([9,1])
            info_holder = st.container()
            tip_col, _ = st.columns([2,1])
            primer_checklist = st.expander('Plates')

            #Help
            with help_col:
                custom_text(size='p', 
                            color=ld.HELP_COLOUR, 
                            text='Step 3: Upload the primer files and add barcodes '+\
                                 'for the first PCR reaction. Then generate and download the picklists.', 
                            align='left',
                            style="italic")
            #Tabs
            with tab_col1:
                primer_tab = create_tabs([("PCR 1", "Components"), ("Generate", "Picklists"), ("Volumes", "")])
            if not primer_tab:
                set_session_state("primer_tab", 1)
                primer_tab = st.session_state['primer_tab']

            with primer_checklist:
                st.subheader('Plate checklist',
                             help='These are the provided plates involved in the reaction. '+\
                                  'You can choose to include / exclude them from the picklist file')
            
            #******** TAB 1: PCR 1 info & adding barcodes and files ********
            if primer_tab == 1:
                if unlocked(exp):
                    tip_col.info('Provide the barcodes for PCR plates and Taq/water plates and '+\
                            'upload primer layouts and volumes to meet the requirements.')
                    
                    #TODO: Checklist doesn't update the components
                    if efs:
                        with primer_checklist:
                            included_DNA_plates, included_PCR_plates, included_taqwater_plates =\
                                        plate_checklist_expander(efs, pcr_stage=pcr_stage)
                            
                        if included_DNA_plates:
                            st.subheader('PCR 1 Components', help='Required plates and volumes for the PCR reaction')
                            dc.display_pcr_components(dna_pids=included_DNA_plates)
                            dc.display_pcr1_components(dna_pids=included_DNA_plates)
                            hline()
                    else:
                        custom_text('h5', '#f63366', "Load Nimbus output files to enable PCR stages")

                    #Barcodes for PCR and taq and water, upload files for primer, custom volumes
                    add_vertical_space(1)
                    st.subheader('Add Barcodes', help='Add barcodes for plates')
                    ld.add_barcodes(key='barcodes_tab1', pcr_stage=pcr_stage)
                    add_vertical_space(1)

                    st.subheader('Upload Files')
                    ld.upload_pcr1_files(key='pcr1_primer1')
                    add_vertical_space(1)

                st.session_state['primer_tab'] = 1
                
            #******** TAB 2: Generate PCR 1 picklists ********
            if primer_tab == 2:
                if unlocked(exp):
                    set_session_state('pcr1 picklist', False)

                    tip_col.info('Ensure you have upload the primer layouts and volumes and provided plate barcodes '+\
                                  'before generating the picklist.')

                    if not efs:
                        st.warning('No DNA plate information available. Have you uploaded Echo input files yet?')
                    else:
                        with primer_checklist:
                            included_DNA_plates,included_PCR_plates,\
                                    included_taqwater_plates = plate_checklist_expander(efs, pcr_stage=pcr_stage)

                            st.session_state['included_DNA_pids'] = included_DNA_plates
                        
                        #show generate button and run generate
                        if included_DNA_plates:
                            get_echo_picklist_btn_pcr1(exp, included_DNA_plates, included_PCR_plates,\
                                                                included_taqwater_plates)
                #show download buttons            
                if st.session_state['pcr1 picklist']:
                    dc.get_echo1_downloads_btns()

                st.session_state['primer_tab'] = 2
            
            #******** TAB 3: Editable volumes ********
            if primer_tab == 3:
                st.info("You can edit the volumes for DNA, primers, indexes, and taq and water, "+
                                "as well as dead volumes, and volume capacity, for each plate type")
                st.subheader("Custom Volumes")
                dc.custom_volumes(exp)
                add_vertical_space(1)
                dc.custom_plate_volumes(exp)
                
                st.session_state['primer_tab'] = 3
            
            #info viewer
            with tab_col2:
                add_vertical_space(2)
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #--------------------------------------------- STAGE 4: PCR 2 Index --------------------------------------------
        if pipeline_stage == 3:
            exp = st.session_state['experiment']
            pcr_stage = 2

            #nimbus fp, echo fp, barcodes not in echo
            nfs, efs, xbcs = exp.get_nimbus_filepaths()
            available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]
            #missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]

            #Container set up
            _,help_col, _ = subsection.columns([1,4,1])
            tab_col1, tab_col2 = subsection.columns([9,1])
            info_holder = st.container()
            tip_col, _ = st.columns([2,1])
            index_checklist = st.expander('Plates')
            title_holder = st.empty()
            pcr_comp_holder = st.container()

            with help_col:
                custom_text(size='p', 
                            color=ld.HELP_COLOUR, 
                            text='Step 4: Upload the index files and add barcodes '+\
                                 'for the second PCR reaction. Then generate and download the picklists.', 
                            align='left',
                            style="italic")
            
            with tab_col1:
                index_tab = create_tabs([("PCR 2", "Components"), ("Generate", "Picklists"), ("Volumes", "")])
            if not index_tab:
                set_session_state('index_tab', 1)
                index_tab = st.session_state['index_tab']
            
            with index_checklist:
                st.subheader('Plate checklist',
                             help='These are the provided plates involved in the reaction. '+\
                                  'You can choose to include / exclude them from the picklist file')

            #******** TAB 1: PCR 2 Components ********
            if index_tab == 1:
                if unlocked(exp):
                    tip_col.info('Provide the barcodes for Taq/water plates and '+\
                                  'upload index layouts and volumes to meet the requirements.')
                    
                    #Provide barcodes for pcr & taq/water, upload index files, adjust volumes
                    st.subheader('Add Barcodes', help='Add barcodes for plates')
                    ld.add_barcodes(key='index_barcodes', pcr_stage=pcr_stage)
                    add_vertical_space(1)

                    st.subheader('Upload Files')
                    ld.upload_pcr2_files(key='pcr2_index1')
                    add_vertical_space(1)

                    st.session_state['index_tab'] = 1
                    
                    if available_nimbus:
                        with index_checklist:
                            included_PCR_plates, included_taqwater_plates,included_index_plates, included_amplicon_plates =\
                                            plate_checklist_expander(available_nimbus, pcr_stage)
                            hline()
                            
                        with pcr_comp_holder:
                            st.subheader('PCR 2 Components')
                            dc.display_pcr_components()
                            dc.display_pcr2_components(pcr_pids=included_PCR_plates, 
                                                       amplicon_pids=included_amplicon_plates)
                            hline()
                            add_vertical_space(1)
                    else:
                        with title_holder:
                            custom_text('h5', '#f63366', "Load Nimbus output files to enable PCR stages")
                        
                st.session_state['index_tab'] = 1

            #******** TAB 2: Generate PCR 2 picklists ********
            if index_tab == 2:
                if unlocked(exp):
                    set_session_state('pcr2 picklist', False)

                    if not available_nimbus:
                        st.warning('No DNA plate information available. Have you uploaded Echo input files yet?')
                    else:
                        with index_checklist:
                            included_PCR_plates, included_taqwater_plates, included_index_plates,\
                                    included_amplicon_plates = plate_checklist_expander(available_nimbus, pcr_stage)
                            hline()
                            
                        get_echo_picklist_btn_pcr2(exp, 
                                                   included_PCR_plates, 
                                                   included_taqwater_plates, 
                                                   included_index_plates, 
                                                   included_amplicon_plates)
                        
                if st.session_state['pcr2 picklist']:
                    dc.get_echo2_download_btns()
                
                st.session_state['index_tab'] = 2

            #******** TAB 3: Editable volumes ********
            if index_tab == 3:
                st.info("You can edit the volumes for DNA, primers, indexes, and taq and water, "+
                                "as well as dead volumes, and volume capacity, for each plate type")
                st.subheader("Custom Volumes")
                dc.custom_volumes(exp)
                add_vertical_space(1)
                dc.custom_plate_volumes(exp)

                st.session_state['index_tab'] = 3        
                
            #Info viewer
            with tab_col2:
                add_vertical_space(2)
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #----------------------------------------------- STAGE 5: Miseq ------------------------------------------------
        if pipeline_stage == 4:
            exp = st.session_state['experiment']
            _,help_col, _ = subsection.columns([1,7,1])
            tab_col1, tab_col2 = subsection.columns([9,1])
            info_holder = st.container()
            tip_col, _ = st.columns([2,1])
            
            with help_col:
                custom_text(size='p', 
                            color=ld.HELP_COLOUR, 
                            text='Step 5: Download the Miseq file. Upload reference sequences if you haven\'t already. '+\
                                 'Once you get the Miseq sequences, you can upload them here as well.',
                            align='left',
                            style="italic")
                add_vertical_space(1)

            with tab_col1:
                miseq_tab = create_tabs([("Download", "Miseq Samplesheet"), ("Upload", "Miseq Sequence Files")])
            if not miseq_tab:
                set_session_state('miseq_tab', 1)
                miseq_tab = st.session_state['miseq_tab']

            add_vertical_space(1)

            #******* TAB 1: Download Miseq Samplesheet. *******
            if miseq_tab == 1:
                if exp.locked:
                    st.warning(f'Experiment {exp.name} locked from further modification')
                
                st.subheader('Download File')
                if exp.get_miseq_samplesheets():
                    get_miseq_download_btn(exp)
                    add_vertical_space(4)
                    hline()
                else:
                    st.warning(f'No MiSeq Samplesheet available for download')

                st.session_state['miseq_tab'] = 1

            #******* TAB 2: Upload custom reference file and sequence files. *******
            if miseq_tab == 2:

                st.subheader('Upload Custom Reference Files')
                ld.upload_reference_sequences('reference_miseq')
                add_vertical_space(2)

                ready_messages = []
                if not exp.check_sequence_upload_ready(ready_messages):
                    for msg in ready_messages:
                        st.error(msg)

                    st.warning('These resources are required for allele calling'+\
                               'and must be present before FASTQs can be uploaded')
                else:
                    ld.upload_miseq_fastqs()

                st.session_state['miseq_tab'] = 2

            #Info viewer
            with tab_col2:
                add_vertical_space(2)
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        #----------------------------------------- STAGE 6: Allele Calling ---------------------------------------------
        if pipeline_stage == 5:
            exp = st.session_state['experiment']
            _,help_col, _ = subsection.columns([1,4,1])
            tab_col1, tab_col2 = subsection.columns([9,1])
            info_holder = st.container()

            with help_col:
                custom_text(size='p', 
                            color=ld.HELP_COLOUR, 
                            text="Step 6: If you haven\'t uploaded reference sequences, or added FASTQ files, "+\
                                 "load them now. Then you can run the allele calling.",
                            align='left',
                            style="italic")
                add_vertical_space(1)

            with tab_col1:
                allele_tab = create_tabs([("Allele Calling", "")])
            if not allele_tab:
                set_session_state('allele_tab', 1)
                allele_tab = st.session_state['allele_tab']

            # Only offer upload in the Miseq pipeline section
            if allele_tab == 1:
                rundir = exp.get_exp_dn()
                seq_ready_messages = []
                call_ready_messages = []
                if not exp.check_sequence_upload_ready(seq_ready_messages):
                    for msg in seq_ready_messages:
                        st.error(msg)
                    st.warning('These resources are required for allele calling and must be present before FASTQs can be uploaded')
                    add_vertical_space(1)
                    st.subheader('Upload Custom Reference File')
                    ld.upload_reference_sequences('reference_allele1')
                    add_vertical_space(1)

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
                        exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                        nocache = st.checkbox("Disable caching: slower but removes any target ambiguity")
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
                add_vertical_space(1)
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
                add_vertical_space(1)
                show_info_viewer_checkbox()
            with info_holder:
                if st.session_state['show_info_viewer']:
                    dc.info_viewer(1)

        st.session_state['pipeline_stage'] = pipeline_stage

if __name__ == '__main__':
    main()



