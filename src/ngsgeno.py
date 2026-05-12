#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__version__ = "2.04.000"

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Core interface script. Run with: streamlit --run ngsgeno.py

Needs load_data.py for GUI functions that are responsible for incorporating data into an
experiment, and display_components.py for functions dedicated to the presentation of GUI
elements
"""
from threading import Thread
from re import S
import re
import select
#from telnetlib import theNULL
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
from copy import deepcopy
import weakref
from io import StringIO

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components

from stutil import add_vertical_space, custom_text, hline, init_state, \
        upper_info, upper_height, lower_info, lower_height, m, mq, set_state, add_css

from experiment import Experiment, EXP_FN, load_experiment
import util
import transaction as trans
import generate
import parse
import match
import ngsmatch
import load_data as ld
import display_components as dc
import info_viewer as iv

from stage_load import stage_load
from stage_nimbus import stage_nimbus
from stage_primers import stage_primers
from stage_index import stage_index
from stage_miseq import stage_miseq
from stage_alleles import stage_alleles
from stage_reports import stage_reports

import extra_streamlit_components as stx

import asyncio
from time import sleep

global unsaved_exp
if 'experiment' in st.session_state and st.session_state['experiment']:
    unsaved_exp = st.session_state['experiment']
else:
    unsaved_exp = None


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
    run_folders = [d[4:] for d in os.listdir('.') if d.startswith('run_') and os.path.isdir(d)]
    #sort by date modified
    sorted_folders = sorted(run_folders, key=lambda d: os.path.getmtime(f'run_{d}'), reverse=True)
    sorted_folders.insert(0, '')
    return sorted_folders


def create_run_folder(newpath):
    """ returns Experiment or None, and a message string. Requires a full path and generates a new experiment file """
    m(f'Attempting to create: {newpath}', level='begin')
    if not os.path.exists(newpath):
        try:
            os.mkdir(newpath)
        except Exception as exc:
            m(f'Could not create new folder: {newpath} {exc}', level='error')
            return None, 'Failed to create new folder: ' + newpath
    else:
        if os.path.exists(os.path.join(newpath, EXP_FN)):
            return None, f'Experiment already exists with this name: {newpath}'
    print(f'Generating experiment: {newpath}')
    if 'experiment' in st.session_state and st.session_state['experiment']:
        st.session_state['experiment'].save()
    exp = Experiment(name=newpath[4:])  # chop off 'run_'
    return exp, ''





def report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog):
    """
    Allows the interface to keep running while a background process (ngsmatch.py) progress is tracked
    """
    launch_progress = 0
    match_progress = 0
    while True:    
        progress_files = list(Path(rundir).glob('match_progress_*'))
        if len(progress_files) == 1:
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
            sleep(0.3)
        if match_progress == 100:
            launch_msg.write('Allele calling task launch progress: Done')
            completion_msg.write('Allele calling task completion progress: Done')
            if 'matching_in_progress' in st.session_state:
                st.session_state['matching_in_progress'] = None
            return


def load_experiment_screen():
    """
    Landing screen
    """
    global unsaved_exp
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
            st.rerun()
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
                    st.rerun()
                else:
                    error_msg = "Invalid experiment file in: " + ch_run_path
    new_folder_col.markdown(f'<p style="color:#FF0000; text-align:center">{error_msg}</p>',\
            unsafe_allow_html=True)


def save_button(exp, key):
    """
    Save button for users
    """
    save = st.button('💾', type = 'primary', help = 'Save current experiment', key = key)
    if save:
        try:
            print(f'Saving experiment {exp.name}', file=sys.stderr, flush=True)
            exp.save()
        except Exception as exc:
            m(f'Saving experiment failed! {exc}', level='error')


def home_button(exp):
    unload_button = st.button('🏠', type='primary', help='Go back and change experiment')
    if unload_button:
        if 'experiment' in st.session_state and st.session_state['experiment']:
            st.session_state['experiment'].save()
        st.session_state['experiment'] = None
        st.rerun()


def save_message(exp, key):
    _, col1, col2 = st.columns([9, 3, 1])
    col1.info('**Remember to save**')
    with col2:
        save_button(exp, key)


def var_select_cb():
    """
    Callback for variant selection
    """
    st.session_state['v_chosen'] = st.session_state.get('variant_select1', None)
    if st.session_state['v_chosen'] is not None:
        st.session_state['ref_chosen'] = st.session_state['v_chosen'].split('//')[0]


def ref_select_cb():
    """
    Callback for reference selection
    """
    st.session_state['ref_chosen'] = st.session_state.get('amplicon_select1', None)


async def report_match_progress(exp, preprocess_prog, match_prog):
    """
    Allows the interface to keep running while a background process (ngsmatch.py) progress is tracked
    """
    if 'match_manager' in st.session_state:
        MM = st.session_state['match_manager']
        while True:
            preprocess_prog.progress(MM.get_preprocess_progress())
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
            preprocess_prog.write('Allele calling task launch progress: '+str(launch_progress)+'%')
            preprocess_prog.progress(launch_progress)
            match_prog.write('Allele calling task completion progress: '+str(match_progress)+'%')
            match_prog.progress(match_progress)

            if launch_progress == 100 and match_progress == 100:
                m('Analysis completed', level='info')
                st.session_state['match_running'] = False
                return
            await asyncio.sleep(1)
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
            m('Analysis completed', level='info')
            st.session_state['match_running'] = False
            return
        await asyncio.sleep(1)

def manage_matching(exp: Experiment, lock_path: Path, targets: Path, primer_assayfam: Path, outfn: Path,
        variants: Path, ncpus: int, mincov: int, minprop: float, exhaustive: bool, debug: bool) -> None:
    """
    Manage the NGS match process
    """
    if os.path.exists(lock_path):
        MM = st.session_state.get('match_manager', None)
        if not MM:
            m('No current match manager found, clearing lock file', level='warning')
            try:
                os.remove(lock_path)
            except Exception as exc:
                m(f'Removing lock file: {exc}', level='error')
            return


        # report on progress - NGSMatch object already exists
        print("Analysis already running", file=sys.stderr)
        exit(2)
    else:
        try:
            with open(lock_path,"wt"):
                rundir = exp.get_exp_dn()
                MM = MatchManager(rundir, targets, primer_assayfam, outfn, variants, ncpus, mincov, minprop,
                        exhaustive, debug)
                init_state('match_manager', MM)
                report_progress(rundir, 0, 0)  # set this up asap
                main(args)
            print('Completed regular')
        except Exception as exc:
            m(f'Error in NGSMatch: {exc}', level='error')

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

    home_col, save_col, _ = info_col.columns(3)
    with home_col:
        home_button(exp)
    with save_col:
        save_button(exp, key = 'logo')

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

    print(os.getcwd())

    add_css()

    init_state('experiment', None)

    if 'experiment' not in st.session_state or st.session_state['experiment'] is None:
        load_experiment_screen()

    #================================================ START EXPERIMENT =================================================
    else:  # main program
        exp = st.session_state['experiment']
        if not hasattr(exp, '_finalizer'):
            exp._finalizer = weakref.finalize(exp, exp.save)

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

        init_state('info_expand', False)
        upper_container = st.container()
        message_container = st.container()
        main_body_container = st.container()
        add_vertical_space(4)
        save_container = st.container()
        add_vertical_space(2)
        #hline()
        lower_container = st.container(border = True)

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
        

        # attempt to parse any files that are set for upload
        parse.process_upload_queue(exp)
        if '_upload_pending' not in exp.uploaded_files:
            exp.uploaded_files['_upload_pending'] = {}

        with save_container:
            save_message(exp, key = 'save1')

        with lower_container:
            success = iv.info_selection("bottom_viewer", 'info_panel3', 'info_panel4',
                    'lower_panel_height', default_view1='Log', default_view2='None',
                    default_height=st.session_state.get('lower_panel_height',350))


        #============================================== STAGE 1: Load data =============================================
        if pipeline_stage == 0:
            stage_load(exp, main_body_container, upper_container, message_container)
            
        #=============================================== STAGE 2: Nimbus ===============================================
        if pipeline_stage == 1:
            stage_nimbus(exp, main_body_container, upper_container, message_container)
            
        #=========================================== STAGE 3: PCR 1 Primers ============================================
        if pipeline_stage == 2:
            stage_primers(exp, main_body_container, upper_container, message_container)
            
        #============================================ STAGE 4: PCR 2 Index =============================================
        if pipeline_stage == 3:
            stage_index(exp, main_body_container, upper_container, message_container)
            
        #=============================================== STAGE 5: Miseq ================================================
        if pipeline_stage == 4:
            stage_miseq(exp, main_body_container, upper_container, message_container)
         
        #=========================================== STAGE 6: Allele Calling ===========================================
        if pipeline_stage == 5:
            stage_alleles(exp, main_body_container, upper_container, message_container)

        #=============================================== STAGE 7: Reports ==============================================
        if pipeline_stage == 6:
            stage_reports(exp, main_body_container, upper_container, message_container)
            
        #================================================ PROCESS UPLOADS ==============================================
        with main_body_container:
            with st.spinner('Processing uploaded files...'):
                parse.process_upload_queue(exp)

        #=============================================== UPPER INFO SECTION ============================================

        # build the upper info viewers after main content so that they reflect any changes
        upper_panels = [v for v in (st.session_state.get('info_panel1', 'None'),
                st.session_state.get('info_panel2', 'None')) if v != 'None']
        if any(upper_panels) and st.session_state.get('show_upper_info_viewer'):
            with upper_container:
                iv.show_info_viewer(upper_panels, st.session_state.get('upper_panel_height',250), 'upper_view_panels')

        #=============================================== LOWER INFO SECTION ============================================

        # info panel displays are updated at the bottom of the script, so that they reflect any changes

        lower_panels = [v for v in (st.session_state.get('info_panel3','None'),
                st.session_state.get('info_panel4','None')) if v != "None"]
        if any(lower_panels):
            with lower_container:
                iv.show_info_viewer(lower_panels, st.session_state.get('lower_panel_height',350), 'lower_view_panels')

        #=============================================== UPDATE MESSAGES ==============================================
        with message_container:
            #dc.display_temporary_messages()
            #dc.display_persistent_messages('main1')
            pass

        
        ### End of main display ###

        #================================================ AUTO SAVE ===================================================
        if 'pipeline' not in st.session_state:
            st.session_state['pipeline_stage'] = pipeline_stage
        elif st.session_state['pipeline_stage'] != pipeline_stage:
            st.session_state['pipeline_stage'] = pipeline_stage


        #================================================ PROGRESS REPORT =============================================
        # must come last as it launches an infinite loop until matching is complete
        if 'matching_in_progress' in st.session_state and st.session_state['matching_in_progress'] and pipeline_stage == 5:
            rundir, launch_msg, launch_prog, completion_msg, match_prog = st.session_state['matching_in_progress']
            report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog)
            

if __name__ == '__main__':
    main()


