#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__version__ = "2.03.004"

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

try:
    import bin.match as match
except ModuleNotFoundError:
    import match

try:
    import bin.ngsmatch as ngsmatch
except ModuleNotFoundError:
    import ngsmatch

#import bin.file_io as file_io
#import bin.db_io as db_io

import extra_streamlit_components as stx
import display_components as dc
import load_data as ld
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


def clean_up_match_files(rundir):
    """
    Remove any existing match progress files to ensure a clean start
    """
    lock_path = Path(rundir+'/ngsgeno_lock')
    if os.path.exists(lock_path):
        try:
            os.remove(lock_path)
        except Exception as exc:
            print(f'Could not remove lock file: {exc}', file=sys.stderr)
            #m(f'Could not remove lock file: {exc}', level='error', dest='noGUI')
    progress_files = list(Path(rundir).glob('match_progress_*'))
    for pf in progress_files:
        try:
            os.remove(pf)
        except Exception as exc:
            #m(f'Could not remove progress file: {exc}', level='error', dest='noGUI')
            print(f'Could not remove progress file: {exc}', file=sys.stderr)


def report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog):
    """
    Allows the interface to keep running while a background process (ngsmatch.py) progress is tracked
    """
    while True:
        launch_progress = 0
        match_progress = 0
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
            return


def run_generate(exp, target_func, *args, **kwargs):
    """
    Run target_func() and ask user to respond to any file/pipeline clashes in the given context, if provided
    Saves the experiment after use
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
    exp.save()
    return success


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

    dc.add_css()

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
        def upper_info_viewer_code(tab_col3, tab_col2, widget_key, default_view1='None', default_view2='None', checked=True):
            with tab_col3:
                add_vertical_space(2)
                dc.show_upper_info_viewer_checkbox(widget_key, value=checked)
            if st.session_state['show_upper_info_viewer']:
                with tab_col2:
                    ignore = dc.info_selection(widget_key+"top_viewer", 'info_panel1', 'info_panel2',
                            'upper_panel_height', default_view1=default_view1, default_view2=default_view2,
                            default_height=st.session_state.get('upper_panel_height',250))
            caller_id = 'display_feedback'
            # display any messages for this widget
            if caller_id in mq:
                for msg, lvl in mq[caller_id]:
                    m(msg, level=lvl, no_log=True)
                sleep(0.3)
                mq[caller_id] = set()

        # attempt to parse any files that are set for upload
        parse.process_upload_queue(exp)
        if '_upload_pending' not in exp.uploaded_files:
            exp.uploaded_files['_upload_pending'] = {}

        with save_container:
            save_message(exp, key = 'save1')

        # info panel displays are updated at the bottom of the script, so that they reflect any changes
        with lower_container:
            success = dc.info_selection("bottom_viewer", 'info_panel3', 'info_panel4',
                    'lower_panel_height', default_view1='Log', default_view2='None',
                    default_height=st.session_state.get('lower_panel_height',350))

        #============================================== STAGE 1: Load data =============================================
        if pipeline_stage == 0:

            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
            with tab_col1:
                #load_data_tab = dc.create_tabs([("Load Samples", "Rodentity or custom samples for the complete pipeline"),
                #        ("Load Consumables", "Additional plates and information"),
                #        ("Load Amplicons","Additional amplicons to be sequenced")])
                load_data_tab = dc.create_tabs([("Load Samples", ""),("Load Consumables", ""),
                        ("Load Amplicons", "")])
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
                        st.markdown('**Upload Rodentity plate files, custom manifests, or 384-plate reruns, '+\
                                'then move to the next tab to upload pipeline consumables and other important files**')
                        with st.expander('Upload Rodentity Ear Punch Plates', expanded=True):
                            ld.load_rodentity_data('rodentity_load1')
                            ld.assign_rodentity_dna_plate('rodentity_load2')
                        hline()
                        ld.load_custom_manifests('custom_load1')
                        hline()
                        ld.load_manifest_384('manifest_384_load1')
                        hline()

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
                        st.markdown('**Upload consumables, PCR files, and custom volumes, '+\
                                'then move to the next tab to upload amplicon files**')
                        ld.load_extra_consumables('consumables_load2')
                        ld.load_pcr1_files('pcr1_load2')
                        ld.load_pcr2_files('pcr2_load2')
                        add_vertical_space(1)
                        st.subheader('Custom Volumes')
                        ld.custom_volumes('custom1')
                        hline()
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_consumables', default_view1='Consumables',
                        default_view2='Samples')

            #------------------------------------ Load ~ TAB 3: Load amplicons ---------------------------------------
            if load_data_tab == 3:
                st.session_state['load_tab'] = 3
                if unlocked(exp):
                    init_state('upload stage', None)
                    with main_body_container:
                        st.subheader('Upload Amplicon Files')
                        st.markdown('**Upload amplicon plates and reference files**')
                        ld.load_amplicons('amp_load3')
                        add_vertical_space(1)
                        ld.load_amplicon_references('ref_load3')
                        hline()
                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_consumables', default_view1='Files',
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
                        _, header_col, _ = st.columns([2,2,1])
                        with header_col:
                            header_col.subheader('Generate Echo Files')
                            dc.set_nimbus_title(exp, nfs, efs)
                        add_vertical_space(2)
                        if not ld.check_assay_file():
                            m('Assay list file (assay to primer mapping) is required to proceed. '+\
                                    'Please upload this to continue', level='error')
                            success = ld.load_assaylist('nimbus1_assaylist')
                            if success:
                                st.rerun()
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
                                                level='error')
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
                        ld.load_echo_inputs('1')

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
                            st.warning('Do you need to upload Echo input files (output from Nimbus)?')

                        checkbox_cont = st.container()
                        hline()
                        add_vertical_space(1)

                        display_cont = st.container()
                        hline()

                        add_vertical_space(1)
                        st.subheader('Add Barcodes', help='Add barcodes for plates')
                        pcr_col, taqwater_col = st.columns(2)

                        st.subheader('Upload Files')
                        ld.load_pcr1_files('pcr1_primer1')
                        add_vertical_space(1)

                        st.subheader('Custom Volumes')
                        ld.custom_volumes('cust_vol')
                        selected_pids = {}
                        with checkbox_cont:
                            checkbox_keys = dc.display_plate_checklist('pmr1_checklist',
                                ['dna','pcr','taqwater1','primer'])

                            selected_pids = dc.collect_plate_checklist(checkbox_keys)
                            if not selected_pids['pcr']:
                                m('No PCR plates selected/available', level='display', dest=('css',), color='red',size='p')

                        with pcr_col:
                            ld.add_pcr_barcodes('pcr_bc_tab1')
                        with taqwater_col:
                            ld.add_taqwater_barcodes('tw_bc_tab1', pcr_stage=pcr_stage)

                        with display_cont:
                            dc.display_pcr1_components(selected_pids)

                            add_vertical_space(1)
                            if selected_pids['dna']:
                                primer_max_vol = util.CAP_VOLS[util.PLATE_TYPES['Echo384']]
                                primer_dead_vol = util.DEAD_VOLS[util.PLATE_TYPES['Echo384']]
                                if not selected_pids['primer']:
                                    st.warning(f'Primers will appear highlighted red if no primer plate files have been loaded')
                                st.write(f'Available volumes equal the measured volume - dead volume ({primer_dead_vol/1000}ul). Max primer volume is {primer_max_vol/1000}ul.')
                                dc.display_primers('pcr_tab1', dna_pids=selected_pids['dna'],
                                        primer_pids=selected_pids['primer'], save_buttons=True)

                        add_vertical_space(1)

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
                        if not selected_pids['pcr'] and not selected_pids['amplicon']:
                            m('No PCR or amplicon plates selected/added yet', level='display', dest=('css',), color='red',size='p')
                        hline()
                        #dc.display_pcr_common_components(selected_pids)
                        dc.display_pcr1_components(selected_pids)
                        hline()

                        if selected_pids['dna']:
                            if exp.check_ready_pcr1(selected_pids, caller_id=caller_id):
                                _,button_col,cb_col,_ = st.columns([5, 2, 2, 4])
                                echo_picklist_go = button_col.button('Generate Echo Picklists',
                                            key='echo_pcr1_go_button',
                                            type='primary')
                                cb_force = cb_col.checkbox('Force picklist generation', key='force_pcr1_picklist')
                                if echo_picklist_go:
                                    success = run_generate(exp, exp.generate_echo_PCR1_picklists,
                                            selected_pids, force=cb_force, caller_id=caller_id)
                                    if not success:
                                        st.error('Picklist generation failed. Please see the log')
                                    else:
                                        exp.add_pcr_wells(exp, selected_pids['pcr'], selected_pids['dna'])

                        for msg,lvl in mq[caller_id]:
                            m(msg, lvl)
                        mq[caller_id] = set()

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

                        checkbox_cont = st.container()
                        hline()
                        add_vertical_space(1)

                        display_cont = st.container()
                        hline()
                        add_vertical_space(1)

                        st.info('Indexing (PCR 2) requires index index plates, taq/water plates, '+\
                                'and either Echo plates (prepared by the Nimbus) or amplicon plates')

                        st.subheader('Add Barcodes', help='Add barcodes for plates')
                        pcr_col, taqwater_col = st.columns(2)

                        st.subheader('Upload Files')
                        ld.load_pcr2_files('pcr2_index1')
                        add_vertical_space(1)

                        st.subheader('Custom Volumes')
                        ld.custom_volumes(exp)

                        with checkbox_cont:
                            checkbox_keys = dc.display_plate_checklist('idx_checklist1',
                                    ['pcr','taqwater2','amplicon','index'])

                            selected_pids = dc.collect_plate_checklist(checkbox_keys)
                            if not selected_pids['pcr'] and not selected_pids['amplicon']:
                                m('No PCR or amplicon plates selected', level='display',
                                        dest=('css',), color='red',size='p')

                        with display_cont:
                            dc.display_pcr2_components(selected_pids)

                        with pcr_col:
                            ld.add_pcr_barcodes('pcr_bc_tab2')
                        with taqwater_col:
                            ld.add_taqwater_barcodes('tw_bc_tab2', pcr_stage=pcr_stage)
                        add_vertical_space(1)

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
                        if not selected_pids['pcr'] and not selected_pids['amplicon']:
                            m('No PCR or amplicon plates selected',
                                    level='display', dest=('css',), color='red',size='p')
                        hline()
                        print(f'{selected_pids=}', flush=True)
                        #dc.display_pcr_common_components(selected_pids)
                        dc.display_pcr2_components(selected_pids)
                        print(f'{selected_pids=}', flush=True)
                        hline()
                        add_vertical_space(1)

                        if selected_pids['pcr']:
                            success = exp.check_ready_pcr2(selected_pids, caller_id=caller_id)
                            if success:
                                do_generate = True
                        elif selected_pids['amplicon']:
                            success = exp.check_ready_pcr2(selected_pids, caller_id=caller_id, amplicon_only=True)
                            if success:
                                do_generate = True
                        else:
                            m('Either PCR plates or Amplicon plates must exist for indexing to begin',
                                    level='display')

                        for msg,lvl in mq[caller_id]:
                            m(msg, lvl)
                        mq[caller_id] = set()

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
                                m('Picklist generation failed. Please see the log', level='display')

                        if generate.pcr2_picklists_exist(exp):
                            dc.show_echo2_outputs()

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_index2', default_view1='Files',
                        default_view2='Files')

        #=============================================== STAGE 5: Miseq ================================================
        if pipeline_stage == 4:
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
            info_holder = st.container()

            with tab_col1:
                with tab_col1:
                    #miseq_tab = dc.create_tabs([("Download", "Miseq Samplesheet"), ("Upload", "Miseq Sequence Files")])
                    miseq_tab = dc.create_tabs([("Download", "Miseq Samplesheet")])
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

                    #ld.load_rodentity_references('reference_miseq1')
                    if exp.get_miseq_samplesheets():
                        dc.get_miseq_download_btn(exp)
                        add_vertical_space(4)

                    else:
                        st.warning(f'No MiSeq Samplesheet available for download')

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_miseq1', default_view1='Files',
                        default_view2='Plates')

            #------------------------ Miseq ~ TAB 2: Upload Reference File and Miseq Sequences -------------------------
            if miseq_tab == 2:
                st.session_state['miseq_tab'] = 2
                with main_body_container:
                    st.subheader('Upload Custom Reference Files')
                    ld.load_rodentity_references('reference_miseq2')
                    add_vertical_space(2)

                    caller_id = 'seq_upload_ready'
                    success = exp.check_sequence_upload_ready(caller_id)
                    if caller_id in mq:
                        for msg, lvl in mq[caller_id]:
                            m(msg, level=lvl, no_log=True)
                        sleep(0.3)
                    mq[caller_id] = set()
                    if not success:
                        st.warning('Resources are required for allele calling and must be present before FASTQs can be uploaded')
                    else:
                        ld.load_miseq_fastqs('miseq_tab1')


                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_miseq2', default_view1='Files',
                        default_view2='Plates')

        #=========================================== STAGE 6: Allele Calling ===========================================
        if pipeline_stage == 5:
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])

            with tab_col1:
                allele_tab = dc.create_tabs([("Allele Calling", "For Rodentity Mice"),("Amplicon Calling", "For custom amplicons")])
            if not allele_tab:
                init_state('allele_tab', 1)
                allele_tab = st.session_state['allele_tab']

            # Only offer upload in the Miseq pipeline section
            #-------------------------------- Allele ~ TAB 1: Allele calling --------------------------------
            if allele_tab == 1:
                st.session_state['allele_tab'] = 1
                with main_body_container:
                    rundir = exp.get_exp_dn()
                    caller_id = 'pre-execute-gt-analysis'
                    success = exp.check_sequence_upload_ready(caller_id)
                    for msg,lvl in mq[caller_id]:
                        m(msg, level=lvl, no_log=True)
                    mq[caller_id] = set()
                    ld.display_fastqs('disp_fastq1')
                    if not success:
                        st.warning('Resources are required for allele calling:')
                        st.subheader('Upload reference sequences')
                        ld.load_rodentity_references('reference_allele1')
                        #ld.load_miseq_fastqs('miseq_tab1a')
                    else:
                        #ld.load_miseq_fastqs('miseq_tab1b')
                        # check whether output files already exist and are opened elsewhere
                        targets_fn = exp.get_exp_fn('target.fa')
                        primers_fn = exp.get_exp_fn('primers.csv')
                        results_fn = exp.get_exp_fn('results.csv')
                        matchlog_fn = exp.get_exp_fn('match.log')
                        for fn in [targets_fn, primers_fn, results_fn, matchlog_fn]:
                            if Path(fn).exists():
                                try:
                                    os.rename(fn, fn.replace('.','_tmp_swap.'))
                                    os.rename(fn.replace('.','_tmp_swap.'), fn)
                                except PermissionError:
                                    st.error(f'{fn} appears to be in use. Please close this file before continuing')

                        with st.form('allele_calling_form', clear_on_submit=True):
                            all_fns = [fp for fp in exp.uploaded_files.keys() if fp not in {'_upload_queue','_upload_pending'}]
                            ref_fns = [fp for fp in all_fns if exp.uploaded_files[fp]['purpose'] \
                                    in ['rodentity_reference','custom_reference']]
                            ref = st.selectbox('Select references to match against', options=ref_fns)
                            #num_unique_seq = st.number_input("Number of unique sequences per work unit", value=1)
                            cpus_avail = max(1, os.cpu_count()-2)
                            num_cpus = st.number_input(\
                                    label=f"Number of processes to run simultaneously, default: {cpus_avail}",\
                                            value=cpus_avail, min_value=1)
                            margin = st.number_input(label="Require lengths of read sequences and target "+\
                                    "references to be proportionally similar by this amount. Value must be between 0.0 and 1.0 "+\
                                    "default 0.9", format='%f',min_value=0.0, step=0.05,value=0.9)
                            identity = st.number_input(label="Proportion of identity required for inexact match "+\
                                    ", default 0.9. Must be between 0.0 and 1.0",
                                    format='%f',min_value=0.0, max_value=1.0, value=0.9, step=0.05)
                            mincov = st.number_input(label="Do not match unique sequences with less than this "+\
                                    "many reads coverage, default 5", format='%i',min_value=0, step=1,value=5)
                            minprop = st.number_input(label="Do not match unique sequences with less than this "+\
                                    "proportion of the reads seen for the most observed (expected) allele, default 0.1. Must be between 0.0 and 1.0",
                                    format='%f',min_value=0.0, max_value=1.0, value=0.1, step=0.05)
                            inexact_mode = st.checkbox("Enable inexact matching")
                            exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                            debug_mode = st.checkbox('Turn on debugging for allele calling')
                            do_matching = st.form_submit_button("Run allele calling")

                        if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                            match_fn = exp.get_exp_fn('match_progress_100_100')
                            if Path(match_fn).exists():
                                clean_up_match_files()
                                st.rerun()

                        elif do_matching:
                            success = generate.generate_targets(exp, ref, caller_id=caller_id)
                            if not success:
                                m('failed to save reference sequences to target file', level='critical')
                                sleep(0.5)
                            success = generate.generate_primer_assayfams(exp, caller_id=caller_id)
                            if not success:
                                m('failed to save primers and assay families to file', level='critical')
                                sleep(0.5)
                            else:
                                matching_prog = os.path.join('bin','ngsmatch.py')
                                cmd_str = f'{sys.executable} {matching_prog} ' +\
                                         f'--ncpus {num_cpus} --rundir {rundir} '+\
                                         f'--margin {margin} --identity {identity} --mincov {mincov} '+\
                                         f'--minprop {minprop}'
                                if inexact_mode:
                                    cmd_str += ' --inexact'
                                if exhaustive_mode:
                                    cmd_str += ' --exhaustive'
                                if debug_mode:
                                    cmd_str += ' --debug'
                                m(f'{cmd_str}', level='info')
                                st.write(f'Calling {cmd_str}')

                                launch_msg = st.empty()
                                launch_prog = st.progress(0)
                                completion_msg = st.empty()
                                match_prog = st.progress(0)
                                st.session_state['matching_in_progress'] = (rundir, launch_msg, launch_prog,
                                        completion_msg, match_prog)
                                if sys.platform == "win32":
                                    subprocess.Popen(cmd_str.split(' '),
                                            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS,
                                            stdin=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL,
                                            close_fds=True,
                                    )
                                else:
                                    subprocess.Popen(cmd_str.split(' '),
                                            start_new_session=True,
                                            stdin=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL,
                                            close_fds=True,
                                    )
                                # cp = subprocess.Popen(cmd_str.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                # st.session_state['matching_process'] = cp

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_allele1', default_view1='Files',
                        default_view2='Plates')

            #-------------------------------- Allele ~ TAB 2: Amplicon calling --------------------------------
            elif allele_tab == 2:
                st.session_state['allele_tab'] = 2
                with main_body_container:
                    rundir = exp.get_exp_dn()
                    caller_id = 'pre-execute-amp-analysis'
                    st.subheader('Amplicon Calling')
                    st.info('This section is for custom amplicon calling. Amplicon sequences are always inexact matched against their expected primer to find variants. Please choose the amplicon '+\
                            'reference you wish to use or upload a new one')
                    ld.load_amplicon_references('amplicon_ref1')

                    targets_fn = exp.get_exp_fn('amplicon_targets.fa')
                    results_fn = exp.get_exp_fn('amplicon_results.csv')
                    matchlog_fn = exp.get_exp_fn('amplicon_match.log')
                    # check whether output files already exist and are opened elsewhere
                    for fn in [targets_fn, results_fn, matchlog_fn]:
                        if Path(fn).exists():
                            try:
                                os.rename(fn, fn.replace('.','_tmp_swap.'))
                                os.rename(fn.replace('.','_tmp_swap.'), fn)
                            except PermissionError:
                                st.error(f'{fn} appears to be in use. Please close this file before continuing')
                    # choose amplicon plate(s) and reference file(s) to run
                    plate_col, ref_col = st.columns(2)
                    with plate_col:
                        st.subheader('Select amplicon plates for allele calling')
                        checkbox_keys = dc.display_plate_checklist('amplicon_checklist', ['amplicon'], default_value=False)
                        selected_pids = dc.collect_plate_checklist(checkbox_keys)
                        if not selected_pids['amplicon']:
                            m('No amplicon plates selected', level='display', dest=('css',), color='red',size='p')
                        else:
                            st.write(f'Selected amplicon plates: {", ".join([util.unguard_pbc(gpid, silent=True) for gpid in selected_pids["amplicon"]])}')
                    with ref_col:
                        st.subheader('Select amplicon references to match against')
                        checkbox_keys = dc.display_file_checklist('amplicon_reference_checklist', ['amplicon_reference'], default_value=False)
                        selected_refs = dc.collect_file_checklist(checkbox_keys)
                        if not selected_refs['amplicon_reference']:
                            m('No amplicon references selected', level='display', dest=('css',), color='red',size='p')
                        else:
                            st.write(f'Selected amplicon references: {", ".join(dc.fns_from_checklist(selected_refs))}')

                    with st.form('amplicon_calling_form', clear_on_submit=True):
                        cpus_avail = max(1, os.cpu_count()-2)
                        num_cpus = st.number_input(\
                                label=f"Number of processes to run simultaneously, default: {cpus_avail}",\
                                        value=cpus_avail, min_value=1)
                        mincov = st.number_input(label="Do not match unique sequences with less than this "+\
                                "many reads coverage, default 5", format='%i',min_value=0, step=1,value=5)
                        minprop = st.number_input(label="Do not match unique sequences with less than this "+\
                                "proportion of the reads seen for the most observed (expected) allele, default 0.1. Must be between 0.0 and 1.0",
                                format='%f',min_value=0.0, max_value=1.0, value=0.1, step=0.05)
                        margin = st.number_input(label="Require lengths of read sequences and target "+\
                                "references to be proportionally similar by this amount. Value must be between 0.0 and 1.0 "+\
                                "default 0.9", format='%f',min_value=0.0, step=0.05,value=0.9)
                        identity = st.number_input(label="Proportion of identity required for inexact match "+\
                                ", default 0.9. Must be between 0.0 and 1.0",
                                format='%f',min_value=0.0, max_value=1.0, value=0.9, step=0.05)
                        exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                        debug_mode = st.checkbox('Turn on debugging for allele calling')
                        do_matching = st.form_submit_button("Run amplicon calling")

                    if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                        st.info('Analysis in progress')
                    elif do_matching:
                        success = generate.generate_targets(exp, selected_refs['amplicon_reference'])
                        if not success:
                            m('failed to save reference sequences to target file', level='critical')
                            sleep(0.5)
                        else:
                            matching_prog = Path('bin/ngsmatch.py')
                            cmd_str = f'{sys.executable} {matching_prog} --ncpus {num_cpus} --rundir {rundir} '+\
                                    f'--mincov {mincov} --minprop {minprop}'
                            if exhaustive_mode:
                                cmd_str += ' --exhaustive'
                            if debug_mode:
                                cmd_str += ' --debug'
                            if selected_pids['amplicon']:
                                for pid in selected_pids['amplicon']:
                                    cmd_str += f' --amplicon {",".join(selected_pids["amplicon"])}'
                                cmd_str += f' --targets amplicon_targets.fa'
                                m(f'{cmd_str}', level='info')
                                st.write(f'Calling {cmd_str}')
                                launch_msg = st.empty()
                                launch_prog = st.progress(0)
                                completion_msg = st.empty()
                                match_prog = st.progress(0)
                                st.session_state['matching_in_progress'] = (rundir, launch_msg, launch_prog,
                                        completion_msg, match_prog)
                                if sys.platform == "win32":
                                    subprocess.Popen(cmd_str.split(' '),
                                            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS,
                                            stdin=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL,
                                            close_fds=True,
                                    )
                                else:
                                    subprocess.Popen(cmd_str.split(' '),
                                            start_new_session=True,
                                            stdin=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL,
                                            close_fds=True,
                                    )
                                #cp = subprocess.Popen(cmd_str.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                #st.session_state['matching_process'] = cp

                # ** Info viewer **
                upper_info_viewer_code(tab_col3, tab_col2, 'upper_allele2', default_view1='Files',
                        default_view2='Plates')

        #=============================================== STAGE 7: Reports ==============================================
        if pipeline_stage == 6:
            results_fp = exp.get_exp_fn('results.csv')
            tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
            with tab_col1:
                allele_tab = dc.create_tabs([("Sequencing results", ""),("Amplicon results", "")])
            if not allele_tab:
                init_state('results_tab', 1)
                allele_tab = st.session_state['results_tab']

            with main_body_container:
                if allele_tab == 1:
                    if not Path(results_fp).exists():
                        m('**No allele calling results available**', level='display', dest=('mkdn'))
                    else:
                        if not Path(exp.get_exp_fn('amplicon_targets.fa')).exists():
                            m('**No amplicon targets available**', level='display', dest=('mkdn'))
                        else:
                            # Three displays for Rodentity, custom, and other results
                            key = 'rodentity_results_display'
                            dc.show_results_display('rodentity', key, caller_id=key)
                            key = 'custom_results_display'
                            dc.show_results_display('custom', key, caller_id=key)
                            key = 'other_results_display'
                            dc.show_results_display('other', key, caller_id=key)

                elif allele_tab == 2:
                    if not Path(exp.get_exp_fn('amplicon_results.csv')).exists():
                        m('**No amplicon calling results available**', level='display', dest=('mkdn'))
                    else:
                        if not Path(exp.get_exp_fn('amplicon_targets.fa')).exists():
                            m('**No amplicon targets available**', level='display', dest=('mkdn'))
                        else:
                            key = 'amplicon_results_display'
                            dc.show_results_display('amplicon', key, caller_id=key)
                            
            # ** Info viewer **
            upper_info_viewer_code(tab_col3, tab_col2, 'upper_report1', default_view1='Files',
                    default_view2='Plates')

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
            dc.display_persistent_messages('main1')

        parse.process_upload_queue(exp)
        ### End of main display ###

        #================================================ AUTO SAVE ===================================================
        if 'pipeline' not in st.session_state:
            st.session_state['pipeline_stage'] = pipeline_stage
        elif st.session_state['pipeline_stage'] != pipeline_stage:
            st.session_state['pipeline_stage'] = pipeline_stage


        #================================================ PROGRESS REPORT =============================================
        if 'matching_in_progress' in st.session_state and st.session_state['matching_in_progress']:
            rundir, launch_msg, launch_prog, completion_msg, match_prog = st.session_state['matching_in_progress']
            report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog)
            

if __name__ == '__main__':
    main()


