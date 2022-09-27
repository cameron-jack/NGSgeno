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
from pathlib import PurePath
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components

from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode
from bin.experiment import Experiment, EXP_FN, load_experiment
from bin.util import output_error
from bin.util import CAP_VOLS, DEAD_VOLS
import bin.file_io as file_io
import bin.db_io as db_io
from bin.makehtml import generate_heatmap_html


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
            output_error(exc, msg='Could not create new folder: ' + newpath)
            return None, 'Failed to create new folder: ' + newpath
    else:
        if os.path.exists(os.path.join(newpath, EXP_FN)):
            return None, 'Experiment already exists with this name: ' + newpath
    print('Generating experiment: ', newpath)
    exp = Experiment(name=newpath.lstrip('run_'))
    return exp, ''


def folder_sb():
    """
    Home page
    Sidebar
    Selectbox for loading or creating new experiment
    """
    #New experiment

    ftab1, ftab2 = st.sidebar.tabs(["New Folder", "Existing Folder"])

    add_run_folder = ftab1.text_input('Create new run folder:')
    create_run_folder_button = ftab1.button(label='Create', key='create_run_folder_button')

    try:
        existing_run_folders = get_run_folders()
    except Exception as exc:
        output_error(exc, msg='Cannot locate NGSgeno folder')
        return            

    run_folder = ftab2.selectbox("Select a run folder to open", existing_run_folders)

    if add_run_folder and create_run_folder_button:
        add_run_folder_str = 'run_' + add_run_folder
        exp, msg = create_run_folder(add_run_folder_str)
        if exp:
            exp.save()                                                                      
            st.session_state['experiment'] = exp
            st.experimental_rerun()
        else:
            if 'already exists' in msg:
                st.markdown('<p style="color:#FF0000">Folder name already exists</p>', unsafe_allow_html=True)
            else:
                st.markdown('<p style="color:#FF0000">Fatal path error: ' + msg + '</p>', unsafe_allow_html=True)

    if run_folder:
        if st.session_state['experiment'] == None or st.session_state['experiment'].name != run_folder:
            
            ch_run_path = 'run_' + run_folder
            if os.path.exists(ch_run_path):
                exp = load_experiment(ch_run_path)
                
                #print(dir(exp))
                if not exp:
                    st.markdown('<p style="color:#FF0000">Could not load experiment from: ' + ch_run_path + '</p>', unsafe_allow_html=True)
                elif ch_run_path.endswith(exp.name):
                    # success!
                    st.session_state['experiment'] = exp
                    st.write(exp)
                    st.experimental_rerun()
                else:
                    st.markdown('<p style="color:#FF0000">Invalid experiment file in: ' + ch_run_path + '</p>', unsafe_allow_html=True)
        ch_exp_folder = None
