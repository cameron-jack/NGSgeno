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
from pathlib import PurePath
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO
from shutil import copy2

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
from bin.ngsmatch import match_alleles

import load_data as ld
import display_components as dc
from pipeline_stage import pipe_stages

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

def st_radio_horizontal(*args, **kwargs):
    """Trick to have horizontal st radio to simulate tabs"""
    col, _ = st.columns([3,1])
    with col:
        st.write('<style> div[data-testid=column] > div > div > div > div.stRadio > div{flex-direction: row}</style>', unsafe_allow_html=True)
        return st.radio(*args, **kwargs)

def clear_widget(*keys):
    for k in keys:
        if k in st.session_state:
            st.session_state[k] = ''

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

def pipeline_sb():
    """
    Stages of the pipeline as a drop down box in the sidebar
    """
    pipeline_stage = st.sidebar.selectbox("Pipeline Stage", ("Load Data", "1. Nimbus Plates", "2. Echo Primers (PCR 1)", "3. Echo Indexing (PCR 2)", 
                                                             "4. Miseq", "5. Allele Calling", "6. Genotyping", "7. Review"))
    if pipeline_stage == "Load Data":
        st.session_state['pipe_stage'] = 1

    if pipeline_stage == "1. Nimbus Plates":
        st.session_state['pipe_stage'] = 2

    if pipeline_stage == "2. Echo Primers (PCR 1)":
        st.session_state['pipe_stage'] = 3

    if pipeline_stage == "3. Echo Indexing (PCR 2)":
        st.session_state['pipe_stage'] = 4

    if pipeline_stage == "4. Miseq":
        st.session_state['pipe_stage'] = 5

    if pipeline_stage == "5. Allele Calling":
        st.session_state['pipe_stage'] = 6

    if pipeline_stage == "6. Genotyping":
        st.session_state['pipe_stage'] = 7
    
    if pipeline_stage == "7. Review":
        st.session_state['pipe_stage'] = 8

def folder_sb():
    """
    Selectbox for loading existing experiment or input text for
    creating new experiment in the sidebar
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

    error_msg_area = st.sidebar.empty()
    error_msg = None

    if add_run_folder and create_run_folder_button:
        add_run_folder_str = 'run_' + add_run_folder
        exp, msg = create_run_folder(add_run_folder_str)
        if exp:
            exp.save()                                                                      
            st.session_state['experiment'] = exp
            add_run_folder = None
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
                
                if not exp:
                    error_msg = "Could not load experiment from: "+ ch_run_path
                elif ch_run_path.endswith(exp.name):
                    # success!
                    st.session_state['experiment'] = exp
                    st.experimental_rerun()
                else:
                    error_msg = "Invalid experiment file in: " + ch_run_path
    
        ch_exp_folder = None

    if error_msg:
        error_msg_area.markdown(f'<p style="color:#FF0000; text-align:center">{error_msg}</p>', unsafe_allow_html=True)

def main():
    st.set_page_config(page_icon="ngsg_icon.png", page_title="NGS Genotyping Pipeline", initial_sidebar_state="expanded", layout="wide")
    # Remove whitespace from the top of the page and sidebar
    st.markdown(
        """
        <style>
               .css-18e3th9 {
                    padding-top: 0rem;
                    padding-bottom: 10rem;
                    padding-left: 5rem;
                    padding-right: 5rem;
                }
                header {visibility: hidden;}
                footer {visibility: hidden;}
        </style>
        """, unsafe_allow_html=True)
    padding = 0
    st.markdown(f""" <style>
    .reportview-container .main .block-container{{
        padding-top: {padding}rem;
        padding-right: {padding}rem;
        padding-left: {padding}rem;
        padding-bottom: {padding}rem;
    }} </style> """, unsafe_allow_html=True)
    #This is a workaround from https://www.codegrepper.com/code-examples/python/streamlit+sidebar+width
    st.markdown(
        """
        <style>
        div.block-container {padding-top:0rem;}
        [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
            width: 400px;
            top: 0rem;
            padding-top: 0rem;
            height: 100vh;
        }
        [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
            width: 400px;
            top: 0rem;
            padding-top: 0rem;
            margin-left: -400px;
        }
         
        </style>
        """,
        unsafe_allow_html=True)
    
    sidebar_title = st.sidebar.title('')
    sidebar_title.write("NGS Genotyping Pipeline")
    st.sidebar.image('ngsg_explorer.png')
    main_container = st.container()

    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')
     
    if 'experiment' in st.session_state:
        if 'navigation' not in st.session_state:
            st.session_state['navigation'] = 'load'

        if st.session_state['experiment']:
            st.session_state['navigation'] = 'pipe'

            exp = st.session_state['experiment']
            st.sidebar.header('Current experiment: '+ exp.name)

            pipeline_sb()
            title='Experiment '+ st.session_state['experiment'].name
            main_container.markdown(f'<h4 style="color:#BED2D6">{title}</h2>', unsafe_allow_html=True)

            if 'pipe_stage' not in st.session_state:
                st.session_state['pipe_stage'] = 1
            
            if 'lock' not in st.session_state:
               st.session_state['lock'] = False
            if 'unlock' not in st.session_state:
               st.session_state['unlock'] = False
            if exp.locked:
               #unlock_button = st.sidebar.button('Unlock experiment', help='Allow modification of experiment again')
               #if unlock_button and not st.session_state['unlock']:
               #    st.session_state['unlock'] = True
               pass
            else:
               #lock_button = st.sidebar.button('Lock experiment', help='Freezes the current experiment, preventing further modification')
               #if lock_button and not st.session_state['lock']:
               #    st.session_state['lock'] = True
               pass
            if st.session_state['unlock']:
               exp.unlock()
               st.session_state['unlock'] = False
            if st.session_state['lock']:
               exp.lock()
               st.session_state['lock'] = False

            #if 'nuked' not in st.session_state:
            #    st.session_state['nuked'] = False
            #nuke_button = st.sidebar.button('Delete experiment', help='Hides the current experiment from the user interface. Manual retrieval is required')
            #if nuke_button and not st.session_state['nuked']:
            #    st.session_state['nuked'] = True
            #if st.session_state['nuked']:
            #    st.sidebar.write('Currently not enabled, but in future will hide experiment')
            #    st.session
            else: # pipeline stages

            #    pipeline_sb()

            ##Pipeline Stages

                if 'pipe_stage' not in st.session_state or st.session_state['pipe_stage'] == 1:
                    
                    st.sidebar.write('')
                    st.sidebar.write('')
                    st.sidebar.markdown('<p style="color:#8AB1BD">Or change experiment</p>', unsafe_allow_html=True)
                    #gotta change logic so that this load still
                    folder_sb()
                    #dc.data_table(key='Test')
                    #something in data table
                    pipe_stages(1)
            
                #Nimbus
                if st.session_state['pipe_stage'] == 2:
                    pipe_stages(2)
                    
                #Echo primers   
                if st.session_state['pipe_stage'] == 3:
                    pipe_stages(3)

                if st.session_state['pipe_stage'] == 4:  # echo barcodes
                    pipe_stages(4)

                if st.session_state['pipe_stage'] == 5:  # Miseq
                    pipe_stages(5)

                if st.session_state['pipe_stage'] == 6:  # Genotyping
                    pipe_stages(6)

                if st.session_state['pipe_stage'] == 7:  # Review
                    pipe_stages(7)


    if not st.session_state['experiment'] or (st.session_state['experiment'] and st.session_state['navigation'] == 'load'):
        #st.session_state['experiment'] = None
        main_container.write('')
        main_container.write('')
        main_container.markdown(\
                    '<h2 style="text-align:center;color:#154c79">'\
                            'Please load existing experiment or create a new one</h2>',\
                                        unsafe_allow_html=True)

        folder_sb()


    #hvar = """  <script>
    #                var elements = window.parent.document.querySelectorAll('.streamlit-expanderHeader');
    #                elements.forEach(function (element) {
    #                    element.style.fontSize = 'large';
    #                    });         
    #            </script>"""
    #components.html(hvar, height=0, width=0)

    

    #st.markdown(
    #    """
    #<script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    #<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    #<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    #""",
    #    unsafe_allow_html=True,
    #)


if __name__ == '__main__':
    main()
