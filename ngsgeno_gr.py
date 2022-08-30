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
#from load_css import local_css

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



def aggrid_interactive_table(df: pd.DataFrame, grid_height: int=250):
    """Creates an st-aggrid interactive table based on a dataframe.

    Args:
        df (pd.DataFrame]): Source dataframe

    Returns:
        dict: The selected row
    """
    if 'view_box_size' in st.session_state:
        grid_height = st.session_state['view_box_size']
    options = GridOptionsBuilder.from_dataframe(
        df, enableRowGroup=True, enableValue=True, enablePivot=True
    )

    options.configure_side_bar()

    options.configure_selection("multiple", use_checkbox=False, rowMultiSelectWithClick=False, suppressRowDeselection=False)
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        height=grid_height,
        gridOptions=options.build(),
        # _available_themes = ['streamlit','light','dark', 'blue', 'fresh','material']
        theme="light",
        #update_mode=GridUpdateMode.MODEL_CHANGED,
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        allow_unsafe_jscode=True,
        rowSelection='multiple',
        selection_mode='multiple',
        rowMultiSelectWithClick=True,
    )

    return selection

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


def st_radio_horizontal(*args, **kwargs):
    """Trick to have horizontal st radio to simulate tabs"""
    col, _ = st.columns([3,1])
    with col:
        st.write('<style> div[data-testid=column] > div > div > div > div.stRadio > div{flex-direction: row}</style>', unsafe_allow_html=True)
        return st.radio(*args, **kwargs)


def delete_entries(exp, entries):
    """ removes data from an experiment, called by buttons """
    print(f"Debug: (delete_entries) {exp.name=} {entries=}")
    if exp is None:
        return
    if entries: # list
        for i,entry in enumerate(entries):
            for key in entry:
                #print(key)
                if 'PID' in key and entry[key]:
                    if not file_io.is_guarded_pbc(entry[key]):
                        entries[i][key] = file_io.guard_pbc(entry[key])
        print(f"delete_entries {entries=}")
        exp.remove_entries(entries)
        exp.save()
    print('delete_entries ran')


def display_reaction_stats(assay_usage, show_general=True, show_primers=True):
    """
    Required reaction wells, PCR plats, water volumes, taq volumes, number of wells, uses, volume needed and volume required per primer
    """
    if assay_usage:
        primer_vols, primer_taq_vol, primer_water_vol, barcode_taq_vol, barcode_water_vol =\
                st.session_state['experiment'].get_volumes_required(assay_usage=assay_usage)
        reactions = sum([v for v in assay_usage.values()])
        barcodes_remain, barcodes_avail, barcode_vol_okay = st.session_state['experiment'].get_barcode_remaining_available_volume(assay_usage=assay_usage)

        #for k,p in enumerate(sorted(primer_names)):
        #    if p != '':
        #        r = (k*2)%6  # horizontal offset
        #        uses_columns[r+0].write(p)
        #        uses_columns[r+1].write(assay_usage.get(p,0))
        

    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, barcode_taq_vol, barcode_water_vol =\
            0,{},0,0,0,0 
        barcodes_remain, barcodes_avail, barcode_vol_okay =\
                st.session_state['experiment'].get_barcode_remaining_available_volume()
    primer_avail_counts, primer_avail_vols = st.session_state['experiment'].get_primers_avail()



    if show_general:
        with st.expander('Required Componenents for PCR Reactions', expanded=False):
            taq_water_avail_PCR1 = st.session_state['experiment'].get_taq_water_avail()
            taq_water_avail_PCR2 = st.session_state['experiment'].get_taq_water_avail()

            #Ask Cam about this - filtered vs unfiltered
            #Sometimes doesn't show up..
            filter_columns = st.columns([1,5,1,5,1])
            if 'assay_filter' not in st.session_state:
                st.session_state['assay_filter'] = 1

                filter_assays = filter_columns[1].button('Filter assays', key='assay_filter_button', disabled=True)
                unfilter_assays = filter_columns[3].button('Unfilter assays', key='assay_unfilter_button')


                if unfilter_assays:
                    st.session_state['assay_filter'] = 0

                if filter_assays:
                    st.session_state['assay_filter'] = 1

           # if st.session_state['assay_filter'] == 1:


            req_cols = st.columns([6, 4, 6, 4])

            req_cols[0].markdown('<h6 style="color:#B1298D">Required reaction wells</h6>', unsafe_allow_html=True)
            req_cols[1].write(str(reactions))
            req_cols[2].markdown('<h6 style="color:#B1298D">Required PCR plates</h6>', unsafe_allow_html=True)
            req_cols[3].write(str(ceil(reactions/384)))
            st.write('')

            #PCR1 componenents - available amounts need taq water plates, still need to implement
            st.markdown('<h6 style="color:#B1298D">PCR 1</h6>', unsafe_allow_html=True)
            if taq_water_avail_PCR1[0] == 0 or taq_water_avail_PCR1[1] == 0:
                st.markdown('<p style="color:#B92A5D"> Upload taq and water plates to fulfill PCR 1 requirements </p>', unsafe_allow_html=True)
            PCR1_cols = st.columns([6, 4, 6, 4])
            PCR1_cols[0].markdown('**Required water volume**', unsafe_allow_html=True)
            PCR1_cols[1].write(str(primer_water_vol/1000)+' μl')
            PCR1_cols[2].markdown('**Available water volume**')
            PCR1_cols[3].write(str(taq_water_avail_PCR1[1]/1000)+' μl')
            PCR1_cols[0].markdown('**Required taq volume**')
            PCR1_cols[1].write(str(primer_taq_vol/1000)+' μl')
            PCR1_cols[2].markdown('**Available taq volume**')
            PCR1_cols[3].write(str(taq_water_avail_PCR1[0]/1000)+' μl')
            st.write('')
            #PCR 2 components - avail amounts plus barcodes?
            st.markdown('<h6 style="color:#B1298D">PCR 2</h6>', unsafe_allow_html=True)
            PCR2_cols = st.columns([6, 4, 6, 4])
            PCR2_cols[0].markdown('**Required water volume**')
            PCR2_cols[1].write(str(barcode_water_vol/1000)+' μl')
            PCR2_cols[2].markdown('**Available water volume**')
            PCR2_cols[3].write(str(taq_water_avail_PCR2[1]/1000)+' μl')
            PCR2_cols[0].markdown('**Required taq volume**')
            PCR2_cols[1].write(str(barcode_taq_vol/1000)+' μl')
            PCR2_cols[2].markdown('**Available taq volume**')
            PCR2_cols[3].write(str(taq_water_avail_PCR2[0]/1000)+' μl')
            
            #st.write('')
            #st.markdown('<h6 style="color:#B1298D">Barcode Pairs</h6>', unsafe_allow_html=True)
            bar_cols = st.columns([6, 4, 6, 4])
            bar_cols[0].markdown('**Barcode Pairs Remaining**')
            bar_cols[1].write(str(barcodes_remain))
            bar_cols[2].markdown('**Barcode Pairs Available**')
            bar_cols[3].write(str(barcodes_avail))
            bar_cols[0].markdown('**Sufficient Barcode Volumes**')
            if barcode_vol_okay:
                bar_cols[1].markdown('<p style="color:green">True</p>', unsafe_allow_html=True)
            else:
                bar_cols[1].markdown('<p style="color:red">False</p>', unsafe_allow_html=True)

    #Check all the calculations with Cam
    if show_primers:
        #Probably should write a function for this
        with st.expander('Primer volumes', expanded=False):
            
            ptab1, ptab2, ptab3, ptab4 = st.tabs(["Wells", "Uses", "Volume μL", "Available μL"])

            #Number of wells
            with ptab1:
                uses_columns = st.columns([3, 3, 3, 3, 3, 3])
                primer_names = set(primer_vols.keys())
                uses_columns[0].markdown('**Primer**')
                uses_columns[2].markdown('**Primer**')
                uses_columns[4].markdown('**Primer**')
                uses_columns[1].markdown('**Wells**')
                uses_columns[3].markdown('**Wells**')
                uses_columns[5].markdown('**Wells**')
                for k,p in enumerate(sorted(primer_names)):
                    if p != '':
                        r = (k*2)%6  # horizontal offset
                        uses_columns[r+0].write(p)
                        need_vol = primer_vols.get(p,0)
                        #this should probably be in experiment.py
                        per_use_vol = (CAP_VOLS['384PP_AQ_BP'] - DEAD_VOLS['384PP_AQ_BP'])
                        num_wells = ceil((need_vol)/per_use_vol)
                        uses_columns[r+1].write(num_wells)



            with ptab2:
                uses_columns = st.columns([3, 3, 3, 3, 3, 3])
                primer_names = set(primer_vols.keys())
                uses_columns[0].markdown('**Primer**')
                uses_columns[2].markdown('**Primer**')
                uses_columns[4].markdown('**Primer**')
                uses_columns[1].markdown('**Uses**')
                uses_columns[3].markdown('**Uses**')
                uses_columns[5].markdown('**Uses**')
                for k,p in enumerate(sorted(primer_names)):
                    if p != '':
                        r = (k*2)%6  # horizontal offset
                        uses_columns[r+0].write(p)
                        uses_columns[r+1].write(assay_usage.get(p,0))

            #Amount of volume for each primer
            with ptab3:
                #9 columns (8 in python terms)
                volume_columns = st.columns([3, 2, 2, 3, 2, 2, 3, 2, 2])
                primer_names = set(primer_vols.keys())
                volume_columns[0].markdown('**Primer**')
                volume_columns[3].markdown('**Primer**')
                volume_columns[6].markdown('**Primer**')
                volume_columns[1].markdown('**Volume**')
                volume_columns[4].markdown('**Volume**')
                volume_columns[7].markdown('**Volume**')
                volume_columns[2].markdown('**Available**')
                volume_columns[5].markdown('**Available**')
                volume_columns[8].markdown('**Available**')

                for k,p in enumerate(sorted(primer_names)):
                    r = (k*3)%9 #horizontal offset
                    if p != '':
                        volume_columns[r+0].write(p)
                        volume_columns[r+1].write(primer_vols.get(p,0)/1000)
                        volume_columns[r+2].write(primer_avail_vols.get(p,0)/1000)

            #Volume available for each primer
            with ptab4:
                avail_columns = st.columns([3, 3, 3, 3, 3, 3])
                primer_names = set(primer_vols.keys())
                avail_columns[0].markdown('**Primer**')
                avail_columns[2].markdown('**Primer**')
                avail_columns[4].markdown('**Primer**')
                avail_columns[1].markdown('**Available μL**')
                avail_columns[3].markdown('**Available μL**')
                avail_columns[5].markdown('**Available μL**')
                for k,p in enumerate(sorted(primer_names)):
                    r = (k*2)%6 #horizontal offset
                    if p != '':
                        avail_columns[r+0].write(p)
                        avail_columns[r+1].write(primer_avail_vols.get(p,0)/1000)
                       


def clear_widget(*keys):
    for k in keys:
        if k in st.session_state:
            st.session_state[k] = ''



def folder_sb():
    """
    Home page

    Selectbox for loading or creating new experiment
    """
    #New experiment

    ftab1, ftab2 = st.sidebar.tabs(["New Folder", "Existing Folder"])
    if 'experiment' not in st.session_state or not st.session_state['experiment']:
        with ftab1:
            add_run_folder = st.text_input('Create new run folder:')
            create_run_folder_button = st.button(label='Create', key='create_run_folder_button')
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

    #Existing experiments
    #if not st.session_state['experiment'] or (st.session_state['experiment']):
        with ftab2:
            try:
                existing_run_folders = get_run_folders()
            except Exception as exc:
                output_error(exc, msg='Cannot locate NGSgeno folder')
                return
            run_folder = st.selectbox("Select a run folder to open", existing_run_folders)
            if run_folder:
                if not st.session_state['experiment'] or st.session_state['experiment'].name != run_folder:   
                    ch_run_path = 'run_' + run_folder
                    if os.path.exists(ch_run_path):
                        exp = load_experiment(ch_run_path)
                        #print(dir(exp))
                        if not exp:
                            st.markdown('<p style="color:#FF0000">Could not load experiment from: ' + ch_run_path + '</p>', unsafe_allow_html=True)
                        elif ch_run_path.endswith(exp.name):
                            # success!
                            st.session_state['experiment'] = exp
                            st.experimental_rerun()
                        else:
                            st.markdown('<p style="color:#FF0000">Invalid experiment file in: ' + ch_run_path + '</p>', unsafe_allow_html=True)
                ch_exp_folder = None

def pipeline_sb():
    """
    Pipeline page
    Drop down box for all of the option for the stages of the pipeline
    """

    pipeline_stage = st.sidebar.selectbox("Pipeline Stage", ("1. Nimbus Files", "2. Echo Primers (PCR 1)", "3. Echo Barcodes (PCR 2)", "4. Miseq", "5. Genotyping", "6. Review"))


    if pipeline_stage == "1. Nimbus Files":
        st.session_state['pipe_stage'] = 2


    if pipeline_stage == "2. Echo Primers (PCR 1)":
        st.session_state['pipe_stage'] = 3

    if pipeline_stage == "3. Echo Barcodes (PCR 2)":
        st.session_state['pipe_stage'] = 4

    if pipeline_stage == "4. Miseq":
        st.session_state['pipe_stage'] = 5

    if pipeline_stage == "5. Genotyping":
        st.session_state['pipe_stage'] = 6

    if pipeline_stage == "6. Review":
        st.session_state['pipe_stage'] = 7


def sidebar_display(page):
    """
    Both pages

    All of the sidebar images, text, buttons and text inputs
    Args:
        page(int): indicate which page / .py file the website is on
    """

    st.sidebar.title('NGS Genotyping Pipeline')

    st.sidebar.image('ngsg_explorer.png')

    if 'experiment' in st.session_state and st.session_state['experiment']:
        st.sidebar.header('Current experiment: '+ st.session_state['experiment'].name)
        if page == 1:
            st.sidebar.markdown('<p style="color:#6f9db3;font-size:14px">Continue to Pipeline once data is loaded</p>', unsafe_allow_html=True)
            st.sidebar.write('')
            st.sidebar.markdown('<p style="color:#6f9db3;font-size:14px">Or change experiment:</p>', unsafe_allow_html=True)
            #local_css('style.css')

        print("Starting NGS Genotyping Explorer")
    
    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')

    #Home page 
    if page == 1:
        folder_sb()

    #Pipeline
    if page == 2:

        if st.session_state['experiment']:               
            exp = st.session_state['experiment']
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

                pipeline_sb()


def header_selectbox():
    """
    Home page

    Different view options for the data when viewing it
    """
    data_option = st.selectbox('View Data options',
                                ('Summary', 'Edit table', 'Plate view', 'Explore Assays', 'View log'))

    st.write('Selected:', data_option)

    if 'view option' not in st.session_state:
        st.session_state['view option'] = 1
    elif data_option=='Summary' and not st.session_state['view option'] == 1:
        st.session_state['view option'] = 1
    elif data_option=='Edit table' and not st.session_state['view option'] == 2:
        st.session_state['view option'] = 2
    elif data_option=='Plate view' and not st.session_state['view option'] == 3:
        st.session_state['view option'] = 3

    elif data_option=='Explore assays' and not st.session_state['view option'] == 4:
        st.session_state['view option'] = 4

    elif data_option=='View log' and not st.session_state['view option'] == 5:
        st.session_state['view option'] = 5


def data_table(table_options, show_usage, show_primers):
    """
    Home page (ish)

    Shows the loaded data in a table.
    Args:
        show_usage(bool): If True, show the usage data as an expanded widget. Otherwise, don't show.
        show_primers(bool): If True, show the primer data, otherwise don't show it.
    """
    if table_options:
        header_selectbox()
    else:
        st.session_state['view option'] = 1

    #Experiment summary - table
    if st.session_state['view option'] == 1:           
        selection = []
        st.session_state['delete_selection'] = False
        st.session_state['nuke'] = ''
        if st.session_state['experiment']:
            #st.write(st.session_state['experiment'].name)
            df = st.session_state['experiment'].inputs_as_dataframe()
            if df is None or not isinstance(df, pd.DataFrame):
                st.write('No data loaded')
            else:
                #print(f"{type(df)=}, {df=})")
                selection = aggrid_interactive_table(df)
                #print(f"{type(selection)=} {selection=}")
                with st.container():
                    if selection['selected_rows']:
                        rows = selection["selected_rows"]
                        lines = '\n'.join(['DNA PID: '+r['DNA PID'] for r in rows if r['DNA PID'] != 'Total'])
                        if lines:
                            st.write("You selected:")
                            st.write(lines)
                            if st.button('DELETE selection', key='delete_button', help='Remove these entries'):
                        #        st.session_state['delete_selection'] = True
                        #if st.session_state['delete_selection']:    
                                del_col1, del_col2, del_col3 = st.columns([2,1,1])
                                with del_col1:
                                    st.warning("Are you sure you wish to permanently remove these entries?")
                                with del_col2:
                                    mistake_button = st.button("No, it was a mistake!",on_click=delete_entries,
                                            args=(st.session_state['experiment'], []))
                                    if mistake_button:
                                        st.experimental_rerun()
                                with del_col3:
                                    delete_button = st.button("Yes, DELETE them!",on_click=delete_entries,
                                            args=(st.session_state['experiment'], selection['selected_rows']))
                                    if delete_button:
                                        st.experimental_rerun()
            assay_usage = st.session_state['experiment'].get_assay_usage()
            display_reaction_stats(assay_usage=assay_usage, show_general=show_usage, show_primers=show_primers)


    elif st.session_state['view option'] == 2:
        # Show plates in table form and let the user edit them to fix minor mistakes
        pass

    #Show a visual of the plates
    elif st.session_state['view option'] == 3:
        plate_input = st.text_input('Plate ID to view', key='plate_input1')
        if plate_input:
            plate_id = file_io.guard_pbc(plate_input)
            if plate_id in st.session_state['experiment'].plate_location_sample:
                heatmap_str = generate_heatmap_html(st.session_state['experiment'], plate_id, scaling=1.2)
                #with open("debug.html", 'wt') as outf:
                #    print(heatmap_str, file=outf)
                components.html(heatmap_str, height=700, scrolling=True)
            else:
                st.markdown('<p style="color:#FF0000">Plate barcode not found in experiment</p>', unsafe_allow_html=True)
        # let user choose plates and then display them as a plate layout with well contents on hover
        # plate_col1, plate_col2 = st.columns(2)
        # for plate in st.session_state['experiment'].plate_location_sample:
        #    break
        # heatmap_str = generate_heatmap_html(plate_id, st.session_state['experiment'], scaling=1)
        # components.html(heatmap_str, height=600, scrolling=True)

    elif st.session_state['view option'] == 4:  # Assay explorer
        # display allowed and denied assays, missing assays, used assays, etc.
        assay_col1, assay_col2, assay_col3, assay_col4, assay_col5 = st.columns(5)
        with assay_col1:
            pass
            #rodentity_assays_avail = st.session_state['experiment'].get_rodentity_assays_avail()
            #st.write('Rodentity assays available')
            #for a in rodentity_assays_avail:
            #    st.write(a)
        with assay_col2:
            pass
        with assay_col3:
            pass
        with assay_col4:
            pass
        with assay_col5:
            pass
                
    elif st.session_state['view option'] == 5:
        # display the experiment log and let the user filter the view in a number of ways
        with st.container():
            df = pd.DataFrame(st.session_state['experiment'].get_log(30))
            if 'view_box_size' in st.session_state:
                height = st.session_state['view_box_size']
            else:
                height = 250
            st.dataframe(df, height=height)

def extra_data():
    """
    Home page: load data

    Upload options for extra custom data including reference files, assay lists, primer plates and barcode plates
    """
    with st.container():
        with st.expander('Upload reference sequences, assay lists, taq and water plates, primer plates, and barcode plates', expanded=False):

            upload_option = st.radio('Upload Options:', ('Custom Reference File', 'Assay List', 'Taq and Water Plates', 'Primer Plates', 'Barcode Plates'))

            if upload_option == 'Custom Reference File':
                uploaded_reference = st.file_uploader('Upload custom reference file', key='ref_uploader', type=['txt','fa','fasta'])
                if uploaded_reference:
                    if 'reference_upload' not in st.session_state or st.session_state['reference_upload'] != uploaded_reference.name:
                        success = st.session_state['experiment'].add_reference(uploaded_reference)
                        st.session_state['reference_upload'] = uploaded_reference.name
                

            if upload_option == 'Assay List':
                uploaded_assaylist = st.file_uploader('Upload assay list', key='assaylist_uploader', type=['txt','csv'])
                if uploaded_assaylist:
                    if 'assaylist_upload' not in st.session_state or st.session_state['assaylist_upload'] != uploaded_assaylist.name:
                        success = st.session_state['experiment'].add_assaylist(uploaded_assaylist)
                        st.session_state['assaylist_upload'] = uploaded_assaylist.name

            if upload_option == 'Primer Plates':
                uploaded_primer_layouts = st.file_uploader('Upload primer plate layouts', 
                        key='primer_uploader', type='csv', accept_multiple_files=True)
                if uploaded_primer_layouts:
                    if 'primer_layouts_upload' not in st.session_state or st.session_state['primer_layouts_upload'] != uploaded_primer_layouts[0].name:
                        success = st.session_state['experiment'].add_primer_layouts(uploaded_primer_layouts)
                        st.session_state['primer_layouts_upload'] = uploaded_primer_layouts[0].name
                uploaded_primer_volumes = st.file_uploader('Upload primer plate volumes', key='primer_vol_uploader', type='csv', accept_multiple_files=True)
                if uploaded_primer_volumes:
                    if 'primer_volumes_upload' not in st.session_state or st.session_state['primer_volumes_upload'] != [upv.name for upv in uploaded_primer_volumes]:
                        success = st.session_state['experiment'].add_primer_volumes(uploaded_primer_volumes)
                        st.session_state['primer_volumes_upload'] = [upv.name for upv in uploaded_primer_volumes]

            if upload_option == 'Barcode Plates':
                uploaded_barcode_layouts = st.file_uploader('Upload i7i5 barcode plate layout', 
                        key='barcode_uploader', type='csv', accept_multiple_files=True)
                if uploaded_barcode_layouts:
                    if 'barcode_layout_upload' not in st.session_state or st.session_state['barcode_layout_upload'] != uploaded_barcode_layouts[0].name:
                        success = st.session_state['experiment'].add_barcode_layouts(uploaded_barcode_layouts)
                        st.session_state['barcode_layout_upload'] = uploaded_barcode_layouts[0].name
                            
                uploaded_barcode_volumes = st.file_uploader('Upload i7i5 barcode plate volumes', key='barcode_vol_uploader',
                        type='csv', accept_multiple_files=True)
                if uploaded_barcode_volumes:
                    if 'barcode_volume_upload' not in st.session_state or \
                            st.session_state['barcode_volume_upload'] !=\
                            [ubv.name for ubv in uploaded_barcode_volumes]:
                        success = st.session_state['experiment'].add_barcode_volumes(uploaded_barcode_volumes)
                        st.session_state['barcode_volume_upload'] = [ubv.name for ubv in uploaded_barcode_volumes]

            #Option to add taq water
            if upload_option == 'Taq and Water Plates':
                uploaded_taqwater_plates = st.file_uploader('Upload taq and water resevoir plates', 
                        key='taq_water_upload', type='csv', accept_multiple_files=True)
                if uploaded_taqwater_plates:
                    if 'taqwater_upload' not in st.session_state or \
                        st.session_state['taqwater_upload'] != uploaded_taqwater_plates[0].name:
                        success = st.session_state['experiment'].add_standard_taqwater_plates(uploaded_taqwater_plates) #not sure if this is right
                        st.session_state['taqwater_upload'] = uploaded_taqwater_plates[0].name #????????



                

def load_data():
    """
    Home page 
    """
    with st.expander('Add data from Rodentity JSON files',expanded=True):
        # C1: upload, C2: four slots plus 'X' remove, C3: dest plate; "accept" button
        
        rod_dp = ''
        plates_to_clear = [False, False, False, False]
        rod_up_disabled = False  # disable uploading more plates until space available
        exp = st.session_state['experiment']
        #print(f"{exp.unassigned_plates=}")
        #if 'prev_rod_upload' in st.session_state:
        #    print(f"{st.session_state['prev_rod_upload']=}")
        #if 'rod_dp_key' in st.session_state:
        #    print(f"{st.session_state['rod_dp_key']=}")
        if all(exp.unassigned_plates[k] for k in exp.unassigned_plates):
            rod_up_disabled = True
                

        rodentity_epp = st.file_uploader('Rodentity JSON file', type='json', disabled=rod_up_disabled)
        if rodentity_epp:
            rodentity_plate_name = file_io.guard_pbc(rodentity_epp.name.rstrip('.json'))
            if 'prev_rod_upload' not in st.session_state or \
                st.session_state['prev_rod_upload'] != rodentity_plate_name:
                if rodentity_epp.name.rstrip('.json') in exp.unassigned_plates:
                    print('We see this already')
                else: 
                    for k in sorted(exp.unassigned_plates):
                        if not exp.unassigned_plates[k]:
                            exp.unassigned_plates[k] = rodentity_plate_name
                            # Save a copy we can edit and reload
                            fn = os.path.join('run_'+exp.name, db_io._eppfn_r(rodentity_plate_name))
                            with open(fn, 'wt') as outf:
                                outf.write(rodentity_epp.getvalue().decode("utf-8"))
                            print('here', k)
                            st.session_state['prev_rod_upload'] = rodentity_plate_name
                            st.experimental_rerun()     

        col_r1, col_r2= st.columns(2)
        with col_r1:
            st.write('Checkboxes required for clearing plate IDs')
                    
            for i,p in enumerate(plates_to_clear):
                plates_to_clear[i] = st.checkbox(f"P{str(i+1)}: {file_io.unguard_pbc(exp.unassigned_plates[i+1], silent=True)}", 
                        help='Click the checkbox to allow a set plate ID to be cleared', key='check'+str(i+1))
        with col_r2:
            accept_disabled = False
            rod_dp = st.text_input('Destination plate barcode', max_chars=20, key='rod_dp_key')
            if rod_dp and any(exp.unassigned_plates):
                rod_dp = file_io.guard_pbc(rod_dp, silent=True)
                if rod_dp in st.session_state['experiment'].dest_sample_plates or \
                        rod_dp in st.session_state['experiment'].plate_location_sample or \
                        rod_dp in st.session_state['experiment'].unassigned_plates:
                    st.markdown('<p style="color:#FF0000">Destination plate barcode already in use: ' +\
                            file_io.unguard_pbc(rod_dp) + '</p>', unsafe_allow_html=True)
            else:
                accept_disabled = True

            accept_rod_button = st.button('Accept', help='Read and confirm plates', 
                    disabled=accept_disabled, key='accept_rod_button_key')
            if accept_rod_button:
                success = exp.add_rodentity_plate_set([exp.unassigned_plates[k] for k in exp.unassigned_plates], rod_dp)
                if not success:
                    st.markdown('<p style="color:#FF0000">Failed to incorporate plate set. Please read the log.</p>', unsafe_allow_html=True)
                    if not success:
                        st.session_state['view option index'] = len(view_headers)-1
                else:
                    st.write('Successfully added plate set')
                    rod_dp = ''
                    exp.unassigned_plates = {1:'',2:'',3:'',4:''}
                    clear_widget(('rod_dp_key','prev_rod_upload', ))
                    plates_to_clear = [False, False, False, False]
                    accept_rod_button = False
                    exp.save()
                    st.experimental_rerun()
            clear_disabled = True
            for i,plate in enumerate(plates_to_clear):
                if plate and exp.unassigned_plates[i+1]:
                    clear_disabled = False
                    break  # only need to enable the button once
            clear_plates_button = st.button('Clear IDs', help='Clear selected Rodentity plate IDs', disabled=clear_disabled)
            cpb_activated = False
            if clear_plates_button:
                for i, plate in enumerate(plates_to_clear):
                    if plate and exp.unassigned_plates[i+1]:
                        if 'prev_rod_upload' in st.session_state:
                            if st.session_state['prev_rod_upload'] == plate:
                                st.session_state['prev_rod_upload'] = ''
                        exp.unassigned_plates[i+1] = ''
                        plates_to_clear[i] = False
                        cpb_activated = True
                clear_plates_button = False
                if cpb_activated:
                    exp.save()
                    st.experimental_rerun()


    with st.container():
        with st.expander('Add data from a custom manifest CSV file', expanded=False):
            col_m1, col_m2 = st.columns(2)
            st.session_state['default_manifest_type'] = 'c'
            with col_m1:
                manifest_choice = st.radio('The default contents of my manifest file (if not declared) are', 
                        ['Custom', 'Rodentity mouse', 'Musterer mouse'])
                if manifest_choice == 'Rodentity mouse':
                    st.session_state['default_manifest_type'] = 'r'
                elif manifest_choice == 'Musterer mouse':
                    st.session_state['default_manifest_type'] = 'm'
            with col_m2:
                uploaded_file = st.file_uploader('', type='csv')
                if uploaded_file:
                    if 'custom_upload' not in st.session_state or st.session_state['custom_upload'] != uploaded_file.name:
                        success = st.session_state['experiment'].add_manifest(uploaded_file, st.session_state['default_manifest_type'])
                        if not success:
                            st.markdown('<p style="color:#FF0000">Failed to incorporate manifest. Please read the log.</p>', unsafe_allow_html=True)
                            st.session_state['view option index'] = len(view_headers) -1
                        st.session_state['custom_upload'] = uploaded_file.name
                        st.experimental_rerun()

   
    with st.container():
        with st.expander('Add data from database'):
            #st.subheader('Add mouse plate set')
            row1col1, row1col2, row1col3, row1col4 = st.columns(4)
            epp1 = ''
            epp2 = ''
            epp3 = ''
            epp4 = ''
            dnap = ''
            with row1col1:
                epp1_str = st.text_input(label='', placeholder='Ear punch plate 1 barcode', key='epp1_input').strip()
                if epp1_str:
                    try:
                        epp1 = file_io.guard_pbc(epp1_str)
                    except Exception as exc:
                        output_error(exc, msg="Ear punch plate 1 barcode")         
            with row1col2:
                epp2_str = st.text_input(label='', placeholder='Ear punch plate 2 barcode', key='epp2_input').strip()
                if epp2_str:
                    try:
                        epp2 = file_io.guard_pbc(epp2_str)
                    except Exception as exc:
                        output_error(exc, msg="Ear punch plate 2 barcode")              
            with row1col3:
                epp3_str = st.text_input(label='', placeholder='Ear punch plate 3 barcode', key='epp3_input').strip()
                if epp3_str:
                    try:
                        epp3 = file_io.guard_pbc(epp3_str)
                    except Exception as exc:
                        output_error(exc, msg="Ear punch plate 3 barcode")          
            with row1col4:
                epp4_str = st.text_input(label='', placeholder='Ear punch plate 4 barcode', key='epp4_input').strip()
                if epp4_str:
                    try:
                        epp4 = file_io.guard_pbc(epp4_str)
                    except Exception as exc:
                        output_error(exc, msg="Ear punch plate 4 barcode")      
                
            row2col1, row2col2, row2col3, row2col4 = st.columns(4)
            with row2col1:
                dnap_str = st.text_input(label='', placeholder='DNA plate barcode', key='dnap_input').strip()
                if dnap_str:
                    try:
                        dnap = file_io.guard_pbc(dnap_str)
                    except Exception as exc:
                        output_error(exc, msg="DNA plate barcode")
            # nothing in row2col2
            with row2col3:
                st.markdown('')
                st.markdown('')
                add_musterer_button = st.button('Search Musterer')
                if add_musterer_button:
                    epps = [epp.strip() for epp in [epp1,epp2,epp3,epp4] if epp.strip() != '']
                    for epp in epps:
                        success = db_io.get_plate_musterer(st.session_state['experiment'].name, epp, True)
                    success = st.session_state['experiment'].add_musterer_plate_set(epps,dnap)
                    if not success:
                        st.markdown('<p style="color:#FF0000">Failed to incorporate manifest. Please read the log.</p>', unsafe_allow_html=True)
                        st.session_state['view option index'] = len(view_headers)-1
                    st.experimental_rerun()

            with row2col4:
                st.markdown('')
                st.markdown('')


def pipe_stages(exp, stage):
    '''
    Pipeline page

    Goes through each stage of the pipeline, depending on the session state.
    Args:
        int stage: st.session_state['pipe_stage']
    '''

    if stage == 1:
        with st.container():
            data_table(options=True, show_usage=True, show_primers=True)



    #Nimbus: 96-wells of sample plates to 384-wells of DNA plates
    if stage == 2:
        with st.container():
            st.markdown('<h2 style="text-align:center;color:#0f6b8e">Nimbus</h2>', unsafe_allow_html=True)
            st.markdown('<h5 style="text-align:center;color:#2BA2D0">Nimbus input files have been generated from your data</h5>', unsafe_allow_html=True)
            #st.write('Combine 96-well sample plates to 384-well DNA plate.')
            st.write('')
            st.write('')
            
            if not st.session_state['experiment'].dest_sample_plates:
                st.markdown('<p style="color:#FF0000">Load data inputs to enable Nimbus input file generation.</p>', unsafe_allow_html=True)
            else:
                nim_tab1, nim_tab2 = st.tabs(["1. Download", "2. Upload"])



                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                with nim_tab1:

                    # do we have any Nimbus inputs to generate + download
                    yet_to_run = len(st.session_state['experiment'].dest_sample_plates) - len(nfs)
                    if yet_to_run and yet_to_run > 0: #and yet_to_run > 0
                        st.write(f"{yet_to_run} 96-well plate sets need Nimbus input file generation")
                        plates_to_run = [dest_plate for dest_plate in exp.dest_sample_plates if all([dest_plate not in nf for nf in nfs])]
                        plates_to_run_str = '\n'.join(plates_to_run)
                        st.write(plates_to_run_str)
                        success = st.session_state['experiment'].generate_nimbus_inputs()
                        if not success:
                            st.markdown('<p style="color:#FF0000">Failed to generate Nimbus files. Please read the log.</p>', unsafe_allow_html=True)
                            st.session_state['view option index'] = len(view_headers)-1
                        else:
                            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                    if nfs:
                        st.markdown('<h5 style="text-align:center">1. Download Nimbus input files</h5>', unsafe_allow_html=True)
                        st.markdown('<p style="text-align:center">Download the 96-well sample plates as input files for Nimbus to run.</p>', unsafe_allow_html=True)
                        st.write('')
                        st.write('')
                        _,dl_col1,dl_col2,dl_col3,dl_col4,_=st.columns([1,9,6,9,6,1])
                        for i,nf in enumerate(nfs):
                            file_name=nf.split('\\')[1]
                            if (i+1) % 2 != 0:
                                dl_col1.markdown(f"<p style='text-align:left;color:#4b778c;padding:5px'>{file_name}</p>", unsafe_allow_html=True)
                                #dl_cols[0].write('')
                                dl_col2.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                        key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                                #dl_cols[2].write('')
                            else:
                                dl_col3.markdown(f"<p style='text-align:left;color:#4b778c;padding:5px''>{file_name}</p>", unsafe_allow_html=True)
                                #dl_cols[0].write('')
                                dl_col4.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                        key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                                #dl_cols[2].write('')
                        
                with nim_tab2:
                    nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                    print (f"{nfs=} {efs=} {xbcs=}")
                    # Nimbus upload
                    if not efs and not xbcs:
                        st.markdown('<h4 style="text-align=center">Awaiting Nimbus output files</h4>', unsafe_allow_html=True)
                    elif efs and not xbcs:
                        st.markdown('<h4 style="text-align:center">All Nimbus outputs now uploaded</h4>', unsafe_allow_html=True) 
                    elif xbcs:
                        st.markdown('<h5 style="text-align:center">2. Upload Nimbus output files</h5>', unsafe_allow_html=True)
                        st.markdown('<p style="text-align:center">Upload the resulting 384-well DNA plate files for Echo.</p>', unsafe_allow_html=True)
                        # list nimbus outputs that we still need
                        st.markdown('<h6 style="color:#8e0f5b">Missing the Nimbus output files:</h6>', unsafe_allow_html=True)
                        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                        missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
                        mn_str = '</br>'.join(missing_nims)
                        st.markdown(f"<p style='color:#87adc7'>{mn_str}</p>", unsafe_allow_html=True)
                        if 'nim_upload' not in st.session_state:
                            st.session_state['nim_upload'] = ''
                        nim_outputs = st.file_uploader(' ', type='csv', 
                                accept_multiple_files=True, help='You can upload more than one file at once')
                        if nim_outputs: # and nim_outputs != st.session_state['nim_upload']:
                            #st.session_state['nim_uploads'] = nim_outputs
                            for nim_output in nim_outputs:
                                #print(f"{nim_output.name=}")
                                fp = st.session_state['experiment'].get_exp_fp(nim_output.name)
                                print(f"Copying {fp} to experiment folder")
                                with open(fp, 'wt') as outf:
                                    #print(nim_output.getvalue().decode("utf-8"))
                                    outf.write(nim_output.getvalue().decode("utf-8").replace('\r\n','\n'))
                        #else:
                        #    st.session_state['nim_uploads'] = ''
                                    
                    
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            if not xbcs:
                st.markdown('<h4 style="color:#FF0000;text-align:center">Ready to run Echo primer stage</h4>', unsafe_allow_html=True)

    if stage == 3:
        #Primer stage
        
        with st.container():
            st.markdown('<h2 style="text-align:center;color:#0f6b8e">Echo Primers (PCR 1)</h2>', unsafe_allow_html=True)
            with st.expander('Primer plate details'):
                pass
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]
            missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
            #All of the nimbus output files haven't been uploaded
            if missing_nims:
                st.markdown('<h5 style="text-align:center;color:#B92A5D">Missing the Nimbus Output Files:</h5>', unsafe_allow_html=True)
                for mn in missing_nims:
                    st.markdown(f'<p style="text-align:center">{mn}</p>', unsafe_allow_html=True)
                st.markdown('<p style="text-align:center;color:#B92A5D">Upload or continue without them</p>', unsafe_allow_html=True)
  
                if 'nim_upload' not in st.session_state:
                        st.session_state['nim_upload'] = ''
                nim_outputs = st.file_uploader('Upload Nimbus output files', type='csv', 
                        accept_multiple_files=True, help='You can upload more than one file at once')
                if nim_outputs: # and nim_outputs != st.session_state['nim_upload']:
                    #st.session_state['nim_uploads'] = nim_outputs
                    for nim_output in nim_outputs:
                        #print(f"{nim_output.name=}")
                        fp = st.session_state['experiment'].get_exp_fp(nim_output.name)
                        print(f"Copying {fp} to experiment folder")
                        with open(fp, 'wt') as outf:
                            #print(nim_output.getvalue().decode("utf-8"))
                            outf.write(nim_output.getvalue().decode("utf-8").replace('\r\n','\n'))

            if available_nimbus:
                st.write('')
                st.write('')
                included_DNA_plates = set()
                st.markdown('<h4 style="color:#0f6b8e;text-align:center">Choose Nimbus plates for Echo (PCR 1 Assay)</h4>', unsafe_allow_html=True)
                st.write('')
                for nim in available_nimbus:
                    plate=nim.split('\\')[1].split('.')[0]
                    inc = st.checkbox(plate, value=False, key='chk_box_'+nim)
                    if inc:
                        included_DNA_plates.add(nim)
                assay_usage = st.session_state['experiment'].get_assay_usage(dna_plate_list=included_DNA_plates)
                st.write('')
                display_reaction_stats(assay_usage)
                extra_data()

