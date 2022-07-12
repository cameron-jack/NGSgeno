#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New streamlit interactive web interface
@last_edit: 2022-05-16
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

def aggrid_interactive_table(df: pd.DataFrame, grid_height: int=250):
    """Creates an st-aggrid interactive table based on a dataframe.

    Args:
        df (pd.DataFrame]): Source dataframe

    Returns:
        dict: The selected row
    """
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
        theme="fresh",
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
    run_folders = [''] + [d[4:] for d in os.listdir(st.session_state['app_path']) if d.startswith('run_') and os.path.isdir(d)]
    return sorted(run_folders)


def create_run_folder(new_run_name):
    """ returns Experiment or None, and a message string. Requires a full path and generates a new experiment file """
    newpath = os.path.join(st.session_state['app_path'], new_run_name)
    print('Attempting to create: ', newpath)
    if not os.path.exists(newpath):
        try:
            os.mkdir(newpath)
        except Exception as exc:
            output_error(exc, msg='Could not create new folder: ' + newpath)
            return None, 'Failed to create new folder: ' + newpath
    else:
        if os.path.exists(os.path.join(newpath, EXP_FN)):
            return None, 'Experiment already exists with this name: ' + new_run_name
    print('Generating experiment: ', new_run_name)
    exp = Experiment(name=new_run_name.lstrip('run_'))
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


def clear_widget(*keys):
    for k in keys:
        if k in st.session_state:
            st.session_state[k] = ''


### Lay out the pipeline interface sections
# 1. Run folder (includes sidebar)
# 2. Inventory (plates and assays + DB lookup buttons + remaining i7i5 barcodes)
# 3. Nimbus progress
# 4. Primer plates + taq/water
# 5. Add spare samples for barcoding + barcode completion status (plates done/not done)
# 6. Miseq data available (notice + instructions if not yet complete)
# 7. Analysis stages
###
 

def main():
    st.set_page_config(page_icon="ngsg_icon.png", page_title="NGS Genotyping Pipeline", initial_sidebar_state="expanded", layout="wide")
    # Remove whitespace from the top of the page and sidebar
    st.markdown("""
        <style>
               .css-18e3th9 {
                    padding-top: 0rem;
                    padding-bottom: 10rem;
                    padding-left: 5rem;
                    padding-right: 5rem;
                }
        </style>
        """, unsafe_allow_html=True)
    #This is a workaround from https://www.codegrepper.com/code-examples/python/streamlit+sidebar+width
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
            width: 450px;
        }
        [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
            width: 500px;
            margin-left: -500px;
        }
        </style>
        """,
        unsafe_allow_html=True)
    
    st.sidebar.title('NGS Genotyping Pipeline')
    st.sidebar.image('ngsg_explorer.png')
    if 'experiment' in st.session_state and st.session_state['experiment']:
        st.sidebar.header('Current experiment: '+ st.session_state['experiment'].name)
    
    print("Starting NGS Genotyping Explorer")
    if 'app_path' not in st.session_state:
        app_path = str(PurePath(os.getcwd()))
        if 'NGSgeno' not in app_path:
            st.write('NGSgeno not in application path, please restart the app')
            st.stop()
        app_path = app_path.split('NGSgeno')[0] + 'NGSgeno'
        st.session_state['app_path'] = app_path
        print('Initialising app path to:', app_path)

    os.chdir(os.path.join(st.session_state['app_path']))
    print("Loading library defaults...")
    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')

    if 'nav_pipe_button' not in st.session_state:
        st.session_state['nav_pipe_button'] = 'nav'

    if st.session_state['experiment']:
        sb_col1, sb_col2, sb_col3, sb_col4, sb_col5, sb_col6 = st.sidebar.columns([1,2,1,1,2,1])
        nav_clicked = sb_col2.button('Navigation', help='Change experiment, or create a new one')
        if nav_clicked:
            st.session_state['nav_pipe_button'] = 'nav'

        pipe_clicked = sb_col5.button('Pipeline', help='Work through pipeline stages')
        if pipe_clicked:
            st.session_state['nav_pipe_button'] = 'pipe'

        if st.session_state['nav_pipe_button'] == 'nav':
            sb_col1.button('\u25B6', key='nav_pipe_col1_button', disabled=True)
            sb_col3.button('\u25C0', key='nav_pipe_col3_button', disabled=True)
            sb_col4.button('.', key='nav_pipe_col4_button', disabled=True)
            sb_col6.button('.', key='nav_pipe_col6_button', disabled=True)
        elif st.session_state['nav_pipe_button'] == 'pipe':
            sb_col4.button('\u25B6', key='nav_pipe_col4_button', disabled=True)
            sb_col6.button('\u25C0', key='nav_pipe_col6_button', disabled=True)
            sb_col1.button('.', key='nav_pipe_col1_button', disabled=True)
            sb_col3.button('.', key='nav_pipe_col3_button', disabled=True)
            
        

    ### TODO: get paths for mouse files in the library folder
    # Preferably load the data once, rather than each time a script is run - maybe need a "refresh library" button?

    # Run build folder navigation, load run folder and experiment
    if not st.session_state['experiment'] or (st.session_state['experiment'] and st.session_state['nav_pipe_button'] == 'nav'):
        add_run_folder = st.sidebar.text_input('Create new run folder:')
        
        create_run_folder_button = st.sidebar.button(label='Create', key='create_run_folder_button')
        if add_run_folder and create_run_folder_button:
            add_run_folder_str = 'run_' + add_run_folder
            exp, msg = create_run_folder(add_run_folder_str)
            if exp:
                exp.save()                                                                      
                st.session_state['experiment'] = exp
                st.experimental_rerun()
            else:
                if 'already exists' in msg:
                    st.sidebar.markdown('<p style="color:#FF0000">Folder name already exists</p>', unsafe_allow_html=True)
                else:
                    st.markdown('<p style="color:#FF0000">Fatal path error: ' + msg + '</p>', unsafe_allow_html=True)

        try:
            existing_run_folders = get_run_folders()
        except Exception as exc:
            output_error(exc, msg='Cannot locate NGSgeno folder')
            st.stop()

        ch_exp_folder = st.sidebar.selectbox("Select a run folder to open", existing_run_folders)
        if ch_exp_folder:
            if not st.session_state['experiment'] or st.session_state['experiment'].name != ch_exp_folder:   
                ch_run_folder = 'run_' + ch_exp_folder
                ch_run_path = os.path.join(st.session_state['app_path'], ch_run_folder)
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
    
    if st.session_state['experiment']:               
        exp = st.session_state['experiment']
        if 'pipe_stage' not in st.session_state:
            st.session_state['pipe_stage'] = 1
            
        if st.session_state['nav_pipe_button'] == 'nav':
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
            _,sb_pipe_col1, sb_pipe_col2, sb_pipe_col3,_ = st.sidebar.columns([1,1,8,1,2])
            #if 'pipe_stage' not in st.session_state:
            #    st.session_state['pipe_stage'] = 1
            #    sb_pipe_col1.button('\u25B6', key='1l', disabled=True)
            #    sb_pipe_col3.button('\u25C0', key='1r', disabled=True)
            # Load data
            pipe_load_data_clicked = sb_pipe_col2.button('1. Load data', key='1c')
            pipe_nimbus_clicked = sb_pipe_col2.button('2. Nimbus 96->384', key='2c')
            pipe_echo_primers_clicked = sb_pipe_col2.button('3. Echo Primers', key='3c')
            pipe_echo_barcode_clicked = sb_pipe_col2.button('4. Echo Barcodes', key='4c')
            pipe_miseq_clicked = sb_pipe_col2.button('5. Miseq', key='5c')
            pipe_genotyping_clicked = sb_pipe_col2.button('6. Genotyping', key='6c')
            pipe_review_clicked = sb_pipe_col2.button('7. Review', key='7c')
            
            if 'pipe_stage' not in st.session_state or pipe_load_data_clicked:
                st.session_state['pipe_stage'] = 1
                sb_pipe_col1.button('\u25B6', key='1l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='1r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='1l', disabled=True)
                sb_pipe_col3.button('.', key='1r', disabled=True)
            
            if pipe_nimbus_clicked:
                st.session_state['pipe_stage'] = 2
                sb_pipe_col1.button('\u25B6', key='2l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='2r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='2l', disabled=True)
                sb_pipe_col3.button('.', key='2r', disabled=True)
            
            if pipe_echo_primers_clicked:
                st.session_state['pipe_stage'] = 3
                sb_pipe_col1.button('\u25B6', key='3l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='3r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='3l', disabled=True)
                sb_pipe_col3.button('.', key='3r', disabled=True)
            
            if pipe_echo_barcode_clicked:
                st.session_state['pipe_stage'] = 4
                sb_pipe_col1.button('\u25B6', key='4l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='4r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='4l', disabled=True)
                sb_pipe_col3.button('.', key='4r', disabled=True)
            
            if pipe_miseq_clicked:
                st.session_state['pipe_stage'] = 5
                sb_pipe_col1.button('\u25B6', key='5l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='5r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='5l', disabled=True)
                sb_pipe_col3.button('.', key='5r', disabled=True)
            
            if pipe_genotyping_clicked:
                st.session_state['pipe_stage'] = 6
                sb_pipe_col1.button('\u25B6', key='6l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='6r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='6l', disabled=True)
                sb_pipe_col3.button('.', key='6r', disabled=True)
            
            if pipe_review_clicked:
                st.session_state['pipe_stage'] = 7
                sb_pipe_col1.button('\u25B6', key='7l', disabled=True)
                sb_pipe_col3.button('\u25C0', key='7r', disabled=True)
            else:
                sb_pipe_col1.button('.', key='7l', disabled=True)
                sb_pipe_col3.button('.', key='7r', disabled=True)

                    
   
    #with st.container():
    #    with open('heatmap.html', 'rt') as f:
    #        heatmap_str = f.read()
    #        components.html(heatmap_str, height=600, scrolling=True)

    if not st.session_state['experiment']:
        with st.container():
            st.write('')
            st.write('')
            st.header('Please load existing experiment or create a new one')
    else:
        with st.container():
            view_sp1, view_col1, view_sp2, view_col2, view_sp3, view_col3, view_sp4, view_col4, \
                    view_sp5, view_col5, view_sp6, _, view_resize_col = st.columns([2,8,2,5,2,5,2,6,2,5,2,2,7])
            view_buttons = [view_col1, view_col2, view_col3, view_col4, view_col5]
            view_spacers = [view_sp1, view_sp2, view_sp3, view_sp4, view_sp5, view_sp6]
            for vb in view_buttons:
                vb.write('.')
            for vs in view_spacers:
                vs.write('.')
            # if 'view_option' not in st.session_state:
            #    st.session_state['pipe_stage'] = 1
            #    sb_pipe_col1.button('\u25B6', key='1l', disabled=True)
            #    sb_pipe_col3.button('\u25C0', key='1r', disabled=True)
            # Load data
            view_experiment_summary = view_col1.button('Experiment Summary', key='1v')
            view_table_editor = view_col2.button('Table Editor', key='2v')
            view_plate_layout = view_col3.button('Plate Layout', key='3v')
            view_assay_explorer = view_col4.button('Assay Explorer', key='4v')
            view_log_viewer = view_col5.button('Log Viewer', key='5v')
            view_resize = view_resize_col.number_input('',min_value=10, max_value=1000, value=250, step=20)
            
            ### Update the visual indicators first
            # We could save a lot of code by just creating new buttons with '.' symbols then overwriting with the
            # indicators we want, but this will slow the interface and may make it "blink" more. Instead we update
            # exactly the elements that need to change.
            if 'view option' not in st.session_state:
                st.session_state['view option'] = 1
                view_sp1.button('\u25B6', key='sp1', disabled=True)
                view_sp2.button('\u25C0', key='sp2', disabled=True)
                for i in range(3, len(view_spacers)+1):
                    view_spacers[i-1].button('.', key='sp'+str(i), disabled=True)
            elif view_experiment_summary and not st.session_state['view option'] == 1:
                old_view = st.session_state['view option']
                st.session_state['view option'] = 1
                view_sp1.button('\u25B6', key='sp1', disabled=True)
                view_sp2.button('\u25C0', key='sp2', disabled=True)
                if old_view == st.session_state['view option'] + 1: # next button
                    view_sp3.button('.', key='sp3', disabled=True)
                else:
                    view_spacers[old_view-1].button('.', key='sp'+str(old_view), disabled=True)
                    view_spacers[old_view].button('.', key='sp'+str(old_view+1), disabled=True)
            elif view_table_editor and not st.session_state['view option'] == 2:
                old_view = st.session_state['view option']
                st.session_state['view option'] = 2
                view_sp2.button('\u25B6', key='sp2', disabled=True)
                view_sp3.button('\u25C0', key='sp3', disabled=True)
                if old_view == st.session_state['view option'] - 1: # old view was previous button
                    view_sp1.button('.', key='sp1', disabled=True)            
                elif old_view == st.session_state['view option'] + 1: # old view was next button
                    view_sp4.button('.', key='sp4', disabled=True)
                else:
                    view_spacers[old_view-1].button('.', key='sp'+str(old_view), disabled=True)
                    view_spacers[old_view].button('.', key='sp'+str(old_view+1), disabled=True)
            elif view_plate_layout and not st.session_state['view option'] == 3:
                old_view = st.session_state['view option']
                st.session_state['view option'] = 3
                view_sp3.button('\u25B6', key='sp3', disabled=True)
                view_sp4.button('\u25C0', key='sp4', disabled=True)
                if old_view == st.session_state['view option'] - 1: # old view was previous button
                    view_sp2.button('.', key='sp2', disabled=True)            
                elif old_view == st.session_state['view option'] + 1: # old view was next button
                    view_sp5.button('.', key='sp4', disabled=True)
                else:
                    view_spacers[old_view-1].button('.', key='sp'+str(old_view), disabled=True)
                    view_spacers[old_view].button('.', key='sp'+str(old_view+1), disabled=True)
            elif view_assay_explorer and not st.session_state['view option'] == 4:
                old_view = st.session_state['view option']
                st.session_state['view option'] = 4
                view_sp3.button('\u25B6', key='sp4', disabled=True)
                view_sp4.button('\u25C0', key='sp5', disabled=True)
                if old_view == st.session_state['view option'] - 1: # old view was previous button
                    view_sp2.button('.', key='sp3', disabled=True)            
                elif old_view == st.session_state['view option'] + 1: # old view was next button
                    view_sp5.button('.', key='sp6', disabled=True)
                else:
                    view_spacers[old_view-1].button('.', key='sp'+str(old_view), disabled=True)
                    view_spacers[old_view].button('.', key='sp'+str(old_view+1), disabled=True)
            elif view_log_viewer and not st.session_state['view option'] == 5:
                old_view = st.session_state['view option']
                st.session_state['view option'] = 5
                view_sp3.button('\u25B6', key='sp5', disabled=True)
                view_sp4.button('\u25C0', key='sp6', disabled=True)
                if old_view == st.session_state['view option'] - 1: # old view was previous button
                    view_sp2.button('.', key='sp4', disabled=True)            
                else:
                    view_spacers[old_view-1].button('.', key='sp'+str(old_view), disabled=True)
                    view_spacers[old_view].button('.', key='sp'+str(old_view+1), disabled=True)

            # now do the logic
            if st.session_state['view option'] == 1:  # Experiment Summary            
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
                    with st.expander('Required components: primers, barcodes, taq and water:', expanded=True):
                        assay_usage = st.session_state['experiment'].get_assay_usage()
                        if assay_usage:
                            primer_vols, primer_taq_vol, primer_water_vol, barcode_taq_vol, barcode_water_vol =\
                                    st.session_state['experiment'].get_volumes_required(assay_usage=assay_usage)
                            reactions = sum([v for v in assay_usage.values()])
                            barcodes_remain, barcodes_avail, barcode_vol_okay = st.session_state['experiment'].get_barcode_remaining_available_volume(assay_usage=assay_usage)
                        else:
                            reactions, primer_vols, primer_taq_vol, primer_water_vol, barcode_taq_vol, barcode_water_vol =\
                                0,{},0,0,0,0 # 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
                            barcodes_remain, barcodes_avail, barcode_vol_okay =\
                                    st.session_state['experiment'].get_barcode_remaining_available_volume()
                        primer_avail_counts, primer_avail_vols = st.session_state['experiment'].get_primers_avail()
                        taq_water_avail_PCR1 = st.session_state['experiment'].get_taq_water_avail(1)
                        taq_water_avail_PCR2 = st.session_state['experiment'].get_taq_water_avail(2)

                        req_cols = st.columns([5,2,5,2,5,2,5,2])
                        req_cols[0].write('Required reaction wells')
                        req_cols[1].write(str(reactions))
                        req_cols[2].write('Required PCR plates')
                        req_cols[3].write(str(ceil(reactions/384)))
                        req_cols[4].write('PCR1 required taq volume')
                        req_cols[5].write(str(primer_taq_vol/1000)+' ul')
                        req_cols[6].write('PCR1 available taq volume')
                        req_cols[7].write(str(taq_water_avail_PCR1[0]/1000)+' ul')
                        req_cols[0].write('PCR1 required water volume')
                        req_cols[1].write(str(primer_water_vol/1000)+' ul')
                        req_cols[2].write('PCR1 available water volume')
                        req_cols[3].write(str(taq_water_avail_PCR1[1]/1000)+' ul')
                        req_cols[0].write('PCR2 required taq volume')
                        req_cols[1].write(str(barcode_taq_vol/1000)+' ul')
                        req_cols[2].write('PCR2 available taq volume')
                        req_cols[3].write(str(taq_water_avail_PCR2[0]/1000)+' ul')
                        req_cols[4].write('PCR2 required water volume')
                        req_cols[5].write(str(barcode_water_vol/1000)+' ul')
                        req_cols[6].write('PCR2 available water volume')
                        req_cols[7].write(str(taq_water_avail_PCR2[1]/1000)+' ul')
                        req_cols[0].write('Barcode pairs remaining')
                        req_cols[1].write(str(barcodes_remain))
                        req_cols[2].write('Barcode pairs available')
                        req_cols[3].write(str(barcodes_avail))
                        req_cols[4].write('Barcode vols sufficient')
                        req_cols[5].write(str(barcode_vol_okay))
                    with st.expander('Required/available primer volumes', expanded=False):
                        primer_names = set(primer_vols.keys())
                        print(f"{primer_names=}")
                        assay_columns = st.columns(12)
                        assay_columns[0].write('Primer')
                        assay_columns[1].write('Uses')
                        assay_columns[2].write('Volume ul')
                        assay_columns[3].write('Avail ul')
                        assay_columns[4].write('Primer')
                        assay_columns[5].write('Uses')
                        assay_columns[6].write('Volume ul')
                        assay_columns[7].write('Avail ul')
                        assay_columns[8].write('Primer')
                        assay_columns[9].write('Uses')
                        assay_columns[10].write('Volume ul')
                        assay_columns[11].write('Avail ul')
                        #for k,p in enumerate(sorted(primer_names)):
                        #    r = k%3  # horizontal offset
                        #    assay_columns[r+0].write(p)
                        #    assay_columns[r+1].write(assay_usage[p])
                        #    assay_columns[r+2].write(primer_vols[p])
                        #    assay_columns[r+3].write(primer_avail_vols[p])


            elif st.session_state['view option'] == 2:
                # show plates in table form and let the user edit them to fix minor mistakes
                pass
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
                    st.dataframe(df)

        with st.container():
            with st.expander('Add reference sequences, assay lists, primer plates, barcode plate', expanded=False):
                misc_upload_col1, misc_upload_col2 = st.columns(2)
                with misc_upload_col1:
                    uploaded_reference = st.file_uploader('Upload custom reference file', key='ref_uploader', type=['txt','fa','fasta'])
                    if uploaded_reference:
                        if 'reference_upload' not in st.session_state or st.session_state['reference_upload'] != uploaded_reference.name:
                            success = st.session_state['experiment'].add_reference(uploaded_reference)
                            st.session_state['reference_upload'] = uploaded_reference.name
                            st.experimental_rerun()
                    # display already uploaded reference file names here
                with misc_upload_col2:
                    uploaded_assaylist = st.file_uploader('Upload assay list', key='assaylist_uploader', type=['txt','csv'])
                    if uploaded_assaylist:
                        if 'assaylist_upload' not in st.session_state or st.session_state['assaylist_upload'] != uploaded_assaylist.name:
                            success = st.session_state['experiment'].add_assaylist(uploaded_assaylist)
                            st.session_state['assaylist_upload'] = uploaded_assaylist.name
                            st.experimental_rerun()
                    # display already uploaded assaylist file names here
                misc_upload_col3, misc_upload_col4 = st.columns(2)
                with misc_upload_col3:
                    uploaded_primer_layouts = st.file_uploader('Upload primer plate layouts', 
                            key='primer_uploader', type='csv', accept_multiple_files=True)
                    if uploaded_primer_layouts:
                        if 'primer_layouts_upload' not in st.session_state or st.session_state['primer_layouts_upload'] != uploaded_primer_layouts[0].name:
                            success = st.session_state['experiment'].add_primer_layouts(uploaded_primer_layouts)
                            st.session_state['primer_layouts_upload'] = uploaded_primer_layouts[0].name
                            st.experimental_rerun()
                    # display already uploaded primer plate file names here
                with misc_upload_col4:
                    uploaded_barcode_layouts = st.file_uploader('Upload i7i5 barcode plate layout', 
                            key='barcode_uploader', type='csv', accept_multiple_files=True)
                    if uploaded_barcode_layouts:
                        if 'barcode_layout_upload' not in st.session_state or st.session_state['barcode_layout_upload'] != uploaded_barcode_layouts[0].name:
                            success = st.session_state['experiment'].add_barcode_layouts(uploaded_barcode_layouts)
                            st.session_state['barcode_layout_upload'] = uploaded_barcode_layouts[0].name
                            st.experimental_rerun()
                    # display already uploaded barcode layout file names here


        if st.session_state['pipe_stage'] == 1:
            with st.container():
                _, title_col, _ = st.columns([1,1,1])
                title_col.subheader('Load sample plate data')
                with st.expander('Add data from Rodentity JSON files',expanded=True):
                    # C1: upload, C2: four slots plus 'X' remove, C3: dest plate; "accept" button
                    col_r1, col_r2, col_r3= st.columns(3)
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
                
                    with col_r1:
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
                    with col_r2:
                        st.write('Checkboxes required for clearing plate IDs')
                    
                        for i,p in enumerate(plates_to_clear):
                            plates_to_clear[i] = st.checkbox(f"P{str(i+1)}: {file_io.unguard_pbc(exp.unassigned_plates[i+1], silent=True)}", 
                                    help='Click the checkbox to allow a set plate ID to be cleared', key='check'+str(i+1))
                    with col_r3:
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
                with st.expander('Add data from a custom manifest CSV file', expanded=True):
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

            # Add data to experiment
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
                        
            
        if st.session_state['pipe_stage'] == 2:
            with st.container():
                _, title_col, _ = st.columns([1,5,1])
                title_col.subheader('Nimbus: combine 96-well sample plates to 384-well DNA plate')
                if not st.session_state['experiment'].dest_sample_plates:
                    st.markdown('<p style="color:#FF0000">Load data inputs to enable Nimbus input file generation.</p>', unsafe_allow_html=True)
                else:
                    nim_col1, nim_col2, nim_col3 = st.columns(3)
                    nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                    with nim_col1:
                        # do we have any Nimbus inputs to generate + download
                        yet_to_run = len(st.session_state['experiment'].dest_sample_plates) - len(nfs)
                        if yet_to_run:
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
                            st.markdown('<h3>Download Nimbus input files</h3>', unsafe_allow_html=True)
                            for i,nf in enumerate(nfs):
                                st.download_button(f"{nf}", open(nf), file_name=PurePath(nf).name, 
                                        key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                        
                    with nim_col2:
                        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                        print (f"{nfs=} {efs=} {xbcs=}")
                        # Nimbus upload
                        if not efs and not xbcs:
                            st.markdown('<h3>Awaiting Nimbus output files</h3>', unsafe_allow_html=True)
                        elif efs and not xbcs:
                            st.markdown('<h3>All Nimbus outputs now uploaded</h3>', unsafe_allow_html=True) 
                        elif xbcs:
                            st.markdown('<h3>Upload Nimbus output files</h3>', unsafe_allow_html=True)
                            if 'nim_upload' not in st.session_state:
                                st.session_state['nim_upload'] = ''
                            nim_outputs = st.file_uploader('Nimbus output files "Echo_384_COC_0001_...csv"', type='csv', 
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
                                    
                    with nim_col3:
                        # list nimbus outputs that we still need
                        st.markdown('<h3>Missing Nimbus outputs:</h3>', unsafe_allow_html=True)
                        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                        missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
                        mn_str = '</br>'.join(missing_nims)
                        st.markdown(f"<p>{mn_str}</p>", unsafe_allow_html=True)
                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                if not xbcs:
                    _, msg_col, _ = st.columns([1,1,1])
                    msg_col.markdown('<h3 style="color:#FF0000">Ready to run Echo primer stage</h3>', unsafe_allow_html=True)
                            
                    
                
        if st.session_state['pipe_stage'] == 3:  # echo primers
            with st.container():
                st.subheader('Primer plates')
                with st.expander('Primer plate details'):
                    pass
        if st.session_state['pipe_stage'] == 4:  # echo barcodes
            pass

        if st.session_state['pipe_stage'] == 5:  # Miseq
            pass

        if st.session_state['pipe_stage'] == 6:  # Genotyping
            pass

        if st.session_state['pipe_stage'] == 7:  # Review
            pass

        

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























