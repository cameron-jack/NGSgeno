#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New streamlit interactive web interface
@last_edit: 2022-08-29
@edit_comment:

Display methods for the main GUI pipeline. Methods in include data_table, display_pcr_componenent, display_pcr_componenent 
as well as aggrid_interactive_table and delete_entries

"""

from ctypes.wintypes import WIN32_FIND_DATAA
from msilib.schema import File
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

def aggrid_interactive_table(df: pd.DataFrame, grid_height: int=250, key: int=1):
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
        theme="streamlit",
        #update_mode=GridUpdateMode.MODEL_CHANGED,
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        key=key,
        allow_unsafe_jscode=True,
        rowSelection='multiple',
        selection_mode='multiple',
        rowMultiSelectWithClick=True,
    )
    return selection


def delete_entries(exp, entries):
    """ removes data from an experiment, called by buttons """
    print(f"Debug: (delete_entries) {exp.name=} {entries=}")
    if exp is None:
        if 'delete_entries' in st.session_state:
            st.session_state['delete_entries'] = []
        return
    if entries: # list
        for i,entry in enumerate(entries):
            for key in entry:
                #print(key)
                if 'PID' in key and entry[key]:
                    if not file_io.is_guarded_pbc(entry[key]):
                        entries[i][key] = file_io.guard_pbc(entry[key])
        print(f"delete_entries {entries=}")
        #exp.delete_plates(entries)
        #exp.save()
    print('delete_entries ran')


def data_table(key, view='Summary'):
    """
    Shows the loaded data in a table.
    Args:
        key (Any): key for the table so that it is a unique widget, if there's more than one in a page
        view ('Summary', 'Edit Table', 'Plate View', 'Explore Assays', 'View Log'): Different options for viewing the data

    Need to fix the delete selection part - clicking no should reset the table view. Also deleting the entry doesn't affect the pipeline further down.
    """

    #Need to rearrange
    table_container = st.container()
    log_container = st.container()
    plate_view_container = st.container()
    explore_assay_container = st.container()
    view_log_container = st.container()

    assay_col1, assay_col2, assay_col3, assay_col4, assay_col5 = explore_assay_container.columns(5)

    exp = st.session_state['experiment']

    #Experiment summary - table

    if view == 'Summary':
        selection = []
        st.session_state['nuke'] = ''

        #st.session_state['delete_selection'] = False

        
        if st.session_state['experiment']:
            #st.write(st.session_state['experiment'].name)
            df = st.session_state['experiment'].inputs_as_dataframe()
            if df is None or not isinstance(df, pd.DataFrame):
                st.write('No data loaded')
            else:
                #print(f"{type(df)=}, {df=})")
                selection = aggrid_interactive_table(df, key=key)
            #print(f"{type(selection)=} {selection=}")
            if selection['selected_rows']:
                # only do the code below if this is a fresh selection
                rows = selection["selected_rows"]
                lines = '\n'.join(['DNA PID: '+r['DNA PID'] for r in rows if r['DNA PID'] != 'Total'])
                if 'previous_delete_group_selection' not in st.session_state:
                    st.session_state['previous_delete_group_selection'] = None
                if lines and st.session_state['previous_delete_group_selection'] != lines:
                    table_container.markdown(f"**You selected {lines}**")
                    del_col1, del_col2, del_col3, _ = table_container.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    delete_button = del_col2.button("Yes",on_click=delete_entries, 
                            args=(st.session_state['experiment'], selection['selected_rows']), key="delete " + str(key), help=f"Delete {lines}")
                    mistake_button = del_col3.button("No",on_click=delete_entries, args=(st.session_state['experiment'], []), 
                            key="keep " + str(key), help=f"Keep {lines}")
                    
                    if mistake_button:
                        st.session_state['previous_delete_group_selection'] = lines
                        lines = None
                        selection['selected_rows'] = None
                        st.experimental_rerun()
                    elif delete_button:
                        selection['selected_rows'] = None
                        lines = None
                        st.session_state['previous_delete_group_selection'] = None
                        st.experimental_rerun()

    elif view == 'Edit Table':
        # Show plates in table form and let the user edit them to fix minor mistakes
        pass

    #Show a visual of the plates
    elif view == 'Plate View':
        plate_ids = []
        for pid in exp.plate_location_sample:
            plate_ids.append(pid.strip('p'))
        #plate_view_container.write(pids)
        plate_selectbox = plate_view_container.selectbox('Plate ID to view', plate_ids, key=key)
        #plate_input = plate_view_container.text_input('Plate ID to view', key='plate_input1')
        if plate_selectbox:
            plate_id = file_io.guard_pbc(plate_selectbox)
            if plate_id in st.session_state['experiment'].plate_location_sample:
                heatmap_str = generate_heatmap_html(st.session_state['experiment'], plate_id, scaling=1.2)
                #with open("debug.html", 'wt') as outf:
                #    print(heatmap_str, file=outf)
                components.html(heatmap_str, height=700, scrolling=True)
            else:
                plate_view_container.markdown('<p style="color:#FF0000">Plate barcode not found in experiment</p>', unsafe_allow_html=True)
        # let user choose plates and then display them as a plate layout with well contents on hover
        # plate_col1, plate_col2 = st.columns(2)
        # for plate in st.session_state['experiment'].plate_location_sample:
        #    break
        # heatmap_str = generate_heatmap_html(plate_id, st.session_state['experiment'], scaling=1)
        # components.html(heatmap_str, height=600, scrolling=True)

    elif view == 'Explore Assays':  # Assay explorer
        # display allowed and denied assays, missing assays, used assays, etc.
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
                
    elif view == 'View Log':
        # display the experiment log and let the user filter the view in a number of ways
        df = pd.DataFrame(st.session_state['experiment'].get_log(30))
        if 'view_box_size' in st.session_state:
            height = st.session_state['view_box_size']
        else:
            height = 250
        log_container.dataframe(df, height=height)


def display_pcr_components(assay_usage, PCR_stage='all'):
    """
    Expander widget that shows the required componenents for each PCR reaction, including wells, PCR plates, taq+water plates, index pairs 
    Args:
        assay_usage: from Experiment.get_assay_usage
        PCR_stage(1, 2, 'all'): 1 = Echo Primer stage, 2 = Echo Indexing, all = Both
    
    """
    AFTER_NIM_WELLS = 384
    
    #pcr_components_exp = st.expander('Required Componenents for PCR Reactions', expanded=False)
    pcr_components_container = st.container()
    _,filter_col1,_,filter_col2, _ = pcr_components_container.columns([1,5,1,5,1])
    col_size = [6, 4, 6, 4]
    req_cols = pcr_components_container.columns(col_size)

    if assay_usage:
        primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
                st.session_state['experiment'].get_volumes_required(assay_usage=assay_usage)
        reactions = sum([v for v in assay_usage.values()])
        index_remain, index_avail, index_vol_capacity = st.session_state['experiment'].get_index_remaining_available_volume(assay_usage=assay_usage)
        
    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            0,{},0,0,0,0 
        index_remain, index_avail, index_vol_capacity =\
                st.session_state['experiment'].get_index_remaining_available_volume()

    taq_avail, water_avail, pids = st.session_state['experiment'].get_taqwater_avail()

    if 'assay_filter' not in st.session_state:
        st.session_state['assay_filter'] = 1

        filter_assays = filter_col1.button('Filter assays', key='assay_filter_button', disabled=True)
        unfilter_assays = filter_col2.button('Unfilter assays', key='assay_unfilter_button')

        if unfilter_assays:
            st.session_state['assay_filter'] = 0

        if filter_assays:
            st.session_state['assay_filter'] = 1

    required_wells = ceil(reactions/AFTER_NIM_WELLS)
    user_supplied_PCR = ', '.join([file_io.unguard_pbc(p, silent=True) for p in st.session_state['experiment'].get_pcr_plates()])
    user_supplied_taqwater = ', '.join([file_io.unguard_pbc(p, silent=True) for p in st.session_state['experiment'].get_taqwater_avail()[2]])

    #Wells and plates
    req_cols[0].markdown('<h6 style="color:#B1298D">Required reaction wells</h6>', unsafe_allow_html=True)
    req_cols[1].write(str(reactions))

    req_cols[2].markdown('<h6 style="color:#B1298D">Required PCR plates</h6>', unsafe_allow_html=True)
    req_cols[3].write(str(required_wells))

    req_cols[0].markdown('<h6 style="color:#000000">User supplied PCR plates</h6>', unsafe_allow_html=True)
    req_cols[1].write(user_supplied_PCR)

    req_cols[2].markdown('<h6 style="color:#000000">User supplied taq+water plates</h6>', unsafe_allow_html=True)
    req_cols[3].write(user_supplied_taqwater)

    pcr_components_container.write('')


    #PCR1 and PCR2 taq and water 
    taq_avail_vol = taq_avail/1000
    water_avail_vol = water_avail/1000

    if PCR_stage == 1 or PCR_stage == 'all':
        #PCR1 componenents
        pcr_components_container.markdown('<h6 style="color:#B1298D">PCR 1</h6>', unsafe_allow_html=True)
        if taq_avail < primer_taq_vol or water_avail < primer_water_vol:
            pcr_components_container.markdown('<p style="color:#B92A5D"> Upload more taq and water plates to fulfill PCR 1 requirements </p>', unsafe_allow_html=True)
    
        PCR1_cols = pcr_components_container.columns(col_size)

        PCR1_cols[0].markdown('**Required water volume**', unsafe_allow_html=True)
        PCR1_cols[1].write(str(primer_water_vol/1000)+' μl')

        PCR1_cols[2].markdown('**Available water volume**')
        PCR1_cols[3].write(str(water_avail_vol)+' μl')

        PCR1_cols[0].markdown('**Required taq volume**')
        PCR1_cols[1].write(str(primer_taq_vol/1000)+' μl')

        PCR1_cols[2].markdown('**Available taq volume**')
        PCR1_cols[3].write(str(taq_avail_vol)+' μl')
        pcr_components_container.write('')

    if PCR_stage == 2 or PCR_stage == 'all':
        #PCR 2 components
        pcr_components_container.markdown('<h6 style="color:#B1298D">PCR 2</h6>', unsafe_allow_html=True)
    
        if taq_avail < (primer_taq_vol + index_taq_vol) or water_avail < (primer_water_vol + index_water_vol):
            pcr_components_container.markdown('<p style="color:#B92A5D"> Upload more taq and water plates to fulfill PCR 2 requirements </p>', unsafe_allow_html=True)
    
        PCR2_cols = pcr_components_container.columns(col_size)

        PCR2_cols[0].markdown('**Required water volume**')
        PCR2_cols[1].write(str(index_water_vol/1000)+' μl')

        PCR2_cols[2].markdown('**Available water volume**')
        PCR2_cols[3].write(str(water_avail_vol)+' μl')

        PCR2_cols[0].markdown('**Required taq volume**')
        PCR2_cols[1].write(str(index_taq_vol/1000)+' μl')

        PCR2_cols[2].markdown('**Available taq volume**')
        PCR2_cols[3].write(str(taq_avail_vol)+' μl')
    
        #Index
        index_cols = pcr_components_container.columns(col_size)

        if index_vol_capacity > (index_avail - index_remain):
            index_pairs_allowed = '<p style="color:green">'+str(index_vol_capacity)+'</p>'
        else:
            index_pairs_allowed = '<p style="color:red">'+str(index_vol_capacity)+'</p>'

        index_cols[0].markdown('**Index Pairs Remaining**')
        index_cols[1].write(str(index_remain))
        index_cols[2].markdown('**Max Possible Index Pairs**')
        index_cols[3].write(str(index_avail))
        index_cols[0].markdown('**Index Pairs Allowed by Volume**')
        index_cols[1].markdown(index_pairs_allowed, unsafe_allow_html=True)


def display_primer_components(assay_usage, show_primers=True):

    primer_components_exp = st.expander('Primer Wells, Uses, and Volumes', expanded=False)
    ptab1, ptab2, ptab3, ptab4 = primer_components_exp.tabs(["Wells", "Uses", "Volume μL", "Available μL"])
    
    #Set up columns
    col_size = 3
    num_cols = 6
    columns = []
    i = 0
    while i < num_cols:
        columns.append(col_size)
        i+=1

    wells_col = ptab1.columns(columns)
    uses_col = ptab2.columns(columns)
    volume_col = ptab3.columns(columns)
    avail_col = ptab4.columns(columns)

    #Values
    if assay_usage:
        primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
                st.session_state['experiment'].get_volumes_required(assay_usage=assay_usage)
    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            0,{},0,0,0,0 

    primer_avail_counts, primer_avail_vols = st.session_state['experiment'].get_primers_avail()
    per_use_vol = CAP_VOLS['384PP_AQ_BP'] - DEAD_VOLS['384PP_AQ_BP']
    primer_names = set(primer_vols.keys())

    for i in range(6):
        if (i+2) % 2 == 0:
            wells_col[i].markdown('**Primer**')
            uses_col[i].markdown('**Primer**')
            volume_col[i].markdown('**Primer**')
            avail_col[i].markdown('**Primer**')
        else:
            wells_col[i].markdown('**Wells**')
            uses_col[i].markdown('**Uses**')
            volume_col[i].markdown('**Volume μL**')
            avail_col[i].markdown('**Available μL**')

    for k,p in enumerate(sorted(primer_names)):
        if p != '':
            r = (k*2)%6  # horizontal offset
            need_vol = primer_vols.get(p,0)
            num_wells = ceil((need_vol)/per_use_vol)

            wells_col[r+0].write(p)
            wells_col[r+1].write(num_wells)

            uses_col[r+0].write(p)
            uses_col[r+1].write(assay_usage.get(p,0))

            volume_col[r+0].write(p)
            volume_col[r+1].write(primer_vols.get(p,0)/1000)

            avail_col[r+0].write(p)
            avail_col[r+1].write(primer_avail_vols.get(p,0)/1000)
