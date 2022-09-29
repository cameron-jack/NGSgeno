#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New streamlit interactive web interface
@last_edit: 2022-08-29
@edit_comment:

Display methods for the main GUI pipeline. Methods in include data_table, display_pcr_componenent,
display_pcr_componenent as well as aggrid_interactive_table and delete_entries
"""

import os
import sys
from pathlib import PurePath, Path
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO
import inspect

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
    #print('aggrid1', st.session_state['experiment'].name)
    
    options.configure_side_bar()
    #st.write(df)

    options.configure_selection("multiple", use_checkbox=False, \
                rowMultiSelectWithClick=False, suppressRowDeselection=False)
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        height=grid_height,
        gridOptions=options.build(),
        # _available_themes = ['streamlit','light','dark', 'blue', 'fresh','material']
        theme="alpine",
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


def data_table(key, options=False, table_option=None):
    """
    Shows the loaded data in a table.
    Args:
        key (Any): key for the table so that it is a unique widget, if there's more than one in a page
        options (bool): Shows different options for table view in a selectbox if True
        table_option (str): If options is False, can choose a table_option instead ('Summary', 'Edit Table', 
        'Plate View', 'Explore Assays', 'View Log)

    Need to fix the delete selection part - clicking no should reset the table view. 
    Also deleting the entry doesn't affect the pipeline further down.
    """

    data_container = st.container()

    if options==True:
        table_option = data_container.selectbox('View Data options',
                                ('Summary', 'Edit Table', 'Plate View', 'Explore Assays', 'View Log'),\
                                             key='select_' + str(key))
    elif options==False and table_option is None:
        table_option = 'Summary'

    data_table_title = data_container.title('')
    data_table_title.markdown(f'<h5 style="text-align:center;color:#7eaecc">{table_option}</h5>',\
        unsafe_allow_html=True)

    #exp = st.session_state['experiment']

    #Experiment summary - table

    if table_option == 'Summary':
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
            # if selection['selected_rows']:
            #     # only do the code below if this is a fresh selection
            #     rows = selection["selected_rows"]
            #     lines = '\n'.join(['DNA PID: '+r['DNA PID'] for r in rows if r['DNA PID'] != 'Total'])
            #     if 'previous_delete_group_selection' not in st.session_state:
            #         st.session_state['previous_delete_group_selection'] = None
            #     if lines and st.session_state['previous_delete_group_selection'] != lines:
            #         data_container.markdown(f"**You selected {lines}**")
            #         del_col1, del_col2, del_col3, _ = data_container.columns([2,1,1,4])
            #         del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
            #         delete_button = del_col2.button("Yes",on_click=delete_entries, 
            #                 args=(st.session_state['experiment'], selection['selected_rows']), key="delete " + str(key), help=f"Delete {lines}")
            #         mistake_button = del_col3.button("No",on_click=delete_entries, args=(st.session_state['experiment'], []), 
            #                 key="keep " + str(key), help=f"Keep {lines}")
                    
            #         if mistake_button:
            #             st.session_state['previous_delete_group_selection'] = lines
            #             lines = None
            #             selection['selected_rows'] = None
            #             st.experimental_rerun()
            #         elif delete_button:
            #             selection['selected_rows'] = None
            #             lines = None
            #             st.session_state['previous_delete_group_selection'] = None
            #             st.experimental_rerun()

    elif table_option == 'Edit Table':
        # Show plates in table form and let the user edit them to fix minor mistakes
        pass

    #Show a visual of the plates
    elif table_option == 'Plate View':
        plate_ids = []
        for pid in st.session_state['experiment'].plate_location_sample:
            plate_ids.append(file_io.unguard(pid))

        #Let user choose which plate to view
        plate_selectbox = data_container.selectbox('Plate ID to view', plate_ids, key=key)
        #plate_input = data_container.text_input('Plate ID to view', key='plate_input1')
        if plate_selectbox:
            plate_id = file_io.guard_pbc(plate_selectbox)
            if plate_id in st.session_state['experiment'].plate_location_sample:
                heatmap_str = generate_heatmap_html(st.session_state['experiment'], plate_id, scaling=1.2)
                #with open("debug.html", 'wt') as outf:
                #    print(heatmap_str, file=outf)
                components.html(heatmap_str, height=700, scrolling=True)
            else:
                plate_barcode_error_msg = "Plate barcode not found in experiment"
        # let user choose plates and then display them as a plate layout with well contents on hover
        # plate_col1, plate_col2 = st.columns(2)
        # for plate in st.session_state['experiment'].plate_location_sample:
        #    break
        # heatmap_str = generate_heatmap_html(plate_id, st.session_state['experiment'], scaling=1)
        # components.html(heatmap_str, height=600, scrolling=True)
        plate_barcode_error = data_container.markdown(f'<p style="color:#FF0000">{plate_barcode_error_msg}</p>',\
                unsafe_allow_html=True)

    elif table_option == 'Explore Assays':
        assay_col1, assay_col2, assay_col3, assay_col4, assay_col5 = data_container.columns(5)
          # Assay explorer
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
                
    elif table_option == 'View Log':
        func = sys._getframe(1).f_code.co_name
        func_line = inspect.getframeinfo(sys._getframe(1)).lineno
        print("Function", func, type(func_line))
        # display the experiment log and let the user filter the view in a number of ways
        df = pd.DataFrame(st.session_state['experiment'].get_log(30),\
                     columns=st.session_state['experiment'].get_log_header())
        print(df)
        if 'view_box_size' in st.session_state:
            height = st.session_state['view_box_size']
        else:
            height = 250
        data_container.dataframe(df, height=height)


def display_pcr_components(assay_usage, PCR_stage=1, show_general=True):
    """
    Expander widget that shows the required componenents for each PCR reaction, 
    including wells, PCR plates, taq+water plates, index pairs 
    Args:
        assay_usage: from Experiment.get_assay_usage
        PCR_stage (1, 2): 1 = Echo Primer stage, 2 = Echo Indexing
        show_general (bool): if True show all of the general information for PCR reactions,
        including user supplied plates.  
    
    """
    AFTER_NIM_WELLS = 384
    pcr_comps_area = st.container()
    
    _,filter_col1,_,filter_col2, _ = pcr_comps_area.columns([1,5,1,5,1])
    col_size = [6, 4, 6, 4]
    
    if assay_usage:
        primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
                st.session_state['experiment'].get_volumes_required(assay_usage=assay_usage)
        reactions = sum([v for v in assay_usage.values()])
        index_remain, index_avail, index_vol_capacity =\
                     st.session_state['experiment'].get_index_remaining_available_volume(assay_usage=assay_usage)
        
    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            0,{},0,0,0,0 
        index_remain, index_avail, index_vol_capacity =\
                st.session_state['experiment'].get_index_remaining_available_volume()

    taq_avail, water_avail, pids = st.session_state['experiment'].get_taqwater_avail()
    
    if show_general:
        if 'assay_filter' not in st.session_state:
            st.session_state['assay_filter'] = 1

            filter_assays = filter_col1.button('Filter assays', key='assay_filter_button', disabled=True)
            unfilter_assays = filter_col2.button('Unfilter assays', key='assay_unfilter_button')

            if unfilter_assays:
                st.session_state['assay_filter'] = 0

            if filter_assays:
                st.session_state['assay_filter'] = 1

        required_wells = ceil(reactions/AFTER_NIM_WELLS)
        user_supplied_PCR = ', '.join([file_io.unguard_pbc(p, silent=True)\
                    for p in st.session_state['experiment'].get_pcr_plates()])
        user_supplied_taqwater = ', '.join([file_io.unguard_pbc(p, silent=True)\
                    for p in st.session_state['experiment'].get_taqwater_avail()[2]])

        comp_warning_area = pcr_comps_area.empty()
        comp_warning = ''

        req_cols = pcr_comps_area.columns(col_size)
        #Wells and plates
        req_cols[0].markdown('**Number of required PCR plates**', unsafe_allow_html=True)
        req_cols[1].write(str(required_wells))
        req_cols[0].markdown('**Number of required reaction wells**', unsafe_allow_html=True)
        req_cols[1].write(str(reactions))
        req_cols[2].markdown('**User supplied PCR plates**', unsafe_allow_html=True)
        if user_supplied_PCR:
            req_cols[3].write(user_supplied_PCR)
        else:
            req_cols[3].write('None')
        req_cols[2].markdown('**User supplied taq/water plates**', unsafe_allow_html=True)
        if user_supplied_PCR:
            req_cols[3].write(user_supplied_taqwater)
        else:
            req_cols[3].write('None')
        pcr_comps_area.write('')

    #PCR1 and PCR2 taq and water 
    taq_avail_vol = taq_avail/1000
    water_avail_vol = water_avail/1000

    pcr_cols = pcr_comps_area.columns(col_size)

    if PCR_stage == 1:
        if taq_avail < primer_taq_vol or water_avail < primer_water_vol:
            comp_warning = 'Provide more taq and water plates to fulfill PCR 1 requirements'

        required_water_vol_str = str(primer_water_vol/1000) + ' μl'
        water_avail_vol_str = str(water_avail_vol)+' μl'
        required_taq_vol_str = str(primer_taq_vol/1000)+' μl'
        avail_taq_vol_str = str(taq_avail_vol)+' μl'

    if PCR_stage == 2:
        if taq_avail < (primer_taq_vol + index_taq_vol) or water_avail < (primer_water_vol + index_water_vol):
            comp_warning = "Upload more taq and water plates to fulfill PCR 2 requirements"
        
        required_water_vol_str = str(index_water_vol/1000)+ ' μl'
        water_avail_vol_str = str(water_avail_vol) + ' μl'
        required_taq_vol_str = str(index_taq_vol/1000)+ ' μl'
        avail_taq_vol_str = str(taq_avail_vol)+ ' μl'

        index_cols = pcr_comps_area.columns(col_size)

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

    pcr_cols[0].markdown('**Required water volume**')
    pcr_cols[1].write(required_water_vol_str)
    pcr_cols[2].markdown('**Available water volume**')
    pcr_cols[3].write(water_avail_vol_str)
    pcr_cols[0].markdown('**Required taq volume**')
    pcr_cols[1].markdown(required_taq_vol_str)
    pcr_cols[2].markdown('**Available taq volume**')
    pcr_cols[3].write(avail_taq_vol_str)

    comp_warning_area.markdown(f'<p style="color:#f63366">{comp_warning}</p>', unsafe_allow_html=True)


def display_primer_components(assay_usage, expander=True):
    if expander:
        primer_components_exp = st.expander('Primer Wells, Uses, and Volumes', expanded=False)
        ptab1, ptab2, ptab3, ptab4 = primer_components_exp.tabs(["Wells", "Uses", "Volume μL", "Available μL"])
    else:
        ptab1, ptab2, ptab3, ptab4 = st.tabs(["Wells", "Uses", "Volume μL", "Available μL"])

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

def st_directory_picker(label='Selected directory:', initial_path=Path(),\
            searched_file_types=['fastq','fastq.gz','fq','fq.gz']):
    """
    Streamlit being JS/AJAX has no ability to select a directory. This is for server paths only.
    Initial code by Aidin Jungo here: https://github.com/aidanjungo/StreamlitDirectoryPicker
    Cannot be used inside an st.column as it produces its own columns.
    """
    if "path" not in st.session_state:
        st.session_state['path'] = initial_path.absolute()

    st.text_input(label, st.session_state['path'])
    subdirectories = [f.stem for f in st.session_state['path'].iterdir() if f.is_dir() \
            and (not f.stem.startswith(".") and not f.stem.startswith("__"))]
    contains_files = [f.name for f in st.session_state['path'].iterdir() if f.is_file() \
            and any([f.name.endswith(sft) for sft in searched_file_types])]

    col1, col2, col3, _ = st.columns([1, 3, 1, 5])
    
    with col1:
        st.markdown("Back")
        if st.button("←") and "path" in st.session_state:
            st.session_state['path'] = st.session_state['path'].parent
            st.experimental_rerun()

    with col2:   
        if subdirectories:
            st.session_state['new_dir'] = st.selectbox("Subdirectories", sorted(subdirectories))
        else:
            st.markdown("#")
            st.markdown("<font color='#FF0000'>No subdir</font>", unsafe_allow_html=True)

    with col3:
        if subdirectories:
            st.markdown("Select")
            if st.button("→") and "path" in st.session_state:
                st.session_state['path'] = Path(st.session_state['path'], st.session_state['new_dir'])
                st.experimental_rerun()
    st.markdown(\
            f'<h5 style="color:#000000">Contains {len(contains_files)} files of type {",".join(searched_file_types)}</h5>',\
                    unsafe_allow_html=True)

    return st.session_state['path']