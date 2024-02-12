#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University

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
from time import sleep

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components
import extra_streamlit_components as stx
from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode

from stutil import custom_text
try:
    from bin.experiment import Experiment, EXP_FN, load_experiment
except ModuleNotFoundError:
    from experiment import Experiment, EXP_FN, load_experiment
try:
    import bin.util as util
except ModuleNotFoundError:
    import util                    
try:
    import bin.db_io as db_io
except ModuleNotFoundError:
    import db_io
try:
    from bin.makehtml import generate_heatmap_html
except ModuleNotFoundError:
    from makehtml import generate_heatmap_html

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
    selection = None
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        height=grid_height,
        gridOptions=options.build(),
        # _available_themes = ['streamlit','light','dark', 'blue', 'fresh','material']
        #theme="alpine",
        theme='streamlit',
        update_mode=GridUpdateMode.MODEL_CHANGED | GridUpdateMode.SELECTION_CHANGED | GridUpdateMode.VALUE_CHANGED,
        key=key,
        reload_data=True,
        allow_unsafe_jscode=True,
        rowSelection='multiple',
        selection_mode='multiple',
        rowMultiSelectWithClick=True,
    )
    return selection


def manage_delete(delbox, category, ids):
    """
    Callback for deletion operations
    """
    exp = st.session_state['experiment']
    if category == 'file':
        for id in ids:
            if id not in exp.uploaded_files:
                continue
            success = exp.del_file_record(id)     
            if success:
                #delbox.write(f'{id} removed')
                st.write(f'{id} removed')
            else:
                #delbox.write(f'{id} could not be removed')
                st.write(f'{id} could not be removed')
            st.session_state['previous_file_delete_selection'] = ids
    elif category == 'plate':
        gids = [util.guard_pbc(pid, silent=True) for pid in ids]
        gids = [gid for gid in gids if gid in exp.plate_location_sample]
        success = exp.delete_plates(gids)
        if success:
            #delbox.write(f'{ids} removed')
            st.write(f'{ids} removed')
        else:
            #delbox.write(f'{ids} could not be removed')
            st.write(f'{ids} could not be removed')
        st.session_state['previous_plate_delete_selection'] = ids
    elif category == 'group':  # from summary
        for row in ids:
            dest_pid = util.guard_pbc(row['DNA PID'], silent=True)
            if dest_pid not in exp.dest_sample_plates:
                exp.log(f"Error: {row['DNA PID']=} doesn't actually exist in the experiment!")
                continue
            sample_pids = exp.dest_sample_plates[dest_pid]
            delete_pids = sample_pids + [dest_pid]
            exp.delete_plates(delete_pids)
        st.session_state['previous_group_delete_selection'] = ids
    exp.save()
    return True


def cancel_delete(category, ids):
    """
    Callback for delete operations
    """
    if category == 'plate':
        st.session_state['previous_plate_delete_selection'] = ids
    elif category == 'file':
        st.session_state['previous_file_delete_selection'] = ids
    elif category == 'group':  # from summary
        st.session_state['previous_group_delete_selection'] = ids
    return True


def display_samples(key, height=250):
    """ display a summary of all loaded DNA, amplicon, and sample plates """
    exp = st.session_state['experiment']
    selection = []
    df = exp.inputs_as_dataframe()
    if df is None or not isinstance(df, pd.DataFrame):
        st.write('No 384-well DNA plate data loaded')
    else:    
        selection = aggrid_interactive_table(df, key=key, grid_height=height)
        if 'selected_rows' in selection and selection['selected_rows']:
            rows = selection["selected_rows"]
            rows = [r for r in rows if 'DNA PID' in r and r['DNA PID'] != 'Total']
            if rows:
                if 'previous_group_delete_selection' not in st.session_state:
                    st.session_state['previous_group_delete_selection'] = None
                if st.session_state['previous_group_delete_selection'] != rows:
                    # only do the code below if this is a fresh selection
                    lines = '\n'.join(['DNA PID: '+r['DNA PID'] for r in rows if r['DNA PID'] != 'Total'])
                    st.markdown(f"**You selected {lines}**")
                    delbox = st.container()
                    del_col1, del_col2, del_col3, _ = st.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    del_col2.button("Yes",on_click=manage_delete, 
                            args=(delbox, 'group',rows), key="delete " + str(key), help=f"Delete {lines}")
                    del_col3.button("No",on_click=cancel_delete, args=('group',rows), 
                            key="keep " + str(key), help=f"Keep {lines}")
        selection = None


def display_consumables(key, height=300):
    """
    summarise_consumables():
        d = {'taqwater_pids_pcr1':[], 'taqwater_pids_pcr2':[], 'taq_vol_pcr1':0, 'taq_vol_pcr2':0,'water_vol_pcr1':0, 
                'water_vol_pcr2':0, 'primer_pids':[], 'primer_count_ngs':0, 'primer_count_custom':0, 'unique_primers':set(), 
                'primer_well_count':0, 'assay_primer_mappings':0, 'reference_files':[], 'unique_references':set(), 
                'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
    """
    exp = st.session_state['experiment']
    # display all the required files whether they are, or are not present
    # Plates: primer, index, taq/water, PCR, references, primer/assaylist, 
    headers = ['Purpose','Barcode/ID','Wells','Type','Entries']
    consumables = exp.summarise_consumables()

    data_rows = []
    for tp in consumables['taqwater_pids_pcr1']:
        data_rows.append(['taq/water (PCR 1)', util.unguard_pbc(tp, silent=True), 6, util.PLATE_TYPES['Echo6'], 6])
    for tp in consumables['taqwater_pids_pcr2']:
        data_rows.append(['taq/water (PCR 2)', util.unguard_pbc(tp, silent=True), 6, util.PLATE_TYPES['Echo6'], 6])
    for pp in consumables['primer_pids']:
        data_rows.append(['primer', util.unguard_pbc(pp, silent=True), 384, util.PLATE_TYPES['Echo384'], 
                len(exp.plate_location_sample[pp]['wells'])])
    for ip in consumables['index_pids']:
        data_rows.append(['index', util.unguard_pbc(ip, silent=True), 384, util.PLATE_TYPES['Echo384'],
                len(exp.plate_location_sample[ip]['wells'])])
    for f in consumables['reference_files']:
        data_rows.append(['reference sequences', f, 0, 'File', len(exp.reference_sequences[f])])
    for f in exp.uploaded_files:
        if exp.uploaded_files[f].get('purpose','') == 'assay_primer_map':
            data_rows.append(['assay-primer mappings', f, 0, 'File', consumables['assay_primer_mappings']])
    plate_df = pd.DataFrame(data_rows, columns=headers)
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        st.write('No plates loaded')
    else:
        if 'view_box_size' in st.session_state:
            height = st.session_state['view_box_size']
        else:
            height = height
        selection = aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'consumables_aggrid')


def display_plates(key, plate_usage, height=300): 
    exp = st.session_state['experiment']
    plate_df = pd.DataFrame(plate_usage, columns=['Plates', 'Num Wells', 'Purpose'])
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        st.write('No plates loaded')
    else:
        if 'view_box_size' in st.session_state:
            height = st.session_state['view_box_size']
        else:
            height = height

        selection = aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'plate_aggrid')
        #print(f'{selection=} {selection["selected_rows"]=} {[row["Plates"] for row in selection["selected_rows"]]=}')
        if selection and 'selected_rows' in selection:
            pids = [row['Plates'] for row in selection['selected_rows']]
            gids = [util.guard_pbc(pid, silent=True) for pid in pids]
            gids = [gid for gid in gids if gid in exp.plate_location_sample]
            if gids:
                #print(f'{pids=}')
                if 'previous_plate_delete_selection' not in st.session_state:
                    st.session_state['previous_plate_delete_selection'] = None
                if st.session_state['previous_plate_delete_selection'] == pids:
                    st.session_state['previous_plate_delete_selection'] = None
                else:
                    if pids != st.session_state['previous_plate_delete_selection']:
                        st.markdown(f"**You selected {pids}**")
                    delbox = st.container() # doesn't work reliably in st1.26
                    del_col1, del_col2, del_col3, _ = delbox.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    del_col2.button("Yes",on_click=manage_delete,
                            args=(delbox,'plate',pids), key="delete " + str(key), help=f"Delete {pids}")
                    del_col3.button("No", on_click=cancel_delete,
                            args=('plate',pids), key="keep " + str(key), help=f"Keep {pids}")
    
def display_pcr_components(dna_pids=None, pcr_pids=None, amplicon_pids=None):
    """
    Expander widget that shows the PCR plate information for both PCR reactions
    Args:
        dna_pids: DNA plate IDs included in the reaction
    """
    exp = st.session_state['experiment']
    DNA_PLATE_WELLS = 384
    
    #put as functions in exp
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids=dna_pids)
    num_reactions = sum([primer_usage[p] for p in primer_usage])

    #num reactinos = wells in amplicon
    #then required = num reactions x pcr
    #PCR
    required_pcr_plates = ceil(num_reactions/DNA_PLATE_WELLS)
    
    user_supplied_pcr = ', '.join([util.unguard_pbc(p, silent=True) for p in pcr_pids]) if pcr_pids else ''
    user_supplied_amplicon = ', '.join([util.unguard_pbc(a, silent=True) for a in amplicon_pids]) if amplicon_pids else ''

    if not (pcr_pids or amplicon_pids):
        user_supplied_pcr = '<p style="color:#FF0000">None</p>'

    num_supplied_pcr = len(exp.get_pcr_pids())

    req_PCR_text = '**Number of required PCR plates**'
    req_PCR_num = str(required_pcr_plates)
    if num_supplied_pcr < required_pcr_plates:
        req_PCR_text = '<p style="color:#FF0000"><b>Number of required PCR plates</b></p>'
        req_PCR_num = f'<p style="color:#FF0000">{str(required_pcr_plates)}</p>'
        
    #Page set up
    pcr_comps_area = st.container()
    pcr_cols = pcr_comps_area.columns([6, 4, 6, 4])
    pcr_cols[0].markdown(req_PCR_text, unsafe_allow_html=True)
    pcr_cols[1].markdown(req_PCR_num, unsafe_allow_html=True)

    pcr_cols[0].markdown('**Worst case required reaction wells**')
    pcr_cols[1].write(str(num_reactions), unsafe_allow_html=True)

    pcr_cols[2].markdown('**User supplied PCR plates**')
    pcr_cols[3].markdown(user_supplied_pcr, unsafe_allow_html=True)

    if amplicon_pids:
        pcr_cols[2].markdown('**User supplied amplicon plates**')
        if not pcr_pids:
            pcr_cols[3].write('')
            pcr_cols[3].write('')
        pcr_cols[3].markdown(user_supplied_amplicon, unsafe_allow_html=True)
    
    #Adding white space
    i=0
    while i < 3:
        pcr_cols[2].write('')
        pcr_cols[3].write('')
        i+=1

    for i in range(4):
        pcr_cols[i].write('')


def display_pcr1_components(dna_pids=None):
    """
    Expander widget that shows the required componenents for each PCR reaction, 
    including wells, PCR plates, taq+water plates, index pairs 
    Args:
        dna_pids: DNA plate IDs provided by user
    """
    exp = st.session_state['experiment']
    ul_conv = 1000
    pcr_stage = 1
    
    #Page set up
    pcr_comps_area = st.container()
    col_size = [6, 4, 6, 4]
    req_cols = pcr_comps_area.columns(col_size)
    pcr_cols = pcr_comps_area.columns(col_size)

    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids=dna_pids)
    num_reactions = sum([primer_usage[p] for p in primer_usage])
    
    #Required taq/water
    primer_taq_vol, primer_water_vol = exp.get_taqwater_req_vols_primers(num_reactions)

    taq_avail, water_avail, pids = exp.get_taqwater_avail(pcr_stage=pcr_stage)
    taq_avail_vol = taq_avail/ul_conv
    water_avail_vol = water_avail/ul_conv            

    #get actual values for volume of taq water plates
    required_water_vol_str = str(primer_water_vol/ul_conv) + ' μl'
    water_avail_vol_str = str(water_avail_vol)+' μl'
    required_taq_vol_str = str(primer_taq_vol/ul_conv) + ' μl'
    avail_taq_vol_str = str(taq_avail_vol)+' μl'

    num_req_taq_water_plates = util.num_req_taq_water_plates(primer_taq_vol, primer_water_vol)

    user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True) \
                            for p in exp.get_taqwater_avail(pcr_stage=pcr_stage)[2]])
    
    num_supplied_taqwater = len(user_supplied_taqwater)
    user_taqwater_text = user_supplied_taqwater if user_supplied_taqwater else '<p style="color:#FF0000">None</p>'
    
    if num_supplied_taqwater < num_req_taq_water_plates:
        req_taqwater_text = '<p style="color:#FF0000"><b>Number of required taq/water plates</b></p>'
        req_taqwater_num = f'<p style="color:#FF0000">{str(num_req_taq_water_plates)}</p>'
    else:
        req_taqwater_text = '**Number of required taq/water plates**'
        req_taqwater_num = str(num_req_taq_water_plates)

    for i in range(4):
        req_cols[i].write('')

    req_cols[0].markdown(req_taqwater_text, unsafe_allow_html=True)
    req_cols[1].markdown(req_taqwater_num, unsafe_allow_html=True)

    req_cols[2].markdown(f'**User supplied taq/water plates (PCR {pcr_stage})**', unsafe_allow_html=True)
    req_cols[3].write(user_taqwater_text, unsafe_allow_html=True)

    pcr_cols[0].markdown('**Required water volume**')
    pcr_cols[1].write(required_water_vol_str, unsafe_allow_html=True)
    pcr_cols[2].markdown('**Available water volume**')
    pcr_cols[3].write(water_avail_vol_str, unsafe_allow_html=True)
    pcr_cols[0].markdown('**Required taq volume**')
    pcr_cols[1].markdown(required_taq_vol_str, unsafe_allow_html=True)
    pcr_cols[2].markdown('**Available taq volume**')
    pcr_cols[3].write(avail_taq_vol_str, unsafe_allow_html=True)

def display_pcr2_components(dna_pids=None, pcr_pids=None, amplicon_pids=None):
    """
    Expander widget that shows the required componenents for PCR 2 reaction (index).
    Args:
        pcr_stage (1, 2): 1 = Echo Primer stage, 2 = Echo Indexing
        dna_pids:
        pcr_pids: 
        amplicon_pids:
    """
    #Need to add info about taq water
    exp = st.session_state['experiment']
    ul_conv = 1000
    pcr_stage = 2

    #Index
    fwd_idx, rev_idx = exp.get_index_avail()
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids=dna_pids)
    #num reactions:
    num_reactions = exp.get_num_reactions(primer_usage, pcr_pids, amplicon_pids)
    print(f'{num_reactions=}')

    # num_reactions = sum([primer_usage[p] for p in primer_usage])
    #index_fwd, index_rev = exp.get_index_avail()
    index_max = len(exp.get_index_pairs_avail())
    index_remain = index_max - num_reactions
    #index_remain, index_max = exp.get_index_reactions(primer_usage, fwd_idx, rev_idx)

    #Taq/water (based on pcr stage)
    index_taq_vol, index_water_vol = exp.get_taqwater_req_vols_index(num_reactions)
    
    user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True)\
                                        for p in exp.get_taqwater_avail(pcr_stage=pcr_stage)[2]])
    
    num_supplied_taqwater = len(user_supplied_taqwater)
    
    taq_avail, water_avail, pids = exp.get_taqwater_avail(pcr_stage=pcr_stage)
    taq_avail_vol = taq_avail/ul_conv
    water_avail_vol = water_avail/ul_conv  
    required_water_vol_str = str(index_water_vol/ul_conv)+ ' μl'
    water_avail_vol_str = str(water_avail_vol) + ' μl'
    required_taq_vol_str = str(index_taq_vol/ul_conv)+ ' μl'
    avail_taq_vol_str = str(taq_avail_vol)+ ' μl'

    num_req_taq_water_plates = util.num_req_taq_water_plates(index_taq_vol, index_water_vol)
    
    user_taqwater_text = user_supplied_taqwater
    if not user_supplied_taqwater:
        user_taqwater_text = '<p style="color:#FF0000">None</p>'

    if num_supplied_taqwater < num_req_taq_water_plates:
        req_taqwater_text = f'<p style="color:#FF0000"><b>Number of required taq/water plates (PCR {pcr_stage})</b></p>'
        req_taqwater_num = f'<p style="color:#FF0000">{str(num_req_taq_water_plates)}</p>'
    else:
        req_taqwater_text = f'**Number of required taq/water plates (PCR {pcr_stage})**'
        req_taqwater_num = str(num_req_taq_water_plates)

    if index_remain >= 0:
        index_pairs_remain = '<p style="color:green">'+str(index_remain)+'</p>'     
    else:
        index_pairs_remain = '<p style="color:red">'+str(index_remain)+'</p>'

    #Page set up
    pcr_comps_area = st.container()
    col_size = [6, 4, 6, 4]
    pcr2_col = pcr_comps_area.columns(col_size)

    pcr2_col[0].markdown(req_taqwater_text, unsafe_allow_html=True)
    pcr2_col[1].markdown(req_taqwater_num, unsafe_allow_html=True)
    pcr2_col[2].markdown(f'**Number of user supplied taq/water plates (PCR {pcr_stage})**',unsafe_allow_html=True)
    pcr2_col[3].markdown(user_taqwater_text, unsafe_allow_html=True)

    pcr2_col[0].markdown('**Required water volume**')
    pcr2_col[1].write(required_water_vol_str)
    pcr2_col[2].markdown('**Available water volume**')
    pcr2_col[3].write(water_avail_vol_str)
    pcr2_col[0].markdown('**Required taq volume**')
    pcr2_col[1].markdown(required_taq_vol_str)
    pcr2_col[2].markdown('**Available taq volume**')
    pcr2_col[3].write(avail_taq_vol_str)

    pcr2_col[0].markdown('**Index Pairs Possible**')
    pcr2_col[1].write(str(index_max))
    pcr2_col[0].markdown('**Index Pairs Remaining**')
    pcr2_col[1].markdown(index_pairs_remain, unsafe_allow_html=True)
    



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

def handle_picklist_download(picklist_type, picklist_paths, error_msgs, file_col, btn_col):
    if not picklist_paths:
        error_msgs.append(f"No {picklist_type} picklist available")
    else:
        for ppp in picklist_paths:
            ppp_fn = Path(ppp).name
            with file_col:
                custom_text('p', '#4b778c', ppp_fn, 'right', padding='5px')
            with btn_col:
                st.download_button(label="Download", 
                                   data=open(ppp, 'rt'), 
                                   file_name=ppp_fn, 
                                   mime='text/csv', 
                                   key=f'{picklist_type}_download_'+ppp_fn)


def get_echo1_download_btns():
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)
    error_msgs = []

    dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR1_picklist_filepaths()
    
    picklist_dict = {'DNA': dna_picklist_paths, 'primer': primer_picklist_paths, 'taq/water': taqwater_picklist_paths}

    for pltype, plpath in picklist_dict.items():
        handle_picklist_download(pltype, plpath, error_msgs, picklist_file_col, picklist_btn_col)

    if error_msgs:
        for msg in error_msgs:
            custom_text('p', '#ff0000', msg, 'center', padding='5px')

def get_echo2_download_btns():
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)    
    error_msgs = []

    index_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR2_picklist_filepaths()

    picklist_dict = {'index': index_picklist_paths, 'taq/water': taqwater_picklist_paths}

    for pltype, plpath in picklist_dict.items():
        handle_picklist_download(pltype, plpath, error_msgs, picklist_file_col, picklist_btn_col)
    
    if error_msgs:
        for msg in error_msgs:
            custom_text('p', '#ff0000', msg, 'center', padding='5px')


def show_echo1_outputs():
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)
    dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR1_picklist_filepaths()
    error_msgs = []
    
    if not dna_picklist_paths:
        error_msgs.append('No DNA picklist available')
    else:
        for dpp in dna_picklist_paths:
            dpp_fn = Path(dpp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{dpp_fn}</p>',\
                    unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(dpp, 'rt'), file_name=dpp_fn, mime='text/csv', key='dna_download_'+dpp_fn)

    if not primer_picklist_paths:
        error_msgs.append('No primer picklist available')
    else:
        for ppp in primer_picklist_paths:
            ppp_fn = Path(ppp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ppp_fn}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(ppp, 'rt'), file_name=ppp_fn, mime='text/csv', key='primer_download_'+dpp_fn)
            
    if not taqwater_picklist_paths:
        error_msgs.append('No taq/water picklist available')
    else:
        for tpp in taqwater_picklist_paths:
            tpp_fn = Path(tpp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_fn}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(tpp, 'rt'), file_name=tpp_fn, mime='text/csv', key='taqwater_download_'+tpp_fn)

    if error_msgs:
        for msg in error_msgs:
            picklist_file_col.markdown(f'<p style="color:#ff0000;text-align:right">{msg}</p>',\
                    unsafe_allow_html=True)


def show_echo2_outputs():
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)
    index_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR2_picklist_filepaths()
    
    error_msgs = []

    if not index_picklist_paths:
        error_msgs.append('No index picklist available')
    else:
        for ipp in index_picklist_paths:
            ipp_fn = Path(ipp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ipp_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                    data=open(ipp, 'rt'), file_name=ipp_fn, mime='text/csv', key='index_download_'+ipp_fn)

    if not taqwater_picklist_paths:
        error_msgs.append('No taq/water picklist available')
    else:
        for tpp in taqwater_picklist_paths:
            tpp_fn = Path(tpp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                        data=open(tpp, 'rt'), file_name=tpp_fn, mime='text/csv',\
                                key='taqwater_download_'+tpp_fn)
    if error_msgs:
        for msg in error_msgs:
            picklist_file_col.markdown(f'<p style="color:#ff0000;text-align:right">{msg}</p>',\
                    unsafe_allow_html=True)

    
def display_status(key, height=300):
    """
    Display the progress in the pipeline for this experiment
    Should use aggrid to display the stages and the changes at each stage
    """
    exp = st.session_state['experiment']
    steps, header = exp.get_stages()
    status_df = pd.DataFrame(steps, columns=header)
    #status_df.reset_index(inplace=True)
    #status_df = status_df.rename(columns = {'index':'Steps', 'pending':'Pending Steps'})
    status_df = aggrid_interactive_table(status_df, grid_height=height, key=str(key)+'status_aggrid')


def view_plates(key, height=500):
    """
    Visual view of the plates in the experument
    """
    exp = st.session_state['experiment']
    plate_ids = []
    for pid in exp.plate_location_sample:
        plate_ids.append(f"{exp.plate_location_sample[pid]['purpose']} plate: {util.unguard(pid, silent=True)}")

    #Let user choose which plate to view
    _, col1, _ = st.columns([2, 1,2])
    plate_selectbox = col1.selectbox('Plate ID to view', plate_ids, key=str(key)+'plate_viewer')
    if plate_selectbox:
        plate_id = util.guard_pbc(plate_selectbox.split(':')[1])
        if plate_id in exp.plate_location_sample:
            heatmap_str = generate_heatmap_html(exp, plate_id, scaling=0.9)
            with open("makehtml.html", 'wt', encoding="utf-8") as outf:
                print(heatmap_str, file=outf)
            components.html(heatmap_str, height=height, scrolling=True)
        else:
            plate_barcode_error_msg = "Plate barcode not found in experiment"
            st.markdown(f'<p style="color:#FF0000">{plate_barcode_error_msg}</p>',
                        unsafe_allow_html=True)


def display_files(key, file_usage, height=250):
    """
    Display info for all files that have so far been uploaded into the experiment
    Give info on name, any plates they contain, whether they are required so far, etc
    """
    exp = st.session_state['experiment']
    file_df = pd.DataFrame.from_dict(file_usage, orient='index')
    if file_df is None or not isinstance(file_df, pd.DataFrame):
        st.write('No plates loaded')
        return
    else:
        if 'view_box_size' in st.session_state:
            height = st.session_state['view_box_size']
        else:
            height = height

    file_df.reset_index(inplace=True)
    file_df = file_df.rename(columns = {'index':'File', 'plates':'Plates', 'purpose':'Purpose'})
    selection = aggrid_interactive_table(file_df, grid_height=height, key=str(key)+'file_aggrid')
    if selection and 'selected_rows' in selection:
        fns = [row['File'] for row in selection['selected_rows']]
        fns = [fn for fn in fns if fn in exp.uploaded_files]
        if fns:
            #print(f'{fns=}')
            if 'previous_file_delete_selection' not in st.session_state:
                st.session_state['previous_file_delete_selection'] = None
            if st.session_state['previous_file_delete_selection'] == fns:
                st.session_state['previous_file_delete_selection'] = None
            else:
                if fns != st.session_state['previous_file_delete_selection']:
                    st.markdown(f"**You selected {fns}**")
                delbox = st.container()
                del_col1, del_col2, del_col3, _ = delbox.columns([2,1,1,4])
                del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                del_col2.button("Yes",on_click=manage_delete,
                        args=(delbox,'file',fns), key="delete " + str(key), help=f"Delete {fns}")
                del_col3.button("No", on_click=cancel_delete,
                        args=('file',fns), key="keep " + str(key), help=f"Keep {fns}")
                selection = None


def display_primers(key, dna_pids=None, primer_pids=None, height=350):
    """
    Alternative display_primer_components using aggrid
    Designed as a component for a generic display widget
    """
    exp = st.session_state['experiment']
    ul_conv = 1000
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids=dna_pids)
    primer_avail_vols_doses = exp.get_available_primer_vols_doses(pmr_pids=primer_pids)
    primer_info_array = []
    for primer in exp.primer_assayfam:
        req_vol = exp.transfer_volumes['PRIMER_VOL']*primer_usage.get(primer,0)  # nl
        req_doses = primer_usage.get(primer,0)
        req_wells = util.num_req_wells(req_vol)
        avail_vol = primer_avail_vols_doses.get(primer,[0,0])[0]  # nl
        avail_doses = primer_avail_vols_doses.get(primer,[0,0])[1]
        avail_wells = util.num_req_wells(avail_vol)
        info = [primer, req_doses, req_vol/ul_conv, req_wells, avail_doses, avail_vol/ul_conv, avail_wells]
        primer_info_array.append(info)

    primer_df = pd.DataFrame(primer_info_array, columns=['Primer', 'Required Doses', 
            'Required Volume (μL)', 'Required Wells', 'Available Doses', 'Available Volume (μL)', 'Available Wells'])
    primer_table = aggrid_interactive_table(primer_df, grid_height=height, key=str(key)+'primer_display')

        
def display_indexes(key, dna_pids=None, height=350):
    """
    Display info for all indexes that have been uploaded into the experiment
    """
    exp = st.session_state['experiment']
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids=dna_pids)
    fwd_idx, rev_idx = exp.get_index_avail()
    indexes = {**fwd_idx, **rev_idx}

    index_df = pd.DataFrame.from_dict(indexes,  orient='index')
    index_df.reset_index(inplace=True)
    index_df = index_df.rename(columns = {'index':'Index', 'well_count':'Wells', 'avail_transfers':'Available transfers', 'avail_vol':'Available Volume (μL)'})
    index_table = aggrid_interactive_table(index_df,grid_height=height,key=str(key)+'index_display')


def display_references(key, height=350):
    """
    Show the amplicon targets/reference sequences that are loaded
    """
    exp = st.session_state['experiment']
    refseq = st.session_state['experiment'].reference_sequences
    for group in refseq:
        dataset = [(ref,refseq[group][ref]) for ref in exp.reference_sequences[group]]
        ref_df = pd.DataFrame(dataset, columns=['Name','Sequence'])
        st.markdown('**'+group+'**')
        aggrid_interactive_table(ref_df, key=str(key)+'seqs'+group, grid_height=height)


def display_log(key, height=250):
    """
    Display the log
    """
    exp = st.session_state['experiment']
    #func = sys._getframe(1).f_code.co_name
    #func_line = inspect.getframeinfo(sys._getframe(1)).lineno
    #print("Function", func, type(func_line))
    # display the experiment log and let the user filter the view in a number of ways
    log_entries = st.session_state['experiment'].get_log(100)
    if len(log_entries) == 0:
        st.write('No entries currently in the log')
    else:
        df = pd.DataFrame(log_entries, columns=exp.get_log_header())
        aggrid_interactive_table(df, grid_height=height, key=str(key)+'logs')

def change_screens_open(screen_name):
    print(f'{screen_name=}')
    print(f'{st.session_state["screens_open"]=}')
    if screen_name in st.session_state['screens_open']:
        st.session_state['screens_open'].remove(screen_name)
    else:
        st.session_state['screens_open'].add(screen_name)

def info_bar(key):
    """
    Define a bar and populate state if the toggle in screens_open.
    """
    if 'screens_open' not in st.session_state:
        st.session_state['screens_open'] = set()
    info_on = False
    if 'info_on' in st.session_state:
        if st.session_state['info_on']:
            info_on = True

    info_cols = st.columns([1,7])

    with info_cols[0]:
        info_on = st.toggle('Info viewer', 
                value=info_on, 
                help='Show extra information toggles')
        
    with info_cols[1]:
        st.divider()

    if not info_on:
        st.session_state['info_on'] = False
        return

    if info_on:
        row1 = st.columns([1,1,1,1,1])
        with row1[0]:  # status
            already_on = False
            if 'status' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='Status', value=already_on, key='status_toggle'+str(key), 
                    on_change=change_screens_open, args=['status'], help='Display status window')

            already_on = False

            if 'primers' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='Primers', value=already_on, key='primers_toggle'+str(key), 
                    on_change=change_screens_open, args=['primers'], help='Display primers window')
            
        with row1[1]:  # files
            already_on = False
            if 'files' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='files', value=already_on, key='files_toggle'+str(key), 
                    on_change=change_screens_open, args=['files'], help='Display files window')
            already_on = False

            if 'indexes' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='Indexes', value=already_on, key='index_toggle'+str(key), 
                    on_change=change_screens_open, args=['indexes'], help='Display indexes window')

        with row1[2]:  # plates
            already_on = False

            if 'plates' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='Plates', value=already_on, key='plates_toggle'+str(key), 
                    on_change=change_screens_open, args=['plates'], help='Display plates window')

            already_on = False

            if 'references' in st.session_state['screens_open']:
                already_on = True

            on = st.toggle(label='References', value=already_on, key='references_toggle'+str(key), 
                    on_change=change_screens_open, args=['references'], help='Display status window')

        with row1[3]:  # plate viewer
            already_on = False

            if 'plate_viewer' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='Plate viewer', value=already_on, key='plate_view_toggle'+str(key), 
                    on_change=change_screens_open, args=['plate_viewer'], help='Display plate viewer window')

            already_on = False

            if 'log' in st.session_state['screens_open']:
                already_on = True
            on = st.toggle(label='Log', value=already_on, key='log_toggle'+str(key), 
                    on_change=change_screens_open, args=['log'], help='Display log window')

        with row1[4]:
            view_height = st.number_input('Set display height', min_value=50, max_value=700, 
                    value=350, step=25, help="Size of display grid", key=str(key))

        st.divider()


def info_viewer(key, dna_pids=None, pcr_pids=None, primer_pids=None, index_pids=None, amp_pids=None, taq_pids=None):
    """
    Container for displaying module info functions, each of which provides a dataframe for display in an aggrid.
    Because aggrid allows selection, each module can also handle a standard set of operations (such as delete).
    """
    exp = st.session_state['experiment']
    if 'info_expand' not in st.session_state:
        st.session_state['info_expand'] = False
    
    #info_expander = st.expander('Info Panel', expanded=st.session_state['info_expand'])
    container = st.container()
        
    col1,col2 = st.columns([7,1])
    with col2:
        view_height = st.number_input('Set display height', min_value=50, max_value=700, 
                value=350, step=25, help="Size of display grid", key=str(key))
    with col1:
        view_tab = stx.tab_bar(data=[
            stx.TabBarItemData(id=1, title="Status", description=""),
            stx.TabBarItemData(id=2, title="Files", description=""),
            stx.TabBarItemData(id=3, title="Plates", description=""),
            stx.TabBarItemData(id=4, title="Plate Viewer", description=""),
            stx.TabBarItemData(id=5, title="Primers", description=""),
            stx.TabBarItemData(id=6, title="Indexes", description=""),
            #stx.TabBarItemData(id=7, title="Reference sequences", description=""),
            stx.TabBarItemData(id=7, title="Log", description="")
        ], return_type=int, default=1)

    if view_tab == 1:
        # Status tab should tell us where we are up to in the pipeline and what's happened so far
        with container:
            display_status(key, height=view_height)
                
    if view_tab == 2:
        #view_expander = container.expander(label='All uploaded files', expanded=False)
        file_usage = exp.get_file_usage()
        with container:
            display_files(key, file_usage, height=view_height)


    if view_tab == 3:
        plate_usage = exp.get_plate_usage()
        with container:
            display_plates(key, plate_usage, height=view_height)

    if view_tab == 4:
        with container:
            view_height = 500
            view_plates(key, height=view_height)


    if view_tab == 5:
        with container:
            display_primers(key, dna_pids=dna_pids, height=view_height)
        
    if view_tab == 6:
        with container:
            display_indexes(key, dna_pids=dna_pids, height=view_height)


    #if view_tab == 7:
    #    with container:
    #        display_references(key, height=view_height)


    if view_tab == 7:
        with container:
            display_log(key, height=view_height)
    

def plate_checklist_pcr1(exp):
    """
    Allows the selection/deselection all plates involved in the PCR1 (primer) reaction stage
    DNA plates, PCR plates, taq/water plates, primer plates
    """
    pcr_stage = 1
    checklist_col = st.columns(4)
    included_DNA_plates = set()
    included_PCR_plates = set()
    included_taqwater_plates = set()
    
    pcr_plate_title = "PCR Plates"
    taqwater_plate_title = "Taq/Water Plates"
    dna_plate_title = "DNA Plates"
    
    #nimbus fp, echo fp, barcodes not in echo
    nfs, efs, xbcs = exp.get_nimbus_filepaths()

    #print(f'{efs=}')

    #missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
    #available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]
   
    checklist_col[0].markdown(f'**{dna_plate_title}**')
    for nim in efs:
        echo_filename=Path(nim).stem
        #print(echo_filename)
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
    
    #print(included_DNA_plates)
        
    return included_DNA_plates, included_PCR_plates, included_taqwater_plates


def plate_checklist_pcr2(exp):
    """
    Allows the selection/deselection all plates involved in the PCR2 (indexing) reaction stage
    DNA plates, PCR plates, taq/water plates, primer plates
    """
    pcr_stage = 2
    checklist_col = st.columns(4)
    included_PCR_plates = set()
    included_index_plates = set()
    included_taqwater_plates = set()
    included_amplicon_plates = set()
    #could make a for loop

    index_plate_title = "Index Plates"
    amplicon_plate_title = "Amplicon Plates"
    pcr_plate_title = "PCR Plates"
    taqwater_plate_title = "Taq/Water Plates"

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
 
    


