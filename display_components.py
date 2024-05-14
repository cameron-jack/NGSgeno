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
from tabnanny import check
import jsonpickle
from time import sleep
import datetime

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components
import extra_streamlit_components as stx
from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode

from stutil import custom_text, add_vertical_space, add_pm, m, init_state, mq
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
    options = GridOptionsBuilder.from_dataframe(
        df, enableRowGroup=True, enableValue=True, enablePivot=True,
    )

    if 'Message' in df.columns:
        options.configure_column(field = 'Message', width = 800)
    
    if 'Func line' in df.columns:
        options.configure_column(field = 'Func line', width = 70)
    
    options.configure_side_bar()

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
        fit_columns_on_grid_load=True
    )
    return selection


def manage_delete_cb(caller_id, category, ids):
    """
    Callback for deletion operations
    args:
        caller_id (str): the name of a message queue for a particular display widget
        category (str): file/plate/group - the resource type to be deleted
        ids (list): a list of items to be deleted
    """
    exp = st.session_state['experiment']
    successful_ids = []
    failed_ids = []
    init_state('message_queues',dict())
    if caller_id not in st.session_state['message_queues']:
        st.session_state['message_queues'][caller_id] = []
        
    if category == 'file':
        for id in ids:
            if id not in exp.uploaded_files:
                continue
            success = exp.del_file_record(id)     
            if success:
                successful_ids.append(id)
                m(f'{id} removed', level='info', dest=('log','console'))
            else:
                failed_ids.append(id)
                m(f'{id} could not be removed', level='warning', dest=('log','console'))
            st.session_state['previous_file_delete_selection'] = ids
    elif category == 'plate':
        for pid in ids:
            gid = util.guard_pbc(pid, silent=True)
            if gid not in exp.plate_location_sample:
                failed_ids.append(gid)
            else:
                success = exp.delete_plate(gid)
                if success:
                    successful_ids.append(gid)
                else:
                    failed_ids.append(gid)
        st.session_state['previous_plate_delete_selection'] = ids
    elif category == 'group':  # from summary
        for row in ids:
            dest_pid = util.guard_pbc(row['DNA/amplicon PID'], silent=True)
            if dest_pid in exp.plate_location_sample and exp.plate_location_sample[dest_pid]['purpose'] == 'amplicon':
                success = exp.delete_plate(dest_pid)
                if success:
                    successful_ids.append(dest_pid)
                else:
                    failed_ids.append(dest_pid)
            elif dest_pid not in exp.dest_sample_plates:
                m(f"{dest_pid=} doesn't actually exist in the experiment!",
                        level='error', dest=('log','toast','debug'))
                failed_ids.append(dest_pid)
            else:
                sample_pids = exp.dest_sample_plates[dest_pid]
                delete_pids = sample_pids + [dest_pid]
                for pid in delete_pids:
                    success = exp.delete_plate(pid)
                    if success:
                        successful_ids.append(pid)
                    else:
                        failed_ids.append(pid)
        st.session_state['previous_group_delete_selection'] = ids
    # set up messages
    for sid in successful_ids:
        m(f'{sid} removed', level='info', dest=('log','console','noGUI'))
        st.session_state['message_queues'][caller_id].append((f'{sid}','success'))
    for fid in failed_ids:
        m(f'{id} could not be removed', level='failure', dest=('log','console'))
        st.session_state['message_queues'][caller_id].append((f'{fid}','failure'))
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


# def display_temporary_messages():
#     """ 
#     Display any temporary user alerts
#     Messages are tuples of message, level. Where level: info/warning/error/success
#     """
#     if st.session_state['messages_temp']:
#         for message, level in st.session_state['messages']:
#             if level is None:
#                 st.markdown(message)
#             elif level.lower() == 'info':
#                 st.info(message)
#             elif level.lower() == 'warning':
#                 st.warning(message)
#             elif level.lower() == 'error':
#                 st.error(message)
#             elif level.lower() == 'success':
#                 st.success(message)
#             else:
#                 st.markdown(message)
#     # display only once
#     st.session_state['messages_temp'] = []
    

def display_persistent_messages(key):
    """
    Display any persistent user alerts, and allow the user to choose which ones to clear
    Messages are tuples of message, level. Where level: info/warning/errror/success
    A key is required to ensure the form and checkboxes are unqiue
    """
    if st.session_state['messages_persist']:
        with st.form(key):
            m('**System messages**', dest=('mkdn',))
            check_col, message_col = st.columns([1,10])
            for message, level in st.session_state['messages_persist']:
                with check_col:
                    st.checkbox(message,key=(key,message,level), label_visibility="collapsed")
                with message_col:
                    if level is None:
                        st.markdown(message)
                    elif level.lower() == 'info':
                        st.info(message)
                    elif level.lower() == 'warning':
                        st.warning(message)
                    elif level.lower() == 'error':
                        st.error(message)
                    elif level.lower() == 'success':
                        st.success(message)
                    else:
                        st.markdown(message)
            submitted = st.form_submit_button("Clear ticked items")
        if submitted:
            kept_messages = []
            for message, level in st.session_state['messages_persist']:
                checkbox_name = (key, message, level)
                if checkbox_name not in st.session_state or not st.session_state[checkbox_name]:
                     kept_messages.append((message,level))
            st.session_state['messages_persist'] = kept_messages              
    

def display_samples(key, height=250):
    """ display a summary of all loaded DNA, amplicon, and sample plates """
    exp = st.session_state['experiment']
    caller_id = 'display_samples'
    selection = []
    df = exp.inputs_as_dataframe()
    if df is None or not isinstance(df, pd.DataFrame):
        m('**No 384-well DNA plate data loaded**', dest=('mkdn',))
    else:    
        selection = aggrid_interactive_table(df, key=key, grid_height=height)
        if 'selected_rows' in selection and selection['selected_rows']:
            rows = selection["selected_rows"]
            rows = [r for r in rows if 'DNA/amplicon PID' in r and r['DNA/amplicon PID'] != 'Total']
            if rows:
                if 'previous_group_delete_selection' not in st.session_state:
                    st.session_state['previous_group_delete_selection'] = None
                if st.session_state['previous_group_delete_selection'] != rows:
                    #st.session_state['previous_group_delete_selection'] = rows
                    # only do the code below if this is a fresh selection
                    lines = '\n'.join(['DNA/amplicon PID: '+r['DNA/amplicon PID'] for r in rows if r['DNA/amplicon PID'] != 'Total'])
                    m(f"**You selected {lines}**", dest=('mkdn',))
                    delbox = st.container()
                    del_col1, del_col2, del_col3, _ = st.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    del_col2.button("Yes",on_click=manage_delete_cb, 
                            args=(caller_id, 'group',rows), key="delete " + str(key), help=f"Delete {lines}")
                    del_col3.button("No",on_click=cancel_delete, args=('group',rows), 
                            key="keep " + str(key), help=f"Keep {lines}")
    selection = None
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []

    init_state('message_queues', dict())
    if caller_id in st.session_state['message_queues']:
        for msg,lvl in st.session_state['message_queues'][caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        st.session_state['message_queues'][caller_id] = []
        

def display_consumables(key, height=300):
    """
    summarise_consumables():
        d = {'taqwater_pids_pcr1':[], 'taqwater_pids_pcr2':[], 'taq_vol_pcr1':0, 'taq_vol_pcr2':0,'water_vol_pcr1':0, 
                'water_vol_pcr2':0, 'primer_pids':[], 'primer_count_ngs':0, 'primer_count_custom':0, 'unique_primers':set(), 
                'primer_well_count':0, 'assay_primer_mappings':0, 'reference_files':[], 'unique_references':set(), 
                'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
    """
    exp = st.session_state['experiment']
    caller_id = 'display_consumables'
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
        m('No plates loaded')
    else:
        selection = aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'consumables_aggrid')
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def display_plates(key, plate_usage, height=300): 
    exp = st.session_state['experiment']
    caller_id = 'display_plates'
    plate_df = pd.DataFrame(plate_usage, columns=['Plates', 'Num Wells', 'Purpose'])
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        m('No plates loaded')
    else:
        selection = aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'plate_aggrid')
        if selection and 'selected_rows' in selection:
            pids = [row['Plates'] for row in selection['selected_rows']]
            gids = [util.guard_pbc(pid, silent=True) for pid in pids]
            gids = [gid for gid in gids if gid in exp.plate_location_sample]
            if gids:
                if 'previous_plate_delete_selection' not in st.session_state:
                    st.session_state['previous_plate_delete_selection'] = None
                if st.session_state['previous_plate_delete_selection'] == pids:
                    st.session_state['previous_plate_delete_selection'] = None
                else:
                    if pids != st.session_state['previous_plate_delete_selection']:
                        m(f"**You selected {pids}**", dest=('mkdn'))
                    delbox = st.container() # doesn't work reliably in st1.26
                    del_col1, del_col2, del_col3, _ = delbox.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    del_col2.button("Yes",on_click=manage_delete_cb,
                            args=(caller_id,'plate',pids), key="delete " + str(key), help=f"Delete {pids}")
                    del_col3.button("No", on_click=cancel_delete,
                            args=('plate',pids), key="keep " + str(key), help=f"Keep {pids}")
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def display_pcr_common_components(selected_pids):
    """
    Expander widget that shows the PCR plate information for both PCR reactions
    Args:
        selected_pids (dict): '
    """
    exp = st.session_state['experiment']
    caller_id = 'display_pcr_common_components'
    PCR_PLATE_WELLS = 384
    dna_pids = selected_pids['dna']
    pcr_pids = selected_pids['pcr']
    amplicon_pids = selected_pids['amplicon']
    index_pids = selected_pids['index']

    if dna_pids:
        dna_pids = [util.guard_pbc(dp, silent=True) for dp in dna_pids]
    else:
        dna_pids = []
        
    if pcr_pids:
        pcr_pids = [util.guard_pbc(pp, silent=True) for pp in pcr_pids]
    else:
        pcr_pids = []
        
    if amplicon_pids:
        amplicon_pids = [util.guard_pbc(ap, silent=True) for ap in amplicon_pids]
    else:
        amplicon_pids = []
    
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids)
    dna_reactions = sum([primer_usage[p] for p in primer_usage])
    
    #PCR
    required_pcr_plates = ceil(dna_reactions/PCR_PLATE_WELLS)

    num_supplied_pcr = 0
    supplied_pcr_txt = '<p style="color:#FF0000">None</p>'
    if len(pcr_pids) > 0:
        supplied_pcr_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in pcr_pids])
        num_supplied_pcr = len(pcr_pids)
        
    # amplicons
    amplicon_pid_txt = 'None'
    if len(amplicon_pids) > 0:
        amplicon_pid_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in amplicon_pids])
    
    req_PCR_text = '**Number of required PCR plates**'
    req_PCR_num = str(required_pcr_plates)
    if num_supplied_pcr < required_pcr_plates:
        req_PCR_text = '<p style="color:#FF0000"><b>Number of required PCR plates</b></p>'
        req_PCR_num = f'<p style="color:#FF0000">{str(required_pcr_plates)}</p>'
    
    #Page set up
    pcr_comps_area = st.container()
    col_size = [6, 4, 6, 4]
    pcr_cols = pcr_comps_area.columns(col_size)

    pcr_cols[0].markdown(req_PCR_text, unsafe_allow_html=True)
    pcr_cols[1].markdown(req_PCR_num, unsafe_allow_html=True)

    pcr_cols[0].markdown('**Worst case required reaction wells**')
    pcr_cols[1].write(str(dna_reactions), unsafe_allow_html=True)

    pcr_cols[2].markdown('**User supplied PCR plates**')
    pcr_cols[3].markdown(supplied_pcr_txt, unsafe_allow_html=True)
    
    pcr_cols[2].markdown('**User supplied amplicon plates**')
    pcr_cols[3].markdown(amplicon_pid_txt)
    
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def display_pcr1_components(selected_pids):
    """
    Display panel that shows the required componenents for each PCR reaction, 
    including wells, PCR plates, taq+water plates
    Args:
        selected_pids (dict): 'dna','pcr','taqwater1'
    """
    exp = st.session_state['experiment']
    caller_id = 'display_pcr1_components'
    ul_conv = 1000
    pcr_stage = 1
    PCR_PLATE_WELLS = 384
    
    dna_pids = selected_pids['dna']
    pcr_pids = selected_pids['pcr']
    taqwater1_pids = selected_pids['taqwater1']

    if dna_pids:
        dna_pids = [util.guard_pbc(dp, silent=True) for dp in dna_pids]
    else:
        dna_pids = []
        
    if pcr_pids:
        pcr_pids = [util.guard_pbc(pp, silent=True) for pp in pcr_pids]
    else:
        pcr_pids = []
        
    #Page set up
    pcr_comps_area = st.container()
    col_size = [6, 4, 6, 4]
    req_cols = pcr_comps_area.columns(col_size)
    pcr_cols = pcr_comps_area.columns(col_size)

    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids)
    num_reactions = sum([primer_usage[p] for p in primer_usage])

     #PCR
    required_pcr_plates = ceil(num_reactions/PCR_PLATE_WELLS)

    #Required taq/water
    primer_taq_vol, primer_water_vol = exp.get_taqwater_volumes_primer(num_reactions)

    taq_avail, water_avail, pids = exp.get_taqwater_avail(taqwater_bcs=taqwater1_pids)
    taq_avail_vol = taq_avail/ul_conv
    water_avail_vol = water_avail/ul_conv            

    #get actual values for volume of taq water plates
    required_water_vol_str = str(primer_water_vol/ul_conv) + ' μl'
    water_avail_vol_str = str(water_avail_vol)+' μl'
    required_taq_vol_str = str(primer_taq_vol/ul_conv) + ' μl'
    avail_taq_vol_str = str(taq_avail_vol)+' μl'

    num_req_taq_water_plates = util.num_req_taq_water_plates(primer_taq_vol, primer_water_vol)

    user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True) \
                            for p in taqwater1_pids])
    
    num_supplied_taqwater = len(user_supplied_taqwater)
    user_taqwater_text = user_supplied_taqwater if user_supplied_taqwater else '<p style="color:#FF0000">None</p>'
    
    if num_supplied_taqwater < num_req_taq_water_plates:
        req_taqwater_text = '<p style="color:#FF0000"><b>Number of required taq/water plates</b></p>'
        req_taqwater_num = f'<p style="color:#FF0000">{str(num_req_taq_water_plates)}</p>'
    else:
        req_taqwater_text = '**Number of required taq/water plates**'
        req_taqwater_num = str(num_req_taq_water_plates)


    num_supplied_pcr = 0
    supplied_pcr_txt = '<p style="color:#FF0000">None</p>'
    if len(pcr_pids) > 0:
        supplied_pcr_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in pcr_pids])
        num_supplied_pcr = len(pcr_pids)

    req_PCR_text = '**Number of required PCR plates**'
    req_PCR_num = str(required_pcr_plates)
    if num_supplied_pcr < required_pcr_plates:
        req_PCR_text = '<p style="color:#FF0000"><b>Number of required PCR plates</b></p>'
        req_PCR_num = f'<p style="color:#FF0000">{str(required_pcr_plates)}</p>'



    for i in range(4):
        req_cols[i].write('')


    req_cols[0].markdown(req_PCR_text, unsafe_allow_html=True)
    req_cols[1].markdown(req_PCR_num, unsafe_allow_html=True)

    req_cols[2].markdown('**User supplied PCR plates**')
    req_cols[3].markdown(supplied_pcr_txt, unsafe_allow_html=True)

    req_cols[0].markdown('**Worst case required reaction wells**')
    req_cols[1].write(str(num_reactions), unsafe_allow_html=True)

    with req_cols[2]:
        add_vertical_space(3)
    with req_cols[3]:
        add_vertical_space(3)
    
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
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def display_pcr2_components(selected_pids):
    """
    Expander widget that shows the required componenents for PCR 2 reaction (index).
    Args:
        selected_pids (dict): 'dna','pcr','amplicon','taqwater2','index'
        pcr_stage (1, 2): 1 = Echo Primer stage, 2 = Echo Indexing
    """
    #Need to add info about taq water
    exp = st.session_state['experiment']
    caller_id = 'display_pcr2_components'
    ul_conv = 1000
    pcr_stage = 2
    dna_pids = selected_pids['dna']
    pcr_pids = selected_pids['pcr']
    amplicon_pids = selected_pids['amplicon']
    taqwater2_pids = selected_pids['taqwater2']
    index_pids = selected_pids['index']

    if pcr_pids:
        pcr_pids = [util.guard_pbc(pp, silent=True) for pp in pcr_pids]
    else:
        pcr_pids = []
    
    if amplicon_pids:
        amplicon_pids = [util.guard_pbc(ap, silent=True) for ap in amplicon_pids]
    else:
        amplicon_pids = []
    

    num_reactions = exp.get_num_reactions(pcr_pids = pcr_pids, amplicon_pids = amplicon_pids)
    index_max = len(exp.get_index_pairs_avail(index_pids))
    index_remain = index_max - num_reactions

    #Taq/water (based on pcr stage)
    index_taq_vol, index_water_vol = exp.get_taqwater_req_vols_index(num_reactions)
    
    #user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True)\
    #                                    for p in exp.get_taqwater_avail(pcr_stage=pcr_stage)[2]])
    if taqwater2_pids is None:
        user_supplied_taqwater = ''
    else:
        user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in taqwater2_pids])    

    num_supplied_taqwater = len(user_supplied_taqwater)
                
    taq_avail, water_avail, pids = exp.get_taqwater_avail(taqwater_bcs=taqwater2_pids)
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
    
    supplied_pcr_txt = '<p style="color:#FF0000">None</p>'
    if len(pcr_pids) > 0:
        supplied_pcr_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in pcr_pids])
        
    amplicon_pid_txt = 'None'
    if len(amplicon_pids) > 0:
        amplicon_pid_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in amplicon_pids])


    #Page set up
    pcr_comps_area = st.container()
    col_size = [6, 4, 6, 4]
    pcr2_col = pcr_comps_area.columns(col_size)

    pcr2_col[0].markdown('**User supplied PCR plates**')
    pcr2_col[1].markdown(supplied_pcr_txt, unsafe_allow_html=True)
    pcr2_col[2].markdown('**User supplied amplicon plates**')
    pcr2_col[3].markdown(amplicon_pid_txt, unsafe_allow_html=True)

    pcr2_col[0].markdown(req_taqwater_text, unsafe_allow_html=True)
    pcr2_col[1].markdown(req_taqwater_num, unsafe_allow_html=True)
    pcr2_col[2].markdown(f'**User supplied taq/water plates (PCR {pcr_stage})**',unsafe_allow_html=True)
    pcr2_col[3].markdown(user_taqwater_text, unsafe_allow_html=True)

    pcr2_col[0].markdown('**Required water volume**')
    pcr2_col[1].write(required_water_vol_str)
    pcr2_col[2].markdown('**Available water volume**')
    pcr2_col[3].write(water_avail_vol_str)
    pcr2_col[0].markdown('**Required taq volume**')
    pcr2_col[1].markdown(required_taq_vol_str)
    pcr2_col[2].markdown('**Available taq volume**')
    pcr2_col[3].write(avail_taq_vol_str)

    pcr2_col[0].markdown('**Index Pairs Available**')
    pcr2_col[1].write(str(index_max))
    pcr2_col[0].markdown('**Index Pairs Remaining**')
    pcr2_col[1].markdown(index_pairs_remain, unsafe_allow_html=True)
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []
    

def st_directory_picker(label='Selected directory:', initial_path=Path(),\
            searched_file_types=['fastq','fastq.gz','fq','fq.gz']):
    """
    Streamlit being JS/AJAX has no ability to select a directory. This is for server paths only.
    Initial code by Aidin Jungo here: https://github.com/aidanjungo/StreamlitDirectoryPicker
    """
    caller_id = 'st_directory_picker'
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
            st.rerun()

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
                st.rerun()
    st.markdown(\
            f'<h5 style="color:#000000">Contains {len(contains_files)} files of type {",".join(searched_file_types)}</h5>',\
                    unsafe_allow_html=True)
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []

    return st.session_state['path']


def handle_picklist_download(picklist_type, picklist_paths, error_msgs, file_col, btn_col):
    if not picklist_paths:
        error_msgs.append(f"No {picklist_type} picklist available")
    else:
        for ppp in picklist_paths:
            ppp_fn = Path(ppp).name
            with file_col:
                custom_text('p', '#4b778c', ppp_fn, 'right', padding='5px', display=True)
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
            custom_text('p', '#ff0000', msg, 'center', padding='5px', display=True)
            

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
            custom_text('p', '#ff0000', msg, 'center', padding='5px', display=True)


def show_echo1_outputs():
    exp = st.session_state['experiment']
    caller_id = 'show_echo1_outputs'
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
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def show_echo2_outputs():
    exp = st.session_state['experiment']
    caller_id = 'show_echo2_outputs'
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
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []

    
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
    # _, col1, _ = st.columns([1,2,1])
    plate_selectbox = st.selectbox('Plate ID to view', plate_ids, key=str(key)+'plate_viewer')
    if plate_selectbox:
        plate_id = util.guard_pbc(plate_selectbox.split(':')[1])
        if plate_id in exp.plate_location_sample:
            if exp.plate_location_sample[plate_id]['purpose'] == 'amplicon':
                exp.log(f'DEBUG: {exp.get_plate(plate_id)=}')
            jsonpickle_plate = jsonpickle.encode(exp.get_plate(plate_id), keys=True)
            heatmap_str = generate_heatmap_html(jsonpickle_plate, plate_id, scaling=0.9)
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
    caller_id = 'display_files'
    file_df = pd.DataFrame.from_dict(file_usage, orient='index')
    if file_df is None or not isinstance(file_df, pd.DataFrame):
        st.write('No plates loaded')
        return
   
    file_df.reset_index(inplace=True)
    #file_df = file_df.rename(columns = {'index':'File', 'plates':'Plates', 'purpose':'Purpose'})
    file_df = file_df.rename(columns = {'index':'File', 'date modified':'Date modified', 'purpose':'Purpose'})
    if not file_df.empty:
        file_df.sort_values(by='Date modified', inplace=True, ascending=False)

    selection = aggrid_interactive_table(file_df, grid_height=height, key=str(key)+'file_aggrid')
    if selection and 'selected_rows' in selection:
        fns = [row['File'] for row in selection['selected_rows']]
        fns = [fn for fn in fns if fn in exp.uploaded_files]
        if fns:
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
                del_col2.button("Yes",on_click=manage_delete_cb,
                        args=(caller_id,'file',fns), key="delete " + str(key), help=f"Delete {fns}")
                del_col3.button("No", on_click=cancel_delete,
                        args=('file',fns), key="keep " + str(key), help=f"Keep {fns}")
                selection = None
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl)
        sleep(0.3)
        mq[caller_id] = []


def display_primers(key, dna_pids=None, primer_pids=None, height=350):
    """
    Alternative display_primer_components using aggrid
    Designed as a component for a generic display widget
    """
    exp = st.session_state['experiment']
    ul_conv = 1000
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids)
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
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids)
    fwd_idx, rev_idx = exp.get_index_avail()
    indexes = {**fwd_idx, **rev_idx}

    index_df = pd.DataFrame.from_dict(indexes,  orient='index')
    index_df.reset_index(inplace=True)
    index_df = index_df.rename(columns = {'index':'Index', 
                                          'well_count':'Wells', 
                                          'avail_transfers':'Available transfers', 
                                          'avail_vol':'Available Volume (μL)'})
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
    log_entries = st.session_state['experiment'].get_log(100)
    if len(log_entries) == 0:
        st.write('No entries currently in the log')
    else:
        df = pd.DataFrame(log_entries, columns=exp.get_log_header())
        df = df.drop(['Calling function', 'Call line'], axis = 1)
        aggrid_interactive_table(df, grid_height=height, key=str(key)+'logs')

# def change_screens_open(screen_name):
#     """
#     DEPRECATED - retain for potential future projects
#     Required support function for info_bar()
#     """
#     if screen_name in st.session_state['screens_open']:
#         st.session_state['screens_open'].remove(screen_name)
#     else:
#         st.session_state['screens_open'].add(screen_name)


# def info_bar(key):
#     """
#     DEPRECATED - retain for potential future projects
#     Define a bar and populate state if the toggle in screens_open.
#     """
#     if 'screens_open' not in st.session_state:
#         st.session_state['screens_open'] = set()
#     info_on = False
#     if 'info_on' in st.session_state:
#         if st.session_state['info_on']:
#             info_on = True

#     info_cols = st.columns([1,7])

#     with info_cols[0]:
#         info_on = st.toggle('Info viewer', 
#                 value=info_on, 
#                 help='Show extra information toggles')
        
#     with info_cols[1]:
#         st.divider()

#     if not info_on:
#         st.session_state['info_on'] = False
#         return

#     if info_on:
#         row1 = st.columns([1,1,1,1,1])
#         with row1[0]:  # status
#             already_on = False
#             if 'status' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='Status', value=already_on, key='status_toggle'+str(key), 
#                     on_change=change_screens_open, args=['status'], help='Display status window')

#             already_on = False

#             if 'primers' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='Primers', value=already_on, key='primers_toggle'+str(key), 
#                     on_change=change_screens_open, args=['primers'], help='Display primers window')
            
#         with row1[1]:  # files
#             already_on = False
#             if 'files' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='files', value=already_on, key='files_toggle'+str(key), 
#                     on_change=change_screens_open, args=['files'], help='Display files window')
#             already_on = False

#             if 'indexes' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='Indexes', value=already_on, key='index_toggle'+str(key), 
#                     on_change=change_screens_open, args=['indexes'], help='Display indexes window')

#         with row1[2]:  # plates
#             already_on = False

#             if 'plates' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='Plates', value=already_on, key='plates_toggle'+str(key), 
#                     on_change=change_screens_open, args=['plates'], help='Display plates window')

#             already_on = False

#             if 'references' in st.session_state['screens_open']:
#                 already_on = True

#             on = st.toggle(label='References', value=already_on, key='references_toggle'+str(key), 
#                     on_change=change_screens_open, args=['references'], help='Display status window')

#         with row1[3]:  # plate viewer
#             already_on = False

#             if 'plate_viewer' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='Plate viewer', value=already_on, key='plate_view_toggle'+str(key), 
#                     on_change=change_screens_open, args=['plate_viewer'], help='Display plate viewer window')

#             already_on = False

#             if 'log' in st.session_state['screens_open']:
#                 already_on = True
#             on = st.toggle(label='Log', value=already_on, key='log_toggle'+str(key), 
#                     on_change=change_screens_open, args=['log'], help='Display log window')

#         with row1[4]:
#             view_height = st.number_input('Set display height', min_value=50, max_value=700, 
#                     value=350, step=25, help="Size of display grid", key=str(key))

#         st.divider()


# def info_viewer_old(key, dna_pids=None, pcr_pids=None, primer_pids=None, index_pids=None, amp_pids=None, taq_pids=None):
#     """
#     Container for displaying module info functions, each of which provides a dataframe for display in an aggrid.
#     Because aggrid allows selection, each module can also handle a standard set of operations (such as delete).
#     DEPRECATED: Maintain as an alternative display system
#     """
#     exp = st.session_state['experiment']
#     if 'info_expand' not in st.session_state:
#         st.session_state['info_expand'] = False
    
#     #info_expander = st.expander('Info Panel', expanded=st.session_state['info_expand'])
#     container = st.container()
        
#     col1,col2 = st.columns([7,1])
#     with col2:
#         view_height = st.number_input('Set display height', min_value=50, max_value=700, 
#                 value=350, step=25, help="Size of display grid", key=str(key))
#     with col1:
#         view_tab = stx.tab_bar(data=[
#             stx.TabBarItemData(id=1, title="Status", description=""),
#             stx.TabBarItemData(id=2, title="Files", description=""),
#             stx.TabBarItemData(id=3, title="Plates", description=""),
#             stx.TabBarItemData(id=4, title="Plate Viewer", description=""),
#             stx.TabBarItemData(id=5, title="Primers", description=""),
#             stx.TabBarItemData(id=6, title="Indexes", description=""),
#             #stx.TabBarItemData(id=7, title="Reference sequences", description=""),
#             stx.TabBarItemData(id=7, title="Log", description="")
#         ], return_type=int, default=1)

#     if view_tab == 1:
#         # Status tab should tell us where we are up to in the pipeline and what's happened so far
#         with container:
#             display_status(key, height=view_height)
                
#     if view_tab == 2:
#         #view_expander = container.expander(label='All uploaded files', expanded=False)
#         file_usage = exp.get_file_usage()
#         with container:
#             display_files(key, file_usage, height=view_height)


#     if view_tab == 3:
#         plate_usage = exp.get_plate_usage()
#         with container:
#             display_plates(key, plate_usage, height=view_height)

#     if view_tab == 4:
#         with container:
#             view_height = 500
#             view_plates(key, height=view_height)

#     if view_tab == 5:
#         with container:
#             display_primers(key, dna_pids=dna_pids, height=view_height)
        
#     if view_tab == 6:
#         with container:
#             display_indexes(key, dna_pids=dna_pids, height=view_height)

#     if view_tab == 7:
#         with container:
#             display_log(key, height=view_height)


def info_viewer(selection, key, dna_pids=None, view_height=350):
    exp = st.session_state['experiment']

    if selection == 'Samples':
        display_samples(key=key+selection, height=view_height)
    
    if selection == 'Consumables':
        display_consumables(key=key+selection, height=view_height)

    if selection == "Status":
        # Status tab should tell us where we are up to in the pipeline and what's happened so far
            display_status(key=key+selection, height=view_height)
                
    if selection == "Files":
        file_usage = exp.get_file_usage()
        display_files(key+selection, file_usage, height=view_height)

    if selection == "Plates":
        plate_usage = exp.get_plate_usage()
        display_plates(key+selection, plate_usage, height=view_height)

    if selection == "Plate Viewer":
        view_height = 500
        view_plates(key+selection, height=view_height)

    if selection == "Primers":
        display_primers(key+selection, dna_pids=dna_pids, height=view_height)
        
    if selection == "Indexes":
        display_indexes(key+selection, dna_pids=dna_pids, height=view_height)

    if selection == "Log":
        display_log(key+selection, height=view_height)
    

def set_state(key, value):
    """ Callback function for display elements """
    st.session_state[key] = value


def set_selection(set_key, widget_key):
    """ Callback function for selectbox display elements """
    if widget_key in st.session_state:
        st.session_state[set_key] = st.session_state[widget_key]


def info_selection(key, view1_key, view2_key, height_key, default_view1="None", 
        default_view2="None", default_height=250):
    """
    Container for displaying module info functions, each of which provides a dataframe for display in an aggrid.
    Because aggrid allows selection, each module can also handle a standard set of operations (such as delete).
    Function handles a single case and extra keys are needed for correct naming and identification
    key - the general key applied to this set of widgets
    view1_key - the key for selection1 lookups in other code
    view2_key - the key for selection2 lookups in other code
    height_key - the key for height lookups in other code
    defaults are returned, otherwise these are set via callbacks
    """
    key=str(key)
    
    options = ["None","Samples", "Consumables", "Status", "Files", "Plates", "Plate Viewer", 
            "Primers", "Indexes", "Log"]
    
    disp_col1, disp_col2, height_col = st.columns([4,4,2])
    
    with disp_col1:
        select1_key = key+"_select1"
        selection1 = st.selectbox("Choose info to view", options=options, placeholder='',
                index=options.index(default_view1), on_change=set_selection, 
                args=[view1_key, select1_key], key=select1_key)
        set_state(view1_key, selection1)
        
    with disp_col2:
        select2_key = key+"_select2"    
        selection2 = st.selectbox("Choose info to view", options=options, placeholder='',
                index=options.index(default_view2), on_change=set_selection,
                args=[view2_key, select2_key], key=select2_key)
        set_state(view2_key, selection2)
        
    with height_col:
        height_widget_key = key+"_height"
        view_height = st.number_input('Set display height', min_value=50, max_value=700, 
                value=default_height, step=25, help="Size of display grid", on_change=set_selection,
                args=[height_key, height_widget_key], key=height_widget_key)
        set_state(height_key, view_height)        

    return True


def show_info_viewer(selection, height, groupkey):
    """ Display any number of info view windows side-by-side. Select must be an iterable container """
    if len(selection) != 0:
        columns = st.columns(len(selection))
        for i in range(len(selection)):
            with columns[i]:
                info_viewer(selection[i], str(groupkey)+selection[i]+str(i), view_height=height)


def display_plate_checklist(widget_key:str, inc_plate_types:list) -> list:
    """
    Display a plate checklist, with each included category getting its own column
    Returns a list of all checkbox keys for later lookup
    By choosing which types to select we can customise this for PCR1 or PCR2
    args:
        widget_key (str): a unique id for this widget
        inc_plate_types (list[str]): a list of plate types to include ['dna','pcr','primer',
                'index','taqwater1','taqwater2','amplicon']
    return:
        a list of checkbox ID strings
    """
    exp = st.session_state['experiment']
    checkbox_keys = []
    checklist_cols = st.columns(len(inc_plate_types))
    for i, ipt in enumerate(inc_plate_types): 
        if ipt == 'dna':
            checklist_cols[i].markdown('**DNA Plates**')
            for dp in exp.get_dna_pids(echo_ready=True):
                cb_name = f'{str(widget_key)}_plate_checkbox_dna_{dp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(dp, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
        if ipt == 'pcr':
            checklist_cols[i].markdown('**PCR Plates**')
            for pp in exp.get_pcr_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_pcr_{pp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(pp, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
        if ipt == 'taqwater1':
            checklist_cols[i].markdown('**Taq/Water Plates (PCR 1)**')
            for tp in exp.get_taqwater_pids(pcr_stage=1):
                cb_name = f'{str(widget_key)}_plate_checkbox_taqwater1_{tp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(tp, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
        if ipt == 'taqwater2':
            checklist_cols[i].markdown('**Taq/Water Plates (PCR 2)**')
            for tp in exp.get_taqwater_pids(pcr_stage=2):
                cb_name = f'{str(widget_key)}_plate_checkbox_taqwater2_{tp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(tp, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
        if ipt == 'amplicon':
            checklist_cols[i].markdown('**Amplicon Plates**')
            for ap in exp.get_amplicon_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_amplicon_{ap}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(ap, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
        if ipt == 'primer':
            checklist_cols[i].markdown('**Primer Plates**')
            for pp in exp.get_primer_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_primer_{pp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(pp, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
        if ipt == 'index':
            checklist_cols[i].markdown('**Index Plates**')
            for ip in exp.get_index_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_index_{ip}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(ip, silent=True), 
                        key=cb_name, value=True)
                checkbox_keys.append(cb_name)
    return checkbox_keys
    

def collect_plate_checklist(checkbox_keys):
    """
    Return all PIDs that have been selected by the given checkboxes
    args:
        checkbox_keys (list[str])
    """
    selected_pids = {'dna':[],'pcr':[],'taqwater1':[], 'taqwater2':[],'amplicon':[],'primer':[],'index':[]}
    for cb in checkbox_keys:
        if st.session_state[cb]:
            key_parts = cb.split('_')
            if key_parts[-2] == 'dna':
                selected_pids['dna'].append(key_parts[-1])
            elif key_parts[-2] == 'pcr':
                selected_pids['pcr'].append(key_parts[-1])
            elif key_parts[-2] == 'taqwater1':
                selected_pids['taqwater1'].append(key_parts[-1])
            elif key_parts[-2] == 'taqwater2':
                selected_pids['taqwater2'].append(key_parts[-1])
            elif key_parts[-2] == 'amplicon':
                selected_pids['amplicon'].append(key_parts[-1])
            elif key_parts[-2] == 'primer':
                selected_pids['primer'].append(key_parts[-1])
            elif key_parts[-2] == 'index':
                selected_pids['index'].append(key_parts[-1])
            else:
                m(f'Plate selection checkbox key {cb} of unknown type', level='critical', dest=('log','debug','noGUI'))
    return selected_pids
 

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


def show_upper_info_viewer_checkbox(widget_key, value=True): #, default_panel1='None', default_panel2='None'):
    """
    Allows the user to turn the upper info viewer panel on and off
    The bottom info viewer display is always on
    Allows the default display panels to be set (page specific content)
    args:
        widget_key (str): a unique identifier to prevent same-page clashes
        value (bool): the default value of the enabling checkbox
    """
    if 'show_info_viewer' not in st.session_state:
        st.session_state['show_upper_info_viewer'] = True
    
    if st.checkbox('Info Viewer', value=value, key=widget_key):
        st.session_state['show_upper_info_viewer'] = True
        # if default_panel1:
        #     st.session_state['info_panel1'] == default_panel1
        # if default_panel2:
        #     st.session_state['info_panel2'] == default_panel2
    else:
        st.session_state['show_upper_info_viewer'] = False


def set_nimbus_title(exp, efs, nfs):
    """
    *Stage 2: Nimbus*
    Title for nimbus stage
    Args:
        exp (st.session_state['experiment'])
        efs (str): file path to echo files
        nfs (str): file path for nimbus files
    """
    title = ''
    colour = '#f63366'
    #first stage sample files haven't been loaded
    if not st.session_state['experiment'].dest_sample_plates:
        title = "Load data inputs to enable Nimbus input file generation."
    else:
        # do we have any Nimbus inputs to generate + download
        echo_files_exist = len(efs) == len(nfs) and len(efs) != 0
        yet_to_run = len(exp.dest_sample_plates) - len(nfs)

        if echo_files_exist:
            title = 'All Echo inputs received.'
            colour = '#83b3c9'
        if yet_to_run > 0:
            title = f'For {str(yet_to_run)} 96-well plate set(s)'
            colour = '#83b3c9'
    
    m(title, dest=('css',), size='h5', color=colour, align='left')
        

def get_echo_download_buttons(nfs):
    """
    *Stage 2: Nimbus*
    Generates the echo file download buttons
    Args:
        nfs (str): nimbus file paths
    """

    if len(nfs) < 5:
        _,dl_col1,dl_col2,_= st.columns([6,3,2,6])
        for i, nf in enumerate(nfs):
            nimbus_fn=Path(nf).name
            
            with dl_col1:
                custom_text("p", "#4b778c", nimbus_fn, "left", display=True)
                add_vertical_space(1)

            dl_col2.download_button("Download ", 
                                    open(nf), 
                                    file_name=nimbus_fn, 
                                    key='nimbus_input_1_'+str(i), 
                                    help=f"Download Nimbus input file {nf}")


    else:
        _,dl_col1,dl_col2,dl_col3,dl_col4,_= st.columns([3,6,6,6,6,3])

        for i,nf in enumerate(nfs):
            nimbus_fn=Path(nf).name

            if (i+1) % 2 != 0:
                with dl_col1:
                    
                    custom_text("p", "#4b778c", nimbus_fn, "left", display=True)
                    add_vertical_space(1)

                dl_col2.download_button("Download ", 
                                        open(nf), 
                                        file_name=nimbus_fn, 
                                        key='nimbus_input_2_'+str(i), 
                                        help=f"Download Nimbus input file {nf}")
        
            else:
                with dl_col3:
                    custom_text("p", "#4b778c", nimbus_fn, "left", display=True)
                    add_vertical_space(1)
            
                dl_col4.download_button("Download ", 
                                        open(nf), file_name=nimbus_fn,\
                                        key='nimbus_input'+str(i), 
                                        help=f"Download Nimbus input file {nf}")
            

def get_miseq_download_btn(exp):
    """
    *Stage 5: Miseq*
    Args:
        exp (st.session_state['experiment])
    """
    add_vertical_space(2)
    _, miseq_col1,_, miseq_col2, _ =  st.columns([5,2,1,3,3])
    for fp in exp.get_miseq_samplesheets():
        
        fp_name = str(Path(fp).name)
        with miseq_col1:
            add_vertical_space(1)
            custom_text(size='h5', color='#cf3276', text=fp_name, align='right', display=True)
        
        with miseq_col2:
            add_vertical_space(1)
            download_miseq = st.download_button(label='Download', 
                                                data=open(fp, 'rt'), 
                                                file_name=fp_name, 
                                                mime='text/csv', 
                                                key='dnld_samplesheet_'+str(fp), 
                                                type='secondary')



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
    
    margins_css = """
    <style>
        .appview-container .main .block-container {
            padding-top: 1rem;
            padding-bottom: 1rem;
            padding-left: 2rem;
            padding-right: 2rem;
            }

    </style>
    """
    st.markdown(margins_css, unsafe_allow_html=True)

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
    


