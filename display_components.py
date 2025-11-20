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
import re
from pathlib import PurePath, Path
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO
import inspect
from tabnanny import check
import jsonpickle
from time import sleep
import datetime
import jsonpickle
import warnings

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import streamlit as st
import streamlit.components.v1 as components
import extra_streamlit_components as stx
from st_aggrid import AgGrid, GridOptionsBuilder, JsCode
from st_aggrid.shared import GridUpdateMode

from stutil import custom_text, add_vertical_space, hline, m, init_state, mq

try:
    from bin.match import reconstruct_sequence, run_msa
except ModuleNotFoundError:
    from match import build_aligned_pair, reconstruct_sequence, run_msa
try:
    from bin.generate import generate_pdf
except ModuleNotFoundError:
    from generate import generate_pdf
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


def aggrid_interactive_table(df: pd.DataFrame, grid_height: int=250, hidden: list=[],editable: list=[],
        filtering: bool=False, key: int=1):
    """Creates an st-aggrid interactive table based on a dataframe.

    Args:
        df (pd.DataFrame]): Source dataframe
        grid_height (int, optional): Height of the grid in pixels. Defaults to 250.
        hidden (list, optional): List of column names to hide. Defaults to [].
        editable (list, optional): List of column names to make editable. Defaults to [].
        filtering (bool, optional): Whether to enable filtering column for alleles. Defaults to False.
        key (int, optional): A unique key for the table. Defaults to 1.
    Returns:
        dict: The selected row

    Special behaviour for log tables to colour significant events and set better column ordering
    and for available primers/indexes
    """
    cols = df.columns.tolist()
    if 'Level' in cols:
        #print(cols, flush=True)
        cols = [cols[0], cols[3], cols[4], cols[1], cols[2]]
        df = df[cols]

    options = GridOptionsBuilder.from_dataframe(
        df, enableRowGroup=True, enableValue=True, enablePivot=True)

    # params.data['report'] = kept_indices.map(i => raw_names[i]).join(';');
    if filtering and 'otherCount' in cols and 'mergeCount' in cols and \
            'filtProportion' in cols:
        filt_js = JsCode("""
            function(params) {
                let raw_counts = params.data.otherCount.split(';').map(Number);
                let raw_names = params.data.otherName.split(';');
                let kept_indices = new Array();
                for (let i = 0; i < raw_counts.length; i++) {
                    if (raw_counts[i] > Number(params.data.filtProportion) * Number(params.data.mergeCount)) {
                        kept_indices.push(i);
                    }
                }
                params.data.report = kept_indices.map(i => raw_names[i]).join(';');
                return kept_indices.map(i => raw_names[i]).join(';');
            };""")


        options.configure_column('report', valueGetter=filt_js)

    if 'Available Doses' in cols and 'Required Doses' in cols:
        cell_js = JsCode("""
            function(params) {
                // mark inadequate amounts
                if (params.data['Required Doses'] > params.data['Available Doses']) {
                     return {'color': 'white', 'backgroundColor': 'red'}
                } else if (params.data['Required Doses'] < params.data['Available Doses']) {
                     return {'color': 'black', 'backgroundColor': 'white'}
                }
            };""")
        options.configure_column('Primer', cellStyle=cell_js)

    if 'Level' in cols:
        cell_js = JsCode("""
            function(params) {
            // different styles for each row
                if (params.value === 'Error') {
                    //mark Error cells as red
                    return {'color': 'white', 'backgroundColor': 'red'}
                } else if (params.value === 'Critical') {
                    //mark Critical cells as purple
                    return {'color': 'white', 'backgroundColor': 'purple'}
                } else if (params.value === 'Warning') {
                    //mark Warning cells as yellow
                    return {'color': 'black', 'backgroundColor': 'yellow'}
                } else {
                    return {'color': 'black', 'backgroundColor': 'white'}
                }
            };""")
        options.configure_column('Level', cellStyle=cell_js)

    if hidden:
        for hid in hidden:
            if hid in df.columns:
                # hide the column
                options.configure_column(field=hid, hide=True)
            else:
                m(f"Column {hid} not found in dataframe", level='Critical', dest='nogui')

    if editable:
        for ed in editable:
            if ed in df.columns:
                # make editable
                options.configure_column(field=ed, editable=True)
            else:
                m(f"Column {ed} not found in dataframe", level='Critical', dest='nogui')

    if 'Message' in df.columns:
        options.configure_column(field = 'Message', width = 800)

    if 'Func line' in df.columns:
        options.configure_column(field = 'Func line', width = 70)

    options.configure_side_bar()

    options.configure_selection("single", use_checkbox=False, \
                rowMultiSelectWithClick=False, suppressRowDeselection=True)

    selection = None
    #df = df.astype(str)
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
        fit_columns_on_grid_load=True,
        rowSelection='single'
    )
    #rowSelection='multiple',
    #    selection_mode='multiple',
    #    rowMultiSelectWithClick=True,
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

    if category == 'file':
        for id in ids:
            if id not in exp.uploaded_files:
                continue
            success = exp.del_file_record(id)
            if success:
                successful_ids.append(id)
                m(f'{id} removed', level='info')
            else:
                failed_ids.append(id)
                m(f'{id} could not be removed', level='warning')
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
        for pid in ids:
            dest_pid = util.guard_pbc(pid, silent=True)
            if dest_pid in exp.plate_location_sample and exp.plate_location_sample[dest_pid]['purpose'] == 'amplicon':
                success = exp.delete_plate(dest_pid)
                if success:
                    successful_ids.append(dest_pid)
                else:
                    failed_ids.append(dest_pid)
            elif dest_pid not in exp.dest_sample_plates:
                m(f"{dest_pid=} doesn't actually exist in the experiment!", level='error')
                failed_ids.append(dest_pid)
            else:
                sample_pids = exp.dest_sample_plates[dest_pid]
                delete_pids = sample_pids + [dest_pid]
                for dpid in delete_pids:
                    success = exp.delete_plate(dpid)
                    if success:
                        successful_ids.append(dpid)
                    else:
                        failed_ids.append(dpid)
        st.session_state['previous_group_delete_selection'] = ids
    # set up messages
    for sid in successful_ids:
        m(f'{sid} removed', level='success', caller_id=caller_id)
    for fid in failed_ids:
        m(f'{fid} could not be removed', level='error', caller_id=caller_id)
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
    DEPRECATED
    Display any persistent user alerts, and allow the user to choose which ones to clear
    Messages are tuples of message, level. Where level: info/warning/errror/success
    A key is required to ensure the form and checkboxes are unqiue
    """
    if st.session_state['messages_persist']:
        with st.form(key):
            st.markdown('**System messages**')
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


def display_samples(key, height=250, caller_id=None):
    """
    Info bar display a summary of all loaded DNA, amplicon, and sample plates
    Use only the provided caller_id
    """
    exp = st.session_state['experiment']
    selection = []
    df = exp.inputs_as_dataframe()
    if df is None or not isinstance(df, pd.DataFrame):
        st.markdown('**No 384-well DNA plate data loaded**')
    else:
        selection = aggrid_interactive_table(df, key=key, grid_height=height)
        if selection is not None:
            if 'selected_rows' in selection and selection['selected_rows'] is not None:
                rows = [r for r in selection["selected_rows"].get('DNA/amplicon PID') if r != 'Total']
                if rows:
                    if 'previous_group_delete_selection' not in st.session_state:
                        st.session_state['previous_group_delete_selection'] = None
                    if st.session_state['previous_group_delete_selection'] != rows:
                        # only do the code below if this is a fresh selection
                        lines = '\n'.join(['DNA/amplicon PID: '+r for r in rows])
                        st.markdown(f"**You selected {lines}**")
                        del_col1, del_col2, del_col3, _ = st.columns([2,1,1,4])
                        del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                        del_col2.button("Yes",on_click=manage_delete_cb,
                                args=(caller_id, 'group',rows), key="delete " + str(key), help=f"Delete {lines}")
                        del_col3.button("No",on_click=cancel_delete, args=('group',rows),
                                key="keep " + str(key), help=f"Keep {lines}")
    selection = None


def display_consumables(key, height=300, caller_id=None):
    """
    info bar display of "consumables" non-sample plate info
    Use only the provided caller_id
    summarise_consumables():
        d = {'taqwater_pids_pcr1':[], 'taqwater_pids_pcr2':[], 'taq_vol_pcr1':0, 'taq_vol_pcr2':0,'water_vol_pcr1':0,
                'water_vol_pcr2':0, 'primer_pids':[], 'primer_count_ngs':0, 'primer_count_custom':0, 'unique_primers':set(),
                'primer_well_count':0, 'assay_primer_mappings':0, 'rodentity_reference_files':[],
                'custom_reference_files':[], 'amplicon_reference_files':[],
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
    for f in consumables['rodentity_reference_files']:
        data_rows.append(['Rodentity references', f, 0, 'File', exp.count_reference_sequences('rodentity_references', silent=True)])
    for f in consumables['custom_reference_files']:
        data_rows.append(['Custom references', f, 0, 'File', exp.count_reference_sequences('custom_references', silent=True)])
    for f in consumables['amplicon_reference_files']:
        data_rows.append(['Amplicon references', f, 0, 'File', exp.count_reference_sequences('amplicon_references', silent=True)])
    for f in exp.uploaded_files:
        if exp.uploaded_files[f].get('purpose','') == 'assay_primer_map':
            data_rows.append(['assay-primer mappings', f, 0, 'File', consumables['assay_primer_mappings']])
    plate_df = pd.DataFrame(data_rows, columns=headers)
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        st.write('No plates loaded')
    else:
        selection = aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'consumables_aggrid')


def display_plates(key, plate_usage, height=300, caller_id=None):
    """
    Info bar display of plates
    Use only the provided caller_id
    """
    exp = st.session_state['experiment']
    caller_id = 'display_plates'
    plate_df = pd.DataFrame(plate_usage, columns=['Plates', 'Num Wells', 'Purpose'])
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        st.write('No plates loaded')
    else:
        selection = aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'plate_aggrid')
        if selection is not None and 'selected_rows' in selection:
            if selection['selected_rows'] is not None:
                pids = [pid for pid in selection['selected_rows'].get('Plates') if pid is not None]
                gids = [util.guard_pbc(pid, silent=True) for pid in pids]
                gids = [gid for gid in gids if gid in exp.plate_location_sample]
                if gids:
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
                        del_col2.button("Yes",on_click=manage_delete_cb,
                                args=(caller_id,'plate',pids), key="delete " + str(key), help=f"Delete {pids}")
                        del_col3.button("No", on_click=cancel_delete,
                                args=('plate',pids), key="keep " + str(key), help=f"Keep {pids}")


def display_pcr1_components(selected_pids, caller_id=None):
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

    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids, caller_id=caller_id)
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
    user_taqwater_text = user_supplied_taqwater if user_supplied_taqwater else ':red[None]'

    taq_water_needed = max(num_req_taq_water_plates - num_supplied_taqwater, 0)
    if taq_water_needed > 0:
        req_taqwater_text = ':red[**Remaining taq/water plates needed**]'
        req_taqwater_num = f':red[{str(taq_water_needed)}]'
    else:
        req_taqwater_text = '**Remaining taq/water plates needed**'
        req_taqwater_num = str(taq_water_needed)


    num_supplied_pcr = 0
    supplied_pcr_txt = ':red[None]'
    if len(pcr_pids) > 0:
        supplied_pcr_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in pcr_pids])
        num_supplied_pcr = len(pcr_pids)

    pcr_plates_needed = max(required_pcr_plates - num_supplied_pcr,0)
    req_PCR_text = '**Remaining PCR plates needed**'
    req_PCR_num = str(pcr_plates_needed)
    if pcr_plates_needed:
        req_PCR_text = ':red[**Remaining PCR plates needed**]'
        req_PCR_num = f':red[{pcr_plates_needed}]'

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
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()


def display_pcr2_components(selected_pids, caller_id=None):
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
        user_taqwater_text = ':red[None]'

    taq_water_needed = max(num_req_taq_water_plates - num_supplied_taqwater, 0)
    if taq_water_needed > 0:
        req_taqwater_text = f':red[**Remaining taq/water plates needed (PCR {pcr_stage})**]'
        req_taqwater_num = f':red[{str(taq_water_needed)}]'
    else:
        req_taqwater_text = f'**Remaining taq/water plates needed (PCR {pcr_stage})**'
        req_taqwater_num = str(taq_water_needed)

    if index_remain >= 0:
        index_pairs_remain = f':green[{str(index_remain)}]'
    else:
        index_pairs_remain = f':red[{str(index_remain)}]'

    supplied_pcr_txt = ':red[None]'
    if len(pcr_pids) > 0:
        supplied_pcr_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in pcr_pids])

    amplicon_pid_txt = 'None'
    if len(amplicon_pids) > 0:
        amplicon_pid_txt = ', '.join([util.unguard_pbc(p, silent=True)\
                for p in amplicon_pids])


    #Page set up
    pcr_comps_area = st.container()
    col_size = [4, 2, 6, 6]
    pcr2_col = pcr_comps_area.columns(col_size)

    pcr2_col[0].markdown('**Index Pairs Available**')
    pcr2_col[1].write(str(index_max))
    pcr2_col[2].markdown('**Index Pairs Remaining**')
    pcr2_col[3].markdown(index_pairs_remain)

    for i in range(4):
        pcr2_col[i].write('')

    pcr2_col[0].markdown(req_taqwater_text)
    pcr2_col[1].markdown(req_taqwater_num)
    pcr2_col[2].markdown(f'**User supplied taq/water plates (PCR {pcr_stage})**')
    pcr2_col[3].markdown(user_taqwater_text)

    pcr2_col[0].markdown('**Required water volume**')
    pcr2_col[1].write(required_water_vol_str)
    pcr2_col[2].markdown('**Available water volume**')
    pcr2_col[3].write(water_avail_vol_str)
    pcr2_col[0].markdown('**Required taq volume**')
    pcr2_col[1].markdown(required_taq_vol_str)
    pcr2_col[2].markdown('**Available taq volume**')
    pcr2_col[3].write(avail_taq_vol_str)

    for i in range(4):
        pcr2_col[i].write('')

    pcr2_col[0].markdown('**User supplied amplicon plates**')
    pcr2_col[1].markdown(amplicon_pid_txt)
    pcr2_col[2].markdown('**User supplied PCR plates**')
    pcr2_col[3].markdown(supplied_pcr_txt)



    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()


def st_directory_picker(label='Selected directory:', initial_path=Path(),\
            searched_file_types=['fastq','fastq.gz','fq','fq.gz'], caller_id=None):
    """
    Streamlit being JS/AJAX has no ability to select a directory. This is for server paths only.
    Initial code by Aidin Jungo here: https://github.com/aidanjungo/StreamlitDirectoryPicker
    DEPRECATED - too ugly
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
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()

    return st.session_state['path']


def handle_picklist_download(picklist_type, picklist_paths, file_col, btn_col, caller_id=None):
    if not picklist_paths:
        m(f"No {picklist_type} picklist available", level='error', no_log=True, caller_id=caller_id)
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


def get_echo1_download_btns(caller_id=None):
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)

    dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR1_picklist_filepaths()

    picklist_dict = {'DNA': dna_picklist_paths, 'primer': primer_picklist_paths, 'taq/water': taqwater_picklist_paths}

    for pltype, plpath in picklist_dict.items():
        handle_picklist_download(pltype, plpath, picklist_file_col, picklist_btn_col, caller_id=caller_id)


def get_echo2_download_btns(caller_id=None):
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)

    index_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR2_picklist_filepaths()

    picklist_dict = {'index': index_picklist_paths, 'taq/water': taqwater_picklist_paths}

    for pltype, plpath in picklist_dict.items():
        handle_picklist_download(pltype, plpath, picklist_file_col, picklist_btn_col, caller_id=caller_id)


def show_echo1_outputs(caller_id=None):
    exp = st.session_state['experiment']
    caller_id = 'show_echo1_outputs'
    picklist_file_col, picklist_btn_col = st.columns(2)
    dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR1_picklist_filepaths()

    if not dna_picklist_paths:
        m('No DNA picklist available', level='error', no_log=True, caller_id=caller_id)
    else:
        for dpp in dna_picklist_paths:
            dpp_fn = Path(dpp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{dpp_fn}</p>',\
                    unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",
                    data=open(dpp, 'rt'), file_name=dpp_fn, mime='text/csv', key='dna_download_'+dpp_fn)

    if not primer_picklist_paths:
        m('No primer picklist available', level='error', no_log=True, caller_id=caller_id)
    else:
        for ppp in primer_picklist_paths:
            ppp_fn = Path(ppp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ppp_fn}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",
                    data=open(ppp, 'rt'), file_name=ppp_fn, mime='text/csv', key='primer_download_'+dpp_fn)

    if not taqwater_picklist_paths:
        m('No taq/water picklist available', level='error', no_log=True, caller_id=caller_id)
    else:
        for tpp in taqwater_picklist_paths:
            tpp_fn = Path(tpp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_fn}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",
                    data=open(tpp, 'rt'), file_name=tpp_fn, mime='text/csv', key='taqwater_download_'+tpp_fn)

    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()


def show_echo2_outputs(caller_id=None):
    exp = st.session_state['experiment']
    caller_id = 'show_echo2_outputs'
    picklist_file_col, picklist_btn_col = st.columns(2)
    index_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR2_picklist_filepaths()

    if not index_picklist_paths:
        m('No index picklist available', level='error', no_log=True, caller_id=caller_id)
    else:
        for ipp in index_picklist_paths:
            ipp_fn = Path(ipp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ipp_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                    data=open(ipp, 'rt'), file_name=ipp_fn, mime='text/csv', key='index_download_'+ipp_fn)

    if not taqwater_picklist_paths:
        m('No taq/water picklist available', level='error', no_log=True, caller_id=caller_id)
    else:
        for tpp in taqwater_picklist_paths:
            tpp_fn = Path(tpp).name
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                        data=open(tpp, 'rt'), file_name=tpp_fn, mime='text/csv',\
                                key='taqwater_download_'+tpp_fn)
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()


def display_status(key, height=300, caller_id='display_feedback'):
    """
    Display the progress in the pipeline for this experiment
    Should use aggrid to display the stages and the changes at each stage
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    steps, header = exp.get_stages()
    status_df = pd.DataFrame(steps, columns=header)
    #status_df.reset_index(inplace=True)
    #status_df = status_df.rename(columns = {'index':'Steps', 'pending':'Pending Steps'})
    status_df = aggrid_interactive_table(status_df, grid_height=height, key=str(key)+'status_aggrid')


def view_plates(key, height=500, caller_id='display_feedback'):
    """
    Visual view of the plates in the experument
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    plate_ids = []
    for pid in exp.plate_location_sample:
        plate_ids.append(f"{exp.plate_location_sample[pid]['purpose']} plate: {util.unguard(pid, silent=True)}")

    #Let user choose which plate to view
    # _, col1, _ = st.columns([1,2,1])
    plate_selectbox = st.selectbox('Plate ID to view', plate_ids, key=str(key)+'plate_viewer')
    if plate_selectbox is not None:
        plate_id = util.guard_pbc(plate_selectbox.split(':')[1])
        if plate_id in exp.plate_location_sample:
            jsonpickle_plate = jsonpickle.encode(exp.get_plate(plate_id), keys=True)
            heatmap_str = generate_heatmap_html(jsonpickle_plate, plate_id, scaling=0.9)
            with open("makehtml.html", 'wt', encoding="utf-8") as outf:
                print(heatmap_str, file=outf)
            components.html(heatmap_str, height=height, scrolling=True)
        else:
            plate_barcode_error_msg = "Plate barcode not found in experiment"
            st.markdown(f'<p style="color:#FF0000">{plate_barcode_error_msg}</p>',
                        unsafe_allow_html=True)


def display_files(key, file_usage, height=250, caller_id='display_feedback'):
    """
    Display info for all files that have so far been uploaded into the experiment
    Give info on name, any plates they contain, whether they are required so far, etc
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
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
    if selection is not None:
        if 'selected_rows' in selection and selection['selected_rows'] is not None:
            fns = [fp for fp in selection['selected_rows'].get('File')]
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


def display_primers(key, dna_pids=None, primer_pids=None, height=350, save_buttons=False, caller_id='display_feedback'):
    """
    Alternative display_primer_components using aggrid
    Designed as a component for a generic display widget and as a standalone with optional save to file
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    ul_conv = 1000
    if not dna_pids:
        dna_pids = exp.get_dna_pids(caller_id=caller_id)
    if not primer_pids:
        primer_pids = exp.get_primer_pids(caller_id=caller_id)
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids, caller_id=caller_id)
    primer_avail_vols_doses = exp.get_available_primer_vols_doses(pmr_pids=primer_pids, caller_id=caller_id)
    primer_wells = exp.get_available_primer_wells(pmr_pids=primer_pids, caller_id=caller_id)
    pmr_pids_wells = {pmr:{} for pmr in primer_wells}
    primer_info_array = []
    for primer in exp.primer_assayfam:
        req_vol = exp.transfer_volumes['PRIMER_VOL']*primer_usage.get(primer,0)  # nl
        req_doses = primer_usage.get(primer,0)
        req_wells = util.num_req_wells(req_vol)
        avail_vol = primer_avail_vols_doses.get(primer,[0,0])[0]  # nl
        avail_doses = primer_avail_vols_doses.get(primer,[0,0])[1]
        avail_wells = len([pmr_info for pmr_info in primer_wells.get(primer, [])])
        if req_wells == 0 and avail_wells == 0:
            continue
        if primer not in primer_wells:
            pos_line = ''
        else:
            for ppid, well, vol, dose in primer_wells[primer]:
                if ppid not in pmr_pids_wells[primer]:
                    pmr_pids_wells[primer][ppid] = []
                pmr_pids_wells[primer][ppid].append(well)
            ppid_well_lines = [f"{util.unguard_pbc(ppid, silent=True)}:\
                    {','.join(pmr_pids_wells[primer][ppid])}" for ppid in pmr_pids_wells[primer]]
            pos_line = ' '.join(ppid_well_lines)
        info = [primer, req_doses, req_vol/ul_conv, req_wells, avail_doses, avail_vol/ul_conv,
                avail_wells, pos_line]
        primer_info_array.append(info)

    primer_df = pd.DataFrame(primer_info_array, columns=['Primer', 'Required Doses',
            'Required Volume (μL)', 'Required Wells', 'Available Doses', 'Available Volume (μL)', 'Available Wells', 'Positions'])
    primer_table = aggrid_interactive_table(primer_df, grid_height=height, key=str(key)+'primer_display')
    primer_csv = primer_df.to_csv(index=False) #.encode('utf-8')
    primer_list_fn = exp.get_exp_fn('primer_list.csv', caller_id=caller_id)
    if 'primer_csv' not in st.session_state:
        st.session_state['primer_csv'] = 0
    if st.session_state['primer_csv'] != primer_csv:
        st.session_state['primer_csv'] = primer_csv
        primer_df.to_csv(path_or_buf=primer_list_fn, index=False) #.encode('utf-8')
    st.write(f'Primer table written to file: {primer_list_fn}')


def display_indexes(key, dna_pids=None, height=350, caller_id='display_feedback'):
    """
    Display info for all indexes that have been uploaded into the experiment
    Messages are displayed by the component holding these widgets
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


def choose_reference_files(key, height=350, subset='all', default=True, caller_id=None):
    """
    Allow multiselect of reference files to be used in the experiment
    Args:
        key (str): key for the component
        height (int): height of the table
        subset (str): 'all','amplicon','custom','rodentity' to show all or only the selected purpose
        default (bool): whether to include all by default
        caller_id (str): ID for the caller, used for messages
    Notes: Not used, instead use dc.display_amplicon_file_checklist()
    """
    exp = st.session_state['experiment']
    fns = [fn for fn, purp in exp.reference_sequences if purp == subset or subset == 'all']
    if default is True:
        options = st.multiselect(
                "Which amplicon sequence reference files to use?",
                fns, default=fns,
        )
    else:
        options = st.multiselect(
                "Which amplicon sequence reference files to use?",
                fns, default=[],
        )
    return options


def display_references(key, height=350, subset='all', caller_id='display_feedback'):
    """
    Show the amplicon targets/reference sequences that are loaded
    Args:
        key (str): key for the component
        height (int): height of the table
        subset (str): 'all','amplicon','custom','rodentity' to show all or only the selected purpose
        caller_id (str): ID for the caller, used for messages
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    dataset = []
    for fn,purp in exp.reference_sequences:
        if subset == 'amplicon' and purp != 'amplicon_reference':
            dataset.extend(exp.reference_sequences[(fn,purp)])
        elif subset == 'custom' and purp != 'custom_reference':
            # custom references are not used in the pipeline, we should report a mistake
            dataset.extend(exp.reference_sequences[(fn,purp)])
            m(f'Custom reference sequences are not used in the pipeline, please check the calling code',level='critical')
        elif subset == 'rodentity' and purp != 'rodentity_reference':
            dataset.extend(exp.reference_sequences[(fn,purp)])
        elif subset == 'all':
            dataset.extend(exp.reference_sequences[(fn,purp)])
    ref_df = pd.DataFrame(dataset, columns=['Name','Sequence'])
    st.markdown('**'+group+'**')
    aggrid_interactive_table(ref_df, key=str(key)+'seqs'+group, grid_height=height)


def display_log(key, height=250, caller_id='display_feedback'):
    """
    Display the log
    Messages are displayed by the component holding these widgets
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    log_entries = st.session_state['experiment'].get_log(100)
    if len(log_entries) == 0:
        st.write('No entries currently in the log')
    else:
        df = pd.DataFrame(log_entries, columns=exp.get_log_header())
        df = df.drop(['Calling function', 'Call line'], axis = 1)
        aggrid_interactive_table(df, grid_height=height, key=str(key)+'logs')


def info_viewer(selection, key, dna_pids=None, view_height=350):
    exp = st.session_state['experiment']
    caller_id = 'info_viewer_'+str(key)
    if selection == 'Samples':
        display_samples(key=key+selection, height=view_height, caller_id=caller_id)

    if selection == 'Consumables':
        display_consumables(key=key+selection, height=view_height, caller_id=caller_id)

    if selection == "Status":
        # Status tab should tell us where we are up to in the pipeline and what's happened so far
            display_status(key=key+selection, height=view_height, caller_id=caller_id)

    if selection == "Files":
        file_usage = exp.get_file_usage()
        display_files(key+selection, file_usage, height=view_height, caller_id=caller_id)

    if selection == "Plates":
        plate_usage = exp.get_plate_usage()
        display_plates(key+selection, plate_usage, height=view_height, caller_id=caller_id)

    if selection == "Plate Viewer":
        view_height = 500
        view_plates(key+selection, height=view_height, caller_id=caller_id)

    if selection == "Primers":
        display_primers(key+selection, dna_pids=dna_pids, height=view_height, caller_id=caller_id)

    if selection == "Indexes":
        display_indexes(key+selection, dna_pids=dna_pids, height=view_height, caller_id=caller_id)

    if selection == "Log":
        display_log(key+selection, height=view_height, caller_id=caller_id)


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
                with st.spinner('loading info...'):
                    info_viewer(selection[i], str(groupkey)+selection[i]+str(i), view_height=height)


def display_file_checklist(widget_key:str, inc_file_purposes:list, default_value=True) -> list:
    """
    Display a filename checklist, with each included category getting its own column
    Returns a list of all checkbox keys for later lookup
    args:
        widget_key (str): a unique id for this widget
        inc_file_purposes (list[str]): a list of file purposes to include
                see parse.process_upload for the list of types
                purposes: ['amplicon','DNA','pcr','rodentity_sample','custom_sample','primer_layout','primer_volume',
                'index_layout','index_volume','primer_assay_map','rodentity_reference','amplicon_reference','taq_water']
        default_value (bool): whether to set the checkboxes to True by default
    return:
        a list of checkbox ID strings
    """
    exp = st.session_state['experiment']
    checkbox_keys = []
    checklist_cols = st.columns(len(inc_file_purposes))
    all_fns = [fp for fp in exp.uploaded_files.keys() if fp not in {'_upload_queue','_upload_pending'}]

    for i, ift in enumerate(inc_file_purposes):
        ift_txt = ift.replace('_',' ').capitalize()
        checklist_cols[i].markdown(f'**{ift_txt} Files**')
        for fp in [fp for fp in all_fns if exp.uploaded_files[fp]['purpose'] == ift]:
            cb_name = f'{str(widget_key)}_file_checkbox_^_{ift}_^_{fp}'
            val = checklist_cols[i].checkbox(util.fn_from_path(fp), key=cb_name, value=default_value)
            checkbox_keys.append(cb_name)
    #print(f"Checkbox keys for {widget_key}: {checkbox_keys}", file=sys.stderr)
    return checkbox_keys


def collect_file_checklist(checkbox_keys:list) -> dict:
    """
    Return all filenames that have been selected by the given checkboxes
    args:
        checkbox_keys (list[str])
    """
    selected_files = {'amplicon':[],'DNA':[],'pcr':[],'rodentity_sample':[],
            'custom_sample':[],'primer_layout':[],'primer_volume':[], 'index_layout':[],
            'index_volume':[],'primer_assay_map':[],'rodentity_reference':[],
            'amplicon_reference':[],'taq_water':[]}
    for cb in checkbox_keys:
        if st.session_state[cb]:
            key_parts = cb.split('_^_')
            if len(key_parts) != 3:
                continue  # skip malformed keys
            purpose = key_parts[-2]
            fp = key_parts[-1]
            if purpose in selected_files:
                selected_files[purpose].append(fp)
            else:
                m(f'File selection checkbox key {cb} of unknown type', level='critical', dest=('log','debug','noGUI'))
    return selected_files


def fns_from_checklist(selected_files:dict) -> list:
    """
    Return a list of filenames from the selected files dictionary
    args:
        selected_files (dict): a dictionary of file purposes and lists of filenames
    return:
        a list of filenames
    """
    fps = []  # filepaths
    for f_type in selected_files:
        fps.extend(selected_files[f_type])
    return [util.fn_from_path(fp) for fp in fps]


def display_plate_checklist(widget_key:str, inc_plate_types:list, default_value=True) -> list:
    """
    Display a plate checklist, with each included category getting its own column
    Returns a list of all checkbox keys for later lookup
    By choosing which types to select we can customise this for PCR1 or PCR2
    args:
        widget_key (str): a unique id for this widget
        inc_plate_types (list[str]): a list of plate types to include ['dna','pcr','primer',
                'index','taqwater1','taqwater2','amplicon']
        default_value (bool): whether to set the checkboxes to True by default
    return:
        a list of checkbox ID strings
    """
    exp = st.session_state['experiment']
    checkbox_keys = []
    checklist_cols = st.columns(len(inc_plate_types))
    for i, ipt in enumerate(inc_plate_types):
        if ipt == 'dna':
            checklist_cols[i].markdown('**DNA Plates**')
            for dp in exp.get_dna_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_dna_{dp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(dp, silent=True),
                        key=cb_name, value=default_value)
                checkbox_keys.append(cb_name)
        if ipt == 'pcr':
            checklist_cols[i].markdown('**PCR Plates**')
            for pp in exp.get_pcr_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_pcr_{pp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(pp, silent=True),
                        key=cb_name, value=default_value)
                checkbox_keys.append(cb_name)
        if ipt == 'taqwater1':
            checklist_cols[i].markdown('**Taq/Water Plates (PCR 1)**')
            for tp in exp.get_taqwater_pids(pcr_stage=1):
                cb_name = f'{str(widget_key)}_plate_checkbox_taqwater1_{tp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(tp, silent=True),
                        key=cb_name, value=default_value)
                checkbox_keys.append(cb_name)
        if ipt == 'taqwater2':
            checklist_cols[i].markdown('**Taq/Water Plates (PCR 2)**')
            for tp in exp.get_taqwater_pids(pcr_stage=2):
                cb_name = f'{str(widget_key)}_plate_checkbox_taqwater2_{tp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(tp, silent=True),
                        key=cb_name, value=default_value)
                checkbox_keys.append(cb_name)
        if ipt == 'amplicon':
            checklist_cols[i].markdown('**Amplicon Plates**')
            for ap in exp.get_amplicon_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_amplicon_{ap}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(ap, silent=True),
                        key=cb_name, value=default_value)
                checkbox_keys.append(cb_name)
        if ipt == 'primer':
            checklist_cols[i].markdown('**Primer Plates**')
            for pp in exp.get_primer_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_primer_{pp}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(pp, silent=True),
                        key=cb_name, value=default_value)
                checkbox_keys.append(cb_name)
        if ipt == 'index':
            checklist_cols[i].markdown('**Index Plates**')
            for ip in exp.get_index_pids():
                cb_name = f'{str(widget_key)}_plate_checkbox_index_{ip}'
                val = checklist_cols[i].checkbox(util.unguard_pbc(ip, silent=True),
                        key=cb_name, value=default_value)
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
            plate_type = key_parts[-2]  # e.g. 'dna', 'pcr', 'taqwater1', etc.
            plate_id = key_parts[-1]  # e.g. 'plate123'
            if plate_type in selected_pids:
                selected_pids[plate_type].append(plate_id)
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

    m(title, level='display', dest=('css',), size='h5', color=colour, align='left')


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
                    data=open(fp, 'rt'), file_name=fp_name,
                    mime='text/csv', key='dnld_samplesheet_'+str(fp),
                    type='secondary')


def show_results_table(hdr, well_data, rtype, filter_dict, filter_param_fn, key):
    """
    build an interactive table for amplicon results
    Args:
        hdr (list): header for the amplicon data
        well_data (list of lists): table of data, one row per well
        rtype (str): type of results, in ['rodentity','custom','other','amplicon']
        filter_dict (dict): dictionary of existing filter parameters
        filter_param_fn (str): file path to the filter parameters file
        key (str): key (unique id) for the view component
    Returns:
        list of chosen variants
    """
    # hide these columns but keep them available for the user
    hidden_cols = ['sex','sampleNo','sampleName','strain','clientName','alleleSymbol',
                'alleleKey','assayKey','assays','assayFamilies','primerPlate','primerWell',
                'pcrPlate','pcrWell','index_plate','i7bc','i7well','i7name','i5bc','i5well',
                'i5name','cleanCount','dnaPlate','dnaWell']
    if not well_data:
        st.info('All wells reported')
        return variant_info

    well_dataframe = pd.DataFrame(well_data, columns=hdr)
    well_table = aggrid_interactive_table(well_dataframe, key=key+'_key',
            hidden=hidden_cols, editable=['filtProportion'])
    if well_table and 'selected_rows' in well_table and well_table['selected_rows'] is not None:
        selected_rows = well_table['selected_rows'].to_dict(orient='records')
        st.session_state[f'{rtype}_chosen_vars'] = get_variants_from_table(selected_rows, filter_dict, filter_param_fn)


def get_variants_from_table(well_table, filter_dict, filter_param_fn):
    """
    Get the chosen variant from the interactive table
    args:
        well_table is a list of dicts, one per selected row
        filter_dict (dict): dictionary of existing filter parameters
        filter_param_fn (str): file path to the filter parameters file
    """
    exp = st.session_state['experiment']
    variant_info = []
    first_variant = well_table[0]
    for field in ['samplePlate','sampleWell','sampleBarcode','primer','mergeCount','seqCount','otherCount','otherName','filtProportion']:
        if field in first_variant:
            variant_info.append(first_variant[field])
        else:
            variant_info.append('')
    samplePlate = variant_info[0]
    sampleWell = variant_info[1]
    sampleBarcode = variant_info[2]
    filtProportion = variant_info[-1]
    new_id = f'{samplePlate}\t{sampleWell}\t{sampleBarcode}'
    if new_id not in filter_dict or filter_dict[new_id] != filtProportion:
        filter_dict[new_id] = filtProportion
        with open(filter_param_fn, 'wt') as f:
            for id in filter_dict:
                sp, sw, sb = id.split('\t')
                f.write(f'{sp}\t{sw}\t{sb}\t{filter_dict[id]}\n')
            if new_id not in filter_dict:
                f.write(f'{samplePlate}\t{sampleWell}\t{sampleBarcode}\t{filtProportion}\n')
    return variant_info


def make_summary_line_and_format(var_seqs, wt_seq):
    """
    Make a summary line for the chosen variant as well as a collection of formatting instructions
    """
    #print(f'make_summary_line_and_format {var_seqs=} {wt_seq=}', file=sys.stderr)
    summary_line = []
    formatting = []
    seq_array = list(var_seqs.values())
    for i in range(len(seq_array[0])):
        char_set = set([seq_array[j][i] for j in range(len(seq_array)) if seq_array[j][i] == '-' or seq_array[j][i] != wt_seq[i]])
        indel = False
        if '-' in char_set:
            indel = True
            char_set.remove('-')
        if len(char_set) == 2:
            if 'A' and 'G' in char_set:
                summary_line.append('R')
            elif 'C' and 'T' in char_set:
                summary_line.append('Y')
            elif 'G' and 'C' in char_set:
                summary_line.append('S')
            elif 'A' and 'T' in char_set:
                summary_line.append('W')
            elif 'G' and 'T' in char_set:
                summary_line.append('K')
            elif 'A' and 'C' in char_set:
                summary_line.append('M')
        elif len(char_set) == 3:
            if 'C' in char_set and 'G' in char_set and 'T' in char_set:
                summary_line.append('B')
            elif 'A' in char_set and 'G' in char_set and 'T' in char_set:
                summary_line.append('D')
            elif 'A' in char_set and 'C' in char_set and 'T' in char_set:
                summary_line.append('H')
            elif 'A' in char_set and 'C' in char_set and 'G' in char_set:
                summary_line.append('V')
        elif len(char_set) == 4:
            summary_line.append('N')
        elif len(char_set) == 1:
            summary_line.append(char_set.pop())
        elif indel:
            summary_line.append('-')
        else:
            summary_line.append('*')

    return ''.join(summary_line), formatting


def show_alignments(chosen_vars, tmp_fn, pdf_fn, target_fn, key):
    """
    Show the alignments for the chosen variants
    Run alignments against the reference sequence for the chosen primer, found in amplicon_targets.fa
    Show whether the amplicon has been reported or not using amplicons_reported.csv
    """
    exp = st.session_state['experiment']
    amplicons_reported = set()  # plate\twell\tbarcode\tprimer
    if not chosen_vars:
        st.empty()
        return None
    primer_chosen = chosen_vars[3]
    try:
        merged_counts = int(chosen_vars[4])
    except ValueError:
        merged_counts = 0
    try:
        wt_counts = int(chosen_vars[5])
    except ValueError:
        wt_counts = 0
    try:
        min_prop = float(chosen_vars[-1])
    except ValueError:
        min_prop = 0.15
    ref_seqs = {}
    with open(exp.get_exp_fn(target_fn), 'rt') as rfn:
        for line in rfn:
            if line.startswith('>'):
                ref_name = line[1:].strip()
                ref_seqs[ref_name] = ''
            else:
                ref_seqs[ref_name] += line.strip()
    ref_chosen = None
    if primer_chosen in ref_seqs:
        ref_chosen = primer_chosen
    else:
        for r in ref_seqs:
            if primer_chosen in r and 'wt' in r.lower():
                ref_chosen = r
                break
    if not ref_chosen:
        return None
    
    #print(f'{chosen_vars=} {ref_chosen=}', file=sys.stderr)
    var_list = []
    var_list_counts = {}
    for cv,cc in zip(chosen_vars[-2].split(';'),chosen_vars[-3].split(';')):
        if '//' in cv and ref_chosen in cv:
            if int(cc) >= (merged_counts * min_prop):
                var_list.append(cv.split('//')[1])
                var_list_counts[cv.split('//')[1]] = int(cc)

    print(var_list, file=sys.stderr)
    if not var_list:
        #st.warning('No variant sequences found for this amplicon')
        return None

    #vars_chosen = {cv.split('//')[1]:int(cc) for cv,cc in zip(chosen_vars[-2].split(';'),chosen_vars[-3].split(';')) if ref_chosen in cv and '//' in cv and int(cc)>=(merged_counts*min_prop)}
    id = f'{chosen_vars[0]}\t{chosen_vars[1]}\t{chosen_vars[2]}\t{chosen_vars[3]}'
    #print(f'show_alignments() {ref_chosen=} {ref_seqs[ref_chosen]=}, {var_list=}', file=sys.stderr)
    ref_seq_id = (ref_chosen, ref_seqs[ref_chosen])
    aligned = run_msa((ref_chosen,ref_seqs[ref_chosen]),var_list)
    var_df = pd.DataFrame({
        'Report': [True for name in aligned.names],
        'Variant': [name for name in aligned.names],
        'Counts': [var_list_counts[name] if name != ref_chosen else wt_counts for name in aligned.names],
        'Sequence': [str(aligned.get_gapped_seq(name)) for name in aligned.names]
    })
    # st.markdown(
    #     """
    #     <style>
    #     /* Apply monospace to all text within the main content area */
    #     [class^=dvn-scroller] {
    #         font-family: Courier New", Courier, monospace;
    #     }
    #     </style>
    #     """,
    #     unsafe_allow_html=True
    # )

    # Apply the style to the DataFrame
    var_de = st.data_editor(var_df, key='alignment_editor_'+key, hide_index=True,
        column_config={
            'Report': st.column_config.CheckboxColumn('Report', width=40,
                    help='Check to include this variant in the PDF report', default=True),
            'Variant': st.column_config.TextColumn('Variant', width=150, disabled=True,
                    help='Name of the variant sequence'),
            'Counts': st.column_config.NumberColumn('Counts', width=80, disabled=True,
                    help='Number of reads supporting this variant'),
            'Sequence': st.column_config.TextColumn('Aligned sequence', width=1000, disabled=True,
                    help='The aligned sequence of this variant')
        },
        disabled=['Variant','Counts','Sequence'])

    prefix = st.text_input(label='Enter any descriptive notes you wish to prefix to this alignment')
    make_report_entry = st.button('Make report entry', key='build_pdfs_'+chosen_vars[0]+'_'+chosen_vars[1]+'_'+chosen_vars[2])
    if make_report_entry:
        if var_de is not None:
            aligned_stuff = [(n,f'{s}') for n,s,c in zip(var_de['Variant'], var_de['Sequence'], var_de['Report']) if c]
        #valid_keys = [cb.split('_seq_checkbox')[0] for cb in checkbox_keys if cb in st.session_state and st.session_state[cb]]
        #if valid_keys:
            #if not prefix:
            #    prefix = ''
            #aligned_stuff = {n:f'{s}' for n,s in aligned_names_seqs.items() if n in valid_keys}
            #aligned_stuff = [(n,f'{s}') for n,s in aligned_names_seqs.items() if n in valid_keys]
            if aligned_stuff:
                formatting = {'font':'Courier New', 'fontsize':10, 'lineheight':1.2, 'colour_changes':True}
                formatting = [(k,v) for k,v in formatting.items()]
                page = (prefix, tuple(aligned_stuff), tuple(formatting))
                if Path(tmp_fn).exists():
                    with open(tmp_fn, 'rt') as f:
                        pages = jsonpickle.decode(f.read(), keys=True, handle_readonly=True)
                        #print(type(pages), pages)
                        pages = set(pages)
                        pages.add(page)
                else:
                    pages = set()
                    pages.add(page)

                json_pages = jsonpickle.encode(list(pages), indent=4, keys=True, warn=True, handle_readonly=True)
                with open(tmp_fn, 'wt') as fout:
                    print(json_pages, file=fout)

                generate_pdf(pdf_fn, pages)
                st.success(f'Report entry made in {pdf_fn}')
                st.success(f'Adding {id} to reported amplicons list')
                return id
    
    return None


def build_paired_sequence_views(ref_seq, var_anno, display_width=120, colour_all=False, colour_changes=True):
    """
    Build the views for the reference sequence, variant sequence, link and summary
    Args:
        ref_seq (str): reference sequence
        var_anno (str): the sequence annotation chosen by the user
        display_width (int): width of the display in characters, wrap outputs to this width
        colour_all (bool): False, whether to colour every position in the display
        colour_changes (bool): True, whether to colour the changes of variable sites only
    Returns:
        outputs (list of tuples): each tuple contains four strings:
            ref_view, link_view, var_view, summary_view (str): the views for the reference sequences,
            the link between the reference and variable sequences,
            the variable sequence, and a summary of the variable sequence
    Notes:
        This function builds on the variant sequence generated in ngsmatch.get_variant_seq()
        It also adjusts the reference sequence so that bases remain in matching positions
    """
    if ref_seq is None or var_anno is None:
        return [('', '', '', '')]

    var_seq, mod_seq = build_aligned_pair(var_anno, ref_seq)

    # link view shows direct matches between sequences
    link_view = ['|' if c1 == c2 else ' ' for c1, c2 in zip(mod_seq, var_seq)]
    summary_view = []  # summary of differences
    for m,v in zip(mod_seq, var_seq):
        if m == v:
            summary_view.append('*')
        elif m == '-':
            summary_view.append('+')
        elif v == '-':
            summary_view.append('-')
        elif v != m:
            summary_view.append(v)
        else:
            summary_view.append(' ')
    summary_view = ''.join(summary_view)
    outputs = []
    link_chrs = ''.join(link_view)
    summary_chrs = ''.join(summary_view)
    for i in range(0, len(var_seq), display_width):
        outputs.append((mod_seq[i:i+display_width],
                link_chrs[i:i+display_width],
                var_seq[i:i+display_width],
                summary_chrs[i:i+display_width]))
    return outputs


def show_results_display(results_type, key, caller_id=None):
    """
    Show the results display for the chosen variants
    Args:
        hdr (list): header for the amplicon data
        amplicon_data (list of lists): the amplicon data
        filter_dict (dict): dictionary of per-amplicon filter parameters
        key (str): key (unique id) for the view component
    """
    exp = st.session_state['experiment']
    if results_type not in ['rodentity','custom','other','amplicon']:
        m(f'Unknown results type {results_type} in show_results_display()', caller_id=caller_id, level='critical')
        return
    rtype = results_type
    perform_reset = False
    pdf_fn = exp.get_exp_fn(f'{rtype}_alignments.pdf')
    tmp_fn = exp.get_exp_fn(f'{rtype}_alignments.json')
    reported_fn = exp.get_exp_fn(f'{rtype}_reported.txt')
    filter_param_fn = exp.get_exp_fn(f'{rtype}_table_filter_params.txt')
    if results_type == 'amplicon':
        results_fn = exp.get_exp_fn('amplicon_results.csv')
        target_fn = 'amplicon_targets.fa'
    else:
        results_fn = exp.get_exp_fn('results.csv')
        target_fn = 'targets.fa'
    # if rtype == 'amplicon':
    #     hdr, data, filter_dict, data_reported, data_unreported = exp.gather_amplicon_results(reported_amps_fn)
    # else:
    hdr, data, filter_dict, data_reported, data_unreported = exp.gather_results(rtype, reported_fn, results_fn)
    unreported_container = st.container(key=f'{rtype}_unreported_container')
    alignment_container = st.container(key=f'{rtype}_alignment_container')
    reported_container = st.container(key=f'{rtype}_reported_container')

    with unreported_container:
        if data_unreported:
            st.write(f'Unreported {rtype}:')
            show_results_table(hdr, data_unreported, rtype, filter_dict, filter_param_fn, f'unreported_{rtype}')
    with reported_container:
        if data_reported:
            st.write(f'Reported {rtype}:')
            show_results_table(hdr, data_reported, rtype, filter_dict, filter_param_fn, f'reported_{rtype}')
            wipe_report = st.button(f'Reset {rtype} report', key=f'reset_{rtype}_report_button')
            if wipe_report:
                success = True
                if Path(reported_fn).exists():
                    try:
                        Path(reported_fn).unlink()
                    except Exception as exc:
                        st.failure(f'Could not delete {reported_fn}, perhaps you have it open? {exc}')
                        success = False
                if Path(tmp_fn).exists():
                    try:
                        Path(tmp_fn).unlink()
                    except Exception as exc:
                        st.failure(f'Could not delete {tmp_fn}, perhaps you have it open {exc}')
                        success = False
                if Path(pdf_fn).exists():
                    try:
                        Path(pdf_fn).unlink()
                    except Exception as exc:
                        st.failure(f'Could not delete {pdf_fn}, perhaps you have it open? {exc}')
                        success = False
                if success:
                    st.success(f'{rtype} report reset to empty')
                    perform_reset = True

    with alignment_container:
        if f'{rtype}_chosen_vars' in st.session_state and st.session_state[f'{rtype}_chosen_vars']:
            #well_variants = get_variants_from_table(st.session_state[f'{rtype}_chosen_vars'], filter_dict, filter_param_fn)
            #print('Ready to do alignments', st.session_state[f'{rtype}_chosen_vars'], file=sys.stderr)
            with st.spinner('Generating alignments...'):
                st.session_state[f'{rtype}_chosen_id'] = show_alignments(st.session_state[f'{rtype}_chosen_vars'], tmp_fn, pdf_fn, target_fn, f'{rtype}_alignment_viewer')
        else:
            with st.spinner('Generating alignments...'):
                st.session_state[f'{rtype}_chosen_id'] = show_alignments(None, tmp_fn, pdf_fn, target_fn, f'{rtype}_alignment_viewer')
        if f'{rtype}_chosen_id' in st.session_state and st.session_state[f'{rtype}_chosen_id']:
            entries_reported = set()
            if Path(reported_fn).exists():
                with open(reported_fn, 'rt') as f:
                    for line in f:
                        entries_reported.add(line.strip())
            entries_reported.add(st.session_state[f'{rtype}_chosen_id'])
            with open(reported_fn, 'wt') as fout:
                for er in entries_reported:
                    print(er, file=fout)
            perform_reset = True
        else:
            st.empty()



    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()

    if perform_reset:
        st.session_state[f'{rtype}_chosen_vars'] = None
        st.session_state[f'{rtype}_chosen_id'] = None
        sleep(0.4)
        st.rerun()



def add_css():
    # Disabled since we no longer need monospaced font everywhere
    # with open( "style.css" ) as css:
    #     st.markdown( f'<style>{css.read()}</style>' , unsafe_allow_html= True)
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



