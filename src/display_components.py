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

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import streamlit as st
import extra_streamlit_components as stx
from st_aggrid import AgGrid, GridOptionsBuilder, JsCode
from st_aggrid.shared import GridUpdateMode

from stutil import custom_text, add_vertical_space, hline, m, init_state, mq
from match import build_aligned_pair, reconstruct_sequence, run_msa
from generate import generate_pdf
from experiment import Experiment, EXP_FN, load_experiment
import util
import db_io


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


def set_state(key, value):
    """ Callback function for display elements """
    st.session_state[key] = value


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






