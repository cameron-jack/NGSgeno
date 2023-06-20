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

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components
import extra_streamlit_components as stx
from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode
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
                delbox.write(f'{id} removed')
            else:
                delbox.write(f'{id} could not be removed')
            st.session_state['previous_file_delete_selection'] = ids
    elif category == 'plate':
        gids = [util.guard_pbc(pid, silent=True) for pid in ids if pid in exp.plate_location_sample]
        success = exp.delete_plates(gids)
        if success:
            delbox.write(f'{ids} removed')
        else:
            delbox.write(f'{ids} could not be removed')
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
    """ display a summary of all loaded DNA and sample plates """
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
        d = {'taqwater_pids':[], 'taq_vol':0, 'water_vol':0, 'primer_pids':[], 'primer_count_ngs':0, 
                'primer_count_custom':0, 'unique_primers':set(), 'primer_well_count':0, 'assay_primer_mappings':0,
                'reference_files':[], 'unique_references':set(), 'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
    """
    exp = st.session_state['experiment']
    # display all the required files whether they are, or are not present
    # Plates: primer, index, taq/water, PCR, references, primer/assaylist, 
    headers = ['Purpose','Barcode/ID','Wells','Type','Entries']
    consumables = exp.summarise_consumables()

    data_rows = []
    for tp in consumables['taqwater_pids']:
        data_rows.append(['taq/water', util.unguard_pbc(tp, silent=True), 6, util.PLATE_TYPES['Echo6'], 6])
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
        selection = aggrid_interactive_table(plate_df, grid_height=height, key='plate_aggrid')


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

        selection = aggrid_interactive_table(plate_df, grid_height=height, key='plate_aggrid')
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
                    delbox = st.container()
                    del_col1, del_col2, del_col3, _ = delbox.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    del_col2.button("Yes",on_click=manage_delete,
                            args=(delbox,'plate',pids), key="delete " + str(key), help=f"Delete {pids}")
                    del_col3.button("No", on_click=cancel_delete,
                            args=('plate',pids), key="keep " + str(key), help=f"Keep {pids}")
    

def display_pcr_components(PCR_stage=1, dna_pids=None, show_general=True):
    """
    Expander widget that shows the required componenents for each PCR reaction, 
    including wells, PCR plates, taq+water plates, index pairs 
    Args:
        assay_usage: from Experiment.get_assay_usage
        PCR_stage (1, 2): 1 = Echo Primer stage, 2 = Echo Indexing
        show_general (bool): if True show all of the general information for PCR reactions,
        including user supplied plates.  
    
    """
    exp = st.session_state['experiment']
    DNA_PLATE_WELLS = 384
    ul_conv = 1000
    pcr_comps_area = st.container()
    
    _,filter_col1,_,filter_col2, _ = pcr_comps_area.columns([1,5,1,5,1])
    col_size = [6, 4, 6, 4]

    fwd_idx, rev_idx = exp.get_index_avail()
    print(f'{dna_pids=}')
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids=dna_pids)
    num_reactions = sum([primer_usage[p] for p in primer_usage])
    # print(f'{num_reactions=}')
    # print(f'{primer_usage=}')
    primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            exp.get_taq_water_volumes_required(num_reactions)
    num_req_pcr_taq_water_plates = util.num_req_taq_water_plates(primer_taq_vol, primer_water_vol)
    num_req_index_taq_water_plates = util.num_req_taq_water_plates(index_taq_vol, index_water_vol)
    num_req_total_taq_water_plates = util.num_req_taq_water_plates(primer_taq_vol+index_taq_vol, 
            primer_water_vol+index_water_vol)
        
    index_remain, index_max = exp.get_index_reactions(primer_usage, fwd_idx, rev_idx)
        
    taq_avail, water_avail, pids = exp.get_taqwater_avail()
    
    if show_general:
        required_PCR_plates = ceil(num_reactions/DNA_PLATE_WELLS)
        user_supplied_PCR = ', '.join([util.unguard_pbc(p, silent=True)\
                    for p in exp.get_pcr_pids()])
        user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True)\
                    for p in exp.get_taqwater_avail()[2]])

        comp_warning_area = pcr_comps_area.container()

        num_supplied_PCR = len(exp.get_pcr_pids())
        num_supplied_taqwater = len(user_supplied_taqwater)
        if num_supplied_PCR < required_PCR_plates:
            comp_warning_area.warning(f'Please add PCR plates. {required_PCR_plates-num_supplied_PCR} PCR plates are required')
        if num_supplied_taqwater < num_req_total_taq_water_plates:
            comp_warning_area.warning(f'Please add taq/water plates. {num_req_total_taq_water_plates} taq/water plates required')

        req_cols = pcr_comps_area.columns(col_size)
        #Wells and plates
        if num_supplied_PCR < required_PCR_plates:
            req_cols[0].markdown('<p style="color:#FF0000"><b>Number of required PCR plates</b></p>', unsafe_allow_html=True)
            req_cols[1].markdown(f'<p style="color:#FF0000">{str(required_PCR_plates)}</p>', unsafe_allow_html=True)
        else:
            req_cols[0].markdown('**Number of required PCR plates**')
            req_cols[1].write(str(required_PCR_plates))
        req_cols[0].markdown('**Worst case required reaction wells**')
        req_cols[1].write(str(num_reactions))
        req_cols[2].markdown('**User supplied PCR plates**')
        if user_supplied_PCR:
            req_cols[3].write(user_supplied_PCR)
        else:
            req_cols[3].markdown('<p style="color:#FF0000">None</p>', unsafe_allow_html=True)
        req_cols[2].markdown('**User supplied taq/water plates**')
        if user_supplied_taqwater:
            req_cols[3].write(user_supplied_taqwater)
        else:
            req_cols[3].markdown('<p style="color:#FF0000">None</p>', unsafe_allow_html=True)
        for i in range(4):
            req_cols[i].write('')

    #PCR1 and PCR2 taq and water 
    taq_avail_vol = taq_avail/ul_conv
    water_avail_vol = water_avail/ul_conv                  

    pcr_cols = pcr_comps_area.columns(col_size)

    if PCR_stage == 1:
        #get actual values for volume of taq water plates
        required_water_vol_str = str(primer_water_vol/ul_conv) + ' μl'
        water_avail_vol_str = str(water_avail_vol)+' μl'
        required_taq_vol_str = str(primer_taq_vol/ul_conv) + ' μl'
        avail_taq_vol_str = str(taq_avail_vol)+' μl'

    if PCR_stage == 2:
        required_water_vol_str = str(index_water_vol/ul_conv)+ ' μl'
        water_avail_vol_str = str(water_avail_vol) + ' μl'
        required_taq_vol_str = str(index_taq_vol/ul_conv)+ ' μl'
        avail_taq_vol_str = str(taq_avail_vol)+ ' μl'

        index_cols = pcr_comps_area.columns(col_size)
        if index_remain >= 0:
            index_pairs_allowed = '<p style="color:green">'+str(index_remain)+'</p>'
        else:
            index_pairs_allowed = '<p style="color:red">'+str(index_remain)+'</p>'

        index_cols[0].markdown('**Index Pairs Remaining**')
        index_cols[1].write(str(index_remain))
        index_cols[2].markdown('**Max Possible Index Pairs**')
        index_cols[3].write(str(index_max))
        index_cols[0].markdown('**Index Pairs Allowed by Volume**')
        index_cols[1].markdown(index_pairs_allowed, unsafe_allow_html=True)

    if num_supplied_taqwater < num_req_total_taq_water_plates:
        req_cols[0].markdown('<p style="color:#FF0000"><b>Number of required taq/water plates</b></p>', unsafe_allow_html=True)
        req_cols[1].markdown(f'<p style="color:#FF0000">{str(num_req_total_taq_water_plates)}</p>', unsafe_allow_html=True)
    else:
        req_cols[0].markdown('**Number of required taq/water plates**')
        req_cols[1].write(str(num_req_total_taq_water_plates))
    pcr_cols[0].markdown('**Required water volume**')
    pcr_cols[1].write(required_water_vol_str)
    pcr_cols[2].markdown('**Available water volume**')
    pcr_cols[3].write(water_avail_vol_str)
    pcr_cols[0].markdown('**Required taq volume**')
    pcr_cols[1].markdown(required_taq_vol_str)
    pcr_cols[2].markdown('**Available taq volume**')
    pcr_cols[3].write(avail_taq_vol_str)


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

            with open("debug.html", 'wt') as outf:
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
    selection = aggrid_interactive_table(file_df, grid_height=height, key='file_aggrid')
    if selection and 'selected_rows' in selection:
        fns = [row['File'] for row in selection['selected_rows']]
        fns = [fn for fn in fns if fn in exp.uploaded_files]
        if fns:
            print(f'{fns=}')
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
                value=350, step=25, help="Size of display grid", key=key)
    with col1:
        view_tab = stx.tab_bar(data=[
            stx.TabBarItemData(id=1, title="Status", description=""),
            stx.TabBarItemData(id=2, title="Files", description=""),
            stx.TabBarItemData(id=3, title="Plates", description=""),
            stx.TabBarItemData(id=4, title="Plate Viewer", description=""),
            stx.TabBarItemData(id=5, title="Primers", description=""),
            stx.TabBarItemData(id=6, title="Indexes", description=""),
            stx.TabBarItemData(id=7, title="Reference sequences", description=""),
            stx.TabBarItemData(id=8, title="Log", description="")
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


    if view_tab == 7:
        with container:
            display_references(key, height=view_height)


    if view_tab == 8:
        with container:
            display_log(key, height=view_height)
    




