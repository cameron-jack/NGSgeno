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
try:
    from bin.experiment import Experiment, EXP_FN, load_experiment
except ModuleNotFoundError:
    from experiment import Experiment, EXP_FN, load_experiment
try:
    import bin.util as util
except ModuleNoteFoundError:
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
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        height=grid_height,
        gridOptions=options.build(),
        # _available_themes = ['streamlit','light','dark', 'blue', 'fresh','material']
        theme="fresh",
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
                if 'PID' in key and entry[key]:
                    if not util.is_guarded_pbc(entry[key]):
                        entries[i][key] = util.guard_pbc(entry[key])
        exp.delete_plates(entries)
        exp.save()
    #print('delete_entries ran')


def data_table(key, options=False, table_option=None):
    """
    Shows the loaded data in a table.
    Args:
        key (Any): key for the table so that it is a unique widget, if there's more than one in a page
        options (bool): Shows different options for table view in a selectbox if True
        table_option (str): If options is False, can choose a table_option instead ('Summary', 'Edit Table', 
        'Plate View', 'Explore Assays', 'View Log)

    Need to fix the delete selection part - clicking no should reset the table view. 
    """

    data_container = st.container()
    exp = st.session_state['experiment']
    if options==True:
        table_option = data_container.selectbox('View Data options',
                ('Loaded Samples', 'Loaded Consumables', 'Edit Table', 'Plate View', 'Explore Assays', 
                'Reference Sequences', 'View Log'), key='select_' + str(key))
    elif options==False and table_option is None:
        table_option = 'Loaded Samples'

    data_table_title = data_container.title('')
    data_table_title.markdown(f'<h5 style="text-align:center;color:#7eaecc">{table_option}</h5>',\
            unsafe_allow_html=True)

    #exp = st.session_state['experiment']

    #Experiment summary - table

    if table_option == 'Loaded Samples':
        if exp:
            df = exp.inputs_as_dataframe()
            if df is None or not isinstance(df, pd.DataFrame):
                st.write('No data loaded')
            else:
                #print(f"{type(df)=}, {df=})")
                selection = aggrid_interactive_table(df, key=key)
                #print(f"{type(selection)=} {selection=}")
                if 'selected_rows' in selection and selection['selected_rows']:
                    # only do the code below if this is a fresh selection
                    rows = selection["selected_rows"]
                    lines = '\n'.join(['DNA PID: '+r['DNA PID'] for r in rows if r['DNA PID'] != 'Total'])
                    if 'previous_delete_group_selection' not in st.session_state:
                        st.session_state['previous_delete_group_selection'] = None
                    if lines and st.session_state['previous_delete_group_selection'] != lines:
                        data_container.markdown(f"**You selected {lines}**")
                        del_col1, del_col2, del_col3, _ = data_container.columns([2,1,1,4])
                        del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                        delete_button = del_col2.button("Yes",on_click=delete_entries, 
                                args=(exp, selection['selected_rows']), key="delete " + str(key), help=f"Delete {lines}")
                        mistake_button = del_col3.button("No",on_click=delete_entries, args=(exp, []), 
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

    elif table_option == 'Loaded Consumables':
        # d = {'taqwater_pids':[], 'taq_vol':0, 'water_vol':0, 'primer_pids':[], 'primer_count_ngs':0, 
        #        'primer_count_custom':0, 'unique_primers':set(), 'primer_well_count':0, 'assay_primer_mappings':0,
        #        'reference_files':[], 'unique_references':set(), 'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
        consumables = exp.summarise_consumables()
        # display pids and files
        list_col1, list_col2 = data_container.columns([1,4])
        list_col1.markdown('**Taq+water plates**')
        list_col2.write(', '.join([util.unguard_pbc(pid, silent=True) for pid in consumables['taqwater_pids']]))
        list_col1.markdown('**Primer plates**')
        list_col2.write(', '.join([util.unguard_pbc(pid, silent=True) for pid in consumables['primer_pids']]))
        list_col1.markdown('**Index plates**')
        list_col2.write(', '.join([util.unguard_pbc(pid, silent=True) for pid in consumables['index_pids']]))
        list_col1.markdown('**Reference files**')
        list_col2.write(', '.join(consumables['reference_files']))
        numeric_cols = data_container.columns([2,1,2,1,2,1])
        numeric_cols[0].markdown('**Available taq (uL)**')
        numeric_cols[1].write(consumables['taq_vol']/1000)
        numeric_cols[2].markdown('**Available water (uL)**')
        numeric_cols[3].write(consumables['water_vol']/1000)
        numeric_cols[4].markdown('**Assay to primer mappings**')
        numeric_cols[5].write(consumables['assay_primer_mappings'])
        numeric_cols[0].markdown('**Unique primers**')
        numeric_cols[1].write(len(consumables['unique_primers']))
        numeric_cols[2].markdown('**Unique forward indexes**')
        numeric_cols[3].write(len(consumables['unique_i7s']))
        numeric_cols[4].markdown('**Unique reverse indexes**')
        numeric_cols[5].write(len(consumables['unique_i5s']))
        numeric_cols[0].markdown('**Unique reference sequences**')
        numeric_cols[1].write(len(consumables['unique_references']))

    elif table_option == 'Edit Table':
        # Show plates in table form and let the user edit them to fix minor mistakes
        pass

    #Show a visual of the plates
    elif table_option == 'Plate View':
        plate_ids = []
        for pid in exp.plate_location_sample:
            plate_ids.append(util.unguard(pid))

        #Let user choose which plate to view
        plate_selectbox = data_container.selectbox('Plate ID to view', plate_ids, key=key)
        if plate_selectbox:
            plate_id = util.guard_pbc(plate_selectbox)
            if plate_id in exp.plate_location_sample:
                heatmap_str = generate_heatmap_html(exp, plate_id, scaling=1.2)
                #with open("debug.html", 'wt') as outf:
                #    print(heatmap_str, file=outf)
                components.html(heatmap_str, height=700, scrolling=True)
            else:
                plate_barcode_error_msg = "Plate barcode not found in experiment"
                data_container.markdown(f'<p style="color:#FF0000">{plate_barcode_error_msg}</p>',\
                unsafe_allow_html=True)
        # let user choose plates and then display them as a plate layout with well contents on hover
        # plate_col1, plate_col2 = st.columns(2)
        # for plate in st.session_state['experiment'].plate_location_sample:
        #    break
        # heatmap_str = generate_heatmap_html(plate_id, st.session_state['experiment'], scaling=1)
        # components.html(heatmap_str, height=600, scrolling=True)
        

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

    elif table_option == 'Reference Sequences':
        #print(f"{st.session_state['experiment'].reference_sequences=}", file=sys.stderr)
        refseq = st.session_state['experiment'].reference_sequences
        with data_container:
            for group in refseq:
                dataset = [[id,refseq[group][id]] for id in refseq[group]]
                ref_df = pd.DataFrame(dataset, columns=['Name','Sequence'])
                st.markdown('**'+group+'**')
                aggrid_interactive_table(ref_df, key='seqs'+group)
                
    elif table_option == 'View Log':
        func = sys._getframe(1).f_code.co_name
        func_line = inspect.getframeinfo(sys._getframe(1)).lineno
        print("Function", func, type(func_line))
        # display the experiment log and let the user filter the view in a number of ways
        log_entries = st.session_state['experiment'].get_log(100)
        if len(log_entries) == 0:
            st.write('No entries currently in the log')
        else:
            df = pd.DataFrame(log_entries, columns=st.session_state['experiment'].get_log_header())
            #print(df)
            if 'view_box_size' in st.session_state:
                height = st.session_state['view_box_size']
            else:
                height = 250
            aggrid_interactive_table(df, key='logs')
            #data_container.dataframe(df, height=height)


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
        user_supplied_PCR = ', '.join([util.unguard_pbc(p, silent=True)\
                    for p in st.session_state['experiment'].get_pcr_pids()])
        user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True)\
                    for p in st.session_state['experiment'].get_taqwater_avail()[2]])

        comp_warning_area = pcr_comps_area.empty()
        comp_warning = ''

        req_cols = pcr_comps_area.columns(col_size)
        #Wells and plates
        req_cols[0].markdown('**Number of required PCR plates**')
        req_cols[1].write(str(required_wells))
        req_cols[0].markdown('**Number of required reaction wells**')
        req_cols[1].write(str(reactions))
        req_cols[2].markdown('**User supplied PCR plates**')
        if user_supplied_PCR:
            req_cols[3].write(user_supplied_PCR)
        else:
            req_cols[3].markdown('<p style="color:#FF0000">None</p>', unsafe_allow_html=True)
        req_cols[2].markdown('**User supplied taq/water plates**', unsafe_allow_html=True)
        if user_supplied_PCR:
            req_cols[3].write(user_supplied_taqwater)
        else:
            req_cols[3].markdown('<p style="color:#FF0000">None</p>', unsafe_allow_html=True)
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
    exp  = st.session_state['experiment']
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
                exp.get_volumes_required(assay_usage=assay_usage)
    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            0,{},0,0,0,0 

    primer_avail_counts, primer_avail_vols = exp.get_primers_avail()
    usable_well_vol = util.CAP_VOLS['384PP_AQ_BP'] - util.DEAD_VOLS['384PP_AQ_BP']
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
            num_wells = ceil((need_vol)/usable_well_vol)

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


def show_echo1_outputs():
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)
    dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths = exp.get_echo_PCR1_picklist_filepaths()
    error_msgs = []
    
    if not dna_picklist_paths:
        error_msgs.append('No DNA picklist available')
    else:
        for dpp in dna_picklist_paths:
            dpp_file = str(Path(dpp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{dpp_file}</p>',\
                    unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(dpp, 'rt'), file_name=dpp, mime='text/csv', key='dna_download_'+dpp)

    if not primer_picklist_paths:
        error_msgs.append('No primer picklist available')
    else:
        for ppp in primer_picklist_paths:
            ppp_file = str(Path(ppp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ppp_file}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(ppp, 'rt'), file_name=ppp, mime='text/csv', key='primer_download_'+dpp)
            
    if not taqwater_picklist_paths:
        error_msgs.append('No taq/water picklist available')
    else:
        for tpp in taqwater_picklist_paths:
            tpp_file = str(Path(tpp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_file}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(tpp, 'rt'), file_name=tpp, mime='text/csv', key='taqwater_download_'+tpp)

    if error_msgs:
        for msg in error_msgs:
            picklist_file_col.markdown(f'<p style="color:#ff0000;text-align:right">{msg}</p>',\
                    unsafe_allow_html=True)

def show_echo2_outputs():
    exp = st.session_state['experiment']
    picklist_file_col, picklist_btn_col = st.columns(2)
    index_picklist_paths, taqwater_picklist_paths, amplicon_picklist_paths = exp.get_echo_PCR2_picklist_filepaths()
    
    error_msgs = []

    # amplicon_picklist_paths is a list
    if not index_picklist_paths:
        error_msgs.append('No index picklist available')
    else:
        for ipp in index_picklist_paths:
            ipp_file = str(Path(ipp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ipp_file}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                    data=open(ipp, 'rt'), file_name=ipp, mime='text/csv',\
                            key='index_download_'+ipp)

    if not amplicon_picklist_paths:
        error_msgs.append('No amplicon picklist available')               
    else:
        for app in amplicon_picklist_paths:
            app_file = str(Path(app).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{app_file}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                    data=open(app, 'rt'), file_name=app, mime='text/csv',\
                            key='amplicon_download_'+app)

    if not taqwater_picklist_paths:
        error_msgs.append('No taq/water picklist available')

    else:
        for tpp in taqwater_picklist_paths:
            tpp_file = str(Path(tpp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_file}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                        data=open(tpp, 'rt'), file_name=tpp, mime='text/csv',\
                                key='taqwater_download_'+tpp)
    if error_msgs:
        for msg in error_msgs:
            picklist_file_col.markdown(f'<p style="color:#ff0000;text-align:right">{msg}</p>',\
                    unsafe_allow_html=True)
    