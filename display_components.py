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


#def delete_entries(exp, entries):
#    """ removes data from an experiment, called by buttons """
#    print(f"Debug: (delete_entries) {exp.name=} {entries=}")
#    if exp is None:
#        if 'delete_entries' in st.session_state:
#            st.session_state['delete_entries'] = []
#        return
#    if entries: # list
#        pids = []
#        for i,entry in enumerate(entries):
#            for key in entry:
#                if 'PID' in key and len(entry[key].strip()) != 0:
#                    pids.append(util.guard_pbc(entry[key], silent=True))
#        print(f'Following pids selected for deletions {pids=}', file=sys.stderr)
#        exp.delete_plates(pids)
#        exp.save()
#    #print('delete_entries ran')


def data_table(key, options=False, table_option=None, height=350):
    """
    Shows the loaded data in a table.
    Args:
        key (Any): key for the table so that it is a unique widget, if there's more than one in a page
        options (bool): Shows different options for table view in a selectbox if True
        table_option (str): If options is False, can choose a table_option instead ('Summary', 'Edit Table', 
        'Plate View', 'Explore Assays', 'View Log)

    Need to fix the delete selection part - clicking no should reset the table view and clicking yes creates error
    """
    exp=st.session_state['experiment']
    data_container = st.container()
    selectbox_col,_,_,_ = data_container.columns(4)
    exp = st.session_state['experiment']
    if options==True:
        table_option = selectbox_col.selectbox('View Data options',
                                ('Loaded Samples', 'Edit Table', 'Plate View', 'Explore Assays', 'Reference Sequences', 'View Log'),\
                                             key='select_' + str(key))
    elif options==False and table_option is None:
        table_option = 'Summary'

    data_table_title = data_container.title('')
    data_table_title.markdown(f'<h5 style="text-align:center;color:#7eaecc">{table_option}</h5>',\
            unsafe_allow_html=True)

    #exp = st.session_state['experiment']

    #Experiment summary - table

    if table_option == 'Loaded Samples':
        selection = []
        st.session_state['nuke'] = ''

        if exp:
            df = exp.inputs_as_dataframe()
            if df is None or not isinstance(df, pd.DataFrame):
                st.write('No data loaded')
            else:
                #print(f"{type(df)=}, {df=})")
                selection = aggrid_interactive_table(df, key=key, grid_height=height)
                #print(f"{type(selection)=} {selection['selected_rows']=}")
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
                        delete_button = del_col2.button("Yes",on_click=exp.remove_entries, 
                                args=rows, key="delete " + str(key), help=f"Delete {lines}")
                        mistake_button = del_col3.button("No",on_click=exp.remove_entries, args=rows, 
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
                heatmap_str = generate_heatmap_html(exp, plate_id, scaling=1.5)

                with open("debug.html", 'wt') as outf:
                    print(heatmap_str, file=outf)
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

def display_plates():
    plate_array = []
    exp = st.session_state['experiment']
    for keys, values in exp.plate_location_sample.items():
        plate_array.append([util.unguard(keys), len(values['wells']), values['purpose']])
    
    plate_df = pd.DataFrame(plate_array, columns=['Plates', 'Num Wells', 'Purpose']) 
    plate_table = aggrid_interactive_table(plate_df, grid_height=350, key='plate_aggrid')

def display_pcr_components(assay_usage, PCR_stage=1, show_general=True, filtered=True):
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
    AFTER_NIM_WELLS = 384
    pcr_comps_area = st.container()
    
    _,filter_col1,_,filter_col2, _ = pcr_comps_area.columns([1,5,1,5,1])
    col_size = [6, 4, 6, 4]
    comp_warning = ''

    fwd_idx, rev_idx, warning_idxs = exp.get_index_avail()

    if assay_usage:
        primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
                exp.get_volumes_required(assay_usage=assay_usage, filtered=True)
        reactions = sum([v for v in assay_usage.values()])
        
        index_remain, index_avail, index_capacity =\
                exp.get_index_remaining_available_volume(assay_usage=assay_usage, fwd_idx=fwd_idx, rev_idx=rev_idx)
    
    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            0,{},0,0,0,0 
        index_remain, index_avail, index_vol_capacity =\
                st.session_state['experiment'].get_index_remaining_available_volume()
        
    taq_avail, water_avail, pids = st.session_state['experiment'].get_taqwater_avail()
    
    if show_general:
        #if 'assay_filter' not in st.session_state:
        #    st.session_state['assay_filter'] = 1

        #    filter_assays = filter_col1.button('Filter assays', key='assay_filter_button', disabled=True)
        #    unfilter_assays = filter_col2.button('Unfilter assays', key='assay_unfilter_button')

        #    if unfilter_assays:
        #        st.session_state['assay_filter'] = 0

        #    if filter_assays:
        #        st.session_state['assay_filter'] = 1

        required_wells = ceil(reactions/AFTER_NIM_WELLS)
        user_supplied_PCR = ', '.join([util.unguard_pbc(p, silent=True)\
                    for p in exp.get_pcr_pids()])
        user_supplied_taqwater = ', '.join([util.unguard_pbc(p, silent=True)\
                    for p in exp.get_taqwater_avail()[2]])

        comp_warning_area = pcr_comps_area.empty()

        num_supplied_PCR = len(exp.get_pcr_pids())
        if num_supplied_PCR < required_wells:
            comp_warning += f'{required_wells-num_supplied_PCR} PCR plate(s) '

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
            req_cols[3].markdown('<p style="color:#FF0000">None</p>', unsafe_allow_html=True)
        pcr_comps_area.write('')

    #PCR1 and PCR2 taq and water 
    taq_avail_vol = taq_avail/1000
    water_avail_vol = water_avail/1000

    #num_taqwater_plates = ceil(primer_taq_vol / taq_avail_vol)

    pcr_cols = pcr_comps_area.columns(col_size)

    if PCR_stage == 1:
        if taq_avail < primer_taq_vol or water_avail < primer_water_vol:
            if comp_warning != '':
                comp_warning += ' and '
            
            comp_warning += f' {ceil((primer_taq_vol/1000)/7650)} taq and water plate(s) '
    #get actual values for volume of taq water plates
        

        required_water_vol_str = str(primer_water_vol/1000) + ' μl'
        water_avail_vol_str = str(water_avail_vol)+' μl'
        required_taq_vol_str = str(primer_taq_vol/1000)+' μl'
        avail_taq_vol_str = str(taq_avail_vol)+' μl'

    if PCR_stage == 2:
        if taq_avail < (primer_taq_vol + index_taq_vol) or water_avail < (primer_water_vol + index_water_vol):
            if comp_warning != '':
                comp_warning +=  ' and'
 
            comp_warning += ' taq and water plate(s)'
        
        required_water_vol_str = str(index_water_vol/1000)+ ' μl'
        water_avail_vol_str = str(water_avail_vol) + ' μl'
        required_taq_vol_str = str(index_taq_vol/1000)+ ' μl'
        avail_taq_vol_str = str(taq_avail_vol)+ ' μl'

        index_cols = pcr_comps_area.columns(col_size)

        if index_capacity > (index_avail - index_remain):
            index_pairs_allowed = '<p style="color:green">'+str(index_capacity)+'</p>'
        else:
            index_pairs_allowed = '<p style="color:red">'+str(index_capacity)+'</p>'

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

    if comp_warning:
        msg = 'Provide: ' + comp_warning + ' to fulfill PCR requirements'
    else:
        msg = None
        
    return msg

    # comp_warning_area.markdown(f'<p style="color:#f63366">{comp_warning}</p>', unsafe_allow_html=True)


#def display_primer_components(assay_usage, expander=True):
#    #Values
#    if assay_usage:
#        primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
#                st.session_state['experiment'].get_volumes_required(assay_usage=assay_usage)
#    else:
#        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
#            0,{},0,0,0,0 

#    primer_avail_counts, primer_avail_vols = st.session_state['experiment'].get_primers_avail()
#    per_use_vol = util.CAP_VOLS['384PP_AQ_BP'] - util.DEAD_VOLS['384PP_AQ_BP']
#    primer_names = set(primer_vols.keys())

#    primer_array = []
#    for p in primer_names:
#        need_vol = primer_vols.get(p,0)
#        num_wells = ceil(need_vol/per_use_vol)

#        primer_array.append([p, num_wells, assay_usage.get(p,0), primer_vols.get(p,0)/100,\
#                    primer_avail_vols.get(p,0)/1000])
    
#    primer_df = pd.DataFrame(primer_array, columns=['Primer', 'Num Wells','Uses', 'Volume(μL)', 'Available Volume(μL)'])
    
#    primer_table = aggrid_interactive_table(primer_df, grid_height=350)
    
#    if expander:
#        primer_components_exp = st.expander('Primer Wells, Uses, and Volumes', expanded=False)
#        ptab1, ptab2, ptab3 = primer_components_exp.tabs(["Wells", "Uses", "Volume"])
#    else:
#        ptab1, ptab2, ptab3 = st.tabs(["Wells", "Uses", "Volume"])

#    #Set up columns
#    col_size = 3
#    num_cols = 6
#    columns = []
#    i = 0
#    while i < num_cols:
#        columns.append(col_size)
#        i+=1

    
#    wells_col = ptab1.columns(columns)
#    uses_col = ptab2.columns(columns)
#    volume_col = ptab3.columns([2, 2, 2, 2, 2, 2, 2, 2,2])

   

    

#    for i in range(6):
#        if (i+2) % 2 == 0:
#            wells_col[i].markdown('**Primer**')
#            uses_col[i].markdown('**Primer**')
#        else:
#            wells_col[i].markdown('**Wells**')
#            uses_col[i].markdown('**Uses**')
            
#    for i in range(9):
#        if i in [0, 3, 6]:
#            volume_col[i].markdown('**Primer**')
#        elif i in [1, 4, 7]:
#            volume_col[i].markdown('**Volume μL**')
#        elif i in [2, 5, 8]:
#            volume_col[i].markdown('**Available μL**')

#    for k,p in enumerate(sorted(primer_names)):
#        if p != '':
#            r = (k*2)%6  # horizontal offset
#            need_vol = primer_vols.get(p,0)
#            num_wells = ceil((need_vol)/per_use_vol)

#            wells_col[r+0].write(p)
#            wells_col[r+1].write(num_wells)

#            uses_col[r+0].write(p)
#            uses_col[r+1].write(assay_usage.get(p,0))

#            k=(k*3)%9
#            volume_col[k+0].write(p)
#            if primer_vols.get(p,0) > primer_avail_vols.get(p,0):
#                colour = 'red'
#            else:
#                colour = 'green'
#            volume_col[k+1].markdown(f'<p style="colour:{colour}">{primer_vols.get(p,0)/1000}</p>')
#            volume_col[k+2].write(primer_avail_vols.get(p,0)/1000)

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
            dpp_fn = str(Path(dpp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{dpp_fn}</p>',\
                    unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(dpp, 'rt'), file_name=dpp_fn, mime='text/csv', key='dna_download_'+dpp)

    if not primer_picklist_paths:
        error_msgs.append('No primer picklist available')
    else:
        for ppp in primer_picklist_paths:
            ppp_fn = str(Path(ppp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ppp_fn}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(ppp, 'rt'), file_name=ppp_fn, mime='text/csv', key='primer_download_'+dpp)
            
    if not taqwater_picklist_paths:
        error_msgs.append('No taq/water picklist available')
    else:
        for tpp in taqwater_picklist_paths:
            tpp_fn = str(Path(tpp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_fn}</p>',\
                        unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", 
                    data=open(tpp, 'rt'), file_name=tpp_fn, mime='text/csv', key='taqwater_download_'+tpp)

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
            ipp_fn = str(Path(ipp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{ipp_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                    data=open(ipp, 'rt'), file_name=ipp_fn, mime='text/csv', key='index_download_'+ipp)

    if not amplicon_picklist_paths:
        error_msgs.append('No amplicon picklist available')               
    else:
        for app in amplicon_picklist_paths:
            app_fn = str(Path(app).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{app_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download", \
                    data=open(app, 'rt'), file_name=app_fn, mime='text/csv', key='amplicon_download_'+app)

    if not taqwater_picklist_paths:
        error_msgs.append('No taq/water picklist available')

    else:
        for tpp in taqwater_picklist_paths:
            tpp_fn = str(Path(tpp).name)
            picklist_file_col.markdown(\
                        f'<p style="text-align:right;color:#4b778c;padding:5px">{tpp_fn}</p>',\
                                unsafe_allow_html=True)
            picklist_btn_col.download_button(label=f"Download",\
                        data=open(tpp, 'rt'), file_name=tpp_fn, mime='text/csv',\
                                key='taqwater_download_'+tpp)
    if error_msgs:
        for msg in error_msgs:
            picklist_file_col.markdown(f'<p style="color:#ff0000;text-align:right">{msg}</p>',\
                    unsafe_allow_html=True)
    
def display_status(exp, height=350):
    """
    Display the progress in the pipeline for this experiment
    Should use aggrid to display the stages and the changes at each stage
    """
    steps, header = exp.get_stages()
    status_df = pd.DataFrame(steps, columns=header)
    #status_df.reset_index(inplace=True)
    #status_df = status_df.rename(columns = {'index':'Steps', 'pending':'Pending Steps'})
    status_df = aggrid_interactive_table(status_df, grid_height=height, key='status_aggrid')

def display_primers(exp, assay_usage, height=350):
    """
    Alternative display_primer_components using aggrid
    Designed as a component for a generic display widget
    """
    if assay_usage:
        primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
                exp.get_volumes_required(assay_usage=assay_usage)
    else:
        reactions, primer_vols, primer_taq_vol, primer_water_vol, index_taq_vol, index_water_vol =\
            0,{},0,0,0,0
    
    primer_avail_counts, primer_avail_vols = st.session_state['experiment'].get_primers_avail()
    per_use_vol = util.CAP_VOLS['384PP_AQ_BP'] - util.DEAD_VOLS['384PP_AQ_BP']
    primer_names = set(primer_vols.keys())
    warning_primers = ''
    primer_array = []
    for p in primer_names:
        need_vol = primer_vols.get(p,0)
        num_wells = ceil(need_vol/per_use_vol)

        if primer_vols.get(p, 0)/100 > primer_avail_vols.get(p,0)/1000:
            warning_primers += p + ', '

        primer_array.append([p, num_wells, assay_usage.get(p,0), primer_vols.get(p,0)/100,\
                    primer_avail_vols.get(p,0)/1000])
    
    # if warning_primers != '':
    #     st.warning('The following primers do not have enough volume: ' + warning_primers[:-1])
    primer_df = pd.DataFrame(primer_array, columns=['Primer', 'Num Wells','Uses', 'Required Volume(μL)', 'Available Volume(μL)'])
    primer_table = aggrid_interactive_table(primer_df, grid_height=height, key='primer_display')
    
def display_plates(exp, plate_usage, height=350):
    """
    Simple plate summary using aggrid.
    To be used as a generic viewing component
    """
    #st.write(exp.get_musterer_pids())
    plate_df = pd.DataFrame(plate_usage, columns=['Plate', 'Num Wells', 'Purpose'])
    plate_table = aggrid_interactive_table(plate_df, grid_height=height, key='plate_aggrid')

def view_plates(exp):
    """
    Visual view of the plates in the experument
    """
    plate_ids = []
    for pid in exp.plate_location_sample:
        plate_ids.append(f"{exp.plate_location_sample[pid]['purpose']} plate: {util.unguard(pid)}")

    #Let user choose which plate to view
    _, col1, _ = st.columns([2, 1,2])
    plate_selectbox = col1.selectbox('Plate ID to view', plate_ids, key='plate_viewer')
    if plate_selectbox:
        plate_id = util.guard_pbc(plate_selectbox.split(':')[1])
        if plate_id in exp.plate_location_sample:
            heatmap_str = generate_heatmap_html(exp, plate_id, scaling=1.3)

            #with open("debug.html", 'wt') as outf:
            #    print(heatmap_str, file=outf)
            components.html(heatmap_str, height=700, scrolling=True)
        else:
            plate_barcode_error_msg = "Plate barcode not found in experiment"
            st.markdown(f'<p style="color:#FF0000">{plate_barcode_error_msg}</p>',\
            unsafe_allow_html=True)


def display_files(exp, file_usage, height=350):
    """
    Display info for all files that have so far been uploaded into the experiment
    Give info on name, any plates they contain, whether they are required so far, etc
    """
    file_df = pd.DataFrame.from_dict(file_usage, orient='index')
    file_df.reset_index(inplace=True)
    file_df = file_df.rename(columns = {'index':'File', 'plates':'Plates', 'purpose':'Purpose'})
    file_table = aggrid_interactive_table(file_df, grid_height=height, key='file_aggrid')

def display_indexes(exp, assay_usage, height=350):
    """
    Display info for all indexes that have been uploaded into the experiment
    """
    fwd_idx, rev_idx, warning_idxs = exp.get_index_avail()
    indexes = {**fwd_idx, **rev_idx}
    for index in indexes:
        indexes[index]['avail_vol'] = sum(indexes[index]['avail_vol'])/1000
    # if fwd_idx_vols:
    #     idx_remaining, max_idx_pairs, reaction_vol_capacity, fwd_idx_reactions =\
    #                 exp.get_index_remaining_available_volume(assay_usage, fwd_idx_vols, rev_idx_vols)
    #     print(idx_remaining, max_idx_pairs, reaction_vol_capacity)
    #     per_use_vol = util.CAP_VOLS['384PP_AQ_BP'] - util.DEAD_VOLS['384PP_AQ_BP']
    #     fwd_idx_names = fwd_idx_vols.keys()
        
    #     for p in fwd_idx_names:
    #         need_vol = idx_vols.get(p,0)
    #         num_wells = ceil(need_vol/per_use_vol)
    #         index_array.append([p, num_wells, idx_vols.get(p,0)/100,\
    #                     fwd_idx_vols.get(p,0)[0]/1000])
    #     rev_idx_names = rev_idx_vols.keys()
    #     for p in rev_idx_names:
    #         need_vol = idx_vols.get(p,0)
    #         num_wells = ceil(need_vol/per_use_vol)
    #         index_array.append([p, num_wells, idx_vols.get(p,0)/100,\
    #                     rev_idx_vols.get(p,0)[0]/1000])
    if warning_idxs != '':
        st.warning("The following indexes do not have enough volume: " + warning_idxs[:-2])

    index_df = pd.DataFrame.from_dict(indexes,  orient='index')
    index_df.reset_index(inplace=True)
    index_df = index_df.rename(columns = {'index':'Index', 'count':'Wells', 'req_vol':'Required Volume (μL)', 'avail_vol':'Available Volume (μL)'})
    index_table = aggrid_interactive_table(index_df,grid_height=height,key='index_display')

def display_log(exp, height=250):
    """
    Display the log
    """
    func = sys._getframe(1).f_code.co_name
    func_line = inspect.getframeinfo(sys._getframe(1)).lineno
    print("Function", func, type(func_line))
    # display the experiment log and let the user filter the view in a number of ways
    log_entries = st.session_state['experiment'].get_log(100)
    if len(log_entries) == 0:
        st.write('No entries currently in the log')
    else:
        df = pd.DataFrame(log_entries, columns=exp.get_log_header())
        aggrid_interactive_table(df, grid_height=height, key='logs')

def info_viewer(key):
    """
    Container for displaying module info functions, each of which provides a dataframe for display in an aggrid.
    Because aggrid allows selection, each module can also handle a standard set of operations (such as delete).
    """
    exp = st.session_state['experiment']
    info_expander = st.expander('Info Panel')
    container = info_expander.container()
        
    col1,col2,col3 = info_expander.columns([1 ,7,1])
    with col1:
        st.write('')
        st.markdown('<h5 style="color:#7eafc4">Display panel</h5>', unsafe_allow_html=True)
    with col3:
        view_height = st.number_input('Set display height', min_value=50, max_value=700, value=350, step=25, help="Size of display grid", key=key)
    with col2:
        view_tab = stx.tab_bar(data=[
            stx.TabBarItemData(id=1, title="Status", description=""),
            stx.TabBarItemData(id=2, title="Files", description=""),
            stx.TabBarItemData(id=3, title="Plates", description=""),
            stx.TabBarItemData(id=4, title="Plate Viewer", description=""),
            stx.TabBarItemData(id=5, title="Primers", description=""),
            stx.TabBarItemData(id=6, title="Indexes", description=""),
            stx.TabBarItemData(id=7, title="Log", description="")
        ], return_type=int, default=1)

    
    if view_tab == 1:
        # Status tab should tell us where we are up to in the pipeline and what's happened so far
        with container:
            #with st.expander('Pipeline status', expanded=False):
                display_status(exp)
                
    if view_tab == 2:
        #view_expander = container.expander(label='All uploaded files', expanded=False)
        file_usage = exp.get_file_usage()
        with container:
                display_files(exp, file_usage, height=view_height)


    if view_tab == 3:
        plate_usage = exp.get_plate_usage()
        with container:
                display_plates(exp, plate_usage, height=view_height)


    if view_tab == 4:
        with container:
            view_plates(exp)


    if view_tab == 5:
        assay_usage = exp.get_assay_usage()
        with container:
            display_primers(exp, assay_usage, height=view_height)


    if view_tab == 6:
        assay_usage = exp.get_assay_usage()
        with container:
            display_indexes(exp, assay_usage, height=view_height)


    if view_tab == 7:
        with container:
            display_log(exp, height=view_height)

