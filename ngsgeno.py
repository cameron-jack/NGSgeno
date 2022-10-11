# from asyncio.windows_utils import pipe
# from distutils.command.upload import upload
import os
from ssl import SSLSession
import sys
from pathlib import PurePath
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO
from xml.etree.ElementInclude import include
import subprocess

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
#from bin.ngsmatch import match_alleles

import extra_streamlit_components as stx
from display_components import data_table, display_pcr_components, display_primer_components
from load_data import *
import asyncio

def get_run_folders():
    """ return an alphabetically sorted list of run folders, without their run_ prefix """
        
    # get run folder names without 'run_'
    run_folders = [''] + [d[4:] for d in os.listdir('.') if d.startswith('run_') and os.path.isdir(d)]
    return sorted(run_folders)

def create_run_folder(newpath):
    """ returns Experiment or None, and a message string. Requires a full path and generates a new experiment file """
    print('Attempting to create: ', newpath)
    if not os.path.exists(newpath):
        try:
            os.mkdir(newpath)
        except Exception as exc:
            output_error(exc, msg='Could not create new folder: ' + newpath)
            return None, 'Failed to create new folder: ' + newpath
    else:
        if os.path.exists(os.path.join(newpath, EXP_FN)):
            return None, 'Experiment already exists with this name: ' + newpath
    print('Generating experiment: ', newpath)
    exp = Experiment(name=newpath.lstrip('run_'))
    return exp, ''

#generalise
def generate_picklist(plates, pcr_stage=1):
    picklist_file_col, picklist_btn_col = st.columns(2)
    if pcr_stage == 1:
        DNA_plates = plates[0]
        PCR_plates = plates[1]
        taqwater_plates = plates[2]
        error_msgs = []
        
        success = st.session_state['experiment'].generate_echo_PCR1_picklist_interface(DNA_plates,\
                PCR_plates, taqwater_plates)
                        
        if success:
            dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths =\
                st.session_state['experiment'].get_echo_PCR1_picklist_filepaths()

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
    
    if pcr_stage == 2:                      
        PCR_plates = plates[0]
        taqwater_plates = plates[1]
        index_plates = plates[2]
        amplicon_plates = plates[3]
        error_msgs = []

        success = st.session_state['experiment'].generate_echo_PCR2_picklist_interface(\
                                PCR_plates, index_plates, taqwater_plates,\
                                            amplicon_plates)
        if success:
            index_picklist_paths, taqwater_picklist_paths, amplicon_picklist_paths =\
                    st.session_state['experiment'].get_echo_PCR2_picklist_filepaths()
            
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

def plate_checklist_expander(available_nimbus, pcr_stage=1):
    included_PCR_plates = set()
    included_taqwater_plates = set()
    checklist_col = st.columns(4)
    pcr_plate_title = "PCR Plates"
    taqwater_plate_title = "Taq/Water Plates"

    if pcr_stage == 1:
        included_DNA_plates = set()
        dna_plate_title = "DNA Plates"
        
        checklist_col[0].markdown(f'**{dna_plate_title}**')
        for nim in available_nimbus:
            echo_filename=Path(nim).stem
            inc_dna = checklist_col[0].checkbox(echo_filename, value=True, key='chk_box_dna_'+nim)
            if inc_dna:
                included_DNA_plates.add(echo_filename.split('_')[-2])
    
        checklist_col[1].markdown(f'**{pcr_plate_title}**')
        for pcr_pid in st.session_state['experiment'].get_pcr_plates():
            inc_pcr = checklist_col[1].checkbox(file_io.unguard_pbc(pcr_pid, silent=True),\
                            value=True, key='chk_box_pcr_'+pcr_pid)
            if inc_pcr:
                included_PCR_plates.add(pcr_pid)

        checklist_col[2].markdown(f'**{taqwater_plate_title}**')
        for taqwater_pid in st.session_state['experiment'].get_taqwater_plates():
            inc_taqwater = checklist_col[2].checkbox(file_io.unguard_pbc(taqwater_pid, silent=True), 
                        value=True, key='chk_box_taqwater_'+taqwater_pid)
            if inc_taqwater:
                included_taqwater_plates.add(taqwater_pid)
        
        return included_DNA_plates, included_PCR_plates, included_taqwater_plates

    if pcr_stage == 2:
        included_index_plates = set()
        included_amplicon_plates = set()
        #could make a for loop

        index_plate_title = "Index Plates"
        amplicon_plate_title = "Amplicon Plates"

        checklist_col[0].markdown(f'**{pcr_plate_title}**')
        for pcr_pid in st.session_state['experiment'].get_pcr_plates():
            inc_pcr = checklist_col[0].checkbox(file_io.unguard_pbc(pcr_pid, silent=True),
                    value=True, key='chk_box_pcr_'+pcr_pid)
            if inc_pcr:
                included_PCR_plates.add(pcr_pid)

        checklist_col[1].markdown(f'**{taqwater_plate_title}**')
        for taqwater_pid in st.session_state['experiment'].get_taqwater_plates():
            inc_taqwater = checklist_col[1].checkbox(file_io.unguard_pbc(taqwater_pid, silent=True), 
                    value=True, key='chk_box_taqwater_'+taqwater_pid)
            if inc_taqwater:
                included_taqwater_plates.add(taqwater_pid)

        checklist_col[2].markdown(f'**{index_plate_title}**')
        for index_pid in st.session_state['experiment'].get_index_plates():
            inc_index = checklist_col[2].checkbox(file_io.unguard_pbc(index_pid, silent=True),
                                                value=True, key='chk_box_index_'+index_pid)
            if inc_index:
                included_index_plates.add(index_pid)
        
        checklist_col[3].markdown(f'**{amplicon_plate_title}**')
        # for amplicon_pid in st.session_state['experiment'].get_amplicon_plates():
        #         amplicon_index = checklist_col[3].checkbox(file_io.unguard_pbc(amplicon_pid, silent=True), 
        #                                                      value=True, key='chk_box_amplicon_'+amplicon_pid)
        #         if amplicon_index:
        #             included_amplicon_plates.add(amplicon_pid)

        return included_PCR_plates, included_taqwater_plates, included_index_plates, included_amplicon_plates
 
    
async def report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog):
    while True:
        progress_files = list(Path(rundir).glob('match_progress_*'))
        if len(progress_files) == 0:
            launch_progress = 0
            match_progress = 0
        elif len(progress_files) == 1:
            launch_progress = int(str(progress_files[0]).split('_')[-2])
            match_progress = int(str(progress_files[0]).split('_')[-1])
        
        launch_msg.write('Allele calling task launch progress: '+str(launch_progress)+'%')
        launch_prog.progress(launch_progress)
        completion_msg.write('Allele calling task completion progress: '+str(match_progress)+'%')
        match_prog.progress(match_progress)

        if launch_progress == 100 and match_progress == 100:
            st.write('Analysis completed')
            st.session_state['match_running'] = False
            return
        await asyncio.sleep(1)

def main():
    st.set_page_config(
        page_title="NGS Genotyping",
        page_icon="ngsg_icon.png",
        layout="wide"
    )

    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')

    if 'folder' not in st.session_state:
        st.session_state['folder'] = None
    #Folder inputs
    title = "NGS Genotyping Pipeline"
    _,title_column,new_folder_col, ex_folder_col = st.columns([1,4,4,4])
    title_column.markdown(f'<h2>{title}</h2>', unsafe_allow_html=True)

    add_run_folder = new_folder_col.text_input('Create new run folder:')
    #create_run_folder_button = ftab1.button(label='Create', key='create_run_folder_button')

    try:
        existing_run_folders = get_run_folders()
    except Exception as exc:
        output_error(exc, msg='Cannot locate NGSgeno folder')
        return

    run_folder = ex_folder_col.selectbox("Select a run folder to open", existing_run_folders)
    error_msg=''

    #error message comes up now after entering a new folder name - need to fix
    if add_run_folder:
        add_run_folder_str = 'run_' + add_run_folder
        exp, msg = create_run_folder(add_run_folder_str)
        if exp:
            exp.save()                                                                      
            st.session_state['experiment'] = exp
            st.experimental_rerun()
        else:
            if 'already exists' in msg:
                error_msg = "Folder name already exists"
            else:
                error_msg = "Fatal path error: " + msg
        
        new_folder_col.markdown(f'<p style="color:#FF0000; text-align:center">{error_msg}</p>',\
                unsafe_allow_html=True)

    if 'stage' not in st.session_state:
            st.session_state['stage'] = None

    if run_folder:
        if st.session_state['experiment'] == None or st.session_state['experiment'].name != run_folder:
            
            ch_run_path = 'run_' + run_folder
            if os.path.exists(ch_run_path):
                exp = load_experiment(ch_run_path)
                st.session_state['folder'] = 'existing'
                #print(dir(exp))
                if not exp:
                    error_msg = "Could not load experiment from: "+ ch_run_path
                elif ch_run_path.endswith(exp.name):
                    # success!
                    st.session_state['experiment'] = exp
                    st.experimental_rerun()
                else:
                    error_msg = "Invalid experiment file in: " + ch_run_path
        
        new_folder_col.markdown(f'<p style="color:#FF0000; text-align:center">{error_msg}</p>',\
                unsafe_allow_html=True)

    if st.session_state['experiment']:
        exp = st.session_state['experiment']
        
        pipeline_stages=["Load", "Nimbus", "Primers", "Index", "Miseq", "Alleles", "Genotyping", "Review"]
        pipeline_stage = stx.stepper_bar(steps=pipeline_stages, lock_sequence=False)

        stage_list = ['load', 'nimbus', 'primer', 'index','miseq', 'allele', 'genotyping', 'review']

        if st.session_state['folder']:
            for i in range(6):
                if st.session_state['stage'] == stage_list[i]:
                    pipeline_stage = i
            
        #Load data
        if pipeline_stage == 0:
            load_data_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="Load Data", description=""),
                stx.TabBarItemData(id=2, title="View Data", description="")
            ], return_type=int)

            if not load_data_tab:
                if 'load_tab' not in st.session_state:
                    st.session_state['load_tab'] = 1
                load_data_tab = st.session_state['load_tab']

            #view data
            if load_data_tab == 1:
                load_rodentity_data()
                load_database_data()
                load_custom_csv()
                st.session_state['load_tab'] = 1

            #load data
            if load_data_tab == 2:
                assay_usage = st.session_state['experiment'].get_assay_usage()
                data_table(key=load_data_tab, options=True)
                display_primer_components(assay_usage, expander=True)
                st.session_state['load_tab'] = 2
        
        #Nimbus
        if pipeline_stage == 1:
            
            nim_tab1_err = ''
            nimbus_title = ''

            nimbus_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="Download", description="Nimbus input files"),
                stx.TabBarItemData(id=2, title="Upload", description="Nimbus output files"),
                stx.TabBarItemData(id=3, title="View Data", description="")
            ], return_type=int)

            if not nimbus_tab:
                if 'nimbus_tab' not in st.session_state:
                    st.session_state['nimbus_tab'] = 1
                nimbus_tab = st.session_state['nimbus_tab']

            if not st.session_state['experiment'].dest_sample_plates:
                nimbus_title = "Load data inputs to enable Nimbus input file generation."
            else: 
                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                # do we have any Nimbus inputs to generate + download
                yet_to_run = len(st.session_state['experiment'].dest_sample_plates) - len(nfs)
            
                if yet_to_run and yet_to_run > 0: 
                    plates_to_run = [dest_plate for dest_plate in exp.dest_sample_plates\
                                 if all([dest_plate not in nf for nf in nfs])]
                    plates_to_run_str = '\n'.join(plates_to_run)
                    success = st.session_state['experiment'].generate_nimbus_inputs()
                    if not success:
                        nim_tab1_err = "Failed to generate Nimbus files. Please read the log."
                    else:
                        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()

                    nim_tab2_title = str(yet_to_run) + " 96-well plate sets need Nimbus input file generation"

                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            
            #download nimbus
            if nimbus_tab == 1:
                st.markdown(f'<h5 style="text-align:center;color:#f63366">{nimbus_title}</h5>',\
                         unsafe_allow_html=True)
                st.write('')

                if nim_tab1_err:
                    st.write(f'<p style="color:#FF0000">{nim_tab1_err}</p>', unsafe_allow_html=True)

                _,dl_col1,dl_col2,dl_col3,dl_col4,_= st.columns([1,9,6,9,6,1])
                
                #Generate file names to download + download buttons
                print(f"{nfs=} {efs=} {xbcs=}")
                for i,nf in enumerate(nfs):
                    nimbus_fn=Path(nf).name

                    if (i+1) % 2 != 0:
                        dl_col1.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{nimbus_fn}</p>',\
                                 unsafe_allow_html=True)

                        dl_col2.download_button("Download ", open(nf), file_name=nimbus_fn, 
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                
                    else:
                        dl_col3.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{nimbus_fn}</p>',\
                                 unsafe_allow_html=True)
                 
                        dl_col4.download_button("Download ", open(nf), file_name=nimbus_fn,\
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}") 
                st.session_state['nimbus_tab'] = 1

            #upload nimbus
            if nimbus_tab == 2:
                nim_tab2_title_area = st.empty()
                nim_tab2_title = ''
                nim_upload_area = st.container()
                files_str = None

                if not efs and not xbcs:
                    nim_tab2_title = "Load data inputs to enable Nimbus input file generation."
                elif efs and not xbcs:
                    nim_tab2_title = 'All Nimbus outputs now uploaded'
                    uploaded_nims = [Path(ef).name for ef in efs]
                    files_str = '</br>'.join(uploaded_nims)
                    st.markdown(f'<p style="text-align:center;color:#17754d">{files_str}</p>', unsafe_allow_html=True)
                elif xbcs:
                    nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                    missing_nims = ['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs]
                    files_str = '</br>'.join(missing_nims)

                    nim_outputs = st.file_uploader('Echo_384_COC_0001_....csv', type='csv', 
                            accept_multiple_files=True, help='You can upload more than one file at once')

                    if nim_outputs: # and nim_outputs != st.session_state['nim_upload']:
                        #st.session_state['nim_uploads'] = nim_outputs
                        for nim_output in nim_outputs:
                            #print(f"{nim_output.name=}")
                            fp = st.session_state['experiment'].get_exp_fp(nim_output.name)
                            st.session_state['experiment'].log(f"Info: copying {fp} to experiment folder")
                            with open(fp, 'wt') as outf:
                                #print(nim_output.getvalue().decode("utf-8"))
                                outf.write(nim_output.getvalue().decode("utf-8").replace('\r\n','\n'))
        
                nim_tab2_title_area.markdown(f'<h5 style="text-align:center;color:#f63366">{nim_tab2_title}</h5>',\
                                unsafe_allow_html=True)
                st.session_state['nimbus_tab'] = 2

            #view data
            if nimbus_tab == 3:
                data_table(key=load_data_tab, options=True)
                st.session_state['nimbus_tab'] = 3

        #Primer PCR
        if pipeline_stage == 2:
            primer_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="PCR 1", description="Components"),
                stx.TabBarItemData(id=2, title="Provide", description="Plates"),
                stx.TabBarItemData(id=3, title="Generate", description="Picklists")
            ], return_type=int)

            if not primer_tab:
                if 'primer_tab' not in st.session_state:
                    st.session_state['primer_tab'] = 1
                primer_tab = st.session_state['primer_tab']
            
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            missing_nims = ['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs]

            if primer_tab == 1:
                primer_checklist_exp = st.expander('Plate Checklist', expanded=False)
                if efs:
                    with primer_checklist_exp:
                        included_DNA_plates, included_PCR_plates, included_taqwater_plates =\
                                    plate_checklist_expander(efs, pcr_stage=1)
                    
                    if included_DNA_plates:
                        assay_usage = st.session_state['experiment'].get_assay_usage(dna_plate_list=included_DNA_plates)
                        #assay_usage = st.session_state['experiment'].get_assay_usage()
                        display_pcr_components(assay_usage, 1)
                else:
                    no_nimbus_msg = "Load Nimbus output files to enable PCR stages"
                    st.markdown(f'<h5 style="text-align:center;color:#f63366">{no_nimbus_msg}</h5',\
                            unsafe_allow_html=True)
                st.session_state['primer_tab'] = 1
                
            #provide plates
            if primer_tab == 2:
                primer_checklist_exp = st.expander('Plate Checklist', expanded=False)
                if efs:
                    upload_pcr_files(1)
                    with primer_checklist_exp:
                        included_DNA_plates, included_PCR_plates, included_taqwater_plates =\
                                plate_checklist_expander(efs, pcr_stage=1)
                st.session_state['primer_tab'] = 2

            #generate picklists
            if primer_tab == 3:
                primer_checklist_exp = st.expander('Plate Checklist', expanded=False)
                if efs:
                    with primer_checklist_exp:
                        included_DNA_plates, included_PCR_plates, included_taqwater_plates =\
                                plate_checklist_expander(efs,pcr_stage=1)

                    if included_DNA_plates:
                        if st.session_state['experiment'].check_ready_echo1(included_DNA_plates,\
                                    included_PCR_plates, included_taqwater_plates):

                            plates = [included_DNA_plates, included_PCR_plates, included_taqwater_plates]
                            _,picklist_button_col,_ = st.columns([2, 2, 1])

                            echo_picklist_go = picklist_button_col.button('Generate Echo Picklist',\
                                    key='echo_pcr1_go_button')

                            picklist_button_col.write('')

                            if echo_picklist_go:
                                generate_picklist(plates, pcr_stage=1)
                    
                    else:
                        no_dna_msg = "Include at least one DNA plate to carry on with the pipeline"
                        st.markdown(f'<h5 style="text-align:center;color:#f63366">{no_dna_msg}</h5',\
                            unsafe_allow_html=True)
                st.session_state['primer_tab'] = 3

        #Index PCR
        if pipeline_stage == 3:
            index_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="PCR 2", description="Components"),
                stx.TabBarItemData(id=2, title="Provide", description="Plates"),
                stx.TabBarItemData(id=3, title="Generate", description="Picklists")
            ], return_type=int)

            if not index_tab:
                if 'index_tab' not in st.session_state:
                    st.session_state['index_tab'] = 1
                index_tab = st.session_state['index_tab']

            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
            available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]
            #PCR components
            if index_tab == 1:
                index_checklist_exp = st.expander('Plate Checklist', expanded=False)
                if available_nimbus:
                    with index_checklist_exp:
                        included_PCR_plates, included_taqwater_plates, included_index_plates,\
                                    included_amplicon_plates =\
                                    plate_checklist_expander(available_nimbus, pcr_stage=2)
                
                    assay_usage = st.session_state['experiment'].get_assay_usage()
                    display_pcr_components(assay_usage, 2)
                else:
                    no_nimbus_msg = "Load Nimbus output files to enable PCR stages"
                    st.markdown(f'<h5 style="text-align:center;color:#f63366">{no_nimbus_msg}</h5',\
                            unsafe_allow_html=True)
                st.session_state['index_tab'] = 1

            #plate info
            if index_tab == 2:
                index_checklist_exp = st.expander('Plate Checklist', expanded=False)
                if available_nimbus:
                    with index_checklist_exp:
                        included_PCR_plates, included_taqwater_plates, included_index_plates,\
                                    included_amplicon_plates =\
                                    plate_checklist_expander(available_nimbus, pcr_stage=2)

                    upload_pcr_files(2)
                st.session_state['index_tab'] = 2

            #generate picklist
            picklist_err = ''
            if index_tab == 3:
                index_checklist_exp = st.expander('Plate Checklist', expanded=False)
                if available_nimbus:
                    with index_checklist_exp:
                        included_PCR_plates, included_taqwater_plates, included_index_plates,\
                                    included_amplicon_plates =\
                                    plate_checklist_expander(available_nimbus, pcr_stage=2)

                    if st.session_state['experiment'].check_ready_echo2(included_PCR_plates,\
                                included_taqwater_plates, included_index_plates):

                        _,picklist_button_col,_ = st.columns([2, 2, 1])

                        plates = [included_PCR_plates, included_taqwater_plates, included_index_plates,\
                                    included_amplicon_plates]

                        echo_picklist_go = picklist_button_col.button('Generate Echo Picklists',\
                                    key='echo_pcr2_go_button')

                        picklist_button_col.write('')

                        if echo_picklist_go:
                            generate_picklist(plates, pcr_stage=2)

                    else:
                        picklist_err = "Include at lease one PCR plate, "+\
                            "one index plate and one taq/water plate to carry on with the pipeline"
                
                        st.markdown(f'<h5 style="color:#f63366;text-align:center">{picklist_err}</h5>',\
                                unsafe_allow_html=True)
                st.session_state['index_tab'] = 3
            
        #Miseq
        if pipeline_stage == 4:
            miseq_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="Download", description="Samplesheet"),
                stx.TabBarItemData(id=2, title="Upload", description="Sequence Files"),
            ], default=1, return_type=int)

            if not miseq_tab:
                if 'miseq_tab' not in st.session_state:
                    st.session_state['miseq_tab'] = 1
                miseq_tab = st.session_state['miseq_tab']

            st.write('')
            exp = st.session_state['experiment']
            if miseq_tab == 1:
                _, miseq_col1, miseq_col2, _ =  st.columns([2,1,1,2])

                if exp.get_miseq_samplesheets():
                    for fp in exp.get_miseq_samplesheets():
                        fp_name = str(Path(fp).name)
                        miseq_col1.markdown(f'<strong style="color:#486e7a">{fp_name}</strong>', unsafe_allow_html=True)
                    
                        download_miseq = miseq_col2.button(label='Download', key='dnld_samplesheet_'+str(fp))
                else:
                    no_miseq_msg = ""
                st.session_state['miseq_tab'] = 1

            if miseq_tab == 2:
                upload_miseq_fastqs(exp)
                st.session_state['miseq_tab'] = 2

        #Allele Calling
        if pipeline_stage == 5:
            exp = st.session_state['experiment']

            allele_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="Upload", description="Sequence Files"),
                stx.TabBarItemData(id=2, title="Allele Calling", description=""),
                stx.TabBarItemData(id=3, title="View Data", description=""),
            ], key='allele_tab_bar' , return_type=int)

            if not allele_tab:
                if 'allele_tab' not in st.session_state:
                    st.session_state['allele_tab'] = 1
                allele_tab = st.session_state['allele_tab']

            if allele_tab == 1:
                st.write('')
                upload_miseq_fastqs(exp)
                st.session_state['allele_tab'] = 1

            if allele_tab == 2:
                rundir = exp.get_exp_dir()

                num_unique_seq = st.number_input("Number of unique sequences per work unit", value=1)
                num_cpus = st.number_input(\
                            label="Number of processes to run simultaneously (defaults to # of CPUs)",\
                                    value=os.cpu_count())
                exhaustive_mode = st.checkbox(\
                            "Exhaustive mode: try to match every sequence, no matter how few counts")
                
                progress_files = list(Path(rundir).glob('match_progress_*'))
                if len(progress_files) == 1:
                    st.session_state['match_running'] = True
                else:
                    st.session_state['match_running'] = False

                disable_calling = False
                if Path(exp.get_exp_fp('ngsgeno_lock')).exists():
                    st.markdown('**Analysis in progress**')
                    disable_calling = True
                
                do_matching = st.button("Run allele calling", disabled=disable_calling)

                if do_matching:
                    success = exp.generate_targets()
                    if not success:
                        exp.log('Critical: failed to save reference sequences to target file')
                    else:
                        matching_prog = os.path.join('bin','ngsmatch.py')
                        cmd_str = f'python {matching_prog} --ncpus {num_cpus} --chunk_size {num_unique_seq} --rundir {rundir}'
                        if exhaustive_mode:
                            cmd_str += ' --exhaustive'                     
                        exp.log(f'Info: {cmd_str}')
                        cp = subprocess.run(cmd_str.split(' '))
                        st.session_state['match_running'] = True
                        exp.log(f'Debug: {cp}')
                    do_matching = None

                if 'match_running' in st.session_state and st.session_state['match_running']:
                    launch_msg = st.empty()
                    launch_prog = st.progress(0)
                    completion_msg = st.empty()
                    match_prog = st.progress(0)
                    asyncio.run(report_progress(rundir, launch_msg, launch_prog, completion_msg, match_prog))
                    
                st.session_state['allele_tab'] = 2

            if allele_tab == 3:
                data_table(key='allele_calling', options=True)
                st.session_state['allele_tab'] = 3
            
        
        #Genotyping
        if pipeline_stage == 6:
            pass

        #Review
        if pipeline_stage == 7:
            pass
        
        st.session_state['stage'] = stage_list[pipeline_stage]
        st.session_state['folder'] = None

if __name__ == '__main__':
    main()


