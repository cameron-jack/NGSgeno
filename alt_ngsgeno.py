from distutils.command.upload import upload
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

import extra_streamlit_components as stx
from display_components import data_table, display_pcr_components
from load_data import load_custom_csv, load_database_data, load_rodentity_data, upload_pcr_files

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

def generate_picklist(DNA_plates, PCR_plates, taqwater_plates):
    error_msgs = {}
    dna_picklist_path = primer_picklist_path = taqwater_picklist_path = None
    if DNA_plates:
        assay_usage = st.session_state['experiment'].get_assay_usage(dna_plate_list=DNA_plates)

        if st.session_state['experiment'].check_ready_echo1(DNA_plates, PCR_plates, taqwater_plates):

            success = st.session_state['experiment'].generate_echo_PCR1_picklist_interface(DNA_plates,\
                    PCR_plates, taqwater_plates)
                            
            if success:
                dna_picklist_path, primer_picklist_path, taqwater_picklist_path =\
                    st.session_state['experiment'].get_echo_PCR1_picklist_filepaths()

                if not dna_picklist_path:
                    error_msgs['dna'] = 'No DNA picklist available'
                if not primer_picklist_path:
                    error_msgs['primer'] = 'No primer picklist available'
                if not taqwater_picklist_path:
                    error_msgs['taqwater'] = 'No taq/water picklist available'
            else:
                error_msgs['no_dna'] = 'Include at least one DNA plate to generate picklist'
    
    return dna_picklist_path, primer_picklist_path, taqwater_picklist_path, error_msgs
            

def main():
    st.set_page_config(
        page_title="NGS Genotyping",
        page_icon="ngsg_icon.png",
        layout="wide"
    )

    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')
    #Folder inputs
    _,title_column,new_folder_col, ex_folder_col = st.columns([1,4,4,4])
    title_column.markdown('<h2>NGS Genotyping Pipeline</h2>', unsafe_allow_html=True)

    add_run_folder = new_folder_col.text_input('Create new run folder:')
    #create_run_folder_button = ftab1.button(label='Create', key='create_run_folder_button')

    try:
        existing_run_folders = get_run_folders()
    except Exception as exc:
        output_error(exc, msg='Cannot locate NGSgeno folder')
        return

    run_folder = ex_folder_col.selectbox("Select a run folder to open", existing_run_folders)

    error_msg_area = new_folder_col.empty()
    error_msg=None

    if add_run_folder:
        add_run_folder_str = 'run_' + add_run_folder
        exp, msg = create_run_folder(add_run_folder_str)
        if exp:
            exp.save()                                                                      
            st.session_state['experiment'] = exp
            add_run_folder = None
            st.experimental_rerun()
        else:
            if 'already exists' in msg:
                error_msg = "Folder name already exists"
            else:
                error_msg = "Fatal path error: " + msg

    if run_folder:
        if st.session_state['experiment'] == None or st.session_state['experiment'].name != run_folder:
            
            ch_run_path = 'run_' + run_folder
            if os.path.exists(ch_run_path):
                exp = load_experiment(ch_run_path)
                #print(dir(exp))
                if not exp:
                    error_msg = "Could not load experiment from: "+ ch_run_path
                elif ch_run_path.endswith(exp.name):
                    # success!
                    st.session_state['experiment'] = exp
                    #st.write(exp)
                    st.experimental_rerun()
                else:
                    error_msg = "Invalid experiment file in: " + ch_run_path
    
        ch_exp_folder = None

    if error_msg:
        error_msg_area.markdown(f'<p style="color:#FF0000; text-align:center">{error_msg}</p>',\
                 unsafe_allow_html=True)


    
    
    if st.session_state['experiment']:
        exp = st.session_state['experiment']

        pipeline_stages=["Load", "Nimbus", "Primers", "Index", "Miseq", "Alleles", "Genotyping", "Review"]
        pipeline_stage = stx.stepper_bar(steps=pipeline_stages, lock_sequence=False)

        #Load data
        if pipeline_stage == 0:

            load_data_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="Load Data", description=""),
                stx.TabBarItemData(id=2, title="View Data", description="")
            ], default=1)

            if load_data_tab == '1':
                load_rodentity_data()
                load_database_data()
                load_custom_csv()

            if load_data_tab == '2':
                data_table(key=1, options=True)
        
        #Nimbus
        if pipeline_stage == 1:
            
            nim_tab1_err = ''
            nim_tab2_title = ''
            nimbus_title = ''

            nimbus_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="Download", description=""),
                stx.TabBarItemData(id=2, title="Upload", description=""),
                stx.TabBarItemData(id=3, title="View Data", description="")
            ], default=1)

            if not st.session_state['experiment'].dest_sample_plates:
                nimbus_title = "Load data inputs to enable Nimbus input file generation."
            
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

            if nfs:
                nimbus_title = "Nimbus Input Files"

            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            
            if nimbus_tab == '1':

                st.markdown(f'<h5 style="text-align:center;color:#f63366">{nimbus_title}</h5>',\
                         unsafe_allow_html=True)
                st.write('')

                if nim_tab1_err:
                    st.write(f'<p style="color:#FF0000">{nim_tab1_err}</p>', unsafe_allow_html=True)

                _,dl_col1,dl_col2,dl_col3,dl_col4,_= st.columns([1,9,6,9,6,1])
                
                #Generate file names to download + download buttons
                for i,nf in enumerate(nfs):
                    file_name=nf.split('/')[1]

                    if (i+1) % 2 != 0:
                        dl_col1.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{file_name}</p>',\
                                 unsafe_allow_html=True)

                        dl_col2.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                
                    else:
                        dl_col3.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{file_name}</p>',\
                                 unsafe_allow_html=True)
                 
                        dl_col4.download_button("Download ", open(nf), file_name=PurePath(nf).name,\
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")    

            if nimbus_tab == '2':
                nim_tab2_title_area = st.empty()
                nim_upload_area = st.container()
                files_str = None

                if not efs and not xbcs:
                    nim_tab2_title = "Load data inputs to enable Nimbus input file generation."
                elif efs and not xbcs:
                    nim_tab2_title = 'All Nimbus outputs now uploaded'
                    uploaded_nims = [ef.split('/')[1] for ef in efs]
                    files_str = '</br>'.join(uploaded_nims)
                    # for file in efs:
                    #     file_name = file.split('/')[1]
                    #     st.markdown(f'<p style="text-align:center;color:#17754d">{file_name}</p>', unsafe_allow_html=True)
                elif xbcs:
                    nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                    nim_tab2_title = "Upload Nimbus output files"
                    missing_nims = ['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs]
                    files_str = '</br>'.join(missing_nims)

                    nim_outputs = st.file_uploader(' ', type='csv', 
                            accept_multiple_files=True, help='You can upload more than one file at once')

                    if nim_outputs: # and nim_outputs != st.session_state['nim_upload']:
                        #st.session_state['nim_uploads'] = nim_outputs
                        for nim_output in nim_outputs:
                            #print(f"{nim_output.name=}")
                            fp = st.session_state['experiment'].get_exp_fp(nim_output.name)
                            print(f"Copying {fp} to experiment folder")
                            with open(fp, 'wt') as outf:
                                #print(nim_output.getvalue().decode("utf-8"))
                                outf.write(nim_output.getvalue().decode("utf-8").replace('\r\n','\n'))
        
    
                nim_tab2_title_area.markdown(f'<h5 style="text-align:center;color:#f63366">{nim_tab2_title}</h5>',\
                                unsafe_allow_html=True)
                if files_str:
                    nim_upload_area.markdown(f'<p style="text-align:center;color:#4b778c">{files_str}</p>',\
                            unsafe_allow_html=True)


            if nimbus_tab == '3':
                data_table(key='nimbus', options=True)
        
        #Primer PCR
        if pipeline_stage == 2:
            primer_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="PCR Components", description=""),
                stx.TabBarItemData(id=2, title="Upload", description=""),
                stx.TabBarItemData(id=3, title="Generate Picklists", description="")
            ], default=1)

            dna_picklist_path = primer_picklist_path = taqwater_picklist_path = None
            error_msgs = {}

            if primer_tab == '1':
                primer_checklist_exp = st.expander('Plate Checklist', expanded=False)
                echo1_dna_col, echo1_pcr_col, echo1_taqwater_col = primer_checklist_exp.columns(3)
                echo1_dna_col.markdown('**DNA Plates**')
                echo1_pcr_col.markdown('**PCR Plates**')
                echo1_taqwater_col.markdown('**Taq/Water Plates**')

                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
                available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]
                
                if available_nimbus:
                    included_DNA_plates = set()
                    included_PCR_plates = set()
                    included_taqwater_plates = set()

                    for nim in available_nimbus:
                        echo_filename=nim.split('/')[1].split('.')[0]
                        inc_dna = echo1_dna_col.checkbox(echo_filename, value=True, key='chk_box_dna_'+nim)
                        if inc_dna:
                            included_DNA_plates.add(echo_filename.split('_')[-2])
                    
                    for pcr_pid in st.session_state['experiment'].get_pcr_plates():
                        inc_pcr = echo1_pcr_col.checkbox(file_io.unguard_pbc(pcr_pid, silent=True),\
                                     value=True, key='chk_box_pcr_'+pcr_pid)
                        if inc_pcr:
                            included_PCR_plates.add(pcr_pid)

                    for taqwater_pid in st.session_state['experiment'].get_taqwater_plates():
                        inc_taqwater = echo1_taqwater_col.checkbox(file_io.unguard_pbc(taqwater_pid, silent=True), 
                                    value=True, key='chk_box_taqwater_'+taqwater_pid)
                        if inc_taqwater:
                            included_taqwater_plates.add(taqwater_pid)
                    
                    if included_DNA_plates:
                        dna_picklist_path, primer_picklist_path, taqwater_picklist_path, error_msgs = \
                            generate_picklist(included_DNA_plates, included_PCR_plates, included_PCR_plates)
                    
                assay_usage = st.session_state['experiment'].get_assay_usage()
                display_pcr_components(assay_usage, 1)


            if primer_tab == '2':
                upload_pcr_files(1)

            if primer_tab == '3':
                _,picklist_button_col,_ = st.columns([2, 2, 1])
                echo1_dna_download_col, echo1_primer_download_col, echo1_taqwater_download_col = st.columns(3)

                echo_picklist_go = picklist_button_col.button('Generate Echo Picklist',\
                                         key='echo_pcr1_go_button')

                if echo_picklist_go:
                    #more tweaking to be done

                    if 'no_dna' in error_msgs:
                        st.markdown(f'{error_msgs.no_dna}', unsafe_allow_html=True)

                    if dna_picklist_path:
                        echo1_dna_download_col.download_button(label=f"Download {dna_picklist_path}",\
                                            data=open(dna_picklist_path, 'rt'), file_name=dna_picklist_path,\
                                                    mime='text/csv')
                    elif 'dna' in error_msgs:
                        echo1_dna_download_col.markdown(f'{error_msgs.dna}', unsafe_allow_html=True)

                    if primer_picklist_path:
                        echo1_primer_download_col.download_button(label=f"Download {primer_picklist_path}", 
                                            data=open(primer_picklist_path, 'rt'), file_name=primer_picklist_path,\
                                                    mime='text/csv')
                    elif 'primer' in error_msgs:
                        echo1_primer_download_col.markdown(f'{error_msgs.primer}', unsafe_allow_html=True)

                    if taqwater_picklist_path:
                        echo1_taqwater_download_col.download_button(label=f"Download {taqwater_picklist_path}", 
                                            data=open(taqwater_picklist_path, 'rt'), file_name=taqwater_picklist_path,\
                                                    mime='text/csv')
                    elif 'taqwater' in error_msgs:
                        echo1_taqwater_download_col.markdown(f'{error_msgs.taqwater}', unsafe_allow_html=True)


        #Index PCR
        if pipeline_stage == 3:

            index_tab = stx.tab_bar(data=[
                stx.TabBarItemData(id=1, title="PCR Components", description=""),
                stx.TabBarItemData(id=2, title="Upload", description=""),
                stx.TabBarItemData(id=3, title="Generate Picklists", description="")
            ], default=1)

            dna_picklist_path = primer_picklist_path = taqwater_picklist_path = error_msgs = None

            if index_tab == '1':
                index_checklist_exp = st.expander('Plate Checklist', expanded=False)
                echo2_pcr_col, echo2_index_col, echo2_taqwater_col, echo2_amplicon_col = index_checklist_exp.columns(4)
                echo2_pcr_col.markdown('**PCR Plates**')
                echo2_index_col.markdown('**Index Plates**')
                echo2_taqwater_col.markdown('**Taq/Water Plates**')
                echo2_amplicon_col.markdown('**Amplicon Plates**')

                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
                available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]

                if available_nimbus:
                    included_PCR_plates = set()
                    included_taqwater_plates = set()
                    included_index_plates = set()
                    included_amplicon_plates = set()

                    for pcr_pid in st.session_state['experiment'].get_pcr_plates():
                        inc_pcr = echo2_pcr_col.checkbox(file_io.unguard_pbc(pcr_pid, silent=True),
                                value=True, key='chk_box_pcr_'+pcr_pid)
                        if inc_pcr:
                            included_PCR_plates.add(pcr_pid)

                    for taqwater_pid in st.session_state['experiment'].get_taqwater_plates():
                        inc_taqwater = echo2_taqwater_col.checkbox(file_io.unguard_pbc(taqwater_pid, silent=True), 
                                value=True, key='chk_box_taqwater_'+taqwater_pid)
                        if inc_taqwater:
                            included_taqwater_plates.add(taqwater_pid)

                    for index_pid in st.session_state['experiment'].get_taqwater_plates():
                        inc_index = echo2_index_col.checkbox(file_io.unguard_pbc(index_pid, silent=True),
                                                            value=True, key='chk_box_index_'+index_pid)
                        if inc_index:
                            included_index_plates.add(index_pid)

            if index_tab == '2':
                upload_pcr_files(2)
            

        #miseq
        if pipeline_stage == 4:
            title = 'Download firle for Miseq'
            st.markdown(f'{title}')

        #Allele Calling
        if pipeline_stage == 5:
            pass

if __name__ == '__main__':
    main()


