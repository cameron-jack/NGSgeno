"""
@created: 20/09
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: Seperate file for pipeline
@last_edit: 2022-09-20
@edit_comment:

Function for each stage of the pipeline for the UI
"""

#from asyncio.windows_utils import pipe
from asyncio import selector_events
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
from bin.ngsmatch import match_alleles

import load_data as ld
from display_components import data_table, display_pcr_components, display_primer_components

credits="""
@created: March 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: New interactive web application interface
@last_edit: 
@edit_comment: 

Web application style interface using Streamlit. Present information as both a dashboard and a workflow in one page.

The GUI interacts with a single Experiment object at one time. Methods are called on this to activate pipeline
functionality. The Experiment then deals directly with the pipeline logic.
"""

def pipe_stages(stage):
    '''
    Goes through each stage of the pipeline, depending on the session state.
    Args:
        int stage: st.session_state['pipe_stage']
    '''

    pipe_container = st.container()
    pipe_title = pipe_container.title('')
    pipe_subtitle = pipe_container.empty()
    pipe_subtitle = pipe_container.empty()
    
    
    if stage == 1:
        
        pipe_title.markdown('<h2 style="text-align:center;color:#0f6b8e">Load Data</h2>',\
                 unsafe_allow_html=True)
        load_data_tab, view_data_tab = pipe_container.tabs(['Load Data', 'View Data'])

        assay_usage = st.session_state['experiment'].get_assay_usage()

        with load_data_tab:            
            ld.load_rodentity_data()
            ld.load_custom_csv()
            ld.load_database_data()
            #data_table(key='1a')
            # pcr_files_exp = st.expander('Upload PCR files')
            # with pcr_files_exp:
            #     ld.upload_pcr_files()

        with view_data_tab:
            data_table(key='1b', options=True)
            pcr_components_exp = st.expander('Required Components for PCR Reactions',\
                         expanded=False)
            with pcr_components_exp:
                display_pcr_components(assay_usage)

            display_primer_components(assay_usage)

    #Nimbus: 96-wells of sample plates to 384-wells of DNA plates
    if stage == 2:
        pipe_title.markdown('<h2 style="text-align:center;color:#0f6b8e">Nimbus</h2>', unsafe_allow_html=True)
            
        if not st.session_state['experiment'].dest_sample_plates:
            pipe_subtitle.markdown('<h5 style="text-align:center;color:#FF0000">'+\
                    'Load data inputs to enable Nimbus input file generation.</h5>', unsafe_allow_html=True)
        else:
            pipe_subtitle.markdown('<h5 style="text-align:center;color:#2BA2D0">'+\
                    'Nimbus input files have been generated from your data</h5>', unsafe_allow_html=True)
            
            nim_tab1, nim_tab2, view_data_tab = pipe_container.tabs(["Download", "Upload", "View Data"])
            nim_tab1_title = nim_tab1.empty()
            nim_tab2_title = nim_tab2.empty()
            
            #nfs = Nimbus file paths, efs = Echo file paths, xbcs = barcodes    
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()

            # do we have any Nimbus inputs to generate + download
            yet_to_run = len(st.session_state['experiment'].dest_sample_plates) - len(nfs)
            if yet_to_run and yet_to_run > 0: 
                plates_to_run = [dest_plate for dest_plate in exp.dest_sample_plates if all([dest_plate not in nf for nf in nfs])]
                plates_to_run_str = '\n'.join(plates_to_run)
                success = st.session_state['experiment'].generate_nimbus_inputs()
                if not success:
                    nim_tab1.markdown('<p style="color:#FF0000">Failed to generate Nimbus files. Please read the log.</p>',\
                             unsafe_allow_html=True)
                else:
                    nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()

                #nim_tab1.write(plates_to_run_str)
                nim_tab2.markdown(f'<p style="text-align:center">{yet_to_run} 96-well plate sets need Nimbus input file generation</p>',\
                            unsafe_allow_html=True)

            if nfs:
                nim_tab1_title.markdown('<h4 style="color:#0f6b8e;text-align:center">Download Nimbus Files</h4>',\
                         unsafe_allow_html=True)
                nim_tab1.markdown('<p style="text-align:center">Download the 96-well sample plates as input files for Nimbus to run.</p>',\
                             unsafe_allow_html=True)
                nim_tab1.write('')
                nim_tab1.write('')
                _,dl_col1,dl_col2,dl_col3,dl_col4,_=nim_tab1.columns([1,9,6,9,6,1])
                
                #Generate file names to download + download buttons
                for i,nf in enumerate(nfs):
                    file_name=nf.split('/')[1]

                    if (i+1) % 2 != 0:
                        dl_col1.markdown(f"<p style='text-align:left;color:#4b778c;padding:5px'>{file_name}</p>",\
                                 unsafe_allow_html=True)

                        dl_col2.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                
                    else:
                        dl_col3.markdown(f"<p style='text-align:left;color:#4b778c;padding:5px''>{file_name}</p>",\
                                 unsafe_allow_html=True)
                 
                        dl_col4.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
              
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            print (f"{nfs=} {efs=} {xbcs=}")
            # Nimbus upload
            if not efs and not xbcs:
                nim_tab2_title.markdown('<h4 style="color:#0f6b8e;text-align=center">Awaiting Nimbus output files</h4>',\
                         unsafe_allow_html=True)
            elif efs and not xbcs:
                pipe_subtitle.markdown('<h4 style="color:#FF0000;text-align:center">Ready to run Echo primer stage</h4>',\
                         unsafe_allow_html=True)
                nim_tab2_title.markdown('<h4 style="color:#0f6b8e;text-align:center">All Nimbus outputs now uploaded</h4>',\
                         unsafe_allow_html=True)
                nim_tab2.markdown('')
                _,file_col1, file_col2,_ = nim_tab2.columns([1,2,2,1])
                for i, ef in enumerate(efs):
                    file_name = ef.split('/')[1]
                    if i % 2 == 0:
                        file_col2.markdown(f'<p style="text-align:center;color:#17754d;padding:5px">{file_name}</p>',\
                                 unsafe_allow_html=True)
                    else:
                        file_col1.markdown(f'<p style="text-align:center;color:#17754d; padding:5px">{file_name}</p>',\
                                 unsafe_allow_html=True)
                    #nim_tab2.markdown(f'<p style="text-align:center;color:#17754d">{file_name}</p>', unsafe_allow_html=True)
            elif xbcs:
                nim_tab2_title.markdown('<h4 style="color:#0f6b8e;text-align:center">Upload Nimbus Output Files</h4>',\
                         unsafe_allow_html=True)
                #nim_tab2.markdown('<p style="text-align:center">Upload the resulting 384-well DNA plate files for Echo.</p>', unsafe_allow_html=True)
                # list nimbus outputs that we still need
                #nim_tab2.markdown('<h6 style="color:#8e0f5b">Nimbus output files to upload:</h6>', unsafe_allow_html=True)
                nim_tab2_missing = nim_tab2.empty() 

                if 'nim_upload' not in st.session_state:
                    st.session_state['nim_upload'] = ''
                nim_outputs = nim_tab2.file_uploader(' ', type='csv', 
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


                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                missing_nims = ['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs]
                mn_str = '</br>'.join(missing_nims)

                nim_tab2_missing.markdown(f"<p style='color:#87adc7;text-align:center'>{mn_str}</p>",\
                         unsafe_allow_html=True)
                #else:
                #    st.session_state['nim_uploads'] = ''
            with view_data_tab:
                data_table(key='2a', options=True)                       
                    
        # nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
        # if not xbcs:
            

    if stage == 3:
        #Primer stage
        pipe_title.markdown('<h2 style="text-align:center;color:#0f6b8e">Echo Primers</h2>',\
                 unsafe_allow_html=True)
        s3_checklist_exp = pipe_container.expander('Checklist for Plates', expanded=False)
        echo1_dna_col, echo1_pcr_col, echo1_taqwater_col = s3_checklist_exp.columns(3)
        pipe_container.write('')
        pipe_container.write('')
        comp_tab, pcr1_upload_tab, pcr1_picklist_tab, view_data_tab = \
                    pipe_container.tabs(["PCR Components", "Provide Plates", "Generate Picklists", "View Data"])
        comp_tab.markdown('<h4 style="color:#0f6b8e;text-align:left">Components</h4>',\
                 unsafe_allow_html=True)
        pcr1_upload_tab.markdown('<h4 style="color:#0f6b8e;text-align:left">Provide Plate Information</h4>',\
                 unsafe_allow_html=True)
        pcr1_upload_tab.write('Provide primer plate layout and volumes, taq and water barcodes and PCR plate barcodes.')
        pcr1_upload_tab.write('')
        pcr1_picklist_tab.markdown('<h4 style="color:#0f6b8e;text-align:center">Picklists for Echo</h4>',\
                 unsafe_allow_html=True)
        
        with view_data_tab:
            data_table(key = '3a', options=True)

        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
        missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
        available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]   

        if available_nimbus:
            included_DNA_plates = set()
            included_PCR_plates = set()
            included_taqwater_plates = set()

            assay_usage = st.session_state['experiment'].get_assay_usage()

            # Arrange available DNA, PCR, and taq+water plates in columns
            
            #pipe_container.markdown('<h6 style="color:#0f6b8e;text-align:center">Components for PCR 1</h6>', unsafe_allow_html=True)
            
            with pcr1_upload_tab:
                ld.upload_pcr_files(1)
            
            echo1_dna_col.markdown('**DNA Plates**')
            echo1_pcr_col.markdown('**PCR Plates**')
            echo1_taqwater_col.markdown('**Taq+Water Plates**')
            #Glitches when clicking - fix
            for nim in available_nimbus:
                #Naming system - should look better
                echo_filename=nim.split('/')[1].split('.')[0]
                inc_dna = echo1_dna_col.checkbox(echo_filename, value=True, key='chk_box_dna_'+nim)
                if inc_dna:
                    included_DNA_plates.add(echo_filename.split('_')[-2])
                

            for pcr_pid in st.session_state['experiment'].get_pcr_plates():
                inc_pcr = echo1_pcr_col.checkbox(file_io.unguard_pbc(pcr_pid, silent=True), value=True, key='chk_box_pcr_'+pcr_pid)
                if inc_pcr:
                    included_PCR_plates.add(pcr_pid)

            for taqwater_pid in st.session_state['experiment'].get_taqwater_plates():
                inc_taqwater = echo1_taqwater_col.checkbox(file_io.unguard_pbc(taqwater_pid, silent=True), 
                        value=True, key='chk_box_taqwater_'+taqwater_pid)
                if inc_taqwater:
                    included_taqwater_plates.add(taqwater_pid)

            #st.write(f"{included_DNA_plates=}")
            if included_DNA_plates:
                pipe_container.write('')
                
                assay_usage = st.session_state['experiment'].get_assay_usage(dna_plate_list=included_DNA_plates)


                with comp_tab:
                    display_pcr_components(assay_usage, PCR_stage=1)
                    display_primer_components(assay_usage)
                
                if st.session_state['experiment'].check_ready_echo1(included_DNA_plates, included_PCR_plates, included_taqwater_plates):

                    pcr1_picklist_tab.markdown('<h5 style="color:#95B7C3;text-align:center">Ready to run Echo PCR 1</h5>',\
                             unsafe_allow_html=True)
                    pcr1_picklist_tab.write('')
                    _,picklist_button_col,_ = pcr1_picklist_tab.columns([3, 2, 2])

                    echo_picklist_go = picklist_button_col.button('Generate Echo Picklist', key='echo_pcr1_go_button')
                    pcr1_picklist_tab.write('')
                    pcr1_picklist_tab.write('')
                    if echo_picklist_go:
                        success = st.session_state['experiment'].generate_echo_PCR1_picklist_interface(included_DNA_plates,
                                included_PCR_plates, included_taqwater_plates)
                        if success:
                            dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths =\
                                    st.session_state['experiment'].get_echo_PCR1_picklist_filepaths()
                            _, echo1_picklist_file_col, echo1_picklist_download_col, _ = pcr1_picklist_tab.columns([1,2,1,1])
                            if not dna_picklist_paths:
                                echo1_picklist_file_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No DNA picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for dpp in dna_picklist_paths:
                                    dpp_file = dpp.split('/')[1]
                                    echo1_picklist_file_col.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{dpp_file}</p>', \
                                            unsafe_allow_html=True)
                                    echo1_picklist_download_col.download_button(label=f"Download", 
                                            data=open(dpp, 'rt'), file_name=dpp, mime='text/csv', key='dna_download_'+dpp)
                            if not primer_picklist_paths:
                                echo1_picklist_file_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No primer picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for ppp in primer_picklist_paths:
                                    ppp_file = ppp.split('/')[1]
                                    echo1_picklist_file_col.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{ppp_file}</p>',\
                                             unsafe_allow_html=True)
                                    echo1_picklist_download_col.download_button(label=f"Download", 
                                            data=open(ppp, 'rt'), file_name=ppp, mime='text/csv', key='primer_download_'+dpp)
                            if not taqwater_picklist_paths:
                                echo1_picklist_file_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No taq+water picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for tpp in taqwater_picklist_paths:
                                    tpp_file = tpp.split('/')[1]
                                    echo1_picklist_file_col.markdown(f'<p style="text-align:left;color:#4b778c;padding:5px">{tpp_file}</p>',\
                                             unsafe_allow_html=True)
                                    echo1_picklist_download_col.download_button(label=f"Download", 
                                            data=open(tpp, 'rt'), file_name=tpp, mime='text/csv', key='taqwater_download_'+tpp)

                #display_reaction_stats(assay_usage)

                        else:
                            pcr1_picklist_tab.markdown('<h5 style="color:#B72572;text-align:center">Please include at least one '+\
                                    'DNA plate (Echo file) to carry on with pipeline</h5>', unsafe_allow_html=True)


        
        #Missing Nimbus Output Files - add this warning in previous stage? Reduces space...
        #if missing_nims:
        #    pipe_container.markdown('<h6 style="text-align:left;color:#B92A5D">Missing the Nimbus Output Files:</h6>', unsafe_allow_html=True)
        #    for mn in missing_nims:
        #        pipe_container.markdown(f'<p style="text-align:left">{mn}</p>', unsafe_allow_html=True)
        #    pipe_container.markdown('<p style="text-align:left;color:#B92A5D">Upload now or continue without them</p>', unsafe_allow_html=True)
        #    pipe_container.write('')
  
        #    if 'nim_upload' not in st.session_state:
        #            st.session_state['nim_upload'] = ''
        #    nim_outputs = pipe_container.file_uploader('Upload Nimbus output files', type='csv', 
        #            accept_multiple_files=True, help='You can upload more than one file at once')

        #    if nim_outputs: # and nim_outputs != st.session_state['nim_upload']:
        #        #st.session_state['nim_uploads'] = nim_outputs
        #        for nim_output in nim_outputs:
        #            #print(f"{nim_output.name=}")
        #            fp = st.session_state['experiment'].get_exp_fp(nim_output.name)
        #            print(f"Copying {fp} to experiment folder")
        #            with open(fp, 'wt') as outf:
        #                #print(nim_output.getvalue().decode("utf-8"))
        #                outf.write(nim_output.getvalue().decode("utf-8").replace('\r\n','\n'))


    if stage == 4:
        pipe_title.markdown('<h2 style="text-align:center;color:#0f6b8e">Echo Index</h2>', unsafe_allow_html=True)
        s4_checklist_expander = pipe_container.expander('Checklist for Plates', expanded=False)
        echo2_pcr_col, echo2_index_col, echo2_taqwater_col, echo2_amplicon_col = s4_checklist_expander.columns(4)

        echo2_pcr_col.markdown('**PCR Plates**')
        echo2_index_col.markdown('**Index Plates**')
        echo2_taqwater_col.markdown('**Taq/Water Plates**')
        echo2_amplicon_col.markdown('**Amplicon Plates**')

        pcr2_comp_tab, pcr2_upload_tab, pcr2_download_tab, pcr2_view_data_tab = \
                pipe_container.tabs(["PCR Components", "Provide Plates", "Generate Picklists", "View Data"])

        pcr2_comp_tab.markdown('<h4 style="color:#0f6b8e;text-align:left">Components</h4>',\
                 unsafe_allow_html=True)
        pcr2_upload_tab.markdown('<h4 style="color:#0f6b8e;text-align:left">Provide Plate Information</h4>',\
                 unsafe_allow_html=True)
        pcr2_download_tab.markdown('<h4 style="color:#0f6b8e;text-align:center">Generate the Picklists for Indexing (Echo stage 2)</h4>',\
                 unsafe_allow_html=True)

        with pcr2_view_data_tab:
            data_table(key='4a', options=True)

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

            for index_pid in st.session_state['experiment'].get_index_plates():
                inc_index = echo2_index_col.checkbox(file_io.unguard_pbc(index_pid, silent=True),
                                                    value=True, key='chk_box_index_'+index_pid)
                if inc_index:
                    included_index_plates.add(index_pid)

            for amplicon_pid in st.session_state['experiment'].get_amplicon_plates():
                amplicon_index = echo2_amplicon_col.checkbox(file_io.unguard_pbc(amplicon_pid, silent=True), 
                                                             value=True, key='chk_box_amplicon_'+amplicon_pid)
                if amplicon_index:
                    included_amplicon_plates.add(amplicon_pid)

        

            assay_usage = st.session_state['experiment'].get_assay_usage()

            with pcr2_comp_tab:
                display_pcr_components(assay_usage, PCR_stage=2)
                #display_primer_components(assay_usage)

            with pcr2_upload_tab:
                ld.upload_pcr_files(2)

            with pcr2_download_tab:
                if st.session_state['experiment'].check_ready_echo2(included_PCR_plates, included_taqwater_plates,\
                             included_index_plates):
                    pcr2_download_tab.markdown('<h5 style="color:#95B7C3;text-align:center">Ready to run index '+\
                             'picklist generation (Echo PCR 2)</h5>', unsafe_allow_html=True)
                    _,picklist_button_col,_ = pcr2_download_tab.columns([2, 2, 1])

                    echo_picklist_go = picklist_button_col.button('Generate Echo Picklists', key='echo_pcr2_go_button')

                    if echo_picklist_go:  # pcr_plates, index_plates, taq_water_plates, amplicon_plates=None
                        success = st.session_state['experiment'].generate_echo_PCR2_picklist_interface(included_PCR_plates, 
                                included_index_plates, included_taqwater_plates, included_amplicon_plates)
                        if success:
                            index_picklist_paths, taqwater_picklist_paths, amplicon_picklist_paths =\
                                    st.session_state['experiment'].get_echo_PCR2_picklist_filepaths()
                            # amplicon_picklist_paths is a list
                            echo2_index_download_col, echo2_taqwater_download_col, echo2_amplicon_download_col = st.columns(3)
                            if not index_picklist_paths:
                                echo2_index_download_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No index picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for ipp in index_picklist_paths:
                                    echo2_index_download_col.download_button(label=f"Download {ipp}", data=open(ipp, 'rt'), 
                                            file_name=ipp, mime='text/csv', key='index_download_'+ipp)
                            if not amplicon_picklist_paths:
                                echo2_amplicon_download_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No amplicon picklists available</h4>', unsafe_allow_html=True)
                            else:
                                for app in amplicon_picklist_paths:
                                    echo2_amplicon_download_col.download_button(label=f"Download {app}", data=open(app, 'rt'),
                                            file_name=app, mime='text/csv', key='amplicon_download_'+app)
                            if not taqwater_picklist_paths:
                                echo2_taqwater_download_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No taq+water picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for tpp in taqwater_picklist_paths:
                                    echo2_taqwater_download_col.download_button(label=f"Download {tpp}", data=open(tpp, 'rt'),
                                            file_name=tpp, mime='text/csv', key='taqwater_download_'+tpp)

                #display_reaction_stats(assay_usage)

                else:
                    pcr1_picklist_tab.markdown('<h5 style="color:#B72572;text-align:center">Please include at least one PCR plate, '+\
                            'one index plate and one taq+water plate to carry on with pipeline</h5>', unsafe_allow_html=True)
            
        
    if stage == 5:
        # Miseq
        st.title('Miseq')
        st.markdown('<h5 style="color:#2BA2D0;text-align:center">Your samples are reading for sequencing</h5>',\
                 unsafe_allow_html=True)
        st.write('')
        st.write('')
        s5_download_tab, s5_seq_tab = st.tabs(["Miseq samplesheet", "Upload sequence files"])
        s5_download_tab.markdown('<h4 style="color:#0f6b8e;text-align:center">Download Miseq samplesheet</h4>',\
                 unsafe_allow_html=True)
        s5_download_tab.write('')
        _, miseq_col1, miseq_col2, _ =  s5_download_tab.columns([2,1,1,2])
        exp = st.session_state['experiment']
        for fp in exp.get_miseq_samplesheets():
            miseq_col1.markdown(f'<strong style="color:#486e7a">{fp}</strong>', unsafe_allow_html=True)
            download_miseq = miseq_col2.button(label='Download', key='dnld_samplesheet_'+str(fp))

        with s5_seq_tab:
            ld.upload_miseq_fastqs(exp)


    if stage == 6:
        # Allele calling
        s4_title = st.title('Allele calling')
        exp = st.session_state['experiment']
        #s6_seq_tab, s6_alleles_tab = st.tabs(["Upload sequence files", "Call alleles"])

        if 'test_status' in st.session_state and st.session_state['test_status'] == 'upload':
            s6_seq_tab, s6_alleles_tab, s6_view_tab = st.tabs(["Upload sequence files", "Call alleles", "View"])
            with s6_seq_tab:
                st.write('')
                ld.upload_miseq_fastqs(exp)
                st.session_state['test_status'] = None
        elif 'test_status' in st.session_state and st.session_state['test_status'] == 'call':
            with s6_alleles_tab:
                num_unique_seq = st.number_input("Number of unique sequences per work unit", value=1)
                num_cpus = st.number_input(label="Number of processes to run simultaneously (defaults to # of CPUs)", value=os.cpu_count())
                exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                clear_lock = st.checkbox("Clear process lock")
                do_matching = st.button("Run allele calling")
                if do_matching:
                    match_alleles(exp, ncpus=num_cpus, chunk_size=num_unique_seq, exhaustive=exhaustive_mode,
                            stagefile=exp.get_exp_fp("Stage3.csv"), force_restart=clear_lock)
                    do_matching = None
                st.session_state['test_status'] = None
        elif 'test_status' in st.session_state and st.session_state['test_status'] == 'view':
            with s6_view_tab:
                st.session_state['test_status'] = 'view'
                table_option = st.selectbox('View Data options',
                        ('Summary', 'Edit Table', 'Plate View', 'Explore Assays', 'View Log'))

                st.markdown(f'<h5 style="text-align:center;color:#2BA2D0"> Selected: {table_option} </h5>', unsafe_allow_html=True)
                data_table(key=2, view=table_option)
                st.session_state['test_status'] = None
        else:
            st.session_state['test_status'] = None
            s6_seq_tab, s6_alleles_tab, s6_view_tab = st.tabs(["Upload sequence files", "Call alleles", "View data"])
            with s6_seq_tab:
                st.session_state['test_status'] = 'upload'
                st.write('')
                ld.upload_miseq_fastqs(exp)
                st.session_state['test_status'] = None

            with s6_alleles_tab:
                st.session_state['test_status'] = 'call'
                num_unique_seq = st.text_input("Number of unique sequences per work unit", "1")
                num_cpus = st.number_input(label="Number of processes to run simultaneously (defaults to # of CPUs)", value=os.cpu_count())
                exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                clear_lock = st.checkbox("Clear process lock")
                do_matching = st.button("Run allele calling")
                if do_matching:
                    match_alleles(exp, ncpus=num_cpus, chunk_size=num_unique_seq, exhaustive=exhaustive_mode,
                            stagefile=exp.get_exp_fp("Stage3.csv"), force_restart=clear_lock)
                    do_matching = None
                st.session_state['test_status'] = None

            with s6_view_tab:
                st.session_state['test_status'] = 'view'
                table_option = st.selectbox('View Data options',
                        ('Summary', 'Edit Table', 'Plate View', 'Explore Assays', 'View Log'))

                st.markdown(f'<h5 style="text-align:center;color:#2BA2D0"> Selected: {table_option} </h5>', unsafe_allow_html=True)
                data_table(key=2)
                st.session_state['test_status'] = None
