#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New streamlit interactive web interface
@last_edit: 2022-08-29
@edit_comment:

A web based interactive GUI with Streamlit. Plate and sample barcodes here are unguarded.
In all other code they must be guarded. We guard them here before we send them to any external function.
"""

import os
from pathlib import PurePath
from math import fabs, factorial, floor, ceil
from io import StringIO
from shutil import copy2

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
import display_components as dc
import sidebar_components as sb


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




def st_radio_horizontal(*args, **kwargs):
    """Trick to have horizontal st radio to simulate tabs"""
    col, _ = st.columns([3,1])
    with col:
        st.write('<style> div[data-testid=column] > div > div > div > div.stRadio > div{flex-direction: row}</style>', unsafe_allow_html=True)
        return st.radio(*args, **kwargs)




def clear_widget(*keys):
    for k in keys:
        if k in st.session_state:
            st.session_state[k] = ''


def pipeline_sb():
    """
    Pipeline page
    Sidebar
    Drop down box for all of the option for the stages of the pipeline
    """
    pipeline_stage = st.sidebar.selectbox("Pipeline Stage", ("Load Data", "1. Nimbus Plates", "2. Echo Primers (PCR 1)", "3. Echo Indexing (PCR 2)", 
                                                             "4. Miseq", "5. Genotyping", "6. Review"))
    if pipeline_stage == "Load Data":
        st.session_state['pipe_stage'] = 1

    if pipeline_stage == "1. Nimbus Plates":
        st.session_state['pipe_stage'] = 2

    if pipeline_stage == "2. Echo Primers (PCR 1)":
        st.session_state['pipe_stage'] = 3

    if pipeline_stage == "3. Echo Indexing (PCR 2)":
        st.session_state['pipe_stage'] = 4

    if pipeline_stage == "4. Miseq":
        st.session_state['pipe_stage'] = 5

    if pipeline_stage == "5. Genotyping":
        st.session_state['pipe_stage'] = 6

    if pipeline_stage == "6. Review":
        st.session_state['pipe_stage'] = 7


def pipe_stages(exp, stage):
    '''
    Pipeline page

    Goes through each stage of the pipeline, depending on the session state.
    Args:
        int stage: st.session_state['pipe_stage']
    '''

    s1_st = st.container()
    s2_st = st.container()
    #s3_st = st.container()
    s4_st = st.container()
    
    
    if stage == 1:
        assay_usage = exp.get_assay_usage()
        data_tab1, data_tab2 = s1_st.tabs(['Load Data', 'View Data'])


        with data_tab1:
            st.markdown('<h5 style="text-align:center;color:#2BA2D0">Table Summary</h5>', unsafe_allow_html=True)
            dc.data_table(key=1)
            ld.load_rodentity_data()
            ld.load_custom_csv()
            ld.load_database_data()
            pcr_files_exp = st.expander('Upload PCR files')
            with pcr_files_exp:
                ld.upload_pcr_files()

        with data_tab2:
            table_option = st.selectbox('View Data options',
                                ('Summary', 'Edit Table', 'Plate View', 'Explore Assays', 'View Log'))

            st.markdown(f'<h5 style="text-align:center;color:#2BA2D0"> Selected: {table_option} </h5>', unsafe_allow_html=True)
            dc.data_table(key=2)

            pcr_components_exp = st.expander('Required Components for PCR Reactions', expanded=False)
            with pcr_components_exp:
                dc.display_pcr_components(assay_usage)

            dc.display_primer_components(assay_usage)

    #Nimbus: 96-wells of sample plates to 384-wells of DNA plates
    if stage == 2:
        s2_st.markdown('<h2 style="text-align:center;color:#0f6b8e">Nimbus</h2>', unsafe_allow_html=True)
            
        if not st.session_state['experiment'].dest_sample_plates:
            s2_st.markdown('<p style="color:#FF0000">Load data inputs to enable Nimbus input file generation.</p>', unsafe_allow_html=True)
        else:
            s2_st.markdown('<h5 style="text-align:center;color:#2BA2D0">Nimbus input files have been generated from your data</h5>', unsafe_allow_html=True)
            #st.write('Combine 96-well sample plates to 384-well DNA plate.')
            s2_st.write('')
            s2_st.write('')
            nim_tab1, nim_tab2, view_data_tab = s2_st.tabs(["Download", "Upload", "View Data"])
            
            
            #Nimbus files, Echo files?, xbcs    
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()

            # do we have any Nimbus inputs to generate + download
            yet_to_run = len(st.session_state['experiment'].dest_sample_plates) - len(nfs)
            if yet_to_run and yet_to_run > 0: 
                plates_to_run = [dest_plate for dest_plate in exp.dest_sample_plates if all([dest_plate not in nf for nf in nfs])]
                plates_to_run_str = '\n'.join(plates_to_run)
                success = st.session_state['experiment'].generate_nimbus_inputs()
                if not success:
                    nim_tab1.markdown('<p style="color:#FF0000">Failed to generate Nimbus files. Please read the log.</p>', unsafe_allow_html=True)
                else:
                    nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()

                #nim_tab1.write(plates_to_run_str)
                nim_tab2.write(f"{yet_to_run} 96-well plate sets need Nimbus input file generation")

            if nfs:
                nim_tab1.markdown('<h5 style="text-align:center">1. Download Nimbus input files</h5>', unsafe_allow_html=True)
                nim_tab1.markdown('<p style="text-align:center">Download the 96-well sample plates as input files for Nimbus to run.</p>', unsafe_allow_html=True)
                nim_tab1.write('')
                nim_tab1.write('')
                _,dl_col1,dl_col2,dl_col3,dl_col4,_=nim_tab1.columns([1,9,6,9,6,1])
                
                #Generate file names to download + download buttons
                for i,nf in enumerate(nfs):
                    file_name=nf.split('\\')[1]

                    if (i+1) % 2 != 0:
                        dl_col1.markdown(f"<p style='text-align:left;color:#4b778c;padding:5px'>{file_name}</p>", unsafe_allow_html=True)

                        dl_col2.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
                
                    else:
                        dl_col3.markdown(f"<p style='text-align:left;color:#4b778c;padding:5px''>{file_name}</p>", unsafe_allow_html=True)
                 
                        dl_col4.download_button("Download ", open(nf), file_name=PurePath(nf).name, 
                                key='nimbus_input'+str(i), help=f"Download Nimbus input file {nf}")
              
                        
            #?
            nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
            print (f"{nfs=} {efs=} {xbcs=}")
            # Nimbus upload
            if not efs and not xbcs:
                nim_tab2.markdown('<h4 style="text-align=center">Awaiting Nimbus output files</h4>', unsafe_allow_html=True)
            elif efs and not xbcs:
                nim_tab2.markdown('<h4 style="text-align:center">All Nimbus outputs now uploaded</h4>', unsafe_allow_html=True) 
            elif xbcs:
                nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
                missing_nims = ['Echo_384_COC_0001_'+xbc+'_0.csv' for xbc in xbcs]
                mn_str = '</br>'.join(missing_nims)

                nim_tab2.markdown('<h5 style="text-align:center">2. Upload Nimbus output files</h5>', unsafe_allow_html=True)
                nim_tab2.markdown('<p style="text-align:center">Upload the resulting 384-well DNA plate files for Echo.</p>', unsafe_allow_html=True)
                # list nimbus outputs that we still need
                nim_tab2.markdown('<h6 style="color:#8e0f5b">Nimbus output files to upload:</h6>', unsafe_allow_html=True)
                nim_tab2.markdown(f"<p style='color:#87adc7'>{mn_str}</p>", unsafe_allow_html=True)

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
                #else:
                #    st.session_state['nim_uploads'] = ''
                                    
                    
        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
        if not xbcs:
            s2_st.markdown('<h4 style="color:#FF0000;text-align:center">Ready to run Echo primer stage</h4>', unsafe_allow_html=True)


        with view_data_tab:
            st.markdown(f'<h5 style="text-align:center;color:#2BA2D0"> Table Summary </h5>', unsafe_allow_html=True)
            dc.data_table(key=3)

    if stage == 3:
        #Primer stage
        s3_st = st.container()
        s3_title = s3_st.title('')
        s3_checklist_exp = s3_st.expander('Checklist for Plates', expanded=False)
        echo1_dna_col, echo1_pcr_col, echo1_taqwater_col = s3_checklist_exp.columns(3)
        s3_st.write('')
        s3_st.write('')
        comp_tab, pcr1_upload_tab, pcr1_download_tab = s3_st.tabs(["PCR Components", "Provide Plates", "Download Picklists"])
        pcr1_comp1_tab_title = comp_tab.title('')
        pcr1_upload_tab_title = pcr1_upload_tab.title('')
        pcr1_upload_cont = pcr1_upload_tab.container()
        pcr1_download_tab_title = pcr1_download_tab.title('')
        
        nfs, efs, xbcs = st.session_state['experiment'].get_nimbus_filepaths()
        missing_nims = ['Echo_384_COC_0001_'+xbc+'_01.csv' for xbc in xbcs]
        available_nimbus = ['Echo_384_COC_0001_'+ef+'_01.csv' for ef in efs]   

        if available_nimbus:
            included_DNA_plates = set()
            included_PCR_plates = set()
            included_taqwater_plates = set()

            assay_usage = st.session_state['experiment'].get_assay_usage()

            # Arrange available DNA, PCR, and taq+water plates in columns
            
            #s3_st.markdown('<h6 style="color:#0f6b8e;text-align:center">Components for PCR 1</h6>', unsafe_allow_html=True)
            
            with pcr1_upload_tab:
                ld.upload_pcr_files(1)
             
            #want this to be below the rest of the info but not sure if that's possible
            #comp_tab.markdown('<h6 style="color:#0f6b8e;text-align:center">Checklist of plates for Echo (PCR 1 Assay)</h6>', unsafe_allow_html=True)
            #comp_tab.write('')
            #choose_plates_exp = comp_tab.expander('Plate Selection', expanded=False)
            
            echo1_dna_col.markdown('**DNA Plates**')
            echo1_pcr_col.markdown('**PCR Plates**')
            echo1_taqwater_col.markdown('**Taq+Water Plates**')
            #Glitches when clicking - fix
            for nim in available_nimbus:
                #Naming system - should look better
                echo_filename=nim.split('\\')[1].split('.')[0]
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


            #components
            


            #upload buttons


            #download buttons

            #st.write(f"{included_DNA_plates=}")
            if included_DNA_plates:
                s3_st.write('')
                
                assay_usage = st.session_state['experiment'].get_assay_usage(dna_plate_list=included_DNA_plates)


                with comp_tab:
                    dc.display_pcr_components(assay_usage, PCR_stage=1)
                    dc.display_primer_components(assay_usage)
                
                if st.session_state['experiment'].check_ready_echo1(included_DNA_plates, included_PCR_plates, included_taqwater_plates):

                    pcr1_download_tab.markdown('<h5 style="color:#95B7C3;text-align:center">Ready to run Echo PCR 1</h5>', unsafe_allow_html=True)
                    _,picklist_button_col,_ = pcr1_download_tab.columns([2, 2, 1])

                    echo_picklist_go = picklist_button_col.button('Generate Echo Picklist', key='echo_pcr1_go_button')

                    if echo_picklist_go:
                        success = st.session_state['experiment'].generate_echo_PCR1_picklist_interface(included_DNA_plates,
                                included_PCR_plates, included_taqwater_plates)
                        if success:
                            dna_picklist_paths, primer_picklist_paths, taqwater_picklist_paths =\
                                    st.session_state['experiment'].get_echo_PCR1_picklist_filepaths()
                            echo1_dna_download_col, echo1_primer_download_col, echo1_taqwater_download_col = st.columns(3)
                            if not dna_picklist_paths:
                                echo1_dna_download_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No DNA picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for dpp in dna_picklist_paths:
                                    echo1_dna_download_col.download_button(label=f"Download {dpp}", 
                                            data=open(dpp, 'rt'), file_name=dpp, mime='text/csv', key='dna_download_'+dpp)
                            if not primer_picklist_paths:
                                echo1_primer_download_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No primer picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for ppp in primer_picklist_paths:
                                    echo1_primer_download_col.download_button(label=f"Download {ppp}", 
                                            data=open(ppp, 'rt'), file_name=ppp, mime='text/csv', key='primer_download_'+dpp)
                            if not taqwater_picklist_paths:
                                echo1_taqwater_download_col.markdown('<h4 style="color:#ff0000;text-align:center">'+\
                                        'No taq+water picklist available</h4>', unsafe_allow_html=True)
                            else:
                                for tpp in taqwater_picklist_paths:
                                    echo1_taqwater_download_col.download_button(label=f"Download {tpp}", 
                                            data=open(tpp, 'rt'), file_name=tpp, mime='text/csv', key='taqwater_download_'+tpp)

                #display_reaction_stats(assay_usage)

                else:
                    pcr1_download_tab.markdown('<h5 style="color:#B72572;text-align:center">Please include at least one '+\
                            'DNA plate (Echo file) to carry on with pipeline</h5>', unsafe_allow_html=True)
            #extra_data()

        s3_title.markdown('<h2 style="text-align:center;color:#0f6b8e">Echo Primers</h2>', unsafe_allow_html=True)
        #comp_tab, pcr1_upload_tab, pcr1_download_tab = s3_st.tabs(["PCR Components", "Upload Files", "Download Picklists"])
        pcr1_comp1_tab_title.markdown('<h4 style="color:#0f6b8e;text-align:left">Components</h4>', unsafe_allow_html=True)
        pcr1_upload_tab_title.markdown('<h4 style="color:#0f6b8e;text-align:left">Provide Plate Information</h4>', unsafe_allow_html=True)

        pcr1_upload_cont.write('Provide primer plate layout and volumes, taq and water barcodes and PCR plate barcodes.')
        pcr1_upload_cont.write('')
        #pcr1_upload_tab.markdown('<h4 style="color:#0f6b8e;text-align:center">Files and barcodes are required for taq and water plates, primer plates</h4>', unsafe_allow_html=True)
        pcr1_download_tab_title.markdown('<h4 style="color:#0f6b8e;text-align:center">Download the Picklists for Echo</h4>', unsafe_allow_html=True)
        #primer_plate_detail_exp = s3_st.expander('Primer Plate Details')

        
        #Missing Nimbus Output Files - add this warning in previous stage? Reduces space...
        #if missing_nims:
        #    s3_st.markdown('<h6 style="text-align:left;color:#B92A5D">Missing the Nimbus Output Files:</h6>', unsafe_allow_html=True)
        #    for mn in missing_nims:
        #        s3_st.markdown(f'<p style="text-align:left">{mn}</p>', unsafe_allow_html=True)
        #    s3_st.markdown('<p style="text-align:left;color:#B92A5D">Upload now or continue without them</p>', unsafe_allow_html=True)
        #    s3_st.write('')
  
        #    if 'nim_upload' not in st.session_state:
        #            st.session_state['nim_upload'] = ''
        #    nim_outputs = s3_st.file_uploader('Upload Nimbus output files', type='csv', 
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
        s4_title = s4_st.title('')
        s4_checklist_expander = s4_st.expander('Checklist for Plates', expanded=False)
        echo2_pcr_col, echo2_index_col, echo2_taqwater_col, echo2_amplicon_col = s4_checklist_expander.columns(4)

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



        

        pcr2_comp_tab, pcr2_upload_tab, pcr2_download_tab = s4_st.tabs(["PCR 2 Components", "Upload Files", "Download Picklists"])

        pcr2_comp_tab.markdown('<h4 style="color:#0f6b8e;text-align:left">Components</h4>', unsafe_allow_html=True)
        pcr2_upload_tab.markdown('<h4 style="color:#0f6b8e;text-align:left">Provide Plate Information</h4>', unsafe_allow_html=True)
        pcr2_download_tab.markdown('<h4 style="color:#0f6b8e;text-align:center">Download the Picklists for Indexing (Echo stage 2)</h4>', unsafe_allow_html=True)

        assay_usage = st.session_state['experiment'].get_assay_usage()

        with pcr2_comp_tab:
            dc.display_pcr_components(assay_usage, PCR_stage=2)
            #dc.display_primer_components(assay_usage)

        with pcr2_upload_tab:
            ld.upload_pcr_files(2)

        with pcr2_download_tab:
            if st.session_state['experiment'].check_ready_echo2(included_PCR_plates, included_taqwater_plates, included_index_plates):
                pcr2_download_tab.markdown('<h5 style="color:#95B7C3;text-align:center">Ready to run index picklist generation (Echo PCR 2)</h5>', unsafe_allow_html=True)
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
                pcr1_download_tab.markdown('<h5 style="color:#B72572;text-align:center">Please include at least one PCR plate, '+\
                        'one index plate and one taq+water plate to carry on with pipeline</h5>', unsafe_allow_html=True)
            


        s4_title.markdown('<h2 style="text-align:center;color:#0f6b8e">Echo Index</h2>', unsafe_allow_html=True)
        
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
        for fp in exp.get_miseq_samplesheets():
            miseq_col1.markdown(f'<strong style="color:#486e7a">{fp}</strong>', unsafe_allow_html=True)
            download_miseq = miseq_col2.button(label='Download', key='dnld_samplesheet_'+str(fp))

        with s5_seq_tab:
            ld.upload_miseq_fastqs(exp)

    if stage == 6:
        # Allele calling
        s4_title = s4_st.title('Allele calling')
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
                dc.data_table(key=2, view=table_option)
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
                dc.data_table(key=2, view=table_option)
                st.session_state['test_status'] = None


def main():
    st.set_page_config(page_icon="ngsg_icon.png", page_title="NGS Genotyping Pipeline", initial_sidebar_state="expanded", layout="wide")
    # Remove whitespace from the top of the page and sidebar
    st.markdown(
        """
        <style>
               .css-18e3th9 {
                    padding-top: 0rem;
                    padding-bottom: 10rem;
                    padding-left: 5rem;
                    padding-right: 5rem;
                }
                header {visibility: hidden;}
                footer {visibility: hidden;}
        </style>
        """, unsafe_allow_html=True)
    padding = 0
    st.markdown(f""" <style>
    .reportview-container .main .block-container{{
        padding-top: {padding}rem;
        padding-right: {padding}rem;
        padding-left: {padding}rem;
        padding-bottom: {padding}rem;
    }} </style> """, unsafe_allow_html=True)
    #This is a workaround from https://www.codegrepper.com/code-examples/python/streamlit+sidebar+width
    st.markdown(
        """
        <style>
        div.block-container {padding-top:0rem;}
        [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
            width: 400px;
            top: 0rem;
            padding-top: 0rem;
            height: 100vh;
        }
        [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
            width: 400px;
            top: 0rem;
            padding-top: 0rem;
            margin-left: -400px;
        }
         
        </style>
        """,
        unsafe_allow_html=True)
    
    st.sidebar.title('NGS Genotyping Pipeline')
    st.sidebar.image('ngsg_explorer.png')
    new_container = st.container()

    if 'experiment' not in st.session_state:
        st.session_state['experiment'] = None
        print('Current experiment set to clear')
     
    if 'experiment' in st.session_state:
        if 'navigation' not in st.session_state:
            st.session_state['navigation'] = 'load'

        if st.session_state['experiment']:
            st.session_state['navigation'] = 'pipe'

            exp = st.session_state['experiment']
            st.sidebar.header('Current experiment: '+ exp.name)

            pipeline_sb()
            title='Experiment '+ st.session_state['experiment'].name
            new_container.markdown(f'<h4 style="color:#BED2D6">{title}</h2>', unsafe_allow_html=True)

            if 'pipe_stage' not in st.session_state:
                st.session_state['pipe_stage'] = 1
            
            #if 'lock' not in st.session_state:
            #    st.session_state['lock'] = False
            #if 'unlock' not in st.session_state:
            #    st.session_state['unlock'] = False
            #if exp.locked:
            #    #unlock_button = st.sidebar.button('Unlock experiment', help='Allow modification of experiment again')
            #    #if unlock_button and not st.session_state['unlock']:
            #    #    st.session_state['unlock'] = True
            #    pass
            #else:
            #    #lock_button = st.sidebar.button('Lock experiment', help='Freezes the current experiment, preventing further modification')
            #    #if lock_button and not st.session_state['lock']:
            #    #    st.session_state['lock'] = True
            #    pass
            #if st.session_state['unlock']:
            #    exp.unlock()
            #    st.session_state['unlock'] = False
            #if st.session_state['lock']:
            #    exp.lock()
            #    st.session_state['lock'] = False

            ##if 'nuked' not in st.session_state:
            ##    st.session_state['nuked'] = False
            ##nuke_button = st.sidebar.button('Delete experiment', help='Hides the current experiment from the user interface. Manual retrieval is required')
            ##if nuke_button and not st.session_state['nuked']:
            ##    st.session_state['nuked'] = True
            ##if st.session_state['nuked']:
            ##    st.sidebar.write('Currently not enabled, but in future will hide experiment')
            ##    st.session
            #else: # pipeline stages

            #    pipeline_sb()

            ##Pipeline Stages

            if 'pipe_stage' not in st.session_state or st.session_state['pipe_stage'] == 1:
                pipe_stages(exp, 1)
                st.sidebar.write('')
                st.sidebar.write('')
                st.sidebar.markdown('<p style="color:#8AB1BD">Or change experiment</p>', unsafe_allow_html=True)

                #gotta change logic so that this load still
                sb.folder_sb()
        
            #Nimbus
            if st.session_state['pipe_stage'] == 2:
                pipe_stages(exp, 2)
                   
            #Echo primers   
            if st.session_state['pipe_stage'] == 3:
                pipe_stages(exp, 3)

  
            if st.session_state['pipe_stage'] == 4:  # echo barcodes
                pipe_stages(exp, 4)

            if st.session_state['pipe_stage'] == 5:  # Miseq
                pipe_stages(exp, 5)

            if st.session_state['pipe_stage'] == 6:  # Genotyping
                pipe_stages(exp, 6)

            if st.session_state['pipe_stage'] == 7:  # Review
                pipe_stages(exp, 7)


    if not st.session_state['experiment'] or (st.session_state['experiment'] and st.session_state['navigation'] == 'load'):
        #st.session_state['experiment'] = None
        new_container.write('')
        new_container.write('')
        new_container.markdown('<h2 style="text-align:center;color:#154c79">Please load existing experiment or create a new one</h2>', unsafe_allow_html=True)

        sb.folder_sb()


    
       
    #folder_sb()



        

    #hvar = """  <script>
    #                var elements = window.parent.document.querySelectorAll('.streamlit-expanderHeader');
    #                elements.forEach(function (element) {
    #                    element.style.fontSize = 'large';
    #                    });         
    #            </script>"""
    #components.html(hvar, height=0, width=0)

    

    #st.markdown(
    #    """
    #<script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    #<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    #<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    #""",
    #    unsafe_allow_html=True,
    #)


if __name__ == '__main__':
    main()