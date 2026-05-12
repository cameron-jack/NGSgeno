import streamlit as st
import load_data as ld
import info_viewer as iv
from stutil import init_state, m, add_vertical_space, hline, unlocked, mq, set_state
import display_components as dc
from generate import run_generate

def stage_miseq(exp, main_body_container, upper_container, message_container):
    tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    info_holder = st.container()

    with tab_col1:
        with tab_col1:
            #miseq_tab = dc.create_tabs([("Download", "Miseq Samplesheet"), ("Upload", "Miseq Sequence Files")])
            miseq_tab = dc.create_tabs([("Download", "Miseq Samplesheet")])
        if not miseq_tab:
            init_state('miseq_tab', 1)
            miseq_tab = st.session_state['miseq_tab']

    add_vertical_space(1)

    #-------------------------------- Miseq ~ TAB 1: Download Miseq Samplesheet --------------------------------
    if miseq_tab == 1:
        st.session_state['miseq_tab'] = 1
        with main_body_container:
            _,header_col,_ = st.columns([2,2,1])

            if exp.locked:
                st.warning(f'Experiment {exp.name} locked from further modification')

            with header_col:
                st.subheader('Download MiSeq File')
                add_vertical_space(1)
            #hline()

            #ld.load_rodentity_references('reference_miseq1')
            if exp.get_miseq_samplesheets():
                dc.get_miseq_download_btn(exp)
                add_vertical_space(4)

            else:
                st.warning(f'No MiSeq Samplesheet available for download')

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_miseq1', default_view1='Files',
                default_view2='Plates')

    #------------------------ Miseq ~ TAB 2: Upload Reference File and Miseq Sequences -------------------------
    if miseq_tab == 2:
        st.session_state['miseq_tab'] = 2
        with main_body_container:
            st.subheader('Upload Custom Reference Files')
            ld.load_rodentity_references('reference_miseq2')
            add_vertical_space(2)

            caller_id = 'seq_upload_ready'
            success = exp.check_sequence_upload_ready(caller_id)
            if caller_id in mq:
                for msg, lvl in mq[caller_id]:
                    m(msg, level=lvl, no_log=True)
                sleep(0.3)
            mq[caller_id] = set()
            if not success:
                st.warning('Resources are required for allele calling and must be present before FASTQs can be uploaded')
            else:
                ld.load_miseq_fastqs('miseq_tab1')


        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_miseq2', default_view1='Files',
                default_view2='Plates')
