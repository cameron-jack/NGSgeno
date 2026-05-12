import streamlit as st
import load_data as ld
import info_viewer as iv
import display_components as dc
from stutil import init_state, m, mq, add_vertical_space, hline, unlocked
from generate import run_generate

def stage_nimbus(exp, main_body_container, upper_container, message_container):
    tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    with tab_col1:
        nimbus_tab = dc.create_tabs([("Download", "Nimbus input files"),("Upload", "Echo input files")])
    if not nimbus_tab:
        init_state('nimbus_tab', 1)
        nimbus_tab = st.session_state['nimbus_tab']

    nfs, efs, xbcs = exp.get_nimbus_filepaths()

    #------------------------------------ Nimbus ~ TAB 1: Download Nimbus --------------------------------------
    if nimbus_tab == 1:
        st.session_state['nimbus_tab'] = 1
        if unlocked(exp):
            with main_body_container:
                _, header_col, _ = st.columns([2,2,1])
                with header_col:
                    header_col.subheader('Generate Echo Files')
                    dc.set_nimbus_title(exp, nfs, efs)
                add_vertical_space(2)
                if not ld.check_assay_file():
                    m('Assay list file (assay to primer mapping) is required to proceed. '+\
                            'Please upload this to continue', level='error')
                    success = ld.load_assaylist('nimbus1_assaylist')
                    if success:
                        st.rerun()
                else:
                    #Generate files button
                    _, btn_col,_ = st.columns([2,2,1])
                    if st.session_state['experiment'].dest_sample_plates:
                        with btn_col:
                            run_gen_nimbus = st.button('Generate Nimbus input files', type="primary")

                        #Generate files
                        if run_gen_nimbus:
                            success = run_generate(exp, exp.generate_nimbus_inputs)
                            if not success:
                                m('Failed to generate the Nimbus files. See the log for details',
                                        level='error')
                            else:
                                add_vertical_space(2)
                                nfs, efs, xbcs = exp.get_nimbus_filepaths()

                add_vertical_space(3)
                dc.get_echo_download_buttons(nfs)

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_nimbus1', default_view1='Status',
                default_view2='Files')

    #---------------------------------- Nimbus ~ TAB 2: Upload echo input files --------------------------------
    if nimbus_tab == 2:
        st.session_state['nimbus_tab'] = 2
        if unlocked(exp):
            with main_body_container:
                _, header_col, _ = st.columns([2,2,1])

                header_col.subheader('Upload Echo Input Files')
                ld.load_echo_inputs('1')

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_nimbus2', default_view1='Status',
                default_view2='Files')
