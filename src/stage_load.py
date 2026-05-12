import streamlit as st
import load_data as ld
import info_viewer as iv
import display_components as dc
from stutil import init_state, m, add_vertical_space, hline, unlocked, mq, set_state

def stage_load(exp, main_body_container, upper_container, message_container):
    """
    The function `stage_load` in Python defines a multi-tab interface for uploading sample data,
    consumables, and amplicons in a pipeline.
    
    Args:
      exp: An Experiment object (the currently loaded run)
      main_body_container: The `main_body_container` parameter is a container where the main content
      of the stage is displayed. This container is where the different tabs for loading samples, 
      consumables, and amplicons are rendered based on the selected tab in the user interface.
      upper_container: The `upper_container` parameter in the `stage_load` function is used to define
    the layout structure for the upper section of the user interface. It is divided into three
    columns (`tab_col1`, `tab_col2`, `tab_col3`) where different components or information can be
    displayed
      message_container: The `message_container` parameter in the `stage_load` function is used to
    specify the container where messages or notifications related to the loading stage will be
    displayed. This container can be used to show information, warnings, errors, or any other messages
    that are relevant to the user during the loading process.
    """
    tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    with tab_col1:
        #load_data_tab = dc.create_tabs([("Load Samples", "Rodentity or custom samples for the complete pipeline"),
        #        ("Load Consumables", "Additional plates and information"),
        #        ("Load Amplicons","Additional amplicons to be sequenced")])
        load_data_tab = dc.create_tabs([("Load Samples", ""),("Load Consumables", ""),
                ("Load Amplicons", "")])
    if not load_data_tab:
        init_state('load_tab', 1)
        load_data_tab = st.session_state['load_tab']
    #------------------------------------ Load ~ TAB 1: Load sample data  --------------------------------------
    if load_data_tab == 1:
        st.session_state['load_tab'] = 1
        if unlocked(exp):
            init_state('run queue', [])
            with main_body_container:
                st.subheader('Upload Sample Files')
                st.markdown('**Upload Rodentity plate files, custom manifests, or 384-plate reruns, '+\
                        'then move to the next tab to upload pipeline consumables and other important files**')
                with st.expander('Upload Rodentity Ear Punch Plates', expanded=True):
                    ld.load_rodentity_data('rodentity_load1')
                    ld.assign_rodentity_dna_plate('rodentity_load2')
                hline()
                ld.load_custom_manifests('custom_load1')
                hline()
                ld.load_manifest_384('manifest_384_load1')
                hline()

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_load', default_view1='Samples',
                default_view2='Files')

    #------------------------------------ Load ~ TAB 2: Load consumables ---------------------------------------
    if load_data_tab == 2:
        st.session_state['load_tab'] = 2
        if unlocked(exp):
            init_state('upload stage', None)
            with main_body_container:
                st.subheader('Upload Consumables')
                st.markdown('**Upload consumables, PCR files, and custom volumes, '+\
                        'then move to the next tab to upload amplicon files**')
                ld.load_extra_consumables('consumables_load2')
                ld.load_pcr1_files('pcr1_load2')
                ld.load_pcr2_files('pcr2_load2')
                add_vertical_space(1)
                st.subheader('Custom Volumes')
                ld.custom_volumes('custom1')
                hline()
        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_consumables', default_view1='Consumables',
                default_view2='Samples')

    #------------------------------------ Load ~ TAB 3: Load amplicons ---------------------------------------
    if load_data_tab == 3:
        st.session_state['load_tab'] = 3
        if unlocked(exp):
            init_state('upload stage', None)
            with main_body_container:
                st.subheader('Upload Amplicon Files')
                st.markdown('**Upload amplicon plates and reference files**')
                ld.load_amplicons('amp_load3')
                add_vertical_space(1)
                ld.load_amplicon_references('ref_load3')
                hline()
        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_consumables', default_view1='Files',
                default_view2='Samples')
