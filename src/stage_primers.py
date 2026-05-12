import streamlit as st
import load_data as ld
import info_viewer as iv
from stutil import init_state, m, mq, add_vertical_space, hline, unlocked
import display_components as dc
from generate import run_generate, pcr1_picklists_exist
import util

def stage_primers(exp, main_body_container, upper_container, message_container):
    st.session_state['assay_filter'] = True
    pcr_stage = 1

    tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])

    #Tabs
    with tab_col1:
        primer_tab = dc.create_tabs([("PCR 1", "Components"), ("Generate", "Picklists")])
    if not primer_tab:
        init_state("primer_tab", 1)
        primer_tab = st.session_state['primer_tab']

    #nimbus fp, echo fp, barcodesnot in echo
    nfs, efs, xbcs = exp.get_nimbus_filepaths()
    missing_nims = ['Echo_384_COC_0001_'+util.unguard(xbc, silent=True)+'_0.csv' for xbc in xbcs]

    #------------------------------------ Primers ~ TAB 1: PCR 1 Components ------------------------------------
    if primer_tab == 1:
        st.session_state['primer_tab'] = 1
        if unlocked:
            with main_body_container:
                st.subheader('Primer (PCR 1) Components')
                if efs:
                    st.info('Provide the barcodes for PCR plates and Taq/water plates and '+\
                            'upload primer layouts and volumes here, then move to the *Generate Picklists* tab')
                else:
                    st.warning('Do you need to upload Echo input files (output from Nimbus)?')

                checkbox_cont = st.container()
                hline()
                add_vertical_space(1)

                display_cont = st.container()
                hline()

                add_vertical_space(1)
                st.subheader('Add Barcodes', help='Add barcodes for plates')
                pcr_col, taqwater_col = st.columns(2)

                st.subheader('Upload Files')
                ld.load_pcr1_files('pcr1_primer1')
                add_vertical_space(1)

                st.subheader('Custom Volumes')
                ld.custom_volumes('cust_vol')
                selected_pids = {}
                with checkbox_cont:
                    checkbox_keys = dc.display_plate_checklist('pmr1_checklist',
                        ['dna','pcr','taqwater1','primer'])

                    selected_pids = dc.collect_plate_checklist(checkbox_keys)
                    if not selected_pids['pcr']:
                        m('No PCR plates selected/available', level='display', dest=('css',), color='red',size='p')

                with pcr_col:
                    ld.add_pcr_barcodes('pcr_bc_tab1')
                with taqwater_col:
                    ld.add_taqwater_barcodes('tw_bc_tab1', pcr_stage=pcr_stage)

                with display_cont:
                    dc.display_pcr1_components(selected_pids)

                    add_vertical_space(1)
                    if selected_pids['dna']:
                        primer_max_vol = util.CAP_VOLS[util.PLATE_TYPES['Echo384']]
                        primer_dead_vol = util.DEAD_VOLS[util.PLATE_TYPES['Echo384']]
                        if not selected_pids['primer']:
                            st.warning(f'Primers will appear highlighted red if no primer plate files have been loaded')
                        st.write(f'Available volumes equal the measured volume - dead volume ({primer_dead_vol/1000}ul). Max primer volume is {primer_max_vol/1000}ul.')
                        iv.infoview_primers('pcr_tab1', dna_pids=selected_pids['dna'],
                                primer_pids=selected_pids['primer'], save_buttons=True)

                add_vertical_space(1)

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_pcr1', default_view1='Primers',
                default_view2='Files')

    #-------------------------------- Primers ~ TAB 2: Generate PCR 1 picklists --------------------------------
    if primer_tab == 2:
        st.session_state['primer_tab'] = 2
        if unlocked:
            with main_body_container:
                caller_id = 'pcr1 generate'
                st.subheader('Generate Primer (PCR 1) Echo Picklists')
                st.info('Select the resources you wish to include for primer PCR and then click on the '+\
                        '**Generate Echo Picklists** button below, to create Echo picklist files')
                checkbox_keys = dc.display_plate_checklist('pmr1_checklist',
                        ['dna','pcr','taqwater1','primer'])
                selected_pids = dc.collect_plate_checklist(checkbox_keys)
                if not selected_pids['pcr'] and not selected_pids['amplicon']:
                    m('No PCR or amplicon plates selected/added yet', level='display', dest=('css',), color='red',size='p')
                hline()
                dc.display_pcr1_components(selected_pids)
                hline()

                if selected_pids['dna']:
                    if exp.check_ready_pcr1(selected_pids, caller_id=caller_id):
                        _,button_col,cb_col,_ = st.columns([5, 2, 2, 4])
                        echo_picklist_go = button_col.button('Generate Echo Picklists',
                                    key='echo_pcr1_go_button',
                                    type='primary')
                        cb_force = cb_col.checkbox('Force picklist generation', key='force_pcr1_picklist')
                        if echo_picklist_go:
                            success = run_generate(exp, exp.generate_echo_PCR1_picklists,
                                    selected_pids, force=cb_force, caller_id=caller_id)
                            if not success:
                                st.error('Picklist generation failed. Please see the log')
                            else:
                                exp.add_pcr_wells(exp, selected_pids['pcr'], selected_pids['dna'])

                for msg,lvl in mq[caller_id]:
                    m(msg, lvl)
                mq[caller_id] = set()

                if pcr1_picklists_exist(exp):
                    dc.get_echo1_download_btns()

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_pcr2', default_view1='Primers',
                default_view2='Consumables')
