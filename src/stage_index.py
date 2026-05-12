import streamlit as st
import load_data as ld
import info_viewer as iv
from stutil import init_state, m, add_vertical_space, hline, unlocked, set_state, mq
import display_components as dc
from generate import run_generate, pcr2_picklists_exist
import util

def stage_index(exp, main_body_container, upper_container, message_container):
    pcr_stage = 2
    tab_col1, tab_col2,tab_col3 = upper_container.columns([5,5,1])

    #Tab setup
    with tab_col1:
        index_tab = dc.create_tabs([("PCR 2", "Components"), ("Generate", "Picklists")])
    if not index_tab:
        init_state('index_tab', 1)
        index_tab = st.session_state['index_tab']

    #------------------------------------- Index ~ TAB 1: PCR 2 Components -------------------------------------
    if index_tab == 1:
        st.session_state['index_tab'] = 1
        if unlocked(exp):
            with main_body_container:
                st.subheader('Indexing (PCR 2) Components')
                st.info('Provide the resources needed to perform a sufficient number of indexing reactions '+\
                        'for your experiment, then move to the *Generate Picklists* tab')

                checkbox_cont = st.container()
                hline()
                add_vertical_space(1)

                display_cont = st.container()
                hline()
                add_vertical_space(1)

                st.info('Indexing (PCR 2) requires index index plates, taq/water plates, '+\
                        'and either Echo plates (prepared by the Nimbus) or amplicon plates')

                st.subheader('Add Barcodes', help='Add barcodes for plates')
                pcr_col, taqwater_col = st.columns(2)

                st.subheader('Upload Files')
                ld.load_pcr2_files('pcr2_index1')
                add_vertical_space(1)

                st.subheader('Custom Volumes')
                ld.custom_volumes(exp)

                with checkbox_cont:
                    checkbox_keys = dc.display_plate_checklist('idx_checklist1',
                            ['pcr','taqwater2','amplicon','index'])

                    selected_pids = dc.collect_plate_checklist(checkbox_keys)
                    if not selected_pids['pcr'] and not selected_pids['amplicon']:
                        m('No PCR or amplicon plates selected', level='display',
                                dest=('css',), color='red',size='p')

                with display_cont:
                    dc.display_pcr2_components(selected_pids)

                with pcr_col:
                    ld.add_pcr_barcodes('pcr_bc_tab2')
                with taqwater_col:
                    ld.add_taqwater_barcodes('tw_bc_tab2', pcr_stage=pcr_stage)
                add_vertical_space(1)

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_index1', default_view1='Indexes',
                default_view2='Consumables')

    #--------------------------------- Index ~ TAB 2: Generate PCR 2 picklists ---------------------------------
    if index_tab == 2:
        st.session_state['index_tab'] = 2
        init_state('run_index', False)
        if unlocked(exp):
            with main_body_container:
                caller_id = 'pcr2 generate'
                st.subheader('Generate Index (PCR 2) Echo Picklists')
                st.info('Select the resources you wish to include for indexing PCR and then click on the '+\
                        '**Generate Echo Picklists** button below, to create Echo picklist files')
                do_generate = False
                checkbox_keys = dc.display_plate_checklist('idx_checklist2',
                        ['pcr','taqwater2','amplicon','index'])
                selected_pids = dc.collect_plate_checklist(checkbox_keys)
                if not selected_pids['pcr'] and not selected_pids['amplicon']:
                    m('No PCR or amplicon plates selected',
                            level='display', dest=('css',), color='red',size='p')
                hline()
                print(f'{selected_pids=}', flush=True)
                dc.display_pcr2_components(selected_pids)
                print(f'{selected_pids=}', flush=True)
                hline()
                add_vertical_space(1)

                if selected_pids['pcr']:
                    success = exp.check_ready_pcr2(selected_pids, caller_id=caller_id)
                    if success:
                        do_generate = True
                elif selected_pids['amplicon']:
                    success = exp.check_ready_pcr2(selected_pids, caller_id=caller_id, amplicon_only=True)
                    if success:
                        do_generate = True
                else:
                    m('Either PCR plates or Amplicon plates must exist for indexing to begin',
                            level='display')

                for msg,lvl in mq[caller_id]:
                    m(msg, lvl)
                mq[caller_id] = set()

                if do_generate:
                    _,picklist_button_col,_ = st.columns([2, 2, 1])
                    echo_picklist_go = picklist_button_col.button('Generate Echo Picklists',\
                            key='echo_pcr2_go_button')
                    picklist_button_col.write('')
                    if echo_picklist_go:
                        if selected_pids['pcr']:
                            st.session_state['run_index'] = True
                        elif selected_pids['amplicon'] and not selected_pids['pcr']:
                            _, amp1, amp2, amp3, _ = st.columns([3, 3, 1, 1, 3])
                            amp1.warning('Create picklists with only amplicons?')
                            amp2.button('Yes', on_click=set_state, args=('run_index', True))
                            amp3.button('No', on_click=set_state, args=('run_index', False))

                if st.session_state['run_index']:
                    success = run_generate(exp, exp.generate_echo_PCR2_picklists,
                            selected_pids)
                    set_state('run_index', False)
                    if not success:
                        m('Picklist generation failed. Please see the log', level='display')

                if pcr2_picklists_exist(exp):
                    dc.show_echo2_outputs()

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_index2', default_view1='Files',
                default_view2='Files')
