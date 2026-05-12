import streamlit as st
import load_data as ld
import info_viewer as iv
from stutil import init_state, m, add_vertical_space, hline, unlocked, mq
import display_components as dc
import display_alignments as da
from pathlib import Path


def stage_reports(exp, main_body_container, upper_container, message_container):
    results_fp = exp.get_exp_fn('results.csv')
    tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    with tab_col1:
        allele_tab = dc.create_tabs([("Sequencing results", ""),("Amplicon results", "")])
    if not allele_tab:
        init_state('results_tab', 1)
        allele_tab = st.session_state['results_tab']

    with main_body_container:
        if allele_tab == 1:
            if not Path(results_fp).exists():
                m('**No allele calling results available**', level='display', dest=('mkdn'))
            else:
                if not Path(exp.get_exp_fn('amplicon_targets.fa')).exists():
                    m('**No amplicon targets available**', level='display', dest=('mkdn'))
                else:
                    # Three displays for Rodentity, custom, and other results
                    key = 'rodentity_results_display'
                    da.show_results_display('rodentity', key, caller_id=key)
                    key = 'custom_results_display'
                    da.show_results_display('custom', key, caller_id=key)
                    key = 'other_results_display'
                    da.show_results_display('other', key, caller_id=key)

        elif allele_tab == 2:
            if not Path(exp.get_exp_fn('amplicon_results.csv')).exists():
                m('**No amplicon calling results available**', level='display', dest=('mkdn'))
            else:
                if not Path(exp.get_exp_fn('amplicon_targets.fa')).exists():
                    m('**No amplicon targets available**', level='display', dest=('mkdn'))
                else:
                    key = 'amplicon_results_display'
                    da.show_results_display('amplicon', key, caller_id=key)
                    
    # ** Info viewer **
    iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_report1', default_view1='Files',
            default_view2='Plates')
