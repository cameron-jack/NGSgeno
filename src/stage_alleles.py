import streamlit as st
import load_data as ld
import info_viewer as iv
from stutil import init_state, m, mq, add_vertical_space, hline, unlocked
import display_components as dc
from generate import run_generate, generate_targets, generate_primer_assayfams
from util import unguard_pbc
import os, sys, subprocess
from pathlib import Path


def clean_up_match_files(rundir):
    """
    Remove any existing match progress files to ensure a clean start
    """
    lock_path = Path(rundir+'/ngsgeno_lock')
    if os.path.exists(lock_path):
        try:
            os.remove(lock_path)
        except Exception as exc:
            print(f'Could not remove lock file: {exc}', file=sys.stderr)
            #m(f'Could not remove lock file: {exc}', level='error', dest='noGUI')
    progress_files = list(Path(rundir).glob('match_progress_*'))
    for pf in progress_files:
        try:
            os.remove(pf)
        except Exception as exc:
            #m(f'Could not remove progress file: {exc}', level='error', dest='noGUI')
            print(f'Could not remove progress file: {exc}', file=sys.stderr)


def stage_alleles(exp, main_body_container, upper_container, message_container):
    tab_col1, tab_col2, tab_col3 = upper_container.columns([5,5,1])
    with tab_col1:
        allele_tab = dc.create_tabs([("Allele Calling", "For Rodentity Mice"),("Amplicon Calling", "For custom amplicons")])
    if not allele_tab:
        init_state('allele_tab', 1)
        allele_tab = st.session_state['allele_tab']

    # Only offer upload in the Miseq pipeline section
    #-------------------------------- Allele ~ TAB 1: Allele calling --------------------------------
    if allele_tab == 1:
        st.session_state['allele_tab'] = 1
        with main_body_container:
            rundir = exp.get_exp_dn()
            caller_id = 'pre-execute-gt-analysis'
            success = exp.check_sequence_upload_ready(caller_id)
            for msg,lvl in mq[caller_id]:
                m(msg, level=lvl, no_log=True)
            mq[caller_id] = set()
            ld.display_fastqs('disp_fastq1')
            if not success:
                st.warning('Resources are required for allele calling:')
                st.subheader('Upload reference sequences')
                ld.load_rodentity_references('reference_allele1')
                #ld.load_miseq_fastqs('miseq_tab1a')
            else:
                #ld.load_miseq_fastqs('miseq_tab1b')
                # check whether output files already exist and are opened elsewhere
                targets_fn = exp.get_exp_fn('target.fa')
                primers_fn = exp.get_exp_fn('primers.csv')
                results_fn = exp.get_exp_fn('results.csv')
                matchlog_fn = exp.get_exp_fn('match.log')
                for fn in [targets_fn, primers_fn, results_fn, matchlog_fn]:
                    if Path(fn).exists():
                        try:
                            os.rename(fn, fn.replace('.','_tmp_swap.'))
                            os.rename(fn.replace('.','_tmp_swap.'), fn)
                        except PermissionError:
                            st.error(f'{fn} appears to be in use. Please close this file before continuing')

                with st.form('allele_calling_form', clear_on_submit=True):
                    all_fns = [fp for fp in exp.uploaded_files.keys() if fp not in {'_upload_queue','_upload_pending'}]
                    ref_fns = [fp for fp in all_fns if exp.uploaded_files[fp]['purpose'] \
                            in ['rodentity_reference','custom_reference']]
                    ref = st.selectbox('Select references to match against', options=ref_fns)
                    #num_unique_seq = st.number_input("Number of unique sequences per work unit", value=1)
                    cpus_avail = max(1, os.cpu_count()-2)
                    num_cpus = st.number_input(\
                            label=f"Number of processes to run simultaneously, default: {cpus_avail}",\
                                    value=cpus_avail, min_value=1)
                    margin = st.number_input(label="Require lengths of read sequences and target "+\
                            "references to be proportionally similar by this amount. Value must be between 0.0 and 1.0 "+\
                            "default 0.9", format='%f',min_value=0.0, step=0.05,value=0.9)
                    identity = st.number_input(label="Proportion of identity required for inexact match "+\
                            ", default 0.9. Must be between 0.0 and 1.0",
                            format='%f',min_value=0.0, max_value=1.0, value=0.9, step=0.05)
                    mincov = st.number_input(label="Do not match unique sequences with less than this "+\
                            "many reads coverage, default 5", format='%i',min_value=0, step=1,value=5)
                    minprop = st.number_input(label="Do not match unique sequences with less than this "+\
                            "proportion of the reads seen for the most observed (expected) allele, default 0.1. Must be between 0.0 and 1.0",
                            format='%f',min_value=0.0, max_value=1.0, value=0.1, step=0.05)
                    inexact_mode = st.checkbox("Enable inexact matching")
                    exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                    debug_mode = st.checkbox('Turn on debugging for allele calling')
                    do_matching = st.form_submit_button("Run allele calling")

                if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                    match_fn = exp.get_exp_fn('match_progress_100_100')
                    if Path(match_fn).exists():
                        clean_up_match_files()
                        st.rerun()

                elif do_matching:
                    success = generate_targets(exp, ref, caller_id=caller_id)
                    if not success:
                        m('failed to save reference sequences to target file', level='critical')
                        sleep(0.5)
                    success = generate_primer_assayfams(exp, caller_id=caller_id)
                    if not success:
                        m('failed to save primers and assay families to file', level='critical')
                        sleep(0.5)
                    else:
                        matching_prog = os.path.join('src','ngsmatch.py')
                        cmd_str = f'{sys.executable} {matching_prog} ' +\
                                    f'--ncpus {num_cpus} --rundir {rundir} '+\
                                    f'--margin {margin} --identity {identity} --mincov {mincov} '+\
                                    f'--minprop {minprop}'
                        if inexact_mode:
                            cmd_str += ' --inexact'
                        if exhaustive_mode:
                            cmd_str += ' --exhaustive'
                        if debug_mode:
                            cmd_str += ' --debug'
                        m(f'{cmd_str}', level='info')
                        st.write(f'Calling {cmd_str}')

                        launch_msg = st.empty()
                        launch_prog = st.progress(0)
                        completion_msg = st.empty()
                        match_prog = st.progress(0)
                        st.session_state['matching_in_progress'] = (rundir, launch_msg, launch_prog,
                                completion_msg, match_prog)
                        if sys.platform == "win32":
                            subprocess.Popen(cmd_str.split(' '),
                                    creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS,
                                    stdin=subprocess.DEVNULL,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    close_fds=True,
                            )
                        else:
                            subprocess.Popen(cmd_str.split(' '),
                                    start_new_session=True,
                                    stdin=subprocess.DEVNULL,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    close_fds=True,
                            )
                        # cp = subprocess.Popen(cmd_str.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        # st.session_state['matching_process'] = cp

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_allele1', default_view1='Files',
                default_view2='Plates')

    #-------------------------------- Allele ~ TAB 2: Amplicon calling --------------------------------
    elif allele_tab == 2:
        st.session_state['allele_tab'] = 2
        with main_body_container:
            rundir = exp.get_exp_dn()
            caller_id = 'pre-execute-amp-analysis'
            st.subheader('Amplicon Calling')
            st.info('This section is for custom amplicon calling. Amplicon sequences are always inexact matched against their expected primer to find variants. Please choose the amplicon '+\
                    'reference you wish to use or upload a new one')
            ld.load_amplicon_references('amplicon_ref1')

            targets_fn = exp.get_exp_fn('amplicon_targets.fa')
            results_fn = exp.get_exp_fn('amplicon_results.csv')
            matchlog_fn = exp.get_exp_fn('amplicon_match.log')
            # check whether output files already exist and are opened elsewhere
            for fn in [targets_fn, results_fn, matchlog_fn]:
                if Path(fn).exists():
                    try:
                        os.rename(fn, fn.replace('.','_tmp_swap.'))
                        os.rename(fn.replace('.','_tmp_swap.'), fn)
                    except PermissionError:
                        st.error(f'{fn} appears to be in use. Please close this file before continuing')
            # choose amplicon plate(s) and reference file(s) to run
            plate_col, ref_col = st.columns(2)
            with plate_col:
                st.subheader('Select amplicon plates for allele calling')
                checkbox_keys = dc.display_plate_checklist('amplicon_checklist', ['amplicon'], default_value=False)
                selected_pids = dc.collect_plate_checklist(checkbox_keys)
                if not selected_pids['amplicon']:
                    m('No amplicon plates selected', level='display', dest=('css',), color='red',size='p')
                else:
                    st.write(f'Selected amplicon plates: {", ".join([unguard_pbc(gpid, silent=True) for gpid in selected_pids["amplicon"]])}')
            with ref_col:
                st.subheader('Select amplicon references to match against')
                checkbox_keys = dc.display_file_checklist('amplicon_reference_checklist', ['amplicon_reference'], default_value=False)
                selected_refs = dc.collect_file_checklist(checkbox_keys)
                if not selected_refs['amplicon_reference']:
                    m('No amplicon references selected', level='display', dest=('css',), color='red',size='p')
                else:
                    st.write(f'Selected amplicon references: {", ".join(dc.fns_from_checklist(selected_refs))}')

            with st.form('amplicon_calling_form', clear_on_submit=True):
                cpus_avail = max(1, os.cpu_count()-2)
                num_cpus = st.number_input(\
                        label=f"Number of processes to run simultaneously, default: {cpus_avail}",\
                                value=cpus_avail, min_value=1)
                margin = st.number_input(label="Require lengths of read sequences and target "+\
                        "references to be proportionally similar by this amount. Value must be between 0.0 and 1.0 "+\
                        "default 0.9", format='%f',min_value=0.0, step=0.05,value=0.9)
                identity = st.number_input(label="Proportion of identity required for inexact match "+\
                        ", default 0.9. Must be between 0.0 and 1.0",
                        format='%f',min_value=0.0, max_value=1.0, value=0.9, step=0.05)
                mincov = st.number_input(label="Do not match unique sequences with less than this "+\
                        "many reads coverage, default 5", format='%i',min_value=0, step=1,value=5)
                minprop = st.number_input(label="Do not match unique sequences with less than this "+\
                        "proportion of the reads seen for the most observed (expected) allele, default 0.1. Must be between 0.0 and 1.0",
                        format='%f',min_value=0.0, max_value=1.0, value=0.1, step=0.05)
                exhaustive_mode = st.checkbox("Exhaustive mode: try to match every sequence, no matter how few counts")
                debug_mode = st.checkbox('Turn on debugging for allele calling')
                do_matching = st.form_submit_button("Run amplicon calling")

            if Path(exp.get_exp_fn('ngsgeno_lock')).exists():
                st.info('Analysis in progress')
            elif do_matching:
                success = generate_targets(exp, selected_refs['amplicon_reference'])
                if not success:
                    m('failed to save reference sequences to target file', level='critical')
                    sleep(0.5)
                else:
                    matching_prog = Path('src/ngsmatch.py')
                    cmd_str = f'{sys.executable} {matching_prog} --ncpus {num_cpus} --rundir {rundir} '+\
                            f'--mincov {mincov} --minprop {minprop}'
                    if exhaustive_mode:
                        cmd_str += ' --exhaustive'
                    if debug_mode:
                        cmd_str += ' --debug'
                    if selected_pids['amplicon']:
                        for pid in selected_pids['amplicon']:
                            cmd_str += f' --amplicons {",".join(selected_pids["amplicon"])}'
                        cmd_str += f' --targets amplicon_targets.fa'
                        m(f'{cmd_str}', level='info')
                        st.write(f'Calling {cmd_str}')
                        launch_msg = st.empty()
                        launch_prog = st.progress(0)
                        completion_msg = st.empty()
                        match_prog = st.progress(0)
                        st.session_state['matching_in_progress'] = (rundir, launch_msg, launch_prog,
                                completion_msg, match_prog)
                        if sys.platform == "win32":
                            subprocess.Popen(cmd_str.split(' '),
                                    creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS,
                                    stdin=subprocess.DEVNULL,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    close_fds=True,
                            )
                        else:
                            subprocess.Popen(cmd_str.split(' '),
                                    start_new_session=True,
                                    stdin=subprocess.DEVNULL,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    close_fds=True,
                            )
                        #cp = subprocess.Popen(cmd_str.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        #st.session_state['matching_process'] = cp

        # ** Info viewer **
        iv.upper_info_viewer_panel(tab_col3, tab_col2, 'upper_allele2', default_view1='Files',
                default_view2='Plates')
