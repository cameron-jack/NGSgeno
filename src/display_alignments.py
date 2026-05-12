import streamlit as st
import jsonpickle
import warnings
import pandas as pd
import sys
import re
from pathlib import Path
from time import sleep
from stutil import m, mq
import display_components as dc

from generate import generate_pdf

from cogent3 import get_app, make_unaligned_seqs, make_tree
from cogent3.align.progressive import tree_align


def show_results_table(hdr, well_data, rtype, filter_dict, filter_param_fn, key):
    """
    build an interactive table for amplicon results
    Args:
        hdr (list): header for the amplicon data
        well_data (list of lists): table of data, one row per well
        rtype (str): type of results, in ['rodentity','custom','other','amplicon']
        filter_dict (dict): dictionary of existing filter parameters
        filter_param_fn (str): file path to the filter parameters file
        key (str): key (unique id) for the view component
    Returns:
        list of chosen variants
    """
    # hide these columns but keep them available for the user
    hidden_cols = ['sex','sampleNo','sampleName','strain','clientName','alleleSymbol',
                'alleleKey','assayKey','assays','assayFamilies','primerPlate','primerWell',
                'pcrPlate','pcrWell','index_plate','i7bc','i7well','i7name','i5bc','i5well',
                'i5name','cleanCount','dnaPlate','dnaWell']
    if not well_data:
        st.info('All wells reported')
        return variant_info

    well_dataframe = pd.DataFrame(well_data, columns=hdr)
    well_table = dc.aggrid_interactive_table(well_dataframe, key=key+'_key',
            hidden=hidden_cols, editable=['filtProportion'])
    if well_table and 'selected_rows' in well_table and well_table['selected_rows'] is not None:
        selected_rows = well_table['selected_rows'].to_dict(orient='records')
        st.session_state[f'{rtype}_chosen_vars'] = get_variants_from_table(selected_rows, filter_dict, filter_param_fn)


def compare_var_to_ref(ref_seq, var_anno, display_width=120, colour_all=False, colour_changes=True):
    """
    Build the views for the reference sequences, sequence, variable and summary
    Args:
        ref_seq (str): dictionary of reference sequences
        var_anno (str): the sequence annotation chosen by the user
        display_width (int): width of the display in characters, wrap outputs to this width
        colour_all (bool): False, whether to colour every position in the display
        colour_changes (bool): True, whether to colour the changes of variable sites only
    Returns:
        outputs (list of tuples): each tuple contains four strings:
            ref_view, link_view, var_view, summary_view (str): the views for the reference sequences,
            the link between the reference and variable sequences, 
            the variable sequence, and a summary of the variable sequence
    Notes:
        This function builds on the variant sequence generated in ngsmatch.get_variant_seq()
        It also adjusts the reference sequence so that bases remain in matching positions
    """
    if ref_seq is None or var_anno is None:
        return [('', '', '', '')]

    parts = re.split(r'(\d+)', var_anno)
    rev_parts = parts[::-1]  # work from the end to the start
    var_seq = ref_seq  # variant sequence to be built up
    mod_seq = ref_seq  # modified sequence to be built up
    link_view = []  # link between reference and variable sequences
    summary_view = []  # summary of the variable sequence
    # need to go through changes in reverse order to avoid length changes from affecting position
    for i, p in enumerate(rev_parts):
        if i % 2 == 0:
            if p == '':
                break  # we've reached the end
            change = p
            if len(parts) <= i+1:
                print(f'Error: no matching change for position: {pos} in {parts} from {var_name}', flush=True)
                return ''
            try:
                pos = int(rev_parts[i+1]) -1  # 1-based
            except Exception as e:
                print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
                return ''
            #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
            if '+' in change:
                var_seq = var_seq[:pos] + change[1:] + var_seq[pos:]
                mod_seq = mod_seq[:pos] + '-'*len(change[1:]) + mod_seq[pos:]
            if '-' in change:
                var_seq = var_seq[:pos] + '-'*len(change[1:]) + var_seq[pos+len(change)-1:]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                var_seq = var_seq[:pos] + repl_bases +var_seq[pos+len(repl_bases):]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
    link_view = ['|' if c1 == c2 else ' ' for c1, c2 in zip(mod_seq, var_seq)]
    summary_view = []
    for m,v in zip(mod_seq, var_seq):
        if m == v:
            summary_view.append('*')
        elif m == '-':
            summary_view.append('+')
        elif v == '-':
            summary_view.append('-')
        elif v != m:
            summary_view.append(v)
        else:
            summary_view.append(' ')
    summary_view = ''.join(summary_view)
    outputs = []
    link_chrs = ''.join(link_view)
    summary_chrs = ''.join(summary_view)
    for i in range(0, len(var_seq), display_width):
        outputs.append((mod_seq[i:i+display_width],
                link_chrs[i:i+display_width],
                var_seq[i:i+display_width],
                summary_chrs[i:i+display_width]))
    return outputs


def reconstruct_sequence(var_anno, ref_seq):
    """
    Recreate the full variant sequence based on the variant annotations and the reference sequence
    args:
        var_anno (str): variant annotations
        ref_seq (str): the reference sequence variants are based on
    returns:
        variant_sequence (str): the aligned variant sequence
        mod_ref_seq (str): the matching aligned
    """
    if ref_seq is None or var_anno is None:
        return ''

    parts = re.split(r'(\d+)', var_anno)
    rev_parts = parts[::-1]  # work from the end to the start
    var_seq = ref_seq  # variant sequence to be built up
    mod_seq = ref_seq  # modified sequence to be built up
    # need to go through changes in reverse order to avoid length changes from affecting position
    for i, p in enumerate(rev_parts):
        #print(f'reconstruct_sequence() {i=} {p=}', flush=True)
        if i % 2 == 0:
            if p == '':
                break  # we've reached the end
            change = p
            if len(parts) <= i+1:
                 print(f'Error: no matching change for position: {pos} in {parts}', flush=True)
                 return ''
            try:
                pos = int(rev_parts[i+1]) -1  # 1-based
            except Exception as e:
                print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
                return ''
            #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
            if '+' in change:
                var_seq = var_seq[:pos] + change[1:] + var_seq[pos:]
                mod_seq = mod_seq[:pos] + '-'*len(change[1:]) + mod_seq[pos:]
            if '-' in change:
                var_seq = var_seq[:pos] + '-'*len(change[1:]) + var_seq[pos+len(change)-1:]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                var_seq = var_seq[:pos] + repl_bases +var_seq[pos+len(repl_bases):]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
    return var_seq


def build_aligned_pair(var_anno, ref_seq):
    """
    Recreate the full variant sequence based on the variant annotations and the reference sequence
    args:
        var_anno (str): variant annotations
        ref_seq (str): the reference sequence variants are based on
    returns:
        var_seq (str): the aligned variant sequence
        mod_seq (str): the matching aligned reference sequence
    """
    if ref_seq is None or var_anno is None:
        return '',''

    parts = re.split(r'(\d+)', var_anno)
    rev_parts = parts[::-1]  # work from the end to the start
    var_seq = ref_seq  # variant sequence to be built up
    mod_seq = ref_seq  # modified sequence to be built up
    # need to go through changes in reverse order to avoid length changes from affecting position
    for i, p in enumerate(rev_parts):
        if i % 2 == 0:
            if p == '':
                break  # we've reached the end
            change = p
            if len(parts) <= i+1:
                print(f'Error: no matching change for position: {pos} in {parts} from {var_name}', flush=True)
                return ''
            try:
                pos = int(rev_parts[i+1]) -1  # 1-based
            except Exception as e:
                print(f'Error: could not convert position {rev_parts[i+1]=} to integer', flush=True)
                return ''
            #print(f'{rev_parts=} {i=} {p=} {pos=} {new_seq=}',flush=True)
            if '+' in change:
                var_seq = var_seq[:pos] + change[1:] + var_seq[pos:]
                mod_seq = mod_seq[:pos] + '-'*len(change[1:]) + mod_seq[pos:]
            if '-' in change:
                var_seq = var_seq[:pos] + '-'*len(change[1:]) + var_seq[pos+len(change)-1:]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
            if '/' in change:
                repl_bases = change.split('/')[1] 
                var_seq = var_seq[:pos] + repl_bases +var_seq[pos+len(repl_bases):]
                mod_seq = mod_seq[:pos] + mod_seq[pos:]
    return var_seq, mod_seq


def run_msa(ref_id_seq, var_list):
    """
    Reconstruct full variant sequences and then perform multiple sequence alignment
    args:
        ref_id_seq (tuple(str,str)): tuple of reference id and reference seq
        var_list (list[str]): a list of variant annotations
    returns:
        aligned (cogent3 alignment object)
    """
    # need to reconstruct the full sequences from the annotations
    reconstructed_seqs = {ref_id_seq[0]:ref_id_seq[1].replace('-','')}  # start with the reference sequence
    for var_anno in var_list:
        reconstructed_seqs[var_anno] = reconstruct_sequence(var_anno, ref_id_seq[1]).replace('-','')  # remove gaps for alignment
    no_tree = make_tree(tip_names=list(reconstructed_seqs.keys()), underscore_unmunge=True)
    unaligned_seqs = make_unaligned_seqs(reconstructed_seqs, moltype='dna')
    with warnings.catch_warnings(action="ignore"):
        aln, tree = tree_align("HKY85", unaligned_seqs, tree=no_tree, show_progress=False)
    return aln


def get_variants_from_table(well_table, filter_dict, filter_param_fn):
    """
    Get the chosen variant from the interactive table
    args:
        well_table is a list of dicts, one per selected row
        filter_dict (dict): dictionary of existing filter parameters
        filter_param_fn (str): file path to the filter parameters file
    """
    exp = st.session_state['experiment']
    variant_info = []
    first_variant = well_table[0]
    for field in ['samplePlate','sampleWell','sampleBarcode','primer','mergeCount','seqCount','otherCount','otherName','filtProportion']:
        if field in first_variant:
            variant_info.append(first_variant[field])
        else:
            variant_info.append('')
    samplePlate = variant_info[0]
    sampleWell = variant_info[1]
    sampleBarcode = variant_info[2]
    filtProportion = variant_info[-1]
    new_id = f'{samplePlate}\t{sampleWell}\t{sampleBarcode}'
    if new_id not in filter_dict or filter_dict[new_id] != filtProportion:
        filter_dict[new_id] = filtProportion
        with open(filter_param_fn, 'wt') as f:
            for id in filter_dict:
                sp, sw, sb = id.split('\t')
                f.write(f'{sp}\t{sw}\t{sb}\t{filter_dict[id]}\n')
            if new_id not in filter_dict:
                f.write(f'{samplePlate}\t{sampleWell}\t{sampleBarcode}\t{filtProportion}\n')
    return variant_info


def make_summary_line_and_format(var_seqs, wt_seq):
    """
    Make a summary line for the chosen variant as well as a collection of formatting instructions
    """
    #print(f'make_summary_line_and_format {var_seqs=} {wt_seq=}', file=sys.stderr)
    summary_line = []
    formatting = []
    seq_array = list(var_seqs.values())
    for i in range(len(seq_array[0])):
        char_set = set([seq_array[j][i] for j in range(len(seq_array)) if seq_array[j][i] == '-' or seq_array[j][i] != wt_seq[i]])
        indel = False
        if '-' in char_set:
            indel = True
            char_set.remove('-')
        if len(char_set) == 2:
            if 'A' and 'G' in char_set:
                summary_line.append('R')
            elif 'C' and 'T' in char_set:
                summary_line.append('Y')
            elif 'G' and 'C' in char_set:
                summary_line.append('S')
            elif 'A' and 'T' in char_set:
                summary_line.append('W')
            elif 'G' and 'T' in char_set:
                summary_line.append('K')
            elif 'A' and 'C' in char_set:
                summary_line.append('M')
        elif len(char_set) == 3:
            if 'C' in char_set and 'G' in char_set and 'T' in char_set:
                summary_line.append('B')
            elif 'A' in char_set and 'G' in char_set and 'T' in char_set:
                summary_line.append('D')
            elif 'A' in char_set and 'C' in char_set and 'T' in char_set:
                summary_line.append('H')
            elif 'A' in char_set and 'C' in char_set and 'G' in char_set:
                summary_line.append('V')
        elif len(char_set) == 4:
            summary_line.append('N')
        elif len(char_set) == 1:
            summary_line.append(char_set.pop())
        elif indel:
            summary_line.append('-')
        else:
            summary_line.append('*')

    return ''.join(summary_line), formatting


def show_alignments(chosen_vars, tmp_fn, pdf_fn, target_fn, key):
    """
    Show the alignments for the chosen variants
    Run alignments against the reference sequence for the chosen primer, found in amplicon_targets.fa
    Show whether the amplicon has been reported or not using amplicon_reported.csv
    """
    exp = st.session_state['experiment']
    amplicons_reported = set()  # plate\twell\tbarcode\tprimer
    if not chosen_vars:
        st.empty()
        return None
    primer_chosen = chosen_vars[3]
    try:
        merged_counts = int(chosen_vars[4])
    except ValueError:
        merged_counts = 0
    try:
        wt_counts = int(chosen_vars[5])
    except ValueError:
        wt_counts = 0
    try:
        min_prop = float(chosen_vars[-1])
    except ValueError:
        min_prop = 0.15
    ref_seqs = {}
    with open(exp.get_exp_fn(target_fn), 'rt') as rfn:
        for line in rfn:
            if line.startswith('>'):
                ref_name = line[1:].strip()
                ref_seqs[ref_name] = ''
            else:
                ref_seqs[ref_name] += line.strip()
    ref_chosen = None
    if primer_chosen in ref_seqs:
        ref_chosen = primer_chosen
    else:
        for r in ref_seqs:
            if primer_chosen in r and 'wt' in r.lower():
                ref_chosen = r
                break
    if not ref_chosen:
        return None
    
    #print(f'{chosen_vars=} {ref_chosen=}', file=sys.stderr)
    var_list = []
    var_list_counts = {}
    for cv,cc in zip(chosen_vars[-2].split(';'),chosen_vars[-3].split(';')):
        if '//' in cv and ref_chosen in cv:
            if int(cc) >= (merged_counts * min_prop):
                var_list.append(cv.split('//')[1])
                var_list_counts[cv.split('//')[1]] = int(cc)

    #print(var_list, file=sys.stderr)
    if not var_list:
        #st.warning('No variant sequences found for this amplicon')
        return None

    #vars_chosen = {cv.split('//')[1]:int(cc) for cv,cc in zip(chosen_vars[-2].split(';'),chosen_vars[-3].split(';')) if ref_chosen in cv and '//' in cv and int(cc)>=(merged_counts*min_prop)}
    id = f'{chosen_vars[0]}\t{chosen_vars[1]}\t{chosen_vars[2]}\t{chosen_vars[3]}'
    #print(f'show_alignments() {ref_chosen=} {ref_seqs[ref_chosen]=}, {var_list=}', file=sys.stderr)
    ref_seq_id = (ref_chosen, ref_seqs[ref_chosen])
    aligned = run_msa((ref_chosen,ref_seqs[ref_chosen]),var_list)
    var_df = pd.DataFrame({
        'Report': [True for name in aligned.names],
        'Variant': [name for name in aligned.names],
        'Counts': [var_list_counts[name] if name != ref_chosen else wt_counts for name in aligned.names],
        'Sequence': [str(aligned.get_gapped_seq(name)) for name in aligned.names]
    })
    # st.markdown(
    #     """
    #     <style>
    #     /* Apply monospace to all text within the main content area */
    #     [class^=dvn-scroller] {
    #         font-family: Courier New", Courier, monospace;
    #     }
    #     </style>
    #     """,
    #     unsafe_allow_html=True
    # )

    # Apply the style to the DataFrame
    var_de = st.data_editor(var_df, key='alignment_editor_'+key, hide_index=True,
        column_config={
            'Report': st.column_config.CheckboxColumn('Report', width=40,
                    help='Check to include this variant in the PDF report', default=True),
            'Variant': st.column_config.TextColumn('Variant', width=150, disabled=True,
                    help='Name of the variant sequence'),
            'Counts': st.column_config.NumberColumn('Counts', width=80, disabled=True,
                    help='Number of reads supporting this variant'),
            'Sequence': st.column_config.TextColumn('Aligned sequence', width=1000, disabled=True,
                    help='The aligned sequence of this variant')
        },
        disabled=['Variant','Counts','Sequence'])

    prefix = st.text_input(label='Enter any descriptive notes you wish to prefix to this alignment')
    make_report_entry = st.button('Make report entry', key='build_pdfs_'+chosen_vars[0]+'_'+chosen_vars[1]+'_'+chosen_vars[2])
    if make_report_entry:
        if var_de is not None:
            aligned_stuff = [(n,f'{s}') for n,s,c in zip(var_de['Variant'], var_de['Sequence'], var_de['Report']) if c]
        #valid_keys = [cb.split('_seq_checkbox')[0] for cb in checkbox_keys if cb in st.session_state and st.session_state[cb]]
        #if valid_keys:
            #if not prefix:
            #    prefix = ''
            #aligned_stuff = {n:f'{s}' for n,s in aligned_names_seqs.items() if n in valid_keys}
            #aligned_stuff = [(n,f'{s}') for n,s in aligned_names_seqs.items() if n in valid_keys]
            if aligned_stuff:
                formatting = {'font':'Courier New', 'fontsize':10, 'lineheight':1.2, 'colour_changes':True}
                formatting = [(k,v) for k,v in formatting.items()]
                page = (prefix, tuple(aligned_stuff), tuple(formatting))
                if Path(tmp_fn).exists():
                    with open(tmp_fn, 'rt') as f:
                        pages = jsonpickle.decode(f.read(), keys=True, handle_readonly=True)
                        #print(type(pages), pages)
                        pages = set(pages)
                        pages.add(page)
                else:
                    pages = set()
                    pages.add(page)

                json_pages = jsonpickle.encode(list(pages), indent=4, keys=True, warn=True, handle_readonly=True)
                with open(tmp_fn, 'wt') as fout:
                    print(json_pages, file=fout)

                generate_pdf(pdf_fn, pages)
                st.success(f'Report entry made in {pdf_fn}')
                st.success(f'Adding {id} to reported amplicons list')
                return id
    
    return None


def build_paired_sequence_views(ref_seq, var_anno, display_width=120, colour_all=False, colour_changes=True):
    """
    Build the views for the reference sequence, variant sequence, link and summary
    Args:
        ref_seq (str): reference sequence
        var_anno (str): the sequence annotation chosen by the user
        display_width (int): width of the display in characters, wrap outputs to this width
        colour_all (bool): False, whether to colour every position in the display
        colour_changes (bool): True, whether to colour the changes of variable sites only
    Returns:
        outputs (list of tuples): each tuple contains four strings:
            ref_view, link_view, var_view, summary_view (str): the views for the reference sequences,
            the link between the reference and variable sequences,
            the variable sequence, and a summary of the variable sequence
    Notes:
        This function builds on the variant sequence generated in ngsmatch.get_variant_seq()
        It also adjusts the reference sequence so that bases remain in matching positions
    """
    if ref_seq is None or var_anno is None:
        return [('', '', '', '')]

    var_seq, mod_seq = build_aligned_pair(var_anno, ref_seq)

    # link view shows direct matches between sequences
    link_view = ['|' if c1 == c2 else ' ' for c1, c2 in zip(mod_seq, var_seq)]
    summary_view = []  # summary of differences
    for m,v in zip(mod_seq, var_seq):
        if m == v:
            summary_view.append('*')
        elif m == '-':
            summary_view.append('+')
        elif v == '-':
            summary_view.append('-')
        elif v != m:
            summary_view.append(v)
        else:
            summary_view.append(' ')
    summary_view = ''.join(summary_view)
    outputs = []
    link_chrs = ''.join(link_view)
    summary_chrs = ''.join(summary_view)
    for i in range(0, len(var_seq), display_width):
        outputs.append((mod_seq[i:i+display_width],
                link_chrs[i:i+display_width],
                var_seq[i:i+display_width],
                summary_chrs[i:i+display_width]))
    return outputs


def show_results_display(results_type, key, caller_id=None):
    """
    Show the results display for the chosen variants
    Args:
        hdr (list): header for the amplicon data
        amplicon_data (list of lists): the amplicon data
        filter_dict (dict): dictionary of per-amplicon filter parameters
        key (str): key (unique id) for the view component
    """
    exp = st.session_state['experiment']
    if results_type not in ['rodentity','custom','other','amplicon']:
        m(f'Unknown results type {results_type} in show_results_display()', caller_id=caller_id, level='critical')
        return
    rtype = results_type
    perform_reset = False
    pdf_fn = exp.get_exp_fn(f'{rtype}_alignments.pdf')
    tmp_fn = exp.get_exp_fn(f'{rtype}_alignments.json')
    reported_fn = exp.get_exp_fn(f'{rtype}_reported.txt')
    filter_param_fn = exp.get_exp_fn(f'{rtype}_table_filter_params.txt')
    if results_type == 'amplicon':
        results_fn = exp.get_exp_fn('amplicon_results.csv')
        target_fn = 'amplicon_targets.fa'
    else:
        results_fn = exp.get_exp_fn('results.csv')
        target_fn = 'targets.fa'
    # if rtype == 'amplicon':
    #     hdr, data, filter_dict, data_reported, data_unreported = exp.gather_amplicon_results(reported_amps_fn)
    # else:
    hdr, data, filter_dict, data_reported, data_unreported = exp.gather_results(rtype, reported_fn, results_fn)
    unreported_container = st.container(key=f'{rtype}_unreported_container')
    alignment_container = st.container(key=f'{rtype}_alignment_container')
    reported_container = st.container(key=f'{rtype}_reported_container')

    with unreported_container:
        if data_unreported:
            st.write(f'Unreported {rtype}:')
            show_results_table(hdr, data_unreported, rtype, filter_dict, filter_param_fn, f'unreported_{rtype}')
    

    with alignment_container:
        if f'{rtype}_chosen_vars' in st.session_state and st.session_state[f'{rtype}_chosen_vars']:
            #well_variants = get_variants_from_table(st.session_state[f'{rtype}_chosen_vars'], filter_dict, filter_param_fn)
            #print('Ready to do alignments', st.session_state[f'{rtype}_chosen_vars'], file=sys.stderr)
            with st.spinner('Generating alignments...'):
                st.session_state[f'{rtype}_chosen_id'] = show_alignments(st.session_state[f'{rtype}_chosen_vars'], tmp_fn, pdf_fn, target_fn, f'{rtype}_alignment_viewer')
        else:
            with st.spinner('Generating alignments...'):
                st.session_state[f'{rtype}_chosen_id'] = show_alignments(None, tmp_fn, pdf_fn, target_fn, f'{rtype}_alignment_viewer')
        if f'{rtype}_chosen_id' in st.session_state and st.session_state[f'{rtype}_chosen_id']:
            entries_reported = set()
            if Path(reported_fn).exists():
                with open(reported_fn, 'rt') as f:
                    for line in f:
                        entries_reported.add(line.strip())
            entries_reported.add(st.session_state[f'{rtype}_chosen_id'])
            with open(reported_fn, 'wt') as fout:
                for er in entries_reported:
                    print(er, file=fout)
            perform_reset = True
        else:
            st.empty()

    hdr, data, filter_dict, data_reported, data_unreported = exp.gather_results(rtype, reported_fn, results_fn)
    print(f'{data_reported=}', file=sys.stderr)
    with reported_container:
        if data_reported:
            st.write(f'Reported {rtype}:')
            show_results_table(hdr, data_reported, rtype, filter_dict, filter_param_fn, f'reported_{rtype}')
            wipe_report = st.button(f'Reset {rtype} report', key=f'reset_{rtype}_report_button')
            if wipe_report:
                success = True
                if Path(reported_fn).exists():
                    try:
                        Path(reported_fn).unlink()
                    except Exception as exc:
                        st.failure(f'Could not delete {reported_fn}, perhaps you have it open? {exc}')
                        success = False
                if Path(tmp_fn).exists():
                    try:
                        Path(tmp_fn).unlink()
                    except Exception as exc:
                        st.failure(f'Could not delete {tmp_fn}, perhaps you have it open {exc}')
                        success = False
                if Path(pdf_fn).exists():
                    try:
                        Path(pdf_fn).unlink()
                    except Exception as exc:
                        st.failure(f'Could not delete {pdf_fn}, perhaps you have it open? {exc}')
                        success = False
                if success:
                    st.success(f'{rtype} report reset to empty')
                    perform_reset = True


    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()

    if perform_reset:
        st.session_state[f'{rtype}_chosen_vars'] = None
        st.session_state[f'{rtype}_chosen_id'] = None
        sleep(0.4)
        st.rerun()