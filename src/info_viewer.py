import display_components as dc
import streamlit as st
import jsonpickle
import pandas as pd
import util
from stutil import set_state, mq, m, sleep, add_vertical_space, hline, unlocked
from makehtml import generate_heatmap_html
import streamlit.components.v1 as components

def manage_delete_cb(caller_id, category, ids):
    """
    Callback for deletion operations
    args:
        caller_id (str): the name of a message queue for a particular display widget
        category (str): file/plate/group - the resource type to be deleted
        ids (list): a list of items to be deleted
    """
    exp = st.session_state['experiment']
    successful_ids = []
    failed_ids = []

    if category == 'file':
        for id in ids:
            if id not in exp.uploaded_files:
                continue
            success = exp.del_file_record(id)
            if success:
                successful_ids.append(id)
                m(f'{id} removed', level='info')
            else:
                failed_ids.append(id)
                m(f'{id} could not be removed', level='warning')
            st.session_state['previous_file_delete_selection'] = ids
    elif category == 'plate':
        for pid in ids:
            gid = util.guard_pbc(pid, silent=True)
            if gid not in exp.plate_location_sample:
                failed_ids.append(gid)
            else:
                success = exp.delete_plate(gid)
                if success:
                    successful_ids.append(gid)
                else:
                    failed_ids.append(gid)
        st.session_state['previous_plate_delete_selection'] = ids
    elif category == 'group':  # from summary
        for pid in ids:
            dest_pid = util.guard_pbc(pid, silent=True)
            if dest_pid in exp.plate_location_sample and exp.plate_location_sample[dest_pid]['purpose'] == 'amplicon':
                success = exp.delete_plate(dest_pid)
                if success:
                    successful_ids.append(dest_pid)
                else:
                    failed_ids.append(dest_pid)
            elif dest_pid not in exp.dest_sample_plates:
                m(f"{dest_pid=} doesn't actually exist in the experiment!", level='error')
                failed_ids.append(dest_pid)
            else:
                sample_pids = exp.dest_sample_plates[dest_pid]
                delete_pids = sample_pids + [dest_pid]
                for dpid in delete_pids:
                    success = exp.delete_plate(dpid)
                    if success:
                        successful_ids.append(dpid)
                    else:
                        failed_ids.append(dpid)
        st.session_state['previous_group_delete_selection'] = ids
    # set up messages
    for sid in successful_ids:
        m(f'{sid} removed', level='success', caller_id=caller_id)
    for fid in failed_ids:
        m(f'{fid} could not be removed', level='error', caller_id=caller_id)
    return True


def cancel_delete(category, ids):
    """
    Callback for delete operations
    """
    if category == 'plate':
        st.session_state['previous_plate_delete_selection'] = ids
    elif category == 'file':
        st.session_state['previous_file_delete_selection'] = ids
    elif category == 'group':  # from summary
        st.session_state['previous_group_delete_selection'] = ids
    return True


def set_selection(set_key, widget_key):
    """ Callback function for selectbox display elements """
    if widget_key in st.session_state:
        st.session_state[set_key] = st.session_state[widget_key]


def info_viewer(selection, key, dna_pids=None, view_height=350):
    """
    Module for displaying info viewers, which are aggrid tables with information about the experiment. 
    These are displayed in a horizontal layout and can be used to view information about the samples, 
    consumables, files, plates, primers, indexes and log. 
    Each viewer has a selection box to choose what information to display and a number input to set 
    the height of the display grid. The viewers are displayed side-by-side if multiple are selected. 
    The viewers are updated when the selection or height is changed.
    """
    exp = st.session_state['experiment']
    caller_id = 'info_viewer_'+str(key)
    if selection == 'Samples':
        infoview_samples(key=key+selection, height=view_height, caller_id=caller_id)

    if selection == 'Consumables':
        infoview_consumables(key=key+selection, height=view_height, caller_id=caller_id)

    if selection == "Status":
        # Status tab should tell us where we are up to in the pipeline and what's happened so far
        infoview_status(key=key+selection, height=view_height, caller_id=caller_id)

    if selection == "Files":
        file_usage = exp.get_file_usage()
        infoview_files(key+selection, file_usage, height=view_height, caller_id=caller_id)

    if selection == "Plates":
        plate_usage = exp.get_plate_usage()
        infoview_plates(key+selection, plate_usage, height=view_height, caller_id=caller_id)

    if selection == "Plate Viewer":
        view_height = 500
        infoview_platelayout(key+selection, height=view_height, caller_id=caller_id)

    if selection == "Primers":
        infoview_primers(key+selection, dna_pids=dna_pids, height=view_height, caller_id=caller_id)

    if selection == "Indexes":
        infoview_indexes(key+selection, dna_pids=dna_pids, height=view_height, caller_id=caller_id)

    if selection == "Log":
        infoview_log(key+selection, height=view_height, caller_id=caller_id)


def upper_info_viewer_panel(tab_col3, tab_col2, widget_key, default_view1='None', default_view2='None', checked=True):
    with tab_col3:
        add_vertical_space(2)
        show_upper_info_viewer_checkbox(widget_key, value=checked)
    if st.session_state['show_upper_info_viewer']:
        with tab_col2:
            ignore = info_selection(widget_key+"top_viewer", 'info_panel1', 'info_panel2',
                    'upper_panel_height', default_view1=default_view1, default_view2=default_view2,
                    default_height=st.session_state.get('upper_panel_height',250))
    caller_id = 'display_feedback'
    # display any messages for this widget
    if caller_id in mq:
        for msg, lvl in mq[caller_id]:
            m(msg, level=lvl, no_log=True)
        sleep(0.3)
        mq[caller_id] = set()


def show_upper_info_viewer_checkbox(widget_key, value=True): #, default_panel1='None', default_panel2='None'):
    """
    Allows the user to turn the upper info viewer panel on and off
    The bottom info viewer display is always on
    Allows the default display panels to be set (page specific content)
    args:
        widget_key (str): a unique identifier to prevent same-page clashes
        value (bool): the default value of the enabling checkbox
    """
    if 'show_info_viewer' not in st.session_state:
        st.session_state['show_upper_info_viewer'] = True

    if st.checkbox('Info Viewer', value=value, key=widget_key):
        st.session_state['show_upper_info_viewer'] = True
        # if default_panel1:
        #     st.session_state['info_panel1'] == default_panel1
        # if default_panel2:
        #     st.session_state['info_panel2'] == default_panel2
    else:
        st.session_state['show_upper_info_viewer'] = False


def info_selection(key, view1_key, view2_key, height_key, default_view1="None",
        default_view2="None", default_height=250):
    """
    Container for displaying module info functions, each of which provides a dataframe for display in an aggrid.
    Because aggrid allows selection, each module can also handle a standard set of operations (such as delete).
    Function handles a single case and extra keys are needed for correct naming and identification
    key - the general key applied to this set of widgets
    view1_key - the key for selection1 lookups in other code
    view2_key - the key for selection2 lookups in other code
    height_key - the key for height lookups in other code
    defaults are returned, otherwise these are set via callbacks
    """
    key=str(key)

    options = ["None","Samples", "Consumables", "Status", "Files", "Plates", "Plate Viewer",
            "Primers", "Indexes", "Log"]

    disp_col1, disp_col2, height_col = st.columns([4,4,2])

    with disp_col1:
        select1_key = key+"_select1"
        selection1 = st.selectbox("Choose info to view", options=options, placeholder='',
                index=options.index(default_view1), on_change=set_selection,
                args=[view1_key, select1_key], key=select1_key)
        set_state(view1_key, selection1)

    with disp_col2:
        select2_key = key+"_select2"
        selection2 = st.selectbox("Choose info to view", options=options, placeholder='',
                index=options.index(default_view2), on_change=set_selection,
                args=[view2_key, select2_key], key=select2_key)
        set_state(view2_key, selection2)

    with height_col:
        height_widget_key = key+"_height"
        view_height = st.number_input('Set display height', min_value=50, max_value=700,
                value=default_height, step=25, help="Size of display grid", on_change=set_selection,
                args=[height_key, height_widget_key], key=height_widget_key)
        set_state(height_key, view_height)

    return True


def show_info_viewer(selection, height, groupkey):
    """ Display any number of info view windows side-by-side. Select must be an iterable container """
    if len(selection) != 0:
        columns = st.columns(len(selection))
        for i in range(len(selection)):
            with columns[i]:
                with st.spinner('loading info...'):
                    info_viewer(selection[i], str(groupkey)+selection[i]+str(i), view_height=height)

### Info viewer functions - these take care of displaying the dataframes in the info viewers, and any associated messages
# These can used without the infoviewer frame or selection functions
# references, log, plate layout, files, samples, consumables, status, primers, indexes, plates

def infoview_references(key, height=350, subset='all', caller_id='display_feedback'):
    """
    Show the amplicon targets/reference sequences that are loaded
    Args:
        key (str): key for the component
        height (int): height of the table
        subset (str): 'all','amplicon','custom','rodentity' to show all or only the selected purpose
        caller_id (str): ID for the caller, used for messages
    Messages are displayed by the component holding these widgets

    NOTE: Currently unused
    """
    exp = st.session_state['experiment']
    dataset = []
    for fn,purp in exp.reference_sequences:
        if subset == 'amplicon' and purp != 'amplicon_reference':
            dataset.extend(exp.reference_sequences[(fn,purp)])
        elif subset == 'custom' and purp != 'custom_reference':
            # custom references are not used in the pipeline, we should report a mistake
            dataset.extend(exp.reference_sequences[(fn,purp)])
            m(f'Custom reference sequences are not used in the pipeline, please check the calling code',level='critical')
        elif subset == 'rodentity' and purp != 'rodentity_reference':
            dataset.extend(exp.reference_sequences[(fn,purp)])
        elif subset == 'all':
            dataset.extend(exp.reference_sequences[(fn,purp)])
    ref_df = pd.DataFrame(dataset, columns=['Name','Sequence'])
    st.markdown('**'+group+'**')
    dc.aggrid_interactive_table(ref_df, key=str(key)+'seqs'+group, grid_height=height)


def infoview_log(key, height=250, caller_id='display_feedback'):
    """
    Display the log
    Messages are displayed by the component holding these widgets
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    log_entries = st.session_state['experiment'].get_log(100)
    if len(log_entries) == 0:
        st.write('No entries currently in the log')
    else:
        df = pd.DataFrame(log_entries, columns=exp.get_log_header())
        df = df.drop(['Calling function', 'Call line'], axis = 1)
        dc.aggrid_interactive_table(df, grid_height=height, key=str(key)+'logs')


def infoview_platelayout(key, height=500, caller_id='display_feedback'):
    """
    Visual view of the plates in the experiment
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    plate_ids = []
    for pid in exp.plate_location_sample:
        plate_ids.append(f"{exp.plate_location_sample[pid]['purpose']} plate: {util.unguard(pid, silent=True)}")

    #Let user choose which plate to view
    # _, col1, _ = st.columns([1,2,1])
    plate_selectbox = st.selectbox('Plate ID to view', plate_ids, key=str(key)+'plate_viewer')
    if plate_selectbox is not None:
        plate_id = util.guard_pbc(plate_selectbox.split(':')[1])
        if plate_id in exp.plate_location_sample:
            jsonpickle_plate = jsonpickle.encode(exp.get_plate(plate_id), keys=True)
            heatmap_str = generate_heatmap_html(jsonpickle_plate, plate_id, scaling=0.9)
            with open("makehtml.html", 'wt', encoding="utf-8") as outf:
                print(heatmap_str, file=outf)
            components.html(heatmap_str, height=height, scrolling=True)
        else:
            plate_barcode_error_msg = "Plate barcode not found in experiment"
            st.markdown(f'<p style="color:#FF0000">{plate_barcode_error_msg}</p>',
                        unsafe_allow_html=True)


def infoview_plates(key, plate_usage, height=300, caller_id=None):
    """
    Info bar display of plates
    Use only the provided caller_id
    """
    exp = st.session_state['experiment']
    caller_id = 'display_plates'
    plate_df = pd.DataFrame(plate_usage, columns=['Plates', 'Num Wells', 'Purpose'])
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        st.write('No plates loaded')
    else:
        selection = dc.aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'plate_aggrid')
        if selection is not None and 'selected_rows' in selection:
            if selection['selected_rows'] is not None:
                pids = [pid for pid in selection['selected_rows'].get('Plates') if pid is not None]
                gids = [util.guard_pbc(pid, silent=True) for pid in pids]
                gids = [gid for gid in gids if gid in exp.plate_location_sample]
                if gids:
                    if 'previous_plate_delete_selection' not in st.session_state:
                        st.session_state['previous_plate_delete_selection'] = None
                    if st.session_state['previous_plate_delete_selection'] == pids:
                        st.session_state['previous_plate_delete_selection'] = None
                    else:
                        if pids != st.session_state['previous_plate_delete_selection']:
                            st.markdown(f"**You selected {pids}**")
                        delbox = st.container() # doesn't work reliably in st1.26
                        del_col1, del_col2, del_col3, _ = delbox.columns([2,1,1,4])
                        del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                        del_col2.button("Yes",on_click=manage_delete_cb,
                                args=(caller_id,'plate',pids), key="delete " + str(key), help=f"Delete {pids}")
                        del_col3.button("No", on_click=cancel_delete,
                                args=('plate',pids), key="keep " + str(key), help=f"Keep {pids}")

        
def infoview_files(key, file_usage, height=250, caller_id='display_feedback'):
    """
    Display info for all files that have so far been uploaded into the experiment
    Give info on name, any plates they contain, whether they are required so far, etc
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    file_df = pd.DataFrame.from_dict(file_usage, orient='index')
    if file_df is None or not isinstance(file_df, pd.DataFrame):
        st.write('No plates loaded')
        return

    file_df.reset_index(inplace=True)
    #file_df = file_df.rename(columns = {'index':'File', 'plates':'Plates', 'purpose':'Purpose'})
    file_df = file_df.rename(columns = {'index':'File', 'date modified':'Date modified', 'purpose':'Purpose'})
    if not file_df.empty:
        file_df.sort_values(by='Date modified', inplace=True, ascending=False)

    selection = dc.aggrid_interactive_table(file_df, grid_height=height, key=str(key)+'file_aggrid')
    if selection is not None:
        if 'selected_rows' in selection and selection['selected_rows'] is not None:
            fns = [fp for fp in selection['selected_rows'].get('File')]
            if fns:
                if 'previous_file_delete_selection' not in st.session_state:
                    st.session_state['previous_file_delete_selection'] = None
                if st.session_state['previous_file_delete_selection'] == fns:
                    st.session_state['previous_file_delete_selection'] = None
                else:
                    if fns != st.session_state['previous_file_delete_selection']:
                        st.markdown(f"**You selected {fns}**")
                    delbox = st.container()
                    del_col1, del_col2, del_col3, _ = delbox.columns([2,1,1,4])
                    del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                    del_col2.button("Yes",on_click=manage_delete_cb,
                            args=(caller_id,'file',fns), key="delete " + str(key), help=f"Delete {fns}")
                    del_col3.button("No", on_click=cancel_delete,
                            args=('file',fns), key="keep " + str(key), help=f"Keep {fns}")
                    selection = None


def infoview_status(key, height=300, caller_id='display_feedback'):
    """
    Display the progress in the pipeline for this experiment
    Should use aggrid to display the stages and the changes at each stage
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    steps, header = exp.get_stages()
    status_df = pd.DataFrame(steps, columns=header)
    #status_df.reset_index(inplace=True)
    #status_df = status_df.rename(columns = {'index':'Steps', 'pending':'Pending Steps'})
    status_df = dc.aggrid_interactive_table(status_df, grid_height=height, key=str(key)+'status_aggrid')


def infoview_consumables(key, height=300, caller_id=None):
    """
    info bar display of "consumables" non-sample plate info
    Use only the provided caller_id
    summarise_consumables():
        d = {'taqwater_pids_pcr1':[], 'taqwater_pids_pcr2':[], 'taq_vol_pcr1':0, 'taq_vol_pcr2':0,'water_vol_pcr1':0,
                'water_vol_pcr2':0, 'primer_pids':[], 'primer_count_ngs':0, 'primer_count_custom':0, 'unique_primers':set(),
                'primer_well_count':0, 'assay_primer_mappings':0, 'rodentity_reference_files':[],
                'custom_reference_files':[], 'amplicon_reference_files':[],
                'index_pids':[], 'unique_i7s':set(), 'unique_i5s':set()}
    """
    exp = st.session_state['experiment']
    caller_id = 'display_consumables'
    # display all the required files whether they are, or are not present
    # Plates: primer, index, taq/water, PCR, references, primer/assaylist,
    headers = ['Purpose','Barcode/ID','Wells','Type','Entries']
    consumables = exp.summarise_consumables()

    data_rows = []
    for tp in consumables['taqwater_pids_pcr1']:
        data_rows.append(['taq/water (PCR 1)', util.unguard_pbc(tp, silent=True), 6, util.PLATE_TYPES['Echo6'], 6])
    for tp in consumables['taqwater_pids_pcr2']:
        data_rows.append(['taq/water (PCR 2)', util.unguard_pbc(tp, silent=True), 6, util.PLATE_TYPES['Echo6'], 6])
    for pp in consumables['primer_pids']:
        data_rows.append(['primer', util.unguard_pbc(pp, silent=True), 384, util.PLATE_TYPES['Echo384'],
                len(exp.plate_location_sample[pp]['wells'])])
    for ip in consumables['index_pids']:
        data_rows.append(['index', util.unguard_pbc(ip, silent=True), 384, util.PLATE_TYPES['Echo384'],
                len(exp.plate_location_sample[ip]['wells'])])
    for f in consumables['rodentity_reference_files']:
        data_rows.append(['Rodentity references', f, 0, 'File', exp.count_reference_sequences('rodentity_references', silent=True)])
    for f in consumables['custom_reference_files']:
        data_rows.append(['Custom references', f, 0, 'File', exp.count_reference_sequences('custom_references', silent=True)])
    for f in consumables['amplicon_reference_files']:
        data_rows.append(['Amplicon references', f, 0, 'File', exp.count_reference_sequences('amplicon_references', silent=True)])
    for f in exp.uploaded_files:
        if exp.uploaded_files[f].get('purpose','') == 'assay_primer_map':
            data_rows.append(['assay-primer mappings', f, 0, 'File', consumables['assay_primer_mappings']])
    plate_df = pd.DataFrame(data_rows, columns=headers)
    if plate_df is None or not isinstance(plate_df, pd.DataFrame):
        st.write('No plates loaded')
    else:
        selection = dc.aggrid_interactive_table(plate_df, grid_height=height, key=str(key)+'consumables_aggrid')


def infoview_samples(key, height=250, caller_id=None):
    """
    Info bar display a summary of all loaded DNA, amplicon, and sample plates
    Use only the provided caller_id
    """
    exp = st.session_state['experiment']
    selection = []
    df = exp.inputs_as_dataframe()
    if df is None or not isinstance(df, pd.DataFrame):
        st.markdown('**No 384-well DNA plate data loaded**')
    else:
        selection = dc.aggrid_interactive_table(df, key=key, grid_height=height)
        if selection is not None:
            if 'selected_rows' in selection and selection['selected_rows'] is not None:
                rows = [r for r in selection["selected_rows"].get('DNA/amplicon PID') if r != 'Total']
                if rows:
                    if 'previous_group_delete_selection' not in st.session_state:
                        st.session_state['previous_group_delete_selection'] = None
                    if st.session_state['previous_group_delete_selection'] != rows:
                        # only do the code below if this is a fresh selection
                        lines = '\n'.join(['DNA/amplicon PID: '+r for r in rows])
                        st.markdown(f"**You selected {lines}**")
                        del_col1, del_col2, del_col3, _ = st.columns([2,1,1,4])
                        del_col1.markdown('<p style="color:#A01751">Delete selection?</p>', unsafe_allow_html=True)
                        del_col2.button("Yes",on_click=manage_delete_cb,
                                args=(caller_id, 'group',rows), key="delete " + str(key), help=f"Delete {lines}")
                        del_col3.button("No",on_click=cancel_delete, args=('group',rows),
                                key="keep " + str(key), help=f"Keep {lines}")
    selection = None


def infoview_primers(key, dna_pids=None, primer_pids=None, height=350, save_buttons=False, caller_id='display_feedback'):
    """
    Alternative display_primer_components using aggrid
    Designed as a component for a generic display widget and as a standalone with optional save to file
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    ul_conv = 1000
    if not dna_pids:
        dna_pids = exp.get_dna_pids(caller_id=caller_id)
    if not primer_pids:
        primer_pids = exp.get_primer_pids(caller_id=caller_id)
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids, caller_id=caller_id)
    primer_avail_vols_doses = exp.get_available_primer_vols_doses(pmr_pids=primer_pids, caller_id=caller_id)
    primer_wells = exp.get_available_primer_wells(pmr_pids=primer_pids, caller_id=caller_id)
    pmr_pids_wells = {pmr:{} for pmr in primer_wells}
    primer_info_array = []
    for primer in exp.primer_assayfam:
        req_vol = exp.transfer_volumes['PRIMER_VOL']*primer_usage.get(primer,0)  # nl
        req_doses = primer_usage.get(primer,0)
        req_wells = util.num_req_wells(req_vol)
        avail_vol = primer_avail_vols_doses.get(primer,[0,0])[0]  # nl
        avail_doses = primer_avail_vols_doses.get(primer,[0,0])[1]
        avail_wells = len([pmr_info for pmr_info in primer_wells.get(primer, [])])
        if req_wells == 0 and avail_wells == 0:
            continue
        if primer not in primer_wells:
            pos_line = ''
        else:
            for ppid, well, vol, dose in primer_wells[primer]:
                if ppid not in pmr_pids_wells[primer]:
                    pmr_pids_wells[primer][ppid] = []
                pmr_pids_wells[primer][ppid].append(well)
            ppid_well_lines = [f"{util.unguard_pbc(ppid, silent=True)}:\
                    {','.join(pmr_pids_wells[primer][ppid])}" for ppid in pmr_pids_wells[primer]]
            pos_line = ' '.join(ppid_well_lines)
        info = [primer, req_doses, req_vol/ul_conv, req_wells, avail_doses, avail_vol/ul_conv,
                avail_wells, pos_line]
        primer_info_array.append(info)

    primer_df = pd.DataFrame(primer_info_array, columns=['Primer', 'Required Doses',
            'Required Volume (μL)', 'Required Wells', 'Available Doses', 'Available Volume (μL)', 'Available Wells', 'Positions'])
    primer_table = dc.aggrid_interactive_table(primer_df, grid_height=height, key=str(key)+'primer_display')
    primer_csv = primer_df.to_csv(index=False) #.encode('utf-8')
    primer_list_fn = exp.get_exp_fn('primer_list.csv', caller_id=caller_id)
    if 'primer_csv' not in st.session_state:
        st.session_state['primer_csv'] = 0
    if st.session_state['primer_csv'] != primer_csv:
        st.session_state['primer_csv'] = primer_csv
        primer_df.to_csv(path_or_buf=primer_list_fn, index=False) #.encode('utf-8')
    st.write(f'Primer table written to file: {primer_list_fn}')


def infoview_indexes(key, dna_pids=None, height=350, caller_id='display_feedback'):
    """
    Display info for all indexes that have been uploaded into the experiment
    Messages are displayed by the component holding these widgets
    """
    exp = st.session_state['experiment']
    assay_usage, primer_usage = exp.get_assay_primer_usage(dna_pids)
    fwd_idx, rev_idx = exp.get_index_avail()
    indexes = {**fwd_idx, **rev_idx}

    index_df = pd.DataFrame.from_dict(indexes,  orient='index')
    index_df.reset_index(inplace=True)
    index_df = index_df.rename(columns = {'index':'Index',
                                          'well_count':'Wells',
                                          'avail_transfers':'Available transfers',
                                          'avail_vol':'Available Volume (μL)'})
    index_table = dc.aggrid_interactive_table(index_df,grid_height=height,key=str(key)+'index_display')