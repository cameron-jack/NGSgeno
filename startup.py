"""
@created: 3 Jan 2023
@author: Gabi Ryan, gabrielle.ryan@anu.edu.au

Initial starting interface for NGSgeno that checks for remote updates and the allows user
to accept any new changes as well as revert back to old versions.
"""

import streamlit as st
import os 
import sys
from subprocess import check_output, CalledProcessError, STDOUT
import subprocess

def change_ver(selection):
    reset_result = subprocess.call(["git", "reset",  "--hard", f"{selection}"])

def getstatusoutput(cmd):
    try:
        data = check_output(cmd, shell=True, universal_newlines=True, stderr=STDOUT)
        status = 0
    except CalledProcessError as ex:
        data = ex.output
        status = ex.returncode
    if data[-1:] == '\n':
        data = data[:-1]
    return status, data

def main():
  st.set_page_config(
        page_title="NGS Genotyping",
        page_icon="ngsg_icon.png"
        )
  #Page setup:
  _, version_col, _, select_col, change_col = st.columns([1, 4, 2, 3, 2])
  select_col.markdown(
        '<p style="color:#83b3c9; font-size: 90%">' +
        'Change versions:</p>',
        unsafe_allow_html=True
        )
  change_col.write('')
  change_col.write('')
  change_col.write('')
  st.image('ngsg_explorer.png')
  _, msg_col, _ = st.columns([1, 8, 1])
  update_msg = "The pipeline is up to date"
  update_color = '#5bc777'
  
  #Fetch data:

  #current tag
  current_status, current_ver = getstatusoutput("git describe")
  #all tags sorted by date
  all_status, all_ver = getstatusoutput("git tag --sort=-creatordate")
  latest_ver = all_ver.split("\n")[0]

  all_versions = []
  for ver in all_ver.split("\n"):
    all_versions.append(ver)
  
  #get the remote info from the git repo and then get status for any external updates
  fetch_status, fetch_data = getstatusoutput("git fetch")
  update_status, update_ver = getstatusoutput("git status")
  update_version = update_ver.split("\n")[1]

  #Update the display:

  version_col.markdown(
        '<p style="color:#4e839c; font-size: 110%; text-align:left">' +
        f'Current version:<br>{current_ver}</p>',
        unsafe_allow_html=True
        )

  #select previous tags
  selection_ver = select_col.selectbox(
        options=all_versions, label=' ', 
        label_visibility="collapsed"
        )
  
  change_ver_button = change_col.button(
        'Change', 
        on_click=change_ver, 
        args=[selection_ver]
        )

  #if there's a remote update to the repo
  if update_version != "Your branch is up to date with 'origin/master'." and latest_ver == current_ver:
    update_msg = "There is a new update to the pipeline"
    update_color = '#3d2ffa'

    st.markdown(
          '<p style="text-align:center;color:#5c8cfa;font-size:120%">' +
          'Do you want to download this update?</p>',
          unsafe_allow_html=True
          )    
    
    st.write('')
    _,update_col, no_update_col,_ = st.columns([4, 3, 3, 3])

    clicked_update = update_col.button(
          label='Yes, update', 
          help='Update the pipeline and run'
          )
    clicked_no_update = no_update_col.button(
          label='No, don\'t update', 
          help='Run the pipeline on the current version'
          )

    if clicked_update:
      update_result = subprocess.call(["git", "pull"])
      run_result = subprocess.call(["streamlit", "run", 'ngsgeno.py'])

    if clicked_no_update:
      run_result = subprocess.call(["streamlit", "run", 'ngsgeno.py'])

  #no updates
  else:
    _, run_col, _ = st.columns([2, 2, 1])
    run_pipeline = run_col.button('Start NGSgeno', type="primary")
    if run_pipeline:
      run_result = subprocess.call(["streamlit", "run", 'ngsgeno.py'])

  msg_col.markdown(
      f'<p style="text-align:center;color:{update_color};font-size:150%">{update_msg}</p>',
      unsafe_allow_html=True
      )


if __name__ == '__main__':
    main()
