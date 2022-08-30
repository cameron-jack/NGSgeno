#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.16
@version_comment: New streamlit interactive web interface
@last_edit: 2022-08-29
@edit_comment:

A web based interactive GUI with Streamlit. Plate and sample barcodes here are unguarded.
In all other code they must be guarded. We guard them here before we send them to any external function.
"""

import os
import sys
from pathlib import PurePath
import itertools
from math import fabs, factorial, floor, ceil
from io import StringIO

import pandas as pd

import streamlit as st
import streamlit.components.v1 as components

from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode
from bin.experiment import Experiment, EXP_FN, load_experiment
from bin.util import output_error
from bin.util import CAP_VOLS, DEAD_VOLS
import bin.file_io as file_io
import bin.db_io as db_io
from bin.makehtml import generate_heatmap_html
#from load_css import local_css
from ngsgeno_gr import *

credits="""
@created: March 2022
@author: Cameron Jack, ANU Bioinformatics Consultancy, 2019-2021
@version: 0.16
@version_comment: New interactive web application interface
@last_edit: 29/08/2022 by Gabi Ryan 
@edit_comment: Separate pages for the Home and Pipeline, both reference back to ngsgeno (or ngsgeno_gr) for functions.

Web application style interface using Streamlit. Present information as both a dashboard and a workflow in one page.

The GUI interacts with a single Experiment object at one time. Methods are called on this to activate pipeline
functionality. The Experiment then deals directly with the pipeline logic.
"""


### Lay out the pipeline interface sections
# 1. Run folder (includes sidebar)
# 2. Inventory (plates and assays + DB lookup buttons + remaining i7i5 barcodes)
# 3. Nimbus progress
# 4. Primer plates + taq/water
# 5. Add spare samples for barcoding + barcode completion status (plates done/not done)
# 6. Miseq data available (notice + instructions if not yet complete)
# 7. Analysis stages
###
 

def main():
    st.set_page_config(page_icon="ngsg_icon.png", page_title="NGS Genotyping Pipeline", initial_sidebar_state="expanded", layout="wide")
    # Remove whitespace from the top of the page and sidebar
    st.markdown(
        """
        <style>
               .css-18e3th9 {
                    padding-top: 0rem;
                    padding-bottom: 10rem;
                    padding-left: 5rem;
                    padding-right: 5rem;
                }
                header {visibility: hidden;}
                footer {visibility: hidden;}
        </style>
        """, unsafe_allow_html=True)
    padding = 0
    st.markdown(f""" <style>
    .reportview-container .main .block-container{{
        padding-top: {padding}rem;
        padding-right: {padding}rem;
        padding-left: {padding}rem;
        padding-bottom: {padding}rem;
    }} </style> """, unsafe_allow_html=True)
    #This is a workaround from https://www.codegrepper.com/code-examples/python/streamlit+sidebar+width
    st.markdown(
        """
        <style>
        div.block-container {padding-top:0rem;}
        [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
            width: 400px;
            top: 0rem;
            padding-top: 0rem;
            height: 100vh;
        }
        [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
            width: 400px;
            top: 0rem;
            padding-top: 0rem;
            margin-left: -400px;
        }
         
        </style>
        """,
        unsafe_allow_html=True)


    #TO DO:
    #Taq and water plates. 
    #Load data once rather than each time a script is run - "refresh library" button?

    
    #Sidebar components for home page
    sidebar_display(page=1)

    if not st.session_state['experiment']:
        with st.container():
            st.write('')
            st.write('')
            st.markdown('<h2 style="text-align:center;color:#154c79">Please load existing experiment or create a new one</h2>', unsafe_allow_html=True)
    else:
        title='Experiment '+ st.session_state['experiment'].name
        st.markdown(f'<h2 style="color:#2B4864">{title}</h2>', unsafe_allow_html=True)

        #Options for interacting with the data:
        data_tab1,data_tab2=st.tabs(['Load', 'View'])

        with data_tab2:
                st.markdown('<h2 style="text-align:center;color:#0f6b8e">View Data</h2>', unsafe_allow_html=True)
                data_table(table_options=True,show_usage=True, show_primers=True)
                extra_data()

     
        with data_tab1:
                st.markdown('<h2 style="text-align:center;color:#0f6b8e">Load Sample Data</h2>', unsafe_allow_html=True)
                load_data()
          
                    
        

    #hvar = """  <script>
    #                var elements = window.parent.document.querySelectorAll('.streamlit-expanderHeader');
    #                elements.forEach(function (element) {
    #                    element.style.fontSize = 'large';
    #                    });         
    #            </script>"""
    #components.html(hvar, height=0, width=0)

    

    #st.markdown(
    #    """
    #<script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    #<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    #<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    #""",
    #    unsafe_allow_html=True,
    #)


if __name__ == '__main__':
    main()























