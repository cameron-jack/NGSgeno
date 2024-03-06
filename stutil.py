#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Utility functions for streamlit
"""

import streamlit as st

upper_info = "upper_info_viewer"
upper_height = "upper_info_height"
lower_info = "lower_info_viewer"
lower_height = "lower_info_height"

def init_state(key, value):
    """
    Initialises a session state with the value if the key is not in session state.
    Args:
        key (str): for the st.session_state dictionary
        value (str or None): value for the key in st.session_state
    """
    if key not in st.session_state:
        st.session_state[key] = value


def add_vertical_space(num_lines: int):
  """
  Add vertical space to your container
  Args:
      num_lines (int): number of lines to add vertical space
  """
  for line in range(num_lines):
    st.write("")

def hline():
   """
   Add a line across the page
   """
   st.write('***')


def custom_text(size, color, text, align="center", style="normal", padding='0px'):
    """
    Centred customised streamlit text
    Args:
        size (str): style size, either h1-6 or p
        color (str): can use hex # for a specific colour
        text (str): string to display
    Return:
        str: Text with the size and colour implemented for html markdown.
    """
    custom_css = f'<{size} style="text-align:{align};'+\
                 f'padding:{padding};'+\
                 f'font-style:{style};'+\
                 f'color:{color}">'+\
                 f'{text}</{size}>'

    st.markdown(custom_css, unsafe_allow_html=True)


def custom_button(color, text):
    """
    Customised streamlit button
    Args:
        color (str): colour of button (can have hex)
        text (str): text for button
    """
    st.markdown(f"""
    <style>
    .custom-button {{
        background-color: {color};
        color: white;
        padding: 0.25rem 0.75rem;
        margin: 8px 0;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        font-weight: 400;
        width: fit-content;
        height: auto
    }}
    .custom-button:hover {{
        opacity: 0.8;
    }}
    </style>
    <button class="custom-button">{text}</button>
    """, unsafe_allow_html=True)