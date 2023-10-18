#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Utility functions for streamlit
"""

import streamlit as st

def add_vertical_space(num_lines: int):
  """
  Add vertical space to your container
  Args:
      num_lines (int): number of lines to add vertical space
  """
  for line in range(num_lines):
    st.write("")


def custom_text(size, color, text, align="center", style="normal"):
    """
    Centred text
    Args:
        size (str): style size, either h1-6 or p
        color (str): can use hex # for a specific colour
        text (str): string to display
    Return:
        str: Text with the size and colour implemented for html markdown.
    """
    custom_css = f'<{size} style="text-align:{align};'+\
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