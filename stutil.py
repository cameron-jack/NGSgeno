#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 1 May 2022
@author: Gabrielle Ryan, Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Utility functions for streamlit
"""

import streamlit as st
import sys
from time import sleep
from collections import defaultdict
import inspect

upper_info = "upper_info_viewer"
upper_height = "upper_info_height"
lower_info = "lower_info_viewer"
lower_height = "lower_info_height"

# universal message queue. Use this rather than st.session_state['message_queues']
mq = defaultdict(list)

def init_state(key, value):
    """
    Initialises a session state with the value if the key is not in session state.
    Args:
        key (str): for the st.session_state dictionary
        value (str or None): value for the key in st.session_state
    """
    if key not in st.session_state:
        st.session_state[key] = value


def set_state(key, value):
    st.session_state[key] = value


def flip_state(key):
    """ Use for dynamic checkboxes """
    if key in st.session_state:
        if st.session_state[key]:
            st.session_state[key] = False
        else:
            st.session_state[key] = True
    else:
        st.session_state[key] = True
        

def do_tm(message, level=None):
    if level is None:
        st.markdown(message)
    elif level.lower() == 'info':
        st.info(message)
    elif level.lower() == 'warning':
        st.warning(message)
    elif level.lower() == 'error':
        st.error(message)
    elif level.lower() == 'success':
        st.success(message)
    else:
        st.markdown(message)
        

def add_tm(message, level=None):
    """ DEPRECATED. adds a temporary message to the st.session_state['messages_temp'] list """
    if 'messages_temp' not in st.session_state:
        st.session_state['messages_temp'] = []
    st.session_state['messages_temp'].append((message, level))
    

def add_pm(message, level=None):
    """ Adds a persistent message to the st.session_state['messages_persist'] list """
    # will we ever use this?
    if 'messages_persist' not in st.session_state:
        st.session_state['messages_persist'] = []
    st.session_state['messages_persist'].append((message, level))


def m(message, level, dest=None, caller_id=None,
        size='p', color='black', align="center", font_style="normal", padding='0px', log_debug=False, no_log=False,
        call_func=None, call_line=None):
    """
    Standardised message interface. Removes the need for multiple function calls for the same message.
    Args:
        message (str): message text
        level is required {display, debug, info, warning, error, critical, begin, end, success, failure}
          display only goes to the screen
          failure is for expected non-error failures (user is doing something wrong)
        dest (None/tuple('log','persist','toast','css','mkdn','console','debug','noGUI','no_mkdn','silent','text')): display in these streams
          no_mkdn: just print on screen in regular text
          silent: no on-screen display
          text: don't use fancy formatting or display
        style (None/{'size':'p', 'color':'black', 'align':'center','style':'normal','padding':'0px'}):
            A dictionary of css stylings that will be applied if the css is a destination
        caller_id (str): if provided, store this message for later display
    #     log (bool): if True, send to exp.log()
    #     persist (bool): if True, make this a persistent message AS WELL as an immediate message
    #     css (bool): if True, apply css stylings and display as markdown
    #     mkdn (bool): True to show that messsage text is already markdown
    #     console (bool): if True, write to stdout as well
    #     debug (bool): if True, write to stderr and flush
    #     noscreen (bool): if True, do not write to the screen
    # The following args are passed to custom_text():
    #     size (str): style size, either h1-6 or p
    #     color (str): can use hex # for a specific colour
    #     text (str): string to display
        log_debug (False): don't log debug messages unless True
        no_log (False): don't log if True. Used by queued display messages
        call_func (str): name of calling function (used by queue)
        call_line (str): integer line number
    
    'debug' level messages are not logged unless log_debug is True
    'begin' and 'end' level message are not displayed to screen
    Automatically add dests('log', 'toast' and 'debug') for any warning, failure, error or critical message! (No dest required)
    """
    if dest is None:
        dest = {}
    dest = set([d.lower() for d in dest])
    level = level.lower()
    
    # do not add debug to the log
    if level == 'debug':
        dest.add('debug')
        if log_debug:
            dest.add('log')
    elif level == 'display':
        pass
    else:
        dest.add('log')
    
    if level == 'begin' or level == 'end':
        dest.add('noGUI')
        
    if level in set(['warning', 'error', 'critical', 'failure']):
        if 'debug' not in dest: 
            dest.add('debug')
        if 'toast' not in dest:
            dest.add('toast')
            
    try:
        func = sys._getframe(1).f_code.co_name
    except Exception as exc:
        print(f'Critical: {exc}', file=sys.stderr, flush=True)
    try:
        func_line = inspect.getframeinfo(sys._getframe(1)).lineno
    except Exception as exc:
        print(f'Critical: {exc}', file=sys.stderr, flush=True)
    try:
        call_func = sys._getframe(2).f_code.co_name
    except Exception as exc:
        print(f'Critical: {exc}', file=sys.stderr, flush=True)
    try:
        call_line = inspect.getframeinfo(sys._getframe(2)).lineno
    except Exception as exc:
        print(f'Critical: {exc}', file=sys.stderr, flush=True)
        
    if level != 'display':
        full_message = level.capitalize() + ': ' + message
            
    if 'console' in dest:
        print(f'{func=} {func_line} {call_func=} {call_line} {full_message}')
    if 'debug' in dest:
        print(f'{func=} {func_line} {call_func=} {call_line} {full_message}', file=sys.stderr, flush=True)
        
    if 'experiment' not in st.session_state or not st.session_state['experiment']:
        if 'log' in dest:
            st.error('Experiment not yet loaded, cannot log messages!')
            st.toast('Experiment not yet loaded, cannot log messages!')
            return
            
    # Most widgets will not handle markdown
    if level != 'display' and 'css' in dest:
        st.warning('Streamlit level displays are not compatible with css stylings')
        print('Streamlit level displays are not compatible with css stylings', file=sys.stderr)
    
    if 'log' in dest and not no_log:  # we were duplicating log entries from queued messages
        if 'mkdn' in dest:
            st.warning('Log cannot render markdown')
        if 'experiment' in st.session_state and st.session_state['experiment']:
            exp = st.session_state['experiment']
            if not level:
                level = 'info'
            
            exp.log(message, level=level, func=func, func_line=func_line, call_func=call_func, call_line=call_line)

    print(message, flush=True)
    if 'noGUI' in dest:
        print('do not write to GUI', flush=True)
        return

    if 'toast' in dest:
        st.toast(message)
            
    if 'css' in dest:
        # if style is None:
        #     style = {}
        # if 'size' not in style:
        #     style['size'] = 'p'
        # if 'color' not in style:
        #     style['color'] = 'black'
        # if 'align' not in style:
        #     style['align'] = 'center'
        # if 'style' not in style:
        #     style['style'] = 'normal'
        # if 'padding' not in style:
        #     style['padding'] = '0px'
        
        style = {}
        style['size'] = size
        style['color'] = color
        style['align'] = align
        style['style'] = font_style
        style['padding'] = padding
        message = custom_text(style['size'], style['color'], message, align=style['align'], style=style['style'], padding=style['padding'], display = False)    
        
    if caller_id:
        # init_state('message_queues', dict())
        # if caller_id not in st.session_state['message_queues']:
        #     st.session_state['message_queues'][caller_id] = []
        # st.session_state['message_queues'][caller_id].append((message, level))
        # return        
        mq[caller_id].append((message,level))
        return

    if 'persist' in dest:
        add_pm(message, level=level)

    if 'css' in dest:
        st.markdown(message, unsafe_allow_html=True)
        return
    
    if 'mkdn' in dest and not level:
        st.markdown(message)
        return
    
    if message.lower().startswith('critical:') or message.lower().startswith('failure:'):
        level = 'error'
    elif message.lower().startswith('success:') or message.lower().startswith('begin:') or\
            message.lower().startswith('end:'):
        level = 'info'

    if not level:
        level = 'info'

    if level.lower() in set(['info', 'success', 'begin', 'end']) and not 'no_mkdn' in dest:
        st.info(message)
    elif level.lower() == 'warning' and not 'no_mkdn' in dest:
        st.warning(message)
    elif level.lower() in set(['error', 'critical','failure']) and not 'no_mkdn' in dest:
        st.error(message)
    else:
        st.write(message)
    sleep(0.3)
            

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


def custom_text(size, color, text, align="center", style="normal", padding='0px', display=True):
    """
    Centred customised streamlit text
    Args:
        size (str): style size, either h1-6 or p
        color (str): can use hex # for a specific colour
        text (str): string to display
        display (bool): if True, display the markdown immediately, else return the string
    Return:
        str: Text with the size and colour implemented for html markdown.
    """
    custom_css = f'<{size} style="text-align:{align};'+\
                 f'padding:{padding};'+\
                 f'font-style:{style};'+\
                 f'color:{color}">'+\
                 f'{text}</{size}>'
    if display:
        st.markdown(custom_css, unsafe_allow_html=True)
    else:
        return custom_css
    

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