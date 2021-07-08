#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: July 2021
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.8
@version_comment:
@last_edit:
@edit_comment:

Common code for connecting to the Musterer database web services. Should be extended to Rodentity when that becomes available.
"""

import os
import re
import csv
import glob
import json
import cgi
import cgitb
import requests
import subprocess
import collections
import sys
from configparser import ConfigParser as CP
from pathlib import PurePath

def _eppfn(bc):
    return "Plate{0}-EP.json".format(bc) # ear-punch plate (filename) format

mustfix = re.compile(r'(,\s*"(mouseBarcode|plateId)":\s*)(\d+)')

def getPlate_app(plateID, cache=True):
    """ Application version access Musterer ear-punch plate webservice for user ngs.genotyping
    Sets app.status and mb.showerror """
    global app
    config = CP()
    config.read(PurePath('../bin/config.ini'))
    url = config['earpunch']['url']
    username = config['earpunch']['username']
    passwd = config['earpunch']['passwd']
    service = 'exportPlate'
    fn = _eppfn(plateID)
    if cache and os.path.isfile(fn): # if cached in local file
        app.setstatus('reading cached', fn)
        with open(fn) as src:
            res = json.load(src)
    else:    
        param = 'plateId='+plateID
        url += "/{0}?".format(service)+param
        r = requests.get(url, auth=(username, passwd), verify=False)

        assert r.headers['content-type']=="application/json;charset=UTF-8"
        if not r.ok:
            app.setstatus("Return status =", r.status_code)
            app.setstatus("when accessing URL:", url)
            return None
        assert r.ok
        if not r.status_code==200:
           app.setstatus("Request for {service} failed. Return code={status}".format(service=service, status=r.status_code))
           return None
        # fix strings ... Musterer sends numbers but we really want strings - best to fix them here.
        global mustfix
        myjson = mustfix.sub(r'\1"\3"', r.text) # put number inside quotes
        res = json.loads(myjson)
        if len(res)==0:
            msg = "Musterer result for "+plateID+" is empty."
            app.setstatus(msg)
            mb.showerror("Musterer error", msg)
            return None
        if cache and len(res)==1: # save it in cache file
            with open(fn, "wt") as dst:
                dst.write(myjson)
            app.setstatus('cached', fn)
    if len(res)!=1 or str(res[0]['plateId'])!=plateID:
        px = ' '.join(str(p['plateId']) for p in res)
        msg = "Looking for plate {} but found plates {}.".format(plateID, px)
        app.setstatus(msg)
        mb.showerror("Musterer error.", msg)
        return None
    assert res[0]['plateId']==plateID # check that we got the right PlateID
    print('\n',res[0],'\n', file=sys.stderr)
    return res[0]  # return one plate from the result list


def getPlate(plateID, cache=True):
    """ Web-app version, access Musterer ear-punch plate webservice for user ngs.genotyping """
    service = 'exportPlate'
    config = CP()
    config.read(PurePath('../bin/config.ini'))
    url = config['earpunch']['url']
    username = config['earpunch']['username']
    passwd = config['earpunch']['passwd']
    fn = _eppfn(plateID)
    errs = []
    if cache and os.path.isfile(fn): # if cached in local file
        with open(fn) as src:
            res = json.load(src)
    else:    
        param = 'plateId='+plateID
        url += "/{0}?".format(service)+param
        r = requests.get(url, auth=(username, passwd), verify=False)

        assert r.headers['content-type']=="application/json;charset=UTF-8"
        if not r.ok:
            errs.append("Request for URL={} returned status={}".format(url, r.status_code))
            return None, errs
        assert r.ok
        if not r.status_code==200:
           errs.append("Request for Musterer {service} returned code={status}".format(service=service, status=r.status_code))
           return None, errs
        # fix strings ... Musterer sends numbers but we really want strings - best to fix them here.
        global mustfix
        myjson = mustfix.sub(r'\1"\3"', r.text) # put number inside quotes
        res = json.loads(myjson)
        if len(res)==0:
            msg = "Musterer result for ear-punch plate <code>"+plateID+"</code> is empty."
            errs.append(msg)
            return None, errs
        if cache and len(res)==1: # save it in cache file
            with open(fn, "wt") as dst:
                dst.write(myjson)
    if len(res)!=1 or str(res[0]['plateId'])!=plateID:
        px = ' '.join(str(p['plateId']) for p in res)
        msg = "Looking for plate {} but found plates {}.".format(plateID, px)
        errs.append(msg)
        return None, errs
    assert str(res[0]['plateId'])==plateID # check that we got the right PlateID
    print('\n',res[0],'\n', file=sys.stderr)
    return res[0], errs  # return one plate from the result list


def simplify_mouse_JSON(json_data):
    """
    Musterer mouse lookup web service returns JSON. Pick out the info we need and return:
    mouse_data = {}  # {mouseBarcode: {sex:M/F, strain:strainName, assays_names_values: {assayName:assayValue},
                     # sire: mouseID, dams:[mouseID]}
    from:
    [
       {
          "barcode":"105254500035",
          "strainId":52545,
          "strainName":"ASD552:Cahill:3::N3F4",
          "mouseId":35,
          "sex":"F",
          "bornDate":"09/04/2020",
          "mouseState":"Alive",
          "dams":[
             {
                "mouseId":14,
                "barcode":105139600013,
                "strainName":"ASD552:Cahill:3::N3F3",
                "strainId":51396
        },
         {
            "mouseId":14,
            "barcode":105139600018,
            "strainName":"ASD552:Cahill:3::N3F3",
            "strainId":51396
    }

    ],
          "sire":{
             "mouseId":14,
             "barcode":105139600014,
             "strainName":"ASD552:Cahill:3::N3F3",
             "strainId":51396
    },
          "assays":[
             {
                "assayMethod":"Amplifluor",
                "assayValueOptions":["wt/wt","mut/wt","mut/mut"],
                "assayId":2134,
                "assayName":"Ikbkb-V203I",
                "assayValue":"mut/wt"
    }
    ]
    }
    ]

    """
    mice = defaultdict(None)
    #print(json_data[0], file=sys.stderr)
    for bunch in json_data:
        for mouse in bunch:
            barcode = str(mouse.get('barcode'))
            sire = mouse.get('sire')
            dams = mouse.get('dams')
            assays = mouse.get('assays')
            m = {'sex' : mouse.get('sex'),
                'strain' : mouse.get('strainName')}
            if sire:
                m['sire_barcode'] = str(sire.get('barcode'))
            if dams:
                m['dams_barcodes'] = [str(d.get('barcode')) for d in dams]
            if assays:
                m['assay_names_values'] = {a.get('assayName'):a.get('assayValue') for a in assays}
                m['assay_value_options'] = {a.get('assayName'):a.get('assayValueOptions') for a in assays}
            mice[barcode] = m
    return mice


def get_musterer_mouse_info(mouse_barcodes, debug=False):
    """
        Contact the Musterer web service for mouse parental barcodes and Assay/Observables
        inputs: list of mouse barcodes
        intermediate: Send JSON blocks to simplify_mouse_JSON(json_data)
        output: dictionary of simplified mouse data
            mouse_data = {}  # {mouseID: {sex:M/F, strain:strainName, assays: [], sire: mouseID, dams:[mouseID]}
        Notes: No caching of parental information. Not sure how many mice we can ask for info in one transaction.
        https://musterer.apf.edu.au/musterer2/webservice/getMouseGenotypes?barcode=105254500035&barcode=105254500037
    """
    barcode_limit = 50 # max chars in GET request is 2048
    barcode_groups = [mouse_barcodes[i:i+barcode_limit] for i in range(0, len(mouse_barcodes), barcode_limit)]
    #print(barcode_groups)
    results = []
    config = CP()
    config.read(PurePath('../bin/config.ini'))
    url = config['mice']['url']
    username = config['mice']['username']
    passwd = config['mice']['passwd']
    for bg in barcode_groups:
        #print(bg)
        url += '?barcode='
        url += '&barcode='.join(bg)
        try:
            r = requests.get(url, auth=(username, passwd), verify=False)
        except requests.exceptions.ConnectionError:
            print('No connection to Musterer!', file=sys.stderr)
            return None

        assert r.headers['content-type']=="application/json;charset=UTF-8", r.headers['content-type']
        if not r.ok:
            #app.setstatus("Return status =", r.status_code)
            #app.setstatus("when accessing URL:", url)
            return None
        assert r.ok
        if not r.status_code==200:
           #app.setstatus("Request for {service} failed. Return code={status}".format(service=service, status=r.status_code))
           print("Request for {service} failed. Return code={status}".format(service=service, status=r.status_code), file=sys.stderr)
           return None
        res = r.json()
        if len(res)==0:
            msg = "Musterer returned no results for barcodes:" + ','.join(bg)
            #app.setstatus(msg)
            #mb.showerror("Musterer error", msg)
            if debug:
                print(msg, file=sys.stderr)
        results.append(res)
    return simplify_mouse_JSON(results)


def upload_plate(records):
    """
    Upload results using this web service:
    https://test.musterer.apf.edu.au/musterer2/genotyping/populateResults
    Need to test if this successfully uploads just one plate or several
    """
    config = CP()
    config.read(PurePath('../bin/config.ini'))
    pass