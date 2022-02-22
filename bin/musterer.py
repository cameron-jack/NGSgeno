#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: July 2021
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.15
@version_comment: Error handling updated
@last_edit: 2022-02-16
@edit_comment: Moved config path to global variable and location of file moved to library folder. Desktop application support removed.

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
from collections import defaultdict
import sys
from configparser import ConfigParser as CP
from pathlib import PurePath
import file_io  # guarding and unguarding of barcodes
from util import CONFIG_PATH, output_error

def _eppfn(bc):
    return "Plate-{0}-EP.json".format(bc) # ear-punch plate (filename) format

def _eppfn_old(bc):
    """ older files which are unguarded """
    return "Plate{0}-EP.json".format(file_io.unguard_pbc(bc)) # ear-punch plate (filename) format

mustfix = re.compile(r'(,\s*"(mouseBarcode|plateId)":\s*)(\d+)')


def get_plate(plateID, cache=True):
    """ Web-app version, access Musterer ear-punch plate webservice for user ngs.genotyping """
    try:
        service = 'exportPlate'
        config = CP()
        config.read(CONFIG_PATH)
        url = config['earpunch']['url']
        username = config['earpunch']['username']
        passwd = config['earpunch']['passwd']
        if not file_io.is_guarded_pbc(plateID):
            bc = file_io.guard_pbc(plateID)

        fn = _eppfn(plateID)
        errs = []
        if cache and os.path.isfile(fn): # if cached in local file
            print('Reading cached data:',fn, file=sys.stderr)
            with open(fn) as src:
                res = json.load(src)
        elif cache and os.path.isfile(_eppfn_old(plateID)):  # unguarded file name (older)
            fn = _eppfn_old(plateID)
            print('Reading cached data:',fn, file=sys.stderr)
            with open(fn) as src:
                res = json.load(src)
        else:  
            param = 'plateId='+file_io.unguard_pbc(plateID)
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
        if len(res)!=1 or str(res[0]['plateId'])!=file_io.unguard_pbc(plateID):
            px = ' '.join(str(p['plateId']) for p in res)
            msg = "Looking for plate {} but found plates {}.".format(file_io.unguard_pbc(plateID), px)
            errs.append(msg)
            return None, errs
        assert str(res[0]['plateId'])==file_io.unguard_pbc(plateID) # check that we got the right PlateID
        print(f"Returning {len(res[0]['wells'])} plate-well entries from Musterer\n", file=sys.stderr)
        return res[0], errs  # return one plate from the result list
    except Exception as exc:
        output_error(exc, msg='Error in musterer.get_plate')


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
    try:
        mice = defaultdict(None)
        #print(json_data[0], file=sys.stderr)
        for bunch in json_data:
            for mouse in bunch:
                barcode = str(mouse.get('barcode'))
                if file_io.is_guarded_mbc(barcode):
                    barcode = file_io.unguard_mbc(barcode)  # we need to look up unguarded barcodes
                sire = mouse.get('sire')
                dams = mouse.get('dams')
                assays = mouse.get('assays')
                m = defaultdict(None)
                m['sex'] = mouse.get('sex')
                m['strain'] = mouse.get('strainName')
                m['mouseId'] = mouse.get('mouseId')
                if sire:
                    m['sire_barcode'] = file_io.guard_mbc(str(sire.get('barcode')))
                if dams:
                    m['dams_barcodes'] = [file_io.guard_mbc(str(d.get('barcode'))) for d in dams]
                if assays:
                    m['assay_names_values'] = defaultdict(str)
                    m['assay_value_options'] = defaultdict(list)
                    for a in assays:
                        if 'assayMethod' in a and 'assayName' in a:
                            m['assay_names_values'][a.get('assayMethod')+'::'+a.get('assayName')] = a.get('assayValue')
                            m['assay_value_options'][a.get('assayMethod')+'::'+a.get('assayName')] = a.get('assayValueOptions')
                mice[file_io.guard_mbc(barcode)] = m
        return mice
    except Exception as exc:
        output_error(exc, msg='Error in musterer.simplify_mouse_JSON')


def get_musterer_mouse_info(mouse_barcodes, debug=False):
    """
        Contact the Musterer web service for mouse parental barcodes and Assay/Observables
        inputs: list of mouse barcodes
        intermediate: Send JSON blocks to simplify_mouse_JSON(json_data)
        output: dictionary of simplified mouse data
            mouse_data = {}  # {mouse_barcode: {sex:M/F, strain:strainName, assays: [], sire: mouse_barcode, dams:[mouse_barcode]}
        Notes: No caching of parental information. Not sure how many mice we can ask for info in one transaction.
        https://musterer.apf.edu.au/musterer2/webservice/getMouseGenotypes?barcode=105254500035&barcode=105254500037
    """
    try:
        # take off the guards for DB lookups
        mouse_barcodes = [file_io.unguard(mbc) if file_io.is_guarded(mbc) else mbc for mbc in mouse_barcodes]
        barcode_limit = 50 # max chars in GET request is 2048
        barcode_groups = [mouse_barcodes[i:i+barcode_limit] for i in range(0, len(mouse_barcodes), barcode_limit)]
        #print(barcode_groups)
        results = []
        config = CP()
        config.read(CONFIG_PATH)
        username = config['mice']['username']
        passwd = config['mice']['passwd']
        for bg in barcode_groups:
            #print(bg)
            url = config['mice']['url'] + '?barcode=' + '&barcode='.join(bg)
            try:
                #print(url, username, passwd, file=sys.stderr)
                r = requests.get(url, auth=(username, passwd), verify=False)
            except requests.exceptions.ConnectionError:
                print('No connection to Musterer!', file=sys.stderr)
                return None

            #assert r.headers['content-type']=="application/json;charset=UTF-8", r.headers['content-type']
            #assert r.headers['content-type']=="text/html;charset=UTF-8", r.headers['content-type']
            if not r.ok:
                #app.setstatus("Return status =", r.status_code)
                #app.setstatus("when accessing URL:", url)
                print('Return status:', r.status_code, 'when accessing URL:', url, file=sys.stderr)
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
                print(msg, file=sys.stderr)
            results.append(res)
        simple_mice = simplify_mouse_JSON(results)
        print('Returned', len(simple_mice), 'from', len(mouse_barcodes), 'lookups', file=sys.stderr)
        return simplify_mouse_JSON(results)
    except Exception as exc:
        output_error(exc, msg='Error in musterer.get_musterer_mouse_info')


def upload_plate(records):
    """
    Upload results using this web service:
    https://test.musterer.apf.edu.au/musterer2/genotyping/populateResults
    Need to test if this successfully uploads just one plate or several
    """
    config = CP()
    config.read(CONFIG_PATH)
    pass