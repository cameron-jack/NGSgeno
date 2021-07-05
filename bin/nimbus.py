#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: Oct 2019
@author: Bob Buckley & Cameron Jack, ANU Bioinformatics Consultancy, JCSMR, Australian National University
@version: 0.7
@version_comment: column header changes
@last_edit: 2021-07-02
@edit_comment: epwell->EPwell, epplate->EPplate, Sample No->mouseID in reporting

GUI for Nimbus files for the NGSgenotyping project.
Note - The GUI code is no longer used - the nimbus() function is called by cgi-nimbus2.py
in the new HTML interface.

Input is 1-4 96-well sample (ear-punch) plate bar codes ... and a destination 384-well plate bar code.

The program uses either a Musterer webservice for each ear-punch plate, or has 
the required info passed in from CSV format and converted to JSON.
The webservice returns JSON containing a list of wells and mouse information 
(IDs, barcodes, expected "observables", ...).

The output format is CSV
The output asks the Nimbus Robot to tranfer a sample from each well of 
the sample plates to a well in a 384-well DNA plate. 
This is the first stage of the NGS Genotyping pipeline.
"""

# can we validate barcode entry widgets as GUI focus leaves them, and colour their background accordingly?
# activate go button when ear-punch fields are valid

import os
import sys
import csv
import json
import tkinter as tk
#import tkinter.ttk as ttk
import tkinter.filedialog as fd
import tkinter.messagebox as mb  
import requests
from musterer import getPlate_app

#import xlrd
#import xlwt

import primercheck

assayMethod = "NGS" 

app = None # global!




def nimbus(tgtid, platesdata):
    """
    Process Musterer punch-plate data to produce a workfile for the BRF's Nimbus robot.
    See nimbus_custom() for custom data version.
    
    It needs the name of an output file and ear-punch data for each plate.
    """
    
    # use common code to read library files
    global assayMethod
    pl = primercheck.PrimerLookup()
    def getassays(winf):
        "combine old and new assay family names from Musterer"
        ngsassays = [n['assayName'] for n in winf.get('assays', []) if n['assayMethod']==assayMethod]
        gx = [n['assayName'] for n in winf.get('assays', []) if n['assayMethod']!=assayMethod and n['assayName'] not in pl.isna]
        mas = sorted(frozenset(x for xs in (ngsassays, gx) for x in xs))
        genold = (mt for ma in gx for mt in pl.mmdict.get(ma, [ma]))
        return mas, sorted(frozenset(x.split('_',1)[0] for xs in (genold, ngsassays) for x in xs))
    # basic validation
    tgtfn = 'Nimbus'+tgtid+'.csv'
    fnstg = 'Stage1-P'+tgtid+'.csv'
    # rn - row number in Nimbus picklist file - Nimbus needs this
    # wn - count the number of wells needed in the DNA plate
    # plist - list of primer families needed for the mice (samples) in the plate
    # unk - unknown assay family names; names that don't map to a primer family name
    rn, wn, plist, unk = 0, 0, [], []
    with open(tgtfn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
        dstpl  = csv.writer(dstfd1, dialect='unix')
        # Note: apparently, the Nimbus is extremely particular about column names.
        dstpl.writerow(['Sample no', 'Plate barcode', 'Well', 'Mouse barcode'])
        dststg = csv.writer(dstfd2, dialect='unix')
        stghdr = ['mouseID', 'EPplate', 'EPwell', 'mouseBarcode', 'strainName', 'mouseAssays', 'assayFamilies',
                 'sex', 'parentLitter', 'parentSequence', 'parentStrain', 'parentDate']
        dststg.writerow(stghdr)
        for no, pbc, shx in platesdata: # plate no, barcode, platedata (from JSON format)
            # print('Plate', pbc, 'data:', ', '.join(k+': '+str(v) for k,v in shx.items() if k!=' ))
            # Nope! assert str(shx['barcode']).endswith(pbc)
            if not shx:
                continue
            assert 'wells' in shx
            assert all(k in w for w in shx['wells'] for k in ('mouse', 'wellLocation'))
            wells = dict((w['wellLocation'], w['mouse']) for w in shx['wells'])
            # print("Sheet", shx.name) # use status bar
            for c in sorted(frozenset(int(w[1:]) for w in wells.keys())):
                # rows are two blocks of 4 rows in this order - doesn't actually matter
                for r in 'ACEGBDFH':
                    nwix = r+str(c) # well index
                    nwixm = r+str(c).zfill(2) # Musterer well ID
                    # nimbus needs non-empty Mouse Barcode so use '0' for 'empty'
                    mbc = str(wells[nwixm]["mouseBarcode"]) if nwixm in wells else '0'
                    rn += 1
                    rd = [rn, pbc, nwix, mbc]
                    dstpl.writerow(rd)
                    if nwixm in wells:
                        winf = wells[nwixm]
                        mas, mafs = getassays(wells[nwixm])
                        unk += [maf for maf in mafs if maf not in pl.pfdict]
                        if mafs:
                            xtras = {"mouseAssays": ';'.join(mas), "assayFamilies": ';'.join(mafs)}
                            xtra = [xtras[k] if k in xtras else winf.get(k, '') for k in stghdr[4:]]
                            dststg.writerow(rd+xtra)
                            wn += 1
                            plist += mafs
    # return values for use in cgi-nimbus2.py code
    return wn, sorted(frozenset(plist)), sorted(frozenset(unk))


def nimbus_custom(tgtid, platesdata, custom_assay_fn=None, custom_primer_fn=None):
    """
    Process custom sample-plate data, to produce a workfile for the BRF's Nimbus robot.
    
    It needs the name of an output file and sample data for each plate.
    """
    # use common code to read library files
    global assayMethod
    pl = primercheck.PrimerLookup(fnmm=custom_assay_fn, fnpp=custom_primer_fn)
    def getassays(winf):
        "combine old and new assay family names from Musterer"
        ngsassays = [n['assayName'] for n in winf.get('assays', []) if n['assayMethod']==assayMethod]
        gx = [n['assayName'] for n in winf.get('assays', []) if n['assayMethod']!=assayMethod and n['assayName'] not in pl.isna]
        mas = sorted(frozenset(x for xs in (ngsassays, gx) for x in xs))
        genold = (mt for ma in gx for mt in pl.mmdict.get(ma, [ma]))
        return mas, sorted(frozenset(x.split('_',1)[0] for xs in (genold, ngsassays) for x in xs))
    # basic validation
    tgtfn = 'Nimbus'+tgtid+'.csv'
    fnstg = 'Stage1-P'+tgtid+'.csv'
    # rn - row number in Nimbus picklist file - Nimbus needs this
    # wn - count the number of wells needed in the DNA plate
    # plist - list of primer families needed for the samples in the plate
    # unk - unknown assay family names; names that don't map to a primer family name
    rn, wn, plist, unk = 0, 0, [], []
    with open(tgtfn, "wt", newline='') as dstfd1, open(fnstg, "wt", newline='') as dstfd2:
        dstpl  = csv.writer(dstfd1, dialect='unix')
        # Note: apparently, the Nimbus is extremely particular about column names.
        dstpl.writerow(['Sample no', 'Plate barcode', 'Well', 'Sample barcode'])
        dststg = csv.writer(dstfd2, dialect='unix')
        stghdr = ['SampleNo', 'EPplate', 'EPwell', 'sampleBarcode', 'assays', 'assayFamilies']
        dststg.writerow(stghdr)
        for no, pbc, shx in platesdata: # plate no, barcode, platedata (from JSON format)
            #print('Plate', pbc, 'data:', ', '.join(k+': '+str(v) for k,v in shx.items() if k!=''), file=sys.stderr)
            # Nope! assert str(shx['barcode']).endswith(pbc)
            if not shx:
                print('No plate data found for plate', no, pbc, file=sys.stderr)
                continue
            assert 'wells' in shx
            assert all(k in w for w in shx['wells'] for k in ('organism', 'wellLocation'))
            wells = dict((w['wellLocation'], w['organism']) for w in shx['wells'])
            #print(wells.keys(), file=sys.stderr)
            # print("Sheet", shx.name) # use status bar
            for c in sorted(frozenset(int(w[1:]) for w in wells.keys())):
                # rows are two blocks of 4 rows in this order - doesn't actually matter
                for r in 'ACEGBDFH':
                    nwix = r+str(c) # well index
                    #TODO remove this restriction of wells looking like A01, etc
                    nwixm = r+str(c).zfill(2) # Sample plate well ID 
                    # nimbus needs non-empty sample barcode so use '0' for 'empty'
                    mbc = str(wells[nwixm]["sampleBarcode"]) if nwixm in wells else '0'
                    rn += 1
                    rd = [rn, pbc, nwix, mbc]
                    dstpl.writerow(rd)
                    if nwixm in wells:
                        winf = wells[nwixm]
                        mas, mafs = getassays(wells[nwixm])
                        unk += [maf for maf in mafs if maf not in pl.pfdict]
                        if mafs:
                            xtras = {"assays": ';'.join(mas), "assayFamilies": ';'.join(mafs)}
                            xtra = [xtras[k] if k in xtras else winf.get(k, '') for k in stghdr[4:]]
                            dststg.writerow(rd+xtra)
                            wn += 1
                            plist += mafs
    # return values for use in cgi-nimbus2.py code
    return wn, sorted(frozenset(plist)), sorted(frozenset(unk))


class NGSapp(tk.Frame):
    
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.master = master
        
        # work directory
        wdr=2
        l0 = tk.Label(self, text='Directory:')
        l0.grid(row=wdr, column=0, padx=2, pady=10, sticky='E')
        #self.wdv = tk.StringVar(self)
        self.wd = tk.Entry(self, width=60)
        self.wd.insert('end', os.getcwd())
        self.wd.grid(row=wdr, column=1, columnspan=3, pady=20)
        wb = tk.Button(self, text="Choose", command=self.getdir, width=8)
        wb.grid(row=wdr, column=4, padx=2, pady=5, sticky="W")
        
        def nextEnt(evt):
            "behave a bit like Tab"
            ent = self.focus_get()
            if ent.get():
                # move on if not empty
                ent.tk_focusNext().focus_set()
            return
        
        # punch plates
        ppr = 3
        ep_ent = []
        for i in range(4):
            lab = tk.Label(self, text="Ear-punch Plate No. {}:".format(i+1))
            lab.grid(row=ppr+i, column=0, padx=(5, 0), sticky=tk.E)
            ent = tk.Entry(self, relief="sunken", width=6)
            ent.grid(row=ppr+i, column=1, padx=(0, 5), sticky=tk.W)
            ent.bind('<Return>', nextEnt)
            ep_ent.append(ent)
        self.eps = ep_ent
            
        # sample/DNA plate
        dpr = 8
        l2 = tk.Label(self, text="barcode for 384-well DNA/sample plate:")
        l2.grid(row=dpr, column=0, columnspan=2, pady=(30, 0))
        #self.dev = tk.StringVar(self)
        self.de = tk.Entry(self, width=8)
        self.de.bind('<Return>', nextEnt)
        self.de.grid(row=dpr+1, column=0, columnspan=2, pady=5)
        
        # Go & Reset buttons
        qb = tk.Button(self, text="Go", command=self.gobutton, width=8)
        qb.bind('<Return>', self.gobutton)
        qb.grid(row=dpr+3, column=0, columnspan=2, pady=10)
        qc = tk.Button(self, text="Reset", command=self.reset, width=8)
        qc.grid(row = dpr+3, column=2, pady=10)
        return
    
    def reset(self):
        for ent in self.eps+[self.de]:
            ent.delete("0", tk.END)
        self.master.clearstatus()
        return
        
    def getdir(self):
        "set working directory for the application"
        dn = fd.askdirectory(title="Choose work directory.")
        os.chdir(dn)
        # maybe search for existing plate files and populate plateID textbox?
        self.wd.delete("0", tk.END)
        self.wd.insert(tk.END, dn)
        # print("changed to directory", dn)
        return
    
    def plateids(self):
        "list of plate IDs (barcodes) as strings."
        # should check formats and for duplicates ... drop empties ...
        return [x for x in ((i, ent.get().strip()) for i, ent in enumerate(self.eps)) if x[1]]
        # return [x for x in self.xx.get("1.0", tk.END).strip().split() if x]      
    
    def gobutton(self, evt=None):
        """
        This does most of the processing.
        Check that the required fields are completed.
        Fetch Plate files from Musterer database.
        Process Plate files to produce an input file for the Nimbus robot.
        """
        global app
        wdn = self.wd.get()
        assert wdn==os.getcwd()

        # fetch the plate data
        pnos = self.plateids()
        if pnos and not len(pnos)<=4: 
            mb.showwarning(title="Warn!", message="Please provide 1-4 Plate IDs.")
            return
        pxs = [(n, pno, getPlate_app(pno)) for n, pno in pnos]
        if not pxs:
            mb.showerror(message="No valid ear-punch plate IDs.")
            return
        if not all(x[2] for x in pxs):
            app.setstatus("invalid ear-punch plates: "+' '.join(x[1] for x in pxs if not x[2]))
            return
        
        tgtid = self.de.get()
        if not tgtid:
            app.setstatus("Missing target plate ID!")
            # change to ask about overwrite ...
            mb.showwarning(title="Warn!", message="Please provide Target Plate ID.")
            return
        app.setstatus("Target Plate ID: "+tgtid)
        res = nimbus(tgtid, pxs)
        app.setstatus("Created Plate", tgtid, '- click exit when done.')
        return res
           
import base64
        
class MyApp(tk.Frame):
    # embed the logo image in this file so no separate file is needed
    logo = base64.b64decode('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')
    
    def __init__(self, master=None):
        global app
        app = self
        tk.Frame.__init__(self, master, bg="lightgreen")
        #self.master = master
        self.pack()
        
        # logo
        # ix = tk.PhotoImage(file="ABCLogo.png")
        ix = tk.PhotoImage(data=self.logo)
        zzz=tk.Label(self, image=ix)
        zzz.image=ix # this is required or the image is garbage collected
        zzz.grid(row=0, columnspan=2, column=0, pady=20)
        
        ngsapp=NGSapp(self)
        ngsapp.grid(row=1, column=0, columnspan=2, padx=5, pady=(5, 0), ipadx=10)        

        # extra
        epr = 5
        cb = tk.Button(self, text="Exit", command=self.cancelbutton, width=8)
        cb.grid(row=epr, column=0, padx=5, pady=20)        

    
        brandfont=("Times", 10, "italic")
#        brandstyle = ttk.Style()
#        brandstyle.configure("Brand.TLabel", font=brandfont, foreground="#22D")
#        brandstyle.configure("Brand.TLabel", background="#FF8")
#        brand=ttk.Label(self, text="Created by ANU Bioinformatics Consultancy.", style="Brand.TLabel")
        brand=tk.Label(self, text="Created by ANU Bioinformatics Consultancy.", font=brandfont, bg="#FF8", fg="#22D")
        #brand.pack(side=tk.RIGHT, anchor=tk.E, padx=5, pady=5)
        brand.grid(row=epr+1, column=1)
        
        statfont=("Times", 12)
        self.status = tk.StringVar()
#        brandstyle.configure("Stat.TLabel", font=statfont)
#        statlab = ttk.Label(self, relief=tk.RIDGE, anchor=tk.W, width=60, textvar=self.status, style="Stat.TLabel")
        statlab = tk.Label(self, relief=tk.RIDGE, anchor=tk.W, width=50, 
                           textvar=self.status, font=statfont, bg='#EFE')
        #statlab.pack(side=tk.LEFT, anchor=tk.W, padx=5, pady=2)
        statlab.grid(row=epr+1, column=0)
        self.setstatus("Initialising ...")
        self.after(500, self.clearstatus)
        return
    
    # should have background colours for warnings, info ...
    # use these instead of print() or mb.messagebox()
    def setstatus(self, *args):
        "display status in the status field (print is not visible)"
        self.status.set("  "+' '.join(map(str, args))) #extra spaces since ipadx option doesn't seem to work
        self.update()
        return
    
    def clearstatus(self):
        "set status to empty"
        self.status.set('')
        self.update()
        return

    def cancelbutton(self):
        "exit program ..."
        self.setstatus("Exit ...")
        self.after(300, self.quit) # quit after a slight delay
        return

if __name__=="__main__":
    root = tk.Tk()
    root.title("Nimbus setup")
    #root.geometry("800x300")
    iconfiles = [fn for fn in (os.path.join(d, 'abc.ico') for d in ('.', '../library', os.path.dirname(sys.argv[0]))) if os.path.isfile(fn)]
    root.wm_iconbitmap(iconfiles[0] if iconfiles else '')
    MyApp(root)
    root.mainloop()
