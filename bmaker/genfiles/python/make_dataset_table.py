#!/usr/bin/env python

from ROOT import TChain, TFile, TTree, TH1D
import os, sys, subprocess
import pprint
import glob
import json
import string
import time
import math

mcdir = '/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed'
dadir = '/cms2r0/babymaker/babies/2016_06_14/data/unskimmed'

mcset = set([i.split('_RunIISpring')[0].split('fullbaby_')[-1] for i in glob.glob(mcdir+"/*.root")])
daset = set([i.strip(i.split('_')[-1]).rstrip('_').split('fullbaby_')[-1] for i in glob.glob(dadir+"/*.root")])

def table_header(cols=3):
    header = "\\documentclass{article}\n" 
    header += "\\usepackage{amsmath,graphicx,rotating}\n" 
    header += "\\usepackage[landscape]{geometry}\n" 
    header += "\\thispagestyle{empty}\n" 
    header += "\\begin{document}\n" 
    header += "\\begin{table}\n" 
    header += "\\centering\n" 
    header += "\\resizebox{\\textwidth}{!}{\n" 
    header += "\n\\begin{tabular}[tbp!]{ l " 
    for i in range(cols): header += "r" 
    header += "}\\hline\\hline\n"
    return header

def table_footer():
    footer = "\\hline\\hline\n\\end{tabular}\n\n"
    footer += "}\n"
    footer += "\\end{table}\n"
    footer += "\\end{document}\n"
    return footer


darows = [] 
for ds in sorted(daset):
    darows.append('\\multicolumn{1}{l}{'+ds.replace('_','\\_')+'}\\\\\n')

mcrows = []
for ds in sorted(mcset):
    c = TChain("tree")
    c.Add(mcdir+"/*"+ds+"*.root")
    ntot = c.GetEntries()
    f = TFile(glob.glob(mcdir+"/*"+ds+"*.root")[0],"READ")
    g = f.Get("treeglobal")
    xsec = -1
    for event in g:
        xsec = event.xsec
        break
    lumi = ntot/xsec/1000.
    print ds, ntot, xsec, lumi
    mcrows.append(' &'.join([ds.replace('_','\\_'),'{:,}'.format(ntot),'{:,.2f}'.format(lumi)])+' \\\\\n')

texfile = "dataset_table.tex"
tex = open(texfile,"w")
tex.write(table_header(3))
for i in darows: tex.write(i)
tex.write('\\hline\\\\\n')
for i in mcrows: tex.write(i)
tex.write(table_footer())
tex.close()

os.system("cat "+texfile)
