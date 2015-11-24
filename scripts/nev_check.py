#!/usr/bin/env python
import os, sys, re
import glob
import string
from pprint import pprint
import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t","--timestamp")
args = parser.parse_args()

if (args.timestamp):
  timestamp = args.timestamp
else:
  sys.exit("Please provide a timestamp either as, e.g. 151019_011440")


bdir = os.getcwd()

if ('babymaker' not in bdir.split("/").pop()):
  sys.exit("Execute from babymaker directory")

outdir = os.path.join(bdir,'out',timestamp)
print outdir
if not os.path.exists(outdir):
  sys.exit("Can't find out directory %s" %outdir)

#one line to get list of unique dataset names in outdir  
datasets = set([x.split("_mf")[0].split("/").pop() for x in glob.glob(outdir+"/*.root")])

mismatch = 0
match = 0
totalsuccess = 0
totaltotal = 0 #awesomename
for dset in datasets:
    print "Let's take a look at " + dset
    ch = ROOT.TChain("tree")
    chglob = ROOT.TChain("treeglobal")
    if("T1tttt" in dset): # for now just checking mass skims
        ch.Add(outdir+"/"+dset+"_mf*mlsp*.root")
        chglob.Add(outdir+"/"+dset+"_mf*mlsp*.root")
    else:
        ch.Add(outdir+"/"+dset+"_mf*.root")
        chglob.Add(outdir+"/"+dset+"_mf*.root")
        
    #make sure all root files agree how many events there should be
    nev_list = []
    for entry in chglob:
        nev_list.append(entry.nev_sample)
    if len(set(nev_list)) > 1:
        print "incompatible root files combined"
    
    nsuccess = ch.GetEntries()
    ntotal = nev_list[0]
    totaltotal += ntotal
    totalsuccess +=nsuccess

    if ntotal == nsuccess:
        print "COMPLETE: "+str(nsuccess)+" / "+str(ntotal)+" events made it. Great job."
        match+=1
    elif ntotal > nsuccess:
        print "NOT COMPLETE:"+str(nsuccess)+" / "+str(ntotal)+" events made it. Looks like you can't go play in the mountains after all."
        mismatch+=1

    elif ntotal < nsuccess:
        print "OVER COMPLETE:"+str(nsuccess)+" / "+str(ntotal)+" events made it. Looks like you can't go play in the mountains after all."
        mismatch+=1    

    print ""
    
print "Attempted "+str(totaltotal)+" events. Completed "+str(totalsuccess)
print "Completed "+str(match) +"/"+str(match+mismatch)+" datasets."
    
    
    
    
