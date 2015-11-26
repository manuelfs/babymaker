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
for dset in datasets:
    print "Let's take a look at " + dset
    #get list of all original batch jobs (no skim)
    batches = glob.glob(outdir+"/"+dset+"_mf?_batch?.root")
    for b in batches:
        # for each original job, look for resubmits
        if len(glob.glob(b.split(".root")[0]+"_rs*.root")) > 0:
            print "delete "+b
            #get names of skims of original output,too
            to_remove = glob.glob(b.split(".root")[0]+"_mgluino*.root")
            os.rename(b,outdir+"/to_rm/"+b.split("/").pop()) #stash original
            for rm in to_remove:
                os.rename(rm,outdir+"/to_rm/"+rm.split("/").pop()) #stash skims of original
    
