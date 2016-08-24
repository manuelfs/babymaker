#!/usr/bin/env python

###### Script to send jobs to merge ntuples
import os, sys, subprocess
import pprint
import glob
import json
import string
import time
import pprint

# Setting folders
datasets = 'txt/alldata.txt'
infolder  = '/net/cms2/cms2r0/babymaker/babies/2016_08_10/data/unprocessed/'
outfolder = '/net/cms2/cms2r0/babymaker/babies/2016_08_10/data/unskimmed/'
jsonfile = '../../data/json/golden_Cert_271036-277933_13TeV_PromptReco_Collisions16_NoL1T.json'

runs_file = 1 # Number of runs in each ntuple


runs = []
with open(jsonfile) as jfile:
  for line in jfile:
    for word in line.split():
      if '"' in word: 
        word = word.split('"')[1]
        runs.append(word)

# Dividing runs into sets of "runs_file" elements
runs = [runs[i:i+runs_file] for i in xrange(0, len(runs), runs_file)]

# Sending jobs for each set of runs
os.system("JobSetup.csh")
for run in runs:
  cmd = "JobSubmit.csh ./run/wrapper.sh ./run/combine_datasets.exe -i "+infolder+" -o "+outfolder+" -f "+datasets+" -b "+run[0]+" -e "+run[-1]
  # print cmd
  os.system(cmd)

sys.exit(0)
