#!/usr/bin/env python

###### Script that adds the era to data root files that have been split by run (and have the run number in the name)
import os, sys, subprocess
import pprint
import glob
import argparse

## Parsing input arguments
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-i", "--infolder", help="Folder to find root files in", 
                    default="/net/cms29/cms29r0/babymaker/babies/2016_11_21/data/unskimmed/")
parser.add_argument("-o", "--outfolder", help="Folder to write files to", 
                    default="/net/cms29/cms29r0/babymaker/babies/2016_11_21/data/eras/")
args = parser.parse_args()
args.outfolder = args.outfolder+"/"

eras = [["RunB", [272007, 275376]], ["RunC", [275657, 276283]], ["RunD", [276315, 276811]], ["RunE", [276831, 277420]], 
        ["RunF", [277772, 278808]], ["RunG", [278820, 280385]], ["RunH", [280919, 28038500]]] ## Added zeros to runH 

## Create output folder
if not os.path.exists(args.outfolder):
  print "\nCreating "+args.outfolder
  os.system("mkdir -p "+args.outfolder)

noera_runs = []
nfiles = 0
files = glob.glob(args.infolder+'/*.root')
for file in files:
  outfile = file.replace(args.infolder, args.outfolder)

  ## Parsing run from file name
  run = file.split('runs')[-1]
  run = int(run.split('.root')[0])

  ## Finding era for run
  found_era = False
  for era in eras:
    if run >= era[1][0] and run <= era[1][1]: 
      outfile = outfile.replace("baby_", "baby_"+era[0]+"_")
      found_era = True
      break
  if not found_era: noera_runs.append(run)

  ## Moving file
  cmd = "mv "+file+" "+outfile
  os.system(cmd)
  nfiles += 1


## Final printouts
if len(noera_runs)>0:
  print "\nNot found era for runs "
  print noera_runs
print "\nCopied and renamed "+str(nfiles)+" files into "+args.outfolder+"\n\n"

