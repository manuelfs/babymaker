#!/usr/bin/env python

###### Script to send jobs to merge ntuples
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

# Setting folders
folder    = '/net/cms27/cms27r0/babymaker/2016_04_29/mc/'
skims = ['abcd', 'baseline', 'sys_abcd', '1lht500met200']
skims = ['abcd']



os.system("JobSetup.csh")
for skim in skims:
  # Checking if input folder exists and creating output folder
  infolder = folder+"/skim_"+skim+"/"
  if not os.path.exists(infolder):
    print "Skim "+infolder+"does not exist. Skipping"
    continue
  outfolder = folder+"/merged_"+skim+"/"
  if not os.path.exists(outfolder):
    os.system("mkdir -p "+outfolder)

  ## Finding tags for each dataset
  infiles = set()
  for file in glob.glob(infolder+'/*.root'):
    tag = file.split("RunII")[0]
    if ("TTJets_Tune" not in tag) and ("DYJetsToLL_M-50_Tune" not in tag): tag = tag.split("Tune")[0]
    tag = tag.split("13TeV")[0]
    tag = tag.split("baby_")[1]
    if "DYJetsToLL_M-50_Tune" in tag: infiles.add("_"+tag)

  # Sending jobs for each tag
  runfolder = outfolder+"run/" 
  if not os.path.exists(runfolder):
    os.system("mkdir -p "+runfolder)
  ijob = 0
  for ifile in infiles:
    ijob += 1
    exename = runfolder+"/merge_ntuples_"+str(ijob)+".sh"
    fexe = open(exename,"w")
    os.system("chmod u+x "+exename)
    fexe.write("#!/bin/bash\n\n")
    fexe.write("./run/merge_ntuples.exe "+infolder+' '+outfolder+' '+ifile+' '+skim+'\n')
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh "+exename
    #print cmd
    #os.system('ls '+folder+'*'+ifile+'*')
    os.system(cmd)
	
  print "\nSubmitted "+str(len(infiles))+" merging jobs. Output goes to "+outfolder+"\n"

sys.exit(0)
