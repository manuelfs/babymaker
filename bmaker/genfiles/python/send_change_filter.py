#!/usr/bin/env python

###### Script to send 50 jobs to the batch changing the CSC beam halo filter in data
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

# Setting folders
infolder  = "/cms2r0/babymaker/babies/2015_11_20/data/singlelep/combined/"
infolder  = "/cms2r0/babymaker/babies/2015_11_20/data/singlelep/combined/skim_1lht500met200/"
outfolder = "out_data/" 
runfolder = outfolder+"run/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

#input datasets
inputfiles = [i for i in os.listdir(infolder) if "root" in i]

os.system("JobSetup.csh")
njobs = 1
files_job = (len(inputfiles)+njobs-1)/njobs
ifile = 0
ijob = 0
for file in inputfiles:
  ifile += 1
  # Creating executable
  if ifile % files_job == 1:
    ijob += 1
    exename = runfolder+"/change_weights_"+str(ijob)+".sh"
    fexe = open(exename,"w")
    os.system("chmod u+x "+exename)
    fexe.write("#!/bin/bash\n\n")
  fexe.write("./run/change_filter.exe -i "+infolder+' -s "*'+file+'*" -o '+outfolder+'\n')
  if ifile % files_job == 0 or ifile == len(inputfiles): 
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh ./"+exename
    #print cmd
    #os.system(cmd)

print "\nSubmitted "+str(ifile)+" files in "+str(ijob)+" jobs. Output goes to "+outfolder+"\n"
sys.exit(0)
