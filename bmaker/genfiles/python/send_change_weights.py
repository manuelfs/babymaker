#!/usr/bin/env python

###### Script to send jobs to the batch changing weights to 20 signal files at a time
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

# Setting folders
infolder  = "/net/cms2/cms2r0/babymaker/babies/2015_11_27/sms/split_sms/"
outfolder = "outsms/" 
runfolder = outfolder+"run/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

#input datasets
inputfiles = [i for i in os.listdir(infolder) if "SMS" in i]

os.system("JobSetup.csh")
files_job = 20
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
  fexe.write("./run/change_weights.exe "+infolder+' "*'+file+'*" '+outfolder+'\n')
  if ifile % files_job == 0 or ifile == len(inputfiles): 
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh ./"+exename
    #print cmd
    os.system(cmd)
    #sys.exit(0)

print "\nSubmitted "+str(ifile)+" files in "+str(ijob)+" jobs. Output goes to "+outfolder+"\n"
sys.exit(0)
