#!/usr/bin/env python

###### Script to send jobs renormalizing all samples in infolder
import os, sys, subprocess
import string
from utilities import *

infolder  = '/net/cms29/cms29r0/babymaker/babies/2016_08_10/mc/unprocessed/'
outfolder = '/net/cms29/cms29r0/babymaker/babies/2016_08_10/mc/unskimmed/'

infolder  = '/net/cms29/cms29r0/babymaker/babies/2016_08_10/T1tttt/unprocessed/'
outfolder = '/net/cms29/cms29r0/babymaker/babies/2016_08_10/T1tttt/unskimmed/'


## Finding tags for each dataset
sortedfiles = findBaseSampleNames(infolder)

# Setting folders
runfolder = outfolder+"run/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

os.system("JobSetup.csh")
ijob = 0
for ifile in sortedfiles:
  ijob += 1
  exename = runfolder+"/change_weights_"+str(ijob)+".sh"
  fexe = open(exename,"w")
  os.system("chmod u+x "+exename)
  fexe.write("#!/bin/bash\n\n")
  execmd = "\n./run/change_weights.exe "+infolder+' "*'+ifile+'*" '+outfolder
  fexe.write(execmd)
  fexe.close()
  cmd = "JobSubmit.csh ./run/wrapper.sh "+exename
  print execmd
  os.system(cmd)
  #print '\n\n'
  #print "./run/change_weights.exe "+infolder+' "*'+ifile+'*" '+outfolder+'\n'
  # os.system('ls '+infolder+'*'+ifile+'*')

print "\nSubmitted "+str(ijob)+" jobs. Output goes to "+outfolder+"\n"
sys.exit(0)
