#!/usr/bin/env python

# in case we decide to use David's batch:
import os, sys, pprint, glob

ntup_date = '2016_02_09'
model = "T2tt"

if model == "T1tttt":
  outname = "baby_SMS-T1tttt_MASS_TAG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9.root"
if model == "T5tttt":
  outname = "baby_SMS-T5ttttDM175_MASS_TAG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1.root"
if model == "T6ttWW":
  outname = "baby_SMS-T6ttWW-mLSP50_MASS_TAG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9-v1.root"
if model == "T2tt":
  outname = "baby_SMS-T2tt_MASS_TAG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9-v1.root"

outfolder = "/net/cms2/cms2r0/babymaker/babies/"+ntup_date+"/to_normalize/"+model 
if not os.path.exists(outfolder):
  os.system("mkdir -p "+outfolder)

# CHANGE ACCORDING TO NEEDS: Only sending a specific mass point for a given date here
dirs = [i for i in glob.glob("/mnt/hadoop/cms/store/user/manuelf/*/*/*/*/") if model in i]

os.system("JobSetup.csh")
for infolder in dirs:
  cmd = "JobSubmit.csh ./run/wrapper.sh ./run/skim_scan_onefile.exe %s %s" % (infolder, outfolder+"/"+outname)
  print cmd+"\n"
  os.system(cmd)


##### For condor submission

# import os, sys, subprocess
# import pprint
# import glob
# import json
# import string
# import time
# import argparse
# 
# cmsswdir = '/net/cms2/cms2r0/babymaker/CMSSW_7_4_14/src'
# script = '/net/cms2/cms2r0/manuelf/cmssw/babies/CMSSW_7_4_14/src/babymaker/bmaker/genfiles/run/skim_scan_onefile.exe'
# 
# ntup_date = '2015_11_27'
# outfolder = "/net/cms2/cms2r0/babymaker/babies/"+ntup_date+"/sms/split_sms3" 
# baby_name_template = "baby_SMS-T1tttt_MASS_TAG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9.root"
# 
# if not os.path.exists(outfolder):
#   os.mkdir(outfolder)
# 
# tmpdir = os.path.join(outfolder, 'tmp')
# if not os.path.exists(tmpdir):
#   os.makedirs(tmpdir)
# print "INFO: Babies will be written to: ", outfolder
# print "INFO: Logs will be written to:   ", tmpdir
# 
# #input datasets
# datasets = [i for i in os.listdir("/mnt/hadoop/cms/store/user/manuelf/") if "SMS" in i]
# dirs = [i for i in glob.glob("/mnt/hadoop/cms/store/user/manuelf/*/*/*/*/") if "SMS" in i]
# 
# # check if there is more then one directory per dataset
# if len(datasets)!=len(dirs): 
#   print("We have a little problem... Found dataset with more than 1 sub-folder and TChain does not take e.g. dir/*/*root")
#   sys.exit(1)
# 
# # Create condor submission cmd file
# cmdfile = tmpdir+'/submit_all.cmd'
# fcmd = open(cmdfile,"w")
# 
# total_jobs = 0
# for infolder in dirs:
#   ids = infolder.split("/")[7]
#   # set up script arguments
#   condor_args = []
#   condor_args.append(infolder)
#   condor_args.append(outfolder+"/"+baby_name_template)
# 
#   # Create executable that will be transfered to the work node by condor
#   exefile = tmpdir+'/'+ids+".sh"
#   fexe = open(exefile,"w")
#   fexe.write("#! /bin/bash\n")
#   fexe.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
#   fexe.write("cd "+cmsswdir+"\n")
#   fexe.write("eval `scramv1 runtime -sh`\n")
#   fexe.write(script+" \\\n"+" \\\n".join(condor_args)+"\n")
#   fexe.close()
#   os.system("chmod u+x "+exefile)
# 
#   fcmd.write("Executable   = "+exefile+"\n")
#   fcmd.write("Universe     = vanilla\n")
#   fcmd.write("Log          = "+tmpdir+ "/"+ids+".log\n")
#   fcmd.write("output       = "+tmpdir+"/"+ids+".out\n")
#   fcmd.write("error        = "+tmpdir+"/"+ids+".err\n")
#   fcmd.write("Notification = never\n")
#   fcmd.write("Queue\n")
#   fcmd.write("\n")
#   total_jobs = total_jobs + 1
# 
# fcmd.close()
# 
# # Submit condor jobs
# #os.system("ssh cms25.physics.ucsb.edu condor_submit "+cmdfile)
# #print "Submitted ", total_jobs, "jobs"

