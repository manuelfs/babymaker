#!/usr/bin/env python

#=======================================================
# From e.g. cms29 
  
# export SCRAM_ARCH=slc6_amd64_gcc491
# source /cvmfs/cms.cern.ch/cmsset_default.sh
  
# # checkout CMSSW 
# cmsrel CMSSW_7_4_6_patch6
# cd CMSSW_7_4_6_patch6/src/
# cmsenv
  
# # checkout and compile CfANtupler 
# git clone git@github.com:manuelfs/babymaker
# cd babymaker
# ./compile.sh
# cd ..
  
# # now you will have to edit generate_crab_cfg.py in the scripts directory 
# # only change is to put the datasets you want in the list at the beginning 
  
# # setup crab 
#  source /cvmfs/cms.cern.ch/crab3/crab.sh 
#  voms-proxy-init --voms cms --valid 168:00
  
# # submit crab job 
# python babymaker/scripts/generate_crab_cfg.py
#========================================================

import das_client as das
import json
import os
import sys

def getNumberOfEvents(dataset):
  query = "file dataset=" + dataset + " | sum(file.nevents)"

  data = das.get_data(query)
  if isinstance(data, basestring):
    dasjson = json.loads(data)
  else:
    dasjson = data
  status  = dasjson.get('status')
  if  status == 'ok':
    data = dasjson.get('data')
    sumevents=0
    for idata in data:
      sumevents+=idata.get('result').get('value')
    return sumevents

doAdam = False 
doAna = False
doRyan = False
doMissingDatasets=False;

doFastSimScans = False

datasets = [[]]

if doMissingDatasets:
    datasets.append(["/SingleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD"])
    
if doRyan:
    datasets.append(["/MET/Run2016B-03Feb2017_ver2-v2/MINIAOD"])
    datasets.append(["/MET/Run2016C-03Feb2017-v1/MINIAOD"])
    datasets.append(["/MET/Run2016D-03Feb2017-v1/MINIAOD"])
    
    datasets.append(["/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD"])
    datasets.append(["/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD"])
    datasets.append(["/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD"])

    datasets.append(["/SingleMuon/Run2016C-03Feb2017-v1/MINIAOD"])
    datasets.append(["/SingleMuon/Run2016D-03Feb2017-v1/MINIAOD"])

    datasets.append(["/JetHT/Run2016B-03Feb2017_ver2-v2/MINIAOD"])
    datasets.append(["/JetHT/Run2016C-03Feb2017-v1/MINIAOD"])
    datasets.append(["/JetHT/Run2016D-03Feb2017-v1/MINIAOD"])

if doAdam:
 
    datasets.append(["/MET/Run2016E-03Feb2017-v1/MINIAOD"])
    datasets.append(["/MET/Run2016F-03Feb2017-v1/MINIAOD"])

    datasets.append(["/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD"])
    datasets.append(["/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD"])
    
    datasets.append(["/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD"])
    datasets.append(["/SingleMuon/Run2016F-03Feb2017-v1/MINIAOD"])
    
    datasets.append(["/JetHT/Run2016E-03Feb2017-v1/MINIAOD"])
    datasets.append(["/JetHT/Run2016F-03Feb2017-v1/MINIAOD"])

if doAna:
    datasets.append(["/MET/Run2016G-03Feb2017-v1/MINIAOD"])
    datasets.append(["/MET/Run2016H-03Feb2017_ver2-v1/MINIAOD"])
    datasets.append(["/MET/Run2016H-03Feb2017_ver3-v1/MINIAOD"])
    
    datasets.append(["/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD"])
    datasets.append(["/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD"])
    datasets.append(["/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD"])
 
    datasets.append(["/SingleMuon/Run2016G-03Feb2017-v1/MINIAOD"])
    datasets.append(["/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD"])
    datasets.append(["/SingleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD"])

    datasets.append(["/JetHT/Run2016G-03Feb2017-v1/MINIAOD"])
    datasets.append(["/JetHT/Run2016H-03Feb2017_ver2-v1/MINIAOD"])
    datasets.append(["/JetHT/Run2016H-03Feb2017_ver3-v1/MINIAOD"])



if doFastSimScans:
  datasets.append(["/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM"])
  datasets.append(["/SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM"])

for ilist in datasets:
  nevents=0
  print ""
  for ids in ilist:
    #nevents += getNumberOfEvents(ids) # Not working currently
    nevents += 1000000 ## Only needed to get w_lumi right on the first go, will fix at renormalization
    #print "running nevents is "+str(nevents)

  for ids in ilist:
    cmssw_base = os.getenv("CMSSW_BASE")
    datasetID = ids.replace('/','',1).replace('/', '_', 1)
    datasetID = datasetID[0:datasetID.find('/')]
    inputfile = cmssw_base + "/src/babymaker/bmaker/python/crab_cfg_template.py"
    outputfile = "crab_cfg_" + datasetID + ".py"

    s = open(inputfile).read()
    s = s.replace('DATASETNAME', ids)
    s = s.replace('NEVENTS', str(nevents))
    if ("Run2016F" in datasetID): # split RunF in two parts due to JECs
      outputfile = "crab_cfg_" + datasetID + "_part1.py"
      with open(outputfile, 'w') as f: f.write(s.replace('RUN_RANGE', 'Run2016F1'))
      os.system("crab submit -c " + outputfile)

      outputfile = "crab_cfg_" + datasetID + "_part2.py"
      with open(outputfile, 'w') as f: f.write(s.replace('RUN_RANGE', 'Run2016F2'))      
      os.system("crab submit -c " + outputfile)
    else:
      with open(outputfile, 'w') as f: f.write(s)
      os.system("crab submit -c " + outputfile)

    
    print "Submitted ",ids
    # sys.exit(0)
