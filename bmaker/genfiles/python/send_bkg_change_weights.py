#!/usr/bin/env python

###### Script to send 40 jobs to the batch changing weights to the signal scan
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

infiles = [
['/net/cms2/cms2r0/babymaker/babies/2015_11_28/to_normalize/bkg/', '_TTJets_SingleLeptFromT_'],
['/net/cms2/cms2r0/babymaker/babies/2015_11_28/to_normalize/bkg/', '_TTJets_SingleLeptFromTbar_'],
['/net/cms2/cms2r0/babymaker/babies/2015_11_28/to_normalize/bkg/', '_TTJets_DiLept_'],
['/mnt/hadoop/cms/store/user/ana/TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1__MINIAODSIM/151127_101155/0000', 'TTJets_HT-600to800_'],
['/mnt/hadoop/cms/store/user/ana/TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_101126/0000', 'TTJets_HT-2500toInf_'],
['/mnt/hadoop/cms/store/user/ana/TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_101112/0000', 'TTJets_HT-1200to2500_'],
['/mnt/hadoop/cms/store/user/ana/TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_101208/0000', 'TTJets_HT-800to1200_'],
['/mnt/hadoop/cms/store/user/ana/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_101317/0000', 'TTJets_TuneCU'],
['/mnt/hadoop/cms/store/user/ana/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_104351/0000', 'ST_s-channel_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_105423/0000', 'QCD_HT1000to1500_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_105518/0000', 'QCD_HT1500to2000_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_105632/0000', 'QCD_HT2000toInf_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_105028/0000', 'QCD_HT200to300_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_105133/0000', 'QCD_HT300to500_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_105233/0000', 'QCD_HT500to700_'],
['/mnt/hadoop/cms/store/user/manuelf/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_105329/0000', 'QCD_HT700to1000_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_112743/0000', 'WZTo2L2Q_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WWToLNuQQ_13TeV-powheg/crab_WWToLNuQQ_13TeV-powheg__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112836/0000', 'WWToLNuQQ_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WWTo2L2Nu_13TeV-powheg/crab_WWTo2L2Nu_13TeV-powheg__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112819/0000', 'WWTo2L2Nu_'],
['/mnt/hadoop/cms/store/user/jaehyeok/ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8/crab_ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1__MINIAODSIM/151127_112858/0000', 'ttHJetTobb_'],
['/mnt/hadoop/cms/store/user/jaehyeok/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_112133/0000', 'DYJetsToLL_M-50_HT-100to200_'],
['/mnt/hadoop/cms/store/user/jaehyeok/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_112152/0000', 'DYJetsToLL_M-50_HT-200to400_'],
['/mnt/hadoop/cms/store/user/jaehyeok/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_112217/0000', 'DYJetsToLL_M-50_HT-400to600_'],
['/mnt/hadoop/cms/store/user/jaehyeok/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2__MINIAODSIM/151127_112314/0000', 'DYJetsToLL_M-50_HT-600toInf_'],
['/mnt/hadoop/cms/store/user/jaehyeok/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112115/0000', 'DYJetsToLL_M-50_Tu'],
['/mnt/hadoop/cms/store/user/jaehyeok/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112447/0000', 'ST_t-channel_antitop_'],
['/mnt/hadoop/cms/store/user/jaehyeok/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112505/0000', 'ST_t-channel_top_4f_'],
['/mnt/hadoop/cms/store/user/jaehyeok/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112358/0000', 'ST_tW_antitop_5f_i'],
['/mnt/hadoop/cms/store/user/jaehyeok/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112423/0000', 'ST_tW_top_5f_'],
['/mnt/hadoop/cms/store/user/jaehyeok/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112522/0000', 'TTGJets_'],
['/mnt/hadoop/cms/store/user/jaehyeok/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1__MINIAODSIM/151127_112558/0000', 'TTTT_'],
['/mnt/hadoop/cms/store/user/jaehyeok/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112615/0000', 'TTWJetsToLNu_'],
['/mnt/hadoop/cms/store/user/jaehyeok/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112633/0000', 'TTWJetsToQQ_'],
['/mnt/hadoop/cms/store/user/jaehyeok/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112650/0000', 'TTZToLLNuNu_'],
['/mnt/hadoop/cms/store/user/jaehyeok/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112709/0000', 'TTZToQQ_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_111959/0000', 'WJetsToLNu_HT-100To200_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112017/0000', 'WJetsToLNu_HT-200To400_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3__MINIAODSIM/151127_112036/0000', 'WJetsToLNu_HT-400To600_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112054/0000', 'WJetsToLNu_HT-600ToInf_'],
['/mnt/hadoop/cms/store/user/jaehyeok/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM/151127_112801/0000', 'WZTo3LNu_'],
]

# Setting folders
outfolder = "out/" 
runfolder = outfolder+"run/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

os.system("JobSetup.csh")
ijob = 0
for ifile in infiles:
  ijob += 1
  if ijob == 1: continue
  exename = runfolder+"/change_weights_"+str(ijob)+".sh"
  fexe = open(exename,"w")
  os.system("chmod u+x "+exename)
  fexe.write("#!/bin/bash\n\n")
  fexe.write("./run/change_weights.exe "+ifile[0]+' "*'+ifile[1]+'*" '+outfolder+'\n')
  fexe.close()
  cmd = "JobSubmit.csh ./run/wrapper.sh ./"+exename
  #print cmd
  os.system(cmd)

print "\nSubmitted "+str(ifile)+" files in "+str(ijob)+" jobs. Output goes to "+outfolder+"\n"
sys.exit(0)
