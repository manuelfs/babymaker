#!/bin/env python
""" run skims for the 13 TeV RPV gluino analysis"""
import os
import time

HTCUT = 1200
BASEDIR = "/homes/cawest/links"
OUTDIR = "/net/cms29/cms29r0/cawest/skims/ht" + str(HTCUT)
CUT = "ht>" + str(HTCUT)
DATASETS = []
DATASETS.append(BASEDIR + "/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DATASETS.append(BASEDIR + "/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DATASETS.append(BASEDIR + "/JetHT_Run2015C_25ns-05Oct2015-v1")
DATASETS.append(BASEDIR + "/JetHT_Run2015D-05Oct2015-v1")
DATASETS.append(BASEDIR + "/JetHT_Run2015D-PromptReco-v4")
DATASETS.append(BASEDIR + "/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DATASETS.append(BASEDIR + "/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1")
DATASETS.append(BASEDIR + "/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DATASETS.append(BASEDIR + "/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1")
DATASETS.append(BASEDIR + "/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DATASETS.append(BASEDIR + "/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1")
DATASETS.append(BASEDIR + "/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1")
DATASETS.append(BASEDIR + "/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1")
DATASETS.append(BASEDIR + "/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1")
DATASETS.append(BASEDIR + "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1")
DATASETS.append(BASEDIR + "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1")
DATASETS.append(BASEDIR + "/TTJets_Mtt-1000toInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
DATASETS.append(BASEDIR + "/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
DATASETS.append(BASEDIR + "/TTJets_TuneCUETP8M1_13TeV-madgraphMLM")
DATASETS.append(BASEDIR + "/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8")
DATASETS.append(BASEDIR + "/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8")
DATASETS.append(BASEDIR + "/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8")
DATASETS.append(BASEDIR + "/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8")
DATASETS.append(BASEDIR + "/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8")
DATASETS.append(BASEDIR + "/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DATASETS.append(BASEDIR + "/ZJetsToQQ_HT600toInf_13TeV-madgraph")
DATASETS.append(BASEDIR + "/ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8")
DATASETS.append(BASEDIR + "/TT_TuneCUETP8M1_13TeV-powheg-pythia8")

for idataset in DATASETS:
    cmd = "ssh cms25 condor_run " + os.environ.get("CMSSW_BASE") + "/src/babymaker/bmaker/genfiles/run/wrapper.sh " +os.environ.get("CMSSW_BASE") + "/src/babymaker/bmaker/genfiles/run/skim_ntuples.exe " + idataset + " " + OUTDIR + " " + CUT + " > /dev/null &"
    print cmd
    # hack to avoid too many ssh connections
    time.sleep(1)
    os.system(cmd)
