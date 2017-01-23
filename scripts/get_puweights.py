#!/usr/bin/env python

import os, math
from ROOT import *

# ----- NEEDED INPUTS --------------------------------------
# Instructions:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData

cmssw = os.getenv("CMSSW_BASE")+"/src/"
# retrieve latest json with per-lumiblock luminosity information
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
# fyi, how it gets derived: https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/660/1.html
lumi_json = cmssw + "babymaker/data/json/perlumiblock_pileup.json"

# obtain certified lumis JSON for dataset we are trying to reweight to
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/
json = cmssw + "babymaker/data/json/golden_Cert_271036-284044_13TeV_PromptReco_Collisions16.json"

# Minbias cross-section, (corrected) announcement for ICHEP16:
# https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/613/2/1/1/1.html
mbxsec, mbxsec_relunc = 69200, 0.046

# get hist of the MC PU profile
gROOT.SetBatch(kTRUE)
mcfile = TChain("tree") 
mcfile.Add("/hadoop/cms/store/user/ana/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_SingleLeptFromT_13TeV-MG-PY8__Moriond17_80XMiniAODv2_v6-v1__MINIAODSIM/170121_073722/0000/fullbaby_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_*")
# mcfile.Add("/net/cms2/cms2r0/babymaker/babies/2017_01_21/mc/unprocessed/fullbaby_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_*")
hmc = TH1D("hmc","hmc",75,0,75)
mcfile.Draw("ntrupv_mean>>hmc","","norm")

# To reweight both the in- and out-of-time pile up
# pileupCalc should be used in mode "true" and the variable to reweight is npvtru_mean:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookExercisesHATS_Pileup_2013
calc_mode = "true" 

# ----------------------------------------------------------

mbxsec_dict = {
     "nom": mbxsec,
      "up": mbxsec*(1+mbxsec_relunc),
    "down": mbxsec*(1-mbxsec_relunc)
}

for imb in mbxsec_dict:
    cmd  = "pileupCalc.py"
    cmd += " -i "+json
    cmd += " --inputLumiJSON "+lumi_json
    cmd += " --calcMode "+calc_mode
    cmd += " --minBiasXsec "+str(mbxsec_dict[imb])
    cmd += " --maxPileupBin "+str(hmc.GetNbinsX())
    cmd += " --numPileupBins "+str(hmc.GetNbinsX())
    cmd += " pileup_"+imb+".root"

    print "Obtaining data pile up distribution for variation:", imb
    print cmd
    os.system(cmd)

    fdata = TFile("pileup_"+imb+".root","READ")
    htmp = fdata.Get("pileup").Clone("pu_"+imb)
    htmp.Scale(1./htmp.Integral())
    htmp.Divide(hmc)
    htmp.SetDirectory(0)
    fdata.Close()
    wgt = [htmp.GetBinContent(i+1) for i in range(htmp.GetNbinsX())]
    if imb=="nom": 
        print "Nominal weights:"
        for j,iwgt in enumerate(wgt):
            print "NPV: "+'{:>3d}'.format(j+1),
            print " Weight: "+'{:>10.3e}'.format(iwgt)

    print "------> Vector for weight_tools:"
    print "w_pu_"+imb,
    print " = vector<double>({"+', '.join('{:.3e}'.format(x) for x in wgt)+"});"
    


