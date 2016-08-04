#!/usr/bin/env python

import os, math
from ROOT import *

do80X = True # false will calculate PU weights for 74X MC w.r.t. 2015 data

# ----- NEEDED INPUTS --------------------------------------
# Instructions:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData

# retrieve latest json with per-lumiblock luminosity information
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
if do80X: lumi_json = "/net/cms2/cms2r0/babymaker/misc/pu/2016/pileup_latest.txt"
else: lumi_json = "/net/cms2/cms2r0/babymaker/misc/pu/2015/pileup_latest.txt"

# obtain certified lumis JSON for dataset we are trying to reweight to
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/
# json = "data/json/golden_Cert_271036-276811_13TeV_PromptReco_Collisions16.json"
if do80X: json = "/net/cms2/cms2r0/babymaker/misc/pu/2016/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt"
else: json = "/net/cms2/cms2r0/babymaker/misc/pu/2015/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt"

# Minbias cross-section, (corrected) announcement for ICHEP16:
# https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/613/2/1/1/1.html
if do80X: mbxsec, mbxsec_relunc = 69200, 0.046
else: mbxsec, mbxsec_relunc = 69000, 0.05 

# get hist of the MC PU profile
mcpu_vals = []
if (do80X):
    # https://github.com/cms-sw/cmssw/blob/2595fa9fffd50a63220711e20385df5e234a7c42/SimGeneral/MixingModule/python/mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi.py
    mcpu_vals = [0.000829312873542, 0.00124276120498, 0.00339329181587, 0.00408224735376, 0.00383036590008, 0.00659159288946, 0.00816022734493, 0.00943640833116, 0.0137777376066, 0.017059392038, 0.0213193035468, 0.0247343174676, 0.0280848773878, 0.0323308476564, 0.0370394341409, 0.0456917721191, 0.0558762890594, 0.0576956187107, 0.0625325287017, 0.0591603758776, 0.0656650815128, 0.0678329011676, 0.0625142146389, 0.0548068448797, 0.0503893295063, 0.040209818868, 0.0374446988111, 0.0299661572042, 0.0272024759921, 0.0219328403791, 0.0179586571619, 0.0142926728247, 0.00839941654725, 0.00522366397213, 0.00224457976761, 0.000779274977993, 0.000197066585944, 7.16031761328e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
else: 
    # https://github.com/cms-sw/cmssw/blob/fd6c90c51357635870727e9e0a85959be7b9463f/SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.py
    mcpu_vals = [4.8551E-07, 1.74806E-06, 3.30868E-06, 1.62972E-05, 4.95667E-05, 0.000606966, 0.003307249, 0.010340741, 0.022852296, 0.041948781, 0.058609363, 0.067475755, 0.072817826, 0.075931405, 0.076782504, 0.076202319, 0.074502547, 0.072355135, 0.069642102, 0.064920999, 0.05725576, 0.047289348, 0.036528446, 0.026376131, 0.017806872, 0.011249422, 0.006643385, 0.003662904, 0.001899681, 0.00095614, 0.00050028, 0.000297353, 0.000208717, 0.000165856, 0.000139974, 0.000120481, 0.000103826, 8.88868E-05, 7.53323E-05, 6.30863E-05, 5.21356E-05, 4.24754E-05, 3.40876E-05, 2.69282E-05, 2.09267E-05, 1.5989E-05, 4.8551E-06, 2.42755E-06, 4.8551E-07, 2.42755E-07, 1.21378E-07, 4.8551E-08]

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

hmc = TH1D("hmc","hmc",len(mcpu_vals),0,len(mcpu_vals))
for i, val in enumerate(mcpu_vals):
    hmc.SetBinContent(i+1, val)
hmc.Scale(1./hmc.Integral())

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
    


