#!/bin/bash

#### Script to quickly run babymaker interactively. Format is:
####  ./scripts/cmsRun.sh <inFile> <nEvents=1000> <outName>"

# The default inFile is 7.4.14 Monte Carlo to be run with CMSSW_7_4_14
# inFile=/afs/cern.ch/user/m/manuelf/work/data/SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2_MINIAODSIM.root
# inFile=/afs/cern.ch/user/m/manuelf/work/data/SMS-T1tttt_mGluino-1500to1525_mLSP-50to1125_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix_MINIAODSIM.root
# inFile=/afs/cern.ch/user/m/manuelf/work/data/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2_74X_mcRun2_asymptotic_v2-v1_MINIAODSIM.root
# inFile=/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/0C22259A-F008-E511-B7A9-003048FFCB84.root
#inFile=/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/50A3A073-3B6C-E511-A997-02163E0144CD.root
# inFile=/hadoop/cms/phedex/store/mc/RunIISpring16MiniAODv1/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/70000/8E926E3B-F304-E611-A352-001E67A3FB9B.root
# inFile=/hadoop/cms/phedex/store/mc/RunIISpring16MiniAODv2/ttZJets_13TeV_madgraphMLM/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/FEF1A376-641C-E611-AA2A-002590E7E00A.root
# inFile=/nfs-7/userdata/ana/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1.root
# inFile=/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/275/375/00000/FEDDBC5C-4339-E611-8740-02163E0136C4.root
# inFile=/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/275/066/00000/6CFE9029-6E34-E611-BA01-02163E011AC1.root
# inFile=/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/274/388/00000/A0B6C0F6-D32B-E611-8A5A-02163E014167.root
#inFile=/nfs-7/userdata/ana/TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1_f2.root
inFile=/hadoop/cms/phedex/store/mc/RunIISpring16MiniAODv2/SMS-T1qqqq_mGluino-1000_mLSP-800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/80000/BAE83A6A-0B38-E611-BF40-0025905B85DA.root
inFile=/nfs-7/userdata/ana/SMS-T1qqqq_mGluino-1400_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2.root
# inFile=/nfs-7/userdata/ana/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1.root
inFile=/net/cms2/cms2r0/babymaker/miniaod/RunIISpring16MiniAODv2__TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__MINIAODSIM__PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1__00000__001AFDCE-C33B-E611-B032-0025905D1C54.root

nEvents="10000"
# outName="fullbaby_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1.root"
# outName="fullbaby20_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1.root"
# outName="fullbaby_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1.root"
# outName="fullbaby_TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1.root"
# outName="fullbaby_SingleMuon_Run2016B_PromptReco-v2_run274388.root"
# outName="fullbaby_SMS-T1qqqq_mGluino-1000_mLSP-800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16.root"
# outName="fullbaby_SMS-T1qqqq_mGluino-1400_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16.root"
# outName="fullbaby_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1.root"
outName="fullbaby_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2_PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic.root"

json=do_not_want_json # If it can't find the file, it doesn't pre-apply the JSON
json="babymaker/data/json/golden_Cert_271036-275125_13TeV_PromptReco_Collisions16.json"

if (( "$#" > 0 ))
then
    inFile=$1
fi
if (( "$#" > 1 ))
then
    nEvents=$2
fi
if (( "$#" > 2 ))
then
    outName="outputFile="$3
fi

if [[ ($outName != *"outputFile"*) ]]
then
    outName="outputFile="$outName
fi
   
if [[ ($inFile == *"store"*) && ($inFile != *"hadoop"*) ]]
then
    cmsRun bmaker/python/bmaker_full_cfg.py inputFiles=$inFile nEvents=$nEvents json=$json $outName 2>&1 | tee logout.log
else 
    cmsRun bmaker/python/bmaker_full_cfg.py inputFiles=file:$inFile nEvents=$nEvents json=$json $outName 2>&1 | tee logout.log
fi

