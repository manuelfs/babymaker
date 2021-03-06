#!/bin/bash

#### Script to quickly run babymaker interactively. Format is:
####  ./scripts/cmsRun.sh <inFile> <nEvents=1000> <outName>"

## Storage location of files at UCSD/lxplus, two examples:
# inFile=/afs/cern.ch/user/m/manuelf/work/data/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2_74X_mcRun2_asymptotic_v2-v1_MINIAODSIM.root
# inFile=/nfs-7/userdata/ana/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1_PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1.root

# inFile=/net/cms2/cms2r0/babymaker/miniaod/RunIISpring16MiniAODv2__TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8__MINIAODSIM__PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1__00000__001AFDCE-C33B-E611-B032-0025905D1C54.root
# outName="fullbaby_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2_PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic.root"

#inFile=/net/cms2/cms2r0/babymaker/miniaod/TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root
#outName="fullbaby_TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root"

# inFile=/net/cms2/cms2r0/babymaker/miniaod/SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2.root
# outName="fullbaby_SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2.root"

# inFile=/net/cms2/cms2r0/babymaker/miniaod/MET_Run2016E_PromptReco-v2_1666866B-394D-E611-93B2-FA163E6A5A26_Run276831.root
# outName="fullbaby_MET_Run2016E_PromptReco-v2_1666866B-394D-E611-93B2-FA163E6A5A26_Run276831.root"

# inFile=/net/cms2/cms2r0/babymaker/miniaod/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root
# outName=fullbaby_SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root

#inFile=/store/data/Run2016F/SingleElectron/MINIAOD/PromptReco-v1/000/277/981/00000/08034EC5-C759-E611-9E7B-FA163E4E25AC.root
#inFile=/store/data/Run2016G/MET/MINIAOD/PromptReco-v1/000/280/385/00000/889C6ADD-A078-E611-9996-02163E011A1B.root
#outName="fullbaby_MET_Run2016G_PromptReco-v1.root"

inFile=/store/data/Run2016C/MET/MINIAOD/23Sep2016-v1/70000/027548F7-B28A-E611-A1FB-6CC2173BBAB0.root
outName="fullbaby_MET_Run2016C_23Sep2016-v1_2_oldfilter.root"

inFile=/store/data/Run2016C/SingleElectron/MINIAOD/23Sep2016-v1/50000/00972D76-E68C-E611-B134-0025905A607E.root
outName="fullbaby_SingleElectron_Run2016C_23Sep2016-v1.root"

# inFile=/home/users/ana/CMSSW_8_0_16/src/data/SMS-TChiHH_HToWWZZTauTau_HToWWZZTauTau_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root
# outName="fullbaby_SMS-TChiHH_HToWWZZTauTau_HToWWZZTauTau_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root"

inFile=/home/users/ana/data/SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root
outName="fullbaby_SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2_PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root"

inFile=/home/users/ana/data/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2.root
outName="fullbaby_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2.root"

inFile=/home/users/ana/data/Run2016H_MET_03Feb2017_ver3-v1.root
outName="fullbaby_Run2016H_MET_03Feb2017_ver3-v1.root"

inFile=/hadoop/cms/phedex/store/data/Run2017D/MET/MINIAOD/17Nov2017-v1/50000/06B11CC2-F9EA-E711-826C-A4BF0101202F.root
outName="fullbaby_Run2017D_MET_17Nov2017-v1.root"

inFile=/eos/cms/store/data/Run2018B/MET/MINIAOD/PromptReco-v1/000/317/320/00000/A2062E5A-0968-E811-8E97-02163E01A0B3.root
outName="fullbaby_Run2018B_MET_PromptReco-v1.root"

# inFile=ZprimeToWWToWlepWhad_narrow_M-3000_TuneCP5_13TeV-madgraph_RunIIFall17MiniAOD_94X_mc2017_realistic_v10-v1.root
# outName="fullbaby_TTJets_TuneCP5_13TeV-madgraph_RunIIFall17MiniAOD.root"

# inFile=/home/users/ana/data/TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2.root
# outName="fullbaby_TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2.root"

# inFile=/home/users/ana/data/TTJets_SingleLeptFromT_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2.root
# outName="TTJets_SingleLeptFromT_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2.root"

nEvents="1000"

json=do_not_want_json # If it can't find the file, it doesn't pre-apply the JSON
# json="babymaker/data/json/golden_Cert_271036-275125_13TeV_PromptReco_Collisions16.json"

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
   
cmsRun bmaker/python/bmaker_full_cfg.py inputFiles=file:$inFile nEvents=$nEvents json=$json $outName 2>&1 | tee logout.log

