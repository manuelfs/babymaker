#!/bin/bash

#### Script to quickly run babymaker interactively. Format is:
####  ./scripts/cmsRun.sh <inFile> <nEvents=1000> <outName>"

# The default inFile is 7.4.14 Monte Carlo to be run with CMSSW_7_4_14
inFile=/afs/cern.ch/user/m/manuelf/work/data/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2_74X_mcRun2_asymptotic_v2-v1_MINIAODSIM.root
# inFile=/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/0C22259A-F008-E511-B7A9-003048FFCB84.root
#inFile=/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/50A3A073-3B6C-E511-A997-02163E0144CD.root

nEvents="1000"
outName=""
json=do_not_want_json # If it can't find the file, it doesn't pre-apply the JSON

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

if [[ $inFile == *"store"* ]]
then
    cmsRun bmaker/python/bmaker_basic_cfg.py inputFiles=$inFile nEvents=$nEvents json=$json $outName 2>&1 | tee logout.log
else 
    cmsRun bmaker/python/bmaker_basic_cfg.py inputFiles=file:$inFile nEvents=$nEvents json=$json $outName 2>&1 | tee logout.log
fi

