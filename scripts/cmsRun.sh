#!/bin/bash

#### Script to quickly run babymaker interactively. Format is:
####  ./scripts/cmsRun.sh <inFile> <nEvents=1000> <outName>"

# The default inFile is 7.4.14 Monte Carlo to be run with CMSSW_7_4_14
inFile=/afs/cern.ch/user/m/manuelf/work/data/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2_74X_mcRun2_asymptotic_v2-v1_MINIAODSIM.root
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

cmsRun bmaker/python/bmaker_basic_cfg.py inputFiles=file:$inFile nEvents=$nEvents json=$json $outName 2>&1 | tee logout.log


