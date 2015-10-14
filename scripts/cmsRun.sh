#!/bin/bash

if (( "$#" < 1 ))
then
    echo
    echo "Format is: ./scripts/cmsRun.sh inFile <nEvents=1000> <outName>"
    echo 
    exit 0
fi

inFile=$1
nEvents="1000"
outName=""
if (( "$#" > 1 ))
then
    nEvents=$2
fi
if (( "$#" > 2 ))
then
    outName="outputFile="$3
fi

cmsRun bmaker/python/bmaker_basic_cfg.py inputFiles=file:$inFile nEvents=$nEvents $outName


