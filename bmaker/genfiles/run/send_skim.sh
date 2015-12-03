#!/bin/bash

if (( "$#" < 3 ))
then
    echo
    echo "Format is: ./run/send_skim.sh <infolder> <outfolder> <njobs> <cuts\"nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj>250\">"
    echo 
    exit 0
fi

infolder=$1
outfolder=$2
njobs=$3

if (( "$#" < 4 ))
then
    cuts="nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj>250"
else
    cuts=$4
fi

nfiles=`ls $infolder/*root | wc -l`
for file in `seq 1 $njobs`;
do
    JobSubmit.csh ./run/wrapper.sh ./run/skim_ntuples.exe $infolder $outfolder $cuts $njobs $file
done
