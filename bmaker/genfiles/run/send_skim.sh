#!/bin/bash

if (( "$#" < 4 ))
then
    echo
    echo "Format is: ./run/send_skim.sh <infolder> <outfolder> <cuts> <njobs>"
    echo 
    exit 0
fi

infolder=$1
outfolder=$2
cuts=$3
njobs=$4

nfiles=`ls $infolder/*root | wc -l`
for file in `seq 1 $njobs`;
do
    #echo JobSubmit.csh ./plot/skim_ntuples.exe $infolder $outfolder \"$cuts\" $njobs $file
    JobSubmit.csh ./plot/skim_ntuples.exe $infolder $outfolder $cuts $njobs $file
done
