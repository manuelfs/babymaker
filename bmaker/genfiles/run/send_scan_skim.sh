#!/bin/bash

if (( "$#" < 3 ))
then
    echo
    echo "Format is: ./run/send_scan_skim.sh <infolder> <outfolder> <njobs> (optional: <additional cuts>)"
    echo 
    exit 0
fi


infolder=$1
outfolder=$2
njobs=$3

if(( "$#" >= 3 ))
then
    cuts=$4
else
    cuts=""
fi


nfiles=`ls $infolder/*root | wc -l`
for file in `seq 1 $njobs`;
do
    #echo JobSubmit.csh ./run/skim_scan.exe $infolder $outfolder \"$cuts\" $njobs $file
    JobSubmit.csh ./run/skim_scan.exe $infolder $outfolder $cuts $njobs $file
done
