#!/bin/bash

infolder="/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed/"
basefolder="/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/skim"


./run/send_skim.sh $infolder ${basefolder}_standard/   5 standard
./run/send_skim.sh $infolder ${basefolder}_met150/     5 met150
./run/send_skim.sh $infolder ${basefolder}_baseline/   5 baseline
./run/send_skim.sh $infolder ${basefolder}_qcd/        5 qcd
./run/send_skim.sh $infolder ${basefolder}_ttisr/      5 ttisr
./run/send_skim.sh $infolder ${basefolder}_zisr/       5 zisr
./run/send_skim.sh $infolder ${basefolder}_wisr/       5 wisr
./run/send_skim.sh $infolder ${basefolder}_wisrht200/  5 wisrht200
