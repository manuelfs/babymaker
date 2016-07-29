#!/bin/bash

### RA4 and ISR skims
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


### RA2/b skims
infolder="/net/cms2/cms2r0/treemaker/2016_07_27/clean/"
basefolder=${infolder}/skim
./run/send_skim.sh $infolder ${basefolder}_ra2_ht300/   1  ra2_ht300   MET
./run/send_skim.sh $infolder ${basefolder}_ra2_ht300/   1  ra2_ht300   SingleMuon
./run/send_skim.sh $infolder ${basefolder}_ra2_eht300/  1  ra2_eht300  SingleElectron
./run/send_skim.sh $infolder ${basefolder}_ra2_zmht200/ 1  ra2_zmht200 SingleElectron
./run/send_skim.sh $infolder ${basefolder}_ra2_zmht200/ 1  ra2_zmht200 SingleMuon
./run/send_skim.sh $infolder ${basefolder}_ra2_qcd/     1  ra2_qcd     JetHT

infolder="/net/cms2/cms2r0/treemaker/2016_07_27/dirty/"
basefolder=${infolder}/skim
./run/send_skim.sh $infolder ${basefolder}_ra2_ht300/  1  ra2_ht300  SingleMuon
./run/send_skim.sh $infolder ${basefolder}_ra2_eht300/ 1  ra2_eht300 SingleElectron

