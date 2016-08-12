#!/bin/bash

### RA4 and ISR skims
infolder="/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed/"
./python/send_skim_ntuples.py $infolder 5 standard 
./python/send_skim_ntuples.py $infolder 5 met150	  
./python/send_skim_ntuples.py $infolder 5 baseline 
./python/send_skim_ntuples.py $infolder 5 qcd	  
./python/send_skim_ntuples.py $infolder 5 ttisr	  
./python/send_skim_ntuples.py $infolder 5 zisr	  
./python/send_skim_ntuples.py $infolder 5 wisr	  
./python/send_skim_ntuples.py $infolder 5 wisrht200

### RA2/b skims
infolder="/net/cms2/cms2r0/treemaker/2016_07_27/clean/"
./python/send_skim_ntuples.py $infolder 1 ra2_ht300   MET	       
./python/send_skim_ntuples.py $infolder 1 ra2_ht300   SingleMuon      
./python/send_skim_ntuples.py $infolder 1 ra2_eht300  SingleElectron  
./python/send_skim_ntuples.py $infolder 1 ra2_zmht200 SingleElectron  
./python/send_skim_ntuples.py $infolder 1 ra2_zmht200 SingleMuon      
./python/send_skim_ntuples.py $infolder 1 ra2_qcd     JetHT

infolder="/net/cms2/cms2r0/treemaker/2016_07_27/dirty/"
./python/send_skim_ntuples.py $infolder 1 ra2_ht300  SingleMuon
./python/send_skim_ntuples.py $infolder 1 ra2_eht300 SingleElectron
