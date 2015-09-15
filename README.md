babymaker
==============

CMSSW code to generate babies (small flat ntuples) from MINIAOD

To generate a file named baby.root, issue the following commands on lxplus:

    cmsrel CMSSW_7_4_6_patch6
    cd CMSSW_7_4_6_patch6/src
    cmsenv
    git clone git@github.com:manuelfs/babymaker
    cd babymaker
    ./compile.sh
    cmsRun cmsRun babymaker/bmaker/python/bmaker_basic_cfg.py

To add new branches to the tree, you specify the type and name in
`babymaker/variables/basic`.