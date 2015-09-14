babymaker
==============

Code to flat ntuples from MINIAOD in CMSSW

Issue the following commands on lxplus:

    cmsrel CMSSW_7_4_6_patch6
    cd CMSSW_7_4_6_patch6/src
    cmsenv
    git clone git@github.com:manuelfs/babymaker
    scram b -j$(getconf _NPROCESSORS_ONLN)
    cmsRun cmsRun babymaker/bmaker/python/bmaker_cfg.py