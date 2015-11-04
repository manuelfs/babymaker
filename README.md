babymaker
==============

CMSSW code to generate babies (small flat ntuples) from MINIAOD

To set up the code and generate a file named `baby.root`, issue the following commands 
on lxplus:

    cmsrel CMSSW_7_4_12_patch4
    cd CMSSW_7_4_12_patch4/src
    cmsenv
    git-cms-addpkg CondFormats/BTauObjects
    git clone git@github.com:manuelfs/babymaker
    cd babymaker
    ./compile.sh
    cmsRun bmaker/python/bmaker_basic_cfg.py inputFiles=/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v3/000/256/675/00000/5E85670E-9C5F-E511-AD41-02163E011E80.root outputFile=baby_Run2015D_SingleMuon.root nEvents=100

The `compile.sh` script first compiles the `babymaker/bmaker/genfiles` folder, which
automatically generates the tree structure (see below), and then issues `scram b`
in the `babymaker` folder. 

#### Adding new branches

To add new branches to the tree, you first create the new branch in
`babymaker/variables/basic` where the type and name are specified.
The code in `babymaker/bmaker/genfiles/src/generate_baby.cxx` automatically generates
the files 

    babymaker/bmaker/interface/baby_base.hh
    babymaker/bmaker/interface/baby_basic.hh
    babymaker/bmaker/src/baby_base.cc
    babymaker/bmaker/src/baby_basic.cc

which have the c++ code that defines the class `baby_basic` with the tree, all the branches,
automatic vector clearing in the `baby_basic::Fill` method, and other useful functions.


The branches are filled in the `bmaker_basic::analyze` method in 
`babymaker/bmaker/plugins/bmaker_basic.cc`. Functions that define physics quantities,
like isolation or electron ID, are defined in `babymaker/bmaker/src/phys_objects.cc`.

#### Submitting jobs to condor in the T3 at UCSB

For now, make sure that `process.maxEvents` is set to the number of events you want
(-1 for all events) in `babymaker/bmaker/python/bmaker_basic_cfg.py`.

Log on to cms18 and set up the code are described above.
Define the datasets you want to run over in `babymaker/scripts/get_flist.py` and run it
from the `babymaker` folder with 

    ./scripts/get_flist.py 

This step finds the paths for the files that are to be run over.

In the `babymaker` folder, create an `out` folder. This is typically a soft link to a place
with lots of disk space, such as `ln -s /net/cms2/cms2r0/user/out/ out`.

Then log on to cms0 and submit the condor jobs from `babymaker` folder with

    ./scripts/sub_cond.py

