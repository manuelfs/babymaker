babymaker
==============

CMSSW code to generate babies (small flat ntuples) from MINIAOD

  * [Code setup](#code-setup)
  * [Adding new branches](#adding-new-branches)
  * [Getting an flist](#getting-an-flist)
  * [Submitting jobs to condor in the T3 at UCSB](#submitting-jobs-to-condor-in-the-T3-at-UCSB)
  * [Conventions in babymaker](#conventions-in-babymaker)


#### Code setup
To set up the code and generate a file named `baby.root`, issue the following commands 
on lxplus:

    cmsrel CMSSW_8_0_20
    cd CMSSW_8_0_20/src
    cmsenv
    git cms-init
    git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate	
    git clone git@github.com:manuelfs/babymaker
    cd babymaker
    ./compile.sh
    ./scripts/cmsRun.sh <inputfile> <nevents=1000> <outname>

The `compile.sh` script first compiles the `babymaker/bmaker/genfiles` folder, which
automatically generates the tree structure (see below), and then issues `scram b`
in the `babymaker` folder. 

#### Adding new branches

To add new branches to the tree, you first create the new branch in
`babymaker/variables/full` where the type and name are specified.
The code in `babymaker/bmaker/genfiles/src/generate_baby.cxx` automatically generates
the files 

    babymaker/bmaker/interface/baby_base.hh
    babymaker/bmaker/interface/baby_full.hh
    babymaker/bmaker/src/baby_base.cc
    babymaker/bmaker/src/baby_full.cc

which have the c++ code that defines the class `baby_full` with the tree, all the branches,
automatic vector clearing in the `baby_full::Fill` method, and other useful functions.

The branches are filled in the `bmaker_full::analyze` method in 
`babymaker/bmaker/plugins/bmaker_full.cc`. Functions that define physics quantities,
like isolation or electron ID, are defined in `babymaker/bmaker/src/*_tools.cc`.


#### Getting an flist

To process entire datasets you need the corresponding flist, or list of individual files in the dataset.
A number of these flists are committed to trandbea/flists, and are downloaded with

    git clone git@github.com:trandbea/flists

This repository must be placed at `${CMSSW_BASE}/src`, not inside `babymaker`.

To obtain an flist for a new dataset, or a number of datasets, you use the following script

    ./scripts/get_flist.py -d <dataset_name>
    ./scripts/get_flist.py -f <file_with_dataset_names>

All flists must be placed in the `${CMSSW_BASE}/src/flists` folder.


#### Submitting jobs to condor in the T3 at UCSB

Log on to one of the compute-0-X machines and set up the code are described above.
In the `babymaker` folder, create an `out` folder. This is typically a soft link to a place
with lots of disk space, such as `ln -s /net/cms2/cms2r0/user/out/ out`.

Then, add the datasets you want to run over in the wishlist of `scripts/sub_cond.py` and submit the condor jobs from
`babymaker` folder with

    ./scripts/sub_cond.py

#### Post-processing of produced babies

See [bmaker/genfiles/README.md](bmaker/genfiles/README.md) for detailed instructions.

#### Conventions in babymaker

In order to homogenize the code and know what to expect, we try to follow these conventions in the  development
of `babymaker`:

 * **Branch names** use **all lower case letters**. If words must be separated, use an underscore, e.g. `met_phi`
 * **File names** use **all lower case letters**. If words must be separated, use an underscore, e.g. `lepton_tools.cc`
 * **Function names** follow the standard **CMSSW convention**, that is, first word all lower case, and first letter 
 of subsequent words in upper case. e.g. `bmaker_full::writeFatJets`
 * **Product names** (e.g. `"slimmedElectrons"`) are **only defined in `bmaker_full::analize` or in 
 `babymaker/bmaker/python/bmaker_full_cfg.py`**. 
 * As much as possible, **physics definitions go in the corresponding `src/*_tools.cc` file**, e.g. lepton ID goes in
 `src/lepton_tools.cc`. They should not be part of `plugins/bmaker_full.cc` so that when/if we move to having various 
 baby definitions in parallel, for which the code is setup, all babies would use common definitions.
