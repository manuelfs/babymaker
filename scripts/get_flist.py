#!/usr/bin/env python

import os, sys 
import glob
import string
import ROOT
from ROOT import TChain

# silence ROOT
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")

# List of datasets to run over, could be either MC or data
# Enter either the dataset path starting with "/store", e.g:
# /store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/
# or the dataset name as found in DAS, e.g:
# /TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM
datasets = []
datasets.append('/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM')
# datasets.append('/store/mc/RunIISpring15DR74/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1')
# datasets.append('/store/mc/RunIISpring15DR74/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9_ext1-v1')
# datasets.append('/store/data/Run2015B/HTMHT/MINIAOD/PromptReco-v1')


# Directory to dump all condor-related logs, schell and cmd files
rundir = "run/"
if not (os.path.exists(os.getcwd()+'/'+rundir)):
    os.mkdir(os.getcwd()+'/'+rundir)

# Parsing where the files can be found depends on whether we run on UCSB or UCSD
host = os.environ.get("HOSTNAME")
if "ucsd" in host: host = "sd"
elif host=="compute-0-1.local" or host=="compute-0-0.local": host = "sb"
elif "ucsb" in host: sys.exit("\033[91mERROR: To allow access to hadoop use cms18 or cms19. \033[0m") 
else: sys.exit("\033[91mERROR: Unknown host: "+host+" Exit. \033[0m")

hadoop = '/mnt/hadoop/cms'
if host=="sd": hadoop = '/hadoop/cms/phedex'

for ds in datasets:
    # parse the dataset name and guess the path on hadoop to create the input file list
    path,dsname,campaign,reco = '','','',''
    if (ds[0:6]=='/store'):
        tags = string.split(ds,'/')
        dsname = tags[4]
        campaign = tags[3]
        reco = tags[6]
        filetype = tags[5]
    else:
        tags = string.split(ds,'/')
        dsname = tags[1]
        campaign = (string.split(tags[2],'-'))[0]
        reco = tags[2][len(campaign)+1:]
        filetype = tags[3]
    if 'PromptReco' in ds: path = '/'.join([hadoop+'/store/data',campaign,dsname,filetype,reco,'*/*/*/*/*root'])
    else: path = '/'.join([hadoop+'/store/mc',campaign,dsname,filetype,reco,'*/*root'])

    filelist = glob.glob(path)
    nfiles = len(filelist)
    if nfiles==0:
        print "\033[93m WARNING: "+ds+" not found! Skip dataset. \033[0m"
        continue
    
    fnm = '_'.join(['flist',dsname,campaign,reco+'.txt'])
    f = open(rundir+'/'+fnm,"w")
    for i,ifile in enumerate(filelist):
        tree = TChain("Events")
        tree.Add(ifile)
        nent = tree.GetEntries()
        ifile = string.replace(ifile,hadoop,'')
        f.write('{:<10}'.format(nent)+ifile+'\n')
    f.close()
