#!/usr/bin/env python

import os, sys 
import glob
import json
import string, pprint
import ROOT
import das_client as das
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d","--datasets")
parser.add_argument("-f","--files")
args = parser.parse_args()

datasets = []
if (args.datasets):
  datasets = args.datasets.split(",")
elif (args.files):
  for fnm in args.files.split(","):
    with open(fnm) as f:
      datasets.extend([line for line in f.read().splitlines() if (len(line)>0 and line[0]=="/")])

print "Processing datasets:"
pprint.pprint(datasets)

# Parsing where the files can be found depends on whether we run on UCSB or UCSD
host = os.environ.get("HOSTNAME")
if "ucsd" in host: host = "sd"
elif host=="compute-0-1.local"or host =="compute-0-5.local" or host=="compute-0-0.local": host = "sb"
elif "ucsb" in host: sys.exit("\033[91mERROR: To allow access to hadoop use cms18 or cms19. \033[0m") 
else: sys.exit("\033[91mERROR: Unknown host: "+host+" Exit. \033[0m")

hadoop = '/mnt/hadoop/cms'
if host=="sd": hadoop = '/hadoop/cms/phedex'

# Directory to dump all condor-related logs, schell and cmd files
flistdir = os.path.join(os.getenv("CMSSW_BASE"),"src/flists/")
if not os.path.exists(flistdir):
  sys.exit("ERROR: flists repository not found.")


for ds in datasets:
  # parse the dataset name and guess the path on hadoop to create the input file list
  path,dsname,campaign,reco = '','','',''
  tags = string.split(ds,'/')
  dsname = tags[1]
  campaign = (string.split(tags[2],'-'))[0]
  reco = tags[2][len(campaign)+1:]
  filetype = tags[3]

  # query DAS; the result is a list of dictionaries; one dict per file
  # each dictionary has the name, size and nevents of a file
  # this is a per-dataset level query, so it's quick
  print "INFO: Query DAS for files in:", '_'.join([dsname,campaign,reco])
  this_fdicts = das.getFilesInfo(ds, wanted_keys = ['name','size','nevents'])
  nfiles = len(this_fdicts)
  if nfiles==0:
    print "\033[93m WARNING: "+ds+" not found! Skip dataset. \033[0m"
    continue
  
  fnm = flistdir+'/'+'_'.join(['flist',dsname,campaign,reco+'.txt'])
  do_chmod = True 
  if os.path.exists(fnm): do_chmod = False

  f = open(fnm,"w")
  nent_local = 0 
  nent = 0
  nfiles_local = 0
  for ifile in this_fdicts:
    runlist = ''
    if "Run2015" in ds: 
      # we need to ask DAS for what runs are in each file if it's reprocessed data
      # for prompt reco, the run number can be parsed from the file path
      if "PromptReco" in ds:
        runlist = ifile['name'].split("/000/").pop().split("/00000/")[0].replace("/","")
      else: 
        # this is a per-file query, so it's slow
        runlist = ','.join([str(irun) for irun in das.getFileRunInfo(ifile['name'])])

    # check if file is available locally
    if os.path.exists(hadoop+ifile['name']):
      if (os.path.getsize(hadoop+ifile['name'])==ifile['size']):
        if runlist=='': f.write("local  " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + '\n')
        else: f.write("local  " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + " " + runlist + '\n')
        nfiles_local = nfiles_local + 1
        nent_local = nent_local + ifile['nevents']
      else:
        print "Encountered partial local file, will use xrootd instead"
        if runlist=='': f.write("xrootd " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + '\n')
        else: f.write("xrootd " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + " " + runlist + '\n')
    else:
      if runlist=='': f.write("xrootd " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + '\n')
      else: f.write("xrootd " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + " " + runlist + '\n')
    nent = nent + ifile['nevents']
  f.write("nEventsTotal: "+'{:<10}'.format(nent)+ '\n')
  f.write("stat events available locally: " + "%s%% %i/%i" % ('{:.0f}'.format(float(nent_local)/float(nent)*100.),nent_local,nent) + '\n')
  f.write("stat files available locally: " + "%s%% %i/%i" % ('{:.0f}'.format(float(nfiles_local)/float(nfiles)*100.),nfiles_local,nfiles) + '\n')
  f.close()
  if (do_chmod): os.chmod(fnm, 0777)
