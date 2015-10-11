#!/usr/bin/env python

import os, sys 
import glob
import string, pprint
import ROOT
import das_client as das

datasets = []
with open("txt/dataset/mc.txt") as fmc:
  datasets = [line for line in fmc.read().splitlines() if (len(line)>0 and line[0]=="/")]

# Parsing where the files can be found depends on whether we run on UCSB or UCSD
host = os.environ.get("HOSTNAME")
if "ucsd" in host: host = "sd"
elif host=="compute-0-1.local" or host=="compute-0-0.local": host = "sb"
elif "ucsb" in host: sys.exit("\033[91mERROR: To allow access to hadoop use cms18 or cms19. \033[0m") 
else: sys.exit("\033[91mERROR: Unknown host: "+host+" Exit. \033[0m")

hadoop = '/mnt/hadoop/cms'
if host=="sd": hadoop = '/hadoop/cms/phedex'

# Directory to dump all condor-related logs, schell and cmd files
flistdir = "run/"
if host=="sd":
    if not (os.path.exists(os.getcwd()+'/'+flistdir)):
        os.mkdir(os.getcwd()+'/'+flistdir)
else:
    flistdir = "/net/cms2/cms2r0/babymaker/flist/"

for ds in datasets:
  # parse the dataset name and guess the path on hadoop to create the input file list
  path,dsname,campaign,reco = '','','',''
  tags = string.split(ds,'/')
  dsname = tags[1]
  campaign = (string.split(tags[2],'-'))[0]
  reco = tags[2][len(campaign)+1:]
  filetype = tags[3]

  # query DAS, result is a list of dictionaries
  # each dictionary has the name, size and nevents of a file
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
    if os.path.exists(hadoop+ifile['name']):
      if (os.path.getsize(hadoop+ifile['name'])==ifile['size']):
        f.write("local  " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + '\n')
        nfiles_local = nfiles_local + 1
        nent_local = nent_local + ifile['nevents']
      else:
        print "Encountered partial local file, will use xrootd instead"
        f.write("xrootd " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + '\n')
    else:
      f.write("xrootd " + '{:<10}'.format(ifile['nevents']) + ifile['name'] + '\n')
    nent = nent + ifile['nevents']
  f.write("nEventsTotal: "+'{:<10}'.format(nent)+ '\n')
  f.write("stat events available locally: " + "%s%% %i/%i" % ('{:.0f}'.format(float(nent_local)/float(nent)*100.),nent_local,nent) + '\n')
  f.write("stat files available locally: " + "%s%% %i/%i" % ('{:.0f}'.format(float(nfiles_local)/float(nfiles)*100.),nfiles_local,nfiles) + '\n')
  f.close()
  if (do_chmod): os.chmod(fnm, 0777)
