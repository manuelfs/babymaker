#!/usr/bin/env python

import os, sys, subprocess
import pprint
import glob
import json
import string
import time
import argparse

#What to submit? Use substrings that would be found in the desired dataset, no wild cards!
# e.g. if we want only 25ns TTJets, use a substring that contains it all:
# "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns"
# if we want all the HT-binned TTJets append:
# "TTJets_HT-"
mc_wishlist = []
mc_wishlist.append("SMS-T1tttt_mGluino-1000to1050_mLSP-1to800");
mc_wishlist.append("SMS-T1tttt_mGluino-1050_mLSP-50to775");
mc_wishlist.append("SMS-T1tttt_mGluino-1100_mLSP-1to775");
mc_wishlist.append("SMS-T1tttt_mGluino-1100to1125_mLSP-700to900");
mc_wishlist.append("SMS-T1tttt_mGluino-1150_mLSP-1to800");
mc_wishlist.append("SMS-T1tttt_mGluino-1150to1175_mLSP-750to925");
mc_wishlist.append("SMS-T1tttt_mGluino-1200_mLSP-1to825");
mc_wishlist.append("SMS-T1tttt_mGluino-1200to1225_mLSP-800to1000");
mc_wishlist.append("SMS-T1tttt_mGluino-1250to1275_mLSP-700to1050");
mc_wishlist.append("SMS-T1tttt_mGluino-1275_mLSP-900to975");
mc_wishlist.append("SMS-T1tttt_mGluino-1300_mLSP-1to1075");
mc_wishlist.append("SMS-T1tttt_mGluino-1300to1325_mLSP-700to1100");
mc_wishlist.append("SMS-T1tttt_mGluino-1325to1350_mLSP-1to1125");
mc_wishlist.append("SMS-T1tttt_mGluino-1350to1375_mLSP-50to1025");
mc_wishlist.append("SMS-T1tttt_mGluino-1400_mLSP-1to1175");
mc_wishlist.append("SMS-T1tttt_mGluino-1400to1425_mLSP-50to1100");
mc_wishlist.append("SMS-T1tttt_mGluino-1425to1450_mLSP-1to1200");
mc_wishlist.append("SMS-T1tttt_mGluino-1450to1475_mLSP-50to1075");
mc_wishlist.append("SMS-T1tttt_mGluino-1475to1500_mLSP-1to1250");
mc_wishlist.append("SMS-T1tttt_mGluino-1500to1525_mLSP-50to1125");
mc_wishlist.append("SMS-T1tttt_mGluino-1525to1550_mLSP-1to1300");
mc_wishlist.append("SMS-T1tttt_mGluino-1600to1650_mLSP-1to1350");
mc_wishlist.append("SMS-T1tttt_mGluino-1650to1700_mLSP-1to1400");
mc_wishlist.append("SMS-T1tttt_mGluino-1700to1750_mLSP-1to1450");
mc_wishlist.append("SMS-T1tttt_mGluino-1800to1850_mLSP-1to1450");
mc_wishlist.append("SMS-T1tttt_mGluino-1850to1900_mLSP-1to1450");
mc_wishlist.append("SMS-T1tttt_mGluino-1900to1950_mLSP-0to1450");
mc_wishlist.append("SMS-T1tttt_mGluino-1950_mLSP-700to950");
mc_wishlist.append("SMS-T1tttt_mGluino-600_mLSP-250to325");
mc_wishlist.append("SMS-T1tttt_mGluino-625_mLSP-275to375");
mc_wishlist.append("SMS-T1tttt_mGluino-625to650_mLSP-200to400");
mc_wishlist.append("SMS-T1tttt_mGluino-650to675_mLSP-250to425");
mc_wishlist.append("SMS-T1tttt_mGluino-675_mLSP-325to450");
mc_wishlist.append("SMS-T1tttt_mGluino-700_mLSP-1to450");
mc_wishlist.append("SMS-T1tttt_mGluino-700to750_mLSP-200to500");
mc_wishlist.append("SMS-T1tttt_mGluino-750to775_mLSP-350to525");
mc_wishlist.append("SMS-T1tttt_mGluino-775_mLSP-475to550");
mc_wishlist.append("SMS-T1tttt_mGluino-800to825_mLSP-1to575");
mc_wishlist.append("SMS-T1tttt_mGluino-825to850_mLSP-200to600");
mc_wishlist.append("SMS-T1tttt_mGluino-850to875_mLSP-450to625");
mc_wishlist.append("SMS-T1tttt_mGluino-875to900_mLSP-1to650");
mc_wishlist.append("SMS-T1tttt_mGluino-950to975_mLSP-350to725");
mc_wishlist.append("SMS-T1tttt_mGluino-975_mLSP-600to750");
#------------ MC ------------------
# mc_wishlist.append("T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
# mc_wishlist.append("T1tttt_mGluino-1200_mLSP-800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")

# mc_wishlist.append("TTJets_HT")
# mc_wishlist.append("TTJets_SingleLep")
# mc_wishlist.append("TTJets_DiLep")
# mc_wishlist.append("TTGJets")
# mc_wishlist.append("WJets")
# mc_wishlist.append("DYJets")
# mc_wishlist.append("DYJetsToLL_M-50_TuneCUETP8M1_13TeV")
# mc_wishlist.append("QCD_")
# mc_wishlist.append("ST_")
# mc_wishlist.append("ttHJetTobb")
# mc_wishlist.append("WWTo2L2Nu")
# mc_wishlist.append("WWToLNuQQ")
# mc_wishlist.append("ggZH_HToBB")

data_wishlist = []
#data_wishlist.append("JetHT")
#data_wishlist.append("HTMHT")
#data_wishlist.append("MET")
#data_wishlist.append("SingleElectron")
#data_wishlist.append("SingleMuon")
#data_wishlist.append("DoubleEG")
#data_wishlist.append("DoubleMuon")

jsonlist = glob.glob("data/json/subgolden_*.json")

# for data get the golden runs
goldruns = {}
for jsonfile in jsonlist:
  jdata = {}
  with open(jsonfile) as jfile:
    jdata = json.load(jfile)
  goldruns[jsonfile] = [int(i) for i in jdata.keys()]

# Maximum number of input MINIAOD files per condor job
maxfiles = int(raw_input('Enter max number of files per job: '))

# These keys, one for mc and one for each data period, are used to split the name of a dataset in two parts, 
# if the substring preceeding the key is the same for two datasets, 
# they are considered extensions of each other and combined
# i.e. the output babies and logs are labeled by 'substring-before-key'+'key'
# the sub-string following the key is dropped and forgotten!
comb_keys = ['RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9','RunIISpring15FSPremix_MCRUN2_74_V9','RunIISpring15MiniAODv2_74X_mcRun2_asymptotic_v2','Run2015D']

# for testing... otherwise set to -1
maxjobs = -1
maxds = -1
maxevents_perjob = -1
istest = 'n'
if (maxjobs!=-1 or maxevents_perjob!=-1 or maxds!=-1): 
  istest = raw_input('Running in test mode with %i jobs, %i events per job, over %i datasets. Enter \'y\' to continue: ' % (maxjobs, maxevents_perjob, maxds))
  if (istest!='y'):
    sys.exit("No jobs submitted. Edit sub_cond.py to exit test mode.")

# Only matters if running on UCSD:
# To run on multiple T2's use, e.g:
whitelist = "T2_US_UCSD,T2_US_WISCONSIN,T2_US_FLORIDA,T2_US_PURDUE,T2_US_NEBRASKA,T2_US_CALTECH"
# whitelist = "T2_US_WISCONSIN,T2_US_FLORIDA,T2_US_PURDUE,T2_US_NEBRASKA,T2_US_CALTECH"
# Need to check which is better, running on remote T2 or xrootd-ing the data...
# To run only at UCSD use:
# whitelist = "T2_US_UCSD"

# Condor set up depends on whether we run on UCSB or UCSD
host = os.environ.get("HOSTNAME")
if "ucsd" in host: host = "sd"
elif ("compute" in host) or ( "physics.ucsb.edu" in host) : host = "sb"
else: sys.exit("\033[91mERROR: Unknown host: "+host+" Exit. \033[0m")
print "INFO: Setting up job submission at",('UCSB.' if host=='sb' else 'UCSD.')
hadoop = '/mnt/hadoop/cms'
if host=="sd": hadoop = '/hadoop/cms/phedex'

# Job submission should be done from the babymaker directory, which is under a valid CMSSW release
codedir = os.getcwd()
if not ("/src/babymaker") in codedir:
  print "\033[91mERROR: Please submit from path consistent with: .../src/babymaker/ \033[0m\n"
  sys.exit(0)

# Need a valid proxy to submit condor jobs at UCSD
# or for fallback at UCSB
proxy,valid = "",""
proc = subprocess.Popen('voms-proxy-info', stdout=subprocess.PIPE)
tmp = proc.stdout.read()
if "Proxy not found" in tmp: 
  sys.exit("\033[91mERROR: Proxy not found. \033[0m")
elif ('timeleft' in tmp) and ('0:00:00' in tmp):
  sys.exit("\033[91mERROR: Proxy expired. \033[0m")
else: 
  for info in tmp.splitlines():
    if ("/tmp/x509" in info): proxy = "/tmp/x509"+(string.split(info,"/tmp/x509"))[1]
    if ("timeleft" in info): valid = "Time left before proxy expires: "+info.split()[-1]

print "INFO: Found proxy path",proxy
print "INFO:",valid

# Default output directory is the "out" sub-directory of the current working directory.
outdir = os.getcwd()+'/out/'
sub_time = time.strftime("%y%m%d_%H%M%S", time.gmtime())
if not (os.path.exists(os.path.realpath(outdir))):
    sys.exit("\033[91m ERROR: Directory "+outdir+" does not exist. Please either create a sym link or dir. \033[0m")
outdir = os.path.join(os.path.realpath(outdir),sub_time)
os.mkdir(outdir)

logdir = os.path.join(os.getcwd(),"logs", sub_time)
if not os.path.exists(logdir):
  os.makedirs(logdir)
print "INFO: Babies will be written to: ", outdir
print "INFO: Logs will be written to:   ", logdir

# this is where the condor submission script and job executable are stored
if not (os.path.exists(os.getcwd()+'/run')):
        os.mkdir(os.getcwd()+'/run')
rundir = os.path.join(os.path.realpath(os.getcwd()+'/run'),sub_time)
os.mkdir(rundir)
print rundir
# read in datasets to run over based on the flist_*txt files
# where to find the flists
flistdir = os.path.join(os.getenv("CMSSW_BASE"),"src/flists/")
if not os.path.exists(flistdir):
  sys.exit("ERROR: flists repository not found.")

files_dict = {}
nent_dict = {}
flists_pd = glob.glob(os.path.join(flistdir,"flist*.txt"))
for fnm in flists_pd:
  if any(wish in fnm for wish in mc_wishlist):
    dsname = ''
    if any(ikey in fnm for ikey in comb_keys):
      for ikey in comb_keys: 
        if ikey in fnm:
          dsname = string.split(string.split(fnm,"flist_").pop(),ikey)[0] + ikey
          break
    else:
      sys.exit("ERROR: None of the combination keys (%s) were found in this flist:%s\n" % (comb_keys,fnm))

    print "INFO: Adding PD: ",fnm.replace("flist_","").replace(".txt","")

    if dsname not in files_dict.keys():
      nent_dict[dsname] = 0
      files_dict[dsname] = []
    with open(fnm) as f: 
      for line in f:
        if ("nEventsTotal" in line): # this is read instead of calculated to make resubmission simpler
          nent_dict[dsname] = nent_dict[dsname] + int(line.split().pop())
        if "/store" not in line: continue
        col = line.split()
        files_dict[dsname].append(col[2])

# form new datasets from the data split into subperiods
for pd in data_wishlist:
  # book the dataset names for all sub-periods in advance
  for json in jsonlist:
    dsname = pd + json.replace('data/json/subgolden','').replace('.json','')
    files_dict[dsname] = []
    nent_dict[dsname] = 0 # not filled for data
  # read flists
  flists_pd = glob.glob(os.path.join(flistdir,"flist_"+pd+"_Run2015D*.txt"))
  for fnm in flists_pd:
    with open(fnm) as f: 
      for line in f:
        if "/store" not in line: continue
        col = line.split()
        runlist = [int(irun) for irun in string.split(col[3],",")]
        for run in runlist:
          for jsonfile in goldruns.keys():
            if run in goldruns[jsonfile]:
              dsname = pd + jsonfile.replace('data/json/subgolden','').replace('.json','')
              if (col[2] not in files_dict[dsname]): # don't add same file twice if it has two runs in this subperiod
                files_dict[dsname].append(col[2])

# If on UCSD prep also tarball
if (host=="sd"): 
  print "INFO: Creating babymaker tarball to transfer to work node..."
  os.system("tar --directory=../ --exclude=\"babymaker/out\" --exclude=\"babymaker/run\" --exclude=\"babymaker/logs\" --exclude=\"bmaker/interface/release.hh\" -c babymaker | xz > ../babymaker.tar.xz")

total_jobs = 0
for ids, ds in enumerate(sorted(files_dict.keys())):
  if (maxds!=-1 and ids>=maxds): break

  #release
  cmssw = "CMSSW_7_4_6_patch6"
  if ("Run2015" in ds) or ("RunIISpring15MiniAODv2" in ds): cmssw = "CMSSW_7_4_14"

  # with the list of files in hand, determine the number of condor jobs     
  nfiles = len(files_dict[ds])
  njobs = maxjobs
  if (maxjobs==-1): njobs = (nfiles/maxfiles) if (nfiles%maxfiles==0) else (nfiles/maxfiles+1)
  
  for job in range(0,njobs):
    # name the baby
    bname = "_".join(["baby",ds,"mf"+str(maxfiles),"batch"+str(job)])
    print("INF0: "+bname)

    # check if job had already succeeded on previous submission
    outpath = os.path.join(outdir,bname+".root")
    if os.path.exists(outpath):
      print "\033[38m WARNING: "+outpath+" already exists. Skip job submission. \033[0m"
      continue

    # list of arguments to the cmsRun job
    condor_args = []
    condor_args.append("nEvents="+str(maxevents_perjob))
    condor_args.append("nEventsSample="+str(nent_dict[ds]))
    condor_args.append("inputFiles=\\\n"+",\\\n".join(files_dict[ds][(job*maxfiles):((job+1)*maxfiles)]))
    if (host=="sb"): condor_args.append("outputFile="+outpath)
    else: condor_args.append("outputFile="+bname+".root")
    condor_args.append("condorSubTime="+sub_time)
    if ("Run2015D" in ds):
      json_name = "data/json/subgolden_Run2015D" + ds.split("Run2015D").pop() + ".json"
      if (json_name not in jsonlist): sys.exit("ERROR: Could not find json!")
      condor_args.append("json=babymaker/"+json_name)

    # Create executable that will be transfered to the work node by condor
    exefile =rundir+"/"+bname+".sh"
    fexe = open(exefile,"w")
    if (host=="sb"):
      fexe.write("#! /bin/bash\n")
      fexe.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
      fexe.write("cd "+codedir+"\n")
      fexe.write("eval `scramv1 runtime -sh`\n")
      fexe.write("export ORIGIN_USER="+os.getenv("USER")+"\n")
      fexe.write("cmsRun bmaker/python/bmaker_basic_cfg.py \\\n"+" \\\n".join(condor_args)+"\n")
    else: 
      fexe.write("#! /bin/bash\n")
      fexe.write("source /code/osgcode/cmssoft/cmsset_default.sh\n")
      fexe.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
      fexe.write("eval `scramv1 project CMSSW "+cmssw+"`\n")
      fexe.write("cd "+cmssw+"/src\n")
      fexe.write("eval `scramv1 runtime -sh`\n")
      fexe.write("export ORIGIN_USER="+os.getenv("USER")+"\n")
      fexe.write("tar -xf ../../babymaker.tar.xz\n")
      fexe.write("cd babymaker\n")
      fexe.write("./compile.sh\n")
      fexe.write("cmsRun bmaker/python/bmaker_basic_cfg.py \\\n"+" \\\n".join(condor_args)+"\n")
      fexe.write("echo \"cmsRun exit code \"$?\n")
      #fexe.write("lcg-cp -b -D srmv2 --vo cms -t 2400 --verbose file:"+bname+".root srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN="+outpath+"\n")
      if "T1tttt" in bname:
        fexe.write("./bmaker/genfiles/run/skim_scan_onefile.exe "+bname+".root\n")
      fexe.write("for i in $(ls *.root); do\n")
      fexe.write("\tlcg-cp -b -D srmv2 --vo cms -t 2400 --verbose file:$i srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN="+outdir+"/$i\n")
      fexe.write("done\n")
      fexe.write("cd ../../..\n")
      fexe.write("rm -rf "+cmssw+"\n")
    fexe.close()
    os.system("chmod u+x "+exefile)

    # Create condor submission cmd file
    cmdfile = rundir+"/"+bname+".cmd"
    print "cmdfile is "+ cmdfile
    fcmd = open(cmdfile,"w")
    if (host=="sb"):
      fcmd.write("Executable   = "+exefile+"\n")
      fcmd.write("Universe     = vanilla\n")
      # send proxy even for local submissions
      # in case fallback is necessary
      fcmd.write("use_x509userproxy = True\n")
      fcmd.write("x509userproxy="+proxy+"\n")
      fcmd.write("Log          = "+logdir+ "/"+bname+".log\n")
      fcmd.write("output       = "+logdir+"/"+bname+".out\n")
      fcmd.write("error        = "+logdir+"/"+bname+".err\n")
      fcmd.write("Notification = never\n")
      fcmd.write("Queue\n")
    else:
      fcmd.write("Universe = grid\n")
      fcmd.write("Grid_Resource = condor cmssubmit-r1.t2.ucsd.edu glidein-collector.t2.ucsd.edu\n")
      fcmd.write("use_x509userproxy = True\n")
      fcmd.write("x509userproxy="+proxy+"\n")
      fcmd.write("+remote_DESIRED_Sites=\""+whitelist+"\"\n")
      fcmd.write("Executable = "+exefile+"\n")
      fcmd.write("Transfer_Executable = True\n")
      fcmd.write("should_transfer_files = YES\n")
      fcmd.write("transfer_input_files = ../babymaker.tar.xz\n")
      fcmd.write("Notification = Never\n")
      fcmd.write("Log          = "+logdir+"/"+bname+".log\n")
      fcmd.write("output       = "+logdir+"/"+bname+".out\n")
      fcmd.write("error        = "+logdir+"/"+bname+".err\n")
      fcmd.write("queue 1\n")
    fcmd.close()
    total_jobs = total_jobs + 1

# Submit condor job
if host=="sb":
  cmd = "ssh cms25.physics.ucsb.edu condor_submit "
else:
  cmd = "condor_submit "
  print "INFO: Submitting", cmdfile

# for the sake of efficiency, submit all jobs at once
if host=="sb":
  os.system("scp " + proxy + " cms25.physics.ucsb.edu:/tmp")
os.system("cat " + rundir + "/baby*.cmd > " + rundir + "/submit_all.cmd")
os.system(cmd + rundir + "/submit_all.cmd")
print "Submitted ", total_jobs, "jobs"
