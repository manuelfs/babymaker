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

## Signal scan, 81M events
mc_wishlist.append("SMS-T1tttt_mGluino")

## TTJets, 170M events
mc_wishlist.append("TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")

## QCD, 85M events 
mc_wishlist.append("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")

## Other MC, 95M events 
mc_wishlist.append("DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WWTo2L2Nu_13TeV-powheg_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WWToLNuQQ_13TeV-powheg_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9")
mc_wishlist.append("ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9_ext3")

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
  tarcmd = "tar --directory=../ --exclude=\"out\" --exclude=\"run\""
  tarcmd += " --exclude=\"logs\" --exclude=\"bmaker/interface/release.hh\""
  tarcmd += " --exclude=\"data/csc_beamhalo_filter/*\""
  tarcmd += " --exclude=\".git\""
  tarcmd += " -c babymaker | xz > ../babymaker.tar.xz"
  os.system(tarcmd)

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
