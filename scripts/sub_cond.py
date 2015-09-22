#!/usr/bin/env python

import os, sys, subprocess
import glob
import string
import time

#What to submit? Use substrings that would be found in the desired dataset
# no mid-word wild cards enabled yet, 
# e.g. if we want only 25ns TTJets, use a substring that contains it all:
# "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns"
wishlist = []
wishlist.append("TTJets")

print 

# Maximum number of input MINIAOD files per condor job
maxfiles = 20

cmssw = "CMSSW_7_4_6_patch6"

# Only matters if running on UCSD:
# To run on multiple T2's use, e.g:
# whitelist = "T2_US_UCSD,T2_US_Wisconsin,T2_US_Nebraska,T2_US_Caltech"
# To run only at UCSD use:
whitelist = "T2_US_UCSD"

# Condor set up depends on whether we run on UCSB or UCSD
host = os.environ.get("HOSTNAME")
if "ucsd" in host: host = "sd"
elif host=="cms0.physics.ucsb.edu": host = "sb"
elif ("ucsb" in host) or ("compute" in host): sys.exit("\033[91mERROR: Submission must be done from cms0. \033[0m")
else: sys.exit("\033[91mERROR: Unknown host: "+host+" Exit. \033[0m")
print "INFO: Setting up job submission at",('UCSB.' if host=='sb' else 'UCSD.')

# Job submission should be done from the babymaker directory, which is under a valid CMSSW release
codedir = os.getcwd()
if not (cmssw+"/src/babymaker") in codedir:
    print "\033[91mERROR: Please submit from path consistent with: <basedir>/"+cmssw+"/src/babymaker/ \033[0m\n"
    sys.exit(0)

# Default output directory is the "out" sub-directory of the current working directory.
outdir = os.getcwd()+'/out/'
sub_time = time.strftime("%y%m%d_%H%M%S", time.gmtime())
if not (os.path.exists(os.path.realpath(outdir))):
    sys.exit("\033[91m ERROR: Directory "+outdir+" does not exist. Please either create a sym link or dir. \033[0m")
outdir = os.path.join(os.path.realpath(outdir),sub_time)
os.mkdir(outdir)
logdir = os.path.join(outdir,"logs")
if not os.path.exists(logdir):
    if (host=='sd') or (host=='sb' and ('hadoop' not in logdir)): 
        os.makedirs(logdir)
    else: 
        print "\033[91m ERROR: Cannot create log directory on hadoop at UCSB. Using ./logs/ instead. \033[0m" 
        logdir = './logs/'
        if not os.path.exists(logdir): os.mkdir(logdir)
print "INFO: Babies will be written to: ", outdir
print "INFO: Logs will be written to:   ", logdir

# read in datasets to run over based on the flist_*txt files present in the run directory
# from there, read the filenames for each dataset and add up the number of events
# firstly, bookkeep datasets that are not extensions
flistdir = "run/" 
# regardless of host, create the ./run/ directory if it doesn't exist since this is where auto-generated submission scripts are written
if not (os.path.exists(os.getcwd()+'/'+flistdir)):
        os.mkdir(os.getcwd()+'/'+flistdir)
# the actual flist is stored centrally if running at UCSB and in ./run/ at UCSD
if host=="sb": flistdir = "/net/cms2/cms2r0/babymaker/flist/"

flists_pd = [x for x in glob.glob(os.path.join(flistdir,"flist*.txt")) if "_ext" not in x]
files_dict = {}
nent_dict = {}
for fnm in flists_pd:
    if any(wish in fnm for wish in wishlist):
        dsname = string.split(string.split(fnm,"flist_").pop(),".txt")[0]
        print "INFO: Adding PD: ",dsname
        with open(fnm) as f: files_dict[dsname] = f.read().splitlines()
        nent_dict[dsname] = int(files_dict[dsname].pop().split().pop())
# now add the extensions...
flists_ext = glob.glob(os.path.join(flistdir,"flist*_ext*.txt"))
for fnm in flists_ext:
    if any(wish in fnm for wish in wishlist):
        dsname = string.split(string.split(fnm,"flist_").pop(),".txt")[0]
        print "INFO: Adding PD extension: ",dsname
        #get name without the "_ext?" to check if we are already bookeeping some files from the original version of that dataset
        split_dsname = string.split(dsname,"_ext")
        orig_dsname = split_dsname[0] + split_dsname[1][1:] # skip one character after "_ext" which is the number of the extension
        if (orig_dsname not in files_dict.keys()): # in case it was an extension of a dataset we are not running on
            with open(fnm) as f: files_dict[orig_dsname] = f.read().splitlines()
            nent_dict[orig_dsname] = int(files_dict[orig_dsname].pop().split().pop())
        else:
            with open(fnm) as f: files_dict[orig_dsname].extend(f.read().splitlines())
            nent_dict[orig_dsname] = nent_dict[orig_dsname] + int(files_dict[orig_dsname].pop().split().pop())
    
# add up number of events in each dataset 
for ds in files_dict.keys():
    for iline in files_dict[ds]:
        files_dict[ds] 

# Need a valid proxy to submit condor jobs at UCSD
proxy,valid = "",""
if host=="sd":
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

# If on UCSD prep also tarball
if (host=="sd"): 
    print "INFO: Creating babymaker tarball to transfer to work node..."
    os.system("tar --directory=../ --exclude=\"babymaker/out\" --exclude=\"babymaker/run\" -c babymaker | xz > ../babymaker.tar.xz")

for ds in files_dict.keys():

    # with the list of files in hand, determine the number of condor jobs     
    nfiles = len(files_dict[ds])
    njobs = (nfiles/maxfiles) if (nfiles%maxfiles==0) else (nfiles/maxfiles+1)

    for job in range(0,njobs):
        # name the baby
        bname = "_".join(["baby",ds,"mf"+str(maxfiles),"batch"+str(job)])

        # check if job had already succeeded on previous submission
        outpath = os.path.join(outdir,bname+".root")
        if os.path.exists(outpath):
            print "\033[38m WARNING: "+outpath+" already exists. Skip job submission. \033[0m"
            continue

        # list of arguments to the cmsRun job
        args = []
        args.append("nEvents=-1")
        args.append("nEventsSample="+str(nent_dict[ds]))
        args.append("inputFiles=\\\n"+",\\\n".join(files_dict[ds][(job*maxfiles):((job+1)*maxfiles)]))
        if (host=="sb"): args.append("outputFile="+outpath)
        else: args.append("outputFile="+bname+".root")

        # Create executable that will be transfered to the work node by condor
        exefile ="run/"+bname+".sh"
        fexe = open(exefile,"w")
        if (host=="sb"):
            fexe.write("#! /bin/bash\n")
            fexe.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            fexe.write("cd "+codedir+"\n")
            fexe.write("eval `scramv1 runtime -sh`\n")
            fexe.write("cmsRun bmaker/python/bmaker_basic_cfg.py \\\n"+" \\\n".join(args)+"\n")
        else: 
            fexe.write("#! /bin/bash\n")
            fexe.write("source /code/osgcode/cmssoft/cmsset_default.sh\n")
            fexe.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
            fexe.write("eval `scramv1 project CMSSW "+cmssw+"`\n")
            fexe.write("cd "+cmssw+"/src\n")
            fexe.write("eval `scramv1 runtime -sh`\n")
            fexe.write("tar -xf ../../babymaker.tar.xz\n")
            fexe.write("cd babymaker\n")
            fexe.write("./compile.sh\n")
            fexe.write("cmsRun bmaker/python/bmaker_basic_cfg.py \\\n"+" \\\n".join(args)+"\n")
            fexe.write("lcg-cp -b -D srmv2 --vo cms -t 2400 --verbose file:"+bname+".root srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN="+outpath+"\n")
            fexe.write("cd ../../..\n")
            fexe.write("rm -rf "+cmssw+"\n")
        fexe.close()
        os.system("chmod u+x "+exefile)

        # Create condor submission cmd file
        cmdfile = "run/"+bname+".cmd"
        fcmd = open(cmdfile,"w")
        if (host=="sb"):
            fcmd.write("Executable   = "+exefile+"\n")
            fcmd.write("Universe     = vanilla\n")
            fcmd.write("Log          = "+logdir+ "/"+bname+".log\n")
            fcmd.write("output       = "+logdir+"/"+bname+".out\n")
            fcmd.write("error        = "+logdir+"/"+bname+".err\n")
            fcmd.write("Notification = never\n")
            fcmd.write("Queue\n")
        else:
            fcmd.write("Universe = grid\n")
            fcmd.write("Grid_Resource = condor cmssubmit-r1.t2.ucsd.edu glidein-collector.t2.ucsd.edu\n")
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

        # Submit condor job
        cmd = "condor_submit " + cmdfile
        print "INFO: Submitting", cmdfile
        os.system(cmd)
