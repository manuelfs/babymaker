import os, sys 
import time
import glob
import string

# Maximum number of files per condor job
maxfiles = 1

# Job submission should be done from the babymaker directory, which is under a valid CMSSW release
codedir = os.getcwd()
if not 'CMSSW_7_4_6_patch6/src/babymaker' in codedir:
    print "\033[91m ERROR: Please submit from path consistent with: <basedir>/CMSSW_7_4_6_patch6/src/babymaker/ \033[0m"
    sys.exit(0)

# Default output directory is the "out" sub-directory of the current working directory. It is created if it does not exist
outdir = os.getcwd()+'/out'
if not (os.path.exists(outdir)):
    os.mkdir(outdir)

# Get list of datasets and corresponding list of files
flists = glob.glob("sub/flist*.txt")
filedict = {}
for fnm in flists:
    dsname = string.replace(string.replace(fnm,"sub/flist_",''),".txt",'')
    with open(fnm) as f:
        filedict[dsname] = f.read().splitlines()

for ds in filedict.keys():

    # with the list of files in hand, determine the number of condor jobs     
    nfiles = len(filedict[ds])
    njobs = (nfiles/maxfiles) if (nfiles%maxfiles==0) else (nfiles/maxfiles+1)
    if njobs>2: njobs=2

    for job in range(0,njobs):
        # name the baby
        bname = '_'.join(['baby',ds,'batch'+str(job)])

        # check if job had already succeeded on previous submission
        outpath = os.path.join(outdir,bname+'.root')
        if os.path.exists(outpath):
            print "\033[93m WARNING: "+outpath+" already exists. Skip job submission. \033[0m"
            continue

        # list of arguments to the cmsRun job
        args = []
        args.append('inputFiles=\\\n'+",\\\n".join(filedict[ds][(job*maxfiles):((job+1)*maxfiles)]))
        args.append('outputFile='+outpath)

        # Create executable that will be transfered to the work node by condor
        exefile ='sub/'+bname+".sh"
        fexe = open(exefile,"w")
        fexe.write("#! /bin/bash\n")
        fexe.write(". /cvmfs/cms.cern.ch/cmsset_default.sh\n")
        fexe.write("cd "+codedir+"\n")
        fexe.write("eval `scramv1 runtime -sh`\n")
        fexe.write("echo $CMSSW_BASE\n")
        fexe.write("ls -tr\n")
        fexe.write("cmsRun bmaker/python/bmaker_basic_cfg.py \\\n"+' \\\n'.join(args)+"\n")
        fexe.close()
        os.system("chmod u+x "+exefile)

        # Create condor submission cmd file
        cmdfile = 'sub/'+bname+'.cmd'
        fcmd = open(cmdfile,"w")
        fcmd.write("Executable   = "+exefile+"\n")
        fcmd.write("Universe     = vanilla\n")
        fcmd.write("Log          = sub/"+bname+".log\n")
        fcmd.write("output       = sub/"+bname+".out\n")
        fcmd.write("error        = sub/"+bname+".err\n")
        fcmd.write("Notification = never\n")
        fcmd.write("Queue\n")
        fcmd.close()

        # Submit condor job
        cmd = "condor_submit " + cmdfile
        print cmd
        os.system(cmd)
