#!/usr/bin/env python
import os, sys, re
import glob
import string
from pprint import pprint
import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p","--logpath")
parser.add_argument("-t","--timestamp")
args = parser.parse_args()

if (args.timestamp):
  timestamp = args.timestamp
elif (args.logpath):
  timestamp = args.logpath.rstrip("/").split("/")[-1]
else:
  sys.exit("Please provide a timestamp either as, e.g. 151019_011440, or a path ending with the timestamp")

onePerJob = False

# determine if job is running at UCSB or UCSD
host = os.environ.get("HOSTNAME")
atUCSB = False
if "compute" in host or "physics.ucsb.edu" in host:
  atUCSB = True

# redirector = "root://cmsxrootd.fnal.gov//"
redirector = "root://cms-xrd-global.cern.ch//"

bdir = os.getcwd()
if ('babymaker' not in bdir.split("/").pop()):
  sys.exit("Execute from babymaker directory")

logdir = os.path.join(bdir,'logs',timestamp)
if not os.path.exists(logdir):
  sys.exit("Can't find log directory %s" %logdir)

arxivdir = os.path.join(bdir,'logs',timestamp,'arxiv')
if not os.path.exists(arxivdir):
  os.mkdir(arxivdir)

#loglist = [x for x in glob.glob(logdir+"/*.log") if "_rs" in x]
loglist = [x for x in glob.glob(logdir+"/*.log")] 
print "Found %i logs" %(len(loglist))

failed = set()
unfinished = set()
xrootd_err = set()

for flog in loglist:
  ferr = flog.rstrip(".log") + ".err"
  fout = flog.rstrip(".log") + ".out"
  bname = flog.split("/").pop().rstrip(".log")
  logfile = open(flog).read()
  if "Job was aborted by the user" in logfile:
    failed.add(bname)
    continue
  if os.path.getsize(fout)==0:
    # failed.add(bname)
    unfinished.add(bname)
  else:
    if "BABYMAKER: Written" not in open(fout).read():
      failed.add(bname)
    errfile = open(ferr).read()
    # transfer not necessary at UCSB
    if "Transfer took" not in errfile and not atUCSB:
      failed.add(bname)
    if "cmsRun exit code 1" in errfile:
      failed.add(bname)
    if "Fatal Exception" in errfile:
      failed.add(bname)
    if "Socket error while handshaking: [FATAL] Auth failed" in errfile:
      xrootd_err.add(bname)

    if bname in xrootd_err and bname not in failed:
      print "xrootd err but success(?): ",bname 

if len(unfinished) > 0 :
  print "--------- Unfinished:"
  pprint(unfinished)
  print "--------- Total unfinished ",len(unfinished),"\n"
if len(failed) > 0 :
  print "--------- Failed:"
  pprint(failed)
  print "Total unfinished ",len(unfinished)
  print "Total failed ",len(failed)
  print "Total with xrootd err",len(xrootd_err)
else :
  if len(unfinished) == 0 : sys.exit("\nCongrats, no jobs failed. You might be able to go out and enjoy the mountains now :o)\n")
  else : sys.exit("\nCongrats, no jobs failed, but still "+str(len(unfinished))+" jobs to go\n")

user_input = raw_input('Resubmit jobs [y/N]?')
if (user_input!='y'):
  sys.exit("Bye.")
else:
  user_input = raw_input('Resubmit with one file per job [y/N]?')
  if (user_input=='y'):
    onePerJob = True

# --- resubmission
total_jobs = 0
for old_baby in failed:
  fexe = os.path.join(logdir.replace("/logs/","/run/"), old_baby+".sh")
  fcmd = os.path.join(logdir.replace("/logs/","/run/"), old_baby+".cmd")
  os.rename(logdir+"/"+old_baby+".log", arxivdir+"/"+old_baby+".log")
  os.rename(logdir+"/"+old_baby+".err", arxivdir+"/"+old_baby+".err")
  os.rename(logdir+"/"+old_baby+".out", arxivdir+"/"+old_baby+".out")
  if not onePerJob:
    old_exe = open(fexe).read()
    new_exe = old_exe.replace("file:/hadoop/cms/phedex", redirector)
    with open(fexe,'w') as f: f.write(new_exe)
    sys_cmd = "condor_submit " + fcmd
    if atUCSB: sys_cmd = "ssh cms25.physics.ucsb.edu condor_submit " + fcmd
    print "INFO: Submitting", fcmd
    os.system(sys_cmd)
    total_jobs = total_jobs + 1
  else:
    tags = old_baby.split("_")
    batch = int(tags[-1].strip("batch"))
    nfiles = int(tags[-2].strip("mf"))

    # --- read the old submission files
    if (not os.path.exists(fexe)) or (not os.path.exists(fcmd)): 
      sys.exit("Cannot find either .sh or .cmd for: %s" % logdir)
    old_exe = open(fexe).readlines()
    old_cmd = open(fcmd).read()

    # --- parse for input files
    inputfiles = [line for line in old_exe if ("/store/" in line and "lcg-cp" not in line)]
    for i, line in enumerate(inputfiles):
      inputfiles[i] = "/store/"+line.split("/store/").pop().split(".root")[0]+".root"
    # pprint(inputfiles)

    # --- assume none are local in case job failed because file was not found
    for i, infile in enumerate(inputfiles):
      new_baby = old_baby + "_rs"+str(i)
      new_exe_lines = [line for line in old_exe if ("/store/" not in line or "lcg-cp" in line)]
      ind = new_exe_lines.index('inputFiles=\\\n')
      new_exe_lines.insert(ind+1, redirector + infile + ' \\\n')
      new_exe = ''.join(new_exe_lines)
      new_exe = new_exe.replace(old_baby, new_baby)
      with open(fexe.replace(old_baby,new_baby),'w') as f1: f1.write(new_exe)
      new_cmd = old_cmd.replace(old_baby, new_baby)
      fnew_cmd = fcmd.replace(old_baby,new_baby)
      with open(fnew_cmd,'w') as f2: f2.write(new_cmd)
      sys_cmd = "condor_submit " + fnew_cmd
      if atUCSB: sys_cmd = "ssh cms25.physics.ucsb.edu condor_submit " + fnew_cmd
      print "INFO: Submitting", fnew_cmd
      os.system(sys_cmd)
      total_jobs = total_jobs + 1

print("Submitted %i jobs." % total_jobs)

# --- check if output files already exist
for old_baby in failed:
  fexe = os.path.join(logdir.replace("/logs/","/run/"), old_baby+".sh")
  old_exe = open(fexe).readlines()
  for line in old_exe:
    if ("SFN=") in line:
      if ("T1tttt" in old_baby):
        outputfile = line.split("SFN=").pop().strip("$i\n") + old_baby + ".root"
      else: 
        outputfile = line.split("SFN=").pop().strip("\n")
      if os.path.exists(outputfile):      
        user_input = raw_input('Remove output file %s corresponding to a failed job [Y/n]?' % outputfile)
        if (user_input!='n'):
          print "Removing file: ", outputfile
          os.remove(outputfile)

          
