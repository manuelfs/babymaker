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
  sys.exist("Please provide a timestamp either as, e.g. 151019_011440, or a path ending with the timestamp")

onePerJob = False

# determine if job is running at UCSB or UCSD
host = os.environ.get("HOSTNAME")
atUCSB = "compute" in host or "physics.ucsb.edu" in host

redirector = "root://cmsxrootd.fnal.gov//"

bdir = os.getcwd()
if ('babymaker' not in bdir.split("/").pop()):
  sys.exit("Execute from babymaker directory")

logdir = os.path.join(bdir,'logs',timestamp)
if not os.path.exists(logdir):
  sys.exit("Can't find log directory %s" %logdir)

loglist = [x for x in glob.glob(logdir+"/*.log") if "_rs" not in x] 
print "Found %i logs" %(len(loglist))

failed = set()
unfinished = set()

for flog in loglist:
  ferr = flog.rstrip(".log") + ".err"
  fout = flog.rstrip(".log") + ".out"
  bname = flog.split("/").pop().rstrip(".log")
  if os.path.getsize(fout)==0 or os.path.getsize(ferr)==0:
    unfinished.add(bname)
  else:
    if "BABYMAKER: Written" not in open(fout).read():
      failed.add(bname)
    # transfer not necessary at UCSB
    if "Transfer took" not in open(ferr).read() and not atUCSB:
      failed.add(bname)

print "--------- Unfinished:"
pprint(unfinished)
print "--------- Failed:"
pprint(failed)

user_input = raw_input('Resubmit jobs [y/n]?')
if (user_input!='y'):
  sys.exit("Bye.")
else:
  user_input = raw_input('Resubmit with one file per job [y/n]?')
  if (user_input=='y'):
    onePerJob = True

# --- check if corrupted output files already exist
badfiles = []
for old_baby in failed:
  fexe = os.path.join(logdir.replace("/logs/","/run/"), old_baby+".sh")
  old_exe = open(fexe).readlines()
  for line in old_exe:
    if ("SFN=") in line:
      outputfile = line.split("SFN=").pop().strip("\n")
      if os.path.exists(outputfile):
        print "Output file exists %s" % outputfile
        badfiles.append(outputfile)
      break

# --- resubmission

sys_cmd = ""
# at UCSB submission is only possible from cms25
if atUCSB:
  sys_cmd+="ssh cms25.physics.ucsb.edu "
total_jobs = 0
for old_baby in failed:
  fexe = os.path.join(logdir.replace("/logs/","/run/"), old_baby+".sh")
  fcmd = os.path.join(logdir.replace("/logs/","/run/"), old_baby+".cmd")

  if not onePerJob:
    old_exe = open(fexe).read()
    new_exe = old_exe.replace("file:/hadoop/cms/phedex", redirector)
    with open(fexe,'w') as f: f.write(new_exe)
    sys_cmd += "condor_submit " + fcmd
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
    if len(inputfiles)!=nfiles:
      sys.exit("Could not parse for the input files")
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
      sys_cmd += "condor_submit " + fnew_cmd
      print "INFO: Submitting", fnew_cmd
      os.system(sys_cmd)
      total_jobs = total_jobs + 1

print("Submitted %i jobs." % total_jobs)

# --- deal with old output
print("Corrupt output files:")
pprint(badfiles)
for badfile in badfiles:
  froot = ROOT.TFile(badfile)
  if froot.IsZombie():
    user_input = raw_input('Remove zombie file %s [y/n]?' % badfile)
    if (user_input=='y'):
      os.remove(badfile)
  else: 
    print "**** ERROR: Inputs resubmitted, but output exist. Please check status and remove manually: ", badfile

          
