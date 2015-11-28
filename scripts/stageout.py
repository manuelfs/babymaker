#!/usr/bin/env python 
import os, glob

# ------------------------------------------------------------------------------------------
# This script can be used to "manually" download babies that were produced ok, but failed 
# to actually be transferred to our Tier 3. It *assumes* that all the failed transfers
# are due to one bad site. But it is ok to use it, even if assumption is wrong... it just
# will not find the files and you will get a message "FAIL" in red.
# If we see that there are often >1 bad sites, this can be modified.
# It seems to take ~ 30-50s per file download
#
# 1. from the web task monitoring interface open one log file of the type "Job" 
#    and find the origin PFN, listed on the line starting with "Stage Out Successful"
#    enter the result here as origin, stripping anything including and after the dataset name
#
# 2. in the same file, find the line starting with: "JOB AD: CRAB_Destination"
#    and check the destination header, likely as the example, but with your CERN username
#
# 3. again looking at the web monitoring interface,
#    get the failed task names (without the "$USER_") and the number of jobs 
#    (regardless of status) enter this in the dictionary named 'tasks' below, enter 
# ------------------------------------------------------------------------------------------

# enter origin and destination headers in place of these examples
origin = 'srm://heplnx204.pp.rl.ac.uk:8443/srm/managerv2?SFN=/pnfs/pp.rl.ac.uk/data/cms/store/temp/user/ana.91d36fddf73016fe56a4674b87cd86f61feae489/'
dest = 'srm://cms25.physics.ucsb.edu:8443/srm/v2/server?SFN=/mnt/hadoop/cms/store/user/ana/'

# failed tasks and the total number of jobs
tasks = {
'crab_SMS-T1tttt_mGluino-600_mLSP-250to325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15FSPremix-MCRUN2_74_V9-v1__MINIAODSIM': 11,
'crab_SMS-T1tttt_mGluino-1150to1175_mLSP-750to925_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15FSPremix-MCRUN2_74_V9-v1__MINIAODSIM': 10,
'crab_SMS-T1tttt_mGluino-1400_mLSP-1to1175_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15FSPremix-MCRUN2_74_V9-v1__MINIAODSIM': 11,
'crab_SMS-T1tttt_mGluino-1425to1450_mLSP-1to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15FSPremix-MCRUN2_74_V9-v1__MINIAODSIM': 11,
'crab_SMS-T1tttt_mGluino-1800to1850_mLSP-1to1450_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15FSPremix-MCRUN2_74_V9-v1__MINIAODSIM': 13,
'crab_TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM': 25,
'crab_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM': 70
}

# In general there should be no need to modify anything below...

for task in tasks.keys():
  print '\n' + 20*'=--' + '\n'
  print "Working on task:", task
  print "Number of jobs:", tasks[task]
  # find what is already on hadoop
  hadoop_base = dest.split('SFN=').pop()
  good_files = glob.glob(hadoop_base + task.replace('crab_','').split('__')[0]+"/" + task + "/*/*/*root")
  # see if there were multiple subdirectories in the time-stamped directory
  runs = set([i.split('/')[-2] for i in good_files]) # put in a set to remove repeating entries
  if (len(runs)>1): print "Found multiple runs ", len(runs)

  # print info on files on hadoop
  if (len(good_files)>0):
    print "Found %i files:" % len(good_files)
    # for i in good_files: print i
  else:
    print "No files found on hadoop"
    continue

  # download files
  ndls = 0
  for run in runs:
    for job in range(1,tasks[task]+1):
      abspath = good_files[0][0:good_files[0].rfind('_')] + '_'+str(job)+'.root'
      relpath = abspath.replace(hadoop_base,'')
      # do we already have the file?
      if abspath in good_files: continue

      cmd = 'lcg-cp --verbose --vo=cms -b -D srmv2 '+ origin + relpath + ' ' + dest + relpath
      exitcode = os.system(cmd)
      if (exitcode==0): ndls = ndls + 1
      else: print "ERROR:: Download failed."

  if (tasks[task] == (len(good_files) + ndls)): 
    print "\033[92m -- SUCCESS:\033[0m", task
  else: 
    print "\033[91m -- FAIL\033[0m", task

