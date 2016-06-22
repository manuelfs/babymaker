#!/usr/bin/env python

###### Script that links missing files an "old" production to a "new" production
import os, sys, subprocess
import pprint
import glob

# Setting folders
oldfolder    = '/net/cms2/cms2r0/babymaker/babies/2016_06_14/data/unskimmed/'
newfolder    = '/net/cms2/cms2r0/babymaker/babies/2016_06_21/data/unskimmed/'

oldfiles = glob.glob(oldfolder+'/*.root')
for file in oldfiles:
  filename = file.split('/')[-1]
  filename = filename.replace(".root", "_link.root")
  os.symlink(file, newfolder+filename)

print "\nWritten links to files in "+oldfolder+" into "+newfolder+"\n"
