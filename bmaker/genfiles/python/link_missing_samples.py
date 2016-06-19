#!/usr/bin/env python

###### Script that links missing files an "old" production to a "new" production
import os, sys, subprocess
import pprint
import glob

# Setting folders
oldfolder    = '/net/cms27/cms27r0/babymaker/2016_04_29/mc/unskimmed/'
newfolder    = '/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed/'
link_dest    = '/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed/'

## Finding tags for each dataset
oldsamples, newsamples = set(), set()
oldfiles = glob.glob(oldfolder+'/*.root')
newfiles = glob.glob(newfolder+'/*.root')
for file in oldfiles+newfiles:
  tag = file.split("RunII")[0]
  if "SMS-" in tag: continue
  if ("TTJets_Tune" not in tag) and ("DYJetsToLL_M-50_Tune" not in tag): tag = tag.split("Tune")[0]
  tag = tag.split("13TeV")[0]
  tag = tag.split("pythia")[0]
  tag = tag.split("baby_")[1]
  tag = tag.split("__")[0]
  if tag[0] != '_': tag = "_"+tag
  if tag[-1] != '_': tag = tag+"_"
  if oldfolder in file: oldsamples.add(tag)
  else: newsamples.add(tag)


missing = sorted(list(oldsamples-newsamples))
if len(missing)>0:
    print "Samples in the old folder that are not found in the new folder:"
    for itag in missing:
        print "\t",itag
        for file in oldfiles:
            if itag in file: 
                print "Linking:",file
                os.symlink(file, link_dest+"/"+file.split("/")[-1])
else:
    print "No missing samples in the new production"

extras = sorted(list(newsamples-oldsamples))
if len(extras)>0:
    print "Samples in the new folder that are not found in the old folder:"
    for i in extras:
        print "\t",i
else:
    print "No extra samples in the new production"

