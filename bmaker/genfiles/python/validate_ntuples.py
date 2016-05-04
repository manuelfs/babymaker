#!/usr/bin/env python

###### Script that compares the yields in old and new ntuples
from ROOT import TChain, TH1D
import os, sys, subprocess
import pprint
import glob
import json
import string
import time
import math

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Setting folders
oldfolder    = '/net/cms2/cms2r0/babymaker/babies/2016_01_11/mc/T1tttt/skim_abcd/'
newfolder    = '/net/cms26/cms26r0/babymaker/2016_04_29/normalized/T1tttt/skim_abcd/'

oldfolder    = '/net/cms2/cms2r0/babymaker/babies/2016_01_11/mc/T1tttt/'
newfolder    = '/net/cms26/cms26r0/babymaker/2016_04_29/normalized/T1tttt/'

newfolder    = '/net/cms27/cms27r0/babymaker/2016_04_29/mc/skim_met100nb2nj4nl0'
oldfolder    = '/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_met100nb2nj4nl0'

oldfolder    = '/net/cms2/cms2r0/babymaker/babies/2015_11_28/mc/skim_1lht500met200'
newfolder    = '/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_1lht500met200'

## Finding tags for each dataset
infiles = set()
for file in glob.glob(newfolder+'/*.root'):
  tag = file.split("RunII")[0]
  if ("TTJets_Tune" not in tag) and ("DYJetsToLL_M-50_Tune" not in tag): tag = tag.split("Tune")[0]
  tag = tag.split("13TeV")[0]
  tag = tag.split("pythia")[0]
  tag = tag.split("baby_")[1]
  tag = tag.split("__")[0]
  if tag[0] != '_': tag = "_"+tag
  if tag[-1] != '_': tag = tag+"_"
  infiles.add(tag)
sortedfiles = list()
for file in infiles:
    sortedfiles.append(file)
sortedfiles = sorted(sortedfiles)

print '\nOLD FOLDER: '+oldfolder
print 'NEW FOLDER: '+newfolder
print '\n{:>40}'.format(' Ntuple name          ')+'{:>16}'.format('Difference')+'{:>17}'.format('Old yield'),
print '{:>17}'.format('New yield')+'{:>17}'.format('Old entries')+'{:>17}'.format('New entries')
print '=' * 128
not_in_old = list()
not_in_new = list()
rows = list()
line = 1
histo = TH1D("histo","",10,0,10)
for ifile in sortedfiles:
  ochain = TChain("tree")
  oldntuples = oldfolder+"/*"+ifile+'*root'
  no = ochain.Add(oldntuples)
  if no == 0:
    not_in_old.append(ifile)
    continue;
  nchain = TChain("tree")
  newntuples = newfolder+"/*"+ifile+'*root'
  nn = nchain.Add(newntuples)
  if nn == 0:
    not_in_new.append(ifile)
    continue;

  no = ochain.Draw("1>>histo","weight","goff")
  oldtot = histo.Integral()
  nn = nchain.Draw("1>>histo","weight","goff")
  newtot = histo.Integral()
    
  if oldtot != 0 : diff = (newtot-oldtot)*100/oldtot
  elif newtot == 0 : diff = 0
  else :  diff = 999

  pretag = ""
  posttag = ""
  ## Appending rows with significant differences for later printing
  if abs(diff) > 150/math.sqrt(no+1) and abs(diff) > 150/math.sqrt(nn+1): 
    rows.append([ifile, abs(diff), oldtot, newtot, no, nn])
    pretag = bcolors.FAIL
    posttag = bcolors.ENDC
  ## Printing all rows
  print pretag+'{:>40}'.format(ifile)+'{:>14.2f}'.format(diff)+' %'+'{:>17.2f}'.format(oldtot),
  print '{:>17.2f}'.format(newtot)+'{:>17}'.format(no)+'{:>17}'.format(nn)+posttag
  if line == 5 : 
    print
    line = 0
  line += 1

## Sorting rows by difference
if len(rows) > 0: 
  print bcolors.FAIL + "\nSamples off by more than 1.5 sigma"+ bcolors.ENDC
  print '\n{:>40}'.format(' Ntuple name          ')+'{:>16}'.format('Difference')+'{:>17}'.format('Old yield'),
  print '{:>17}'.format('New yield')+'{:>17}'.format('Old entries')+'{:>17}'.format('New entries')
  print '=' * 128
rows = sorted(rows, key=lambda rows: rows[1])
line = 1
for row in rows:
  ifile = row[0]
  diff = row[1]
  oldtot = row[2]
  newtot = row[3]
  no = row[4]
  nn = row[5]

  ## Printing rows with significant differences
  print '{:>40}'.format(ifile)+'{:>14.2f}'.format(diff)+' %'+'{:>17.2f}'.format(oldtot),
  print '{:>17.2f}'.format(newtot)+'{:>17}'.format(no)+'{:>17}'.format(nn)
  if line == 5 : 
    print
    line = 0
    line += 1

if len(not_in_old) > 0:
  print bcolors.BOLD + '\nNtuples not found in '+oldfolder+':'+ bcolors.ENDC
  for ntu in not_in_old:
    print '\t'+ntu
if len(not_in_new) > 0:
  print bcolors.BOLD + '\nNtuples not found in '+newfolder+':'+ bcolors.ENDC
  for ntu in not_in_new:
    print '\t'+ntu

print
sys.exit(0)
