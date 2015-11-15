#!/usr/bin/env python

import os, sys 
import json
import glob
import string
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-j","--jsonpath")
args = parser.parse_args()

pb_per_json = 100.
period_name = 'Run2015D'

jdata = {}
with open(args.jsonpath) as jfile:
  jdata = json.load(jfile)

lumifile = args.jsonpath.replace('golden_Cert','lumi').replace('.json','.txt')
if (not os.path.exists('lumi.csv')):
  print "First obtain luminosity for each run with following command (requires brilcalc installation):"
  print "brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i "+args.jsonpath+" -u /pb -o "+lumifile

lumidict = {}
with open(lumifile) as lfile:
  for line in lfile:
    if line[0]=='#': continue
    lumidict[line.split(':')[0]] = line.split(',')[-1]

ijson = 0
lumi = 0
newjson = {}
runs = sorted(jdata.keys())
for run in runs:
  newjson[run] = jdata[run]
  lumi = lumi + float(lumidict[run])
  if (lumi > pb_per_json) or run==runs[-1]:
    subjson = args.jsonpath.split('golden_')[0]+'subgolden_'+period_name+str(ijson)+'.json'
    json.dump(newjson, open(subjson,'w'), sort_keys=True)
    with open(subjson) as fsubj:
      data = fsubj.readlines()
      with open(subjson+'.tmp','w') as tmp_fsubj:
        for line in data:
          tmp_fsubj.write(line.replace(', "',',\n "'))
        tmp_fsubj.write('\n')
    os.rename(subjson+'.tmp', subjson)
    print "Wrote json %s with lumi %.2f" % (subjson, lumi)
    ijson = ijson + 1
    lumi=0
    newjson = {}
