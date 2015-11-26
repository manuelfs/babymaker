#!/usr/bin/env python
import das_client as das
import json
import os
import sys

def getNumberOfEvents(dataset):
	query = "file dataset=" + dataset + " | sum(file.nevents)"

	data = das.get_data(query)
	if isinstance(data, basestring):
		dasjson = json.loads(data)
	else:
		dasjson = data
	status  = dasjson.get('status')
	if  status == 'ok':
		data = dasjson.get('data')
		sumevents=0
		for idata in data:
			sumevents+=idata.get('result').get('value')
		return sumevents

#######################################################

print "####################################################################################"
print "# Warning! This script will not set the correct weight for samples with extensions #" 
print "####################################################################################"

if len(sys.argv) != 2:
	print "Must specify dataset name"
else:
	dataset = sys.argv[1]
	nevents = getNumberOfEvents(dataset)

	cmssw_base = os.getenv("CMSSW_BASE")
	datasetID = dataset.replace('/','',1).replace('/', '_', 1)
	datasetID = datasetID[0:datasetID.find('/')]
	inputfile = cmssw_base + "/src/babymaker/bmaker/python/crab_cfg_template.py"
	outputfile = "crab_cfg_" + datasetID + ".py"

	s = open(inputfile).read()
	s = s.replace('DATASETNAME', dataset)
	s = s.replace('NEVENTS', str(nevents))
	f = open(outputfile, 'w')
	f.write(s)
	f.close()
	print "Dataset " + dataset + " has " + str(nevents) + " events."
	print "Wrote crab configuration file " + outputfile
