#!/usr/bin/env python
import os, sys, time, re
import json
import urllib, urllib2, httplib, cookielib
import pprint
from   optparse import OptionParser
from   types import GeneratorType

# Data retrieval function taken/adapted from the DAS command line tool: https://cmsweb.cern.ch/das/cli
DAS_CLIENT = 'das-client/1.1::python/%s.%s' % sys.version_info[:2]
def get_data(query):
  """Contact DAS server and retrieve data for given DAS query"""
  params  = {'input':query, 'idx':0, 'limit':0}
  path    = '/das/cache'
  url = "https://cmsweb.cern.ch" + path
  client = '%s (%s)' % (DAS_CLIENT, os.environ.get('USER', ''))
  headers = {"Accept": "application/json", "User-Agent": client}
  encoded_data = urllib.urlencode(params, doseq=True)
  url += '?%s' % encoded_data
  req  = urllib2.Request(url=url, headers=headers)
  http_hdlr  = urllib2.HTTPHandler(debuglevel=0)
  proxy_handler  = urllib2.ProxyHandler({})
  cookie_jar     = cookielib.CookieJar()
  cookie_handler = urllib2.HTTPCookieProcessor(cookie_jar)
  opener = urllib2.build_opener(http_hdlr, proxy_handler, cookie_handler)
  fdesc = opener.open(req)
  data = fdesc.read()
  fdesc.close()

  pat = re.compile(r'^[a-z0-9]{32}')
  if  data and isinstance(data, str) and pat.match(data) and len(data) == 32: pid = data
  else: pid = None
  iwtime  = 2  # initial waiting time in seconds
  wtime   = 20 # final waiting time in seconds
  sleep   = iwtime
  time0   = time.time()
  while pid:
    params.update({'pid':data})
    encoded_data = urllib.urlencode(params, doseq=True)
    url  = "https://cmsweb.cern.ch" + path + '?%s' % encoded_data
    req  = urllib2.Request(url=url, headers=headers)
    try:
      fdesc = opener.open(req)
      data = fdesc.read()
      fdesc.close()
    except urllib2.HTTPError as err:
      return {"status":"fail", "reason":str(err)}
    if  data and isinstance(data, str) and pat.match(data) and len(data) == 32: pid = data
    else: pid = None
    time.sleep(sleep)
    if  sleep < wtime: sleep *= 2
    elif sleep == wtime: sleep = iwtime # start new cycle
    else: sleep = wtime
    if  (time.time()-time0) > 300:
      return {"status":"fail", "reason":("client timeout after %s sec" % int(time.time()-time0))}
  jsondict = json.loads(data)
  if  ('status' not in jsondict) or jsondict['status'] != 'ok':
    print('DAS record with bad status or without status field:\n%s' % jsondict)
    sys.exit()
    
  return jsondict

def getFilesInfo(dataset, wanted_keys = ['name','size','nevents'], verbose = False):
  jsondict = get_data('file dataset='+dataset)
  fdicts = []
  for entry in jsondict['data']:
    # pprint.pprint(entry)
    orig_fdict = entry['file']
    skim_fdict = {}
    for key in wanted_keys: 
      for i in range(0, len(orig_fdict)):
        if key in orig_fdict[i].keys():
          skim_fdict[key] = orig_fdict[i][key]
      if verbose: print skim_fdict[key],
    if verbose: print 
    fdicts.append(skim_fdict)

  return fdicts

def getDatasetInfo(dataset, wanted_keys = ['name','size','nevents','nfiles'], verbose = False):
  jsondict_ds = get_data('dataset='+dataset)
  # pprint.pprint(jsondict_ds['data'])
  # what file attributes do we want to keep track of
  orig_dsdict = jsondict_ds['data'][0]['dataset']
  skim_dsdict = {}
  for key in wanted_keys:
    for i in range(0, len(orig_dsdict)):
      if key in orig_dsdict[i].keys():
        skim_dsdict[key] = orig_dsdict[i][key] # don't break the loop, get the last (?) 
    if verbose: print skim_dsdict[key],
  if verbose: print

  return skim_dsdict

# Use a wildcard to retrieve info for multiple datasets
def getDatasetsInfo(dataset, wanted_keys = ['name','size','nevents','nfiles'], verbose = False):
  dsdicts = []
  if '*' in dataset:
    jsondict = get_data('dataset='+dataset) 
    # if using a wildcard, I need to get all the names first, 
    # because if I query for multiple datasets at a time it returns only minimal info for each
    for entry in jsondict['data']:
      dsname = entry['dataset'][0]['name']
      skim_dsdict = getDatasetInfo(dsname, wanted_keys = wanted_keys, verbose = verbose)
      dsdicts.append(skim_dsdict)
  else:
    skim_dsdict = getDatasetInfo(dataset, wanted_keys = wanted_keys, verbose = verbose)
    dsdicts.append(skim_dsdict)
  
  return dsdicts

# test
# getDatasetsInfo("/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_*/RunIISpring15DR74*/MINIAODSIM", verbose = True)
# getDatasetsInfo("/ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1/MINIAODSIM", verbose = True)

