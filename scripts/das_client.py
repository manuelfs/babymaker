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
    print('DAS record with bad status or without status field:\n')
    pprint.pprint(jsondict)
    sys.exit()
    
  return jsondict

def findKeyValue(data, key):
  value = 0
  if isinstance(data, list):
    found_key = False
    for i in range(0, len(data)):
      if key in data[i].keys():
        if found_key and value!=data[i][key]: 
          pprint.pprint(data)
          print "ERROR: Found multiple instances of key \'%s\'." % key
          sys.exit(0)
        else:
          found_key = True
          value = data[i][key]
    if not found_key:   
      # pprint.pprint(data)
      print "WARNING: Returning NULL. Could not find key \'%s\' in list." % key
      value = "NULL"
      # sys.exit(0)       
  elif isinstance(data, dict):
    if key in data.keys():
      value = data[key]
    else:   
      # pprint.pprint(data)
      print "WARNING: Returning NULL. Could not find key \'%s\' in dict." % key
      value = "NULL"
      # sys.exit(0)    
  else:
    # pprint.pprint(data)
    print "WARNING: Returning NULL. Dictionary is of neither type list or dict:"
    value = "NULL"    
  return value

def getFileRunInfo(file, getlumis = False, verbose = False):
  jsondict = get_data('lumi file='+file)
  data = jsondict['data']

  # when looking at data (not mc), we expect a list of dictionaries 
  # because a number of runs got combined into one if it was not prompt reco
  rundict = {}
  for idata in data: 
    orig_lumidict = findKeyValue(idata, 'lumi')
    run_number = findKeyValue(orig_lumidict,'run_number')
    if (run_number=="NULL"): run_number = -1 #expecting an int
    if (verbose): print "Found run %s in file %s." % (run_number, file)
    rundict[run_number] = []
    lumis = findKeyValue(orig_lumidict, 'number')
    if (lumis=="NULL"): 
      rundict[run_number].append(-1)
    else:
      for ll in lumis:
        rundict[run_number].extend([i for i in range(ll[0],ll[1]+1)])

  if verbose: 
    print "Contents of %s" % file
    pprint.pprint(rundict)
  if (getlumis):
    return rundict
  else:
    return rundict.keys()

def getFilesInfo(dataset, wanted_keys = ['name','size','nevents'], verbose = False):
  jsondict = get_data('file dataset='+dataset)
  # pprint.pprint(jsondict)
  fdicts = []
  for entry in jsondict['data']:
    orig_fdict = entry['file']
    # if (orig_fdict['name'] = '/store/data/Run2015D/HTMHT/MINIAOD/PromptReco-v4/000/258/706/00000/90EEB778-5271-E511-B309-02163E014364.root'):
    #   print "asdgdfgadfgafdg"
    skim_fdict = {}
    for key in wanted_keys: 
      skim_fdict[key] = findKeyValue(orig_fdict, key)
      if (key=='size' or key=='nevents') and skim_fdict[key]=="NULL": skim_fdict[key] = -1
      if verbose: print skim_fdict[key],
    if verbose: print 
    fdicts.append(skim_fdict)

  return fdicts

def getDatasetInfo(dataset, wanted_keys = ['name','size','nevents','nfiles'], verbose = False):
  jsondict_ds = get_data('dataset='+dataset)
  # pprint.pprint(jsondict_ds)
  # what file attributes do we want to keep track of
  orig_dsdict = findKeyValue(jsondict_ds['data'],'dataset')
  skim_dsdict = {}
  for key in wanted_keys:
    skim_dsdict[key] = findKeyValue(orig_dsdict, key)
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
      ds_entry = findKeyValue(entry,'dataset')
      dsname = findKeyValue(ds_entry,'name')
      skim_dsdict = getDatasetInfo(dsname, wanted_keys = wanted_keys, verbose = verbose)
      dsdicts.append(skim_dsdict)
  else:
    skim_dsdict = getDatasetInfo(dataset, wanted_keys = wanted_keys, verbose = verbose)
    dsdicts.append(skim_dsdict)
  
  return dsdicts

# test
# answer = getDatasetsInfo("/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_*/RunIISpring15DR74*/MINIAODSIM", verbose = True)
# answer = getDatasetsInfo("/ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1/MINIAODSIM", verbose = True)
# answer = getFilesInfo("/HTMHT/Run2015D-PromptReco-v3/MINIAOD", verbose = True)
# answer = getFilesInfo("/HTMHT/Run2015D-PromptReco-v3/MINIAOD", verbose = True)
# answer = getFileRunInfo("/store/data/Run2015D/HTMHT/MINIAOD/05Oct2015-v1/10000/0A668427-B06F-E511-813E-0025904B12A8.root", lumis = True, verbose = True)

# answer = getDatasetsInfo("/HTMHT/Run2015D*/MINIAOD", verbose = True)
# answer = getDatasetsInfo("/MET/Run2015D*/MINIAOD", verbose = True)
# answer = getDatasetsInfo("/SingleElectron/Run2015D*/MINIAOD", verbose = True)
# answer = getDatasetsInfo("/SingleMuon/Run2015D*/MINIAOD", verbose = True)
# answer = getDatasetsInfo("/JetHT/Run2015D*/MINIAOD", verbose = True)
# answer = getDatasetsInfo("/DoubleEG/Run2015D*/MINIAOD", verbose = True)
# answer = getDatasetsInfo("/DoubleMuon/Run2015D*/MINIAOD", verbose = True)
