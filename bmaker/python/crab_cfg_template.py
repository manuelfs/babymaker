## Template file for CRAB submission. The script generate_crab_config.py 
## replaces the following two lines with the appropriate values
## Do not edit manually!
dataset = 'DATASETNAME'
nevents = NEVENTS

jobname = dataset[1:].replace('/','__')
jobname = jobname.replace(':','___')

datasetID = dataset.replace('/','',1).replace('/', '_', 1)
datasetID = datasetID[0:datasetID.find('/')]

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = jobname
config.General.workArea = 'out_crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'babymaker/bmaker/python/bmaker_basic_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['baby_' + datasetID + '.root']
config.JobType.pyCfgParams = ['nEventsSample=' + str(nevents), 'outputFile=baby_' + datasetID + '.root']

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False # used to be True for cfA production
# config.Data.publishDBS = 'phys03'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'

config.section_("Site")
config.Site.storageSite = 'T3_US_UCSB' #'T3_US_UCSB'
#config.Site.whitelist = ['T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT', 'T3_US_UCSB']
config.Site.blacklist = ['T1_RU_JINR']
# you may want to uncomment this line and force jobs to run in the US
# only a few datasets (mostly very new ones) will not be accessible
