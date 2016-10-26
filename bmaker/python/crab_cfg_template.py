## Template file for CRAB submission. The script generate_crab_config.py 
## replaces the following two lines with the appropriate values
## Do not edit manually!
dataset = 'DATASETNAME'
nevents = NEVENTS

# CRAB3 task names can no longer be greater than 100 characters; need to shorten task name
# Do NOT replace: "PUSpring16Fast" since it is used by the bmaker_full code to decide how to read GenInfo !!!
taskname = dataset[1:].replace('/','__')
taskname = taskname.replace('RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2','MiniAODv2')
taskname = taskname.replace('TuneCUETP8M1_13TeV-madgraphMLM-pythia8','13TeV-MG-PY8')
taskname = taskname.replace('RunIISpring15MiniAODv2-Asympt25ns_74X_mcRun2_asymptotic_v2','MiniAODv2')
taskname = taskname.replace('RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016','80XMiniAODv1')
taskname = taskname.replace('RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2','80XMiniAODv2')
taskname = taskname.replace(':','___')

if(len(taskname)>100): taskname = taskname[0:99]

datasetID = dataset.replace('/','',1).replace('/', '_', 1)
datasetID = datasetID[0:datasetID.find('/')]

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = taskname
config.General.workArea = 'out_crab'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'babymaker/bmaker/python/bmaker_full_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['fullbaby_' + datasetID + '.root']
config.JobType.pyCfgParams = ['nEventsSample=' + str(nevents), 'outputFile=fullbaby_' + datasetID + '.root']

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
if "Run201" in taskname:
    config.Data.splitting = 'LumiBased'
    config.Data.unitsPerJob = 100
    config.Data.lumiMask = 'babymaker/data/json/Cert_271036-283059_minus_271036-277933.json' 
                          

else:
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 5

config.Data.publication = False # used to be True for cfA production
# config.Data.publishDBS = 'phys03'

config.section_("Site")
config.Site.storageSite = 'T2_US_UCSD'
#config.Site.storageSite = 'T3_US_UCSB'
#config.Site.whitelist = ['T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT', 'T3_US_UCSB']
config.Site.blacklist = ['T1_RU_JINR', 'T1_FR_CCIN2P3']
# you may want to uncomment this line and force jobs to run in the US
# only a few datasets (mostly very new ones) will not be accessible
