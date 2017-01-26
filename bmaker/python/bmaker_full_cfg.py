###########################################################
### Configuration file to make full babies from miniAOD
###########################################################
import math, sys
from   os import environ
from   os.path import exists, join, basename

def findFileInPath(theFile):
    for s in environ["CMSSW_SEARCH_PATH"].split(":"):
        attempt = join(s,theFile)
        if exists(attempt):
            return attempt
    print "BABYMAKER: Could not find file "+theFile+". Not pre-applying JSON"
    return None

###### Input parameters parsing 
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('nEventsSample',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Total number of events in dataset for event weight calculation.")
options.register('nEvents',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Number of events to run over.")
options.register('json',
                 'babymaker/data/json/golden_Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt', 
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Path to json starting with babymaker/...")
options.register('condorSubTime',
                 '000000_000000',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Timestamp from condor submission")
options.parseArguments()
outName = options.outputFile
if outName == "output.root": # output filename not set
    rootfile = basename(options.inputFiles[0])
    outName = "fullbaby_"+rootfile.replace("file:","")

doSystematics = True
cmsswRel = environ["CMSSW_BASE"]
if "RunIISpring15DR74" in outName or "RunIISpring15FSPremix" in outName: 
    if cmsswRel.find("CMSSW_7_4_6") == -1: sys.exit("ERROR: Trying to run miniAOD V1 in a new release. Exiting")

## JECs must be undone and reapplied when rerunning b-tagging
## => if doJEC = False, DeepCSV discriminator will not be included
doJEC = True  
if doJEC: jets_label = "updatedPatJetsTransientCorrectedUpdatedJEC"
else: jets_label = "slimmedJets"

# to apply JECs with txt files in babymaker, 
# prefix jecLabel with "onthefly_", e.g. onthefly_Spring16_25nsV6_MC
# systematics will also be calculated using this tag, even if JECs are not re-applied
# N.B. JECs change in the middle of RunF, thereby the Run2016F1 vs Run2016F2 distinction
jecLabel = 'onthefly_Spring16_25nsV6_MC'
if ("Run2016B" in outName) or ("Run2016C" in outName) or ("Run2016D" in outName): 
  jecLabel = 'Summer16_23Sep2016BCDV2_DATA'
elif ("Run2016E" in outName) or ("Run2016F1" in outName):
  jecLabel = 'Summer16_23Sep2016EFV2_DATA'
elif ("Run2016F2" in outName) or ("Run2016G" in outName):
  jecLabel = 'Summer16_23Sep2016GV2_DATA'
elif ("Run2016H" in outName): 
  jecLabel = 'Summer16_23Sep2016HV2_DATA'
elif "RunIISpring16MiniAOD" in outName:
    if "Fast" in outName or "FSPremix" in outName: jecLabel = 'onthefly_Spring16_FastSimV1_MC'
    else: jecLabel = 'Spring16_23Sep2016V2_MC' # for ICHEP MC against re-reco data
elif "RunIISummer16MiniAOD" in outName:
    jecLabel = 'Summer16_23Sep2016V2_MC'

if "FSPremix" in outName or "Fast" in outName: fastsim = True
else: fastsim = False

if "Run2016" in outName:
    isData = True
    # These only used for the official application of JECs
    globalTag = "80X_dataRun2_2016SeptRepro_v6"
    processRECO = "RECO"
    jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
else:
    isData = False
    # These only used for the official application of JECs
    globalTag = "80X_mcRun2_asymptotic_2016_miniAODv2"
    if "RunIISummer16MiniAOD" in outName: globalTag = "80X_mcRun2_asymptotic_2016_TrancheIV_v7"
    processRECO = "PAT"
    jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']

###### Defining Baby process, input and output files 
process = cms.Process("Baby")
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
)
if isData: # Processing only lumis in JSON
    import FWCore.PythonUtilities.LumiList as LumiList
    jsonfile = findFileInPath(options.json)
    process.source.lumisToProcess = LumiList.LumiList(filename = jsonfile).getVLuminosityBlockRange()
    doSystematics = False

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

###### Setting global tag 
## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = globalTag

process.baby_full = cms.EDAnalyzer('bmaker_full',
                                    condor_subtime = cms.string(options.condorSubTime),
                                    outputFile = cms.string(outName),
                                    inputFiles = cms.vstring(options.inputFiles),
                                    json = cms.string(options.json),
                                    jec = cms.string(jecLabel),
                                    met = cms.InputTag("slimmedMETs"),
                                    met_nohf = cms.InputTag("slimmedMETsNoHF"),
                                    jets = cms.InputTag(jets_label),
                                    nEventsSample = cms.uint32(options.nEventsSample),
                                    doMetRebalancing = cms.bool(True),
                                    doSystematics = cms.bool(doSystematics),
                                    addBTagWeights = cms.bool(True),
                                    isFastSim = cms.bool(fastsim),
                                    debugMode = cms.bool(False)
)

###### Setting up number of events, and reporing frequency 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100000


if not fastsim:
    process.load('Configuration.StandardSequences.Services_cff')

    ##_____________Bad charged candidate filter________________||
    process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
    process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
    process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
    process.BadChargedCandidateFilter.taggingMode = cms.bool(True)
    process.BadChargedCandidateFilter.debug = cms.bool(False)

    ##_____________Bad muon filter_____________________________||
    process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
    process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
    process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
    process.BadPFMuonFilter.taggingMode = cms.bool(True)
    process.BadPFMuonFilter.debug = cms.bool(False)

if doJEC:
    process.options = cms.untracked.PSet(
      allowUnscheduled = cms.untracked.bool(True),
      wantSummary = cms.untracked.bool(False)
    )
    ###### Setting sqlite file for the JECs that are in newer global tags 
    ## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string('sqlite_fip:babymaker/data/jec/'+jecLabel+'.db'),
                               toGet   = cms.VPSet(
                                   cms.PSet(
                                       record = cms.string("JetCorrectionsRecord"),
                                       tag    = cms.string("JetCorrectorParametersCollection_"+jecLabel+"_AK4PFchs"),
                                       label  = cms.untracked.string("AK4PFchs")
                                   )
                )
    )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
    ###### Applying JECs and including deepCSV info
    ## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
      process,
      jetSource = cms.InputTag('slimmedJets'),
      labelName = 'UpdatedJEC',
      jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None'),
      btagDiscriminators = ["deepFlavourJetTags:probudsg", 
                            "deepFlavourJetTags:probb", 
                            "deepFlavourJetTags:probc", 
                            "deepFlavourJetTags:probbb", 
                            "deepFlavourJetTags:probcc"]
    )
    process.options.allowUnscheduled = cms.untracked.bool(True)

    process.content = cms.EDAnalyzer("EventContentAnalyzer")

    ###### Apply new JECs to MET
    ## From https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    ## If you only want to re-correct and get the proper uncertainties, no reclustering
    runMetCorAndUncFromMiniAOD(process,
                               isData = isData,
    )

###### Path
process.p = cms.Path(process.baby_full)
