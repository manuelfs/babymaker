###########################################################
### Configuration file to make basic babies from miniAOD
###########################################################
import math
from   os import environ
from   os.path import exists, join

def findFileInPath(theFile):
    for s in environ["CMSSW_SEARCH_PATH"].split(":"):
        attempt = join(s,theFile)
        if exists(attempt):
            return attempt
    print "Could not find file "+theFile
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
                 'txt/json/nonblind_Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.json',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Path to json starting with babymaker/...")
options.parseArguments()
outName = options.outputFile

## This refers to the official JEC methods. To apply on-the-fly, just set jecLabel to something different from miniAOD
doJEC = False  
if doJEC: jets_label = "patJetsReapplyJEC"
else: jets_label = "slimmedJets"
jecLabel = 'miniAOD'
if outName.find("Run2015") != -1:
    isData = True
    jecLabel = 'Summer15_25nsV5_DATA'
    globalTag = "74X_dataRun2_v2"
    processRECO = "RECO"
    jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
else:
    isData = False
    jecLabel = 'Summer15_25nsV2_MC'
    globalTag = "74X_mcRun2_asymptotic_v2"
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

process.baby_basic = cms.EDAnalyzer('bmaker_basic',
                                    outputFile = cms.string(options.outputFile),
                                    inputFiles = cms.vstring(options.inputFiles),
                                    jec = cms.string(jecLabel),
                                    met = cms.InputTag("slimmedMETs"),
                                    met_nohf = cms.InputTag("slimmedMETsNoHF"),
                                    jets = cms.InputTag(jets_label),
                                    nEventsSample = cms.uint32(options.nEventsSample)
)

###### Setting up number of events, and reporing frequency 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

###### Setting global tag 
## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = globalTag


###### HBHE
## HBHE noise filter needs to be recomputed in early 2015 data
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

if doJEC:
    ###### Setting sqlite file for the JECs that are in newer global tags 
    ## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string('sqlite_file:txt/jec/'+jecLabel+'.db'),
                               toGet   = cms.VPSet(
                                   cms.PSet(
                                       record = cms.string("JetCorrectionsRecord"),
                                       tag    = cms.string("JetCorrectorParametersCollection_"+jecLabel+"_AK4PFchs"),
                                       label  = cms.untracked.string("AK4PFchs")
                                   ),
	                cms.PSet(
                            record = cms.string("JetCorrectionsRecord"),
                            tag    = cms.string("JetCorrectorParametersCollection_"+jecLabel+"_AK4PF"),
                            label  = cms.untracked.string("AK4PF")
                        )
	        )
    )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
    ###### Jets 
    ## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJets"),
        levels = jecLevels,
        payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.patJetsReapplyJEC = patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
    )

    ###### Apply new JECs to MET
    ## From https://github.com/cms-met/cmssw/blob/METCorUnc74X/PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py
    process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False)
    )
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process,
                               isData=isData,
                               pfCandColl=cms.InputTag("packedPFCandidates")
                           )
    if isData:
        process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
        process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    # Redefine process for data to RECO
    process.slimmedMETs.t01Variation = cms.InputTag("slimmedMETs","",processRECO)

    ###### Path
    process.p = cms.Path(process.patJetCorrFactorsReapplyJEC*
                         process.patJetsReapplyJEC*
                         process.HBHENoiseFilterResultProducer*
                         process.baby_basic)
else:
    ###### Path
    process.p = cms.Path(process.HBHENoiseFilterResultProducer*
                         process.baby_basic)
