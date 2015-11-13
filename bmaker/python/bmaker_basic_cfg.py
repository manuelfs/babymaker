###########################################################
### Configuration file to make basic babies from miniAOD
###########################################################
import math
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
                 'babymaker/txt/json/golden_Cert_246908-258714_13TeV_PromptReco_Collisions15_25ns_JSON.json',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Path to json starting with babymaker/...")
options.register('condorSubTime',
                 '000000_000000',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Timestamp from condor submission")
options.register('isCRABJob',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Switch to change directory structure")
options.parseArguments()
outName = options.outputFile
if outName == "output.root": # output filename not set
    rootfile = basename(options.inputFiles[0])
    outName = "baby_"+rootfile

doSystematics = True

## This refers to the official JEC methods. To apply on-the-fly, just set jecLabel to something different from miniAOD
doJEC = False  
if doJEC: jets_label = "patJetsReapplyJEC"
else: jets_label = "slimmedJets"

# jecLabel must contain the JEC version; this is needed for uncertainty calculation
# if there is no need to apply JECs then the jecLabel must also contain 'miniAOD' (as well as the version)
jecLabel = 'miniAOD_Summer15_25nsV6_MC' # for 7.4.14 mc, don't apply JEC, but still give the JEC tag because of systematics
if "Run2015D" in outName: jecLabel = 'Summer15_25nsV6_DATA' # for 7.4.12 data
elif "RunIISpring15DR74" in outName: jecLabel = 'Summer15_25nsV6_MC' # for 7.4.6.patch4 mc

if "Run2015" in outName:
    isData = True
    # These only used for the official application of JECs
    globalTag = "74X_dataRun2_v2"
    processRECO = "RECO"
    jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
else:
    isData = False
    # These only used for the official application of JECs
    globalTag = "74X_mcRun2_asymptotic_v2"
    processRECO = "PAT"
    jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']

# The 7.4.14 re-miniAOD already has V5 JECs with the new prescription
cmsswRel = environ["CMSSW_BASE"]
# if cmsswRel.find("CMSSW_7_4_14") != -1: jecLabel = 'miniAOD' 

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
                                    isCRABJob = cms.bool(options.isCRABJob)
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
##___________________________HCAL_Noise_Filter________________________________||
if "FSPremix" in outName: fastsim = True
else: fastsim = False
if not fastsim:
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

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
                         process.HBHENoiseFilterResultProducer* #produces HBHE baseline bools
                         process.baby_basic)
else:
    ###### Path
    if not fastsim:
        process.p = cms.Path(process.HBHENoiseFilterResultProducer* #produces HBHE baseline bools
                             process.baby_basic)
    else:
        process.p = cms.Path(process.baby_basic)
