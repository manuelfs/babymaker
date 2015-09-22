###########################################################
### Configuration file to make basic babies from miniAOD
###########################################################
import math

print "\nThis cfg file applies the latest JECs to MET. It REQUIRES\n  git cms-merge-topic -u cms-met:METCorUnc74X\n"

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
options.parseArguments()

###### Defining Baby process, input and output files 
process = cms.Process("Baby")
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
)
process.baby_basic = cms.EDAnalyzer('bmaker_basic',
                                    outputFile = cms.string(options.outputFile),
                                    met = cms.InputTag("slimmedMETsTypeICorr"),
                                    nEventsSample = cms.uint32(options.nEventsSample)
)

###### Setting up number of events, and reporing frequency 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nEvents) )
ntot = options.nEvents if options.nEvents!=-1 else options.nEventsSample
process.MessageLogger.cerr.FwkReport.reportEvery = min(int(math.pow(10,math.floor(math.log(ntot,10)))),10000)

###### Setting global tag 
## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

###### Setting sqlite file for the JECs that are in newer global tags 
## From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
JECname = 'Summer15_25nsV2_MC'
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                           connect = cms.string('sqlite_file:jec/db/'+JECname+'.db'),
                           toGet   = cms.VPSet(
                               cms.PSet(
                                   record = cms.string("JetCorrectionsRecord"),
                                   tag    = cms.string("JetCorrectorParametersCollection_"+JECname+"_AK4PFchs"),
                                   label  = cms.untracked.string("AK4PFchs")
                               ),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag    = cms.string("JetCorrectorParametersCollection_"+JECname+"_AK4PF"),
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
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute'],
    payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

###### HBHE
## HBHE noise filter needs to be recomputed in early 2015 data
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

###### Apply new JECs to MET 
## From https://github.com/cms-met/cmssw/blob/METCorUnc74X/PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    #wantSummary = cms.untracked.bool(True)
)
isData = False
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=isData,
                           pfCandColl=cms.InputTag("packedPFCandidates"),
                           postfix="TypeICorr"
)
if isData:
    process.patPFMetT1T2CorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2CorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDownTypeICorr.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUpTypeICorr.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")

###### Path
process.p = cms.Path(process.patJetCorrFactorsReapplyJEC*
                     process.patJetsReapplyJEC*
                     process.HBHENoiseFilterResultProducer*
                     process.baby_basic)
