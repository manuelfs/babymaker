import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()

process = cms.Process("Baby")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.baby = cms.EDAnalyzer('bmaker_basic',
    outputFile = cms.string(options.outputFile)
)

# process.TFileService = cms.Service("TFileService",
#                                    fileName = cms.string('baby.root')
#                                    )


process.p = cms.Path(process.baby)
