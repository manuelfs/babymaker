import FWCore.ParameterSet.Config as cms

# from FWCore.ParameterSet.VarParsing import VarParsing
# options = VarParsing ('optparser')
# options.parseArguments()

process = cms.Process("Baby")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/m/manuelf/work/data/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9-v2_MINIAODSIM.root'
    )
)

process.baby = cms.EDAnalyzer('bmaker')

# process.TFileService = cms.Service("TFileService",
#                                    fileName = cms.string('baby.root')
#                                    )


process.p = cms.Path(process.baby)
