import FWCore.ParameterSet.Config as cms

process = cms.Process("CounterReader")

# initialize MessageLogger and output report
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:/Users/demattia/CMSSH/soft/Releases/CMSSW_5_3_8/src/onia2MuMuPAT_MC.root'
                                )
                            )

process.counterReader = cms.EDAnalyzer('CounterReader')

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('counters.root')
                                   )


process.p = cms.Path(process.counterReader)
