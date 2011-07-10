import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingEfficiencyFromCosmics")

# # process.load("FWCore.MessageService.MessageLogger_cfi")
# # initialize MessageLogger and output report
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger.cerr.threshold = 'INFO'
# process.MessageLogger.categories.append('Demo')
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)
# )
# process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# process.TFileService=cms.Service('TFileService',
#                                  fileName=cms.string('EfficiencyAnalyzer.root')
#                                  )

process.demo = cms.EDAnalyzer('EfficiencyAnalyzer',
                              InputFileName = cms.string("EfficiencyCleaned.root"),
                              # InputFileName = cms.string("Efficiency.root"),
                              # InputFileName = cms.string("GenToTrackEfficiency.root"),
                              # InputFileName = cms.string("GenToStandAloneEfficiency.root"),
                              Rebin = cms.uint32(4)
                              )

process.p = cms.Path(process.demo)
