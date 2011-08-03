import FWCore.ParameterSet.Config as cms

process = cms.Process("EfficiencyMerger")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.demo = cms.EDAnalyzer('EfficiencyMerger',
                              InputFileName = cms.vstring("EfficiencyCleaned_1.root", "EfficiencyCleaned_2.root")
                              )

process.p = cms.Path(process.demo)
