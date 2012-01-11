import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingEfficiencyFromCosmics")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.demo = cms.EDAnalyzer('EfficiencyRatioProducer',
                              InputFileNameNumerator = cms.string("efficiency_data/EfficiencyCleaned.root"),
                              InputFileNameDenominator = cms.string("efficiency_sim/EfficiencyCleaned.root"),
							  Rebin = cms.uint32(4)
                              )

process.p = cms.Path(process.demo)
