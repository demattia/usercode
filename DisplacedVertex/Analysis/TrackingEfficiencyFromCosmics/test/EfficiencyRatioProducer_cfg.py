import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingEfficiencyFromCosmics")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.demo = cms.EDAnalyzer('EfficiencyRatioProducer',
                              InputFileNameNumerator = cms.string("Efficiency_Data.root"),
                              InputFileNameDenominator = cms.string("Efficiency_Sim.root")
                              )

process.p = cms.Path(process.demo)
