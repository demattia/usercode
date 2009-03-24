import FWCore.ParameterSet.Config as cms    

process = cms.Process("TEST")
process.add_(cms.Service( "MessageLogger"))
process.source = cms.Source( "EmptySource" )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("RecoLocalTracker.SiStripClusterizer.test.ClusterizerUnitTestsTest_cff")
testDefinition = cms.VPSet() + [ process.TestGroup ]

process.runUnitTests = cms.EDAnalyzer("ClusterizerUnitTesterTest",       TestGroups = testDefinition  )

process.path = cms.Path( process.runUnitTests )
