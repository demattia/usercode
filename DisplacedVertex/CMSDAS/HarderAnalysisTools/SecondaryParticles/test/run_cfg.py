import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('HarderAnalysisTools.SecondaryParticles.secondaryParticles_cfi')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_1/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V1-v2/0003/FC835340-F06B-DE11-A5CB-001D09F2960F.root'
    )
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('histograms.root')
)

process.p = cms.Path( process.secondaryParticles )

# print-out all python configuration parameter information
#print process.dumpPython()
