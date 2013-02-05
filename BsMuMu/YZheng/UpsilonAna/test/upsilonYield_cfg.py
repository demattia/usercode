import FWCore.ParameterSet.Config as cms

process = cms.Process("YIELD")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    'file:/uscms_data/d2/zgecse/EJTerm/skim_1.87pb-1.root'
  )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('upsilonYield.root')
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.goodMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull && pt > 3.5 && abs(eta) < 2.1"),
)

process.goodUpsilons  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goodMuons@+ goodMuons@-"),
    cut = cms.string("7 < mass < 12"),
)

process.upsilonYield = cms.EDAnalyzer("YieldTree",
  UpsilonCandidates = cms.untracked.InputTag("goodUpsilons")
)

process.path = cms.Path(
  process.goodMuons * process.goodUpsilons * process.upsilonYield
)

