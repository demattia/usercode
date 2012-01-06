import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load('HarderAnalysis.DisplacedDileptons.displacedDileptonAnalysis_cfi')

# data sample and sample-specific settings
sampleSignalPID=0
sampleRequireCollision=1
sampleRunE=1
sampleRunMu=1
from DisplacedLeptons.Samples.samples.Debug_sample_cff import *

# electrons
process.eAnalysis = process.analyzeDisplacedDileptons.clone()
process.eAnalysis.signalPDGId = sampleSignalPID
process.eAnalysis.leptonPDGId = 11
process.eAnalysis.decayLengthCut = -1
process.eAnalysis.decayLengthSignificanceCut = 8.0
process.eAnalysis.deltaPhiCut = 0.8
process.eAnalysis.leptonPtCut = 38
process.eAnalysis.minD0Significance = 3
process.eAnalysis.maxNumStandAloneMuons = 0
process.eAnalysis.hltPaths = cms.vstring("HLT_DoublePhoton10_L1R", # Run2010A1 (136035-141881)
                                         "HLT_DoublePhoton15_L1R", # Run2010A2 (141956-144114)
                                         "HLT_DoublePhoton17_L1R", # Run2010B1 (146428-148058)
                                         "HLT_Photon22_SC22HE_L1R_v1", # Run2010B2  (148822-150000)
                                         "HLT_DoublePhoton33_v1",
                                         "HLT_DoublePhoton33_v2",
                                         "HLT_DoublePhoton33_v3",
                                         "HLT_DoublePhoton33_v4",
                                         "HLT_DoublePhoton33_v5",
                                         "HLT_DoublePhoton33_HEVT_v1",
                                         "HLT_DoublePhoton33_HEVT_v2",
                                         "HLT_DoublePhoton33_HEVT_v3",
                                         "HLT_DoublePhoton33_HEVT_v4",
                                         "HLT_DoublePhoton38_HEVT_v1",
                                         "HLT_DoublePhoton38_HEVT_v2",
                                         "HLT_DoublePhoton38_HEVT_v3",
                                         "HLT_DoublePhoton43_HEVT_v1")

# muons
process.muAnalysis = process.analyzeDisplacedDileptons.clone()
process.muAnalysis.signalPDGId = sampleSignalPID
process.muAnalysis.leptonPDGId = 13
process.muAnalysis.decayLengthCut = -1
process.muAnalysis.vetoBackToBack = -0.95
process.muAnalysis.minDeltaRBetweenLeptons = 0.2
process.muAnalysis.leptonPtCut = 25
#process.muAnalysis.maxNumStandAloneMuons = 0  # seems to be a problem
process.muAnalysis.hltPaths = cms.vstring("HLT_L2Mu11",  # Run2010A1 (136035-141881)
                                          "HLT_L2Mu15",  # Run2010A2 (141956-144114)
                                          "HLT_L2Mu25",                   # Run2010B1 (146428-147120)
                                          "HLT_L2DoubleMu20_NoVertex_v1", # Run2010B2 (147121-150000)
                                          "HLT_L2DoubleMu23_NoVertex_v1",
                                          "HLT_L2DoubleMu23_NoVertex_v2",
                                          "HLT_L2DoubleMu23_NoVertex_v3",
                                          "HLT_L2DoubleMu23_NoVertex_v4",
                                          "HLT_L2DoubleMu23_NoVertex_v5",
                                          "HLT_L2DoubleMu23_NoVertex_v6",
                                          "HLT_L2DoubleMu23_NoVertex_v7",
                                          "HLT_L2DoubleMu30_NoVertex_v1",
                                          "HLT_L2DoubleMu30_NoVertex_v2",
                                          "HLT_L2DoubleMu30_NoVertex_v3",
                                          "HLT_L2DoubleMu30_NoVertex_v4")
                                          
# track-based analysis channels
process.eTrackAnalysis = process.eAnalysis.clone()
process.eTrackAnalysis.leptonPDGId = -11
process.eTrackAnalysis.minNumCaloMatches = 2
process.muTrackAnalysis = process.muAnalysis.clone()
process.muTrackAnalysis.leptonPDGId = -13
process.muTrackAnalysis.maxNumStandAloneMuons = 0

# create pseudo-leptons from tracks and standalone muons for track-based channels
process.load("HarderAnalysis.DisplacedDileptons.pseudoLeptonProducer_cfi")

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string(sampleDuplicateCheckMode),
    fileNames = cms.untracked.vstring(samplePatFiles)
)
process.GlobalTag.globaltag = sampleGlobalTag
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('histograms.root')
)

process.load("HarderAnalysis.DisplacedDileptons.goodCollFilter_cff")
process.load("DisplacedLeptons.Samples.isoTrackPrefilter_cfi")
process.prefilterPassTwo = process.isoTrackPrefilter.clone()

# print event tree
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTruth = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(20),
  printVertex = cms.untracked.bool(True),
  src = cms.InputTag("genParticles")
)



process.p = cms.Path()
if sampleRequireCollision: process.p+=process.prefilterPassTwo*process.goodcollFilter
process.p+=process.pseudoLeptonProducer
#if sampleRunE: process.p+=process.eAnalysis
#if sampleRunMu: process.p+=process.muAnalysis
if sampleRunE: process.p+=process.eTrackAnalysis
if sampleRunMu: process.p+=process.muTrackAnalysis
if sampleType=="MC":
    process.p+=process.printTruth
