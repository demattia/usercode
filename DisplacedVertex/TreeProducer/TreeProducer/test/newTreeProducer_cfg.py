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

process.load('TreeProducer.TreeProducer.leptonAnalysis_cfi')

# data sample and sample-specific settings
sampleSignalPID=6000113
sampleRequireCollision=1
sampleRunE=1
sampleRunMu=1

from SampleFiles.Samples.Data_Mu_Run2012A1_cff import *

# electrons
process.eAnalysis = process.analyzeLeptons.clone()
process.eAnalysis.signalPDGId = sampleSignalPID
process.eAnalysis.leptonPDGId = 11
process.eAnalysis.decayLengthCut = -1
process.eAnalysis.decayLengthSignificanceCut = 8.0
process.eAnalysis.deltaPhiCut = 0.8
process.eAnalysis.leptonPtCut = 38
process.eAnalysis.minD0Significance = 3
process.eAnalysis.maxNumStandAloneMuons = 0
process.eAnalysis.hltPaths = cms.vstring("*")

# muons
process.muAnalysis = process.analyzeLeptons.clone()
process.muAnalysis.signalPDGId = sampleSignalPID
process.muAnalysis.leptonPDGId = 13
process.muAnalysis.decayLengthCut = -1
process.muAnalysis.vetoBackToBack = -0.95
process.muAnalysis.minDeltaRBetweenLeptons = 0.2
process.muAnalysis.leptonPtCut = 25
process.muAnalysis.FillHistograms = True
# process.muAnalysis.maxNumStandAloneMuons = 0  # seems to be a problem
# process.muAnalysis.hltPaths = cms.vstring("*")
process.muAnalysis.hltPaths = cms.vstring("HLT_Mu17_Mu8_v1",
                                          "HLT_Mu17_Mu8_v2",
                                          "HLT_Mu17_Mu8_v3",
                                          "HLT_Mu17_Mu8_v4",
                                          "HLT_Mu17_Mu8_v5",
                                          "HLT_Mu17_Mu8_v6",
                                          "HLT_Mu17_Mu8_v7",
                                          "HLT_Mu17_Mu8_v8",
                                          "HLT_Mu17_Mu8_v9",
                                          "HLT_Mu17_Mu8_v10",
                                          "HLT_Mu17_Mu8_v11",
                                          "HLT_Mu17_Mu8_v12",
                                          "HLT_Mu17_Mu8_v13",
                                          "HLT_Mu17_Mu8_v14",
                                          "HLT_Mu17_Mu8_v15",
                                          "HLT_Mu17_Mu8_v16",
                                          "HLT_Mu17_Mu8_v17",
                                          "HLT_Mu17_Mu8_v18",
                                          "HLT_Mu17_Mu8_v19",
                                          "HLT_Mu17_Mu8_v20",
                                          "HLT_Mu17_Mu8_v21",
                                          "HLT_Mu17_Mu8_v22",
                                          "HLT_Mu17_Mu8_v23",
                                          "HLT_Mu17_Mu8_v24",
                                          "HLT_Mu17_Mu8_v25",
                                          "HLT_Mu17_Mu8_v26",
                                          "HLT_Mu17_Mu8_v27",
                                          "HLT_Mu17_Mu8_v28",
                                          "HLT_Mu17_Mu8_v29",
                                          "HLT_Mu17_Mu8_v30",
                                          "HLT_Mu17_Mu8_v31")

# track-based analysis channels
process.eTrackAnalysis = process.eAnalysis.clone()
process.eTrackAnalysis.leptonPDGId = -11
process.eTrackAnalysis.minNumCaloMatches = 2
process.muTrackAnalysis = process.muAnalysis.clone()
process.muTrackAnalysis.leptonPDGId = -13
process.muTrackAnalysis.maxNumStandAloneMuons = 0

# create pseudo-leptons from tracks and standalone muons for track-based channels
process.load("TreeProducer.TreeProducer.pseudoLeptonProducer_cfi")

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string(sampleDuplicateCheckMode),
    fileNames = cms.untracked.vstring("file:PATtuple_161_1_tTu._WZroot")
)
process.GlobalTag.globaltag = sampleGlobalTag
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('histograms.root')
)

process.load("TreeProducer.TreeProducer.goodCollFilter_cff")

process.p = cms.Path()
if sampleRequireCollision: process.p+=process.goodcollFilter
process.p+=process.pseudoLeptonProducer
#if sampleRunE: process.p+=process.eAnalysis
#if sampleRunMu: process.p+=process.muAnalysis
if sampleRunE: process.p+=process.eTrackAnalysis
if sampleRunMu: process.p+=process.muTrackAnalysis
