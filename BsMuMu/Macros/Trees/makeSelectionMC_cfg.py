import FWCore.ParameterSet.Config as cms

process = cms.Process("SELECTION")

from YZheng.UpsilonAna.selection_cff import *
selection(process, GlobalTag="START53_V7G::All", MC=True, SelectionTrigger="hltDoubleMu2BsL3Filtered")

process.source.fileNames = cms.untracked.vstring(
    'file:/uscms_data/d3/demattia/BsMuMu/CMSSW_5_3_2_patch1/src/HeavyFlavorAnalysis/Onia2MuMu/test/onia2MuMuPAT.root'
)   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

