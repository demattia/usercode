import FWCore.ParameterSet.Config as cms

process = cms.Process("SELECTION")

from YZheng.UpsilonAna.selection_cff import *
selection(process, GlobalTag="GR_P_V40::All", MC=False, SelectionTrigger="hltDoubleMu2BsL3Filtered")

process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(
   'file:/uscms/home/zhenhu/onia2MuMuPAT_MC_1_1_4ei.root'
 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
