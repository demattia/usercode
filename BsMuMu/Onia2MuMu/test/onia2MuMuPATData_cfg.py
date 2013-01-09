# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms

process = cms.Process("Onia2MuMuPAT")

from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import *

onia2MuMuPAT(process, GlobalTag="GR_R_42_V14::All", MC=False, HLT="HLT", Filter=True)

process.source.fileNames = cms.untracked.vstring(
      # '/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0007/FE96320B-CFEA-DF11-AF94-E0CB4E553666.root',
      # '/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0007/FE72AF8A-CEEA-DF11-B0B9-E0CB4E29C51D.root',
      # '/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0007/FC9C8211-CFEA-DF11-9275-E0CB4E553669.root',
      # '/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0007/F81FBFC1-CCEA-DF11-AB52-E0CB4E1A119A.root',
      # '/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0007/F6318BF2-CEEA-DF11-A00F-E0CB4E553663.root',
    'file:/tmp/testAODData423.root'
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

