# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms

process = cms.Process("Onia2MuMuPAT")

from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import *

onia2MuMuPAT(process, GlobalTag="START39_V8::All", MC=True, HLT="HLT", Filter=False)

process.source.fileNames = cms.untracked.vstring(
       'file:/tmp/testJpsiMuMuRECO.root'
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

