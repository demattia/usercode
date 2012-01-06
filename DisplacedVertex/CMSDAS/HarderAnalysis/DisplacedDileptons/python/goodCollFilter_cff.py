# extracted from DPGAnalysis/Skims/python/GOODCOLL_filter_cfg.py in CMSSW_3_6_3
# however, had to remove technical trigger bit arithmetic because the information is not in PAT.
# also the requirement for bits 40 and 41 that was in there is now considered bad

import FWCore.ParameterSet.Config as cms


primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
   filter = cms.bool(True)
)


noscraping = cms.EDFilter("FilterOutScraping",
applyfilter = cms.untracked.bool(True),
debugOn = cms.untracked.bool(False),
numtrack = cms.untracked.uint32(10),
thresh = cms.untracked.double(0.25)
)


goodcollFilter=cms.Sequence(primaryVertexFilter+noscraping)
