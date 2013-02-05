import FWCore.ParameterSet.Config as cms

from YZheng.UpsilonAna.selection_cff import SELECTIONCUT

filterUpsilon = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    cut = cms.string(SELECTIONCUT + " && 9.3 < mass < 9.6"),
    filter = cms.bool(True),
    variables = cms.PSet(),
    flags = cms.PSet(),
)

from YZheng.UpsilonAna.selection_cff import detailedDimuonTree as detailedUpsilonTree

eventDisplayUpsilonPath = cms.Path(filterUpsilon * detailedUpsilonTree)

outUpsilon = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("eventDisplayUpsilon.root"),
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventDisplayUpsilonPath'))
)

endUpsilon = cms.EndPath(outUpsilon)

