import FWCore.ParameterSet.Config as cms

from YZheng.UpsilonAna.selection_cff import SELECTIONCUT

filterJpsi = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    cut = cms.string(SELECTIONCUT + " && 3.0 < mass < 3.2"),
    filter = cms.bool(True),
    variables = cms.PSet(),
    flags = cms.PSet(),
)

from YZheng.UpsilonAna.selection_cff import detailedDimuonTree as detailedJpsiTree

eventDisplayJpsiPath = cms.Path(filterJpsi * detailedJpsiTree)

outJpsi = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("eventDisplayJpsi.root"),
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventDisplayJpsiPath'))
)

endJpsi = cms.EndPath(outJpsi)

