import FWCore.ParameterSet.Config as cms

onia2MuMuPAT = cms.EDProducer('Onia2MuMuPAT',
  muons = cms.InputTag("patMuons"),
  beamSpotTag = cms.InputTag("offlineBeamSpot"),
  primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
  higherPuritySelection = cms.string("isGlobalMuon"), ## At least one muon must pass this selection
  lowerPuritySelection  = cms.string("isGlobalMuon"), ## BOTH muons must pass this selection
  dimuonSelection  = cms.string("2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25"), ## The dimuon must pass this selection before vertexing
  addCommonVertex = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
  addMuonlessPrimaryVertex = cms.bool(True), ## Embed the primary vertex re-made from all the tracks except the two muons
  addMCTruth = cms.bool(True),      ## Add the common MC mother of the two muons, if any
  resolvePileUpAmbiguity = cms.bool(True),   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt
  # Apply any preselection cut. These cuts are applied just before the candidates are saved in the collection, at the end of the onia2MuMuPAT.
  addThirdTrack = cms.bool(True),  ## search a third track making a good vertex with the dimuon
  minTrackPt =cms.double(0.5), ## minimum pt of the third track
  trackMass = cms.double(0.4936), ## mass for the track
  diMuMassMin= cms.double(2.8),
  diMuMassMax= cms.double(3.4),
  diMuPlusTrackMassMin = cms.double(4.9),
  diMuPlusTrackMassMax = cms.double(5.9),
  preselection = cms.string("")
)
