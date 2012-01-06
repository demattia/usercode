import FWCore.ParameterSet.Config as cms

pseudoLeptonProducer = cms.EDProducer("PseudoLeptonProducer",
                                      trackSrc = cms.InputTag("trackSel"),
                                      muonSrc = cms.InputTag("standAloneMuons"),
                                      minPt = cms.double(25))

