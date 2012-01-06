import FWCore.ParameterSet.Config as cms

isoTrackPrefilter = cms.EDFilter("IsoTrackPrefilter",
                                 trackSrc = cms.InputTag("trackSel"),
                                 muonSrc = cms.InputTag("standAloneMuons"),
                                 generatorSrc = cms.InputTag("genParticles"),
                                 minPt = cms.double(25),
                                 coneSize = cms.double(0.3),
                                 isolationThreshold = cms.double(99999.0),
                                 signalPDGId = cms.int32(36))
