from Configuration.StandardSequences.Reconstruction_cff import *

# NON-STANDARD RECONSTRUCTION: ELECTRONS
pfTrackElec.useFifthStepForTrackerDrivenGsf = True

# NON-STANDARD RECONSTRUCTION: MUONS
muons.inputCollectionLabels = cms.VInputTag(cms.InputTag("generalTracks"),
                                            cms.InputTag("globalMuons"), 
                                            cms.InputTag("standAloneMuons"))
globalMuons.MuonCollectionLabel = cms.InputTag("standAloneMuons")
globalSETMuons.MuonCollectionLabel = cms.InputTag("standAloneSETMuons")
