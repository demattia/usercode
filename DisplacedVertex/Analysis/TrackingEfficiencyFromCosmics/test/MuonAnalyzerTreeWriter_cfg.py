import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonAnalyzerTreeWriter")

# # initialize MessageLogger and output report
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# # process.MessageLogger.cerr.threshold = 'INFO'
# # process.MessageLogger.categories.append('Demo')
# # process.MessageLogger.cerr.INFO = cms.untracked.PSet(
# #     limit = cms.untracked.int32(-1)
# # )

process.load("FWCore.MessageService.test.Services_cff")
# Here is the configuration of the MessgeLogger Service:
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('Message'),
    Message = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    )
)

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
# process.load('Configuration.StandardSequences.MagneticField_cff')

# Careful, this needs to be changed for the data
process.GlobalTag.globaltag = 'START42_V11::All'

# process.load("MagneticField.Engine.uniformMagneticField_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_11.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_12.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_13.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_14.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_15.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_16.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_17.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_18.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_19.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_20.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_21.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_22.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_23.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_24.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_25.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_26.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_27.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_28.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_29.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_30.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_31.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_32.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_33.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_34.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_35.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_36.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_37.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_38.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_39.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_40.root",
    # "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_41.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_42.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_43.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_44.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_45.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_46.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_47.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_48.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_49.root",
    "castor:/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt10/Reco/reco_50.root"
    )
)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('MuonAnalyzerTreeWriter.root')
                                 )

process.MaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
                                            MaxDPhi = cms.double(1.6),
                                            ComponentName = cms.string('PropagatorWithMaterial'),
                                            Mass = cms.double(0.105),
                                            PropagationDirection = cms.string('alongMomentum'),
                                            useRungeKutta = cms.bool(False),
                                            # If ptMin > 0, uncertainty in reconstructed momentum will be taken into account when estimating rms scattering angle.
                                            # (By default, it is neglected). However, it will also be assumed that the track pt can't be below specified value,
                                            # to prevent this scattering angle becoming too big.
                                            ptMin = cms.double(-1.)
                                            )

process.SteppingHelixPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
                                                    ComponentName = cms.string('SteppingHelixPropagatorAny'),
                                                    NoErrorPropagation = cms.bool(False),
                                                    PropagationDirection = cms.string('anyDirection'),
                                                    useTuningForL2Speed = cms.bool(False),
                                                    useIsYokeFlag = cms.bool(True),
                                                    endcapShiftInZNeg = cms.double(0.0),
                                                    SetVBFPointer = cms.bool(False),
                                                    AssumeNoMaterial = cms.bool(False),
                                                    endcapShiftInZPos = cms.double(0.0),
                                                    useInTeslaFromMagField = cms.bool(False),
                                                    VBFName = cms.string('VolumeBasedMagneticField'),
                                                    useEndcapShiftsInZ = cms.bool(False),
                                                    sendLogWarning = cms.bool(False),
                                                    useMatVolumes = cms.bool(True),
                                                    debug = cms.bool(False),
                                                    #This sort of works but assumes a measurement at propagation origin
                                                    ApplyRadX0Correction = cms.bool(True),
                                                    useMagVolumes = cms.bool(True),
                                                    returnTangentPlane = cms.bool(True)
                                                    )

# Smart propagator with IP
process.smartPropagatorWithIPESProducer = cms.ESProducer("SmartPropagatorWithIPESProducer",
                                                         ComponentName = cms.string('SmartPropagatorWithIP'),
                                                         TrackerPropagator = cms.string('PropagatorWithMaterial'),
                                                         # TrackerPropagator = cms.string('RungeKuttaTrackerPropagator'),
                                                         MuonPropagator = cms.string('SteppingHelixPropagatorAny'),
                                                         PropagationDirection = cms.string('alongMomentum'),
                                                         Epsilon = cms.double(10.0) # the standard one has 5., but uses 10 hardcoded internally...
                                                         )



process.demo = cms.EDAnalyzer('MuonAnalyzerTreeWriter',
                              EffDxyMin = cms.double(0),
                              EffDxyMax = cms.double(100),
                              EffDzMin = cms.double(0),
                              EffDzMax = cms.double(100),
                              EffPtMin = cms.double(0),
                              EffPtMax = cms.double(200),
                              
                              MaxDeltaR = cms.double(1),
                              SimMaxDeltaR = cms.double(1000),
                              DzMinCut = cms.double(0),
                              DzMaxCut = cms.double(10),
                              DxyCut = cms.double(1000000),
                              Chi2Cut = cms.double(1000000),
                              TrackPtCut = cms.double(0),
                              StandAlonePtCut = cms.double(0),
                              HighPurity = cms.bool(True),
                              
                              MatchTwoLegs = cms.bool(True),

                              DeltaDxyCut = cms.double(15), # only if matching two legs
                              DeltaDzCut = cms.double(30),  # only if matching two legs
                              DeltaPtCut = cms.double(1000),  # only if matching two legs
                              DeltaPhiCut = cms.double(1),  # only if matching two legs
                              MinimumValidHits = cms.int32(0),

                              UseMCtruth = cms.bool(True),

                              EffOutputFileName = cms.string("Efficiency.root"),
                              EffCleanedOutputFileName = cms.string("EfficiencyCleaned.root"),
                              GenToStandAloneEffOutputFileName = cms.string("GenToStandAloneEfficiency.root"),
                              GenToTrackEffOutputFileName = cms.string("GenToTrackEfficiency.root"),

                              RecomputeIP = cms.bool(True),

                              SingleLegMuon = cms.bool(False),
                              CountSameSide = cms.bool(True),
                              CountOppoSide = cms.bool(False),
                              PhiRegion = cms.bool(False),
                              PhiMinCut = cms.double(-3.2),
                              PhiMaxCut = cms.double(3.2),

                              MuonCollection = cms.InputTag("standAloneMuons"),
                              # MuonCollection = cms.InputTag("cosmicMuons1Leg"),
                              TrackCollection = cms.InputTag("generalTracks"),

                              UseAllTracks = cms.bool(True), # Do not apply analysis quality cuts if True
                              UseTrackParameters = cms.bool(False), # Use track parameters for matched standalone
                              DxyErrorCut = cms.bool(False),
                              DzErrorCut = cms.bool(False),
                              DxyCutForNoDzCut = cms.double(4),

                              GenInsideTkVol = cms.bool(True)
                              )

process.p = cms.Path(process.demo)
