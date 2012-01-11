import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingEfficiencyFromCosmics")

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
process.GlobalTag.globaltag = 'GR_R_42_V14::All'

# process.load("MagneticField.Engine.uniformMagneticField_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.hltPixelTrackerHVOn = cms.EDFilter( "DetectorStateFilter",
                                              DetectorType = cms.untracked.string( "pixel" )
                                            )
process.hltStripTrackerHVOn = cms.EDFilter( "DetectorStateFilter",
                                              DetectorType = cms.untracked.string( "sistrip" )
                                            )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( *(
    FILELIST
    ) )
)

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('TrackingEfficiencyFromCosmics_NUMBER.root')
                                 )

process.demo = cms.EDAnalyzer('TrackingEfficiencyFromCosmics',
                              EffDxyMin = cms.double(0),
                              EffDxyMax = cms.double(100),
                              EffDzMin = cms.double(0),
                              EffDzMax = cms.double(100),
                              EffPtMin = cms.double(0),
                              EffPtMax = cms.double(200),
                              MaxDeltaR = cms.double(1),
                              SimMaxDeltaR = cms.double(1000),
                              DzMinCut = cms.double(0),
                              DzMaxCut = cms.double(1000000),
                              DxyCut = cms.double(50),
                              Chi2Cut = cms.double(1000000),
                              TrackPtCut = cms.double(25),
                              StandAlonePtCut = cms.double(0),
                              HighPurity = cms.bool(True),
                              MatchTwoLegs = cms.bool(True),
                              DeltaDxyCut = cms.double(15), # only if matching two legs
                              DeltaDzCut = cms.double(30),  # only if matching two legs
                              DeltaPtCut = cms.double(1000),  # only if matching two legs
                              DeltaPhiCut = cms.double(1),  # only if matching two legs
                              MinimumValidHits = cms.int32(0),
                              UseMCtruth = cms.bool(False),
                              EffOutputFileName = cms.string("Efficiency_NUMBER.root"),
                              EffCleanedOutputFileName = cms.string("EfficiencyCleaned_NUMBER.root"),
                              GenToStandAloneEffOutputFileName = cms.string("GenToStandAloneEfficiency_NUMBER.root"),
                              GenToTrackEffOutputFileName = cms.string("GenToTrackEfficiency_NUMBER.root"),
                              RecomputeIP = cms.bool(False),
                              SingleLegMuon = cms.bool(True),
                              # MuonCollection = cms.InputTag("standAloneMuons"),
                              MuonCollection = cms.InputTag("cosmicMuons1Leg"),
                              TrackCollection = cms.InputTag("generalTracks"),
                              UseAllTracks = cms.bool(False),
                              UseTrackParameters = cms.bool(False),
                              DxyErrorCut = cms.bool(False),
                              DzErrorCut = cms.bool(False),
                              DxyCutForNoDzCut = cms.double(4),
                              PhiRegion = cms.bool(False),
                              PhiMinCut = cms.double(-3.2),
                              PhiMaxCut = cms.double(0),
                              CountOppoSide = cms.bool(True),
                              CountSameSide = cms.bool(True)
                              )

process.p = cms.Path(process.demo)
