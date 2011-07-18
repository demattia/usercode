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
process.GlobalTag.globaltag = 'START42_V11::All'

# process.load("MagneticField.Engine.uniformMagneticField_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/home/demattia/Simulation/WithGENInfo/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_1.root',
    'file:/home/demattia/Simulation/CosmicSimulation/Second/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_2.root',
    'file:/home/demattia/Simulation/CosmicSimulation/Third/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_3.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_4.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_5.root',
    'file:/home/demattia/Simulation/CosmicSimulation/12/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/13/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/14/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/15/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/16/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/17/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/18/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/19/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'file:/home/demattia/Simulation/CosmicSimulation/20/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    # 'file:/home/demattia/Simulation/CosmicSimulation/21/reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_25.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_26.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_27.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_28.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_29.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_30.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_31.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_32.root',
    )
)

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('TrackingEfficiencyFromCosmics.root')
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
                              DzCut = cms.double(30),
                              DxyCut = cms.double(50),
                              Chi2Cut = cms.double(1000000),
                              TrackPtCut = cms.double(25),
                              StandAlonePtCut = cms.double(35),
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
                              RecomputeIP = cms.bool(False),
                              SingleLegMuon = cms.bool(True),
                              # MuonCollection = cms.InputTag("standAloneMuons"),
                              MuonCollection = cms.InputTag("cosmicMuons1Leg"),
                              TrackCollection = cms.InputTag("generalTracks"),
                              UseAllTracks = cms.bool(False)
                              )


process.p = cms.Path(process.demo)
