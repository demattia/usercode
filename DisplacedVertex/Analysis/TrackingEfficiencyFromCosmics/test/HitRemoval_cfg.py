import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

# process.load("FWCore.MessageService.MessageLogger_cfi")
# initialize MessageLogger and output report
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger.cerr.threshold = 'INFO'
# process.MessageLogger.categories.append('Demo')
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)
# )
# process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/CMSSW_4_2_2/src/reco_RAW2DIGI_L1Reco_RECO_DQM.root'
    )
)

#####################################################################
## BeamSpot from database (i.e. GlobalTag), needed for Refitter
#####################################################################
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Only needed if using trajectories
import RecoTracker.TrackProducer.TrackRefitters_cff
process.TrackRefitter1 = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone()
process.TrackRefitter1.src = 'generalTracks'

# parameters for HitRemoval
process.load("Analysis.TrackingEfficiencyFromCosmics.TrackerTrackHitFilterMod_cff")
process.TrackerTrackHitFilterMod.src = 'TrackRefitter1'
process.TrackerTrackHitFilterMod.commands = cms.vstring(# "keep PXB","keep PXE",
                                                        "drop PXB", "drop PXE",
                                                        # "keep TIB",
                                                        # "drop TID", "drop TIB", "drop TOB", "drop TEC")
                                                        "keep TIB", "keep TID","keep TOB","keep TEC", "drop TIB 1")
process.TrackerTrackHitFilterMod.minimumHits = 4 # 6
# process.TrackerTrackHitFilterMod.replaceWithInactiveHits = False # True
process.TrackerTrackHitFilterMod.replaceWithInactiveHits = True
process.TrackerTrackHitFilterMod.stripAllInvalidHits = False
process.TrackerTrackHitFilterMod.rejectBadStoNHits = False # True
process.TrackerTrackHitFilterMod.StoNcommands = cms.vstring("ALL 18.0")
process.TrackerTrackHitFilterMod.useTrajectories= False # True
process.TrackerTrackHitFilterMod.rejectLowAngleHits= False # True
process.TrackerTrackHitFilterMod.TrackAngleCut= 0. # 0.17 #~20 degrees (minimum angle)
process.TrackerTrackHitFilterMod.usePixelQualityFlag= False # True
process.TrackerTrackHitFilterMod.SaveOnlyAffectedTracks = cms.untracked.bool(False)
# process.TrackerTrackHitFilterMod.PxlCorrClusterChargeCut=10000.0


# # track producer to be run after the hit filter
# import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff
# process.ctfProducerCustomised = RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff.ctfWithMaterialTracks.clone()
# process.ctfProducerCustomised.src = 'TrackerTrackHitFilterMod'
# ##process.ctfProducerCustomised.beamspot='offlineBeamSpot'
# process.ctfProducerCustomised.TTRHBuilder = 'WithAngleAndTemplate'
# process.ctfProducerCustomised.TrajectoryInEvent = False # True

# track producer to be run after the hit filter
import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff
process.generalTracks = RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff.ctfWithMaterialTracks.clone()
process.generalTracks.src = 'TrackerTrackHitFilterMod'
##process.generalTracks.beamspot='offlineBeamSpot'
process.generalTracks.TTRHBuilder = 'WithAngleAndTemplate'
process.generalTracks.TrajectoryInEvent = False # True

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V14::All'

# Output
process.load('Configuration.EventContent.EventContent_cff')
# process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    # outputCommands = cms.untracked.vstring('keep recoTracks_*_*_*',
    #                                        'keep recoTrackExtras_*_*_*',
    #                                        'keep TrackingRecHitsOwned_*_*_*'),
    outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('hitRemoval.root'),
    # dataset = cms.untracked.PSet(
    #     filterName = cms.untracked.string(''),
    #     dataTier = cms.untracked.string('RECO')
    # )
)

# process.RECOoutput_step = cms.EndPath(process.RECOoutput)

process.p = cms.Path(process.offlineBeamSpot+
                     process.TrackRefitter1+
                     process.TrackerTrackHitFilterMod+
                     # process.ctfProducerCustomised
                     process.generalTracks)

process.ep = cms.EndPath(process.RECOoutput)
