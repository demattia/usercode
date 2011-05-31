import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECOALIGN")
#process.load("MuonAnalysis.MomentumScaleCalibration.Summer08_Upsilon1S_cff")

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(

#   "rfio:/castor/cern.ch/user/c/castello/3E23AEB9-3ED0-DF11-B19A-001A4BA6E1CA.root",
#   "rfio:/castor/cern.ch/user/c/castello/E25D7CCB-D3EA-DF11-8F73-E0CB4E29C507.root",
    "rfio:/castor/cern.ch/user/c/castello/180A4793-D2EA-DF11-B061-E0CB4E553643.root"
    
    )
)

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")
#process.load("Alignment.CommonAlignmentProducer.GlobalPosition_Fake_cff")## originally inside, but NOT needed if GlobaTag is used!
#process.load("MagneticField.Engine.uniformMagneticField_cfi") # use the standard one

## ######### NEEDED for Re-Reco ####################################
##### needed for Refitting/Reconstruction #######################
process.load('Configuration/StandardSequences/MagneticField_cff')
#process.load("Configuration/StandardSequences/MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
#################################################################

from CondCore.DBCommon.CondDBSetup_cfi import *
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

############  DATABASE conditions ###############################
#process.GlobalTag.connect="frontier://FrontierProd/CMS_COND_31X"
process.GlobalTag.globaltag ='MC_38Y_V14::All' #'MC_3XY_V27::All' # MC
#process.GlobalTag.globaltag ='FT_R_38X_V14A::All' # DATA --> FT_R_38X_V14A

from CondCore.DBCommon.CondDBSetup_cfi import *


#### Lorentz Angle ############################################################

## process.poolDBESSource1 = cms.ESSource("PoolDBESSource",
##                                        BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
##                                        DBParameters = cms.PSet(messageLevel = cms.untracked.int32(2),
##                                                                authenticationPath = cms.untracked.string('/path/to/authentication') ),
##                                        timetype = cms.untracked.string('runnumber'),
##                                        connect = cms.string('frontier://PromptProd/CMS_COND_31X_STRIP'),
##                                        toGet = cms.VPSet(cms.PSet(record = cms.string('SiStripLorentzAngleRcd'),
##                                                                   tag = cms.string('SiStripLorentzAngle_GR10_v2_offline') )))
## process.es_prefer_LA = cms.ESPrefer('PoolDBESSource','poolDBESSource1')


## #### Backplane corrections ##################################################################

## process.poolDBESSource2 = cms.ESSource("PoolDBESSource",
##                                        BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
##                                        DBParameters = cms.PSet(messageLevel = cms.untracked.int32(2),
##                                                                authenticationPath = cms.untracked.string('/path/to/authentication') ),
##                                        timetype = cms.untracked.string('runnumber'),connect = cms.string('frontier://PromptProd/CMS_COND_31X_STRIP'),
##                                        toGet = cms.VPSet(cms.PSet(record = cms.string('SiStripConfObjectRcd'),
##                                                                   tag = cms.string('SiStripShiftAndCrosstalk_GR10_v2_offline') )))
## process.es_prefer_BP = cms.ESPrefer('PoolDBESSource','poolDBESSource2')

###  Tracker geometry ##########################################################

from CalibTracker.Configuration.Common.PoolDBESSource_cfi import poolDBESSource
import CalibTracker.Configuration.Common.PoolDBESSource_cfi
process.trackerAlignment =  CalibTracker.Configuration.Common.PoolDBESSource_cfi.poolDBESSource.clone(
                                        connect = cms.string('sqlite_file:IdealPlusSaggitta_xTOByTIB100mu.db'),
                                        timetype = cms.string("runnumber"),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                                                                   tag = cms.string('Alignments')
                                                                   ))
                                        )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")

## process.trackerAlignment = cms.ESSource("PoolDBESSource",
##                                         CondDBSetup,
##                                         toGet = cms.VPSet(cms.PSet(
##                                                 record = cms.string('TrackerAlignmentRcd'),
##                                                 tag = cms.string('TrackerAlignment_GR10_v3_offline')
##                                                 )),
##                                         connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT')                              
##                                         )
## process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")

###  Tracker APE  ################################################################
## process.trackerAPE = cms.ESSource("PoolDBESSource",
##                                         CondDBSetup,
##                                         toGet = cms.VPSet(cms.PSet(
##                                                 record = cms.string('TrackerAlignmentErrorRcd'),
##                                                 tag = cms.string('TrackerAlignmentErrors_GR10_v2_offline')
##                                                 )),
##                                         connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT')                              
##                                         )
## process.es_prefer_trackerAPE = cms.ESPrefer("PoolDBESSource", "trackerAPE")


#### Global position  #############################################################

## process.GlobalPosition = cms.ESSource("PoolDBESSource", CondDBSetup,
##                       connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
##                       toGet= cms.VPSet( cms.PSet( record = cms.string('GlobalPositionRcd'),
##                                                   tag = cms.string('GlobalAlignment_v2_offline')
##                                                  )
##                                        )
## )
## process.es_prefer_GlobalPosition= cms.ESPrefer("PoolDBESSource", "GlobalPosition")

## ## DT Muon Geometry  ############################################################
## process.muonAlignment = cms.ESSource("PoolDBESSource", CondDBSetup,
##                       connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
##                       toGet= cms.VPSet( cms.PSet( record = cms.string('DTAlignmentRcd'),
##                                                   tag = cms.string('DTAlignment_2009_v4_offline')
##                                                  )
##                                        )
## )
## process.es_prefer_muonAlignment= cms.ESPrefer("PoolDBESSource","muonAlignment")

## ## CSC Muon Geometry #################################################
## process.muonEndcapAlignment = cms.ESSource("PoolDBESSource", CondDBSetup,
##                       connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
##                       toGet= cms.VPSet( cms.PSet( record = cms.string('CSCAlignmentRcd'),
##                                                   tag = cms.string('CSCAlignment_2009_v4_offline')
##                                                  )
##                                        )
## )
## process.es_prefer_muonEndcapAlignment= cms.ESPrefer("PoolDBESSource","muonEndcapAlignment")

###### Beam spot ######################################################
## process.beamSpot = cms.ESSource("PoolDBESSource",
##                                         CondDBSetup,
##                                         toGet = cms.VPSet(cms.PSet(
##                                                 record = cms.string('BeamSpotObjectsRcd'),
##                                                 tag = cms.string('BeamSpotObjects_2009_LumiBased_v18_offline')
##                                                 )),
##                                         connect = cms.string('frontier://FrontierProd/CMS_COND_31X_BEAMSPOT')                              
##                                         )
## process.es_prefer_beamSpot = cms.ESPrefer("PoolDBESSource", "beamSpot")


process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

######  Tracker refitter ##############################################
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.src ='generalTracks'
process.TrackRefitter.TrajectoryInEvent = True
process.TrackRefitter.TTRHBuilder = "WithAngleAndTemplate" 

######  Global muon reconstruction #####################################
process.load("RecoMuon.GlobalMuonProducer.globalMuons_cff")
process.globalMuons.TrackerCollectionLabel = 'TrackRefitter'

### Muon Reco ################
#from RecoMuon.MuonIdentification.muons_cfi import *
process.muonReReco = process.muons.clone()
process.muonReReco.inputCollectionLabels = ['TrackRefitter', 'globalMuons', 'standAloneMuons:UpdatedAtVtx']

######################### Zmumu TPG skim ###############################################################
#process.load("DPGAnalysis.Skims.WZMuSkim_cff")

######################### Zmumu GoldenSelection sequence ###############################################
process.load("ElectroWeakAnalysis.ZMuMu.ZMuMuGolden_cfi")
#process.goodGlobalMuons.src = 'muons'
process.goodGlobalMuons.src = "muonReReco"
process.goodGlobalMuons.cut = 'isGlobalMuon = 1 & isTrackerMuon = 1 &  pt > 20 & abs(eta)<2.1 & isolationR03().sumPt<3.0 & globalTrack().hitPattern().numberOfValidTrackerHits>10'
process.goodGlobalMuons.filter = cms.bool(True)
process.zmmCands.filter = cms.bool(True)                                                                            
process.dimuonsFilter.filter = cms.bool(True)

process.dimuonsHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI38X")
process.dimuonsHLTFilter.HLTPaths = ["HLT_Mu9","HLT_Mu11","HLT_IsoMu9"]                              ## RunPeriodB_1 146428-147116  ##,"HLT_DoubleMu3"

########## Output module for Reconstruction ############################
process.zmumuSkim=cms.OutputModule(
   "PoolOutputModule",
   outputCommands=cms.untracked.vstring( 
   #"keep *"
   #"keep *_generator_*_*",
   #"keep *_genParticles_*_*",
   #"keep SimTracks_g4SimHits_*_*",
                                         "keep *_muonReReco_*_*",
                                         "keep *_zmmCands_*_*",
                                         "keep *_muons_*_*",
                                         "keep *_standAloneMuons_*_*",
                                         "keep *_generalTracks_*_*",
                                         "keep *_globalMuons_*_*",
                                         "keep *_TriggerResults_*_*",
                                         "keep *_TrackRefitter_*_*"
                                       ),

   SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    ),
   
   fileName=cms.untracked.string("Zmumu_IdealPlusSaggitta_xTOByTIB100mu.root")
   )

process.p = cms.Path(
    process.offlineBeamSpot*
    process.TrackRefitter*process.globalMuons*process.muonReReco*
    process.ewkZMuMuGoldenSequence
    )
#process.p = cms.Sequence(process.offlineBeamSpot*process.TrackRefitter*process.globalMuons*process.muonReReco)
process.outpath = cms.EndPath(process.zmumuSkim)

#########################################################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3000)
)
