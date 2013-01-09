import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("PhysicsTools.HepMCCandAlgos.allMuonsGenParticlesMatch_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

# Input source
process.source = cms.Source("PoolSource",
                            
    fileNames = cms.untracked.vstring(
    # first data
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/68ED95B5-50E2-DE11-B4C8-001D09F27003.root',
    )
)

process.GlobalTag.globaltag = 'GR09_P_V7::All'

process.MuonAnalysis = cms.EDAnalyzer("Onia2MuMu",
    OutputFileName       = cms.string('jpsi_data09.root'),
    OniaType             = cms.int32(443),
    OniaMaxCat           = cms.int32(6),
    skimOnOniaMaxCat     = cms.bool(False),
    DebugLevel           = cms.int32(0),
    genParticlesLabel    = cms.InputTag("genParticles"),
    # StandAloneMuonsLabel = cms.InputTag("standAloneMuons"),
    GlobalMuonsLabel     = cms.InputTag("globalMuons"),
    MuonsLabel           = cms.InputTag("muons"),
    CaloMuonsLabel       = cms.InputTag("calomuons"),
    BeamSpotLabel        = cms.InputTag("offlineBeamSpot"),
    PrimaryVertexLabel   = cms.InputTag("offlinePrimaryVerticesWithBS"),
    UsePrimaryNoMuons    = cms.bool(True),
    TrackLabel           = cms.InputTag("generalTracks"),
    PhotonLabel          = cms.InputTag("particleFlow"),
    PhotonMinEnergy      = cms.double(2.0),
    triggerEventLabel    = cms.string("hltTriggerSummaryAOD"),
    triggerResultsLabel  = cms.string("TriggerResults"),
    HLTprocessName8e29   = cms.string("HLT"),
    HLTprocessName1e31   = cms.string("HLT"),  
    L1GTReadoutRec       = cms.InputTag("gtDigis"),
    L1MuonLabel          = cms.InputTag("l1extraParticles"),
    PATMuonsLabel        = cms.InputTag("selectedLayer1Muons"),
    # PATPhotonsLabel      = cms.InputTag("selectedLayer1Photons"),
    useKinFit            = cms.bool(False),
    StoreGenFlag         = cms.bool(False),
    StoreHLTFlag         = cms.bool(True),
    StoreL1Flag          = cms.bool(True),
    StoreTrkFlag         = cms.bool(True),
    StorePhotonFlag      = cms.bool(False),
    # StorePFMuonFlag      = cms.bool(True),
    StoreTRKMuonFlag     = cms.bool(True),
    StoreGLBMuonFlag     = cms.bool(True),
    StoreCALMuonFlag     = cms.bool(True),
    StoreBeamSpotFlag    = cms.bool(True),
    StorePriVtxFlag      = cms.bool(True),
    StoreOniaFlag        = cms.bool(True),
    StoreChicFlag        = cms.bool(False),
    StoreBpFlag          = cms.bool(False),                     
    StoreWSOnia          = cms.bool(True),                                  
    StoreOniaRadiation   = cms.bool(False),
    UsingBeamSpot        = cms.bool(False),
    minimumFlag          = cms.bool(False),
    UsingAOD             = cms.bool(False),
    StorePATFlag         = cms.bool(False)
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# different choices of trigger conditions:
# bsc minbias
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41)')
# bsc minbias and veto on beam halo
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
# bsc minbias in coincidence with bptx
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41)')
# bsc minbias in coinidence with bptx and veto on beam halo
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')


# this is for filtering on HLT physics bit
process.hltPhysicsDeclared = cms.EDFilter("HLTHighLevel",
                                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                 HLTPaths = cms.vstring("HLT_PhysicsDeclared"
                                                        ),
                                 eventSetupPathsKey = cms.string(''),
                                 andOr = cms.bool(True),
                                 throw = cms.bool(True)
                                 )


# Beamspot temporary fix
from CondCore.DBCommon.CondDBSetup_cfi import *
process.firstCollBeamspot = cms.ESSource("PoolDBESSource",
                                         CondDBSetup,
                                         connect = cms.string("frontier://PromptProd/CMS_COND_31X_BEAMSPOT"),
                                         # connect = cms.string("sqlite_file:/afs/cern.ch/user/y/yumiceva/public/BeamSpotObjects_2009_v1_express@a35f2218-e25f-11de-9d9b-0018f34695d6.db"),
                                         toGet = cms.VPSet(cms.PSet(record = cms.string("BeamSpotObjectsRcd"),
                                                                    tag =cms.string("BeamSpotObjects_2009_v3_offline"))
                                                           )
                                         )
process.es_prefer_firstCollBeamspot = cms.ESPrefer("PoolDBESSource","firstCollBeamspot")
process.load('RecoVertex.BeamSpotProducer.BeamSpot_cfi')

# filter on lumisections
from HeavyFlavorAnalysis.Onia2MuMu.goodLumiSectionList_cfi import *
process.source.lumisToProcess = goodLumisToProcess

process.ana_step = cms.Path(process.hltLevel1GTSeed + process.hltPhysicsDeclared + process.offlineBeamSpot + process.MuonAnalysis)
#process.ana_step = cms.EndPath(process.MuonAnalysis)

# Schedule definition
process.schedule = cms.Schedule(process.ana_step)

