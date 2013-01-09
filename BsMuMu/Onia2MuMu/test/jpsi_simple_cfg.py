import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
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
       ### test
      #'/store/relval/CMSSW_3_1_2/RelValZMM/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/EECE7AB6-CC78-DE11-805C-0019B9F709A4.root'
 #   #2.36 GeV files:
    'file:/tmp/hwoehri/0E8FD2DA-76DB-DE11-8657-001F296BC5BA.root',
    'file:/tmp/hwoehri/0CE5AD13-73DB-DE11-8B4F-00237DA15C96.root'
    )
)

process.GlobalTag.globaltag = 'DESIGN_3X_V8A::All'

process.MuonAnalysis = cms.EDAnalyzer("Onia2MuMu",
    OutputFileName       = cms.string('JPsi900GeV_test.root'),
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
    PrimaryVertexLabel   = cms.InputTag("offlinePrimaryVertices"),
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
    StoreGenFlag         = cms.bool(True),
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
    UsePrimaryNoMuons    = cms.bool(True),
    StoreOniaFlag        = cms.bool(True),
    StoreChicFlag        = cms.bool(False),
    StoreBpFlag          = cms.bool(False),                     
    StoreWSOnia          = cms.bool(True),                                  
    StoreOniaRadiation   = cms.bool(False),
    UsingBeamSpot        = cms.bool(True),
    minimumFlag          = cms.bool(False),
    UsingAOD             = cms.bool(False),
    StorePATFlag         = cms.bool(False)
)

process.ana_step = cms.EndPath(process.MuonAnalysis)

# Schedule definition
process.schedule = cms.Schedule(process.ana_step)

