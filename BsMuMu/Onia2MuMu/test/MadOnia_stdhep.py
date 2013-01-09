import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/MixingInitialLumiPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('L1TriggerConfig/L1GtConfigProducers/Luminosity/lumi1030.L1Menu2008_2E30_Unprescaled_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('HLTrigger/Configuration/HLT_2E30_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.2 $'),
            annotation = cms.untracked.string('GamGamMM.cfi nevts:10'),
            name = cms.untracked.string('PyReleaseValidation')
        )
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(5000)
        )
process.options = cms.untracked.PSet(
        Rethrow = cms.untracked.vstring('ProductNotFound')
        )
process.genParticles.abortOnUnknownPDGCode = False 

# Input source
process.source = cms.Source("MCFileSource",
                                useExtendedAscii = cms.untracked.bool(False),
                                fileNames = cms.untracked.vstring(
    'file:/tmp/aafke/part1.dat',
    'file:/tmp/aafke/part2.dat',
    'file:/tmp/aafke/part3.dat'
    )
                            )
process.MuonAnalysis = cms.EDAnalyzer("Onia2MuMu",
    OutputFileName       = cms.string('/tmp/aafke/onia2mumu_upsmumu.root'),
    OniaType             = cms.int32(553),
    DebugLevel           = cms.int32(0),
    genParticlesLabel    = cms.InputTag("genParticles"),
    #StandAloneMuonsLabel = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    StandAloneMuonsLabel = cms.InputTag("standAloneMuons"),
    GlobalMuonsLabel     = cms.InputTag("globalMuons"),
    MuonsLabel           = cms.InputTag("muons"),
    CaloMuonsLabel       = cms.InputTag("calomuons"),
    BeamSpotLabel        = cms.InputTag("offlineBeamSpot"),
    PrimaryVertexLabel   = cms.InputTag("offlinePrimaryVertices"),
    TrackLabel           = cms.InputTag("generalTracks"),
    triggerEventLabel    = cms.InputTag("hltTriggerSummaryAOD::HLT"),
    L1GTReadoutRec       = cms.InputTag("gtDigis"),
    HLTriggerResults     = cms.InputTag("TriggerResults::HLT"),
    L1MuonLabel          = cms.InputTag("hltL1extraParticles"),
    L2MuonLabel          = cms.InputTag("hltL2Muons:UpdatedAtVtx"),
    L3MuonLabel          = cms.InputTag("hltL3Muons"),
    PATMuonsLabel        = cms.InputTag("selectedLayer1Muons"),
    StoreGenFlag         = cms.bool(True),
    StoreHLTFlag         = cms.bool(False),
    StoreL1Flag          = cms.bool(False),
    StoreTrkFlag         = cms.bool(False),
    StoreSTAMuonFlag     = cms.bool(False),
    StoreGLBMuonFlag     = cms.bool(False),
    StoreAllMuonFlag     = cms.bool(False),
    StoreBeamSpotFlag    = cms.bool(False),
    StorePriVtxFlag      = cms.bool(False),
    StoreOniaFlag        = cms.bool(False),
    StoreOniaRadiation   = cms.bool(True),
    UsingBeamSpot        = cms.bool(False),
    minimumFlag          = cms.bool(False),
    UsingAOD             = cms.bool(True),
    StorePATFlag         = cms.bool(False)
)


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
                                      outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
                                      fileName = cms.untracked.string('file:/tmp/aafke/madonia.10tev.GEN.root'),
                                      dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-HLTDEBUG'),
                    filterName = cms.untracked.string('')
                ),
                                      SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('')
                )
                                  )

# Other statements
process.GlobalTag.globaltag = 'IDEAL_V9::All'

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen*process.MuonAnalysis)
process.out_step = cms.EndPath(process.output)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.out_step)
process.schedule = cms.Schedule(process.generation_step)
print process.dumpPython()
