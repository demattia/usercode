import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'MC_31X_V3::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("Pythia6PtGun",
    PGunParameters = cms.PSet(
        MinPhi = cms.double(-3.14159265359),
        MinPt = cms.double(0.0),
        ParticleID = cms.vint32(553),
        MaxRapidity = cms.double(2.0),
        MaxPhi = cms.double(3.14159265359),
        MinRapidity = cms.double(-2.0),
        AddAntiParticle = cms.bool(False),
        MaxPt = cms.double(20.0)
    ),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
        pythiaUpsilonDecays = cms.vstring(
            'MDME(1034,1)= 0.025200    ! e- e+',
            'MDME(1035,1)= 0.024800    ! mu- mu+',
            'MDME(1036,1)= 0.026700    ! tau- tau+',
            'MDME(1037,1)= 0.015000    ! d dbar',
            'MDME(1038,1)= 0.045000    ! u ubar',
            'MDME(1039,1)= 0.015000    ! s sbar',
            'MDME(1040,1)= 0.045000    ! c cbar',
            'MDME(1041,1)= 0.774300    ! g g g',
            'MDME(1042,1)= 0.029000    ! gamma g'),
        parameterSets = cms.vstring('pythiaUpsilonDecays')
    )
)
process.ProductionFilterSequence = cms.Sequence(process.generator)

# Dumps events to tree
process.DT = cms.EDAnalyzer("DimuonTree",
  OutputFile = cms.string("upsilonGun.root")
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction * process.DT)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reconstruction_step)

# special treatment in case of production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq

