import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")
# this example configuration offers some minimum
# annotation, to help users get through; please
# don't hesitate to read through the comments
# use MessageLogger to redirect/suppress multiple
# service messages coming from the system
#
# in this config below, we use the replace option to make
# the logger let out messages of severity ERROR (INFO level
# will be suppressed), and we want to limit the number to 10
#
process.load("Configuration.Generator.PythiaUESettings_cfi")
process.load("Configuration.StandardSequences.Generator_cff")
# process.load("Configuration.StandardSequences.VtxSmeared")
process.load("IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.generator = cms.PSet(
    initialSeed = cms.untracked.uint32(93278151),
    engineName = cms.untracked.string('HepJamesRandom')
    )

process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")




# process.load("Configuration.StandardSequences.Mixing")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")




process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
)

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
     ignoreTotal=cms.untracked.int32(0),
     oncePerEventMode = cms.untracked.bool(False)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V15B::All'

# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("Pythia6PtGun",
    maxEventsToPrint = cms.untracked.int32(5),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(True),    
    PGunParameters = cms.PSet(
        ParticleID = cms.vint32(1000023),
                    AddAntiParticle = cms.bool(False),
                    MinPhi = cms.double(-3.14159265359),
                    MaxPhi = cms.double(3.14159265359),
                    MinPt = cms.double(50.0),
                    MaxPt = cms.double(50.0001),
                    MinEta = cms.double(-2.6),
                    MaxEta = cms.double(2.6)
                ),
               PythiaParameters = cms.PSet(
                    pythiaUESettings = cms.vstring(
                'IMSS(1)=11',
                'MSTJ(22)=2     ! Decay those unstable particles',
                'PARJ(71)=100 .  ! for which ctau  100 mm',
                'MDME(1034,1)=0                ! Upsilon -> ee',
                'MDME(1035,1)=1                ! Upsilon -> mumu',
                'MDME(1036,1)=0',
                'MDME(1037,1)=0',
                'MDME(1038,1)=0',
                'MDME(1039,1)=0',
                'MDME(1040,1)=0',
                'MDME(1041,1)=0',
                'MDME(1042,1)=0'
                        ),
        SLHAParameters = cms.vstring(
       'SLHAFILE = slha_Z.out'),
        parameterSets = cms.vstring(
            'pythiaUESettings',
        'SLHAParameters'
        )
    )
)

EvtGenEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_evtgenproducer_*_*')
)


process.FEVT = cms.OutputModule("PoolOutputModule",
                                    process.FEVTSIMEventContent,
                                    fileName = cms.untracked.string('gen_singleNeutralino.root')
                                )

process.options = cms.untracked.PSet(
   Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.pstart = cms.Path(process.generator)
process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
process.p2 = cms.Path(process.pdigi)
process.p3 = cms.Path(process.L1Emulator)
process.p4 = cms.Path(process.DigiToRaw)
process.p5 = cms.Path(process.RawToDigi)
process.p6 = cms.Path(process.reconstruction)
process.outpath = cms.EndPath(process.FEVT)

process.schedule = cms.Schedule(process.pstart,
                                process.p0,
                                process.p1,
                                process.p2,
                                process.p3,
                                process.p4,
                                process.p5,
                                process.p6,

process.outpath
)

process.g4SimHits.Generator.HepMCProductLabel = 'generator'
process.genParticleCandidates.src = 'generator'
process.genParticles.src = 'generator'
process.VtxSmearedCommon.src = 'generator'
process.genParticles.abortOnUnknownPDGCode = False
