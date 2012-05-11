# Auto generated configuration file
# using: 
# Revision: 1.372.2.3 
# Source: /local/reps/CMSSW.admin/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: ZMM_8TeV_cfi --conditions auto:startup -s GEN,SIM --datatier GEN-SIM -n 10 --relval 18000,200 --eventcontent RAWSIM --no_exec --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.RandomNumberGeneratorService.generator = cms.PSet(
    initialSeed = cms.untracked.uint32(RANDOMSEED),
    engineName = cms.untracked.string('HepJamesRandom')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.372.2.3 $'),
    annotation = cms.untracked.string('ZMM_8TeV_cfi nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:step1.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START52_V9::All'

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
                'MDME(1042,1)=0',
                'MDME( 174,1) = 0    !Z decay into d dbar',
                'MDME( 175,1) = 0    !Z decay into u ubar',
                'MDME( 176,1) = 0    !Z decay into s sbar',
                'MDME( 177,1) = 0    !Z decay into c cbar',
                'MDME( 178,1) = 0    !Z decay into b bbar',
                'MDME( 179,1) = 0    !Z decay into t tbar',
                'MDME( 182,1) = 1    !Z decay into e- e+',
                'MDME( 183,1) = 0    !Z decay into nu_e nu_ebar',
                'MDME( 184,1) = 1    !Z decay into mu- mu+',
                'MDME( 185,1) = 0    !Z decay into nu_mu nu_mubar',
                'MDME( 186,1) = 0    !Z decay into tau- tau+',
                'MDME( 187,1) = 0    !Z decay into nu_tau nu_taubar'
                ),
        SLHAParameters = cms.vstring(
       'SLHAFILE = slha_Z.out'),
        parameterSets = cms.vstring(
            'pythiaUESettings',
        'SLHAParameters'
        )
    )
)

process.mumugenfilter = cms.EDFilter("MCParticlePairFilter",
    Status = cms.untracked.vint32(1, 1),
    MinPt = cms.untracked.vdouble(2.5, 2.5),
    MaxEta = cms.untracked.vdouble(2.5, 2.5),
    MinEta = cms.untracked.vdouble(-2.5, -2.5),
    ParticleCharge = cms.untracked.int32(-1),
    ParticleID1 = cms.untracked.vint32(13),
    ParticleID2 = cms.untracked.vint32(13)
)


process.ProductionFilterSequence = cms.Sequence(process.generator+process.mumugenfilter)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

