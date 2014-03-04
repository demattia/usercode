# This file was taken from the RelVal configuration of 4_2_2: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/GeneratorInterface/CosmicMuonGenerator/python/CMSCGENsource_cfi.py?revision=1.10&view=markup

# Auto generated configuration file
# using:
# Revision: 1.303
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v
# with command line options: UndergroundCosmicMu_cfi.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:GRun,RAW2DIGI,L1Reco -n 10 --geometry DB --conditions auto:startup --relval 666000,7400 --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --eventcontent FEVTDEBUG --scenario cosmics --no_exec --python_filename RelValCosmics_START42_V11.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimNOBEAM_cff')
process.load('Configuration.StandardSequences.DigiCosmics_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Change the seed
process.RandomNumberGeneratorService.generator.initialSeed = 123456790


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100000)
        )

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

    )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.1 $'),
            annotation = cms.untracked.string('UndergroundCosmicMu_cfi.py nevts:10'),
            name = cms.untracked.string('PyReleaseValidation')
        )

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
                                               splitLevel = cms.untracked.int32(0),
                                               eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                               outputCommands = process.FEVTDEBUGEventContent.outputCommands,
                                               fileName = cms.untracked.string('UndergroundCosmicMu_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root'),
                                               dataset = cms.untracked.PSet(
            filterName = cms.untracked.string(''),
                    dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-HLTDEBUG')
                ),
                                               SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('generation_step')
                )
                                           )

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START42_V11::All'

process.cosmicInTracker = cms.EDFilter("CosmicGenFilterHelix",
                                           maxZ = cms.double(212.0),
                                           src = cms.InputTag("generator"),
                                           minPt = cms.double(0.0),
                                           charges = cms.vint32(1, -1),
                                           pdgIds = cms.vint32(-13, 13),
                                           minZ = cms.double(-212.0),
                                           radius = cms.double(80.0),
                                           doMonitor = cms.untracked.bool(False),
                                           minP = cms.double(0.0),
                                           propagator = cms.string('SteppingHelixPropagatorAlong')
                                       )


process.generator = cms.EDProducer("CosMuoGenProducer",
                                       AcptAllMu = cms.bool(False),
                                       ClayWidth = cms.double(50000.0),
                                       MinP = cms.double(10.0),
                                       MTCCHalf = cms.bool(False),
                                       MinPhi = cms.double(0.0),
                                       TIFOnly_constant = cms.bool(False),
                                       MaxP = cms.double(3000.0),
                                       NuProdAlt = cms.double(7500000.0),
                                       RhoAir = cms.double(0.001214),
                                       ZDistOfTarget = cms.double(15000.0),
                                       PlugVx = cms.double(0.0),
                                       PlugVz = cms.double(-14000.0),
                                       MultiMuonFileFirstEvent = cms.int32(1),
                                       MinT0 = cms.double(-12.5),
                                       RadiusOfTarget = cms.double(8000.0),
                                       RhoPlug = cms.double(2.5),
                                       TIFOnly_linear = cms.bool(False),
                                       MultiMuonNmin = cms.int32(2),
                                       MaxPhi = cms.double(360.0),
                                       Verbosity = cms.bool(False),
                                       MaxT0 = cms.double(12.5),
                                       ZCentrOfTarget = cms.double(0.0),
                                       ElossScaleFactor = cms.double(1.0),
                                       MinTheta = cms.double(0.0),
                                       MinP_CMS = cms.double(-1.0),
                                       RhoRock = cms.double(2.5),
                                       RhoClay = cms.double(2.3),
                                       TrackerOnly = cms.bool(False),
                                       MaxTheta = cms.double(75.0),
                                       MinEnu = cms.double(10.0),
                                       MaxEnu = cms.double(10000.0),
                                       MultiMuon = cms.bool(False),
                                       MultiMuonFileName = cms.string('CORSIKAmultiMuon.root'),
                                       RhoWall = cms.double(2.5)
                                   )


process.SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
                                                          endcapShiftInZNeg = cms.double(0.0),
                                                          PropagationDirection = cms.string('alongMomentum'),
                                                          useMatVolumes = cms.bool(True),
                                                          useTuningForL2Speed = cms.bool(False),
                                                          useIsYokeFlag = cms.bool(True),
                                                          NoErrorPropagation = cms.bool(False),
                                                          SetVBFPointer = cms.bool(False),
                                                          AssumeNoMaterial = cms.bool(False),
                                                          returnTangentPlane = cms.bool(True),
                                                          useInTeslaFromMagField = cms.bool(False),
                                                          VBFName = cms.string('VolumeBasedMagneticField'),
                                                          useEndcapShiftsInZ = cms.bool(False),
                                                          sendLogWarning = cms.bool(False),
                                                          ComponentName = cms.string('SteppingHelixPropagatorAlong'),
                                                          debug = cms.bool(False),
                                                          ApplyRadX0Correction = cms.bool(True),
                                                          useMagVolumes = cms.bool(True),
                                                          endcapShiftInZPos = cms.double(0.0)
                                                      )


process.ProductionFilterSequence = cms.Sequence(process.generator+process.cosmicInTracker)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.endjob_step,process.FEVTDEBUGoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq
