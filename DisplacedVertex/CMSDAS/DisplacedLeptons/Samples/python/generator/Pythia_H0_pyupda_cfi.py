import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUESettings_cfi import *

#
# This produces an H0, which decays to a pair of long-lived exotic particles
# which then each decay to a pair of leptons
#

generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(10),
    pythiaPylistVerbosity = cms.untracked.int32(3),
    filterEfficiency = cms.untracked.double(1.0),
    comEnergy = cms.double(7000.0),
    UseExternalGenerators = cms.untracked.bool(False),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        pythiaMyParameters = cms.vstring(
          'MSTJ(22)=4    ! Decay ALL unstable particles inside a cylinder (needed for long-lived exotica)',
          'PARJ(73)=2000. ! radius of cylinder',
          'PARJ(74)=4000. ! half length of cylinder',
          'MSEL=0',
          'MSUB(152)=1   ! gg->H0 production',
          'MWID(35)=2    ! Let me set H0 properties'
        ),
        PYUPDAParameters = cms.vstring(
          "PYUPDAFILE = 'DisplacedLeptons/Samples/data/Pythia_H0_pyupda.in'"
        ),
        parameterSets = cms.vstring('pythiaUESettings', 
                                    'pythiaMyParameters',
                                    'PYUPDAParameters')
    )
)
