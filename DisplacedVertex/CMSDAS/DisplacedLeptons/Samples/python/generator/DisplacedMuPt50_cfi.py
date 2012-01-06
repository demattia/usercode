import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

# modeled after Configuration/Generator/python/SingleMuPt100_cfi.py
#
# this file creates single 50 GeV p_t muons,
# with an impact parameter with respect to 0,0,0 by up to 100cm


# particle gun
generator = cms.EDProducer("FlatRandomD0GunProducer",
    PGunParameters = cms.PSet(
        MinPt = cms.double(49.99), # in GeV
        MaxPt = cms.double(50.01),
        MinD0 = cms.double(0.0),   # in mm
        MaxD0 = cms.double(1000.0),
        PartID = cms.vint32(-13),
        MinEta = cms.double(-0.1),
        MaxEta = cms.double(0.1),
        MinPhi = cms.double(-3.14159265359), # in radians
        MaxPhi = cms.double(3.14159265359)

    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('single mu pt 50'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
