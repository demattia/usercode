import FWCore.ParameterSet.Config as cms

process = cms.Process("SkimmingOnia2MuMuPAT")

### standard includes
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

### global tag
process.GlobalTag.globaltag = "START3X_V18::All"

### source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_3_4_2/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V15-v1/0011/D68F77F2-7213-DF11-BAA7-002618943842.root',
       '/store/relval/CMSSW_3_4_2/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V15-v1/0011/7C3D4BFE-7113-DF11-B956-002618943884.root',
       '/store/relval/CMSSW_3_4_2/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V15-v1/0011/300CEE61-7013-DF11-8622-00304867924A.root',
       '/store/relval/CMSSW_3_4_2/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V15-v1/0011/22F71073-B413-DF11-BBBE-0030486792B8.root',
       '/store/relval/CMSSW_3_4_2/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V15-v1/0010/1825DD6D-6A13-DF11-AA3F-001BFCDBD15E.root'
    )
)

### number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff")
process.jpsiRecoSeq = cms.Sequence(process.patMuonSequence * process.onia2MuMuPatGlbTrk)

process.genJPsi = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("pdgId == 443"),
)

def mkh(name,title,expr,bins,min,max) :
    return cms.PSet(
        name = cms.untracked.string(name),
        description  = cms.untracked.string(title),
        plotquantity = cms.untracked.string(expr),
        nbins = cms.untracked.int32(bins),
        min   = cms.untracked.double(min),
        max   = cms.untracked.double(max),
    )


process.plotJPsi = cms.EDAnalyzer("CandViewHistoAnalyzer",
    src = cms.InputTag("genJPsi"),
    histograms = cms.VPSet(
        mkh("pt",  "J/Psi P_{T}",          "pt",       20,  0, 20),    
        mkh("eta", "J/Psi pseudorapidity", "eta",      15, -3,  3),    
        mkh("y",   "J/Psi rapidity",       "rapidity", 15, -3,  3),    
    )
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    highEta = cms.PSet(initialSeed = cms.untracked.uint32(37), engineName = cms.untracked.string('HepJamesRandom')),
    lowEta  = cms.PSet(initialSeed = cms.untracked.uint32(37), engineName = cms.untracked.string('HepJamesRandom')),
    highPt  = cms.PSet(initialSeed = cms.untracked.uint32(37), engineName = cms.untracked.string('HepJamesRandom')),
    lowPt   = cms.PSet(initialSeed = cms.untracked.uint32(37), engineName = cms.untracked.string('HepJamesRandom')),
)

process.highEta = cms.EDFilter("GenSpectrumFilter",
    src = cms.InputTag("genParticles"), 
    pdgIds = cms.vint32(443),
    probability = cms.string("1-0.05*rapidity*rapidity")
)
process.lowEta = process.highEta.clone(probability = "0.6+0.05*rapidity*rapidity")
process.lowPt  = process.highEta.clone(probability = "exp(-0.02*pt)")
process.highPt = process.highEta.clone(probability = "1-0.5*exp(-0.1*pt)")

process.core     = cms.Path(process.genJPsi * process.plotJPsi)

process.plotRecoJPsi = process.plotJPsi.clone(src = cms.InputTag("onia2MuMuPatGlbTrk"))
process.coreReco = cms.Path(process.jpsiRecoSeq * process.plotRecoJPsi)

for X in ['lowEta', 'lowPt', 'highEta', 'highPt']:
    setattr(process,X+'Plots', process.plotJPsi.clone())
    setattr(process,X+'Path',  cms.Path(getattr(process,X) + process.genJPsi * getattr(process,X+'Plots')))
    setattr(process,X+'RecoPlots', process.plotRecoJPsi.clone())
    setattr(process,X+'RecoPath',  cms.Path(getattr(process,X) + process.jpsiRecoSeq * getattr(process,X+'RecoPlots')))

process.TFileService = cms.Service("TFileService", fileName = cms.string("genRecoSpectrum.root"))
