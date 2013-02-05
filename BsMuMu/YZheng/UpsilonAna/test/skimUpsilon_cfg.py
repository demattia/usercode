import FWCore.ParameterSet.Config as cms

process = cms.Process("Skim")

ptMin  = 2.5; # minimum pT to consider tracks and muons
ptTag  = 3.0; # minimum pT for the tag
etaTag = 2.5; # maximum |eta| for the tag
triggerProcess = 'HLT8E29'  # when running on 3.1.X MC
triggerPath    = 'HLT_Mu3'  # Trigger name
massRange  = (7.0, 12.0)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",  
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer09/Upsilon1S/GEN-SIM-RECO/MC_31X_V3-v1/0011/FEE67276-CC89-DE11-959B-0030487722A4.root',
        '/store/mc/Summer09/Upsilon1S/GEN-SIM-RECO/MC_31X_V3-v1/0011/F250BF03-7289-DE11-8EAE-00304889D476.root'
    )
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

## ==== pruner of GenParticles ====
process.genMuons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ",                     # this is the default
        "++keep abs(pdgId) = 13",        # keep muons and their parents
        "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
    )
)

## ==== Just ignore all the too low pt stuff ====
process.goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string("pt > %f" % (ptMin,) ),
)

process.slimAOD = cms.Sequence(
    process.genMuons +
    process.goodTracks 
)

### ==== Make PAT Muons ====
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import allLayer1Muons
process.patMuonsWithoutTrigger = allLayer1Muons.clone(
    # embed the tracks, so we don't have to carry them around
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    # then switch off some features we don't need
    embedPickyMuon = False,
    embedTpfmsMuon = False, 
    isolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(), # no heavy isodeposits
    addGenMatch = False,      # no mc: T&P doesn't take it from here anyway.
)

### ==== Unpack trigger, and match ====
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.patTrigger.onlyStandAlone = True
process.patTrigger.processName    = triggerProcess
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import muonTriggerMatchHLT1MuonIso
process.muonMatchHLTMuX = muonTriggerMatchHLT1MuonIso.clone(
    src = 'patMuonsWithoutTrigger',
    pathNames = [ triggerPath ]
)

### ==== Embed ====
process.patMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag( "patMuonsWithoutTrigger" ),
    matches = cms.VInputTag( "muonMatchHLTMuX" ),
)

### ==== Select ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > %f" % (ptMin,)), 
)

### ==== Sequence ====
process.patMuonSequence = cms.Sequence( 
    process.patMuonsWithoutTrigger *
    process.patTrigger * process.muonMatchHLTMuX * process.patMuonsWithTrigger *
    process.patMuons  
)

## ==== For bare tracks, make candidates assuming the muon mass hypothesis ====
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi");
process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodTracks"),      
    particleType = cms.string("mu+"),
) 

## ==== Tracker muons ====
process.trackerMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isTrackerMuon"),
)

## ==== Golden muons to be used for tags. use PAT ones, so I can check HLT =====
#PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % (triggerPath,);
process.goldenMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
#    cut = cms.string("isGlobalMuon && pt > %f && abs(eta) < %f && %s" % (ptTag, etaTag, PASS_HLT) ), 
    cut = cms.string("isGlobalMuon && pt > %f && abs(eta) < %f" % (ptTag, etaTag) ),
)

process.otherMuons = cms.Sequence(
    process.goldenMuons +
    process.tkTracks +
    process.trackerMuons
)

process.allMuons = cms.Sequence(
    process.patMuonSequence *
    process.otherMuons
)

process.main = cms.Path(
    process.slimAOD  *
    process.allMuons
)

## ==== Filter on events that have at least one Tag and Probe pair or one TM-TM pair
process.upsilonTagTk  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ tkTracks@-"),
    cut = cms.string("%f < mass < %f" % massRange),
)
process.upsilonTagTM  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ trackerMuons@-"),
    cut = cms.string("%f < mass < %f" % massRange),
)
process.upsilonTMTM  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("trackerMuons@+ trackerMuons@-"),
    cut = cms.string("%f < mass < %f" % massRange),
)
process.upsilonTagTkFilter  = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("upsilonTagTk"),
    minNumber = cms.uint32(1),
)
process.upsilonTagTMFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('upsilonTagTM'),
    minNumber = cms.uint32(1),
)
process.upsilonTMTMFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('upsilonTMTM'),
    minNumber = cms.uint32(1),
)
process.Skim_upsilonTagTk  = cms.Path(process.upsilonTagTk  * process.upsilonTagTkFilter )
process.Skim_upsilonTagTM = cms.Path(process.upsilonTagTM * process.upsilonTagTMFilter)
process.Skim_upsilonTMTM = cms.Path(process.upsilonTMTM * process.upsilonTMTMFilter)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skimUpsilon.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_patMuons_*_Skim",
        "keep *_trackerMuons_*_Skim",
        "keep *_goodTracks_*_Skim",
        "keep *_genMuons_*_Skim",
        "keep recoTrackExtras_standAloneMuons_*_*",       ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons_UpdatedAtVtx_*", ## bare standalone muon tracks, using standalone muon momentum (with BS constraint)
        "keep edmTriggerResults_*_*_Skim",          ## to know which kind of skim channel got us the event   
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring(
        "Skim_upsilonTagTk",
        "Skim_upsilonTagTM",
        "Skim_upsilonTMTM"
    ))
)
process.end = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.main,
    process.Skim_upsilonTagTk,
    process.Skim_upsilonTagTM,
    process.Skim_upsilonTMTM,
    process.end
)

