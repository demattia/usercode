import FWCore.ParameterSet.Config as cms

analyzeLeptons = cms.EDAnalyzer("LeptonAnalysis",

  # input collections
  electronSrc  = cms.InputTag("selectedPatElectrons"),
  muonSrc      = cms.InputTag("selectedPatMuons"),
  tauSrc       = cms.InputTag("selectedPatTaus"),
  trackSrc     = cms.InputTag("pseudoLeptonProducer"),
  generatorSrc = cms.InputTag("genParticles"),
  trigger      = cms.InputTag("patTrigger"),
  triggerEvent = cms.InputTag("patTriggerEvent"),
  barrelSuperClusters = cms.InputTag("correctedHybridSuperClusters"),
  endcapSuperClusters = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),

  # general analysis setup

  # define which PDG ID is taken as signal
  # 23 = Z
  # 24 = W+
  # 36 = A0
  # 37 = H+
  # 100022 = CalcHEP heavy photon
  signalPDGId = cms.int32(36),

  # cut study mode: true = keep all candidates throughout entire selection to fill pass/fail plots
  #                 false= reject failed candidates immediately
  cutStudyMode = cms.bool(True),

  # define which channel to look at here (with bogus default to force user input)
  leptonPDGId = cms.int32(0),

  # list of triggers we want to look at (here: initialize to a dummy trigger)
  hltPaths = cms.vstring("HLTriggerFinalPath"),

  # number of leptons for which we require a trigger object match (0,1,2)
  numTrigMatches = cms.int32(2),

  # run range for histograms (here: 2010-2011 data)
  minRunNumber = cms.int32(132440),
  maxRunNumber = cms.int32(180300),

  # general selection cuts non individual leptons
  leptonPtCut = cms.double(10.),
  leptonEtaCut = cms.double(2.4),
  leptonIsolationCut = cms.double(4.0),
  minD0Significance = cms.double(2.0),

  # cosmic veto: reject candidates with tracks with cos(angle)<vetoBackToBack
  vetoBackToBack = cms.double(-999),

  # reject candidates where the deltaR between the leptons is too small?
  minDeltaRBetweenLeptons = cms.double(-999),
                                           
  # how many standalone muons are allowed in construction candidates
  # note that this must be zero for the etrack channel, because otherwise the
  # analysis code will match standalone muons to photon trigger objects!
  maxNumStandAloneMuons = cms.int32(0),

  # for etrack channel: require minimum number of offline supercluster matches
  minNumCaloMatches = cms.int32(0),

  # cuts on the dilepton candidate
  leptonChargeCut = cms.bool(True),
  vertexChi2Cut = cms.double(5),
  deltaPhiCut = cms.double(0.2),
  decayLengthCut = cms.double(0.1),
  decayLengthSignificanceCut = cms.double(5.0),
  maxHitsBeforeVertex = cms.int32(1),

  FillHistograms = cms.bool(True),

  UseMCTruth = cms.bool(True),

  # correction weights (sample-specific and thus meant to be overwritten)

  # Summer11 PU_S4, distribution obtained by averaging the number of
  # interactions in each beam crossing to estimate the true mean
  lumiDistrMC = cms.vdouble(
    0.104109,
    0.0703573,
    0.0698445,
    0.0698254,
    0.0697054,
    0.0697907,
    0.0696751,
    0.0694486,
    0.0680332,
    0.0651044,
    0.0598036,
    0.0527395,
    0.0439513,
    0.0352202,
    0.0266714,
    0.019411,
    0.0133974,
    0.00898536,
    0.0057516,
    0.00351493,
    0.00212087,
    0.00122891,
    0.00070592,
    0.000384744,
    0.000219377
  ),

  lumiDistrData = cms.vdouble(
    14.5417,
    34.7743,
    78.9247,
    126.467,
    159.329,
    167.603,
    152.684,
    123.794,
    90.9462,
    61.3973,
    38.505,
    22.628,
    12.5503,
    6.61051,
    3.32403,
    1.60286,
    0.743920,
    0.333477,
    0.144861,
    0.0611127,
    0.0251102,
    0.0100651,
    0.00394398,
    0.00151354,
    0.000896161
    )
)
