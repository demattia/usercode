import FWCore.ParameterSet.Config as cms

SELECTIONCUT  = " charge = 0"
SELECTIONCUT += " && userFloat('vProb') > 0.001"
SELECTIONCUT += " && abs(rapidity) < 2.4"
SELECTIONCUT += " && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 2."

SELECTIONCUT += " && ((daughter('muon1').pt > 3.5 && abs(daughter('muon1').eta) < 1.6) || (daughter('muon1').pt > 2.5 && 1.6 < abs(daughter('muon1').eta) < 2.4))"
SELECTIONCUT += " && daughter('muon1').isTrackerMuon "
SELECTIONCUT += " && daughter('muon1').muonID('TMOneStationTight')"
SELECTIONCUT += " && daughter('muon1').innerTrack.numberOfValidHits > 11 "
SELECTIONCUT += " && daughter('muon1').innerTrack.hitPattern.pixelLayersWithMeasurement > 0 "
SELECTIONCUT += " && daughter('muon1').innerTrack.normalizedChi2 < 5 "
SELECTIONCUT += " && abs(daughter('muon1').innerTrack.dz) < 25 "
SELECTIONCUT += " && abs(daughter('muon1').dB) < 0.2 "
#SELECTIONCUT += " && daughter('muon1').muonID('TMLastStationAngTight') "

SELECTIONCUT += " && ((daughter('muon2').pt > 3.5 && abs(daughter('muon2').eta) < 1.6) || (daughter('muon2').pt > 2.5 && 1.6 < abs(daughter('muon2').eta) < 2.4))"
SELECTIONCUT += " && daughter('muon2').isTrackerMuon "
SELECTIONCUT += " && daughter('muon2').muonID('TMOneStationTight')"
SELECTIONCUT += " && daughter('muon2').innerTrack.numberOfValidHits > 11 "
SELECTIONCUT += " && daughter('muon2').innerTrack.hitPattern.pixelLayersWithMeasurement > 0 "
SELECTIONCUT += " && daughter('muon2').innerTrack.normalizedChi2 < 5 "
SELECTIONCUT += " && abs(daughter('muon2').innerTrack.dz) < 25 "
SELECTIONCUT += " && abs(daughter('muon2').dB) < 0.2 "
#SELECTIONCUT += " && daughter('muon2').muonID('TMLastStationAngTight') "

SELECTIONCUT += " && (!daughter('muon1').triggerObjectMatchesByFilter('SelectionTrigger').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu0QuarkoniumL3PreFiltered').empty())"
SELECTIONCUT += " && (!daughter('muon2').triggerObjectMatchesByFilter('SelectionTrigger').empty() || !daughter('muon2').triggerObjectMatchesByFilter('hltDoubleMu0QuarkoniumL3PreFiltered').empty())"
#SELECTIONCUT += " && (!daughter('muon1').triggerObjectMatchesByFilter('SelectionTrigger').empty())"
#SELECTIONCUT += " && (!daughter('muon2').triggerObjectMatchesByFilter('SelectionTrigger').empty())"
#SELECTIONCUT += " && deltaR(daughter('muon1').eta, daughter('muon2').eta,daughter('muon1').phi, daughter('muon2').phi)>0.6"
GLOBALCUT = "daughter('muon1').isGlobalMuon && daughter('muon2').isGlobalMuon"
 
IPTree = cms.EDFilter('ProbeTreeProducer',
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    cut = cms.string(""),#SELECTIONCUT + GLOBALCUT),
    sortDescendingBy = cms.string("userFloat('vProb')"),
    maxProbes = cms.int32(1),
    variables = cms.PSet(
        upsPt  = cms.string("pt"),
        upsRapidity   = cms.string("rapidity"),
        invariantMass = cms.string("mass"),
        PDGid = cms.string("userInt('momPDGId')"),

        muPlusPt = cms.string("daughter('muon1').pt"),
        muPlusEta = cms.string("daughter('muon1').eta"),
        muPlusPhi = cms.string("daughter('muon1').phi"),
        muPlusCharge = cms.string("daughter('muon1').charge"),
        muPlusInPosX = cms.string("daughter('muon1').outerTrack.innerPosition.X"),
        muPlusInPosY = cms.string("daughter('muon1').outerTrack.innerPosition.Y"),
        muPlusInPosZ = cms.string("daughter('muon1').outerTrack.innerPosition.Z"),

        muMinusPt = cms.string("daughter('muon2').pt"),
        muMinusEta = cms.string("daughter('muon2').eta"),
        muMinusPhi = cms.string("daughter('muon2').phi"),
        muMinusCharge = cms.string("daughter('muon2').charge"),
        muMinusInPosX = cms.string("daughter('muon2').outerTrack.innerPosition.X"),
        muMinusInPosY = cms.string("daughter('muon2').outerTrack.innerPosition.Y"),
        muMinusInPosZ = cms.string("daughter('muon2').outerTrack.innerPosition.Z"),
    ),
    flags = cms.PSet(),
    ignoreExceptions = cms.bool(True),
    addRunLumiInfo = cms.bool(True),
    filter = cms.bool(True),
)
 
yieldTree = cms.EDFilter('ProbeTreeProducer',
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    cut = cms.string(SELECTIONCUT),
    sortDescendingBy = cms.string("userFloat('vProb')"),
    maxProbes = cms.int32(1),
    variables = cms.PSet(
        upsPt  = cms.string("pt"),
        upsRapidity   = cms.string("rapidity"),
        invariantMass = cms.string("mass"),
        PDGid = cms.string("userInt('momPDGId')"),

        muPlusPt = cms.string("daughter('muon1').pt"),
        muPlusEta = cms.string("daughter('muon1').eta"),
        muPlusPhi = cms.string("daughter('muon1').phi"),
        muPlusCharge = cms.string("daughter('muon1').charge"),
        muPlusInPosX = cms.string("daughter('muon1').outerTrack.innerPosition.X"),
        muPlusInPosY = cms.string("daughter('muon1').outerTrack.innerPosition.Y"),
        muPlusInPosZ = cms.string("daughter('muon1').outerTrack.innerPosition.Z"),

        muMinusPt = cms.string("daughter('muon2').pt"),
        muMinusEta = cms.string("daughter('muon2').eta"),
        muMinusPhi = cms.string("daughter('muon2').phi"),
        muMinusCharge = cms.string("daughter('muon2').charge"),
        muMinusInPosX = cms.string("daughter('muon2').outerTrack.innerPosition.X"),
        muMinusInPosY = cms.string("daughter('muon2').outerTrack.innerPosition.Y"),
        muMinusInPosZ = cms.string("daughter('muon2').outerTrack.innerPosition.Z"),
    ),
    flags = cms.PSet(),
    ignoreExceptions = cms.bool(True),
    addRunLumiInfo = cms.bool(True),
    filter = cms.bool(True),
)

detailedDimuonTree = cms.EDFilter('ProbeTreeProducer',
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        y   = cms.string("rapidity"),
        mass = cms.string("mass"),
        vProb = cms.string("userFloat('vProb')"),
        charge = cms.string("charge"),

        mu1_pt = cms.string("daughter('muon1').pt"),
        mu1_eta = cms.string("daughter('muon1').eta"),
        mu1_phi = cms.string("daughter('muon1').phi"),
        mu1_charge = cms.string("daughter('muon1').charge"),
        mu1_nPixHits = cms.string("daughter('muon1').innerTrack.hitPattern.pixelLayersWithMeasurement"),
        mu1_nTrHits = cms.string("daughter('muon1').innerTrack.numberOfValidHits"),
        mu1_nChi2 = cms.string("daughter('muon1').innerTrack.normalizedChi2"),
        mu1_nMuSegs = cms.string("daughter('muon1').numberOfMatches('SegmentAndTrackArbitration')"),
        mu1_nMuSegsCln = cms.string("daughter('muon1').numberOfMatches('SegmentAndTrackArbitrationCleaned')"),
        mu1_dxy = cms.string("daughter('muon1').innerTrack.dxy"),
        mu1_dB = cms.string("daughter('muon1').dB"),
        mu1_dz = cms.string("daughter('muon1').innerTrack.dz"),
        mu1_R03sumPt = cms.string("daughter('muon1').isolationR03.sumPt"),
        mu1_R03emEt = cms.string("daughter('muon1').isolationR03.emEt"),
        mu1_R03hadEt = cms.string("daughter('muon1').isolationR03.hadEt"),
        mu1_R05sumPt = cms.string("daughter('muon1').isolationR05.sumPt"),
        mu1_R05emEt = cms.string("daughter('muon1').isolationR05.emEt"),
        mu1_R05hadEt = cms.string("daughter('muon1').isolationR05.hadEt"),

        mu2_pt = cms.string("daughter('muon2').pt"),
        mu2_eta = cms.string("daughter('muon2').eta"),
        mu2_phi = cms.string("daughter('muon2').phi"),
        mu2_charge = cms.string("daughter('muon2').charge"),
        mu2_nPixHits = cms.string("daughter('muon2').innerTrack.hitPattern.pixelLayersWithMeasurement"),
        mu2_nTrHits = cms.string("daughter('muon2').innerTrack.numberOfValidHits"),
        mu2_nChi2 = cms.string("daughter('muon2').innerTrack.normalizedChi2"),
        mu2_nMuSegs = cms.string("daughter('muon2').numberOfMatches('SegmentAndTrackArbitration')"),
        mu2_nMuSegsCln = cms.string("daughter('muon2').numberOfMatches('SegmentAndTrackArbitrationCleaned')"),
        mu2_dxy = cms.string("daughter('muon2').innerTrack.dxy"),
        mu2_dB = cms.string("daughter('muon2').dB"),
        mu2_dz = cms.string("daughter('muon2').innerTrack.dz"),
        mu2_R03sumPt = cms.string("daughter('muon2').isolationR03.sumPt"),
        mu2_R03emEt = cms.string("daughter('muon2').isolationR03.emEt"),
        mu2_R03hadEt = cms.string("daughter('muon2').isolationR03.hadEt"),
        mu2_R05sumPt = cms.string("daughter('muon2').isolationR05.sumPt"),
        mu2_R05emEt = cms.string("daughter('muon2').isolationR05.emEt"),
        mu2_R05hadEt = cms.string("daughter('muon2').isolationR05.hadEt"),
    ),
    flags = cms.PSet(
        mu1_Calo = cms.string("daughter('muon1').isCaloMuon"),
        mu1_GM = cms.string("daughter('muon1').isGlobalMuon"),
        mu1_GMPT = cms.string("daughter('muon1').isGlobalMuon && daughter('muon1').muonID('GlobalMuonPromptTight')"),
        mu1_TM = cms.string("daughter('muon1').isTrackerMuon"),
        mu1_TMOST = cms.string("daughter('muon1').isTrackerMuon && daughter('muon1').muonID('TMOneStationTight')"),
        mu1_TMOSAT = cms.string("daughter('muon1').isTrackerMuon && daughter('muon1').muonID('TMOneStationAngTight')"),
        mu1_TMLSAT = cms.string("daughter('muon1').isTrackerMuon && daughter('muon1').muonID('TMLastStationAngTight')"),
        mu1_HLT_Mu3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        mu1_HLT_Mu5 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltSingleMu5L3Filtered5').empty()"),
        mu1_HLT_L1DoubleMuOpen = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
        mu1_HLT_L1DoubleMuOpen_Tight = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty()"),
        mu1_HLT_L2DoubleMu0 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty()"),
        mu1_HLT_DoubleMu0 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty()"),
        mu1_HLT_DoubleMu0_Quarkonium = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu0QuarkoniumL3PreFiltered').empty()"),
        mu1_HLT_DoubleMu3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered').empty()"),

        mu2_Calo = cms.string("daughter('muon2').isCaloMuon"),
        mu2_GM = cms.string("daughter('muon2').isGlobalMuon"),
        mu2_GMPT = cms.string("daughter('muon2').isGlobalMuon && daughter('muon2').muonID('GlobalMuonPromptTight')"),
        mu2_TM = cms.string("daughter('muon2').isTrackerMuon"),
        mu2_TMOST = cms.string("daughter('muon2').isTrackerMuon && daughter('muon1').muonID('TMOneStationTight')"),
        mu2_TMOSAT = cms.string("daughter('muon2').isTrackerMuon && daughter('muon1').muonID('TMOneStationAngTight')"),
        mu2_TMLSAT = cms.string("daughter('muon2').isTrackerMuon && daughter('muon2').muonID('TMLastStationAngTight')"),
        mu2_HLT_Mu3 = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        mu2_HLT_Mu5 = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltSingleMu5L3Filtered5').empty()"),
        mu2_HLT_L1DoubleMuOpen = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
        mu2_HLT_L1DoubleMuOpen_Tight = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty()"),
        mu2_HLT_L2DoubleMu0 = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty()"),
        mu2_HLT_DoubleMu0 = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty()"),
        mu2_HLT_DoubleMu0_Quarkonium = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltDoubleMu0QuarkoniumL3PreFiltered').empty()"),
        mu2_HLT_DoubleMu3 = cms.string("!daughter('muon2').triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered').empty()"),
    ),
    ignoreExceptions = cms.bool(True),
    addRunLumiInfo = cms.bool(True),
    filter = cms.bool(True),
)

tagAndProbe = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables 
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        abseta = cms.string("abs(eta)"),
        L1Q  = cms.string("userInt('muonL1Info:quality')"),
        L1dR = cms.string("userFloat('muonL1Info:deltaR')"),
#        L1dPhi = cms.string("hasUserCand('muonL1Info') && userCand('muonL1Info').pt"),
#        L1dPhi = cms.string("? userCand('muonL1Info').isNonnull ? userCand('muonL1Info').phi : 9"),
        L1bx = cms.string("? userCand('muonL1Info').isNonnull ? userCand('muonL1Info').bx : 3"),
        L1pt = cms.string("? !triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() ? triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).pt : -9 "),
        L1phi = cms.string("? !triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() ? triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).phi : -9 "),
        L1eta = cms.string("? !triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() ? triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).eta : -9 "),
    ),
    flags = cms.PSet(
        # MC truth
        MC = cms.string("genParticleRef(0).isNonnull"),
        MCjpsi = cms.string("genParticleRef(0).isNonnull && genParticleRef(0).motherRef().isNonnull && genParticleRef(0).motherRef().pdgId==443"),
        # Track quality
        TQ2  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 && track.normalizedChi2 < 5"),
        TQ  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 0 && track.normalizedChi2 < 5"),
        # MuonID
        Calo   = cms.string("isCaloMuon"),
        Glb    = cms.string("isGlobalMuon"),
        GlbPT  = cms.string("muonID('GlobalMuonPromptTight')"),
        TM = cms.string("isTrackerMuon"),
        TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
        TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
        TMOSAT = cms.string("muonID('TMOneStationAngTight')"),
        TMOST = cms.string("muonID('TMOneStationTight')"),
        TM2 = cms.string("isTrackerMuon && numberOfMatches('SegmentAndTrackArbitration')>1"),
        # Trigger
        L1SingleMuAll =        cms.string("!triggerObjectMatchesByCollection('hltL1extraParticles').empty()"),
        L1SingleMuOpen =       cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMuOpenL1SingleMu0L1SingleMu3').empty() && userFloat('muonL1Info:deltaR')<0.3"),
        L1SingleMu0 =          cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu0').empty() && userFloat('muonL1Info:deltaR')<0.3"),
        L1SingleMu3 =          cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu3').empty() && userFloat('muonL1Info:deltaR')<0.3"),
        L1SingleMu5 =          cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu5').empty() && userFloat('muonL1Info:deltaR')<0.3"),
        L1SingleMu7 =          cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu7').empty() && userFloat('muonL1Info:deltaR')<0.3"),
        L1SingleMu20 =         cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu20').empty() && userFloat('muonL1Info:deltaR')<0.3"),
        L1MuOpen =             cms.string("!triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty()"),
        L2Mu0 =                cms.string("!triggerObjectMatchesByFilter('hltL2Mu0L2Filtered0').empty()"),
        L3Mu0 =                cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && pt > 0 "),#Filter('hltSingleMu0L3Filtered0').empty()"),
        Mu3 =                  cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        Mu5 =                  cms.string("!triggerObjectMatchesByFilter('hltSingleMu5L3Filtered5').empty()"),
        L1DoubleMuOpen =       cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
        L2DoubleMu0 =          cms.string("!triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty()"),
        DoubleMu0 =            cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty()"),
        DoubleMu0Quarkonium =  cms.string("!triggerObjectMatchesByFilter('hltDoubleMu0QuarkoniumL3PreFiltered').empty()"),
        Mu0 =                  cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty() || !triggerObjectMatchesByFilter('hltDoubleMu0QuarkoniumL3PreFiltered').empty()"),
        DoubleMu3 =            cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered').empty()"),
        Mu0_Track0_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
        Mu0_Track0_Jpsi_TK =   cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                          " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
        Mu3_Track0_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
        Mu3_Track0_Jpsi_TK =   cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                          " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
        Mu3_Track3_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track3JpsiTrackMassFiltered')"),
        Mu3_Track3_Jpsi_TK =   cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                          " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3Track3JpsiTrackMassFiltered')"),
        Mu3_Track5_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track5JpsiTrackMassFiltered')"),
        Mu3_Track5_Jpsi_TK =   cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                          " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3Track5JpsiTrackMassFiltered')"),
        Mu5_Track0_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
        Mu5_Track0_Jpsi_TK =   cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                          " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
        Mu0_TkMu0_Jpsi_MU =    cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TkMuJpsiTkMuMassFiltered')"),
        Mu0_TkMu0_Jpsi_TM =    cms.string("!triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').empty() && "+
                                          " triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasFilterLabel('hltMu0TkMuJpsiTkMuMassFiltered')"),
        # Acceptance definition
        Acc_JPsi = cms.string("(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)"),
    ),
    tagVariables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        abseta = cms.string("abs(eta)"),
        L1Q  = cms.string("userInt('muonL1Info:quality')"),
        L1dR = cms.string("userFloat('muonL1Info:deltaR')"),
        L1pt = cms.string("? !triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() ? triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).pt : -9 "),
        L1phi = cms.string("? !triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() ? triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).phi : -9 "),
        L1eta = cms.string("? !triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() ? triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).eta : -9 "),
    ),
    tagFlags = cms.PSet(
        L1MuOpen =             cms.string("!triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty()"),
        L1DoubleMuOpen =       cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
        L2Mu0 =                cms.string("!triggerObjectMatchesByFilter('hltL2Mu0L2Filtered0').empty()"),
        L3Mu0 = 	       cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && pt > 0"),#Filter('hltSingleMu0L3Filtered0').empty()"),
        Mu3 =                  cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        Mu5 =                  cms.string("!triggerObjectMatchesByFilter('hltSingleMu5L3Filtered5').empty()"),
        L2DoubleMu0 =          cms.string("!triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty()"),
        Mu0_Track0_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
        Mu3_Track0_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
        Mu5_Track0_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
        Mu3_Track3_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track3JpsiTrackMassFiltered')"),
        Mu3_Track5_Jpsi_MU =   cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track5JpsiTrackMassFiltered')"),
        MuX_TrackY_Jpsi_MU =  cms.string("(!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                                 "triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered'))"+
                                                 "||(!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                                 "triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered'))"+
                                                 "||(!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                                 "triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered'))"+
                                                 "||(!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                                 "triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track3JpsiTrackMassFiltered'))"+
                                                 "||(!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                                 "triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track5JpsiTrackMassFiltered'))"),
        Mu5L2Mu0 =             cms.string("!triggerObjectMatchesByFilter('hltMu5L2Mu0L3Filtered5').empty()"),
    ),
    pairVariables = cms.PSet(
        dphiVtxTimesQ = cms.InputTag("nearbyMuonsInfo", "dphiVtxTimesQ"),
        drM1          = cms.InputTag("nearbyMuonsInfo", "drM1"),
        drM2          = cms.InputTag("nearbyMuonsInfo", "drM2"),
        dphiM2        = cms.InputTag("nearbyMuonsInfo", "dphiM2"),
        distM2        = cms.InputTag("nearbyMuonsInfo", "distM2"),
        drStaIn       = cms.InputTag("nearbyMuonsInfo", "drStaIn"),
        dphiStaIn     = cms.InputTag("nearbyMuonsInfo", "dphiStaIn"),
    ),
    pairFlags = cms.PSet(
        badCowboy = cms.string("(daughter(0).charge * deltaPhi(daughter(0).phi, daughter(1).phi))>0 && abs(daughter(0).eta-daughter(1).eta)<0.2 && 1.4<abs(daughter(1).eta) && abs(daughter(1).eta)<1.6"),
#        L1arbitrated = cms.string("daughter(0).triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() || "+
#                                  "daughter(1).triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty() || "+
#                                  "deltaR( daughter(0).triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).eta,"+
#                                  "daughter(0).triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).phi,"+
#                                  "daughter(1).triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).eta,"+
#                                  "daughter(1).triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').at(0).phi)>0.001"),
    ),
    isMC = cms.bool(True),
    tagMatches = cms.InputTag("tagMuonsMCMatch"),
    probeMatches = cms.InputTag("probeMuonsMCMatch"),
    makeMCUnbiasTree = cms.bool(True),
    motherPdgId = cms.vint32(22, 23, 443, 100443, 553, 100553, 200553),
    checkMotherInUnbiasEff = cms.bool(True),
    addRunLumiInfo = cms.bool(True),
    allProbes = cms.InputTag("probeMuons"),
)

def selection(process, GlobalTag="GR_R_38X_V8::All", MC=False, SelectionTrigger="hltL1DoubleMuOpenTightL1Filtered"):
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("Configuration.StandardSequences.Geometry_cff")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.GlobalTag.globaltag = cms.string(GlobalTag)
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.load("MuonAnalysis.TagAndProbe.nearbyMuonsInfo_cfi")
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        # SkipEvent = cms.untracked.vstring('ProductNotFound'),
        fileMode = cms.untracked.string('NOMERGE')
    )
    process.MessageLogger.cerr.FwkReport.reportEvery = 100

    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring()
    )

    process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
        vertexCollection = cms.InputTag('offlinePrimaryVertices'),
        minimumNDOF = cms.uint32(4),
        maxAbsZ = cms.double(25),	
        maxd0 = cms.double(2)	
    )

    process.detailedDimuonTree = detailedDimuonTree

  #  process.myEventIdFilter = cms.EDFilter("EventIdFilter")
  #  process.myEventIdFilterPath = cms.Path(
  #      process.myEventIdFilter
  #  )

    process.detailedDimuonTreePath = cms.Path(
        process.detailedDimuonTree
    )

    CUT = SELECTIONCUT.replace("SelectionTrigger", SelectionTrigger)
    process.yieldPsiTree = yieldTree.clone(
        cut = cms.string(CUT + " && 2.5 < mass < 4.5"),
    )
    process.yieldPsiTreePath = cms.Path(
        process.primaryVertexFilter +
        process.yieldPsiTree
    )

    process.yieldUpsilonTree = yieldTree.clone(
        cut = cms.string(CUT + " && 6 < mass < 15"),
    )
    process.yieldUpsilonTreePath = cms.Path(
        process.primaryVertexFilter +
        process.yieldUpsilonTree
    )

    process.yieldZTree = yieldTree.clone(
        cut = cms.string(CUT + " && 60 < mass < 120"),
    )
    process.yieldZTreePath = cms.Path(
        process.primaryVertexFilter +
        process.yieldZTree
    )

    process.yieldMCTree = yieldTree.clone(
        cut = cms.string( "userInt('momPDGId')!=0"
#((daughter('muon1').pt > 3.5 && abs(daughter('muon1').eta) < 1.6) || (daughter('muon1').pt > 2.5 && 1.6 < abs(daughter('muon1').eta) < 2.4))"+
#                          "&& ((daughter('muon2').pt > 3.5 && abs(daughter('muon2').eta) < 1.6) || (daughter('muon2').pt > 2.5 && 1.6 < abs(daughter('muon2').eta) < 2.4))"+
#                          "&& userInt('momPDGId')!=0"
        ),
    )
    process.yieldMCTreePath = cms.Path(
        process.yieldMCTree
    )


    # Tag and Probe

    tagAndProbe.isMC = MC
    process.tagAndProbe = tagAndProbe


    process.tagAndProbePath = cms.Path(
        process.primaryVertexFilter *
        process.nearbyMuonsInfo *
        process.tagAndProbe
    )

    process.probe = cms.EDFilter("ProbeTreeProducer",
        src = cms.InputTag("probeMuons"),
        cut = cms.string("genParticleRef(0).isNonnull"),
        variables = tagAndProbe.variables,
        flags = tagAndProbe.flags,
        addRunLumiInfo = cms.bool(True),
        filter = cms.bool(True),
    )

    process.probePath = cms.Path(
        process.probe
    )

    process.TFileService = cms.Service("TFileService", fileName = cms.string("selection.root"))

