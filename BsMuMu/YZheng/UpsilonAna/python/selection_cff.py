import FWCore.ParameterSet.Config as cms

# #dimuon
# 
# SELECTIONCUT  = " charge = 0"
# SELECTIONCUT += " && userFloat('vProb') > 0.001"
# SELECTIONCUT += " && abs(rapidity) < 2."
# SELECTIONCUT += " && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 2."
# 
# #each muon
# SELECTIONCUT += " && ((daughter('muon1').pt > 3.5 && abs(daughter('muon1').eta) < 1.6) || (daughter('muon1').pt > 2.5 && 1.6 < abs(daughter('muon1').eta) < 2.4))"
# SELECTIONCUT += " && daughter('muon1').isTrackerMuon "
# SELECTIONCUT += " && daughter('muon1').innerTrack.numberOfValidHits > 11 "
# SELECTIONCUT += " && daughter('muon1').innerTrack.hitPattern.pixelLayersWithMeasurement > 0 "
# SELECTIONCUT += " && daughter('muon1').innerTrack.normalizedChi2 < 5 "
# SELECTIONCUT += " && abs(daughter('muon1').innerTrack.dz) < 25 "
# SELECTIONCUT += " && abs(daughter('muon1').dB) < 0.2 "
# #SELECTIONCUT += " && daughter('muon1').muonID('TMLastStationAngTight') "
# 
# SELECTIONCUT += " && ((daughter('muon2').pt > 3.5 && abs(daughter('muon2').eta) < 1.6) || (daughter('muon2').pt > 2.5 && 1.6 < abs(daughter('muon2').eta) < 2.4))"
# SELECTIONCUT += " && daughter('muon2').isTrackerMuon "
# SELECTIONCUT += " && daughter('muon2').innerTrack.numberOfValidHits > 11 "
# SELECTIONCUT += " && daughter('muon2').innerTrack.hitPattern.pixelLayersWithMeasurement > 0 "
# SELECTIONCUT += " && daughter('muon2').innerTrack.normalizedChi2 < 5 "
# SELECTIONCUT += " && abs(daughter('muon2').innerTrack.dz) < 25 "
# SELECTIONCUT += " && abs(daughter('muon2').dB) < 0.2 "
# #SELECTIONCUT += " && daughter('muon2').muonID('TMLastStationAngTight') "
# 
# # SELECTIONCUT += " && !daughter('muon1').triggerObjectMatchesByFilter('SelectionTrigger').empty() && daughter('muon1').userFloat('muonL1Info:deltaR')<0.3"
# # SELECTIONCUT += " && !daughter('muon2').triggerObjectMatchesByFilter('SelectionTrigger').empty() && daughter('muon2').userFloat('muonL1Info:deltaR')<0.3"
# SELECTIONCUT += " && (!daughter('muon1').triggerObjectMatchesByFilter('SelectionTrigger').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu3BsL3Filtered').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2BarrelBsL3Filtered').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2Dimuon6BsL3Filtered').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs6').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs4').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs345').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs3p545').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs47').empty())"
# SELECTIONCUT += " && (!daughter('muon2').triggerObjectMatchesByFilter('SelectionTrigger').empty() || !daughter('muon2').triggerObjectMatchesByFilter('hltDoubleMu3BsL3Filtered').empty() || !daughter('muon2').triggerObjectMatchesByFilter('hltDoubleMu2BarrelBsL3Filtered').empty() || !daughter('muon2').triggerObjectMatchesByFilter('hltDoubleMu2Dimuon6BsL3Filtered').empty() || !daughter('muon2').triggerObjectMatchesByFilter('hltVertexmumuFilterBs6').empty() || !daughter('muon2').triggerObjectMatchesByFilter('hltVertexmumuFilterBs4').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs345').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs3p545').empty() || !daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs47').empty())"

# yieldTree = cms.EDFilter('ProbeTreeProducer',
#     src = cms.InputTag("onia2MuMuPatTrkTrk"),
#     cut = cms.string(SELECTIONCUT),
#     sortDescendingBy = cms.string("userFloat('vProb')"),
#     maxProbes = cms.int32(1),
#     variables = cms.PSet(
#         upsPt  = cms.string("pt"),
#         upsRapidity   = cms.string("rapidity"),
#         invariantMass = cms.string("mass"),
#         PDGid = cms.string("userInt('momPDGId')"),
# 
#         muPlusPt = cms.string("daughter('muon1').pt"),
#         muPlusEta = cms.string("daughter('muon1').eta"),
#         muPlusPhi = cms.string("daughter('muon1').phi"),
#         muPlusCharge = cms.string("daughter('muon1').charge"),
# 
#         muMinusPt = cms.string("daughter('muon2').pt"),
#         muMinusEta = cms.string("daughter('muon2').eta"),
#         muMinusPhi = cms.string("daughter('muon2').phi"),
#         muMinusCharge = cms.string("daughter('muon2').charge"),
#     ),
#     flags = cms.PSet(),
#     ignoreExceptions = cms.bool(True),
#     addRunLumiInfo = cms.bool(True),
#     filter = cms.bool(True),
# )

detailedDimuonTree = cms.EDFilter('ProbeTreeProducer',
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    variables = cms.PSet(
        pt = cms.string("pt"),
        y  = cms.string("rapidity"),
        mass = cms.string("mass"),
        vProb = cms.string("userFloat('vProb')"),
        pvw8 = cms.string("userFloat('pvw8')"),
        charge = cms.string("charge"),
        isolation = cms.string("userFloat('Isolation')"),
        minDca = cms.string("userFloat('minDca')"),
        ntrk = cms.string("userInt('Ntrk')"),
        countTksOfPV = cms.string("userInt('countTksOfPV')"),
        vertexWeight = cms.string("userFloat('vertexWeight')"),
        sumPTPV = cms.string("userFloat('sumPTPV')"),
        dcaxy = cms.string("userFloat('DCAXY')"),
        dca = cms.string("userFloat('DCA')"),
        ctauPV = cms.string("userFloat('ppdlPV')"),
        ctauErrPV = cms.string("userFloat('ppdlErrPV')"),
        cosAlphaXY = cms.string("userFloat('cosAlphaXY')"),
        cosAlpha3D = cms.string("userFloat('cosAlpha3D')"),
        delta3d = cms.string("userFloat('delta3d')"),
        delta3dErr = cms.string("userFloat('delta3dErr')"),
        PDGid = cms.string("userInt('momPDGId')"),
        NChi2 = cms.string("userFloat('vNChi2')"),
        highPurity1 = cms.string("userInt('highPurity1')"),
        highPurity2 = cms.string("userInt('highPurity2')"),
        l3d = cms.string("userFloat('l3d')"),
        l3dsig = cms.string("userFloat('l3dsig')"),

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
        mu1_HLT_DoubleMu2BsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2BsL3Filtered').empty()"),
        mu1_HLT_DoubleMu3BsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu3BsL3Filtered').empty()"),
        mu1_HLT_DoubleMu2BarrelBsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2BarrelBsL3Filtered').empty()"),
        mu1_HLT_DoubleMu2Dimuon6BsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2Dimuon6BsL3Filtered').empty()"),
        mu1_HLT_VertexmumuFilterBs6 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs6').empty()"),
        mu1_HLT_VertexmumuFilterBs4 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs4').empty()"),
        mu1_HLT_VertexmumuFilterBs345 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs345').empty()"),
        mu1_HLT_VertexmumuFilterBs3p545 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs3p545').empty()"),
        mu1_HLT_VertexmumuFilterBs47 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs47').empty()"),
        
        mu2_Calo = cms.string("daughter('muon2').isCaloMuon"),
        mu2_GM = cms.string("daughter('muon2').isGlobalMuon"),
        mu2_GMPT = cms.string("daughter('muon2').isGlobalMuon && daughter('muon2').muonID('GlobalMuonPromptTight')"),
        mu2_TM = cms.string("daughter('muon2').isTrackerMuon"),
        mu2_TMOST = cms.string("daughter('muon2').isTrackerMuon && daughter('muon1').muonID('TMOneStationTight')"),
        mu2_TMOSAT = cms.string("daughter('muon2').isTrackerMuon && daughter('muon1').muonID('TMOneStationAngTight')"),
        mu2_TMLSAT = cms.string("daughter('muon2').isTrackerMuon && daughter('muon2').muonID('TMLastStationAngTight')"),
        mu2_HLT_DoubleMu2BsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2BsL3Filtered').empty()"),
        mu2_HLT_DoubleMu3BsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu3BsL3Filtered').empty()"),
        mu2_HLT_DoubleMu2BarrelBsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2BarrelBsL3Filtered').empty()"),
        mu2_HLT_DoubleMu2Dimuon6BsL3 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltDoubleMu2Dimuon6BsL3Filtered').empty()"),
        mu2_HLT_VertexmumuFilterBs6 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs6').empty()"),
        mu2_HLT_VertexmumuFilterBs4 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs4').empty()"),
        mu2_HLT_VertexmumuFilterBs345 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs345').empty()"),
        mu2_HLT_VertexmumuFilterBs3p545 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs3p545').empty()"),
        mu2_HLT_VertexmumuFilterBs47 = cms.string("!daughter('muon1').triggerObjectMatchesByFilter('hltVertexmumuFilterBs47').empty()")
    ),
    ignoreExceptions = cms.bool(True),
    addRunLumiInfo = cms.bool(True),
    filter = cms.bool(True),
)

def selection(process, GlobalTag="GR_R_38X_V8::All", MC=False, SelectionTrigger="hltL1DoubleMuOpenTightL1Filtered"):
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("Configuration.StandardSequences.Geometry_cff")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.GlobalTag.globaltag = cms.string(GlobalTag)
    process.load("FWCore.MessageService.MessageLogger_cfi")
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

    process.detailedDimuonTreePath = cms.Path(
        process.detailedDimuonTree
    )

    CUT = SELECTIONCUT.replace("SelectionTrigger", SelectionTrigger)
    process.yieldPsiTree = yieldTree.clone(
        cut = cms.string(CUT + " && 2.5 < mass < 4.5"),
    )

    # process.yieldMCTree = yieldTree.clone(
    #     cut = cms.string( "((daughter('muon1').pt > 3.5 && abs(daughter('muon1').eta) < 1.6) || (daughter('muon1').pt > 2.5 && 1.6 < abs(daughter('muon1').eta) < 2.4))"+
    #                       "&& ((daughter('muon2').pt > 3.5 && abs(daughter('muon2').eta) < 1.6) || (daughter('muon2').pt > 2.5 && 1.6 < abs(daughter('muon2').eta) < 2.4))"+
    #                       "&& userInt('momPDGId')!=0"
    #     ),
    # )
    # process.yieldMCTreePath = cms.Path(
    #     process.yieldMCTree
    # )

    process.TFileService = cms.Service("TFileService", fileName = cms.string("selection_test.root"))

