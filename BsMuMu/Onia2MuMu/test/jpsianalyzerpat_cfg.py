import FWCore.ParameterSet.Config as cms

process = cms.Process("ANAPAT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "GR_R_42_V14::All"
#process.GlobalTag.globaltag = "FT_R_42_V13A::All"
process.GlobalTag.globaltag = "GR_P_V20::All"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring("file:onia2MuMuPAT.root")
)

process.hltMuF = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_DoubleMu0"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)


# filter on good vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
                                           )

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.demo = cms.EDAnalyzer('JPsiAnalyzerPAT',

    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    srcWithCaloMuons = cms.InputTag("onia2MuMuPatGlbCal"),

    writeTree = cms.bool(True),
    treeFileName = cms.string("onia2MuMu_tree.root"),

    writeDataSet = cms.bool(True),                 
    dataSetName = cms.string("dataSet.root"),
    triggersForDataset = cms.vstring("HLT_DoubleMu3_Jpsi_v1"),

    massMin = cms.double(2.6),
    massMax = cms.double(3.5),
    pTBinRanges = cms.vdouble(0.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 40.0),
    etaBinRanges = cms.vdouble(0.0, 1.3, 2.5),
    onlyTheBest = cms.bool(True),		
    applyCuts = cms.bool(True),
    applyExpHitCuts = cms.untracked.bool(False),
    applyDiMuonCuts = cms.untracked.bool(False), #H:                         
    useBeamSpot = cms.bool(False),
    useCaloMuons = cms.untracked.bool(False),
    removeSignalEvents = cms.untracked.bool(False),
    removeTrueMuons = cms.untracked.bool(False),
    storeWrongSign = cms.untracked.bool(False),
    writeOutCandidates = cms.untracked.bool(False),
    massCorrectionMode = cms.int32(3),    # mode 0 no correction,
                                          # mode 1 constant corr,
                                          # mode 2 pt dependent corr,
                                          # mode 3 pt and eta dependent corr
    oniaPDG = cms.int32(443),
    genParticles = cms.InputTag("genMuons"),
    isMC = cms.untracked.bool(False),
    storeAllMCEvents = cms.untracked.bool(True),
    isPromptMC = cms.untracked.bool(True),
                              
    # Configuration for the extrapolation at the muon system 
    propagatorStation1 = cms.PSet(
        useStation2 = cms.bool(False), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # for AOD
        useSimpleGeometry = cms.bool(True), 
    ),
    propagatorStation2 = cms.PSet(
        useStation2 = cms.bool(True), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # for AOD
        useSimpleGeometry = cms.bool(True), 
    ),

    # Configuration of trigger matching                           
    triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
                              
    ####### 2011, 5E32 ########
                              
    HLTBitNames_SingleMu = cms.vstring(),
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_SingleMu = cms.vstring(),

    HLTBitNames_DoubleMu = cms.vstring("HLT_DoubleMu3_Jpsi_v1", # 5E32, v4.2-v5.3
                                       "HLT_DoubleMu3_Jpsi_v2", # 5E32, v6.1-v6.2
                                       "HLT_Dimuon6p5_Jpsi_v1", #5E32, v8.1-v8.3
                                       "HLT_Dimuon6p5_Jpsi_Displaced_v1", #5E32, v8.1-v8.3
                                       "HLT_Dimuon6p5_Barrel_Jpsi_v1", #5E32, v8.1-v8.3
                                       "HLT_DoubleMu3_Quarkonium_v1", #5E32, v4.2-v5.3
                                       "HLT_DoubleMu3_Quarkonium_v2", #5E32, v6.1-v6.2
                                       "HLT_DoubleMu3_LowMass_v1", #5E32, v6.1-v6.2
                                       "HLT_Dimuon6p5_LowMass_v1", #5E32, v8.1-v8.3
                                       "HLT_Dimuon6p5_LowMass_Displaced_v1", #5E32, v8.1-v8.3
                                       "HLT_DoubleMu3_Bs_v1", #5E32, v4.2-v5.3
                                       "HLT_DoubleMu2_Bs_v1", #5E32, v6.1-v6.2
                                       "HLT_DoubleMu2_Bs_v2", #5E32, v8.1-v8.3
                                       "HLT_DoubleMu2_Bs_v3", #1E33, v1.2-v2.4
                                       "HLT_DoubleMu2_Bs_v4", #1E33, v2.5
                                       "HLT_DoubleMu2_Bs_v5", #1.4E33, v1.2
                                       "HLT_Dimuon0_Jpsi_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon0_Jpsi_v2", #1E33, v2.5
                                       "HLT_Dimuon0_Jpsi_v3", #1.4E33, v1.2
                                       "HLT_Dimuon0_Jpsi_v5", #2E33, v1.1-v1.2
                                       "HLT_Dimuon0_Jpsi_NoVertexing_v2", #2E33, v1.1-v1.2
                                       "HLT_Dimuon0_Upsilon_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon0_Upsilon_v2",#1E33, v2.5
                                       "HLT_Dimuon0_Upsilon_v3",#1.4E33, v1.2
                                       "HLT_Dimuon0_Upsilon_v5",#2E33, v1.1-v1.2
                                       "HLT_Dimuon4_Bs_Barrel_v2",#1E33, v1.2-v1.3
                                       "HLT_Dimuon4_Bs_Barrel_v3", #1E33, v2.2-v2.4
                                       "HLT_Dimuon4_Bs_Barrel_v4", #1E33, v2.5
                                       "HLT_Dimuon4_Bs_Barrel_v5",#1.4E33, v1.2
                                       "HLT_Dimuon4_Bs_Barrel_v7",#2E33, v1.1-v1.2
                                       "HLT_Dimuon5_Upsilon_Barrel_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon5_Upsilon_Barrel_v2",#1E33, v2.5
                                       "HLT_Dimuon5_Upsilon_Barrel_v3",#1.4E33, v1.2
                                       "HLT_Dimuon5_Upsilon_Barrel_v5",#2E33, v1.1-v1.2
                                       "HLT_Dimuon6_Bs_v1",#1E33, v1.2-v1.3
                                       "HLT_Dimuon6_Bs_v2",#1E33, v2.2-v2.4
                                       "HLT_Dimuon6_Bs_v3",#1E33, v2.5
                                       "HLT_Dimuon6_Bs_v4",#1.4E33, v1.2
                                       "HLT_Dimuon6_Bs_v6",#2E33, v1.1-v1.2
                                       "HLT_Dimuon7_LowMass_Displaced_v1",#1E33, v1.2-v1.3
                                       "HLT_Dimuon7_LowMass_Displaced_v2", #1E33, v2.2-v2.4
                                       "HLT_Dimuon7_LowMass_Displaced_v3",#1E33, v2.5
                                       "HLT_Dimuon7_LowMass_Displaced_v4",#1.4E33, v1.2
                                       "HLT_DoubleMu4_LowMass_Displaced_v2",#2E33, v1.1-v1.2
                                       "HLT_Dimuon7_Jpsi_Displaced_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon7_Jpsi_Displaced_v2",#1E33, v2.5
                                       "HLT_Dimuon7_Jpsi_Displaced_v3",#1.4E33, v1.2
                                       "HLT_DoubleMu3p5_Jpsi_Displaced_v2", #2E33, v1.1-v1.2
                                       "HLT_Dimuon7_Jpsi_X_Barrel_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon7_Jpsi_X_Barrel_v2",#1E33, v2.5
                                       "HLT_Dimuon7_Jpsi_X_Barrel_v3",#1.4E33, v1.2
                                       "HLT_Dimuon7_Jpsi_X_Barrel_v5",#2E33, v1.1-v1.2
                                       "HLT_Dimuon7_PsiPrime_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon7_PsiPrime_v2",#1E33, v2.5
                                       "HLT_Dimuon7_PsiPrime_v3",#1.4E33, v1.2
                                       "HLT_Dimuon7_PsiPrime_v5",#2E33, v1.1-v1.2
                                       "HLT_Dimuon10_Jpsi_Barrel_v1", #1E33, v1.2-v2.4
                                       "HLT_Dimuon10_Jpsi_Barrel_v2",#1E33, v2.5
                                       "HLT_Dimuon10_Jpsi_Barrel_v3",#1.4E33, v1.2
                                       "HLT_Dimuon10_Jpsi_Barrel_v5",#2E33, v1.1-v1.2
                                       "HLT_Dimuon0_Jpsi_Muon_v1", #1E33, v1.2-v1.3
                                       "HLT_Dimuon0_Jpsi_Muon_v2", #1E33, v2.2-v2.4
                                       "HLT_Dimuon0_Jpsi_Muon_v3",#1E33, v2.5
                                       "HLT_Dimuon0_Jpsi_Muon_v4",#1.4E33, v1.2
                                       "HLT_Dimuon0_Jpsi_Muon_v6", #2E33, v1.1-v1.2
                                       "HLT_Dimuon0_Upsilon_Muon_v1",#1E33, v1.2-v1.3
                                       "HLT_Dimuon0_Upsilon_Muon_v2", #1E33, v2.2-v2.4
                                       "HLT_Dimuon0_Upsilon_Muon_v3",#1E33, v2.5
                                       "HLT_Dimuon0_Upsilon_Muon_v4",#1.4E33, v1.2
                                       "HLT_Dimuon0_Upsilon_Muon_v6",#2E33, v1.1-v1.2
                                       ),
   # ONE FILTER NAME PER PATH    
   HLTLastFilterNames_DoubleMu = cms.vstring("hltDoubleMu3JpsiL3Filtered", #HLT_DoubleMu3_Jpsi_v1
                                             "hltDoubleMu3JpsiL3Filtered", #HLT_DoubleMu3_Jpsi_v2
                                             "hltDimuon6p5JpsiL3Filtered", #HLT_Dimuon6p5_Jpsi_v1
                                             "hltDimuon6p5JpsiDisplacedL3Filtered", #HLT_Dimuon6p5_Jpsi_Displaced_v1
                                             "hltDimuon6p5BarrelJpsiL3Filtered", #HLT_Dimuon6p5_Barrel_Jpsi_v1
                                             "hltDoubleMu3QuarkoniumL3Filtered", #HLT_DoubleMu3_Quarkonium_v1	      
                                             "hltDoubleMu3QuarkoniumL3Filtered", #HLT_DoubleMu3_Quarkonium_v2
                                             "hltDoubleMu3LowMassL3Filtered", #HLT_DoubleMu3_LowMass_v1
                                             "hltDimuon6p5LowMassL3Filtered", #HLT_Dimuon6p5_LowMass_v1
                                             "hltDimuon6p5LowMassL3FilteredDisplaced", #HLT_Dimuon6p5_LowMass_Displaced_v1
                                             "hltDoubleMu3BsL3Filtered", #HLT_DoubleMu3_Bs_v1
                                             "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v1
                                             "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v2
                                             "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v3
                                             "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v4
                                             "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v5
                                             "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v1 was "hltJpsiL3Filtered"
                                             "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v2 was "hltJpsiL3Filtered"
                                             "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v3 was "hltJpsiL3Filtered"
                                             "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v5
                                             "hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v2
                                             "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v1 was "hltUpsilonL3Filtered"
                                             "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v2 was "hltUpsilonL3Filtered"
                                             "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v3 was "hltUpsilonL3Filtered"
                                             "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v5
                                             "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v2
                                             "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v3
                                             "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v4
                                             "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v5
                                             "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v7
                                             "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v1 was "hltBarrelUpsilonL3Filtered"
                                             "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v2 was "hltBarrelUpsilonL3Filtered"
                                             "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v3 was "hltBarrelUpsilonL3Filtered"
                                             "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v5
                                             "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v1
                                             "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v2
                                             "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v3
                                             "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v4
                                             "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v6
                                             "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v1 was "hltLowMassDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v2 was "hltLowMassDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v3 was "hltLowMassDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v4 was "hltLowMassDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterLowMass", #HLT_DoubleMu4_LowMass_Displaced_v2
                                             "hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v1 was "hltJpsiDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v2 was "hltJpsiDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v3 was "hltJpsiDisplacedL3Filtered"
                                             "hltDisplacedmumuFilterJpsi", #HLT_DoubleMu3p5_Jpsi_Displaced_v2
                                             "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v1 was "hltJpsiXBarrelL3Filtered"
                                             "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v2 was "hltJpsiXBarrelL3Filtered"
                                             "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v3 was "hltJpsiXBarrelL3Filtered"
                                             "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v5
                                             "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v1 was "hltPsiPrimeL3Filtered"
                                             "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v2 was "hltPsiPrimeL3Filtered"
                                             "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v3 was "hltPsiPrimeL3Filtered"
                                             "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v5
                                             "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v1 was "hltBarrelJpsiL3Filtered"
                                             "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v2 was "hltBarrelJpsiL3Filtered"
                                             "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v3 was "hltBarrelJpsiL3Filtered"
                                             "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v5
                                             "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v1 was "hltJpsiMuonL3Filtered"
                                             "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v2 was "hltJpsiMuonL3Filtered"
                                             "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v3 was "hltJpsiMuonL3Filtered"
                                             "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v4 was "hltJpsiMuonL3Filtered"
                                             "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v6
                                             "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v1 was "hltUpsilonMuonL3Filtered"
                                             "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v2 was "hltUpsilonMuonL3Filtered"
                                             "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v3 was "hltUpsilonMuonL3Filtered"
                                             "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v4 was "hltUpsilonMuonL3Filtered"
                                             "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v6
                                             ),

    HLTBitNames_MuL2Mu = cms.vstring("HLT_Mu5_L2Mu2_v1", #5E32, v4.2-v5.3
                                     "HLT_Mu5_L2Mu2_Jpsi_v1", #5E32, v4.2-v5.3
                                     "HLT_Mu5_L2Mu2_v2", #5E32, v6.1-v6.2
                                     "HLT_Mu5_L2Mu2_Jpsi_v2", #5E32, v6.1-v6.2
                                     "HLT_Mu5_L2Mu2_Jpsi_v3",#5E32, v8.1-v8.3
                                     "HLT_Mu5_L2Mu2_Jpsi_v4",##1E33, v1.2-v2.4
                                     "HLT_Mu5_L2Mu2_Jpsi_v5",#1E33, v2.5
                                     "HLT_Mu5_L2Mu2_Jpsi_v6",#1.4E33, v1.2
                                     "HLT_Mu5_L2Mu2_Jpsi_v8",#2E33, v1.1-v1.2
                                     ),
    # TWO FILTER NAMES PER PATH (FIRST is L3, SECOND is L2)                           
    HLTLastFilterNames_MuL2Mu = cms.vstring("hltMu5L2Mu0L3Filtered5","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0", 
                                            "hltMu5L2Mu0L3Filtered5","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0",
                                            "hltMu5L2Mu2JpsiTrackMassFiltered","hltMu5L2Mu2L2PreFiltered0",
                                            ),

    HLTBitNames_MuTrack = cms.vstring("HLT_Mu3_Track3_Jpsi_v4", #5E32, v4.2-v5.3
                                      "HLT_Mu7_Track5_Jpsi_v1", #5E32, v4.2-v5.3
                                      "HLT_Mu7_Track7_Jpsi_v1", #5E32, v4.2-v5.3    
                                      "HLT_Mu3_Track3_Jpsi_v5", #5E32, v6.1-v6.2
                                      "HLT_Mu5_Track2_Jpsi_v1", #5E32, v6.1-v6.2
                                      "HLT_Mu7_Track5_Jpsi_v2", #5E32, v6.1-v6.2
                                      "HLT_Mu7_Track7_Jpsi_v2", #5E32, v6.1-v6.2                             
                                      "HLT_Mu5_Track2_Jpsi_v2", #5E32, v8.1-v8.3
                                      "HLT_Mu7_Track7_Jpsi_v3", #5E32, v8.1-v8.3
                                      "HLT_Mu5_Track2_Jpsi_v4", #1E33, v1.2-v2.4
                                      "HLT_Mu7_Track7_Jpsi_v5", #1E33, v1.2-v2.4
                                      "HLT_Mu5_Track2_Jpsi_v5", #1E33, v2.5
                                      "HLT_Mu7_Track7_Jpsi_v6", #1E33, v2.5
                                      "HLT_Mu5_Track2_Jpsi_v6", #1.4E33, v1.2
                                      "HLT_Mu7_Track7_Jpsi_v7", #1.4E33, v1.2

                                      "HLT_Mu5_Track2_Jpsi_v8", #2E33, v1.1-v1.2
                                      "HLT_Mu7_Track7_Jpsi_v9", #2E33, v1.1-v1.2
                                      ), 
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_MuTrack = cms.vstring("hltMu3Track3JpsiTrackMassFiltered",
                                             "hltMu7Track5JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",  
                                             "hltMu3Track3JpsiTrackMassFiltered",
                                             "hltMu5Track2JpsiTrackMassFiltered",
                                             "hltMu7Track5JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",
                                             "hltMu5Track2JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",
                                             "hltMu5Track2JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",
                                             "hltMu5Track2JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",
                                             "hltMu5Track2JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",
                                             "hltMu5Track2JpsiTrackMassFiltered",
                                             "hltMu7Track7JpsiTrackMassFiltered",
                                             ),
                              
    HLTBitNames_MuTkMu = cms.vstring(),
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_MuTkMu = cms.vstring(),

    ####### 2010, ALL ########

##     HLTBitNames_SingleMu = cms.vstring("HLT_Mu3", 
##                                        "HLT_Mu5", 
##                                        "HLT_Mu7", 
##                                        "HLT_Mu9", 
##                                        "HLT_Mu11"),
##     # ONE FILTER NAME PER PATH
##     HLTLastFilterNames_SingleMu = cms.vstring("hltSingleMu3L3Filtered3",	
##                                               "hltSingleMu5L3Filtered5",	
##                                               "hltSingleMu7L3Filtered7",	
##                                               "hltSingleMu9L3Filtered9",	
##                                               "hltSingleMu11L3Filtered11"),

##     HLTBitNames_DoubleMu = cms.vstring("HLT_DoubleMu0",		      
##                                        "HLT_DoubleMu0_Quarkonium_v1", 
##                                        "HLT_DoubleMu0_Quarkonium_LS_v1", 
##                                        "HLT_L1DoubleMuOpen",	      
##                                        "HLT_L1DoubleMuOpen_Tight",    
##                                        "HLT_DoubleMu3"), 
##     # ONE FILTER NAME PER PATH    
##     HLTLastFilterNames_DoubleMu = cms.vstring("hltDiMuonL3PreFiltered0",
##                                               "hltDoubleMu0QuarkoniumL3PreFiltered",
##                                               "hltDoubleMu0QuarkoniumLSL3PreFiltered",
##                                               "hltDoubleMuLevel1PathL1OpenFiltered",
##                                               "hltL1DoubleMuOpenTightL1Filtered",
##                                               "hltDiMuonL3PreFiltered"),

##     HLTBitNames_MuL2Mu = cms.vstring("HLT_Mu5_L2Mu0"), 
##     # TWO FILTER NAMES PER PATH (FIRST is L3, SECOND is L2)                           
##     HLTLastFilterNames_MuL2Mu = cms.vstring("hltMu5L2Mu0L3Filtered5","hltDiMuonL2PreFiltered0"),

##     HLTBitNames_MuTrack = cms.vstring("HLT_Mu0_Track0_Jpsi",
##                                       "HLT_Mu3_Track0_Jpsi",
##                                       "HLT_Mu5_Track0_Jpsi",
##                                       "HLT_Mu5_Track0_Jpsi_v2",
##                                       "HLT_Mu3_Track3_Jpsi",
##                                       "HLT_Mu3_Track3_Jpsi_v2",
##                                       "HLT_Mu3_Track5_Jpsi_v1",
##                                       "HLT_Mu3_Track5_Jpsi_v2"),
##     # ONE FILTER NAME PER PATH
##     HLTLastFilterNames_MuTrack = cms.vstring("hltMu0TrackJpsiTrackMassFiltered",
##                                              "hltMu3TrackJpsiTrackMassFiltered",
##                                              "hltMu5TrackJpsiTrackMassFiltered",
##                                              "hltMu5TrackJpsiTrackMassFiltered",
##                                              "hltMu3Track3JpsiTrackMassFiltered",
##                                              "hltMu3Track3JpsiTrackMassFiltered",
##                                              "hltMu3Track5JpsiTrackMassFiltered",
##                                              "hltMu3Track5JpsiTrackMassFiltered"),
                              
##     HLTBitNames_MuTkMu = cms.vstring("HLT_Mu0_TkMu0_Jpsi",	      
##                                      "HLT_Mu3_TkMu0_Jpsi",	      
##                                      "HLT_Mu5_TkMu0_Jpsi",	           
##                                      "HLT_Mu0_TkMu0_OST_Jpsi",     
##                                      "HLT_Mu3_TkMu0_OST_Jpsi",      
##                                      "HLT_Mu5_TkMu0_OST_Jpsi",      	      
##                                      "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1",
##                                      "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1",
##                                      "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1",       
##                                      "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2",
##                                      "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2",
##                                      "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2",      
##                                      "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3",
##                                      "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3",
##                                      "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3"),
##     # ONE FILTER NAME PER PATH
##     HLTLastFilterNames_MuTkMu = cms.vstring("hltMu0TkMuJpsiTkMuMassFiltered",
##                                             "hltMu3TkMuJpsiTkMuMassFiltered", 
##                                             "hltMu5TkMuJpsiTkMuMassFiltered",
##                                             "hltMu0TkMuJpsiTkMuMassFiltered",
##                                             "hltMu3TkMuJpsiTkMuMassFiltered",
##                                             "hltMu5TkMuJpsiTkMuMassFiltered",
##                                             "hltMu0TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu3TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu5TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu0TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu3TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu5TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu0TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu3TkMuJpsiTkMuMassFilteredTight",
##                                             "hltMu5TkMuJpsiTkMuMassFilteredTight"),
                              
)

## no filter
# process.p = cms.Path(process.demo)

## filter on vertex
process.p = cms.Path(process.primaryVertexFilter*process.demo)

import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
#myLumis = LumiList.LumiList(filename = 'Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys.txt').getCMSSWString().split(',')
myLumis = LumiList.LumiList(filename = 'Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)
## filter on vertex and HLT
# process.p = cms.Path(process.primaryVertexFilter*process.hltMuF*process.demo)
