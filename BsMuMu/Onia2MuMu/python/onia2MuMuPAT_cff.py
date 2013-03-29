import FWCore.ParameterSet.Config as cms

def onia2MuMuPAT(process, GlobalTag, MC=False, HLT='HLT', Filter=True):
    # Setup the process
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        # fileMode = cms.untracked.string('MERGE'),
    )
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 100
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = GlobalTag

    # Drop the DQM stuff on input
    process.source = cms.Source("PoolSource",
        inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*"),
        fileNames = cms.untracked.vstring()
    )

    # Scraping filter
    process.scrapingFilter = cms.EDFilter("FilterOutScraping",
        applyfilter = cms.untracked.bool(True),
        debugOn = cms.untracked.bool(False),
        numtrack = cms.untracked.uint32(10),
        thresh = cms.untracked.double(0.25)
    )

    # Merge muons
    process.mergedMuons = cms.EDProducer("CaloMuonMerger",
        muons     = cms.InputTag("muons"), 
        muonsCut = cms.string(""),
        mergeCaloMuons = cms.bool(False),  ### NEEDED TO RUN ON AOD
        caloMuons = cms.InputTag("calomuons"),
        caloMuonsCut = cms.string(""),
        minCaloCompatibility = cms.double(0.6),
        mergeTracks = cms.bool(True),
        tracks = cms.InputTag("generalTracks"),
        # tracksCut = cms.string("pt > 0.5")
        tracksCut = cms.string("(abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <=2.4 && pt > 1.5)"),
    )

    # Prune generated particles to muons and their parents
    process.genMuons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *  ",                     # this is the default
            "++keep abs(pdgId) = 13",        # keep muons and their parents
            "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
        )
    )

    # Make PAT Muons
    process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
    # with some customization
    if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True
        process.muonMatch.matched = "genMuons"
    changeRecoMuonInput(process, "mergedMuons")
    changeTriggerProcessName(process, HLT)
    switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
    #useL1MatchingWindowForSinglets(process)
    process.muonL1Info.maxDeltaR     = 0.3 
    process.muonL1Info.fallbackToME1 = True
    process.muonMatchHLTL1.maxDeltaR = 0.3
    process.muonMatchHLTL1.fallbackToME1 = True
    process.muonMatchHLTL2.maxDeltaR = 0.3
    process.muonMatchHLTL2.maxDPtRel = 10.0
    process.muonMatchHLTL3.maxDeltaR = 0.1
    process.muonMatchHLTL3.maxDPtRel = 10.0
    process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
    process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
    process.muonMatchHLTTrackMu.maxDeltaR = 0.1
    process.muonMatchHLTTrackMu.maxDPtRel = 10.0

    # Make a sequence
    process.patMuonSequence = cms.Sequence(
        process.scrapingFilter *
        process.mergedMuons *
        process.genMuons *
        process.patMuonsWithTriggerSequence
    )
    if not MC:
        process.patMuonSequence.remove(process.genMuons)

    # Make dimuon candidates
    process.onia2MuMuPatTrkTrk = cms.EDProducer('Onia2MuMuPAT',
        muons = cms.InputTag("patMuonsWithTrigger"),
        beamSpotTag = cms.InputTag("offlineBeamSpot"),
        primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
        # At least one muon must pass this selection
        higherPuritySelection  = cms.string("(isGlobalMuon || isTrackerMuon) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35"),
        # BOTH muons must pass this selection
        lowerPuritySelection  = cms.string("(isGlobalMuon || isTrackerMuon) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35"),
        # The dimuon must pass this selection before vertexing
        dimuonSelection  = cms.string("4.5 < mass && mass < 6.3"),
        # Embed the full reco::Vertex out of the common vertex fit
        addCommonVertex = cms.bool(True),
        # Embed the primary vertex re-made from all the tracks except the two muons
        addMuonlessPrimaryVertex = cms.bool(True),
        # Add the common MC mother of the two muons, if any
        addMCTruth = cms.bool(MC),
        # Order PVs by their vicinity to the J/psi vertex, not by sumPt
        resolvePileUpAmbiguity = cms.bool(True),
        addThirdTrack = cms.bool(True),  ## search a third track making a good vertex with the dimuon
        minTrackPt =cms.double(0.5), ## minimum pt of the third track
        trackMass = cms.double(0.4936), ## mass for the track
        diMuPlusTrackMassMax = cms.double(5.9),
        diMuPlusTrackMassMin = cms.double(4.9),
        diMuMassMax = cms.double(3.2),## dimuon mass range where to search for the third track
        diMuMassMin = cms.double(2.8),
        # Preselection cuts
        preselection = cms.string("mass > 4.5 && mass < 6.3 && pt > 5. && pt < 9999. && daughter('muon1').pt > 4. && daughter('muon1').pt < 999. && daughter('muon2').pt > 4. && daughter('muon2').pt < 999. && userFloat('l3d') < 2. && userFloat('l3dsig') > 0. && userFloat('l3dsig') < 120. && userFloat('vNChi2') < 10. && userFloat('delta3d') < 0.1 && userFloat('delta3d')/userFloat('delta3dErr') < 5. && userFloat('DCA') < 0.1 && acos(userFloat('cosAlpha3D')) < 0.3 && userInt('Ntrk') < 21 && userFloat('minDca') < 0.25 && userFloat('Isolation') > 0."),
    )

    # check if there is at least one (inclusive) tracker+tracker di-muon
    process.onia2MuMuPatTrkTrkFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('onia2MuMuPatTrkTrk'),
        minNumber = cms.uint32(1),
    )

    process.onia2MuMuPatMassFilter = cms.EDFilter("CandViewSelector",
        src = cms.InputTag('onia2MuMuPatTrkTrk'),
        cut = cms.string("mass > 4.5 && mass < 6.3 && pt > 5. && pt < 9999. && daughter('muon1').pt > 4. && daughter('muon1').pt < 999. && daughter('muon2').pt > 4. && daughter('muon2').pt < 999. && userFloat('l3d') < 2. && userFloat('l3dsig') > 0. && userFloat('l3dsig') < 120. && userFloat('vNChi2') < 10. && userFloat('delta3d') < 0.1 && userFloat('delta3d')/userFloat('delta3dErr') < 5. && userFloat('DCA') < 0.1 && acos(userFloat('cosAlpha3D')) < 0.3 && userInt('Ntrk') < 21 && userFloat('minDca') < 0.25 && userFloat('Isolation') > 0."),
        # cut = cms.string("mass > 4.5 && mass < 6.3 && pt > 5. && pt < 9999. && daughter('muon1').pt > 4. && daughter('muon1').pt < 999. && daughter('muon2').pt > 4. && daughter('muon2').pt < 999."),
        filter = cms.bool(True)
    )

    # add tracks
    process.load("PhysicsTools.RecoAlgos.allTracks_cfi")
    process.load("PhysicsTools.RecoAlgos.goodTracks_cfi")
    process.goodTracks.particleType = 'pi+'
    # process.goodTracks.cut = 'charge !=0 && found > 1 && pt > 0.5 && hitPattern.pixelLayersWithMeasurement > 0 && abs(d0) < 3.0 && abs(dz) < 25.0 && chi2/ndof < 4.0'
    process.goodTracks.cut = 'charge !=0 && pt > 0.5'

    process.nEventsTotal = cms.EDProducer("EventCountProducer")
    process.nEventsAfterFilter = cms.EDProducer("EventCountProducer")
    
    # the onia2MuMu path
    process.Onia2MuMuPAT = cms.Path(
        process.nEventsTotal *
        process.patMuonSequence *
        process.onia2MuMuPatTrkTrk *
        process.allTracks *
        process.goodTracks *
        process.onia2MuMuPatTrkTrkFilter *
        process.onia2MuMuPatMassFilter *
        process.nEventsAfterFilter
    )

    # output
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('onia2MuMuPAT.root'),
        outputCommands = cms.untracked.vstring('drop *',
            'keep *_genMuons_*_Onia2MuMuPAT',                      # generated muons and parents
             # 'keep patMuons_patMuonsWithTrigger_*_Onia2MuMuPAT',    # All PAT muos including general tracks and matches to triggers
            'keep *_goodTracks_*_Onia2MuMuPAT',                    # All good tracks (heavy!)  
             # 'keep patCompositeCandidates_*__Onia2MuMuPAT',       # PAT di-muons
            'keep patCompositeCandidates_*_*_Onia2MuMuPAT',        # PAT di-muons and dimuons+track
            # 'keep patMuons_tagMuons__Onia2MuMuPAT',                # tagMuons for efficiency
            # 'keep patMuons_probeMuons__Onia2MuMuPAT',              # probeMuons for efficiency
            # 'keep *_tagMuonsMCMatch__Onia2MuMuPAT',                # tagMuons MC matches for efficiency
            # 'keep *_probeMuonsMCMatch__Onia2MuMuPAT',              # probeMuons MC matches for efficiency
            'keep recoCompositeCandidates_*__Onia2MuMuPAT',        # RECO di-muons, tpPairs for efficiency
            'keep *_offlinePrimaryVertices_*_*',                   # Primary vertices: you want these to compute impact parameters
            'keep *_offlineBeamSpot_*_*',                          # Beam spot: you want this for the same reason                                   
            'keep edmTriggerResults_TriggerResults_*_*',           # HLT info, per path (cheap)
            'keep l1extraL1MuonParticles_l1extraParticles_*_*',    # L1 info (cheap)
            'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',       # Prescale info
            'keep *_l1GtRecord_*_*',                               # Prescale info
            'keep edmMergeableCounter_*_*_*'                       # Count the number of events before and after filters
            #'keep *_*_DiMu_*'
        ),
        # outputCommands = cms.untracked.vstring('keep *'),
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('Onia2MuMuPAT') ) if Filter else cms.untracked.PSet()
    )
    process.e = cms.EndPath(process.out)

