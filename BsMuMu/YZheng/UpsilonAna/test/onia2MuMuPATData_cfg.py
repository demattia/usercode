import FWCore.ParameterSet.Config as cms

process = cms.Process("SkimmingOnia2MuMuPAT")

### standard includes
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10

### global tag
process.GlobalTag.globaltag = "GR_R_35X_V7::All"

process.options = cms.untracked.PSet(
  fileMode = cms.untracked.string('MERGE')
)

### source
process.source = cms.Source("PoolSource",
    inputCommands = cms.untracked.vstring("keep *",
        "drop *_MEtoEDMConverter_*_*"
    ),
    fileNames = cms.untracked.vstring(
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/952/66AD7FB3-5951-DF11-A9D2-00E08178C13D.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/937/AC3EDF7B-2851-DF11-A0EA-0025B3E066A0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/935/320855C6-6D51-DF11-8CAD-003048D45FD4.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/934/F47623B3-5951-DF11-96B7-00E08178C01F.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/933/4A53E9C1-5751-DF11-AB1E-00E08178C13D.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/932/1215AA36-5C51-DF11-A386-00E08178C10B.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/929/465D94C9-2151-DF11-8594-003048635FA8.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/928/98190C57-0852-DF11-9FCC-003048D476B0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/928/801750FB-0C52-DF11-8490-0025B3E05D6E.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/928/68B374E0-6651-DF11-89B4-0025B31E3330.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/927/74A7FF55-0251-DF11-A5EE-002481E0D6A0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/926/EE051C40-0551-DF11-9160-0015170ACB6C.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/925/A6FAD20F-FA50-DF11-82EB-003048D479BE.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/924/389CCCC2-0051-DF11-A829-0025B3E05C32.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/922/AABDB939-FF50-DF11-B54C-00E0817918B9.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/919/F0A6805D-D150-DF11-A522-00E08179185B.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/908/943761D9-C750-DF11-AE75-003048D47A1A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/907/10599E5F-9450-DF11-9551-003048D45FBE.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/899/109B40A7-CC50-DF11-B61B-00E08177F267.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/894/E8A462F7-CD50-DF11-8D9B-003048D45FBE.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/887/F2A66D4B-D150-DF11-98B3-0025B3E06378.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/885/E851089D-D050-DF11-A36B-001A6478AB88.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/885/B477E59E-D250-DF11-A4E4-00E0817918D9.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/881/EE1FF30F-8B50-DF11-9F38-0025B3E065F6.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/881/580B9EBA-8F50-DF11-8EC1-003048D4608C.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/ECB3D2F3-8050-DF11-959B-003048673E84.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/D4115E3C-8050-DF11-B70C-003048673F3A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/A671A91F-C550-DF11-B292-00E08178C037.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/A4A1FB94-CA50-DF11-A409-00E081B08DA7.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/7C2F1B45-BD50-DF11-80CC-00E0817917DF.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/7AC303A1-7D50-DF11-880F-003048D45FD6.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/4E9714E7-7E50-DF11-8493-003048D45F4E.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/877/462AF200-D450-DF11-BC0B-00E08179182F.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/876/9E5492DF-5A50-DF11-9A6A-00E0817917EB.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/876/6CB47034-5950-DF11-9782-00E0817917EB.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/875/E8CC3DB7-0650-DF11-A70E-001A649747B0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/874/CAC6F02D-ED4F-DF11-A3FC-00E081B08D11.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/874/A6FCF8EE-F34F-DF11-9B15-003048D45F62.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/874/925C19D7-EF4F-DF11-9640-00E081B08D11.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/874/66963B5D-0750-DF11-99F0-003048670B64.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/874/4E07272D-E94F-DF11-AB55-002481E14D84.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/873/AA18D331-674F-DF11-A470-003048D47A24.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/872/503374A4-6F4F-DF11-AC91-003048D47742.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/868/2E462132-674F-DF11-A538-003048D47756.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/867/8CEB3FB4-5E4F-DF11-B6BD-0025B31E3C24.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/865/722F9234-674F-DF11-9C5D-003048D4607A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/864/D0D76A6E-624F-DF11-95A0-00E0817917E1.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/863/D8B413EE-684F-DF11-90D1-003048D45FBC.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/853/AC198885-384F-DF11-BA85-001A64789D78.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/849/14F382B4-404F-DF11-8947-00E08178C189.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/848/06D0B3FB-2B4F-DF11-99C1-001A64789D30.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/847/7CE679A0-2A4F-DF11-8104-002481E14FEE.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/836/F6A2F34C-1E4F-DF11-A611-00E0817917CF.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/832/DA4B79A2-104F-DF11-8C39-0025B3E05D3E.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/826/1CF2AB28-1A4F-DF11-9025-0025B3E05BA8.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/824/AC4F6989-2C4F-DF11-93AE-003048D460F8.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/820/20E39E48-2B4F-DF11-AC22-003048D47786.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/819/A8A26DFA-2B4F-DF11-911A-003048D460F8.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/817/624F51ED-2F4F-DF11-95B0-00E08178C189.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/816/92FA1E09-2A4F-DF11-824D-0025B3E06514.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/815/5A47F4EF-AE4E-DF11-B8DE-00E08178C107.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/814/F8BC870B-AB4E-DF11-ABCD-00E08178C0F0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/807/78DE2A09-E14E-DF11-972A-003048D47750.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/805/DA2A6AB4-844E-DF11-8385-00E08178C195.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/790/74B7517A-924E-DF11-8E55-00E08178C0F0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/788/D231E29D-B54E-DF11-9F12-00E08178C0AF.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/784/8E7D1B6D-824E-DF11-807A-00E0817917D1.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/776/126A444F-D44E-DF11-A7CB-00E08178C163.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/769/E693CE47-AE4E-DF11-BAD8-00E08178C107.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/768/C65AEB48-BC4E-DF11-8181-0030483268A8.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/767/F225AB4F-AE4E-DF11-9739-003048D460F4.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/765/CE1FF4FF-C04E-DF11-889B-00E081791763.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/763/6641BD98-AF4E-DF11-B29A-00E08178C185.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/762/60F628BB-E54E-DF11-A242-003048D45F94.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/761/8A392747-994E-DF11-A98A-003048D4DBFC.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/759/5A040528-814E-DF11-8E0D-00E0817918CB.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/758/B88E9C19-8D4E-DF11-842A-00E0817917A3.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/757/784205E5-C14E-DF11-AD7B-00E0817917A3.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/753/C2BF5342-184E-DF11-937C-00E0817917A7.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/750/34C524AC-564E-DF11-B1A7-003048D4773A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/748/48D71B75-394E-DF11-A819-00E081791737.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/747/1CDFAEEF-594E-DF11-B64B-003048673FEA.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/733/FC5928B7-A94E-DF11-8AE8-00E08178C16B.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/710/8A822A30-F14D-DF11-9E6E-00E081791747.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/697/6AFCDD54-AB4D-DF11-8124-003048D45F5A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/696/DC88D4CE-D24D-DF11-940A-003048D4777A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/692/3066E273-EE4D-DF11-8972-00E0817917FD.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/684/7605F8B4-9E4D-DF11-9947-003048D460AC.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/664/D0AE1149-984D-DF11-897D-003048D46008.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/654/3E91B034-814D-DF11-867A-001A64789504.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/643/DA181A39-724D-DF11-AA99-003048D47A78.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/636/8E6739A8-724D-DF11-B3BE-003048D47A78.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/618/981D3A5A-AF4D-DF11-AEED-003048D476B0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/617/AEDC2AB7-A04D-DF11-96C6-001A64789D28.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/616/A6D1011D-AE4D-DF11-B34D-003048D476B0.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/615/2AFD04E2-634D-DF11-AC13-003048636240.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/611/02396A4C-BC4D-DF11-9FD8-003048D45F5A.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/601/16F22F70-604D-DF11-8340-00E08178C14F.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/598/4631000A-9C4C-DF11-9F46-003048D476DC.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/596/0C686E38-074D-DF11-AB02-003048D47A44.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/589/9EA8BCE4-AD4C-DF11-BAD7-003048D4775E.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/588/BEEF9670-C24C-DF11-9866-00E08178C143.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/585/3851BAEB-934C-DF11-94E0-00E08178C0E7.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/580/6C341A0E-A54C-DF11-BC84-00E081791891.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/578/3E82EF19-904C-DF11-90A2-00E08178C021.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/572/0C964262-9F4C-DF11-A2CB-001A64789E60.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/561/9CCBAF1A-9B4C-DF11-BCB7-001A64789E60.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/559/A4B2DD48-AF4C-DF11-8284-00E0817917DB.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/558/204F0C96-844C-DF11-B925-00E081791875.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/556/5C7244EE-854C-DF11-9950-00E081791875.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/548/D86BE709-844C-DF11-A1E6-00E081791875.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/547/F83F0D37-834C-DF11-B83B-00E081791821.root',
       '/store/data/Commissioning10/MinimumBias/USER/v9/000/133/537/3889E423-864C-DF11-A336-00E08178C097.root'
   )
)

### number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#############REAL DATA######################
# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
# set the L1MinBiasVetoBeamHalo bit
process.L1MinBiasVetoBeamHalo = cms.Path(process.hltLevel1GTSeed)

# this is for filtering on HLT physics bit
process.hltPhysicsDeclared = cms.EDFilter("HLTHighLevel",
                                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                 HLTPaths = cms.vstring("HLT_PhysicsDeclared"
                                                        ),
                                 eventSetupPathsKey = cms.string(''),
                                 andOr = cms.bool(True),
                                 throw = cms.bool(False)  # Avoid crashes if the job happens to process runs that don't have the HLT_PhysicsDeclared (e.g. 122294)
                                 )
# copy the PhysicsDeclared bit
process.PhysicsDeclared = cms.Path(process.hltPhysicsDeclared)

# this is for filtering on HLT MinBiasBSC bit
process.hltMinBiasBSC = cms.EDFilter("HLTHighLevel",
                                     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                     HLTPaths = cms.vstring("HLT_MinBiasBSC"
                                                            ),
                                     eventSetupPathsKey = cms.string(''),
                                     andOr = cms.bool(True),
                                     throw = cms.bool(True)
                                     )
# copy the MinBiasBSC bit
process.MinBiasBSC = cms.Path(process.hltMinBiasBSC)

# filter to remove scraping ("monster") events
process.scrapingFilter = cms.EDFilter("FilterOutScraping",
                                      applyfilter = cms.untracked.bool(True),
                                      debugOn = cms.untracked.bool(False),
                                      numtrack = cms.untracked.uint32(10),
                                      thresh = cms.untracked.double(0.25)
                                      )
# set the Scraping bit
process.Scraping = cms.Path(process.scrapingFilter)

# filter on good vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15),
                                           maxd0 = cms.double(2)
                                           )
# set the PrimaryVertex bit
process.PrimaryVertex = cms.Path(process.primaryVertexFilter)

###############################################



### load Onia2MuMu PAT
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff")
# it will define:
#   - patMuonSequence: makes the 'patMuons' for Onia2MuMu
#   - onia2MuMuPat<X>: makes di-muons. <X> can be GlbGlb, GlbTrk, GlbCal, TrkTrk, TrkCal.
#                      categories are inclusive: Cal includes Trk which includes Glb.
#                      just run the loosest one you need.

### skim requiring at least one (inclusive) tracker+tracker di-muon
process.selectedEvents = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('onia2MuMuPatTrkTrk'),
    minNumber = cms.uint32(1),
) 

### path
process.Onia2MuMuPatTrkTrk = cms.Path(
        # process.hltLevel1GTSeed +
        # process.hltPhysicsDeclared +
        # process.hltMinBiasBSC +
        # process.primaryVertexFilter +
        # process.scrapingFilter +
        process.patMuonSequence +     # produce PAT muons for Onia2MuMu (includes merging with CaloMuons)
        process.onia2MuMuPatTrkTrk +  # make J/Psi's (inclusively down to tracker+tracker)
        process.selectedEvents        # select events with J/Psi's
)


### output
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('onia2MuMuPATData.root'),
    outputCommands = cms.untracked.vstring('drop *',
        'keep patCompositeCandidates_*__SkimmingOnia2MuMuPAT', ## PAT di-muons
        'keep patMuons_patMuons__SkimmingOnia2MuMuPAT',        ## All PAT muons (note: not necessary if you use only the di-muons)
        'keep *_offlinePrimaryVertices_*_*',                   ## Primary vertices: you want these to compute impact parameters
        'keep *_offlineBeamSpot_*_*',                          ## Beam spot: you want this for the same reason
        'keep edmTriggerResults_TriggerResults_*_*',           ## HLT info, per path (cheap)
        'keep l1extraL1MuonParticles_l1extraParticles_*_*',    ## L1 info (cheap)
        #'keep *_patTrigger_*_*',                               ## HLT info, per object (BIG. Keep only when debugging trigger match)
    ),
    ## Uncomment to activate skimming!
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('Onia2MuMuPatTrkTrk') )
)
process.e = cms.EndPath(process.out)

##This is real data:
from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import onia2MuMu_isNotMC
onia2MuMu_isNotMC(process)

process.load("YZheng.UpsilonAna.eventDisplayJpsi_cff")
process.load("YZheng.UpsilonAna.eventDisplayUpsilon_cff")
process.load("YZheng.UpsilonAna.selection_cff")
process.TFileService = cms.Service("TFileService", fileName = cms.string("trees.root"))


## If this is Summer09 MC
#process.patTrigger.processName    = "HLT8E29" 

