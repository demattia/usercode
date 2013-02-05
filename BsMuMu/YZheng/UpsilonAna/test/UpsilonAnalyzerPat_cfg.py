import FWCore.ParameterSet.Config as cms

process = cms.Process("ANAPAT")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1616) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_9_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_99_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_98_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_97_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_96_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_95_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_94_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_93_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_92_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_91_1.root',
        '/store/user/yzheng/Upsilon3S/PAT_Summer09_Upsilon3S_8E29/e18a91fddaa26d30b8f9cc8915e89f8f/PAT_ups3s_90_1.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('tree.root')
)

process.demo = cms.EDAnalyzer('UpsilonAnalyzerPAT',
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    srcWithCaloMuons = cms.InputTag("onia2MuMuPatGlbCal"),
    whichsample = cms.int32(4),                          
    pTBinRanges = cms.vdouble(3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 14.0, 18.0, 25.0, 35.0, 60.0),
    etaBinRanges = cms.vdouble(0.0, 1.1, 2.5), 
    onlyTheBest = cms.bool(False),		
    applyCuts = cms.bool(True),# apply qualifty cuts			
    storeEfficiency = cms.bool(False),	
    useBeamSpot = cms.bool(True),
    useCaloMuons = cms.untracked.bool(False),
    removeSignalEvents = cms.untracked.bool(False),
#    TriggerResultsLabel = cms.InputTag("TriggerResults","","HLT")           
    # For Summer09 only 
    TriggerResultsLabel = cms.InputTag("TriggerResults","","HLT8E29")
)


process.p = cms.Path(process.demo)
