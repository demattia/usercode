import FWCore.ParameterSet.Config as cms

process = cms.Process("ANAPAT")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         'file:/scratch/scratch95/z/zheng13/myfile.root',
    ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

# process.options = cms.untracked.PSet(
#    IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
# )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ups1sYield.root')
)

#PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % ("HLT_Mu3",);

#process.triggerMuons = cms.EDFilter("PATMuonRefSelector",
#    src = cms.InputTag("patMuons"),
#    cut = cms.string( "abs(eta) < 2.5 && " + PASS_HLT ),
#    filter = cms.bool(True),
#)

process.demo = cms.EDAnalyzer('JPsiAnalyzerPAT',
    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    srcWithCaloMuons = cms.InputTag("onia2MuMuPatGlbCal"),
    whichsample = cms.int32(1), # =1:Y(1S), 2:Y(2S), 3:Y(3S), 0:ppMuMuX                         
    pTBinRanges = cms.vdouble(3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 14.0, 18.0, 25.0, 35.0, 60.0),
    etaBinRanges = cms.vdouble(0.0, 1.1, 2.5), 
    onlyTheBest = cms.bool(False),		
    applyCuts = cms.bool(True),			
    storeEfficiency = cms.bool(False),	
    useBeamSpot = cms.bool(True),
    useCaloMuons = cms.untracked.bool(False),
    removeSignalEvents = cms.untracked.bool(False),
    TriggerResultsLabel = cms.InputTag("TriggerResults","","HLT8E29")           
    # For Summer09 only 
    # TriggerResultsLabel = cms.InputTag("TriggerResults","","HLT8E29")
)


process.p = cms.Path(process.demo)
