import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingEfficiencyFromCosmics")

# process.load("FWCore.MessageService.MessageLogger_cfi")
# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('TrackingEfficiencyFromCosmics.root')
                                 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/CMSSW_4_2_2/src/reco_RAW2DIGI_L1Reco_RECO_DQM.root'
        # 'file:hitRemoval_0.root'
        # 'file:hitRemoval.root'
    )
)

process.demo = cms.EDAnalyzer('TrackingEfficiencyFromCosmics',
                              MaxDeltaR = cms.double(1000),
                              SimMaxDeltaR = cms.double(1000),
                              UseMCtruth = cms.bool(False)
)


process.p = cms.Path(process.demo)
