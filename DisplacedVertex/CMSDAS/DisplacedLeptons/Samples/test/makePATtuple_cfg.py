# the following is a standard PAT production template to start from
from PhysicsTools.PatAlgos.patTemplate_cfg import *


########################################
# INPUT
########################################

# data sample and sample-specific settings
sampleRequireCollision=1
sampleSignalPID=0
from DisplacedLeptons.Samples.samples.Debug_sample_cff import *
rereco=1
prefilter=1

process.source.fileNames = sampleRecoFiles
process.source.duplicateCheckMode = cms.untracked.string(sampleDuplicateCheckMode)
process.GlobalTag.globaltag = cms.string(sampleGlobalTag)
process.maxEvents.input = cms.untracked.int32(-1)

# select high quality reco tracks
import CommonTools.RecoAlgos.recoTrackRefSelector_cfi
process.trackSel = CommonTools.RecoAlgos.recoTrackRefSelector_cfi.recoTrackRefSelector.clone(
    src = 'generalTracks',
    quality = ['highPurity'],
    ptMin = 1.0,
    tip = 30.0,
    lip = 30.0,
    minHit = 6,
    min3DHit = 2
)

# filtering out events that do not have at least two lepton candidates
if prefilter:
    process.load("DisplacedLeptons.Samples.isoTrackPrefilter_cfi")
    process.isoTrackPrefilter.signalPDGId = sampleSignalPID

# setup MC matching
process.muonMatch.checkCharge = False
process.muonMatch.maxDPtRel = 1000.
process.electronMatch.checkCharge = False

from PhysicsTools.PatAlgos.tools.coreTools import *

# for running on real data
if sampleType!="MC":
    removeMCMatching(process, ['All'])


if rereco:
    ########################################
    # RECO CONFIGURATION: MUONS
    ########################################

    # disable parts of PAT embedding that won't work here
    process.patMuons.embedCaloMETMuonCorrs = cms.bool(False)
    process.patMuons.embedTcMETMuonCorrs = cms.bool(False)

    # redo muons without primary vertex constraint
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.globalMuons.MuonCollectionLabel = cms.InputTag("standAloneMuons")
    process.globalSETMuons.MuonCollectionLabel = cms.InputTag("standAloneSETMuons")
    process.muons.inputCollectionLabels = cms.VInputTag(cms.InputTag("generalTracks"),
                                                        cms.InputTag("globalMuons"), 
                                                        cms.InputTag("standAloneMuons"))

    ########################################
    # RECO CONFIGURATION: ELECTRONS
    ########################################

    # need to rerun tracking. this sequence courtesy of Sam Harper
    process.pfTrackElec.useFifthStepForTrackerDrivenGsf = True


    ########################################
    # PROCESS
    ########################################

    process.mainSequence = cms.Sequence(
        process.siPixelRecHits
        *process.siStripMatchedRecHits
        *process.particleFlowCluster
        *process.offlineBeamSpot
        *process.recopixelvertexing
        *process.ckftracks
        *process.caloTowersRec
        *process.vertexreco
        *process.recoJets
        #*process.muonrecoComplete
        *process.electronGsfTracking
        *process.metreco
        *process.particleFlowReco
        *process.recoPFJets
        *process.recoPFMET
        *process.PFTau
        #*process.ckftracks
        #*process.ecalDrivenElectronSeeds
        #*process.trackerDrivenElectronSeeds
        #*process.electronMergedSeeds
        #*process.electronCkfTrackCandidates
        #*process.electronGsfTracks
        *process.pfElectronTranslator
        *process.ecalDrivenGsfElectronCores
        *process.ecalDrivenGsfElectrons
        *process.gsfElectronCores
        *process.gsfElectrons
        *process.egammaIsolationSequence
        *process.eIdSequence
        *process.globalMuons
        *process.muons            # redo global muons
        *process.muIsolation      # redo muon isolation
        *process.patDefaultSequence
        )

else:
    process.mainSequence = process.patDefaultSequence


if sampleRequireCollision and prefilter:
    process.p = cms.Path(process.trackSel*process.isoTrackPrefilter*process.mainSequence)
else:
    process.p = cms.Path(process.trackSel*process.mainSequence)

########################################
# PAT CONFIGURATION
########################################

# use only AOD input.
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAodDataTable for what is available in AOD
restrictInputToAOD(process)

# enable PAT trigger functionality
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process)
switchOnTriggerStandAlone(process)
process.patTrigger.processName = sampleHLTProcess

# some PAT objects fail to be produced properly, especially when rerecoing. Thus drop them.
removeSpecificPATObjects(process, ['METs','Jets','Taus'])
# not sure why the following statements are required in addition
process.p.remove(process.patTaus)
process.p.remove(process.patJets)
process.p.remove(process.patMETs)

# eliminate cleaning
removeCleaning(process)


########################################
# OUTPUT
########################################

if prefilter:
    process.TFileService = cms.Service("TFileService", 
         fileName = cms.string('prefilter.root')
)

from PhysicsTools.PatAlgos.patEventContent_cff import *

if rereco:
    process.out.fileName = cms.untracked.string('PATtuple_rereco.root')
else:
    process.out.fileName = cms.untracked.string('PATtuple.root')
process.out.outputCommands += ['keep *_patTriggerEvent_*_PAT']
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patTriggerStandAloneEventContent
process.out.outputCommands += patExtraAodEventContent # primary vertices, genParticles etc
process.out.outputCommands += ['keep PileupSummaryInfos_*_*_*']
process.out.outputCommands += ['drop *_towerMaker_*_*']
# keep superclusters for track momentum correction in etrack channel
process.out.outputCommands += ['keep *_correctedHybridSuperClusters_*_*']
process.out.outputCommands += ['keep *_correctedMulti5x5SuperClustersWithPreshower_*_*']
# keep high quality track selection
process.out.outputCommands += ['keep *_trackSel_*_*']
# keep standalone muons without IP constraint
process.out.outputCommands += ['keep *_standAloneMuons__RECO']
# make sure we store the TrackExtra collection
process.out.outputCommands += ['keep recoTrackExtras_generalTracks__*']
