import FWCore.ParameterSet.Config as cms

process = cms.Process("CALIB")
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring('cout')
)

process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(166000),
    lastValue = cms.uint64(166000),
    interval = cms.uint64(1)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

sistripconn = cms.ESProducer("SiStripConnectivity")

#-------------------------------------------------
# Calibration
#-------------------------------------------------
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.GlobalTag.globaltag = "GR_R_44_V4::All"

process.poolDBESSource = cms.ESSource("PoolDBESSource",
                                      BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
                                      DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
    ),
                                      timetype = cms.untracked.string('runnumber'),
                                      connect = cms.string('oracle://cms_orcoff_prod/CMS_COND_31X_STRIP'),
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('SiStripFedCablingRcd'),
    tag = cms.string('SiStripFedCabling_GR10_v1_hlt')
    ),
    )
)

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
