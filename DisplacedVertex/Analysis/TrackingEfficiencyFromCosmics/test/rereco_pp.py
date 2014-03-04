# Auto generated configuration file
# using: 
# Revision: 1.303 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco -s RAW2DIGI,L1Reco,RECO,DQM --data --magField AutoFromDBCurrent --scenario pp --datatier RECO,AOD,DQM --eventcontent RECO,AOD,DQM --customise Configuration/GlobalRuns/reco_TLR_42X.customisePPData --no_exec --python_filename=rereco_pp.py --conditions GR_R_42_V14::All
import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring('reco_DIGI2RAW.root')
    # fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-RECO/START42_V11-v1/0007/CE4277E3-806D-E011-BFD9-00248C55CC97.root'),
    # fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/0ED6D02F-5575-E011-A76D-003048D2C01E.root'),
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/005F3FAC-5875-E011-863B-003048F118DE.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/0ED6D02F-5575-E011-A76D-003048D2C01E.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/127488C8-4E75-E011-B348-003048F024F6.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/1AA42BAC-5875-E011-9C6D-003048F024E0.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/28FF0014-4E75-E011-8C2D-0030487A3232.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/30552A60-4D75-E011-9BFA-0030487C6088.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/34137EE7-5075-E011-9632-0030487CF41E.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/343EC67F-5475-E011-B3CF-003048F118E0.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/34B445F2-4B75-E011-85CF-0030487C2B86.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/52DDA6DD-4975-E011-962E-0030487CF41E.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/7ABD903D-5C75-E011-901B-0030487C60AE.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/7C40F6DD-4975-E011-9D15-0030487CD184.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/821FC758-4675-E011-8D37-003048F118D2.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/8291D760-5275-E011-B812-0030487A3C92.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/8C7C995E-4D75-E011-9786-003048D2BBF0.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/9C26FDB4-4775-E011-BBC6-003048F1BF66.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/A629A5E7-5075-E011-8AD6-0030487CD700.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/AEE76A8F-4A75-E011-9D72-0030487CD704.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/BE1A9459-4675-E011-AF95-0030487A3C9A.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/C64F0C59-4675-E011-8AD7-0030487A1884.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/C8D4C249-6375-E011-ADB5-00304879FA4A.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/CCFD7A9E-5175-E011-B75C-003048D2C092.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/D4D96F12-5A75-E011-A272-0030486733B4.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/E281527F-5475-E011-BDAC-003048F1C420.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/E2875F87-4875-E011-ADB5-003048D37514.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/E4DF60E7-5575-E011-A510-003048F1BF68.root',
                                      'rfio:/castor/cern.ch/cms/store/data/Run2011A/Cosmics/RAW/v1/000/163/853/EC8B74A9-5875-E011-870E-003048F118D2.root'
                                      )
                            )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('reco_RAW2DIGI_L1Reco_RECO_DQM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RECO')
    )
)

process.AODoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.AODEventContent.outputCommands,
    fileName = cms.untracked.string('reco_RAW2DIGI_L1Reco_RECO_DQM_inAOD.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('AOD')
    )
)

process.DQMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.DQMEventContent.outputCommands,
    fileName = cms.untracked.string('reco_RAW2DIGI_L1Reco_RECO_DQM_inDQM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DQM')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'GR_R_42_V14::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)
process.AODoutput_step = cms.EndPath(process.AODoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.dqmoffline_step,process.endjob_step,process.RECOoutput_step,process.AODoutput_step,process.DQMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.GlobalRuns.reco_TLR_42X
from Configuration.GlobalRuns.reco_TLR_42X import customisePPData 

#call to customisation function customisePPData imported from Configuration.GlobalRuns.reco_TLR_42X
process = customisePPData(process)

# End of customisation functions
