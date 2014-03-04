# Auto generated configuration file
# using: 
# Revision: 1.303 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco -s RAW2DIGI,L1Reco,RECO,DQM --mc --magField AutoFromDBCurrent --scenario pp --datatier RECO,AOD,DQM --eventcontent RECO,AOD,DQM --customise Configuration/GlobalRuns/reco_TLR_42X.customisePPMC --no_exec --python_filename=rereco_sim_pp.py --conditions START42_V11::All
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//047F2F11-836D-E011-84DC-00248C55CC9D.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//08583CEB-826D-E011-B16B-002618FDA259.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//08683FF1-826D-E011-A84E-002618943836.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//0C31712A-826D-E011-BAF5-003048678D9A.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//0CC996EE-826D-E011-BAFF-002354EF3BDC.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//0E36AD04-836D-E011-A82A-003048678B18.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//1217CC21-856D-E011-92F0-003048678F1C.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//123F4E07-826D-E011-AAB9-002618943896.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//125F1060-836D-E011-AB16-003048678BE8.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//14ECE803-826D-E011-9901-002618943937.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//1860811B-826D-E011-A737-00261894389D.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//18B84D01-826D-E011-88B5-00304867D446.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//343CF300-836D-E011-B25B-00261894396B.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//3AD7DDE9-816D-E011-AD71-0026189438D7.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//3ECC7701-826D-E011-981F-002618FDA248.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//46F08DFA-826D-E011-BB16-002618943854.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//4A31DDF7-826D-E011-A5FA-001A92810AB6.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//4C7E35EC-826D-E011-9EA8-0026189438B9.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//4EBE2DFE-816D-E011-ABCA-001A92971B9C.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//58DA82F1-816D-E011-876B-003048678AC8.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//5C0C48EB-826D-E011-A105-0026189438C9.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//5C0C7009-836D-E011-9B4A-002354EF3BE0.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//645741E1-816D-E011-95AE-002618943937.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//664B2FF8-816D-E011-B4B2-0018F3D09648.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//68C9D107-826D-E011-919D-002618943921.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//74140716-836D-E011-A72E-001A92971B9C.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//786A6904-836D-E011-94D3-0026189438F2.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//7E2EA4FF-826D-E011-8270-002618943963.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//82A4F71E-826D-E011-ADE1-0030486791AA.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//84B4ECCF-836D-E011-9F8D-003048678F8E.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//8A57A76C-836D-E011-B2B3-002618943856.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//8CA18308-836D-E011-9AC3-00261894390C.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//9219F414-836D-E011-ABFF-003048678A6C.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//929CB51B-836D-E011-A13A-003048678FB8.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//96CF1B1B-826D-E011-B8B3-002618943958.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//96E4600A-846D-E011-9E37-00261894391D.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//98113A08-836D-E011-9C66-002618943836.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//982B5073-836D-E011-B380-0026189438D5.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//9A46D919-826D-E011-9D1C-0026189438AC.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//9AB5CBF8-826D-E011-A021-001A92971B04.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//A4EC6699-836D-E011-8336-001A92810AE6.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//AA3AAE04-826D-E011-81E7-00248C0BE013.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//AE797315-826D-E011-B74C-0026189437FE.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//AE8296F2-826D-E011-9081-002618943986.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//B20B03EA-826D-E011-9F09-0026189437EB.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//B64230EB-826D-E011-9BB8-0026189438F2.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//B67498FC-816D-E011-9C93-002618943866.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//BC1E7BD7-826D-E011-9F3A-002618943958.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//BCE21AF2-816D-E011-A66F-003048678A78.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//C404786E-836D-E011-8933-00261894395A.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//C8E5110D-846D-E011-BA9F-002618FDA237.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//C8F16FE2-856D-E011-9D01-00304867BFAE.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//CAF55707-826D-E011-B76C-003048678A76.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//D43E1BF5-826D-E011-8ACA-0026189437FA.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//DCACEA05-826D-E011-9697-0026189437F0.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//E4B340E7-816D-E011-8187-002618943956.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//E659C604-826D-E011-8447-002618943958.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//EEC87A55-836D-E011-AB1F-003048678BB2.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//EEEC446C-836D-E011-9A7F-0026189438EF.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//F263601D-826D-E011-BFB0-003048678FC4.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//F2B87319-826D-E011-B039-002618943831.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//F86A6308-826D-E011-878B-003048D42D92.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//FA8B3E04-836D-E011-B6D1-0026189438F9.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//FC19610A-846D-E011-8ED1-00261894391D.root',
                                      'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_2/RelValCosmics/GEN-SIM-DIGI-RAW-HLTDEBUG/START42_V11-v1/0008//FC632705-836D-E011-8E8B-0018F3D09702.root'
                                      )
                                      # 'file:047F2F11-836D-E011-84DC-00248C55CC9D.root')
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

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('SIM')
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
process.GlobalTag.globaltag = 'START42_V11::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.AODoutput_step = cms.EndPath(process.AODoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
# process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.dqmoffline_step,process.endjob_step,process.RECOoutput_step,process.AODoutput_step,process.DQMoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.dqmoffline_step,process.endjob_step,process.RECOSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.GlobalRuns.reco_TLR_42X
from Configuration.GlobalRuns.reco_TLR_42X import customisePPMC 

#call to customisation function customisePPMC imported from Configuration.GlobalRuns.reco_TLR_42X
process = customisePPMC(process)

# End of customisation functions
