sampleDataSet = '/ZJetToEE_Pt-120to170_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 107304

sampleXSec = 2.74 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZeeJets120"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_z43.root",
  sampleBaseDir+"/pat/PATtuple_2_1_qfN.root",
  sampleBaseDir+"/pat/PATtuple_3_1_CUC.root",
  sampleBaseDir+"/pat/PATtuple_4_1_uNE.root",
  sampleBaseDir+"/pat/PATtuple_5_1_mJs.root",
  sampleBaseDir+"/pat/PATtuple_6_1_IBw.root",
  sampleBaseDir+"/pat/PATtuple_7_1_oYF.root",
  sampleBaseDir+"/pat/PATtuple_8_1_lzT.root",
  sampleBaseDir+"/pat/PATtuple_9_1_G9B.root",
  sampleBaseDir+"/pat/PATtuple_10_1_x1h.root",
  sampleBaseDir+"/pat/PATtuple_11_1_mio.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunMu = 0
