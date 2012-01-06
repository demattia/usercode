sampleDataSet = '/ZJetToMuMu_Pt-80to120_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 110000

sampleXSec = 9.99 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZmumuJets80"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_wTh.root",
  sampleBaseDir+"/pat/PATtuple_2_1_ie9.root",
  sampleBaseDir+"/pat/PATtuple_3_1_GPf.root",
  sampleBaseDir+"/pat/PATtuple_4_1_LjS.root",
  sampleBaseDir+"/pat/PATtuple_5_1_Lml.root",
  sampleBaseDir+"/pat/PATtuple_6_1_I7y.root",
  sampleBaseDir+"/pat/PATtuple_7_1_zOG.root",
  sampleBaseDir+"/pat/PATtuple_8_1_Vvt.root",
  sampleBaseDir+"/pat/PATtuple_9_1_PeT.root",
  sampleBaseDir+"/pat/PATtuple_10_1_e9X.root",
  sampleBaseDir+"/pat/PATtuple_11_1_hd2.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunE = 0
