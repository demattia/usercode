sampleDataSet = '/ZJetToMuMu_Pt-170to230_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 108506

sampleXSec = 0.722 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZmumuJets170"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_X4w.root",
  sampleBaseDir+"/pat/PATtuple_2_1_2FO.root",
  sampleBaseDir+"/pat/PATtuple_3_1_eIV.root",
  sampleBaseDir+"/pat/PATtuple_4_1_5Vo.root",
  sampleBaseDir+"/pat/PATtuple_5_1_KmH.root",
  sampleBaseDir+"/pat/PATtuple_6_1_JE6.root",
  sampleBaseDir+"/pat/PATtuple_7_1_mXB.root",
  sampleBaseDir+"/pat/PATtuple_8_1_P4E.root",
  sampleBaseDir+"/pat/PATtuple_9_1_Chj.root",
  sampleBaseDir+"/pat/PATtuple_10_1_p2e.root",
  sampleBaseDir+"/pat/PATtuple_11_1_erA.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunE = 0
