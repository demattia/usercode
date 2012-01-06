sampleDataSet = '/ZJetToMuMu_Pt-230to300_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 108360

sampleXSec = 0.194 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZmumuJets230"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_iJK.root",
  sampleBaseDir+"/pat/PATtuple_2_1_ppi.root",
  sampleBaseDir+"/pat/PATtuple_3_1_qWD.root",
  sampleBaseDir+"/pat/PATtuple_4_1_eD1.root",
  sampleBaseDir+"/pat/PATtuple_5_1_mZ8.root",
  sampleBaseDir+"/pat/PATtuple_6_1_0ur.root",
  sampleBaseDir+"/pat/PATtuple_7_1_Mpo.root",
  sampleBaseDir+"/pat/PATtuple_8_1_G7J.root",
  sampleBaseDir+"/pat/PATtuple_9_1_iO2.root",
  sampleBaseDir+"/pat/PATtuple_10_1_W4O.root",
  sampleBaseDir+"/pat/PATtuple_11_1_HDY.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunE = 0
