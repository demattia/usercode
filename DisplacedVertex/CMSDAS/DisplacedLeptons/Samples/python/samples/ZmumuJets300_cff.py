sampleDataSet = '/ZJetToMuMu_Pt-300_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 110000

sampleXSec = 0.0758 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZmumuJets300"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_VTH.root",
  sampleBaseDir+"/pat/PATtuple_2_1_FTs.root",
  sampleBaseDir+"/pat/PATtuple_3_1_y99.root",
  sampleBaseDir+"/pat/PATtuple_4_1_VrQ.root",
  sampleBaseDir+"/pat/PATtuple_5_1_oL2.root",
  sampleBaseDir+"/pat/PATtuple_6_1_xCc.root",
  sampleBaseDir+"/pat/PATtuple_7_1_DFu.root",
  sampleBaseDir+"/pat/PATtuple_8_1_PwI.root",
  sampleBaseDir+"/pat/PATtuple_9_1_x3A.root",
  sampleBaseDir+"/pat/PATtuple_10_1_INw.root",
  sampleBaseDir+"/pat/PATtuple_11_1_6gC.root",
  sampleBaseDir+"/pat/PATtuple_12_1_coc.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunE = 0
