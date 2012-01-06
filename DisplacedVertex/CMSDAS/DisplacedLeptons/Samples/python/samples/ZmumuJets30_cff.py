sampleDataSet = '/ZJetToMuMu_Pt-30to50_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 155752

sampleXSec = 84 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZmumuJets30"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_Hp3.root",
  sampleBaseDir+"/pat/PATtuple_2_1_mTT.root",
  sampleBaseDir+"/pat/PATtuple_3_1_ZBX.root",
  sampleBaseDir+"/pat/PATtuple_4_1_sph.root",
  sampleBaseDir+"/pat/PATtuple_5_1_Dki.root",
  sampleBaseDir+"/pat/PATtuple_6_1_P91.root",
  sampleBaseDir+"/pat/PATtuple_7_1_xPu.root",
  sampleBaseDir+"/pat/PATtuple_8_1_Bf1.root",
  sampleBaseDir+"/pat/PATtuple_9_1_Uch.root",
  sampleBaseDir+"/pat/PATtuple_10_1_6Mc.root",
  sampleBaseDir+"/pat/PATtuple_11_1_ipn.root",
  sampleBaseDir+"/pat/PATtuple_12_1_CCd.root",
  sampleBaseDir+"/pat/PATtuple_13_1_0Uh.root",
  sampleBaseDir+"/pat/PATtuple_14_1_YNr.root",
  sampleBaseDir+"/pat/PATtuple_15_1_kWQ.root",
  sampleBaseDir+"/pat/PATtuple_16_1_8RR.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunE = 0
