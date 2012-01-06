sampleDataSet = '/ZJetToEE_Pt-20to30_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 165000

sampleXSec = 131 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZeeJets20"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_py7.root",
  sampleBaseDir+"/pat/PATtuple_2_1_9dB.root",
  sampleBaseDir+"/pat/PATtuple_3_1_5t2.root",
  sampleBaseDir+"/pat/PATtuple_4_1_Shw.root",
  sampleBaseDir+"/pat/PATtuple_5_1_mVk.root",
  sampleBaseDir+"/pat/PATtuple_6_1_Evz.root",
  sampleBaseDir+"/pat/PATtuple_7_1_W1P.root",
  sampleBaseDir+"/pat/PATtuple_8_1_raT.root",
  sampleBaseDir+"/pat/PATtuple_9_1_RHk.root",
  sampleBaseDir+"/pat/PATtuple_10_1_HVK.root",
  sampleBaseDir+"/pat/PATtuple_11_1_RvK.root",
  sampleBaseDir+"/pat/PATtuple_12_1_QTx.root",
  sampleBaseDir+"/pat/PATtuple_13_1_Euc.root",
  sampleBaseDir+"/pat/PATtuple_14_1_mRW.root",
  sampleBaseDir+"/pat/PATtuple_15_1_Xs8.root",
  sampleBaseDir+"/pat/PATtuple_16_1_ZqE.root",
  sampleBaseDir+"/pat/PATtuple_17_1_RMG.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunMu = 0
