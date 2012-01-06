sampleDataSet = '/ZJetToEE_Pt-30to50_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 165000

sampleXSec = 84 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZeeJets30"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_KBr.root",
  sampleBaseDir+"/pat/PATtuple_2_1_ZLI.root",
  sampleBaseDir+"/pat/PATtuple_3_1_qv8.root",
  sampleBaseDir+"/pat/PATtuple_4_1_LLa.root",
  sampleBaseDir+"/pat/PATtuple_5_1_hpZ.root",
  sampleBaseDir+"/pat/PATtuple_6_1_PVm.root",
  sampleBaseDir+"/pat/PATtuple_7_1_Pv0.root",
  sampleBaseDir+"/pat/PATtuple_8_1_rzT.root",
  sampleBaseDir+"/pat/PATtuple_9_1_0i0.root",
  sampleBaseDir+"/pat/PATtuple_10_1_9D1.root",
  sampleBaseDir+"/pat/PATtuple_11_1_4jZ.root",
  sampleBaseDir+"/pat/PATtuple_12_1_Mhz.root",
  sampleBaseDir+"/pat/PATtuple_13_1_mAe.root",
  sampleBaseDir+"/pat/PATtuple_14_1_HgV.root",
  sampleBaseDir+"/pat/PATtuple_15_1_ivH.root",
  sampleBaseDir+"/pat/PATtuple_16_1_4xc.root",
  sampleBaseDir+"/pat/PATtuple_17_1_lpU.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunMu = 0
