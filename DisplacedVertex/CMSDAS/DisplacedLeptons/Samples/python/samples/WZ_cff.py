sampleDataSet = '/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/Summer11-PU_S4_START42_V11-v1/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 204725

sampleXSec = 0.3394 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/WZ"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_2_EUZ.root",
  sampleBaseDir+"/pat/PATtuple_2_2_84r.root",
  sampleBaseDir+"/pat/PATtuple_3_1_mRi.root",
  sampleBaseDir+"/pat/PATtuple_4_1_YiQ.root",
  sampleBaseDir+"/pat/PATtuple_5_0_sQH.root",
  sampleBaseDir+"/pat/PATtuple_6_1_J9u.root",
  sampleBaseDir+"/pat/PATtuple_7_0_VHO.root",
  sampleBaseDir+"/pat/PATtuple_8_0_sgG.root",
  sampleBaseDir+"/pat/PATtuple_9_0_U93.root",
  sampleBaseDir+"/pat/PATtuple_10_0_UAl.root",
  sampleBaseDir+"/pat/PATtuple_11_0_ZIz.root",
  sampleBaseDir+"/pat/PATtuple_12_0_eZ2.root",
  sampleBaseDir+"/pat/PATtuple_13_0_wFV.root",
  sampleBaseDir+"/pat/PATtuple_14_0_0HF.root",
  sampleBaseDir+"/pat/PATtuple_15_0_xvN.root",
  sampleBaseDir+"/pat/PATtuple_16_0_6pO.root",
  sampleBaseDir+"/pat/PATtuple_17_0_Y7f.root",
  sampleBaseDir+"/pat/PATtuple_18_0_UBu.root",
  sampleBaseDir+"/pat/PATtuple_19_0_GX8.root",
  sampleBaseDir+"/pat/PATtuple_20_0_ZC3.root",
  sampleBaseDir+"/pat/PATtuple_21_0_SS5.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
