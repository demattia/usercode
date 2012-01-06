sampleDataSet = '/ZJetToEE_Pt-50to80_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 106592

sampleXSec = 32.3 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZeeJets50"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_23R.root",
  sampleBaseDir+"/pat/PATtuple_2_1_bWX.root",
  sampleBaseDir+"/pat/PATtuple_3_1_8Ng.root",
  sampleBaseDir+"/pat/PATtuple_4_1_fxS.root",
  sampleBaseDir+"/pat/PATtuple_5_1_uMi.root",
  sampleBaseDir+"/pat/PATtuple_6_1_Fvg.root",
  sampleBaseDir+"/pat/PATtuple_7_1_XRx.root",
  sampleBaseDir+"/pat/PATtuple_8_1_XUG.root",
  sampleBaseDir+"/pat/PATtuple_9_1_ljk.root",
  sampleBaseDir+"/pat/PATtuple_10_1_g0x.root",
  sampleBaseDir+"/pat/PATtuple_11_1_tw3.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunMu = 0
