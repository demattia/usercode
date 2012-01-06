sampleDataSet = '/ZJetToMuMu_Pt-20to30_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 165000

sampleXSec = 131 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZmumuJets20"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_iv7.root",
  sampleBaseDir+"/pat/PATtuple_2_1_f2q.root",
  sampleBaseDir+"/pat/PATtuple_3_1_Rqu.root",
  sampleBaseDir+"/pat/PATtuple_4_1_7lL.root",
  sampleBaseDir+"/pat/PATtuple_5_1_cFK.root",
  sampleBaseDir+"/pat/PATtuple_6_1_0MT.root",
  sampleBaseDir+"/pat/PATtuple_7_1_yCL.root",
  sampleBaseDir+"/pat/PATtuple_8_1_RUo.root",
  sampleBaseDir+"/pat/PATtuple_9_1_nIT.root",
  sampleBaseDir+"/pat/PATtuple_10_1_zs7.root",
  sampleBaseDir+"/pat/PATtuple_11_1_mor.root",
  sampleBaseDir+"/pat/PATtuple_12_1_NQW.root",
  sampleBaseDir+"/pat/PATtuple_13_1_lhL.root",
  sampleBaseDir+"/pat/PATtuple_14_1_gmF.root",
  sampleBaseDir+"/pat/PATtuple_15_1_q0L.root",
  sampleBaseDir+"/pat/PATtuple_16_1_2Rb.root",
  sampleBaseDir+"/pat/PATtuple_17_1_582.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunE = 0
