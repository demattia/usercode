sampleDataSet = '/ZJetToEE_Pt-80to120_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 110000

sampleXSec = 9.99 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZeeJets80"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_XLw.root",
  sampleBaseDir+"/pat/PATtuple_2_1_O3t.root",
  sampleBaseDir+"/pat/PATtuple_3_1_E1K.root",
  sampleBaseDir+"/pat/PATtuple_4_1_p6c.root",
  sampleBaseDir+"/pat/PATtuple_5_1_cEg.root",
  sampleBaseDir+"/pat/PATtuple_6_1_E0t.root",
  sampleBaseDir+"/pat/PATtuple_7_1_bRN.root",
  sampleBaseDir+"/pat/PATtuple_8_1_zY6.root",
  sampleBaseDir+"/pat/PATtuple_9_1_eD9.root",
  sampleBaseDir+"/pat/PATtuple_10_1_PyJ.root",
  sampleBaseDir+"/pat/PATtuple_11_1_XAl.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunMu = 0
