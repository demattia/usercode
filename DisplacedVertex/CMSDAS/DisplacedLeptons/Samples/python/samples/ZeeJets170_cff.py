sampleDataSet = '/ZJetToEE_Pt-170to230_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 110000

sampleXSec = 0.722 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZeeJets170"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_1HD.root",
  sampleBaseDir+"/pat/PATtuple_2_1_r5q.root",
  sampleBaseDir+"/pat/PATtuple_3_1_2t1.root",
  sampleBaseDir+"/pat/PATtuple_4_1_pHO.root",
  sampleBaseDir+"/pat/PATtuple_7_1_75o.root",
  sampleBaseDir+"/pat/PATtuple_8_1_pWT.root",
  sampleBaseDir+"/pat/PATtuple_9_1_eFh.root",
  sampleBaseDir+"/pat/PATtuple_10_1_ri8.root",
  sampleBaseDir+"/pat/PATtuple_11_1_Bqu.root",
  sampleBaseDir+"/pat/PATtuple_12_1_Zr6.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
sampleRunMu = 0
