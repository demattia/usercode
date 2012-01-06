sampleDataSet = '/HTo2LongLivedTo4F_MH-1000_MFF-150_CTau-100_7TeV-pythia6/Summer11-PU_S4_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000
sampleXSec = 2 # this actually means 1 for the lepton channel # pb

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 6000113

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 102636

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/Signal_1000_150F"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_2_2Fr.root",
  sampleBaseDir+"/pat/PATtuple_2_2_rgO.root",
  sampleBaseDir+"/pat/PATtuple_3_3_j3r.root",
  sampleBaseDir+"/pat/PATtuple_4_3_Z3T.root",
  sampleBaseDir+"/pat/PATtuple_5_3_VG8.root",
  sampleBaseDir+"/pat/PATtuple_7_2_s0w.root",
  sampleBaseDir+"/pat/PATtuple_8_2_3Vq.root",
  sampleBaseDir+"/pat/PATtuple_9_3_vfT.root",
  sampleBaseDir+"/pat/PATtuple_10_2_8xw.root",
  sampleBaseDir+"/pat/PATtuple_11_3_HCq.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
