sampleDataSet = '/HTo2LongLivedTo4F_MH-400_MFF-150_CTau-400_7TeV-pythia6/Summer11-PU_S4_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000
sampleXSec = 2 # this actually means 1 for the lepton channel # pb

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 6000113

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 102842

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/Signal_400_150F"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_qn5.root",
  sampleBaseDir+"/pat/PATtuple_2_1_oOP.root",
  sampleBaseDir+"/pat/PATtuple_3_1_NFT.root",
  sampleBaseDir+"/pat/PATtuple_4_1_5gH.root",
  sampleBaseDir+"/pat/PATtuple_5_1_uR7.root",
  sampleBaseDir+"/pat/PATtuple_6_1_pJQ.root",
  sampleBaseDir+"/pat/PATtuple_7_1_dW7.root",
  sampleBaseDir+"/pat/PATtuple_8_1_hES.root",
  sampleBaseDir+"/pat/PATtuple_9_1_k5t.root",
  sampleBaseDir+"/pat/PATtuple_10_1_0lW.root",
  sampleBaseDir+"/pat/PATtuple_11_1_umx.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
