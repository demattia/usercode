sampleDataSet = '/HTo2LongLivedTo4F_MH-1000_MFF-20_CTau-15_7TeV-pythia6/Summer11-PU_S4_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000
sampleXSec = 2 # this actually means 1 for the lepton channel # pb

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 6000113

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 101753

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/Signal_1000_020F"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_BUe.root",
  sampleBaseDir+"/pat/PATtuple_2_1_f7C.root",
  sampleBaseDir+"/pat/PATtuple_3_1_lfU.root",
  sampleBaseDir+"/pat/PATtuple_4_1_hC5.root",
  sampleBaseDir+"/pat/PATtuple_5_1_yxe.root",
  sampleBaseDir+"/pat/PATtuple_6_1_vKm.root",
  sampleBaseDir+"/pat/PATtuple_7_1_rVt.root",
  sampleBaseDir+"/pat/PATtuple_8_1_Mwo.root",
  sampleBaseDir+"/pat/PATtuple_9_1_bRA.root",
  sampleBaseDir+"/pat/PATtuple_10_1_BbI.root",
  sampleBaseDir+"/pat/PATtuple_11_1_qE5.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
