sampleDataSet = '/HTo2LongLivedTo4F_MH-200_MFF-20_CTau-70_7TeV-pythia6/Summer11-PU_S4_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000
sampleXSec = 2 # this actually means 1 for the lepton channel # pb

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 6000113

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 103171

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/Signal_200_020F"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_vyH.root",
  sampleBaseDir+"/pat/PATtuple_2_1_T2M.root",
  sampleBaseDir+"/pat/PATtuple_3_1_076.root",
  sampleBaseDir+"/pat/PATtuple_4_1_hf1.root",
  sampleBaseDir+"/pat/PATtuple_5_1_wj5.root",
  sampleBaseDir+"/pat/PATtuple_6_1_cXY.root",
  sampleBaseDir+"/pat/PATtuple_7_1_T3K.root",
  sampleBaseDir+"/pat/PATtuple_8_1_GBN.root",
  sampleBaseDir+"/pat/PATtuple_9_1_1cw.root",
  sampleBaseDir+"/pat/PATtuple_10_1_u6Z.root",
  sampleBaseDir+"/pat/PATtuple_11_1_fUJ.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
