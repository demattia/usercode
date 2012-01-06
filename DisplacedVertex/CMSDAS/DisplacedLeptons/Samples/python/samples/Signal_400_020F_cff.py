sampleDataSet = '/HTo2LongLivedTo4F_MH-400_MFF-20_CTau-40_7TeV-pythia6/Summer11-PU_S4_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000
sampleXSec = 2 # this actually means 1 for the lepton channel # pb

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 6000113

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 105289

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/Signal_400_020F"

sampleRecoFiles = [ ]

samplePatFiles = [
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
