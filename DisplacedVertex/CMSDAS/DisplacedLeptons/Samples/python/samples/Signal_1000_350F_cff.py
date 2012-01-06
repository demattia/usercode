sampleDataSet = '/HTo2LongLivedTo4F_MH-1000_MFF-350_CTau-350_7TeV-pythia6/Summer11-PU_S4_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000
# cross-section in pb.
# note that for this signal sample, this cross-section is the sum of
# cross-sections for the two different types of long-lived particles
# in this sample. By putting in a cross-section of 10pb, we thus actually
# look at 5pb for the physics we are interested in.
sampleXSec = 2 # this actually means 1 for the lepton channel # pb


# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 6000113

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 106724

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/Signal_1000_350F"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_Hs7.root",
  sampleBaseDir+"/pat/PATtuple_2_1_yiq.root",
  sampleBaseDir+"/pat/PATtuple_3_1_9mQ.root",
  sampleBaseDir+"/pat/PATtuple_4_1_x4F.root",
  sampleBaseDir+"/pat/PATtuple_5_1_Yoo.root",
  sampleBaseDir+"/pat/PATtuple_6_1_Cdl.root",
  sampleBaseDir+"/pat/PATtuple_7_1_RX0.root",
  sampleBaseDir+"/pat/PATtuple_8_1_Rwz.root",
  sampleBaseDir+"/pat/PATtuple_9_1_LAm.root",
  sampleBaseDir+"/pat/PATtuple_10_1_0qt.root",
  sampleBaseDir+"/pat/PATtuple_11_1_D83.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
