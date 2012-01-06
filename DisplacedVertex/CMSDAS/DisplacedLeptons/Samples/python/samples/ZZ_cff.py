sampleDataSet = '/ZZTo2L2Nu_TuneZ2_7TeV_pythia6_tauola/Summer11-PU_S4_START42_V11-v1/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_3_patch3" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 220000

sampleXSec = 0.1927 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/ZZ"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_Kl5.root",
  sampleBaseDir+"/pat/PATtuple_2_1_dWQ.root",
  sampleBaseDir+"/pat/PATtuple_3_0_VmM.root",
  sampleBaseDir+"/pat/PATtuple_4_0_i9p.root",
  sampleBaseDir+"/pat/PATtuple_5_0_ox4.root",
  sampleBaseDir+"/pat/PATtuple_6_0_GhV.root",
  sampleBaseDir+"/pat/PATtuple_7_0_a2P.root",
  sampleBaseDir+"/pat/PATtuple_8_0_Kw5.root",
  sampleBaseDir+"/pat/PATtuple_9_0_tWR.root",
  sampleBaseDir+"/pat/PATtuple_10_0_AIz.root",
  sampleBaseDir+"/pat/PATtuple_11_0_zBN.root",
  sampleBaseDir+"/pat/PATtuple_12_0_AJ6.root",
  sampleBaseDir+"/pat/PATtuple_13_0_9NS.root",
  sampleBaseDir+"/pat/PATtuple_14_0_OXI.root",
  sampleBaseDir+"/pat/PATtuple_15_0_NTk.root",
  sampleBaseDir+"/pat/PATtuple_16_0_T8H.root",
  sampleBaseDir+"/pat/PATtuple_17_0_Zsd.root",
  sampleBaseDir+"/pat/PATtuple_18_0_CoK.root",
  sampleBaseDir+"/pat/PATtuple_19_0_0Tj.root",
  sampleBaseDir+"/pat/PATtuple_20_0_fUo.root",
  sampleBaseDir+"/pat/PATtuple_21_1_ITg.root",
  sampleBaseDir+"/pat/PATtuple_22_1_ZtC.root",
  sampleBaseDir+"/pat/PATtuple_23_1_N6X.root",
  sampleBaseDir+"/pat/PATtuple_24_1_btA.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
