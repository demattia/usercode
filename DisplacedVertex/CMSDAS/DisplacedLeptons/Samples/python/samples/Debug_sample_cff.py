# use this only for samples that are created locally with cmsDriver
sampleGeneratorConfig = "python/generator/Pythia_H0_pyupda_cfi.py"

# global tag needs to be chosen to match release, trigger menu and alignment conditions.
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
sampleGlobalTag = 'START42_V13::All' 
sampleHLTProcess = "*"

# use checkAllFilesOpened whenever possible, and noDuplicateCheck when necessary
sampleDuplicateCheckMode = "checkAllFilesOpened"

# "DATA" or "MC"
sampleType = "MC"

# centre-of-mass energy this sample is for
sampleCMSEnergy = 7000

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 36

# total number of events in sample
sampleNumEvents = 10000

# lumi section numbers in this sample are between offset and offset plus number of gen jobs
sampleLumiNumberOffset = 120009

# release used to generate (if applicable) and reconstruct this sample
sampleRelease = "CMSSW_4_2_7"

sampleBaseDir = "dcap://dcap.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms/store/user/"\
                +"harder/longlived/"+sampleRelease+"/Pythia_HtoAA_150GeV_100mm_emutau"

sampleSimFiles = [
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_Xqg.root"
  ]

sampleRecoFiles = [
    sampleBaseDir+"/stdreco/Pythia_H0_pyupda_cfi_py_RAW2DIGI_RECO_1.root"
    ]

samplePatFiles = [ "file:/opt/ppd/scratch/harder/PATtuple_stdreco.root" ]
