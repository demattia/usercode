# use this only for samples that are created locally with cmsDriver
sampleGeneratorConfig = "python/generator/DisplacedMuPt50_cfi.py"

# global tag needs to be chosen to match release, trigger menu and alignment conditions.
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
sampleGlobalTag = 'START42_V13::All' 
sampleHLTProcess = "*"

# use checkAllFilesOpened whenever possible, and noDuplicateCheck when necessary
sampleDuplicateCheckMode = "checkAllFilesOpened"

# "DATA" or "MC"
sampleType = "MC"
sampleRequireCollision=0
sampleRunE=0

# total number of events in sample
sampleNumEvents = 10000

# lumi section numbers in this sample are between offset and offset plus number of gen jobs
sampleLumiNumberOffset = 2000

# release used to generate (if applicable) and reconstruct this sample
sampleRelease = "CMSSW_4_2_7"

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleRelease+"/DisplacedMu_50GeV"

sampleSimFiles = [
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_18_1_GmF.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_8_1_hg6.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_21_1_Yat.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_30_1_T2m.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_28_1_EuG.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_kpU.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_alE.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_25_1_Z2N.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_40_1_02r.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_22_1_3BK.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_24_1_lWE.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_g5s.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_10_1_fUa.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_1_NcY.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_31_1_0Vi.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_26_1_lm7.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_15_1_G1X.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_29_1_K71.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_34_1_OfK.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_Lhx.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_27_1_fez.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_36_1_Szp.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_32_1_6bK.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_23_1_B3K.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37_1_0Mk.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_38_1_58T.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_35_1_lwf.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_39_1_Q2E.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_33_1_O0p.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_6_1_S3y.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_20_1_uFX.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19_1_85n.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_12_1_6HP.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_3_1_lVc.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11_1_2Oh.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2_1_gNO.root", 
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_17_1_dXt.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_13_1_crO.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_16_1_11j.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_14_1_MIk.root"
    ]

sampleRecoFiles = []
for i in range(40):
    sampleRecoFiles.append(sampleBaseDir\
                          +"/stdreco/DisplacedMuPt50_cfi_py_RAW2DIGI_RECO_"\
                          +str(i+1)+".root")

samplePatFiles = [sampleBaseDir+"/pat/PATtuple_stdreco_1.root"]
