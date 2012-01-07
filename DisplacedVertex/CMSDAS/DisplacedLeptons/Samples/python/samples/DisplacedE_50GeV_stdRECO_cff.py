# use this only for samples that are created locally with cmsDriver
sampleGeneratorConfig = "python/generator/DisplacedEPt50_cfi.py"

# global tag needs to be chosen to match release, trigger menu and alignment conditions.
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
sampleGlobalTag = 'START42_V13::All' 
sampleHLTProcess = "*"

# use checkAllFilesOpened whenever possible, and noDuplicateCheck when necessary
sampleDuplicateCheckMode = "checkAllFilesOpened"

# "DATA" or "MC"
sampleType = "MC"
sampleRequireCollision=0
sampleRunMu=0

# total number of events in sample
sampleNumEvents = 10000

# lumi section numbers in this sample are between offset and offset plus number of gen jobs
sampleLumiNumberOffset = 1000

# release used to generate (if applicable) and reconstruct this sample
sampleRelease = "CMSSW_4_2_7"

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleRelease+"/DisplacedE_50GeV"

sampleSimFiles = [
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_23_1_Thx.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_20_1_HJ8.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_12_1_3Uc.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_15_1_v9k.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_Erd.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2_1_R0n.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_28_1_Lpv.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_33_1_HaS.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_31_1_aET.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_36_1_IHc.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_BxB.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_39_1_1C8.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_35_1_5VY.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_18_1_Mlp.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_34_1_J2e.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_JGd.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_30_1_0mW.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_25_1_LDF.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19_1_O6P.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_13_1_imP.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_40_1_Gfv.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_21_1_O9j.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37_1_ICz.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_8_1_Jvj.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_27_1_mpi.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_38_1_nbb.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_29_1_dUI.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_1_jfL.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_16_1_bOz.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11_1_LYm.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_10_1_nTD.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_U8w.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_26_1_fyr.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_3_1_1ut.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_17_1_VrR.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_14_1_gjB.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_32_1_WEa.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_24_1_iyi.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_22_1_zch.root",
    sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_6_1_sHV.root"
    ]

sampleRecoFiles = []
for i in range(40):
    sampleRecoFiles.append(sampleBaseDir\
                          +"/stdreco/DisplacedEPt50_cfi_py_RAW2DIGI_RECO_"\
                          +str(i+1)+".root")

samplePatFiles = [sampleBaseDir+"/pat/PATtuple_stdreco_1.root"]
