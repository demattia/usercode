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
sampleXSec = 10 # pb

# for signal MC: what is the abs(PDGid) of the long-lived signal particle?
sampleSignalPID = 36

# total number of events in sample
sampleNumEvents = 10000

# lumi section numbers in this sample are between offset and offset plus number of gen jobs
sampleLumiNumberOffset = 12000

# release used to generate (if applicable) and reconstruct this sample
sampleRelease = "CMSSW_4_2_7"

sampleBaseDir = "dcap://dcap.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms/store/user/"\
                +"harder/longlived/"+sampleRelease+"/Pythia_HtoAA_150GeV_100mm_emutau"

sampleSimFiles = [
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_UQk.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2_1_uhL.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_3_1_neG.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_4sX.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_JfU.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_6_1_f0M.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_1_2nE.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_8_1_0b9.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_N31.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_10_1_iFF.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11_1_N3U.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_12_1_xqL.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_13_1_YLz.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_14_1_gGz.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_15_1_4Bq.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_16_1_D7B.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_17_1_pBH.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_18_1_7nW.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19_1_x6B.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_20_1_ieu.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_21_1_oFp.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_22_1_PkP.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_23_1_6ae.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_24_1_Qux.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_25_1_xOj.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_26_1_aWm.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_27_1_2Aj.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_28_1_UZN.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_29_1_uQr.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_30_1_fbS.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_31_1_6BO.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_32_1_vwT.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_33_1_0y3.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_34_1_cyb.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_35_1_rUC.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_36_1_114.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37_1_cYa.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_38_1_yGX.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_39_1_Lqz.root",
  sampleBaseDir+"/gen/Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_40_1_ee9.root"
    ]

sampleRecoFiles = []
for i in range(40):
    sampleRecoFiles.append(sampleBaseDir\
                          +"/stdreco/Pythia_H0_pyupda_cfi_py_RAW2DIGI_RECO_"\
                          +str(i+1)+".root")

samplePatFiles = [ sampleBaseDir+"/pat/PATtuple_rereco.root" ]
