# use this only for samples that are created locally with cmsDriver
# but where the GEN-SIM step is taken care of by a different sample config file
sampleGeneratorReference = "python/generator/DisplacedMuPt50_cfi.py"

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

sampleBaseDir = "dcap://dcap.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms/store/user/"\
                +"harder/longlived/"+sampleRelease+"/DisplacedMu_50GeV"

sampleSimFiles = [
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_e2q.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2_1_qLD.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_3_1_aOY.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_G5T.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_OwN.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_6_1_1nz.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_1_AJI.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_8_1_IqC.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_WHN.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_10_1_zV1.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11_1_mIr.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_12_1_iBD.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_13_1_RgC.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_14_1_USL.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_15_1_Plb.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_16_1_RJO.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_17_1_9DU.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_18_1_Nec.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19_1_6sV.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_20_1_MYF.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_21_1_p6U.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_22_1_iRH.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_23_1_4X6.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_24_1_aS8.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_25_1_XMI.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_26_1_sCK.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_27_1_Wcv.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_28_1_oEv.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_29_1_2VZ.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_30_1_YwW.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_31_1_oHf.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_32_1_tXv.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_33_1_yYj.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_34_1_ZIT.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_35_1_U08.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_36_1_Bwm.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37_1_edj.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_38_1_0A3.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_39_1_D8Z.root",
       sampleBaseDir+"/gen/DisplacedMuPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_40_1_69W.root"
     ]

sampleRecoFiles = []
for i in range(40):
    sampleRecoFiles.append(sampleBaseDir\
                          +"/modreco/DisplacedMuPt50_cfi_py_RAW2DIGI_RECO_"\
                          +str(i+1)+".root")

samplePatFiles = [sampleBaseDir+"/pat/PATtuple_modreco_1.root"]
