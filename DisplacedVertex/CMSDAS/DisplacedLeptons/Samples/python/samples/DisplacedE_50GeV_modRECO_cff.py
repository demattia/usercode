# use this only for samples that are created locally with cmsDriver
# but where the GEN-SIM step is taken care of by a different sample config file
sampleGeneratorReference = "python/generator/DisplacedEPt50_cfi.py"

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

sampleBaseDir = "dcap://dcap.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms/store/user/"\
                +"harder/longlived/"+sampleRelease+"/DisplacedE_50GeV"

sampleSimFiles = [
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_wWu.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2_1_ZOS.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_3_1_swh.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_e1x.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_hxv.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_6_1_4hG.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_1_ZxU.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_8_1_NYD.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_Iqg.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_10_1_1rz.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11_1_SrT.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_12_1_cgc.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_13_1_Fvo.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_14_1_ynF.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_15_1_hO7.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_16_1_47R.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_17_1_H6P.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_18_1_qLh.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19_1_Acc.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_20_1_RiX.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_21_1_RV4.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_22_1_HFF.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_23_1_pUP.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_24_1_rKu.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_25_1_aJE.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_26_1_WDC.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_27_1_L73.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_28_1_DVY.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_29_1_0Jr.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_30_1_89q.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_31_1_3di.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_32_1_24e.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_33_1_LI3.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_34_1_Cgs.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_35_1_UTo.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_36_1_OTV.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37_1_n6D.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_38_1_JcA.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_39_1_BMN.root",
       sampleBaseDir+"/gen/DisplacedEPt50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_40_1_16T.root"
    ]

sampleRecoFiles = []
for i in range(40):
    sampleRecoFiles.append(sampleBaseDir\
                          +"/modreco/DisplacedEPt50_cfi_py_RAW2DIGI_RECO_"\
                          +str(i+1)+".root")

samplePatFiles = [sampleBaseDir+"/pat/PATtuple_modreco_1.root"]
