sampleDataSet = '/DoubleMu/Run2011A-PromptReco-v6/AOD'
sampleNumEvents = 6854601 # according to DBS, as of 03 October 2011

# global tag needs to be chosen to match release, trigger menu and alignment conditions.
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
sampleGlobalTag = 'GR_R_42_V21A::All'
sampleHLTProcess = '*'

# data quality run/lumi section selection
sampleJSON = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions11/7TeV/Prompt/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt"

# restrict run range (mainly to get a sample with consistent trigger configuration)
sampleRunRange = [160000,999999]

# use checkAllFilesOpened whenever possible, and noDuplicateCheck when necessary
sampleDuplicateCheckMode = "checkAllFilesOpened"

# "DATA" or "MC"
sampleType = "DATA"

sampleRelease = "CMSSW_4_2_8_patch3" # original (i.e. RECO file) release,
                              # not the one we plan to process them with

sampleProcessRelease = "CMSSW_4_2_8_patch7" # release used to create new files with

sampleBaseDir = "dcap://dcap.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms/store/user/"\
                +"harder/longlived/"+sampleProcessRelease+"/Data_Mu_Run2011A4"

sampleRecoFiles = []

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_Vvf.root",
  sampleBaseDir+"/pat/PATtuple_2_1_lyh.root",
  sampleBaseDir+"/pat/PATtuple_3_0_VxB.root",
  sampleBaseDir+"/pat/PATtuple_4_1_Wb0.root",
  sampleBaseDir+"/pat/PATtuple_5_1_XEb.root",
  sampleBaseDir+"/pat/PATtuple_6_1_ggh.root",
  sampleBaseDir+"/pat/PATtuple_7_1_F1v.root",
  sampleBaseDir+"/pat/PATtuple_8_1_bj6.root",
  sampleBaseDir+"/pat/PATtuple_9_1_YKO.root",
  sampleBaseDir+"/pat/PATtuple_10_1_i3D.root",
  sampleBaseDir+"/pat/PATtuple_11_1_52j.root",
  sampleBaseDir+"/pat/PATtuple_12_1_918.root",
  sampleBaseDir+"/pat/PATtuple_13_1_TNC.root",
  sampleBaseDir+"/pat/PATtuple_14_0_zgD.root",
  sampleBaseDir+"/pat/PATtuple_15_1_fgD.root",
  sampleBaseDir+"/pat/PATtuple_16_1_rcz.root",
  sampleBaseDir+"/pat/PATtuple_17_0_bCE.root",
  sampleBaseDir+"/pat/PATtuple_18_0_BJD.root",
  sampleBaseDir+"/pat/PATtuple_19_0_ZM8.root",
  sampleBaseDir+"/pat/PATtuple_20_0_tn0.root",
  sampleBaseDir+"/pat/PATtuple_21_0_w4U.root",
  sampleBaseDir+"/pat/PATtuple_22_1_A2u.root",
  sampleBaseDir+"/pat/PATtuple_23_1_NSm.root",
  sampleBaseDir+"/pat/PATtuple_24_1_ehf.root",
  sampleBaseDir+"/pat/PATtuple_25_1_zv2.root",
  sampleBaseDir+"/pat/PATtuple_26_1_jik.root",
  sampleBaseDir+"/pat/PATtuple_27_1_6xS.root",
  sampleBaseDir+"/pat/PATtuple_28_1_j1g.root",
  sampleBaseDir+"/pat/PATtuple_29_1_C3I.root",
  sampleBaseDir+"/pat/PATtuple_30_1_MTZ.root",
  sampleBaseDir+"/pat/PATtuple_31_1_ydg.root",
  sampleBaseDir+"/pat/PATtuple_32_1_DKg.root",
  sampleBaseDir+"/pat/PATtuple_33_1_PYe.root",
  sampleBaseDir+"/pat/PATtuple_34_1_f6t.root",
  sampleBaseDir+"/pat/PATtuple_35_1_OHW.root",
  sampleBaseDir+"/pat/PATtuple_36_1_d1R.root",
  sampleBaseDir+"/pat/PATtuple_37_1_PC9.root",
  sampleBaseDir+"/pat/PATtuple_38_1_KuR.root",
  sampleBaseDir+"/pat/PATtuple_39_1_2f5.root",
  sampleBaseDir+"/pat/PATtuple_40_0_toE.root",
  sampleBaseDir+"/pat/PATtuple_41_0_MnG.root",
  sampleBaseDir+"/pat/PATtuple_42_1_PEy.root",
  sampleBaseDir+"/pat/PATtuple_43_1_pNU.root",
  sampleBaseDir+"/pat/PATtuple_44_1_2i2.root",
  sampleBaseDir+"/pat/PATtuple_45_1_EYe.root",
  sampleBaseDir+"/pat/PATtuple_46_1_LMI.root",
  sampleBaseDir+"/pat/PATtuple_47_1_qoy.root",
  sampleBaseDir+"/pat/PATtuple_48_1_xFm.root",
  sampleBaseDir+"/pat/PATtuple_49_1_DXd.root",
  sampleBaseDir+"/pat/PATtuple_50_1_usg.root",
  sampleBaseDir+"/pat/PATtuple_51_1_xcu.root",
  sampleBaseDir+"/pat/PATtuple_52_1_Sfb.root",
  sampleBaseDir+"/pat/PATtuple_53_1_etm.root",
  sampleBaseDir+"/pat/PATtuple_54_1_2Fq.root",
  sampleBaseDir+"/pat/PATtuple_55_1_Uma.root",
  sampleBaseDir+"/pat/PATtuple_56_1_Y7W.root",
  sampleBaseDir+"/pat/PATtuple_57_1_lLa.root",
  sampleBaseDir+"/pat/PATtuple_58_1_JQ7.root",
  sampleBaseDir+"/pat/PATtuple_59_1_gqo.root",
  sampleBaseDir+"/pat/PATtuple_60_1_5Ez.root",
  sampleBaseDir+"/pat/PATtuple_61_1_LIl.root",
  sampleBaseDir+"/pat/PATtuple_62_1_rE3.root",
  sampleBaseDir+"/pat/PATtuple_63_1_OEf.root",
  sampleBaseDir+"/pat/PATtuple_64_1_pBu.root",
  sampleBaseDir+"/pat/PATtuple_65_1_8Ds.root",
  sampleBaseDir+"/pat/PATtuple_66_1_0PS.root",
  sampleBaseDir+"/pat/PATtuple_67_1_JWp.root",
  sampleBaseDir+"/pat/PATtuple_68_1_SER.root",
  sampleBaseDir+"/pat/PATtuple_69_1_vFI.root",
  sampleBaseDir+"/pat/PATtuple_70_0_hlv.root",
  sampleBaseDir+"/pat/PATtuple_71_1_nFu.root",
  sampleBaseDir+"/pat/PATtuple_72_1_Nvv.root",
  sampleBaseDir+"/pat/PATtuple_73_1_6WW.root",
  sampleBaseDir+"/pat/PATtuple_74_1_vli.root",
  sampleBaseDir+"/pat/PATtuple_75_1_6KI.root",
  sampleBaseDir+"/pat/PATtuple_76_1_bpq.root",
  sampleBaseDir+"/pat/PATtuple_77_1_ZRX.root",
  sampleBaseDir+"/pat/PATtuple_78_1_1XF.root",
  sampleBaseDir+"/pat/PATtuple_79_1_sGN.root",
  sampleBaseDir+"/pat/PATtuple_80_1_o5w.root",
  sampleBaseDir+"/pat/PATtuple_81_1_Hug.root",
  sampleBaseDir+"/pat/PATtuple_82_1_dBO.root",
  sampleBaseDir+"/pat/PATtuple_83_1_nlZ.root",
  sampleBaseDir+"/pat/PATtuple_84_0_B2I.root",
  sampleBaseDir+"/pat/PATtuple_85_1_mVA.root",
  sampleBaseDir+"/pat/PATtuple_86_1_p8U.root",
  sampleBaseDir+"/pat/PATtuple_87_1_bnF.root",
  sampleBaseDir+"/pat/PATtuple_88_1_7L5.root",
  sampleBaseDir+"/pat/PATtuple_89_0_UxV.root",
  sampleBaseDir+"/pat/PATtuple_90_0_w7W.root",
  sampleBaseDir+"/pat/PATtuple_91_1_rV5.root",
  sampleBaseDir+"/pat/PATtuple_92_1_uvz.root",
  sampleBaseDir+"/pat/PATtuple_93_1_KLy.root",
  sampleBaseDir+"/pat/PATtuple_94_1_EQ5.root",
  sampleBaseDir+"/pat/PATtuple_95_1_hLd.root",
  sampleBaseDir+"/pat/PATtuple_96_1_jNO.root",
  sampleBaseDir+"/pat/PATtuple_97_1_No7.root",
  sampleBaseDir+"/pat/PATtuple_98_1_INM.root",
  sampleBaseDir+"/pat/PATtuple_99_1_EEw.root",
  sampleBaseDir+"/pat/PATtuple_100_1_QfY.root",
  sampleBaseDir+"/pat/PATtuple_101_1_Uw3.root",
  sampleBaseDir+"/pat/PATtuple_102_1_yWC.root",
  sampleBaseDir+"/pat/PATtuple_103_1_uoO.root",
  sampleBaseDir+"/pat/PATtuple_104_1_jci.root",
  sampleBaseDir+"/pat/PATtuple_105_1_TyR.root",
  sampleBaseDir+"/pat/PATtuple_106_1_4Q4.root",
  sampleBaseDir+"/pat/PATtuple_107_1_jpr.root",
  sampleBaseDir+"/pat/PATtuple_108_1_VBM.root",
  sampleBaseDir+"/pat/PATtuple_109_1_9YF.root",
  sampleBaseDir+"/pat/PATtuple_110_1_oZW.root",
  sampleBaseDir+"/pat/PATtuple_111_1_HEe.root",
  sampleBaseDir+"/pat/PATtuple_112_1_jV0.root",
  sampleBaseDir+"/pat/PATtuple_113_1_klS.root"
    ]
