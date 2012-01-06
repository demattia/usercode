sampleDataSet = '/DoubleMu/Run2011A-05Aug2011-v1/AOD'
sampleNumEvents = 5824586 # according to DBS, as of 24 Nov 2011

# global tag needs to be chosen to match release, trigger menu and alignment conditions.
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
sampleGlobalTag = 'GR_R_42_V21A::All'
sampleHLTProcess = '*'

# data quality run/lumi section selection
sampleJSON = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions11/7TeV/Reprocessing/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_MuonPhys_v2.txt"

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
                +"harder/longlived/"+sampleProcessRelease+"/Data_Mu_Run2011A3fix"

sampleRecoFiles = []

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_NBY.root",
  sampleBaseDir+"/pat/PATtuple_2_1_xud.root",
  sampleBaseDir+"/pat/PATtuple_3_3_B3t.root",
  sampleBaseDir+"/pat/PATtuple_4_2_tj9.root",
  sampleBaseDir+"/pat/PATtuple_5_2_J3M.root",
  sampleBaseDir+"/pat/PATtuple_6_2_Eci.root",
  sampleBaseDir+"/pat/PATtuple_7_2_GqG.root",
  sampleBaseDir+"/pat/PATtuple_8_1_kY8.root",
  sampleBaseDir+"/pat/PATtuple_9_1_19w.root",
  sampleBaseDir+"/pat/PATtuple_10_1_KXy.root",
  sampleBaseDir+"/pat/PATtuple_11_1_kEq.root",
  sampleBaseDir+"/pat/PATtuple_12_1_Lzg.root",
  sampleBaseDir+"/pat/PATtuple_13_1_cLe.root",
  sampleBaseDir+"/pat/PATtuple_14_1_cTb.root",
  sampleBaseDir+"/pat/PATtuple_15_1_wOC.root",
  sampleBaseDir+"/pat/PATtuple_16_1_kqd.root",
  sampleBaseDir+"/pat/PATtuple_17_1_3Yw.root",
  sampleBaseDir+"/pat/PATtuple_18_1_9Us.root",
  sampleBaseDir+"/pat/PATtuple_19_1_r4E.root",
  sampleBaseDir+"/pat/PATtuple_20_1_QnE.root",
  sampleBaseDir+"/pat/PATtuple_21_1_yu0.root",
  sampleBaseDir+"/pat/PATtuple_22_1_cK8.root",
  sampleBaseDir+"/pat/PATtuple_23_1_yzt.root",
  sampleBaseDir+"/pat/PATtuple_24_1_XC6.root",
  sampleBaseDir+"/pat/PATtuple_25_1_Gso.root",
  sampleBaseDir+"/pat/PATtuple_26_1_Yp0.root",
  sampleBaseDir+"/pat/PATtuple_27_1_Lnr.root",
  sampleBaseDir+"/pat/PATtuple_28_2_f8S.root",
  sampleBaseDir+"/pat/PATtuple_29_1_BhV.root",
  sampleBaseDir+"/pat/PATtuple_30_3_Rqz.root",
  sampleBaseDir+"/pat/PATtuple_31_1_lKz.root",
  sampleBaseDir+"/pat/PATtuple_32_1_UIA.root",
  sampleBaseDir+"/pat/PATtuple_33_1_omI.root",
  sampleBaseDir+"/pat/PATtuple_34_1_faD.root",
  sampleBaseDir+"/pat/PATtuple_35_2_9Zc.root",
  sampleBaseDir+"/pat/PATtuple_36_1_ip0.root",
  sampleBaseDir+"/pat/PATtuple_37_2_AuM.root",
  sampleBaseDir+"/pat/PATtuple_38_1_0Rd.root",
  sampleBaseDir+"/pat/PATtuple_39_1_7oe.root",
  sampleBaseDir+"/pat/PATtuple_40_1_5uu.root",
  sampleBaseDir+"/pat/PATtuple_41_1_tmL.root",
  sampleBaseDir+"/pat/PATtuple_42_2_tLd.root",
  sampleBaseDir+"/pat/PATtuple_43_1_v6k.root",
  sampleBaseDir+"/pat/PATtuple_44_1_Vrq.root",
  sampleBaseDir+"/pat/PATtuple_45_1_9jq.root",
  sampleBaseDir+"/pat/PATtuple_46_1_pEf.root",
  sampleBaseDir+"/pat/PATtuple_47_2_GJw.root",
  sampleBaseDir+"/pat/PATtuple_48_2_SHs.root",
  sampleBaseDir+"/pat/PATtuple_49_2_0mz.root",
  sampleBaseDir+"/pat/PATtuple_50_2_Q1f.root",
  sampleBaseDir+"/pat/PATtuple_51_2_Q1J.root",
  sampleBaseDir+"/pat/PATtuple_52_2_2OP.root",
  sampleBaseDir+"/pat/PATtuple_53_2_17B.root",
  sampleBaseDir+"/pat/PATtuple_54_2_7DM.root",
  sampleBaseDir+"/pat/PATtuple_55_1_fSj.root",
  sampleBaseDir+"/pat/PATtuple_56_2_vve.root",
  sampleBaseDir+"/pat/PATtuple_57_2_trt.root",
  sampleBaseDir+"/pat/PATtuple_58_2_h0J.root",
  sampleBaseDir+"/pat/PATtuple_59_1_KU0.root",
  sampleBaseDir+"/pat/PATtuple_60_1_w3E.root",
  sampleBaseDir+"/pat/PATtuple_61_2_KiC.root",
  sampleBaseDir+"/pat/PATtuple_62_1_5Pu.root",
  sampleBaseDir+"/pat/PATtuple_63_2_Yk5.root",
  sampleBaseDir+"/pat/PATtuple_64_2_Nsr.root",
  sampleBaseDir+"/pat/PATtuple_65_1_A97.root",
  sampleBaseDir+"/pat/PATtuple_66_2_5rk.root",
  sampleBaseDir+"/pat/PATtuple_67_2_zWE.root",
  sampleBaseDir+"/pat/PATtuple_68_2_HcM.root",
  sampleBaseDir+"/pat/PATtuple_69_1_D2U.root",
  sampleBaseDir+"/pat/PATtuple_70_1_khT.root",
  sampleBaseDir+"/pat/PATtuple_71_2_PCY.root",
  sampleBaseDir+"/pat/PATtuple_72_2_ZK4.root",
  sampleBaseDir+"/pat/PATtuple_73_2_T4y.root",
  sampleBaseDir+"/pat/PATtuple_74_2_WMZ.root",
  sampleBaseDir+"/pat/PATtuple_75_1_O2Q.root",
  sampleBaseDir+"/pat/PATtuple_76_1_s1x.root",
  sampleBaseDir+"/pat/PATtuple_77_1_K07.root",
  sampleBaseDir+"/pat/PATtuple_78_1_wdL.root",
  sampleBaseDir+"/pat/PATtuple_79_1_jgr.root",
  sampleBaseDir+"/pat/PATtuple_80_2_0bu.root",
  sampleBaseDir+"/pat/PATtuple_81_1_IT5.root",
  sampleBaseDir+"/pat/PATtuple_82_1_LsL.root",
  sampleBaseDir+"/pat/PATtuple_83_1_ggl.root",
  sampleBaseDir+"/pat/PATtuple_84_1_lLD.root",
  sampleBaseDir+"/pat/PATtuple_85_1_6Tc.root",
  sampleBaseDir+"/pat/PATtuple_86_1_GML.root",
  sampleBaseDir+"/pat/PATtuple_87_1_CRY.root",
  sampleBaseDir+"/pat/PATtuple_88_1_h6I.root",
  sampleBaseDir+"/pat/PATtuple_89_1_YcN.root",
  sampleBaseDir+"/pat/PATtuple_90_2_1Ks.root",
  sampleBaseDir+"/pat/PATtuple_91_1_Bjk.root",
  sampleBaseDir+"/pat/PATtuple_92_2_G7V.root",
  sampleBaseDir+"/pat/PATtuple_93_1_wWj.root"
    ]
