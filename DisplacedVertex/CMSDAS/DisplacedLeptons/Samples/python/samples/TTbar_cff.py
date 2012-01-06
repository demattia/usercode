sampleDataSet = '/TT_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'
sampleCMSEnergy = 7000

sampleRelease = "CMSSW_4_2_2_patch2" # original (i.e. RECO file) release, not the one we plan to process them with
sampleProcessRelease = "CMSSW_4_2_7" # release used to create new files with

sampleNumEvents = 1089625

sampleXSec = 94.0 # pb

# global tag can be extracted from file using edmProvDump filename|grep globaltag
# note however that this is the tag for *further* processing, not the original tag
sampleGlobalTag = 'START42_V13::All'
sampleHLTProcess = '*'

sampleBaseDir = "root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/"+sampleProcessRelease+"/TTbar"

sampleRecoFiles = [ ]

samplePatFiles = [
  sampleBaseDir+"/pat/PATtuple_1_1_4CB.root",
  sampleBaseDir+"/pat/PATtuple_2_1_02a.root",
  sampleBaseDir+"/pat/PATtuple_3_1_Von.root",
  sampleBaseDir+"/pat/PATtuple_4_1_cV3.root",
  sampleBaseDir+"/pat/PATtuple_5_1_LpX.root",
  sampleBaseDir+"/pat/PATtuple_6_1_n4x.root",
  sampleBaseDir+"/pat/PATtuple_7_1_QUD.root",
  sampleBaseDir+"/pat/PATtuple_8_1_oGC.root",
  sampleBaseDir+"/pat/PATtuple_10_1_fbv.root",
  sampleBaseDir+"/pat/PATtuple_11_1_4La.root",
  sampleBaseDir+"/pat/PATtuple_12_1_ylw.root",
  sampleBaseDir+"/pat/PATtuple_13_1_Xm1.root",
  sampleBaseDir+"/pat/PATtuple_14_1_54O.root",
  sampleBaseDir+"/pat/PATtuple_15_1_7GK.root",
  sampleBaseDir+"/pat/PATtuple_16_1_Xmm.root",
  sampleBaseDir+"/pat/PATtuple_17_1_gdL.root",
  sampleBaseDir+"/pat/PATtuple_18_1_vSX.root",
  sampleBaseDir+"/pat/PATtuple_19_1_1bt.root",
  sampleBaseDir+"/pat/PATtuple_20_1_SUv.root",
  sampleBaseDir+"/pat/PATtuple_21_1_ILB.root",
  sampleBaseDir+"/pat/PATtuple_22_1_lbr.root",
  sampleBaseDir+"/pat/PATtuple_23_1_duL.root",
  sampleBaseDir+"/pat/PATtuple_24_1_sBF.root",
  sampleBaseDir+"/pat/PATtuple_25_1_KnU.root",
  sampleBaseDir+"/pat/PATtuple_26_1_xsr.root",
  sampleBaseDir+"/pat/PATtuple_27_1_MKx.root",
  sampleBaseDir+"/pat/PATtuple_28_1_FrP.root",
  sampleBaseDir+"/pat/PATtuple_29_1_mMf.root",
  sampleBaseDir+"/pat/PATtuple_30_1_xF8.root",
  sampleBaseDir+"/pat/PATtuple_31_1_RiP.root",
  sampleBaseDir+"/pat/PATtuple_32_1_sHz.root",
  sampleBaseDir+"/pat/PATtuple_33_1_aQq.root",
  sampleBaseDir+"/pat/PATtuple_34_1_d5W.root",
  sampleBaseDir+"/pat/PATtuple_35_1_ihh.root",
  sampleBaseDir+"/pat/PATtuple_36_1_s2D.root",
  sampleBaseDir+"/pat/PATtuple_37_1_Och.root",
  sampleBaseDir+"/pat/PATtuple_38_1_LEW.root",
  sampleBaseDir+"/pat/PATtuple_39_1_H8G.root",
  sampleBaseDir+"/pat/PATtuple_40_1_wFL.root",
  sampleBaseDir+"/pat/PATtuple_41_1_7by.root",
  sampleBaseDir+"/pat/PATtuple_42_1_shw.root",
  sampleBaseDir+"/pat/PATtuple_43_1_bDv.root",
  sampleBaseDir+"/pat/PATtuple_44_1_AF8.root",
  sampleBaseDir+"/pat/PATtuple_45_1_tLV.root",
  sampleBaseDir+"/pat/PATtuple_46_1_95U.root",
  sampleBaseDir+"/pat/PATtuple_47_1_8NU.root",
  sampleBaseDir+"/pat/PATtuple_48_1_4z2.root",
  sampleBaseDir+"/pat/PATtuple_49_1_Wm7.root",
  sampleBaseDir+"/pat/PATtuple_50_1_pyO.root",
  sampleBaseDir+"/pat/PATtuple_51_1_e5x.root",
  sampleBaseDir+"/pat/PATtuple_52_1_MYf.root",
  sampleBaseDir+"/pat/PATtuple_53_1_07G.root",
  sampleBaseDir+"/pat/PATtuple_54_1_9lY.root",
  sampleBaseDir+"/pat/PATtuple_55_1_Wf0.root",
  sampleBaseDir+"/pat/PATtuple_56_1_OTm.root",
  sampleBaseDir+"/pat/PATtuple_57_1_YlA.root",
  sampleBaseDir+"/pat/PATtuple_58_1_Ged.root",
  sampleBaseDir+"/pat/PATtuple_59_1_mLN.root",
  sampleBaseDir+"/pat/PATtuple_60_1_r7p.root",
  sampleBaseDir+"/pat/PATtuple_61_1_JGM.root",
  sampleBaseDir+"/pat/PATtuple_62_1_zNo.root",
  sampleBaseDir+"/pat/PATtuple_63_1_tRf.root",
  sampleBaseDir+"/pat/PATtuple_64_1_qqO.root",
  sampleBaseDir+"/pat/PATtuple_65_1_vsO.root",
  sampleBaseDir+"/pat/PATtuple_66_1_x45.root",
  sampleBaseDir+"/pat/PATtuple_67_1_5L0.root",
  sampleBaseDir+"/pat/PATtuple_68_1_yu5.root",
  sampleBaseDir+"/pat/PATtuple_69_1_8cP.root",
  sampleBaseDir+"/pat/PATtuple_70_1_0zG.root",
  sampleBaseDir+"/pat/PATtuple_71_1_2kh.root",
  sampleBaseDir+"/pat/PATtuple_72_1_Lmn.root",
  sampleBaseDir+"/pat/PATtuple_73_1_mSl.root",
  sampleBaseDir+"/pat/PATtuple_74_1_A1R.root",
  sampleBaseDir+"/pat/PATtuple_75_1_bxh.root",
  sampleBaseDir+"/pat/PATtuple_76_1_4f0.root",
  sampleBaseDir+"/pat/PATtuple_77_1_0AJ.root",
  sampleBaseDir+"/pat/PATtuple_78_1_74D.root",
  sampleBaseDir+"/pat/PATtuple_79_1_kVK.root",
  sampleBaseDir+"/pat/PATtuple_80_1_cOB.root",
  sampleBaseDir+"/pat/PATtuple_81_1_dMy.root",
  sampleBaseDir+"/pat/PATtuple_82_1_7IO.root",
  sampleBaseDir+"/pat/PATtuple_83_1_ZQ6.root",
  sampleBaseDir+"/pat/PATtuple_84_1_78B.root",
  sampleBaseDir+"/pat/PATtuple_85_1_NGf.root",
  sampleBaseDir+"/pat/PATtuple_86_1_juD.root",
  sampleBaseDir+"/pat/PATtuple_87_1_Jij.root",
  sampleBaseDir+"/pat/PATtuple_88_1_lO2.root",
  sampleBaseDir+"/pat/PATtuple_89_1_TVa.root",
  sampleBaseDir+"/pat/PATtuple_90_1_vVp.root",
  sampleBaseDir+"/pat/PATtuple_91_1_9b4.root",
  sampleBaseDir+"/pat/PATtuple_92_1_SQr.root",
  sampleBaseDir+"/pat/PATtuple_93_1_wne.root",
  sampleBaseDir+"/pat/PATtuple_94_1_2Vs.root",
  sampleBaseDir+"/pat/PATtuple_95_1_1jt.root",
  sampleBaseDir+"/pat/PATtuple_96_1_ch4.root",
  sampleBaseDir+"/pat/PATtuple_97_1_foa.root",
  sampleBaseDir+"/pat/PATtuple_98_1_uJ1.root",
  sampleBaseDir+"/pat/PATtuple_99_1_nmv.root",
  sampleBaseDir+"/pat/PATtuple_100_1_Akl.root",
  sampleBaseDir+"/pat/PATtuple_101_1_sMy.root",
  sampleBaseDir+"/pat/PATtuple_102_1_jsW.root",
  sampleBaseDir+"/pat/PATtuple_103_1_cr6.root",
  sampleBaseDir+"/pat/PATtuple_104_1_baF.root",
  sampleBaseDir+"/pat/PATtuple_105_1_JSr.root",
  sampleBaseDir+"/pat/PATtuple_106_1_8c7.root",
  sampleBaseDir+"/pat/PATtuple_107_1_gXs.root",
  sampleBaseDir+"/pat/PATtuple_108_1_8OT.root",
  sampleBaseDir+"/pat/PATtuple_109_1_9DV.root"
]

sampleDuplicateCheckMode = 'checkAllFilesOpened'

sampleType = "MC"
