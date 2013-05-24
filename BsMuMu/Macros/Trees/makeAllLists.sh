#!/bin/sh

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/BsToMuMu_EtaPtFilter_8TeV-pythia6-evtgen/BsToMuMu_BsFilter_8TeV-pythia6-evtgen_Summer12_DR53X-PU_S10_START53_V7A-v1_Onia2MuMuPAT_v13 BsMC
python removeDuplicates.py BsMC

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012A-13Jul2012-v1_Onia2MuMuPAT_v13 Run2012A
python removeDuplicates.py Run2012A

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012ARecover-06Aug2012-v1_Onia2MuMuPAT_v13 Run2012ARecover
python removeDuplicates.py Run2012ARecover

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012B-13Jul2012-v1_Onia2MuMuPAT_v13 Run2012B
python removeDuplicates.py Run2012B

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012C-24Aug2012-v1_Onia2MuMuPAT_v13 Run2012C1
python removeDuplicates.py Run2012C1

#python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012C-EcalRecover-11Dec2012-v1_Onia2MuMuPAT_v13 Run2012CEcalRecover
#python removeDuplicates.py Run2012CEcalRecover
#
python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012C-PromptReco-v2_Onia2MuMuPAT_v13 Run2012C2
python removeDuplicates.py Run2012C2

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012D-v1_Onia2MuMuPAT_v13 Run2012D
python removeDuplicates.py Run2012D

python makeList.py /pnfs/cms/WAX/11/store/user/zhenhu/MuOnia/MuOnia_Run2012DRereco-v1_Onia2MuMuPAT_v13 Run2012DRereco
python removeDuplicates.py Run2012DRereco
