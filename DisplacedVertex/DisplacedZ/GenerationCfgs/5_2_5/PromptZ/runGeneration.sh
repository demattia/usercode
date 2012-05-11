#!/bin/sh

i=NUMBER
seed=THESEED

cfgDir="/afs/cern.ch/user/d/demattia/scratch0/Batch/DisplacedVertex/PromptZ"
castorDir="/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/PromptZ/5_2_5"
castorGenDir=${castorDir}"/genSim"

echo ${castorGenDir}

cd /afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/DisplacedZ/CMSSW_5_2_5
eval `scramv1 r -sh`
cd -
cat ${cfgDir}/ZMM_8TeV_cfi_GEN_SIM_template.py | sed s/RANDOMSEED/${seed}/ > ZMM_8TeV_cfi_GEN_SIM.py
cmsRun ZMM_8TeV_cfi_GEN_SIM.py
rfcp step1.root ${castorGenDir}/genSim_${i}.root
cmsRun ${cfgDir}/step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.py
cmsRun ${cfgDir}/step3_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM.py
rfcp step3.root ${castorDir}/reco_${i}.root
