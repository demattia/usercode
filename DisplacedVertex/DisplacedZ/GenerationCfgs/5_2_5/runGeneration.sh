#!/bin/sh

i=NUMBER
seed=THESEED

cfgDir="/afs/cern.ch/user/d/demattia/scratch0/Batch/DisplacedVertex/DisplacedZ"
castorDir="/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5"
castorGenDir=${castorDir}"/genSim"

cd /afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/DisplacedZ/CMSSW_5_2_5
eval `scramv1 r -sh`
cd -
cat ${cfgDir}/SingleNeutralinoGun_template.py | sed s/RANDOMSEED/${seed}/ > SingleNeutralinoGun.py
cmsRun SingleNeutralinoGun.py
rfcp step1.root ${castorGenDir}/genSim_${i}.root
cmsRun ${cfgDir}/step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.py
cmsRun ${cfgDir}/step3_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM.py
rfcp step3.root ${castorDir}/reco_${i}.root
