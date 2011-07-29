#!/bin/sh

i=NUMBER
seed=THESEED

cd /afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/CMSSW_4_2_2/src
eval `scramv1 r -sh`
cd -
cat /afs/cern.ch/user/d/demattia/scratch0/Batch/CosmicGeneration_cfg.py | sed s/RANDOMSEED/${seed}/ > CosmicGeneration.py
cmsRun CosmicGeneration.py
rfcp UndergroundCosmicMu_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root /castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/UndergroundCosmicMu_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_${i}.root
cp /afs/cern.ch/user/d/demattia/scratch0/Batch/rereco_sim_pp.py .
cmsRun rereco_sim_pp.py
rfcp reco_RAW2DIGI_L1Reco_RECOSIM_DQM.root /castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_RAW2DIGI_L1Reco_RECOSIM_DQM_${i}.root
