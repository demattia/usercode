#!/bin/sh

i=NUMBER
seed=THESEED

cd /afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/CMSSW_4_2_3/src
eval `scramv1 r -sh`
cd -
cat /afs/cern.ch/user/d/demattia/scratch0/Batch/HLT/hltSetupMC.py | sed s/EVENTSTOSKIP/${seed}/ | sed s/TOTALEVENTS/THEEVENTS/ > hlt.py
cmsRun hlt.py
rfcp outputHLTDQM.root /castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/outputHLTDQM_${i}.root
cp /afs/cern.ch/user/d/demattia/scratch0/Batch/HLT/rereco_sim_pp.py .
cmsRun rereco_sim_pp.py
rfcp reco_HLT.root /castor/cern.ch/user/d/demattia/DisplacedVertex/MC/Cosmics/reco_HLT_${i}.root
