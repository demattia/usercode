#!/bin/csh
cd /afs/cern.ch/user/d/demattia/scratch0/CMSSW_1_7_0/src
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/d/demattia/scratch0/BATCH/QCD_170-230/Fast_qcd_170_230_puhl_INDEX.cfg
rfcp QCD_170-230_INDEX.root /castor/cern.ch/user/d/demattia/FastSim/PUHL/QCD_170-230
