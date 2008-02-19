#!/bin/csh
cd /afs/cern.ch/user/d/demattia/scratch0/CMSSW_1_7_5/src
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/d/demattia/scratch0/BATCH/WithTracks/TTH_120/Fast_TTH_120_puhl_INDEX.cfg
rfcp TTH_120_INDEX.root /castor/cern.ch/user/d/demattia/FastSim/PUHL/WithTracks/TTH_120
