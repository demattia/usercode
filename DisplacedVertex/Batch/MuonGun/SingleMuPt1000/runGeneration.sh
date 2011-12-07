#!/bin/sh

i=NUMBER
seed=THESEED

cfgDir=THEDIR
castorDir=THECASTORDIR

eventsPerJob=EVENTSPERJOB

# cd /afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/CMSSW_4_2_2/src
# eval `scramv1 r -sh`
# cd -
# cat ${cfgDir}/gen_template.py | sed s/RANDOMSEED/${seed}/ | sed s/TOTALEVENTS/${eventsPerJob}/ > gen.py
# cmsRun gen.py
# rfcp raw.root ${castorDir}/raw_${i}.root
# 
# cat ${cfgDir}/reco_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
# cmsRun reco.py
# rfcp reco.root ${castorDir}/Reco/reco_${i}.root
# 
# cat ${cfgDir}/reco_refit_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
# cmsRun reco.py
# rfcp reco.root ${castorDir}/Refit/reco_${i}.root
# 
# cat ${cfgDir}/reco_refit_oneIteration_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
# cmsRun reco.py
# rfcp reco.root ${castorDir}/RefitOneIteration/reco_${i}.root
# 
# cat ${cfgDir}/reco_segmentBased_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
# cmsRun reco.py
# rfcp reco.root ${castorDir}/SegmentBased/reco_${i}.root
# 
# cat ${cfgDir}/reco_refit_segmentBased_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
# cmsRun reco.py
# rfcp reco.root ${castorDir}/RefitSegmentBased/reco_${i}.root
# 
# cat ${cfgDir}/reco_refit_oneIteration_segmentBased_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
# cmsRun reco.py
# rfcp reco.root ${castorDir}/RefitOneIterationSegmentBased/reco_${i}.root

cd /afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/ExtendedPropagatorLimit/CMSSW_4_2_2/src
eval `scramv1 r -sh`
cd -

cat ${cfgDir}/reco_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
cmsRun reco.py
rfcp reco.root ${castorDir}/RecoExtendedPropagatorLimit/reco_${i}.root

cat ${cfgDir}/reco_refit_template.py | sed s@INPUTFILE@castor:${castorDir}/raw_${i}.root@ > reco.py
cmsRun reco.py
rfcp reco.root ${castorDir}/RefitExtendedPropagatorLimit/reco_${i}.root
