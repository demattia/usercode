#!/bin/bash

initialSeed=123456796


eventsPerJob=2000


castorDir=/castor/cern.ch/user/d/demattia/MC/MuonGun/SingleMuPt1000


dir=`pwd`

# rfmkdir ${castorDir}
# rfchmod 775 ${castorDir}
# rfmkdir ${castorDir}/Reco
# rfchmod 775 ${castorDir}/Reco
# rfmkdir ${castorDir}/Refit
# rfchmod 775 ${castorDir}/Refit
# rfmkdir ${castorDir}/RefitOneIteration
# rfchmod 775 ${castorDir}/RefitOneIteration
# rfmkdir ${castorDir}/SegmentBased
# rfchmod 775 ${castorDir}/SegmentBased
# rfmkdir ${castorDir}/RefitSegmentBased
# rfchmod 775 ${castorDir}/RefitSegmentBased
# rfmkdir ${castorDir}/RefitOneIterationSegmentBased
# rfchmod 775 ${castorDir}/RefitOneIterationSegmentBased

rfmkdir ${castorDir}/RecoExtendedPropagatorLimit
rfchmod 775 ${castorDir}/RecoExtendedPropagatorLimit
rfmkdir ${castorDir}/RefitExtendedPropagatorLimit
rfchmod 775 ${castorDir}/RefitExtendedPropagatorLimit

echo $dir
echo $castorDir

# for i in {11..50}
for i in {11..12}
  do
  let seed=$initialSeed+$i
  echo "generation number ${i} with seed = ${seed}"
  cat runGeneration.sh | sed s/NUMBER/${i}/ | sed s/THESEED/${seed}/ | sed s@THEDIR@${dir}@ | sed s@THECASTORDIR@${castorDir}@ | sed s@EVENTSPERJOB@${eventsPerJob}@ > runGeneration_${i}.sh
  bsub -R "pool>10000" -q 1nd -J job_Pt1000_${i} < runGeneration_${i}.sh
done
