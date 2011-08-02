#!/bin/bash

eventsPerJob=10000

for i in {0..200}
  do
  let skipEvents=${i}*${eventsPerJob}
  echo "job number ${i} starting from event number ${skipEvents}"
  cat runGeneration.sh | sed s/NUMBER/${i}/ | sed s/THESEED/${skipEvents}/ | sed s/THEEVENTS/${eventsPerJob}/ > runGeneration_${i}.sh
  bsub -R "pool>10000" -q 1nd -J job_${i} < runGeneration_${i}.sh
done
