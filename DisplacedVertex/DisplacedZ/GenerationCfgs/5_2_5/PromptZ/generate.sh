#!/bin/bash

initialSeed=123456796

for i in {1..10}
  do
  let seed=$initialSeed+$i
  echo "generation number ${i} with seed = ${seed}"
  cat runGeneration.sh | sed s/NUMBER/${i}/ | sed s/THESEED/${seed}/ > runGeneration_${i}.sh
  bsub -R "pool>10000" -q 1nd -J job_${i} < runGeneration_${i}.sh
done
