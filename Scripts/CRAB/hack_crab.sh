#!/bin/bash

# Hacks crab to get around the sandbox limit in 1_8_4 for FastSim.
# Works only with crab 2_1_2, newer versions do not use a directory
# modified by this script.
#
# How to use:
# move the root files in FastSimulation/MaterialEffects/data
# out of the data dir. Create the crab job, use the script and then
# submit the job.


if test -z "$1"; then
  CRABDIR=`ls -drt crab_0_* | tail -1`
else
  CRABDIR=$1
fi

echo "CRABDIR: " $CRABDIR 
cd $CRABDIR/share/.boss_cache/BossTask_1_subDir/ 
tar tf BossArchive_1_g0.tar > _file_list.txt
tar xf BossArchive_1_g0.tar
sed -e 's@# END OF PREPARE AND RUN EXECUTABLE@\n\ cd src/FastSimulation/MaterialEffects/data/ \n\ pwd \n\ ls -lrt \n\ wget -i download.url \n\ cd ../../../../ \n\ echo "### hack"\n@' CMSSW.sh > CMSSW.sh-new 
/bin/mv CMSSW.sh-new CMSSW.sh 
/bin/rm BossArchive_1_g0.tar 
tar cf BossArchive_1_g0.tar `cat _file_list.txt` 
/bin/rm `cat _file_list.txt` 
/bin/rm _file_list.txt 
