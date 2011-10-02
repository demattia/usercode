#!/bin/sh

if [ "$#" -lt 2 ]; then
    echo "need the name of the sample and the starting value for the parallelism cut"
    exit
fi

# exit

rm -f parallelCutEfficiency_${1}.txt
touch parallelCutEfficiency_${1}.txt

startingParallelDiffCut=$2
i=0
numFloat=`echo "(3.4 - $2)*10" | bc -l`
# Convert to int
num=${numFloat/\.*}

while [ "$i" -lt "$num" ]; do
    # echo $i
    parallelDiffCut=`echo "${startingParallelDiffCut} + ${i}/10" | bc -l`
    echo "parallelDiffCut = $parallelDiffCut"

    cat Run_template.C | sed s/SAMPLENAME/${1}/ | sed s/PARALLELISMDIFFCUT/${parallelDiffCut}/ > Run.C
    root -l -b -q Run.C
    root -l -b -q drawCanvases.C
    eff=`cat entriesFile.txt | grep parallelism | awk -F"one hit cut = " '{print $2}' | awk -F"," '{print $1}'`
    echo "$parallelDiffCut $eff" >> parallelCutEfficiency_${1}.txt
    mv entriesFile.txt entriesFile_$i.txt

    let "i += 1"
done
