#!/bin/sh

MainDir="/Users/demattia/DisplacedVertex/OpenHLT/"

declare -a samples=("part2" "MH120MFF50" "MH1000MFF20" "MH400MFF50")
# declare -a samples=("MH120MFF50")
length=${#samples[@]}

for sample in ${samples[@]}; do
    echo "running on $sample"
    root -l -b -q "computeEfficiency.C+(\"${MainDir}/${sample}/openhlt_merge.root\", \"effMap_${sample}.root\", 1.5, 23, 40)"
done



exit