#!/bin/sh

echo
echo Running TDAna on all multiplied samples
echo

# Select the numer of events
#if [ "$1" == "" ]; then
#    numEvents=-1
#else
#    numEvents=$1
#fi

i=0

bins=( 
       MULTI_EXTRA_QCD_30-50 
       MULTI_EXTRA_QCD_50-80 
       MULTI_EXTRA_QCD_80-120 
       MULTI_EXTRA_QCD_120-170 
       MULTI_EXTRA_QCD_170-230 
       MULTI_EXTRA_QCD_230-300 
       MULTI_EXTRA_QCD_300-380 
       MULTI_EXTRA_QCD_380-incl 
        )

# Loop on qcd samples
while [ $i -le 7 ]; do

    cat TDAna_template_mult.cfg | sed s/TYPE/${bins[${i}]}/g > TDAna_${bins[${i}]}_mult.cfg
    cat batch_TDAnaTemplate_mult.csh | sed s/TYPE/${bins[${i}]}/g > batch_TDAna_${bins[${i}]}_mult.csh

    echo Submitting job for ${bins[${i}]}
    bsub -R "pool>40" -q 8nh -J ${bins[${i}]} < batch_TDAna_${bins[${i}]}_mult.csh

    i=$[$i+1]

done
