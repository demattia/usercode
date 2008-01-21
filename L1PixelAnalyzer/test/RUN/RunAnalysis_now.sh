#!/bin/sh

echo
echo Running TDAna on all unweighted samples not multiplied
echo

# Select the numer of events
#if [ "$1" == "" ]; then
#    numEvents=-1
#else
#    numEvents=$1
#fi

i=0

bins=( 
       TTH_120
       TT_0JETS 
       TT_1JETS 
       TT_2JETS 
       TT_3JETS 
       TT_4JETS 
       TTBAR
       SINGLETOP_TQ_TQB_LHC_E
       SINGLETOP_TQ_TQB_LHC_MU
       SINGLETOP_TQ_TQB_LHC_TAU
        )

# Loop on no w samples
while [ $i -le 9 ]; do

    cat TDAna_template_now.cfg | sed s/TYPE/${bins[${i}]}/g > TDAna_${bins[${i}]}_now.cfg
    cat batch_TDAnaTemplate_now.csh | sed s/TYPE/${bins[${i}]}/g > batch_TDAna_${bins[${i}]}_now.csh

    echo Submitting job for ${bins[${i}]}
    bsub -R "pool>40" -q 8nh -J ${bins[${i}]} < batch_TDAna_${bins[${i}]}_now.csh

    i=$[$i+1]

done
