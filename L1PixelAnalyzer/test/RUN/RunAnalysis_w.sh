#!/bin/sh

echo
echo Running TDAna on all weighted samples
echo

# Select the numer of events
#if [ "$1" == "" ]; then
#    numEvents=-1
#else
#    numEvents=$1
#fi

i=0

bins=( 
       QCD_30-50 
       QCD_50-80 
       QCD_80-120 
       QCD_120-170 
       QCD_170-230 
       QCD_230-300 
       QCD_300-380 
       QCD_380-incl 
       EXTRA_QCD_30-50 
       EXTRA_QCD_50-80 
       EXTRA_QCD_80-120 
       EXTRA_QCD_120-170 
       EXTRA_QCD_170-230 
       EXTRA_QCD_230-300 
       EXTRA_QCD_300-380 
       EXTRA_QCD_380-incl 
       ZNUNUJETS_120-170
       ZNUNUJETS_170-230
       W_0JETS
       W_1JETS_0ptw100 
       W_1JETS_100ptw300 
       W_2JETS_0ptw100 
       W_2JETS_100ptw300 
       W_3JETS_0ptw100 
       W_3JETS_100ptw300 
       W_4JETS_0ptw100 
       W_4JETS_100ptw300
       W_5JETS_0ptw100 
       W_5JETS_100ptw300
        )

# Loop on qcd samples
while [ $i -le 28 ]; do

    cat TDAna_template_w.cfg | sed s/TYPE/${bins[${i}]}/g > TDAna_${bins[${i}]}_w.cfg
    cat batch_TDAnaTemplate_w.csh | sed s/TYPE/${bins[${i}]}/g > batch_TDAna_${bins[${i}]}_w.csh

    echo Submitting job for ${bins[${i}]}
    bsub -R "pool>40" -q 8nh -J ${bins[${i}]} < batch_TDAna_${bins[${i}]}_w.csh

    i=$[$i+1]

done
