#!/bin/sh

i=1

while [ $i -le 60 ]; do

    FILENAME=Fast_qcd_170_230_puhl_$i

    #  echo i = $i

    # source seed
    a=$[123456789+$i]
    # Vtxsmeared seed
    b=$[123456789+$i]
    # Famos Pileup
    l=$[13578+$i]
    # Famos SimHits
    c=$[13579+$i]
    # siTrackerGaussianSmearingRecHits
    d=$[24680+$i]
    # caloRecHits
    e=$[654321+$i]
    # paramMuons
    f=$[54525+$i]
    # siPixelDigis
    g=$[135711+$i]
    # siStripDigis
    h=$[13171923+$i]

    cat Fast_qcd_170_230_puhl.cfg | sed s/AAAAAAAA/$a/g | sed s/BBBBBBBB/$b/g | sed s/LLLLLLLL/$l/g | sed s/CCCCCCCC/$c/g | sed s/DDDDDDDD/$d/g | sed s/EEEEEEEE/$e/g | sed s/FFFFFFFF/$f/g | sed s/GGGGGGGG/$g/g | sed s/HHHHHHHH/$h/g | sed s/INDEX/$i/g > Fast_qcd_170_230_puhl_$i.cfg

    #  echo "Running cmsRun Fast_qcd_170_230_puhl_$i.cfg > log_qcd_170_230_$i.txt"
    #  cmsRun Fast_qcd_170_230_puhl_$i.cfg > log_qcd_170_230_$i.txt

    cat batch_qcd_170_230_puhl.csh | sed s/INDEX/$i/g > batch_qcd_170_230_puhl_$i.csh

    # Batch job submission
    echo Submitting qcd_170_230_puhl_$i
    bsub -R "pool>40" -q 8nh -J Fast_qcd_170_230_puhl_$i < batch_qcd_170_230_puhl_$i.csh

    #  bsub -R "pool>10000" -q 1nd -J tth_puhl_$i < batch_tth_puhl_$i.csh

    i=$[$i+1]

done
