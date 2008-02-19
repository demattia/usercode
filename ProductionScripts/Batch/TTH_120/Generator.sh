#!/bin/sh

i=10001

while [ $i -le 10010 ]; do

    FILENAME=Fast_TTH_120_puhl_$i

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

    cat Fast_TTH_120_puhl.cfg | sed s/AAAAAAAA/$a/g | sed s/BBBBBBBB/$b/g | sed s/LLLLLLLL/$l/g | sed s/CCCCCCCC/$c/g | sed s/DDDDDDDD/$d/g | sed s/EEEEEEEE/$e/g | sed s/FFFFFFFF/$f/g | sed s/GGGGGGGG/$g/g | sed s/HHHHHHHH/$h/g | sed s/INDEX/$i/g > Fast_TTH_120_puhl_$i.cfg

    #  echo "Running cmsRun Fast_qcd_120_170_puhl_$i.cfg > log_qcd_120_170_$i.txt"
    #  cmsRun Fast_qcd_120_170_puhl_$i.cfg > log_qcd_120_170_$i.txt

    cat batch_TTH_120_puhl.csh | sed s/INDEX/$i/g > batch_TTH_120_puhl_$i.csh

    # Batch job submission
    echo Submitting TTH_120_puhl_$i
    bsub -R "pool>40" -q 8nh -J Fast_TTH_120_puhl_$i < batch_TTH_120_puhl_$i.csh

    #  bsub -R "pool>10000" -q 1nd -J tth_puhl_$i < batch_tth_puhl_$i.csh

    i=$[$i+1]

done
