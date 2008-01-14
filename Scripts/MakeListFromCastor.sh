#!/bin/sh

DIR=/castor/cern.ch/user/d/demattia/FastSim/PUHL

if [ "$1" == "" ]; then
    echo Please, select the sample,
    echo chose from:
    echo QCD_30-50
    echo QCD_50-80
    echo QCD_80-120
    echo QCD_120-170
    echo QCD_170-230
    echo QCD_230-300
    echo QCD_300-380
    echo QCD_380-incl
    echo TTH_120
    echo EXTRA_QCD_30-50
    echo ...
else

    rfdir ${DIR}/$1/ | awk '{print $9}'> temp_filelist.txt

    if [ -e $1_castor.cff ]; then
        rm $1_castor.cff
    fi

    # Use wc (word count) to evaluate the number of lines
    # -l counts the newlines
    totlines=`cat temp_filelist.txt | wc -l`

    echo "source = PoolSource {" >> $1_castor.cff
    echo "  untracked vstring fileNames = {" >> $1_castor.cff

    i=1

    for line in `cat temp_filelist.txt`; do
        if [ "$i" != "$totlines" ]; then
#            echo line = $line
#            echo $line | awk '{print $8}'
            echo "    \""castor:${DIR}/$1/$line"\"," >> $1_castor.cff
        else
            # Do not put the "," at the end of the last file
            echo "    \""castor:${DIR}/$1/$line"\"" >> $1_castor.cff
        fi
        i=$[$i+1]
    done

    echo "  }" >> $1_castor.cff
    echo "}" >> $1_castor.cff

    rm temp_filelist.txt

fi
