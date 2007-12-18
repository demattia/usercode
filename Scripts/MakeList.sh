#!/bin/sh

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
else

    ls /data/demattia/QCD_FILES/$1/ > temp_filelist.txt

    if [ -e $1.cff ]; then
        rm $1.cff
    fi

    # Use wc (word count) to evaluate the number of lines
    # -l counts the newlines
    totlines=`cat temp_filelist.txt | wc -l`

    echo "source = PoolSource {" >> $1.cff
    echo "  untracked vstring fileNames = {" >> $1.cff

    i=1

    for line in `cat temp_filelist.txt`; do
        if [ "$i" != "$totlines" ]; then
            echo "    \""file:/data/demattia/QCD_FILES/$1/$line"\"," >> $1.cff
        else
            # Do not put the "," at the end of the last file
            echo "    \""file:/data/demattia/QCD_FILES/$1/$line"\"" >> $1.cff
        fi
        i=$[$i+1]
    done

    echo "  }" >> $1.cff
    echo "}" >> $1.cff

    rm temp_filelist.txt

fi
