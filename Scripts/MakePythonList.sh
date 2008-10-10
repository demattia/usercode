#!/bin/sh

# Argument is the directory with full path, including the / at the end

if [ "$1" == "" ]; then
    echo Please provide a directory, full path
    exit
else

    ls $1/ > temp_filelist.txt

    if [ -e local_cff.py ]; then
        rm local_cff.py
    fi

    # Use wc (word count) to evaluate the number of lines
    # -l counts the newlines
    totlines=`cat temp_filelist.txt | wc -l`

    echo "import FWCore.ParameterSet.Config as cms" >> local_cff.py
    echo "">> local_cff.py
    echo "source = cms.Source(\"PoolSource\"," >> local_cff.py
    echo "                    fileNames = cms.untracked.vstring(" >> local_cff.py

    i=1

    for line in `cat temp_filelist.txt`; do
        if [ "$i" != "$totlines" ]; then
            echo "    \""file:$1$line"\"," >> local_cff.py
        else
            # Do not put the "," at the end of the last file
            echo "    \""file:$1$line"\"" >> local_cff.py
        fi
        i=$[$i+1]
    done

    echo "    )" >> local_cff.py
    echo "                    )" >> local_cff.py

    rm temp_filelist.txt

fi
