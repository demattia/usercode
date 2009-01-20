#!/bin/sh

# Argument is the directory with full path, including the / at the end

if [ "$1" == "" ]; then
    echo Please provide a directory, full path
    exit
else

    echo Reminder: the path must be provided with the "/" at the end.

    rfdir $1 | awk '{print $9}'> temp_filelist.txt

    if [ -e castor.cff ]; then
        rm castor.cff
    fi

    # Use wc (word count) to evaluate the number of lines
    # -l counts the newlines
    totlines=`cat temp_filelist.txt | wc -l`

    echo "import FWCore.ParameterSet.Config as cms" >> castor_cff.py
    echo "">> castor_cff.py
    echo "source = cms.Source(\"PoolSource\"," >> castor_cff.py
    echo "                    fileNames = cms.untracked.vstring(" >> castor_cff.py

    i=1

    #dir=`echo $1 | sed s-/castor/cern.ch--g`

    for line in `cat temp_filelist.txt`; do
        if [ "$i" != "$totlines" ]; then
            echo "    \""rfio:$1$line"\"," >> castor_cff.py
        else
            # Do not put the "," at the end of the last file
            echo "    \""rfio:$1$line"\"" >> castor_cff.py
        fi
        i=$[$i+1]
    done

    echo "    )" >> castor_cff.py
    echo "                    )" >> castor_cff.py

    rm temp_filelist.txt

fi
