#!/bin/sh

if [ "$1" == "" ]; then
    echo Please provide a directory, full path
    exit
else

    ls $1/ > temp_filelist.txt

    if [ -e local.cff ]; then
        rm local.cff
    fi

    # Use wc (word count) to evaluate the number of lines
    # -l counts the newlines
    totlines=`cat temp_filelist.txt | wc -l`

    echo "source = PoolSource {" >> local.cff
    echo "  untracked vstring fileNames = {" >> local.cff

    i=1

    for line in `cat temp_filelist.txt`; do
        if [ "$i" != "$totlines" ]; then
            echo "    \""file:$1/$line"\"," >> local.cff
        else
            # Do not put the "," at the end of the last file
            echo "    \""file:$1/$line"\"" >> local.cff
        fi
        i=$[$i+1]
    done

    echo "  }" >> local.cff
    echo "}" >> local.cff

    rm temp_filelist.txt

fi
