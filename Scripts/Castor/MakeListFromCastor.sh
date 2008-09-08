#!/bin/sh

if [ "$1" == "" ]; then
    echo Please, input castor directory containing the root files
else

    rfdir $1 | awk '{print $9}'> temp_filelist.txt

    if [ -e castor.cff ]; then
        rm castor.cff
    fi

    # Use wc (word count) to evaluate the number of lines
    # -l counts the newlines
    totlines=`cat temp_filelist.txt | wc -l`

    echo "source = PoolSource {" >> castor.cff
    echo "  untracked vstring fileNames = {" >> castor.cff

    i=1

    for line in `cat temp_filelist.txt`; do
        if [ "$i" != "$totlines" ]; then
#            echo line = $line
#            echo $line | awk '{print $8}'
            echo "    \""castor:$1/$line"\"," >> castor.cff
        else
            # Do not put the "," at the end of the last file
            echo "    \""castor:$1/$line"\"" >> castor.cff
        fi
        i=$[$i+1]
    done

    echo "  }" >> castor.cff
    echo "}" >> castor.cff

    rm temp_filelist.txt

fi
