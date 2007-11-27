#!/bin/sh

if [ "$1" != "" ]; then
    dir=$1
else
    echo "Please, specify castor directory to copy all active files from"
    exit
fi

# Remove the last character in the string if it is a "/" (bash specific command used)
if [[ $dir == */ ]]; then
    localDir=${dir:0:${#dir}-1}
else
    localDir=$dir
fi
# Take the last part of the string
localDir=` echo $localDir | awk -F/ '{print $NF}'`
echo after, local directory = $localDir

mkdir -p $localDir
cd $localDir

#rfdir $dir | awk '{print $9}' > temp_filelist.txt
stager_qry -M $dir | grep STAGED | awk '{print $1}' > temp_filelist.txt

#echo file list = `cat temp_filelist.txt`

#if [ -s temp_filelist.txt ]; then
#    echo No staged files found in the directory,
#    echo please use the Activator.sh to call the stagein of the files.
#else


    for line in `cat temp_filelist.txt`; do
        echo copying file $line
        rfcp $line .
    done

    rm temp_filelist.txt

#fi
