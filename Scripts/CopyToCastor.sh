#!/bin/sh

if [ "$1" != "" ]; then
    dir=$1
else
    echo "Please, specify castor directory to copy all active files to"
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

for file in `ls $localDir`; do
    echo copying file: ${file}
    rfcp ${localDir}/${file} $dir
done
