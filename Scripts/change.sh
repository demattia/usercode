#!/bin/sh

# Performs a replace of parameter 1 with parameter 2 for the
# files in the specified directory

if [ "$3" == "" ]; then
    echo Please select the directory, use "." for the current dir.
    exit
fi

if [ "$1" == "" ] & [ "$2" == "" ]; then
    echo Both argument 1 and 2 are empty. Please provide two arguments
    exit
fi

echo Changing \"$1\" to \"$2\" in files:

cd $3
files=`ls *`
for file in `echo $files`; do
    echo $file
    tempfile=`echo $file | awk '{print $1 ".bak"}'`
    # echo tempfile = $tempfile
    mv $file $tempfile
    sed s/"$1"/"$2"/g $tempfile > $file
    rm $tempfile
done
