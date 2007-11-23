#!/bin/sh

if [ "$1" != "" ]; then
  a=$1
else
  echo "Please, specify castor directory to activate"
  exit
fi

rfdir $a | awk '{print $9}' > temp_filelist.txt

#echo file list = `cat temp_filelist.txt`

for line in `cat temp_filelist.txt`; do
  echo stagein request for file = $line
  stager_get -M $1/$line
done
