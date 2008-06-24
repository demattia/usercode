#!/bin/sh

# This script opens all the files in a directory. If one of the files is not working it will give an error. Change the script so that it can report the error or do something to address it.
# Requires flite_slc4.C to be in the same dir, in order to avoid all the warnings when opening an edm root file.

if [ "$1" == "" ]; then
    echo please provide a castor directory
else
    echo validating castor files in dir: $1

    rfdir $1 | awk '{print $9}' | while read file; do
 
        #if [ `echo $file | grep -c 6` -ne 0 ]; then
        echo validating file = $1/$file
        root -l -b -q flite_slc4.C rfio:$1/$file
        #fi
        #echo $1/$line
    done
    
fi
