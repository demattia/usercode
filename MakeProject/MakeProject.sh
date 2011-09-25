#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Error, please provide the project include file name"
    exit
fi

echo "Filling the project includes for project dir:"

if [ "$CMSSW_BASE" == "" ]; then
    echo "Error: please run this script after doing cmsenv"
    exit
fi

cd $CMSSW_BASE/src

if [ -e "$1" ]; then
    fileEnd=`echo $1 | awk -F. '{print $2}'`
    if [ "$fileEnd" == "includes" ]; then
	mv "$1" "$1.backup"
	pwd | sed s@/afs/cern.ch@/Users/demattia/afs@g >> $1
	echo "$CMSSW_RELEASE_BASE/src" | sed s@/afs/cern.ch@/Users/demattia/afs@g >> $1
	echo "${ROOTSYS}/include" | sed s@/afs/cern.ch@/Users/demattia/afs@g >> $1
	echo "/afs/cern.ch/cms/slc5_amd64_gcc434/external/boost/1.44.0-cms/include" | sed s@/afs/cern.ch@/Users/demattia/afs@g >> $1
    else
	echo "Error:include file does not end in .includes"
    fi
else
    echo "Error: provided include file name $1 does not exist in $CMSSW_BASE/src"
    exit
fi

cd -
