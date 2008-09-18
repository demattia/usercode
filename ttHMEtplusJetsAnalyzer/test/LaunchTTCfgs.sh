#!/bin/sh

if [ "$1" == "" ]; then
    echo no input
fi

# NF is the last position in the line
fileName=`echo $1 | awk -F/ '{print $NF}'`
dirName=`echo $1 | awk -F$fileName '{print $1}'`
#echo fileName = $fileName
#echo dirName = $dirName

names=( "TT_0JETS" "TT_2JETS" "TT_3JETS" )

i=0

while [ $i -lt 3 ]; do

    #echo i = ${names[$i]}

    newDir=`echo $dirName | sed s/TT_1JETS/${names[$i]}/g`
    newFile=`echo $fileName | sed s/TT_1JETS/${names[$i]}/g`

    echo Creating $newDir$newFile

    mkdir -p $newDir
    cat $1 | sed s/TT_1JETS/${names[$i]}/g > $newDir$newFile

    cmsRun $newDir$newFile > log &

    i=$[$i+1]

done

cmsRun $1 > log &
