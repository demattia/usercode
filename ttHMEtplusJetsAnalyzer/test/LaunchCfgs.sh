#!/bin/sh

if [ "$1" == "" ]; then
    echo no input
fi

# NF is the last position in the line
fileName=`echo $1 | awk -F/ '{print $NF}'`
dirName=`echo $1 | awk -F$fileName '{print $1}'`
#echo fileName = $fileName
#echo dirName = $dirName

names=( "30-50" "50-80" "80-120" "170-230" "230-300" "300-380" "380-incl" )

i=0

while [ $i -lt 7 ]; do

    #echo i = ${names[$i]}

    newDir=`echo $dirName | sed s/120-170/${names[$i]}/g`
    newFile=`echo $fileName | sed s/120-170/${names[$i]}/g`

    echo Creating $newDir$newFile

    mkdir -p $newDir
    cat $1 | sed s/120-170/${names[$i]}/g > $newDir$newFile

    cmsRun $newDir$newFile > log &

    i=$[$i+1]

done

cmsRun $1 > log &
