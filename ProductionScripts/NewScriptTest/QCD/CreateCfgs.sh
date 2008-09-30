#!/bin/sh

fileName="Fast_QCD_120-170.cfg"

names=( "30-50" "50-80" "80-120" "170-230" "230-300" "300-380" "380-incl" )

i=0

while [ $i -lt 7 ]; do

    echo i = ${names[$i]}

    cat $fileName | sed s/120-170/${names[$i]}/g > Fast_QCD_${names[$i]}.cfg

    i=$[$i+1]

done
