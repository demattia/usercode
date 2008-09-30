#!/bin/sh

fileName="Fast_TT_1JETS.cfg"

names=( "TT_0JETS" "TT_2JETS" "TT_3JETS" "TT_4JETS" )

i=0

while [ $i -lt 4 ]; do

    echo i = ${names[$i]}

    cat $fileName | sed s/TT_1JETS/${names[$i]}/g > Fast_${names[$i]}.cfg

    i=$[$i+1]

done
