#!/bin/sh
rm -f index.html
#wget http://cmstac11.cern.ch:8080/analysis/CrabAnalysis/TIFNtupleMaker/FNAL_pre6_v16/ index.html
wget $1 index.html

cat index.html | awk -F\" '{print $6}' | awk -F/ '{print $1}' > temp_list.txt

cat temp_list.txt | while read line; do

  echo Downloading file ${line}.root
  wget -A.ps -nd ${1}/${line}/res/${line}_1.root ./Files

done
