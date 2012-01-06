#!/bin/bash

#
# Data
#
for i in `ls -1 ${CMSSW_BASE}/src/DisplacedLeptons/Samples/python/samples/*cff.py|grep -v SLHC|grep Data_` ; do
    echo starting `basename $i .py`...
    ${CMSSW_BASE}/src/HarderAnalysis/DisplacedDileptons/test/run_analysis.py $i &> `basename $i .py`.analysislog &
    sleep 1
done

#
# MC and Cosmics
#
for i in `ls -1 ${CMSSW_BASE}/src/DisplacedLeptons/Samples/python/samples/*cff.py|grep -v SLHC|grep -v Data_` ; do
    echo starting `basename $i .py`...
    ${CMSSW_BASE}/src/HarderAnalysis/DisplacedDileptons/test/run_analysis.py $i &> `basename $i .py`.analysislog &
    sleep 1
done

