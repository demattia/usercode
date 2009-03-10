#!/bin/sh

echo Preparing a CMSSW_2_2_5 area for MuScleFit

scramv1 project CMSSW CMSSW_2_2_5

echo CMSSW area created

cd CMSSW_2_2_5
# eval `scramv1 r -sh`

cvs co -R MuonAnalysis/MomentumScaleCalibration

echo Checked out MuScleFit code

cvs co -R CondFormats/DataRecord
cvs co CondFormats/RecoMuonObjects
cvs co CondCore/RecoMuonPlugins
cvs co CondFormats/DataRecord

echo Checked out database code
# echo Copying Probs file in test
# rfcp /castor/cern.ch/user/d/demattia/MuScleFit/Probs/Probs_new_1000_CTEQ.root MuonAnalysis/MomentumScaleCalibration/test
echo
echo
echo Check out complete. Compile the area with \"scramv1 b\".
echo
