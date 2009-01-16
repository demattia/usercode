#!/bin/sh

cvs co -R MuonAnalysis/MomentumScaleCalibration

echo Checked out MuScleFit code

cvs co -R CondCore/PhysicsToolsPlugins
cvs co -R CondFormats/PhysicsToolsObjects
cvs co -R CondFormats/DataRecord
cvs co UserCode/DeMattia/MuScleFit/CondFormats/MomentumScaleCalibrationObjects
mv UserCode/DeMattia/MuScleFit/CondFormats/MomentumScaleCalibrationObjects ./CondFormats
cvs co UserCode/DeMattia/MuScleFit/CondCore/MomentumScaleCalibrationPlugins
mv UserCode/DeMattia/MuScleFit/CondCore/MomentumScaleCalibrationPlugins ./CondCore
cvs co UserCode/DeMattia/MuScleFit/CondFormats/DataRecord
mv UserCode/DeMattia/MuScleFit/CondFormats/DataRecord/interface/MuScleFitLikelihoodPdfRcd.h ./CondFormats/DataRecord/interface/
mv UserCode/DeMattia/MuScleFit/CondFormats/DataRecord/interface/MuScleFitScaleRcd.h ./CondFormats/DataRecord/interface/
mv UserCode/DeMattia/MuScleFit/CondFormats/DataRecord/src/MuScleFitLikelihoodPdfRcd.cc ./CondFormats/DataRecord/src/
mv UserCode/DeMattia/MuScleFit/CondFormats/DataRecord/src/MuScleFitScaleRcd.cc ./CondFormats/DataRecord/src/
rm -rf UserCode

echo Checked out database code
echo
echo Check out complete. Compile the area with \"scramv1 b\".
echo
