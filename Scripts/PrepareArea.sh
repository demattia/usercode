#!/bin/sh

echo Preparing CMSSW_1_7_0_pre11 area with files necessary to run on the ntuples
echo

scramv1 project CMSSW CMSSW_1_7_0_pre11
mkdir CMSSW_1_7_0_pre11/src/AnalysisExamples
project CMSSW
# Copy AnalysisObjects
cvs co UserCode/DeMattia/AnalysisObjects
mv UserCode/DeMattia/AnalysisObjects CMSSW_1_7_0_pre11/src/AnalysisExamples
# Copy AnalysisClasses
cvs co UserCode/DeMattia/AnalysisObjects
mv UserCode/DeMattia/AnalysisObjects CMSSW_1_7_0_pre11/src/AnalysisExamples
# Copy L1PixelAnalyzer
cvs co UserCode/DeMattia/L1PixelAnalyzer
mv UserCode/DeMattia/AnalysisObjects CMSSW_1_7_0_pre11/src/AnalysisExamples

# Remove UserCode dir
rm -r UserCode

# Not sure this works, better to do it outside
#echo Compiling the area
#cd CMSSW_1_7_0_pre11/src
#scramv1 b

echo
echo
echo Process complete, please compile the area with \"scramv1 b\"
echo then do \"eval \`scramv1 r -csh\`\"
echo
