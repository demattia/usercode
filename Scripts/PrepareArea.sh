#!/bin/sh

echo Preparing CMSSW_1_8_4 area with files necessary to run on the ntuples
echo

scramv1 project CMSSW CMSSW_1_8_4
mkdir CMSSW_1_8_4/src/AnalysisExamples
project CMSSW
# Copy AnalysisObjects
cvs co UserCode/DeMattia/AnalysisObjects
mv UserCode/DeMattia/AnalysisObjects CMSSW_1_8_4/src/AnalysisExamples
# Copy AnalysisClasses
cvs co UserCode/DeMattia/AnalysisClasses
mv UserCode/DeMattia/AnalysisClasses CMSSW_1_8_4/src/AnalysisExamples
# Copy L1PixelAnalyzer
cvs co UserCode/DeMattia/L1PixelAnalyzer
mv UserCode/DeMattia/L1PixelAnalyzer CMSSW_1_8_4/src/AnalysisExamples
# Copy PixelJets
cvs co UserCode/DeMattia/PixelJet
mv UserCode/DeMattia/PixelJet CMSSW_1_8_4/src/AnalysisExamples
# Copy PixelJetFinder
cvs co UserCode/DeMattia/PixelJetFinder
mv UserCode/DeMattia/PixelJetFinder CMSSW_1_8_4/src/AnalysisExamples
# Copy ttHMEtplusJetsAnalyzer
cvs co UserCode/DeMattia/ttHMEtplusJetsAnalyzer
mv UserCode/DeMattia/ttHMEtplusJetsAnalyzer CMSSW_1_8_4/src/AnalysisExamples

# Remove UserCode dir
rm -r UserCode

# Comment some files that are not up to date
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/TDAna.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/TDAna.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/TDAna.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/TDAna.h_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/TDAnaProducer.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/TDAnaProducer.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/TDAnaProducer.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/TDAnaProducer.h_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/AlgValidator.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/AlgValidator.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/AlgValidator.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/AlgValidator.h_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/MultiplicationFilter.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/MultiplicationFilter.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/MultiplicationFilter.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/MultiplicationFilter.h_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/PTag.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/PTag.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/PTag.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/PTag.h_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/RecoAssociator.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/RecoAssociator.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/RecoAssociator.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/RecoAssociator.h_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/VertexAssoc.cc CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/src/VertexAssoc.cc_backup
mv CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/VertexAssoc.h CMSSW_1_8_4/src/AnalysisExamples/L1PixelAnalyzer/interface/VertexAssoc.h_backup

# Not sure this works, better to do it outside
#echo Compiling the area
#cd CMSSW_1_8_4/src
#scramv1 b

echo
echo
echo Process complete, please compile the area with \"scramv1 b\"
echo then do \"eval \`scramv1 r -csh\`\"
echo
