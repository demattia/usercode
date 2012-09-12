#!/bin/bash

mv masses_data_electrons.txt masses_data_electrons_old.txt
mv masses_data_Muons.txt masses_data_Muons_old.txt
mv masses_backgroundMC_electrons.txt masses_backgroundMC_electrons_old.txt
mv masses_backgroundMC_Muons.txt masses_backgroundMC_Muons_old.txt

# Update cross sections from cff fiels
# python updateWeights.py SampleFiles/Samples/python/

# For PU reweighting and total lumi
# Work out which data periods I am running on and produce pileup_<lepton>.root file
# python combinePileupAndLumiFiles.py

./runTreeAnalyzer.sh

# # # Electrons - Track based
# root -l -b -q 'fullAnalyzer.C++(true)'
root -l -b -q 'fullCombination.C++(true)'
root -l -b -q 'makeAllPlots.C++(true)'

# # # Muons - Track based
# root -l -b -q 'fullAnalyzer.C++'
root -l -b -q 'fullCombination.C++'
root -l -b -q 'makeAllPlots.C++'


 
# # # Muons - Global
# root -l -b -q 'fullAnalyzer.C+(false,false,false)'
# root -l -b -q 'fullCombination.C+(false,false,false)'
# root -l -b -q 'makeAllPlots.C+(false,false,false)'
# 
# # # Muons - Refitted SA
# root -l -b -q 'fullAnalyzer.C+(false,false,true)'
# root -l -b -q 'fullCombination.C+(false,false,true)'
# root -l -b -q 'makeAllPlots.C+(false,false,true)'
