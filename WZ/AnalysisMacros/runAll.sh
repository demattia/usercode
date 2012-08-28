#!/bin/bash

# # Electrons
# root -l -b -q 'fullAnalyzer.C+(true)'
# root -l -b -q 'fullCombination.C+(true)'
# root -l -b -q 'makeAllPlots.C+(true)'

# # Muons
root -l -b -q 'fullAnalyzer.C++'
root -l -b -q 'fullCombination.C+'
root -l -b -q 'makeAllPlots.C+'
