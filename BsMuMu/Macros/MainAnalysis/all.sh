#!/bin/sh

root -l -b -q cutTree_BsMuMu.C\(\"0\",\"0\"\)
root -l -b -q cutTree_BsMuMu.C\(\"0\",\"1\"\)
root -l -b -q cutTree_BsMuMu.C\(\"0\",\"2\"\)

rm -f data_afterCuts_0.root
hadd data_afterCuts_0.root data_afterCuts_0_0.root data_afterCuts_0_1.root data_afterCuts_0_2.root

root -l plot.C\(\"0\"\)

root -l -b -q cutTree_BsMuMu.C\(\"1\",\"0\"\)
root -l -b -q cutTree_BsMuMu.C\(\"1\",\"1\"\)
root -l -b -q cutTree_BsMuMu.C\(\"1\",\"2\"\)

rm -f data_afterCuts_1.root
hadd data_afterCuts_1.root data_afterCuts_1_0.root data_afterCuts_1_1.root data_afterCuts_1_2.root

root -l plot.C\(\"1\"\)