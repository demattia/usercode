#!/bin/sh

rm -rf plots
mkdir plots

root -l -b -q mvas.C++\(\"rootfiles/TMVA_barrel_0.root\",3,kFALSE,0\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_barrel_0.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_barrel_1.root\",3,kFALSE,0\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_barrel_1.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_barrel_2.root\",3,kFALSE,0\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_barrel_2.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_endcaps_0.root\",3,kFALSE,0\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_endcaps_0.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_endcaps_1.root\",3,kFALSE,0\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_endcaps_1.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_endcaps_2.root\",3,kFALSE,0\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_endcaps_2.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_barrel_0.root\",3,kFALSE,1\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_barrel_0_log.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_barrel_1.root\",3,kFALSE,1\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_barrel_1_log.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_barrel_2.root\",3,kFALSE,1\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_barrel_2_log.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_endcaps_0.root\",3,kFALSE,1\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_endcaps_0_log.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_endcaps_1.root\",3,kFALSE,1\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_endcaps_1_log.pdf

root -l -b -q mvas.C++\(\"rootfiles/TMVA_endcaps_2.root\",3,kFALSE,1\)
mv plots/overtrain_BDT.pdf BsMuMuLatex/Figures/bdt/overtrain_BDT_endcaps_2_log.pdf

