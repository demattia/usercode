#!/bin/sh

workDir=${1:-.} 
rootExecutable=${2:-root} 

cd ${workDir}
rm -rf figs

# Data

mkdir figs
rm -rf Data_barrel_figs_OLD
mv Data_barrel_figs Data_barrel_figs_OLD
${rootExecutable} -l -b -q comp.C++\(\"../rootfiles/Barrel_preselection.root\",\"../rootfiles/data_main_barrel.root\"\)
# root -l -b -q comp.C++\(\"Barrel_preselection.root\",\"../data_main_barrel.root\"\)
mv figs Data_barrel_figs
cp Data_barrel_figs/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/Data_barrel_figs/
cat page_template.php | sed s/"TITLE"/"Data Barrel"/g > Data_barrel_figs/page.php

mkdir figs
rm -rf Data_endcaps_figs_OLD
mv Data_endcaps_figs Data_endcaps_figs_OLD
${rootExecutable} -l -b -q comp.C++\(\"../rootfiles/Endcaps_preselection.root\",\"../rootfiles/data_main_endcaps.root\"\)
mv figs Data_endcaps_figs
cp Data_endcaps_figs/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/Data_endcaps_figs/
cat page_template.php | sed s/"TITLE"/"Data Endcaps"/g > Data_endcaps_figs/page.php

# MC

mkdir figs
rm -rf MC_barrel_figs_OLD
mv MC_barrel_figs MC_barrel_figs_OLD
${rootExecutable} -l -b -q comp.C++\(\"../rootfiles/BsMC12_barrel_preselection.root\",\"../rootfiles/MC_main_barrel.root\"\)
mv figs MC_barrel_figs
cp MC_barrel_figs/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/MC_barrel_figs/
cat page_template.php | sed s/"TITLE"/"MC Barrel"/g > MC_barrel_figs/page.php

mkdir figs
rm -rf MC_endcaps_figs_OLD
mv MC_endcaps_figs MC_endcaps_figs_OLD
${rootExecutable} -l -b -q comp.C++\(\"../rootfiles/BsMC12_endcaps_preselection.root\",\"../rootfiles/MC_main_endcaps.root\"\)
mv figs MC_endcaps_figs
cp MC_endcaps_figs/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/MC_endcaps_figs/
cat page_template.php | sed s/"TITLE"/"MC Endcaps"/g > MC_endcaps_figs/page.php

cd -
