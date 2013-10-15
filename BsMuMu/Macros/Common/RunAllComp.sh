#!/bin/sh

workDir=${1:-.} 
rootExecutable=${2:-root} 

cd ${workDir}
rm -rf figs

# Data

mkdir figs
rm -rf Data_barrel_figs_OLD
mv Data_barrel_figs Data_barrel_figs_OLD
${rootExecutable} -l -b -q comp_2h.C++\(\"../rootfiles/Barrel_preselection.root\",\"../trees_main/data_afterCuts_0.root\",0\)
# root -l -b -q comp_2h.C++\(\"Barrel_preselection.root\",\"../data_main_barrel.root\"\)
mv figs Data_barrel_figs
cp -rf Data_barrel_figs ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"Data Barrel"/g > Data_barrel_figs/page.php

mkdir figs
rm -rf Data_endcaps_figs_OLD
mv Data_endcaps_figs Data_endcaps_figs_OLD
${rootExecutable} -l -b -q comp_2h.C++\(\"../rootfiles/Endcaps_preselection.root\",\"../trees_main/data_afterCuts_1.root\",0\)
mv figs Data_endcaps_figs
cp -rf Data_endcaps_figs ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"Data Endcaps"/g > Data_endcaps_figs/page.php

# MC

mkdir figs
rm -rf MC_barrel_figs_OLD
mv MC_barrel_figs MC_barrel_figs_OLD
${rootExecutable} -l -b -q comp_2h.C++\(\"../rootfiles/BsMC12_barrel_preselection.root\",\"../trees_main/MC_afterCuts_0.root\",1\)
mv figs MC_barrel_figs
cp -rf MC_barrel_figs ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"MC Barrel"/g > MC_barrel_figs/page.php

mkdir figs
rm -rf MC_endcaps_figs_OLD
mv MC_endcaps_figs MC_endcaps_figs_OLD
${rootExecutable} -l -b -q comp_2h.C++\(\"../rootfiles/BsMC12_endcaps_preselection.root\",\"../trees_main/MC_afterCuts_1.root\",1\)
mv figs MC_endcaps_figs
cp -rf MC_endcaps_figs ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"MC Endcaps"/g > MC_endcaps_figs/page.php

cd -
