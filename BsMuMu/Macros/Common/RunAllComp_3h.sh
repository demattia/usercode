#!/bin/sh

workDir=${1:-.} 
rootExecutable=${2:-root} 

cd ${workDir}
rm -rf figs_3h

# Data

mkdir figs_3h
rm -rf Data_barrel_figs_3h_OLD
mv Data_barrel_figs_3h Data_barrel_figs_3h_OLD
${rootExecutable} -l -b -q comp_3h.C++\(\"../rootfiles/Barrel_preselection_0.root\",\"../rootfiles/Barrel_preselection_1.root\",\"../rootfiles/Barrel_preselection_2.root\"\)
mv figs_3h Data_barrel_figs_3h
cp -rf Data_barrel_figs_3h ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"Data Barrel"/g > Data_barrel_figs_3h/page.php

mkdir figs_3h
rm -rf Data_endcaps_figs_3h_OLD
mv Data_endcaps_figs_3h Data_endcaps_figs_3h_OLD
${rootExecutable} -l -b -q comp_3h.C++\(\"../rootfiles/Endcaps_preselection_0.root\",\"../rootfiles/Endcaps_preselection_1.root\",\"../rootfiles/Endcaps_preselection_2.root\"\)
mv figs_3h Data_endcaps_figs_3h
cp -rf Data_endcaps_figs_3h/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"Data Endcaps"/g > Data_endcaps_figs_3h/page.php

# MC

mkdir figs_3h
rm -rf MC_barrel_figs_3h_OLD
mv MC_barrel_figs_3h MC_barrel_figs_3h_OLD
${rootExecutable} -l -b -q comp_3h.C++\(\"../rootfiles/BsMC12_barrel_preselection_0.root\",\"../rootfiles/BsMC12_barrel_preselection_1.root\",\"../rootfiles/BsMC12_barrel_preselection_2.root\"\)
mv figs_3h MC_barrel_figs_3h
cp -rf MC_barrel_figs_3h/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"MC Barrel"/g > MC_barrel_figs_3h/page.php

mkdir figs_3h
rm -rf MC_endcaps_figs_3h_OLD
mv MC_endcaps_figs_3h MC_endcaps_figs_3h_OLD
${rootExecutable} -l -b -q comp_3h.C++\(\"../rootfiles/BsMC12_endcaps_preselection_0.root\",\"../rootfiles/BsMC12_endcaps_preselection_1.root\",\"../rootfiles/BsMC12_endcaps_preselection_2.root\"\)
mv figs_3h MC_endcaps_figs_3h
cp -rf MC_endcaps_figs_3h/*.pdf ../BsMuMuLatex/Figures/VariablesComparison/
cat page_template.php | sed s/"TITLE"/"MC Endcaps"/g > MC_endcaps_figs_3h/page.php

cd -
