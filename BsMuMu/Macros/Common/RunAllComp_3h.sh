#!/bin/sh

rm -rf figs_3h

# Data

mkdir figs_3h
mkdir figs_3h/delta3d
rm -rf Data_barrel_figs_3h_OLD
mv Data_barrel_figs_3h Data_barrel_figs_3h_OLD
root -l -b -q comp_3h.C++\(\"Barrel_preselection_0.root\",\"Barrel_preselection_1.root\",\"Barrel_preselection_2.root\"\)
mv figs_3h Data_barrel_figs_3h
cat page_template.php | sed s/"TITLE"/"Data Barrel"/g > Data_barrel_figs_3h/page.php

mkdir figs_3h
mkdir figs_3h/delta3d
rm -rf Data_endcaps_figs_3h_OLD
mv Data_endcaps_figs_3h Data_endcaps_figs_3h_OLD
root -l -b -q comp_3h.C++\(\"Endcaps_preselection_0.root\",\"Endcaps_preselection_1.root\",\"Endcaps_preselection_2.root\"\)
mv figs_3h Data_endcaps_figs_3h
cat page_template.php | sed s/"TITLE"/"Data Endcaps"/g > Data_endcaps_figs_3h/page.php

# MC

mkdir figs_3h
mkdir figs_3h/delta3d
rm -rf MC_barrel_figs_3h_OLD
mv MC_barrel_figs_3h MC_barrel_figs_3h_OLD
root -l -b -q comp_3h.C++\(\"BsMC12_barrel_preselection_0.root\",\"BsMC12_barrel_preselection_1.root\",\"BsMC12_barrel_preselection_2.root\"\)
mv figs_3h MC_barrel_figs_3h
cat page_template.php | sed s/"TITLE"/"MC Barrel"/g > MC_barrel_figs_3h/page.php

mkdir figs_3h
mkdir figs_3h/delta3d
rm -rf MC_endcaps_figs_3h_OLD
mv MC_endcaps_figs_3h MC_endcaps_figs_3h_OLD
root -l -b -q comp_3h.C++\(\"BsMC12_endcaps_preselection_0.root\",\"BsMC12_endcaps_preselection_1.root\",\"BsMC12_endcaps_preselection_2.root\"\)
mv figs_3h MC_endcaps_figs_3h
cat page_template.php | sed s/"TITLE"/"MC Endcaps"/g > MC_endcaps_figs_3h/page.php
