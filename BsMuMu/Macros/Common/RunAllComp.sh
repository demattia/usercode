#!/bin/sh

rm -rf figs

# Data

mkdir figs
rm -rf Data_barrel_figs_OLD
mv Data_barrel_figs Data_barrel_figs_OLD
root -l -b -q comp.C++\(\"Barrel_preselection.root\",\"data_main_barrel.root\"\)
mv figs Data_barrel_figs
cat page_template.php | sed s/"TITLE"/"Data Barrel"/g > Data_barrel_figs/page.php

mkdir figs
rm -rf Data_endcaps_figs_OLD
mv Data_endcaps_figs Data_endcaps_figs_OLD
root -l -b -q comp.C++\(\"Endcaps_preselection.root\",\"data_main_endcaps.root\"\)
mv figs Data_endcaps_figs
cat page_template.php | sed s/"TITLE"/"Data Endcaps"/g > Data_endcaps_figs/page.php

# MC

mkdir figs
rm -rf MC_barrel_figs_OLD
mv MC_barrel_figs MC_barrel_figs_OLD
root -l -b -q comp.C++\(\"BsMC12_barrel_preselection.root\",\"MC_main_barrel.root\"\)
mv figs MC_barrel_figs
cat page_template.php | sed s/"TITLE"/"MC Barrel"/g > MC_barrel_figs/page.php

mkdir figs
rm -rf MC_endcaps_figs_OLD
mv MC_endcaps_figs MC_endcaps_figs_OLD
root -l -b -q comp.C++\(\"BsMC12_endcaps_preselection.root\",\"MC_main_endcaps.root\"\)
mv figs MC_endcaps_figs
cat page_template.php | sed s/"TITLE"/"MC Endcaps"/g > MC_endcaps_figs/page.php
