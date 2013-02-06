#!/bin/env python

import os

inputTrees = ['selection_BsMC12.root',
              # 'selection_BsRunA.root'
              ]
cut_based = "0"
region = ["barrel", "endcaps"]
for tree in inputTrees:
    data = "1"
    # NO SPACES IN THE ROOT COMMAND
    for regionIndex in [0,1]:
        outputTree = tree.split("_")[1].split(".")[0]+'_'+region[regionIndex]+'_preselection.root'
        if tree.find("MC") != -1: data = "0"
        os.system('root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+'\)')

