#!/bin/env python

import os

inputTrees = ['/Users/demattia/TMVA-v4.1.2OLD/test/selection_BsMC12.root',
              '/Users/demattia/TMVA-v4.1.2OLD/test/selection_BsRunA.root',
              '/Users/demattia/TMVA-v4.1.2OLD/test/selection_BsRunB.root',
              '/Users/demattia/TMVA-v4.1.2OLD/test/selection_BsRunC1.root'
              ]
region = ["barrel", "endcaps"]
for tree in inputTrees:
    data = "1"
    cut_based = "0"
    # NO SPACES IN THE ROOT COMMAND
    for regionIndex in [0,1]:
        if cut_based == "0":
            print "Applying preselection cuts"
            outputTree = tree.split("_")[1].split(".")[0]+'_'+region[regionIndex]+'_preselection.root'
        else:
            print "Applying analysis cuts"
            outputTree = tree.split("_")[1].split(".")[0]+'_'+region[regionIndex]+'.root'
        if tree.find("MC") != -1: data = "0"
        blinding = "0"
        if data == "1": blinding = "1"
        # print 'root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+'\)'
        os.system('root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+'\)')

