#!/bin/env python

import os

inputTrees = ['/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012A/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012ARecover/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012B/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012C1/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012CEcalRecover/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012C2/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/Run2012D/selection_test.root',
              '/Users/demattia/TMVA-v4.1.2/test/NewTrees/Trees/BsMC/selection_test.root',
              ]
region = ["barrel", "endcaps"]
for tree in inputTrees:
    data = "1"
    cut_based = "0"
    # NO SPACES IN THE ROOT COMMAND
    for regionIndex in [0,1]:
        if cut_based == "0":
            print "Applying preselection cuts"
            # outputTree = tree.split("_")[1].split(".")[0]+'_'+region[regionIndex]+'_preselection.root'
            outputTree = tree.split("/")[-2]+'_'+region[regionIndex]+'_preselection.root'
            # print outputTree
        else:
            print "Applying analysis cuts"
            # outputTree = tree.split("_")[1].split(".")[0]+'_'+region[regionIndex]+'.root'
            outputTree = tree.split("/")[-2]+'_'+region[regionIndex]+'.root'
            print outputTree
        if tree.find("MC") != -1: data = "0"
        blinding = "0"
        if data == "1": blinding = "1"
        # print 'root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+'\)'
        os.system('root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+'\)')


os.system("rm -f Barrel_preselection.root")
os.system("hadd Barrel_preselection.root Run2012A_barrel_preselection.root Run2012ARecover_barrel_preselection.root Run2012B_barrel_preselection.root Run2012C1_barrel_preselection.root Run2012CEcalRecover_barrel_preselection.root Run2012C2_barrel_preselection.root Run2012D_barrel_preselection.root")
os.system("rm -f Endcaps_preselection.root")
os.system("hadd Endcaps_preselection.root Run2012A_endcaps_preselection.root Run2012ARecover_endcaps_preselection.root Run2012B_endcaps_preselection.root Run2012C1_endcaps_preselection.root Run2012CEcalRecover_endcaps_preselection.root Run2012C2_endcaps_preselection.root Run2012D_endcaps_preselection.root")
os.system("rm -f Run*")
os.system("mv BsMC_barrel_preselection.root BsMC12_barrel_preselection.root")
os.system("mv BsMC_endcaps_preselection.root BsMC12_endcaps_preselection.root")

os.system("root -l -b -q AddMuonID.C+\(\\\"Barrel\\\"\)")
os.system("root -l -b -q AddMuonID.C+\(\\\"Endcaps\\\"\)")

os.system("mv Barrel_preselection_muonID.root Barrel_preselection.root")
os.system("mv Endcaps_preselection_muonID.root Endcaps_preselection.root")

os.system("root -l -b -q AddMuonID.C+\(\\\"BsMC12_barrel\\\"\)")
os.system("root -l -b -q AddMuonID.C+\(\\\"BsMC12_endcaps\\\"\)")

os.system("mv BsMC12_barrel_preselection_muonID.root BsMC12_barrel_preselection.root")
os.system("mv BsMC12_endcaps_preselection_muonID.root BsMC12_endcaps_preselection.root")
