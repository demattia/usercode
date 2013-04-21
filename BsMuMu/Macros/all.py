#!/bin/env python

import os

inputTrees = ['Trees/Run2012A/selection_test.root',
              'Trees/Run2012ARecover/selection_test.root',
              'Trees/Run2012B/selection_test.root',
              'Trees/Run2012C1/selection_test.root',
              'Trees/Run2012DRereco/selection_test.root',
              'Trees/Run2012C2/selection_test.root',
              'Trees/Run2012D/selection_test.root',
              'Trees/BsMC/selection_test.root',
              ]
region = ["barrel", "endcaps"]

# This is the value of the rest of event number % 3
# splitting = -1
# This is the maximum run number considered for data (included)
maxRun = "203002"


def combineSamplesAndAddMuonID(appendName):
    # Combine barrel samples
    os.system("rm -f Barrel"+appendName)
    os.system("hadd Barrel"+appendName+" Run2012A_barrel"+appendName+" Run2012ARecover_barrel"+appendName+" Run2012B_barrel"+appendName+" Run2012C1_barrel"+appendName+" Run2012DRereco_barrel"+appendName+" Run2012C2_barrel"+appendName+" Run2012D_barrel"+appendName)
    # Combine endcaps samples
    os.system("rm -f Endcaps"+appendName)
    os.system("hadd Endcaps"+appendName+" Run2012A_endcaps"+appendName+" Run2012ARecover_endcaps"+appendName+" Run2012B_endcaps"+appendName+" Run2012C1_endcaps"+appendName+" Run2012DRereco_endcaps"+appendName+" Run2012C2_endcaps"+appendName+" Run2012D_endcaps"+appendName)
    os.system("rm -f Run*")
    # Add MVA muon-id
    os.system("root -l -b -q AddMuonID.C+\(\\\"Barrel"+appendName+"\\\"\)")
    os.system("root -l -b -q AddMuonID.C+\(\\\"Endcaps"+appendName+"\\\"\)")
    os.system("mv Barrel"+appendName+"_muonID.root Barrel"+appendName)
    os.system("mv Endcaps"+appendName+"_muonID.root Endcaps"+appendName)

    # Process MC
    if appendName.find("unblinded") == -1:
        os.system("mv BsMC_barrel"+appendName+" BsMC12_barrel"+appendName)
        os.system("mv BsMC_endcaps"+appendName+" BsMC12_endcaps"+appendName)
        os.system("root -l -b -q AddMuonID.C+\(\\\"BsMC12_barrel"+appendName+"\\\"\)")
        os.system("root -l -b -q AddMuonID.C+\(\\\"BsMC12_endcaps"+appendName+"\\\"\)")
        os.system("mv BsMC12_barrel"+appendName+"_muonID.root BsMC12_barrel"+appendName)
        os.system("mv BsMC12_endcaps"+appendName+"_muonID.root BsMC12_endcaps"+appendName)


def applySelectionAndSplit(inputTrees, region, splitting, maxRun, blindData = True):
    splitString =""
    if splitting != -1:
        splitString = "_"+str(splitting)
    
    appendName = "_preselection"+splitString+".root"
    # appendNameMC = "_preselection.root"
    
    
    for tree in inputTrees:
        data = "1"
        cut_based = "0"
        if tree.find("MC") != -1: data = "0"
        # append = appendName
        # if data == "0": append = appendNameMC
        # NO SPACES IN THE ROOT COMMAND
        for regionIndex in [0,1]:
            if cut_based == "0":
                print "Applying preselection cuts"
                outputTree = tree.split("/")[-2]+'_'+region[regionIndex]+appendName
            else:
                print "Applying analysis cuts"
                outputTree = tree.split("/")[-2]+'_'+region[regionIndex]+'.root'
            blinding = "0"
            if data == "1":
                if blindData:
                    blinding = "1"
                else:
                    outputTree = outputTree.split(".")[0]+"_unblinded.root"
            # print 'root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+','+str(splitting)+','+maxRun+'\)'
            os.system('root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+','+str(splitting)+','+maxRun+'\)')

    if blindData:
        print "Processing blinded"
        combineSamplesAndAddMuonID(appendName)
    elif splitting == -1 and appendName.find("MC") == -1:
        print "Processing unblinded " + appendName.split(".")[0]+"_unblinded.root"
        combineSamplesAndAddMuonID(appendName.split(".")[0]+"_unblinded.root")

# -------------------------------
# Run the selection and splitting

for splitting in range(-1,3):
    applySelectionAndSplit(inputTrees, region, splitting, maxRun)

# Produce unblinded datasets. The last parameter is blindData = False
applySelectionAndSplit(inputTrees, region, -1, maxRun, False)
