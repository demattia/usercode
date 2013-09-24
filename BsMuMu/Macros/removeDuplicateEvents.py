#!/bin/env python

import os, sys
import ROOT
import json

f = ROOT.TFile("Trees/Run2012D/selection_test.root", "READ")
tree = f.Get("detailedDimuonTree").Get("probe_tree")
jsonFile = open("Cert_207883-208307_8TeV_16Jan2013ReReco_Collisions12_JSON.txt")
jsonObject = json.load(jsonFile)

def checkRunLumi(inputRun, inputLumi):
    for run, lumis in jsonObject.iteritems():
        if int(inputRun) == int(run):
            print "found run", run
            for lumi in lumis:
                if int(inputLumi) >= int(lumi[0]) and inputLumi <= int(lumi[1]):
                    print "Found lumi", inputLumi
                    return False
    return True

# checkRunLumi(208307, 761)

outputFile = ROOT.TFile("Trees/Run2012D/cleaned.root", "RECREATE")
dir = outputFile.mkdir("detailedDimuonTree")
dir.cd()
newtree = tree.CloneTree(0);


# i = 0
for event in tree:
    if checkRunLumi(event.run, event.lumi):
        newtree.Fill()

outputFile.Close()

#     i += 1
# 
#     if i == 100:
#         break

    # print event.pt
