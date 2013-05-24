#!/bin/env python

import os, sys

rootFiles = []
foundRootFiles = False
outputFile = open("makeSelectionBs_cfg_new.py", "w")
dirName = sys.argv[1]

# For MC
index = 13

# For Data
if dirName.find("BsMC") == -1:
    index = 5


rootFilesDict = {}
for line in open(dirName+"/makeSelectionBs_cfg.py"):
    if line.find(".root") != -1:
        # print line.split('_')[index]
        if line.split('_')[index] in rootFilesDict:
            print "Duplicate file found:", line
        rootFilesDict[line.split('_')[index]] = line

        foundRootFiles = True
    elif foundRootFiles == True:
        for a,rootFile in rootFilesDict.iteritems():
            outputFile.write(rootFile)
        foundRootFiles = False
    if foundRootFiles == False:
        outputFile.write(line)


outputFile.close()
