#!/usr/bin/env python

import os

def file_len(fname):
    i = 0
    with open(fname) as f:
        for line in f:
            i = i + 1
    return i + 1

dirList = os.listdir(".")
samplesDir = "../python/samples/"
sampleFilesList = os.listdir(samplesDir)
# print sampleFilesList
for fname in dirList:
    # Skip empty filelist files
    if fname.startswith("filelist") and os.path.getsize(os.getcwd()+"/"+fname) > 0:
        replacing = False
        sampleName = fname.split(".")[1].split("_pat_")[0]
        print "Updating file list for sample "+sampleName
        cffFileName = sampleName+"_cff.py"
        if cffFileName in sampleFilesList:
            inputCff = open(samplesDir+cffFileName)
            sampleCffName = samplesDir+sampleName+"_cff_new.py"
            outputCff = open(sampleCffName, "w")
            # Loop over the original cff
            for line in inputCff:
                if replacing:
                    if line.find("]") != -1:
                        replacing = False
                    pass
                elif line.startswith("sampleBaseDir"):
                    outputCff.write("sampleBaseDir = \"root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/\"+sampleProcessRelease+\"/"+sampleName+"\"\n")
                elif line.find("harder/longlived/") > 0:
                    pass
                # stop at the samplePatFiles list
                else:
                    outputCff.write(line)
                    if line.startswith("samplePatFiles"):
                        print "replacing"
                        replacing = True
                        # Fill the rest with the new list of files
                        fileListFile = open(fname)
                        numFiles = file_len(fname)
                        counter = 1
                        for line in fileListFile:
                            counter = counter + 1
                            if counter == numFiles:
                                outputCff.write(line.replace(",", "\n]"))
                            else:
                                outputCff.write(line)

            outputCff.close()
            os.system("mv "+sampleCffName+" "+samplesDir+cffFileName)

