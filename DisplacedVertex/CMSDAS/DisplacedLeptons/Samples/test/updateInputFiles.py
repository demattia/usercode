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
    if fname.startswith("filelist"):
        startReplacing = False
        sampleName = fname.split(".")[1].split("_pat_")[0]
        print sampleName
        cffFileName = sampleName+"_cff.py"
        if cffFileName in sampleFilesList:
            print "found"
            inputCff = open(samplesDir+cffFileName)
            sampleCffName = sampleName+"_cff_new.py"
            outputCff = open(sampleCffName, "w")
            # Loop over the original cff
            for line in inputCff:
                if startReplacing:
                    break
                if line.startswith("sampleBaseDir"):
                    outputCff.write("sampleBaseDir = \"root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/\"+sampleProcessRelease+\"/Data_Mu_Run2011A1\"")
                elif line.find("harder/longlived/") > 0:
                    pass
                else:
                    outputCff.write(line)
                # stop at the samplePatFiles list
                if line.startswith("samplePatFiles"):
                    startReplacing = True
            # Fill the rest with the new list of files
            fileListFile = open(fname)
            numFiles = file_len(fname)
            counter = 1
            # print "num = "+str(numFiles)
            for line in fileListFile:
                counter = counter + 1
                if counter == numFiles:
                    # print "replacing"
                    outputCff.write(line.replace(",", "\n]"))
                else:
                    outputCff.write(line)
            outputCff.close()
            # print "counter = "+str(counter) 
            os.system("mv "+samplesDir+sampleCffName+" "+samplesDir+cffFileName)

