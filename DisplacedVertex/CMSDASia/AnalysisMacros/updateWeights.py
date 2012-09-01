import os,sys
from SampleFiles.Tools.AnalysisSample import AnalysisSample

# Take line, work out which sample it is and retrieve weight from sample cff file
def findNewWeight(line,cffDir):
    sampleName=line.split("/")[-1]
    sampleName=sampleName.split("_")[0]
    # Get all info from sample configuration file
    sample=AnalysisSample(cffDir+"/"+sampleName+"_cff.py")
    return sample.sampleXSec


#===============================================================================
# Check command line arguments 
#===============================================================================
if len(sys.argv)<2:
    print "ERROR: need to provide location of sample cff files as argument"
    sys.exit(1)

# Try this first

# Filename 
cffDir=sys.argv[1]

if len(sys.argv)==3:
    if sys.argv[1].find(".py")>0: 
        cffDir=sys.argv[2]
        pass
    pass

        
# Find original FilesAndWeights.h
originalFile=open("FilesAndWeights.h","r")
originalFileLines=originalFile.readlines()
originalFile.close()

newFile=open("FilesAndWeights.h.new","w")


for line in originalFileLines:
    if line.find("Data")>=0 or line.find("Signal")>=0:
        # Weight doesn't change for these
        newFile.write(line)
        pass
    elif line.find("fw[")>=0:
        # Replace weight in this line
        line=line[:line.find("=")]
        weight=findNewWeight(line,cffDir)
        newline=line+"= "+str(weight)+";\n"
        newFile.write(newline)
        pass
    else :
        # Normal line, just write to new file
        newFile.write(line)
        pass
    pass

newFile.close()

os.system("mv FilesAndWeights.h FilesAndWeights.h.old")
os.system("mv FilesAndWeights.h.new FilesAndWeights.h")