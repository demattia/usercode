import os,sys

def combineFiles( type, puFiles ):
    
    command = "hadd -f PileupFiles/pileup_"+type+".root "
    
    for file in puFiles:
        command += file+" "
        pass
    
    print "Running command :",command
    os.system(command)
    pass

def checkFileExists( file ):
    if not os.path.exists(file) :
        print "File :",file,"does not exist"
        sys.exit(1)
        pass

def getLumi( lumiFile ):
    openFile = open(lumiFile, "r")
    lumiLines = openFile.readlines()
    openFile.close()
    for line in lumiLines:
        if line.find("Recorded Lumi")>=0 :
            return float( line.split(":")[-1].strip() ) / 1000000
            pass
        pass
    
    print "ERROR : did not find recorded lumi in",lumiFile
    sys.exit(1)
    
def writeLumi( type, lumi ):
    outputFile = open("LumiFiles/lumi_"+type+".txt","w")
    outputFile.write(str(lumi))
    outputFile.close()
    pass


#===============================================================================
# Read FilesAndWeights.h and work out which data periods I am running on
#===============================================================================
# Find original FilesAndWeights.h
originalFile=open("FilesAndWeights.h","r")
originalFileLines=originalFile.readlines()
originalFile.close()

muonFiles = []
electronFiles = []

muonLumi = 0
electronLumi = 0

for line in originalFileLines:
    # Look for data lines
    if line.find("Data")>=0:
        # Not interested if commented out
        firstPart=line.split("fw")[0]
        if firstPart.find("//")>=0: continue
     
        sampleName = ( line.split("_analysis")[0] ).split("/")[-1]
        puFile = "PileupFiles/pileup_"+sampleName+".root"
        lumiFile = "LumiFiles/lumi_"+sampleName+".txt"

        # Check files exist
        checkFileExists(puFile)
        checkFileExists(lumiFile)
        
        # Is this electrons or muons?
        if sampleName.find("Mu")>=0:
            # Muons
            print "Found muon data sample being considered :",sampleName
            muonFiles.append(puFile)
            
            # get lumi
            muonLumi += getLumi( lumiFile )
            pass
        elif sampleName.find("Photon")>=0:
            # Photons
            print "Found photon data sample being considered :",sampleName
            electronFiles.append(puFile)
            
            # get lumi
            electronLumi += getLumi( lumiFile )
            pass
        else : 
            print "Unkown Data type ??? ",sampleName
            pass
        pass
    pass

#===============================================================================
# Combine pileup files so I have one pileup.root file for electrons/muons
#===============================================================================

combineFiles( "muon", muonFiles )
combineFiles( "electron", electronFiles )

#===============================================================================
# Record total lumi for each channel
# In format that can be easily read in root macro???
#===============================================================================

writeLumi( "muon", muonLumi )
writeLumi( "electron", electronLumi )

print "Total muon lumi :",muonLumi,"/pb"
print "Total electron lumi :",electronLumi,"/pb"
