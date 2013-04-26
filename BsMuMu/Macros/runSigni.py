
import os

expectedEventsBarrel  = 60
estimatedBackgroundBarrel = 28884
expectedEventsEndcaps = 35
estimatedBackgroundEndcaps  = 35392


figuresDir = "BsMuMuLatex/Figures/"
tablesDir =  "BsMuMuLatex/Tables/"
rootExecutable = "/Users/nuno/install/root/bin/root"


#os.system(rootExecutable+" -l -b -q significance.C++\("+str(expectedEventsBarrel)+","+str(estimatedBackgroundBarrel)+",\\\"BDT\\\",\\\"barrel\\\",\\\"merged\\\",0\)")
#os.system(rootExecutable+" -l -b -q significance.C++\("+str(expectedEventsEndcaps)+","+str(estimatedBackgroundEndcaps)+",\\\"BDT\\\",\\\"endcaps\\\",\\\"merged\\\",0\)")

os.system(rootExecutable+" -l -b -q significance.C++\("+str(expectedEventsBarrel)+","+str(estimatedBackgroundBarrel)+",\\\"CutsSA\\\",\\\"barrel\\\",\\\"\\\",0\)")

maxSigFile = open("maxsignificance_CutsSA_barrel.txt")
line = maxSigFile.readline()
print line
print "optimal Cut Barrel S/sqrt(S+B):",  line.split()[0], ", S/sqrt(B):",  line.split()[1],  ", S/sqrt(B)+0.5:",  line.split()[2]  
optimalCutBarrel = line.split()[0]

print "ola\n"

