"""
apply Main BDT on cross-check analysis data input
"""

method="BDT"
Region = ""
cutValue = 0.3

for isample in range(3):
    for ireg in range(2,4):
        if ireg==2:
            Region="Barrel"
            cutValue = 0.3
        elif ireg==3:
            Region="Endcaps"
            cutValue = 0.3
        mainWeightFile = "weights_main/TMVA-"+str(ireg)+"-Events"+str(isample)+"_BDT.weights.xml"
        inputFile = "rootfiles/" + Region + "_preselection_"+str(isample)+".root"
        outputFile = inputFile.split(".")[0] + "_mainBDTapplied.root"
        print isample, Region, mainWeightFile, xcheckInputFile, outputFile
        os.system("root -l -b -q TMVAClassificationApplication_main.C+\(\\\""+inputFile+"\\\",\\\""+outputFile+"\\\",\\\""+mainWeightFile+"\\\","+str(cutValue)+",\\\""+method+"\\\"\)")


"""
TMVA-2-Events0_BDT.weights.xml is for type-0 events in 2012 barrel
TMVA-3-Events0_BDT.weights.xml is for type-0 events in 2013 endcap

TMVA-2-Events1_BDT.weights.xml is for type-1 events in 2012 barrel
TMVA-3-Events1_BDT.weights.xml is for type-1 events in 2013 endcap

"""
