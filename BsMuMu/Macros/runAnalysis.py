#!/bin/env python
from setdirs import *
from mvalib import *

"""
executes the mva selection analysis
the methods implementation proper is in mvalib.py 
the common variables and imports are in setdirs.py
"""

def main():

    print "starting analysis..."
    print "select the analysis step to be executed!"
    ##SELECT what is to be done:

    print "applying preselection, trigger, muonid..."
    # doSelection()

    print "starting mva classification..."
    # doMVATraining("3000") # "0" will train all the events 

    print "doing mva comparisons..."
    doComparisons()
    doComparisonsExtra()

    print "mva significance optimization..."
    # doSignificance()

    print "applying mva selection"
    # doApplication()

    print "drawing mva output and mass"
    # doDrawMVA()
    print "...ending analysis"

def doSelection():
    for isplit in range(-1,3):
        appendName = "_preselection.root"
        if isplit != -1:            
            appendName = appendName.split(".")[0] +"_"+str(isplit)+".root"

        print "processing blinded sample", isplit 
        applySelectionAndSplit(inputTrees, isplit, maxRun, True)   
        combineSamples(appendName)
        addMuonID(appendName)

        # print "processing unblinded sample", isplit 
        applySelectionAndSplit(inputTrees, isplit, maxRun, False)
        appendName = appendName.split(".")[0] +"_unblinded"+".root"
        combineSamples(appendName)
        addMuonID(appendName)


def doMVATraining(eventsToTrain):

    mvastr=""
    for im,method in enumerate(methods):
        if im!=0:mvastr+=','
        mvastr += method

    #run a single classifier
    #runClassification("BDT")

    #OR run all default classifiers
    runClassification(mvastr, eventsToTrain)




optimalCut = { "barrel"  : { "BDT" : [-99,-99,-99], "MLP" : [-99,-99,-99], "CutsSA" : [-99,-99,-99] },
               "endcaps" : { "BDT" : [-99,-99,-99], "MLP" : [-99,-99,-99], "CutsSA" : [-99,-99,-99] } }


def doComparisons():

    for region in regions:

        #make signal vs background distribution comparisons of MVA input variables
        TMVAPlots(region)

        prepareTable(region, "IdTransformation")
        prepareMultiTable(region, "IdTransformation")
        for method in methods:
            prepareMultiTable(region, method)
        saveBDTControlPlots(region)
    
    getNumEvtTable()


def doComparisonsExtra():

    # Variable comparison plots
    os.system("./Common/RunAllComp.sh Common "+rootExecutable)
    os.system("./Common/RunAllComp_3h.sh Common "+rootExecutable)
    os.system("./runMvas.sh")



def doSignificance():

    expectedYield = { "barrel" : { "signal" : estimateSignal("barrel"), "background" : estimateBackground("barrel") }, "endcaps" : { "signal" : estimateSignal("endcaps"), "background" : estimateBackground("endcaps") } }
    print "expectedYields:", expectedYield
    sampleIndex = ["","0","1","2","merged"]

    for method in methods:
        for region in regions:

        ## note: CutsSA CANNOT be merged!, cannot add efficiencies!
        ## TBD: there is no overtraining as such, but could use different sets of cuts using different samples
        ## skip this for now
            
            # merge TMVA classification ouputs for the 3 subsamples
            cmd = rootExecutable+" -l -b -q mergeTMVAs.C+\(\\\""+method+"\\\",\\\""+region+"\\\"\)"
            print cmd
            if not "Cuts" in method:
                os.system(cmd)

                # execute the significance macro
                # main case of interest: merged (bdt,mlp,...), full (not merged, index "") for Cuts
                # extra cases: subsamples, to compare roc's
                for ff in sampleIndex:
                    if "Cuts" in method and ff=="merged": ### for CUTS there is no "merged" file
                        continue
                    cmd = rootExecutable+" -l -b -q significance.C+\("+str(expectedYield[region]["signal"])+","+str(expectedYield[region]["background"])+",\\\""+method+"\\\",\\\""+region+"\\\",\\\""+ff+"\\\",1\)"
                    print cmd
                    os.system(cmd)
                    
        #extra: run signficance also for the separate subsamples, eg to comapre roc's
       ###     cmd = rootExecutable+" -l -b -q significance.C+\("+str(expectedYield[region]["signal"])+","+str(expectedYield[region]["background"])+",\\\""+method+"\\\",\\\""+region+"\\\",\\\""+str(ii)+"\\\",1\)"
       ###     print cmd
       ###     os.system(cmd)
                    
                    # retrieve the optimal cut value which maximizes the significance
                    #   note different figures of merit for estimating the signficance are available 
                    maxSigFile = open(logsDir+"maxsignificance_"+method+"_"+region+".txt")
                    signi = maxSigFile.readline()
        #signi = "a 1 b 2 c 3"
                    for isig in range(3):
                        optimalCut[region][method][isig] = signi.split()[2*isig+1]    
                        
                        

#printMvaCut(optimalCut,regions,methods)

trainedOnIAppliedonJ = ["2","0","1"]
# Applied on 2, trained on 0 (tested on 1)
# Applied on 0, trained on 1 (tested on 2)
# Applied on 1, trained on 2 (tested on 0)

#sample = 1
#region="barrel"
#weightDir = region+str(sample)+"Weights/"
#applyMVA("Barrel_preselection.root", "TMVApp"+"Barrel"+"BDT"+".root", weightDir, "BDT", 0.1)

def doApplication():

    for region in regions:            
        Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
        for method in methods:
            #if method is not "BDT": 
            #    continue
            optCut = optimalCut[region][method][0]
            for sample in range(3):
                weightDir = region+str(sample)+"Weights/"
                applyMVA(Region + "_preselection_"+trainedOnIAppliedonJ[sample]+".root", "TMVApp"+Region+method+str(sample)+".root", weightDir, method, optCut)
                applyMVA("BsMC12_"+region+"_preselection_"+trainedOnIAppliedonJ[sample]+".root", "BsMC12TMVApp"+Region+method+str(sample)+".root", weightDir, method, optCut)

            #applying optimal "merged" cut on the total sample is only for fun
            applyMVA(Region + "_preselection.root", "TMVApp"+Region+method+".root", region+"Weights/", method, optCut)
            applyMVA("BsMC12_"+region+"_preselection.root", "BsMC12TMVApp"+Region+method+".root", region+"Weights/", method, optCut)


def doDrawMVA():

    for region in regions:
        for method in methods:
            if not "Cuts" in method:
                drawAppMVAOutputPlots(region,method,0)
                drawAppMVAOutputPlots(region,method,1)
            drawVariablePlots(region,method)
            drawMassPlot(region,method)


if __name__=="__main__":
   main()



