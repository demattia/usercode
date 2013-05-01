
import os
from ROOT import TROOT
from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import THStack



methods = ["BDT"]
regions = ["barrel"]
sampleIndex = [""]

#methods = ["BDT", "MLP", "CutsSA"]
#regions = ["barrel", "endcaps"]
#sampleIndex = ["","0","1","2"]

figuresDir = "BsMuMuLatex/Figures/"
tablesDir =  "BsMuMuLatex/Tables/"
rootExecutable = "root"



optimalCut = { "barrel"  : { "BDT" : [0.1,0.1,0.1], "MLP" : [0.1,0.1,0.1], "CutsSA" : [0.1,0.1,0.1] },
               "endcaps" : { "BDT" : [0.1,0.1,0.1], "MLP" : [0.1,0.1,0.1], "CutsSA" : [0.1,0.1,-99] } }


def applyMVA(inputFileName, outputFileName, weightDir, method, cutValue):
    os.system(rootExecutable + " -q -l TMVAClassificationApplication.C+\(\\\""+inputFileName+"\\\",\\\""+outputFileName+"\\\",\\\""+weightDir+"\\\","+str(cutValue)+",\\\""+method+"\\\"\)")
#    p = subprocess.Popen([rootExecutable, "-q", "-l", "TMVAClassificationApplication.C+(\""+inputFileName+"\",\""+outputFileName+"\",\""+weightDir+"\","+str(cutValue)+",\""+method+"\")"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
 #   out, err = p.communicate()
    # print out
    # print err


def drawAppMVAOutputPlots(region,method,isMC):

    MCstr = ""
    if isMC:
       MCstr = "BsMC12"
    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    HistName = "ApplicationOutput"+MCstr+Region+method
    CvsName  = "canvasMVA_"+MCstr+method+region
    canvas = TCanvas(CvsName,CvsName)
    stackBDT = THStack(HistName,HistName)
    histos = []
    appFiles = []
    canvas.Draw()

    # note: TMVApp file name structure: "TMVApp"+Region+method+sample+".root"
    #                             "BsMC12TMVApp"+Region+method+sample+".root"
    appFile = TFile(MCstr+"TMVApp"+Region+method+".root")
    histo = appFile.Get("MVA_"+method).Clone(HistName)
    histo.Scale(1/histo.GetEntries())
    stackBDT.Add(histo)

    for sample in range(3):
        appFiles.append(TFile(MCstr+"TMVApp"+Region+method+str(sample)+".root"))
        histos.append(appFiles[sample].Get("MVA_"+method).Clone("HistName"+str(sample)))
        histos[sample].SetLineColor(2*sample)
        histos[sample].Scale(1/histos[sample].GetEntries())
        stackBDT.Add(histos[sample])


    stackBDT.Draw("nostack")
    if isMC:
        applicationBDTLegend = TLegend(0.2,0.7,0.5,0.9,"","brNDC")
        #applicationBDTLegend.SetHeader("Bs MC "+region.split("BsMC")[1])
        applicationBDTLegend.SetHeader("Bs MC "+Region)
    else:
        applicationBDTLegend = TLegend(0.55,0.7,0.85,0.9,"","brNDC")
        applicationBDTLegend.SetHeader(Region)

    applicationBDTLegend.AddEntry(histo, "Full sample", "l")
    applicationBDTLegend.AddEntry(histos[0], "Trained on 0, tested on 1, applied on 2", "l")
    applicationBDTLegend.AddEntry(histos[1], "Trained on 1, tested on 2, applied on 0", "l")
    applicationBDTLegend.AddEntry(histos[2], "Trained on 2, tested on 0, applied on 1", "l")
    applicationBDTLegend.Draw("same")
    applicationBDTLegend.SetFillColor(0)
    applicationBDTLegend.SetLineColor(0)
    canvas.SaveAs(figuresDir+"Application"+method+"Output_"+MCstr+region+".pdf")



#sampleIndex = ["","0","1","2"] ## note not needed (nor precise) to apply the optimal cut to the full "" sample (?..)
trainedOnIAppliedonJ = ["2","0","1"]
# Applied on 2, trained on 0 (tested on 1)
# Applied on 0, trained on 1 (tested on 2)
# Applied on 1, trained on 2 (tested on 0)


redoApplication = False

if redoApplication:

    for region in regions:            
        Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
        for method in methods:
            if method is not "BDT": 
                continue
            optCut = optimalCut[region][method][0]
            for sample in range(3):
                weightDir = region+str(sample)+"Weights/"
                applyMVA(Region + "_preselection_"+trainedOnIAppliedonJ[sample]+".root", "TMVApp"+Region+method+str(sample)+".root", weightDir, method, optCut)
                applyMVA("BsMC12_"+region+"_preselection_"+trainedOnIAppliedonJ[sample]+".root", "BsMC12TMVApp"+Region+method+str(sample)+".root", weightDir, method, optCut)

            #applying optimal "merged" cut on the total sample is only for fun
            applyMVA(Region + "_preselection.root", "TMVApp"+Region+method+".root", region+"Weights/", method, optCut)
            applyMVA("BsMC12_"+region+"_preselection.root", "BsMC12TMVApp"+Region+method+".root", region+"Weights/", method, optCut)


drawAppMVAOutputPlots("barrel","BDT",0)


for region in regions:
    for method in methods:
        drawAppMVAOutputPlots(region,method,0)
        drawAppMVAOutputPlots(region,method,1)

exit(42)


exit(45)


methods = ["BDT", "MLP", "CutsSA"]
regions = ["barrel", "endcaps"]
sampleIndex = ["","0","1","2"]

figuresDir = "BsMuMuLatex/Figures/"
tablesDir =  "BsMuMuLatex/Tables/"
rootExecutable = "root"

expectedEventsBarrel  = 60
estimatedBackgroundBarrel = 28884
expectedEventsEndcaps = 35
estimatedBackgroundEndcaps  = 35392

rootExecutable = "root"
#methods = ["BDT"]
#regions = ["barrel"]
methods = ["BDT", "MLP", "CutsSA"]
regions = ["barrel", "endcaps"]

expectedYield = { "barrel" : { "signal" : expectedEventsBarrel, "background" : estimatedBackgroundBarrel }, "endcaps" : { "signal" : expectedEventsEndcaps, "background" : estimatedBackgroundEndcaps } }

#get list of significance figures-of-merit
def getSignFomName(file):
    signifn = [] 
    maxSigFile = open(file)
    line = maxSigFile.readline()
    for i in range(3):
        print i, 2*i
        signifn.append(line.split()[2*i].replace(":",""))
        maxSigFile.close()
    return signifn

signi_fom = getSignFomName("maxsignificance_BDT_barrel.txt")
#print signi_fom

optimalCut = { "barrel"  : { "BDT" : [-99,-99,-99], "MLP" : [-99,-99,-99], "CutsSA" : [-99,-99,-99] },
               "endcaps" : { "BDT" : [-99,-99,-99], "MLP" : [-99,-99,-99], "CutsSA" : [-99,-99,-99] } }
               
for method in methods:
    for region in regions:

        ## note: CutsSA CANNOT be merged!, cannot add efficiencies!
        ## TBD: there is no overtraining as such, but could use different sets of cuts using different samples
        ## skip this for now

        # merge TMVA classification ouputs for the 3 subsamples
        cmd = rootExecutable+" -l -b -q mergeTMVAs.C\(\\\""+method+"\\\",\\\""+region+"\\\"\)"
        print cmd
        if not "Cuts" in method:
            os.system(cmd)

        # execute the significance macro
        # main case of interest: merged (bdt,mlp,...), full (not merged, index "") for Cuts
        # extra cases: subsamples, to comapre roc's
        findex = ["","0","1","2","merged"]
        for ff in findex:
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
        maxSigFile = open("maxsignificance_"+method+"_"+region+".txt")
        signi = maxSigFile.readline()
        #signi = "a 1 b 2 c 3"
        for isig in range(3):
            optimalCut[region][method][isig] = signi.split()[2*isig+1]    


def printMvaCut(arr):
    print "optimal cuts:",
    for rr in regions:
        print "\n\t",rr,
        for mm in methods:
            print "\n\t\t",mm,
            for ii in range(3):
                print "\t",signi_fom[ii],":% 5.3f" % float(arr[rr][mm][ii]), 

#print optimalCut
printMvaCut(optimalCut)
        

exit(23)


os.system(rootExecutable+" -l -b -q significance.C++\("+str(expectedEventsBarrel)+","+str(estimatedBackgroundBarrel)+",\\\"CutsSA\\\",\\\"barrel\\\",\\\"\\\",0\)")
maxSigFile = open("maxsignificance_CutsSA_barrel.txt")
line = maxSigFile.readline()
print line
print "optimal Cut Barrel S/sqrt(S+B):",  line.split()[1], ", S/sqrt(B):",  line.split()[3],  ", S/sqrt(B)+0.5:",  line.split()[5]  
optimalCutBarrel = line.split()[0]

