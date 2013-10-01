from setdirs import *
import fileinput

"""
apply Main BDT on cross-check and main analyses datasets
"""

"""
note: the main trees do not necessarily have all cuts applied (eg pvw8, and likely others)!!!
"""

method="BDT"
Region = ""
cutValue = -999 

runOnMainData = False
looseBdtCut = 0 #1:no cut, #2:zero
doRunComp = False #run comparisons only 
doPlots = False
doSeqBdt = False 
doSeqCnt = False 

if len(sys.argv) > 1: 
    if sys.argv[1] == "unblind": unblind = True
    else : unblind = False 
if len(sys.argv) > 2:
    if sys.argv[2] == "runOnMainData": runOnMainData = True
    else : runOnMainData = False
if len(sys.argv) > 3:
    looseBdtCut = int(sys.argv[3]) 
if len(sys.argv) > 4:
    if sys.argv[4] == "doRunComp": doRunComp = True
    else : doRunComp = False
if len(sys.argv) > 5:
    if sys.argv[5] == "doPlots": doPlots = True
    else : doPlots = False
if len(sys.argv) > 6:
    if sys.argv[6] == "doSeqBdt": doSeqBdt = True
    else : doSeqBdt = False
if len(sys.argv) > 7:
    if sys.argv[7] == "doSeqCnt": doSeqCnt = True
    else : doSeqCnt = False 


tree_name = "probe_tree"
anaType = "xcheck"
if runOnMainData: 
    tree_name = "events"
    anaType = "main"

with open("setdirs.h", "a") as myfile:
    if runOnMainData: 
        myfile.write("typedef double mytype;")
    else:
        myfile.write("typedef float mytype;")


strname_ = ""
if runOnMainData:
    strname_ = "_mainBDTonMainData"
else:
    strname_ = "_mainBDTonXcheckData"

if   looseBdtCut == 1: strname_ += "-bdtnone"
elif looseBdtCut == 2: strname_ += "-bdtzero"



#fn_save = figuresDird + strname_.replace("_main","main")
#if not unblind:
#    fn_save += "_blind"
#else:
#    fn_save += "_unblind"
#if not os.path.exists(fn_save):
#    print "creating directory:",fn_save
#    os.system("mkdir -p "+fn_save)


strname_ += ".root"

strname__ = strname_
if unblind:
    strname__ = strname_.replace(".root","_unblinded.root")


def runComp():

    dirname = ""
    #if runOnMainData:
    #    dirname += "mainBDTonMainData"
    #else:
    #    dirname += "mainBDTonXcheckData"

    if   looseBdtCut == 1: dirname += "bdtnone"
    elif looseBdtCut == 2: dirname += "bdtzero"

    if not unblind:
        dirname += "_blind"
    else:
        dirname += "_unblind"
    dirname += "/"

    for region in regions:

        dirname_ = figuresDird + region +"_"+ dirname

        if not os.path.exists(dirname_):
            print "creating directory:",dirname_
            os.system("mkdir -p "+dirname_)


        parallelAnalysisFile = "rootfiles/"+region+"_mainBDTonXcheckData"
        mainAnalysisFile     = "rootfiles/"+region+"_mainBDTonMainData"
        dir = "figs/"+region

        if looseBdtCut == 1:
            parallelAnalysisFile += "-bdtnone"
            mainAnalysisFile     += "-bdtnone"
            dir += "-bdtnone/"

        elif looseBdtCut == 2:
            parallelAnalysisFile += "-bdtzero"
            mainAnalysisFile     += "-bdtzero"
            dir += "-bdtzero/"

        parallelAnalysisFile += ".root"
        mainAnalysisFile     += ".root"

        cmd = "root -l -b -q Common/comp.C+\(\\\""+parallelAnalysisFile+"\\\",\\\""+mainAnalysisFile+"\\\",\\\""+dirname_+"\\\"\)"
        os.system(cmd)
        print cmd
        #print dirname_

    #exit(44)

    #print dir
    #print parallelAnalysisFile

#print fn_save
#exit(88)
#plotsDir = plotsDir + "bdtzero/"


def main():


    if doRunComp:
        runComp()
        exit(70)

    for region in regions:

        #print "REGION:",region
        #if region=="barrel": continue

        if doSeqBdt:
            sequential_bdt(region)
        
        elif doSeqCnt:
            sequential_count(region)

        elif doPlots:
            TMVAPlots(region)

    if doSeqBdt or doSeqCnt or doPlots:
        exit(70)

    outputfilelist = ""
    outputfilelistB = ""
    outputfilelistE= ""



    mainBdtCut = { "barrel":0.360, "endcaps":0.38 }

    if looseBdtCut:
        print looseBdtCut
        if looseBdtCut == 1:
            bdtv=-999
        elif  looseBdtCut == 2:
            bdtv=0
        for r in regions:
            mainBdtCut[r]=bdtv

    print "mainBdtCut:", mainBdtCut
    #exit(66)

    #mainBdtCut = { "barrel":0, "endcaps":0 }
    #mainBdtCut = { "barrel":-9999, "endcaps":-9999 }
    #rootfiles/all_mainBDTonMainData_tree.root 
    #sequential_count("barrel")
    #sequential_count("endcaps")
    #exit(77)

    #main ana BDT application
    for isample in range(3):
        for ireg,region in enumerate(regions):
            #if region == "endcaps" or isample>0 : continue
            Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
            cutValue = mainBdtCut[region]

            mainWeightFile = "weights_main/TMVA-"+str(ireg+2)+"-Events"+str(isample)+"_BDT.weights.xml"

            inputFile = rootDir + Region + "_preselection_"+str(isample)+".root"
            if runOnMainData:
                inputFile = "trees_main/data_afterCuts_"+str(ireg)+"_"+str(isample)+".root"

            if unblind:
                inputFile = inputFile.replace(".root","_unblinded.root")

            outputFile = inputFile.split(".")[0].split("/")[1]
            outputFile = rootDir + outputFile + strname_
            
            outputfilelist += " " + outputFile
            if region is "barrel": outputfilelistB += " " + outputFile
            else:                  outputfilelistE += " " + outputFile

            print isample, Region, mainWeightFile, inputFile, outputFile
            os.system("root -l -b -q TMVAClassificationApplication_main.C+\(\\\""+inputFile+"\\\",\\\""+outputFile+"\\\",\\\""+mainWeightFile+"\\\","+str(cutValue)+",\\\""+method+"\\\",\\\""+region+"\\\",\\\""+tree_name+"\\\",\\\""+anaType+"\\\"\)")


    cmd = "hadd -f "+rootDir+"all"+strname__+outputfilelist
    #print "CMD:",cmd
    os.system(cmd)

    os.system(cmd)
    #print cmd
    #exit(66)
    cmd = "hadd -f "+rootDir+"barrel"+strname__+outputfilelistB
    os.system(cmd)
    #os.system(cmd.replace(".root","_tree.root"))
    cmd = "hadd -f "+rootDir+"endcaps"+strname__+outputfilelistE
    os.system(cmd)
    #os.system(cmd.replace(".root","_tree.root"))

    #exit(66)

    # display plots
    #for region in regions:
    #    TMVAPlots(region)

    #display mass plot and event count after each Cnc cut
    #for region in regions:
    #    sequential_count(region)



preSelCuts = "me< 0.2 && pt > 5. && pt < 9999. && m1pt > 4. && m2pt > 4. && m1eta<2.0 && m1eta>-2.0 && m2eta<2.0 && m2eta>-2.0 && fl3d < 2. && fls3d > 0. && fls3d < 200. && chi2dof < 20. && pvip < 0.1 && pvips < 5. && maxdoca < 0.1 && alpha < 1. && closetrk < 21 && docatrk < 2.5 && iso > 0. && pvlip < 1.0 && pvlip>-1.0 && pvw8>0.7"  # && (lxysig > 3.)  && pvlips < 5.0 && pvlips>-5.0 
if runOnMainData:
    preSelCuts += " && flsxy>3. && pvlips<5.0 && pvlips>-5.0"
else:
    preSelCuts += " && lxysig>3. && pvlip/pvlipErr<5.0 && pvlip/pvlipErr>-5.0"

preSelCuts += " && m>4.9 && m<5.9"




#list of main analysis CnC cuts
listOfCntCutsBarrel = [preSelCuts,#"mu1_GM && mu2_GM",
                 "((m1pt>m2pt && m1pt>4.5 && m2pt>4.0)||(m1pt<m2pt && m1pt>4.0 && m2pt>4.5))",
                 "pt>6.5", "pvw8>0.6","pvip<0.008", "pvips<2",  
                 "alpha<0.05", "chi2dof<2.2", "fls3d>13", "iso>0.8", "docatrk>0.015", "closetrk<2" 
                 ]
listOfCntCutsEndcaps = [preSelCuts,#"mu1_GM && mu2_GM",
                       "((m1pt>m2pt && m1pt>4.5 && m2pt>4.2)||(m1pt<m2pt && m1pt>4.2 && m2pt>4.5))",
                 "pt>8.5", "pvw8>0.6","pvip<0.008", "pvips<2",  
                 "alpha<0.03", "chi2dof<1.8", "fls3d>15", "iso>0.8", "docatrk>0.015", "closetrk<2" 
                 ]



#tbd: additional cuts such as muon_eta and pvw8 [uhm the bdt preselection pvw8 cut is larger than the cnc cut!] should be applied as well, specially as they are not applied to the main trees we have 

#barrel in root format:                 "((m1pt>m2pt && m1pt>4.5 && m2pt>4.0)||(m1pt<m2pt && m1pt>4.0 && m2pt>4.5))&&pt>6.5&&pvw8>0.6&&pvip<0.008&&pvips<2&&alpha<0.05&&chi2dof<2.2&&fls3d>13&&iso>0.8&&docatrk>0.015&&closetrk<2" 
# events->Draw("m", "m1pt>4.5 && m2pt>4.0&&pt>6.5&&pvw8>0.6&&pvip<0.008&&pvips<2&&alpha<0.05&&chi2dof<2.2&&fls3d>13&&iso>0.8&&docatrk>0.015&&closetrk<2")



Ncuts = len(listOfCntCutsBarrel)
    
gStyle.SetOptStat(0)


#def sequential_count(region):
#    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
#    fname = rootDir + Region+"_preselection.root"
#    print fname
#    inputFile = TFile(rootDir + Region+"_preselection.root")
#    inputDir = inputFile.Get("probe_tree")

def sequential_count(region):
    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    fname = rootDir + Region+"_preselection.root"
    ireg = 1
    if(region == "barrel"):
            ireg = 0
    if runOnMainData:
        fname = "trees_main/data_afterCuts_"+str(ireg)+".root"
    print fname, tree_name
    inputFile = TFile(fname)
    inputDir = inputFile.Get(tree_name)
    #inputDir.ls()

    cuts = ""
#    cuts += preSelCuts

    icol=0
    massPlot = []
    evtcount = []

    #for icut,cut in enumerate(listOfCntCuts):
    for icut in range(Ncuts):
        if icut>0:
            cuts += " && "
        if region=="barrel":
            cuts += listOfCntCutsBarrel[icut]
        elif region == "endcaps":
            cuts += listOfCntCutsEndcaps[icut]
        else:
            exit(99)
        #if (icut < (Ncuts) -1 ):
        #    cuts += " && "
        #print  "=== CUTS:", cuts
        #continue
        #exit(55)
        
        hname = "hn"+str(icut)
        inputDir.Draw("m>>"+hname,cuts)
        massPlot.append(TROOT.gROOT.FindObject(hname))
        icol += 1
        if icol%10==0: icol+=1
        massPlot[icut].SetLineColor(icol)
        evtcount.append(massPlot[icut].GetEntries())
        print icut,"++"+cuts+"++", hname, evtcount[icut]
    canvas = TCanvas()
    tl = TLatex()
    tl.SetNDC();
    tl.SetTextSize( 0.033 );
    step = 0.05
    yy = 0.75
    tit = region + " "
    if runOnMainData:
        tit+="(main)"
    else:
        tit+="(xcheck)"
    massPlot[0].SetTitle(tit)
    icol=0
    for icut in range(Ncuts):
        if icut==0: massPlot[icut].Draw()
        else: massPlot[icut].Draw("same")
        yy-=step
        icol += 1
        if icol%10==0: icol+=1
        tl.SetTextColor(icol)
        txt = listOfCntCutsBarrel[icut]
        if region not in "barrel":
            txt = listOfCntCutsEndcaps[icut]
        txt=txt.split("<")[0].split(">")[0].replace("(","")+": "
        if icut==0:
            yy+=2*step
            tl.DrawLatex(0.38,yy,"preselection")
            yy-=step
            txt = ""
        txt += str("  %i" % (evtcount[icut]))
        tl.DrawLatex(0.38,yy,txt)
        if icut==0:
            #tl.DrawLatex(0.38,yy,"xxx")
            yy-=step

    fname  = "mass_cuts_cnt_" + region + strname__.split(".")[0]
    canvas.SaveAs(plotsDir+fname+".gif")
    canvas.SaveAs(figuresDird+fname+".pdf")
    canvass = TCanvas()
    massPlot[Ncuts-1].Draw()
    fname  = "mass_cnt_" + region + strname__.split(".")[0]
    canvass.SaveAs(plotsDir+fname+".gif")
    canvass.SaveAs(figuresDird+fname+".pdf")


def TMVAPlots(region):
    inputFile = TFile(rootDir+region+strname__)
    strname = strname__.split(".")[0]
    #names = [k.GetName() for k in inputFile.GetListOfKeys()]
    #print names
    histoList = inputFile.GetListOfKeys()
    for key in histoList:
        #print key.GetName()
        histo = key.ReadObj()
        if isinstance(histo, TH1F):
            name = histo.GetName()
            #name = fullName.split("__")[0]
            #hType = fullName.split("__")[1].split("_")[0]
            canvas = TCanvas(name+"_canvas")
            #stack = THStack(name+"_stack", name)
            histo.Draw()
            histo.SetLineColor(4)
            histo.SetFillColor(4);
            histo.SetFillStyle(3004);
            #leg = TLegend(0.5536913,0.770979,0.9345638,0.9702797,"","brNDC")
            #leg.AddEntry(histo, hType, "f")
            canvas.Draw()
            canvas.SaveAs(plotsDir+name+"_"+region+strname+".gif")
            canvas.SaveAs(figuresDird+name+"_"+region+strname+".pdf")




def sequential_bdt(region):
    #Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    if runOnMainData:
        print "for sequential bdt, select xcheck tree"
        exit(55)
    if unblind:
        fname = rootDir + region+"_mainBDTonXcheckData-bdtnone_unblinded.root"
    else:
        fname = rootDir + region+"_mainBDTonXcheckData-bdtnone.root"
    print "file:",fname," tree:",tree_name
    inputFile = TFile(fname)
    inputDir = inputFile.Get(tree_name)
    #inputDir.ls()
    #exit(66)

    cuts = ""
    cut = 0

    icol=0
    massPlot = []
    evtcount = []
    Ncuts=20

    tl = TLatex()
    tl.SetNDC();
    tl.SetTextSize( 0.033 );
    tl.SetTextColor(15)
    #icol=0
    #txt=""

    step = 0.1
    cut=0
    #for icut,cut in enumerate(listOfCntCuts):
    for icut in range(Ncuts):
        if icut<10: step= 0.02
        elif icut<16 and cut>0: step=0.02
        elif icut<18: step=0.03
        else: step=0.04
        cut = 0.36 - icut*step
        if region == "endcaps":
            cut = 0.38 - icut*step
        #if icut==0 : cut=0.36
        cuts = "mass>4.85&&m<5.95&&bdt_v>"+str(cut)
        hname = str("hn%02d" % icut)
        #print hname, cuts
        #continue

        cvsi = TCanvas("cvs"+hname)

        inputDir.Draw("m>>"+hname,cuts)
        massPlot.append(TROOT.gROOT.FindObject(hname))
        #icol += 1
        #if icol%10==0: icol+=1
        ###massPlot[icut].SetLineColor(icol)
        evtcount.append(massPlot[icut].GetEntries())
        print icut," -- "+cuts+"--", hname, evtcount[icut]

        massPlot[icut].SetTitle("")
        massPlot[icut].Draw()
        print "Yield:", massPlot[icut].Print()

        #yy-=step
        #icol += 1
        #if icol%10==0: icol+=1
        #tl.SetTextColor(13)
        txt = str("bdt>%3.2f: %d" % (cut, evtcount[icut]))
        tl.DrawLatex(0.7,0.85,txt)
        ffname  = "mass_cuts_bdt_" + region + strname__.split(".")[0]+"_"+str("%02d" % icut)
        cvsi.SaveAs(plotsDir+ffname+".gif")        
        dirloc = figuresDird+"bdtSeq/"
        if not os.path.exists(dirloc):
            print "creating directory:",dirloc
        os.system("mkdir -p "+dirloc)

        cvsi.SaveAs(dirloc+ffname+".pdf")

    #exit(99)
    #
    #massPlot[0].SetTitle(region)
    #
#   # massPlot[Ncuts-1].Draw()
    #
    #fname  = "mass_cuts_bdt_" + region + strname__.split(".")[0]
    #canvas.SaveAs(plotsDir+fname+".gif")
    #canvas.SaveAs(figuresDird+fname+".pdf")
    #canvass = TCanvas()
    #massPlot[0].Draw()
    #fname  = "mass_bdt_" + region + strname__.split(".")[0]
    #canvass.SaveAs(plotsDir+fname+".gif")
    #canvass.SaveAs(figuresDird+fname+".pdf")






"""
"""


if __name__=="__main__":
   main()


"""
MAIN TMVA dictionary:
  TMVA-2-Events0_BDT.weights.xml is for type-0 events in 2012 barrel
  TMVA-3-Events0_BDT.weights.xml is for type-0 events in 2013 endcap
  TMVA-2-Events1_BDT.weights.xml is for type-1 events in 2012 barrel
  TMVA-3-Events1_BDT.weights.xml is for type-1 events in 2013 endcap
  in directory linked as: weights_main
  =>   ireg(=0,1) + 2 -> 2 (barrel), 3 (endcaps) 

MAIN input filename dictionary:
  data_afterCuts_0_0.root is barrel (0) for sample type 0
  data_afterCuts_1_1.root is endcap (1) for sample type 1
  in directory linked as: trees_main
  =>   ireg(=0,1) + 2 -> 2 (barrel), 3 (endcaps)
  
XCHECK input filename (no need for dictionary, self explanatory):
  inputFile = rootDir + Region + \"_preselection_\"+str(isample)+\".root\"

"""
