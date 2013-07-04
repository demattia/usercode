from setdirs import *
#from mvalib import *


"""
apply Main BDT on cross-check analysis data input
"""

#rootDir="rootfiles/"

method="BDT"
Region = ""
cutValue = -999

runOnMainData = True

tree_name = "probe_tree"
if runOnMainData: tree_name = "events"

with open("setdirs.h", "a") as myfile:
    if runOnMainData: 
        myfile.write("typedef double mytype;")
    else:
        myfile.write("typedef float mytype;")

#exit(44)

strname_ = ""
if runOnMainData:
    strname_ = "_mainBDTonMainData.root"
else:
    strname_ = "_mainBDTonXcheckData.root"

def main():

    outputfilelist = ""
    outputfilelistB = ""
    outputfilelistE= ""

    mainBdtCut = { "barrel":0.360, "endcaps":0.368 }
    #mainBdtCut = { "barrel":0, "endcaps":0 }

    #rootfiles/all_mainBDTonMainData_tree.root 

    #sequential_count("barrel")
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

            outputFile = inputFile.split(".")[0].split("/")[1]
            outputFile = rootDir + outputFile + strname_

            outputfilelist += " " + outputFile
            if region is "barrel": outputfilelistB += " " + outputFile
            else:                  outputfilelistE += " " + outputFile

            print isample, Region, mainWeightFile, inputFile, outputFile
            os.system("root -l -b -q TMVAClassificationApplication_main.C+\(\\\""+inputFile+"\\\",\\\""+outputFile+"\\\",\\\""+mainWeightFile+"\\\","+str(cutValue)+",\\\""+method+"\\\",\\\""+region+"\\\",\\\""+tree_name+"\\\"\)")

#    exit(55)

    cmd = "hadd -f "+rootDir+"all"+strname_+outputfilelist
    os.system(cmd)
    #cmd=cmd.replace(".root","_tree.root")
    os.system(cmd)
    #print cmd
    #exit(66)
    cmd = "hadd -f "+rootDir+"barrel"+strname_+outputfilelistB
    os.system(cmd)
    #os.system(cmd.replace(".root","_tree.root"))
    cmd = "hadd -f "+rootDir+"endcaps"+strname_+outputfilelistE
    os.system(cmd)
    #os.system(cmd.replace(".root","_tree.root"))

    exit(66)


    # display plots
    for region in regions:
        TMVAPlots(region)

    #display mass plot and event count after each Cnc cut
    for region in regions:
        sequential_count(region)



#list of main analysis CnC cuts
listOfCntCutsBarrel = ["","mu1_GM && mu2_GM",
                 "((m1pt>m2pt && m1pt>4.5 && m2pt>4.0)||(m1pt<m2pt && m1pt>4.0 && m2pt>4.5))",
                 "pt>6.5", "pvw8>0.6","pvip<0.008", "pvips<2",  
                 "alpha<0.05", "chi2dof<2.2", "fls3d>13", "iso>0.8", "docatrk>0.015", "closetrk<2" 
                 ]
listOfCntCutsEndcaps = ["","mu1_GM && mu2_GM",
                       "((m1pt>m2pt && m1pt>4.5 && m2pt>4.2)||(m1pt<m2pt && m1pt>4.2 && m2pt>4.5))",
                 "pt>8.5", "pvw8>0.6","pvip<0.008", "pvips<2",  
                 "alpha<0.03", "chi2dof<1.8", "fls3d>15", "iso>0.8", "docatrk>0.015", "closetrk<2" 
                 ]

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
    icol=0
    massPlot = []
    evtcount = []

    #for icut,cut in enumerate(listOfCntCuts):
    for icut in range(Ncuts):
        if icut>1:
            cuts += " && "
        if region=="barrel":
            cuts += listOfCntCutsBarrel[icut]
        elif region == "endcaps":
            cuts += listOfCntCutsEndcaps[icut]
        else:
            exit(99)
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
    massPlot[0].SetTitle(region)
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
        txt += str("%i" % evtcount[icut])
        tl.DrawLatex(0.38,yy,txt)
        if icut==0:
            tl.DrawLatex(0.38,yy,"")
            yy-=step

    fname  = "mass_cuts_cnt_" + region + strname_.split(".")[0]
    canvas.SaveAs(plotsDir+fname+".gif")
    canvas.SaveAs(figuresDird+fname+".pdf")
    canvass = TCanvas()
    massPlot[Ncuts-1].Draw()
    fname  = "mass_cnt_" + region + strname_.split(".")[0]
    canvass.SaveAs(plotsDir+fname+".gif")
    canvass.SaveAs(figuresDird+fname+".pdf")


def TMVAPlots(region):
    inputFile = TFile(rootDir+region+strname_)
    strname = strname_.split(".")[0]
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
