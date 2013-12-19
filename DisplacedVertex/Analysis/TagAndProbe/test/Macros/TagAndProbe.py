
import ROOT
from utils import *
# import math
from ROOT import RooRealVar, RooFormulaVar, RooVoigtian, RooChebychev, RooArgList, RooArgSet, \
    RooAddPdf, RooDataSet, RooCategory, RooSimultaneous, RooGenericPdf, RooWorkspace
import Workspace
from array import array 
ROOT.gROOT.LoadMacro("Loader.C++")



# --------#


MC = True
NoBkgd = True

# --------#


# Define binning
ptBinsX = [15, 20, 23, 26, 30, 35, 40, 45, 50, 60, 70]
#ptBinsY = [26, 30, 35, 40, 45, 50, 60, 70]
#ptBinsX = [15, 70]
ptBinsY = [20, 10000]

# Define cuts and some useful variables
triggerMatchDeltaR = 0.1
minMass = 10
maxMass = 150
minDeltaR = 0.2
p = Properties(minMass, maxMass, ptBinsX, ptBinsY, triggerMatchDeltaR, NoBkgd, minDeltaR)


# Trigger efficiency for new trigger over old trigger
tagTrigger = "IsoMu24_v"
probeTrigger = "HLT_L2DoubleMu23_NoVertex_v"

# Load the input file
tree = ROOT.TChain("T")
tree.Add("Z_mu_mu_tagAndProbe_MC.root")
#tree.Add("Run_2012B.root")
#tree.Add("Z_mumu_Data_tag_and_probe_Run2012A.root")
# Prepare the workspace
ws = RooWorkspace("ws", "workspace")
Workspace.buildPdf(ws, p)

mass = ws.var("mass")
sample = ws.cat("sample")
simPdf = ws.pdf("simPdf")
efficiency = ws.var("efficiency")
meanB = ws.var("meanB")

# prepare_datasets(ws, p)

# Prepare datasets
datasetAllMap = {}
datasetPassMap = {}
hAllMap = {}
hPassMap = {}
for ptBin1 in range(0, len(ptBinsX)):
    for ptBin2 in range(0, len(ptBinsY)):
        datasetAllMap[(ptBin1, ptBin2)] = RooDataSet(buildName("datasetAll_", ptBin1, ptBin2, ptBinsX, ptBinsY), buildName("datasetAll_", ptBin1, ptBin2, ptBinsX, ptBinsY), RooArgSet(mass))
        datasetPassMap[(ptBin1, ptBin2)] = RooDataSet(buildName("datasetPass_", ptBin1, ptBin2, ptBinsX, ptBinsY), buildName("datasetPass_", ptBin1, ptBin2, ptBinsX, ptBinsY), RooArgSet(mass))
        hAllMap[(ptBin1,ptBin2)] = ROOT.TH1F(buildName("hAll_", ptBin1, ptBin2, ptBinsX, ptBinsY), buildName("All events passing old trigger ", ptBin1, ptBin2, ptBinsX, ptBinsY), 100, 60, 120)
        hPassMap[(ptBin1,ptBin2)] = ROOT.TH1F(buildName("hPass_", ptBin1, ptBin2, ptBinsX, ptBinsY), buildName("All events passing old trigger and new trigger ", ptBin1, ptBin2, ptBinsX, ptBinsY), 100, 60, 120)
#Prepare 1D histogram to get the ratio of number of passing probes over all probes muons

counting_all = ROOT.TH1F("counting_all","counting_all", len(ptBinsX)-1, array('d',ptBinsX))
counting_pass = ROOT.TH1F("counting_pass","counting_pass",len(ptBinsX)-1, array('d',ptBinsX))
counting_eff = ROOT.TH1F("counting_eff_SA_MC_probe_negative","counting_eff", len(ptBinsX)-1,  array('d',ptBinsX))
counting_eff.GetXaxis().SetTitle("Probe muon p_{T} [GeV]")
counting_eff.GetYaxis().SetTitle("Number of Passing Probes/Number of All Probes")

# Event loop
allCandidates = 0
passCandidates = 0

processedEvents = 0

totEvents = 1000000
#totEvents = tree.GetEntries()
progress = 0



for event in tree:
    if processedEvents > totEvents:
        break
    processedEvents += 1

    if int(float(processedEvents)/float(totEvents)*100) % 10 == 0 and int(float(processedEvents)/float(totEvents)*100) > progress:
        progress = 100*processedEvents/totEvents
        progressCounter(progress)
    if processedEvents == totEvents - 1:
        progressCounter(100);


    # Get the trigger outcome for this event
    tagTriggerFired = False
    probeTriggerFired = False
    for name in event.triggerNames:
        if name.find(tagTrigger) != -1: tagTriggerFired = True
    # If none of the two triggers fired we can skip the event
    if not tagTriggerFired :
        continue


    tagTriggerObjects = event.triggerFilterObjectsMap["hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15"]
    probeTriggerObjects = event.triggerFilterObjectsMap["hltL2DoubleMu23NoVertexL2PreFiltered"]
    
    matchedMuonsTagTrigger = []
    passingMuonsProbe = []
    allMuonsProbe = []
#check your input trees...
    for muon in event.muons:
        # Find a matching trigger object in DeltaR
#The boolean variable to check if matching between muon or track objects and trigger objects...
#        min_deltaRTag = 1
#        for genParticle in event.genParticles:
#            deltaRTag = deltaR(muon.phi, muon.eta, genParticle.phi, genParticle.eta) 
#            hDeltaRTag.Fill(deltaRTag)
#            if (deltaRTag < min_deltaRTag): min_deltaRTag = deltaRTag 
#            if deltaR(muon.phi(), muon.eta(), genParticle.phi(), genParticle.eta()) < 0.2: 
#        if min_deltaRTag < 0.2:
        fillTriggerMatchedGlobalMuon(muon, tagTriggerObjects, matchedMuonsTagTrigger, p)
#check input
#replace standalonemuon by standAloneMuon
#replace tracks by standAloneMuons later
#check root file if it is generaltracks...
    
#    for standAloneMuon in event.genParticles:
#        min_deltaRProbe = 1
    for standAloneMuon in event.refittedStandAloneMuons:
#The boolean variable to check if matching between muon or track objects and trigger objects...
    # Find a matching trigger object in DeltaR
#        min_deltaRProbe = 1
#        for genParticle in event.genParticles:
#            deltaRProbe = deltaR(standAloneMuon.phi, standAloneMuon.eta, genParticle.phi, genParticle.eta) 
#            if deltaRProbe < min_deltaRProbe and genParticle.motherid == 23 : min_deltaRProbe = deltaRProbe 
#        if passSelectionStandAlone(standAloneMuon) and  min_deltaRProbe < 0.2 :
#        if passSelectionStandAlone(standAloneMuon) and standAloneMuon.motherid == 23 :
        if passSelectionStandAlone(standAloneMuon) :
            fillTriggerMatchedStandAlone(standAloneMuon, probeTriggerObjects, passingMuonsProbe, p)
            allMuonsProbe.append(standAloneMuon)

#    for track in event.tracks:
    # Find a matching trigger object in DeltaR
#        if passSelection(track):
#            fillTriggerMatchedTrack(track, probeTriggerObjects, passingMuonsProbe, p)
#            allMuonsProbe.append(track)


    fillCandidates_tnp(mass, p, matchedMuonsTagTrigger, passingMuonsProbe, hPassMap, datasetPassMap, counting_pass)
    fillCandidates_tnp(mass, p, matchedMuonsTagTrigger, allMuonsProbe, hAllMap, datasetAllMap, counting_all)

print "all candidates =", hAllMap[1,1].GetEntries()
print "pass candidates =", hPassMap[1,1].GetEntries()

#Define a canvas to see superimposed deltaR distribution...

#deltacanvas = ROOT.TCanvas("deltacanvas","deltacanvas", 800, 400)
#deltacanvas.Divide(2,1)
#deltacanvas.cd(1)
#hDeltaRTag.Draw()
#deltacanvas.cd(2)
#hDeltaRProbe.Draw()
#deltacanvas.Print("deltaRcanvas.png")

counting_eff.Divide(counting_pass, counting_all, 1.0, 1.0)

#counting_canvas_pass = ROOT.TCanvas("counting_canvas_pass","counting_canvas_pass", 800, 400)
#counting_pass.Draw()
#counting_canvas_pass.Print("counting_canvas_pass.png")

#counting_canvas_all = ROOT.TCanvas("counting_canvas_all","counting_canvas_all", 800, 400)
#counting_all.Draw()
#counting_canvas_all.Print("counting_canvas_all.png")

#counting_canvas_eff = ROOT.TCanvas("counting_canvas_eff","counting_canvas_eff", 800, 400)
#counting_eff.Draw()
#counting_canvas_eff.Print("counting_canvas_eff.png")
#Save canvas_eff histo in a root file which updateCanvas_Melih.py will make use of

#output = ROOT.TFile("output_tnp_histos.root","UPDATE")
#counting_eff.Write()
#output.Close()

canvas = ROOT.TCanvas("AllAndPassCanvas", "AllAndPassCanvas", 800, 800)
canvas.Divide(len(ptBinsX),len(ptBinsY))
for ptBin1 in range(0, len(ptBinsX)):
    for ptBin2 in range(0, len(ptBinsY)):
        canvas.cd(find_position(ptBin1, ptBin2, ptBinsX)+1)
        hAllMap[(ptBin1, ptBin2)].Draw()
        hPassMap[(ptBin1, ptBin2)].Draw("same")
        hPassMap[(ptBin1, ptBin2)].SetLineColor(2)

canvas.Print("AllAndPassCanvas.pdf")

# Plot the fit results

def plotResults(ptBin1, ptBin2, combData, canvas2, canvas3):
    frame1 = mass.frame(ROOT.RooFit.Bins(30),ROOT.RooFit.Title("All events"))
    # Plot all data tagged as physics sample
    combData.plotOn(frame1,ROOT.RooFit.Cut("sample==sample::all"))
    # Plot "physics" slice of simultaneous pdf.
    # NBL You _must_ project the sample index category with data using ProjWData
    # as a RooSimultaneous makes no prediction on the shape in the index category
    # and can thus not be integrated
    simPdf.plotOn(frame1, ROOT.RooFit.Slice(sample,"all"), ROOT.RooFit.ProjWData(RooArgSet(sample),combData))
    simPdf.plotOn(frame1, ROOT.RooFit.Slice(sample,"all"), ROOT.RooFit.Components("backgroundAll"), ROOT.RooFit.ProjWData(RooArgSet(sample),combData), ROOT.RooFit.LineStyle(ROOT.kDashed))

    # The same plot for the control sample slice
    frame2 = mass.frame(ROOT.RooFit.Bins(30),ROOT.RooFit.Title("Passing events"))
    combData.plotOn(frame2,ROOT.RooFit.Cut("sample==sample::pass"))
    simPdf.plotOn(frame2,ROOT.RooFit.Slice(sample,"pass"),ROOT.RooFit.ProjWData(RooArgSet(sample),combData))
    simPdf.plotOn(frame2,ROOT.RooFit.Slice(sample,"pass"),ROOT.RooFit.Components("backgroundPass"),ROOT.RooFit.ProjWData(RooArgSet(sample),combData),ROOT.RooFit.LineStyle(ROOT.kDashed))

    canvas2.cd(find_position(ptBin1, ptBin2, ptBinsX)+1)
    frame1.GetYaxis().SetTitleOffset(1.4)
    frame1.Draw()
    frame1.SetName(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame1.SaveAs(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")

    canvas3.cd(find_position(ptBin1, ptBin2, ptBinsX)+1)
    frame2.GetYaxis().SetTitleOffset(1.4)
    frame2.Draw()
    frame2.SetName(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame2.SaveAs(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")


# After filling the datasets, build the combined dataset
combDataMap = {}
frMap = {}

from array import array
hEff = ROOT.TH2D("hEff", "hEff", len(ptBinsX)-1, array('d',ptBinsX), len(ptBinsY)-1, array('d',ptBinsY))
hEff1D = ROOT.TH1D("hEff1D", "hEff1D", len(ptBinsX)-1, array('d',ptBinsX))

# Canvases for fit results
canvas2 = ROOT.TCanvas("RooFitCanvas", "RooFitCanvas", 800, 800)
canvas2.Divide(len(ptBinsX),len(ptBinsY))
canvas3 = ROOT.TCanvas("RooFitCanvasPass", "RooFitCanvasPass", 800, 800)
canvas3.Divide(len(ptBinsX),len(ptBinsY))

# Construct combined dataset in (x,sample) and perform simultaneous fit
# Skip the last, overflow, bin from fitting to save time. It does not appear in the final plots.
# Note that the following code assumes an extra bin to build the name of the output file.
for ptBin1 in range(0, len(ptBinsX)-1):
    for ptBin2 in range(0, len(ptBinsY)-1):
        combDataMap[(ptBin1, ptBin2)] = RooDataSet("combData"+"_"+str(ptBin1)+"_"+str(ptBin2),"combined data "+str(ptBin1)+"_"+str(ptBin2),
                                                   RooArgSet(mass),ROOT.RooFit.Index(sample),
                                                   ROOT.RooFit.Import("all",datasetAllMap[(ptBin1, ptBin2)]),
                                                   ROOT.RooFit.Import("pass",datasetPassMap[(ptBin1, ptBin2)]))
        #meanB.setRange(ptBinsX[ptBin1]+ptBinsY[ptBin2], ptBinsX[ptBin1+1]+ptBinsY[ptBin2+1])
        #meanB.setVal((ptBinsX[ptBin1]+ptBinsY[ptBin2]+ptBinsX[ptBin1+1]+ptBinsY[ptBin2+1])/2)
        frMap[(ptBin1, ptBin2)] = simPdf.fitTo(combDataMap[(ptBin1,ptBin2)], ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.NumCPU(4))
        # Use this for minos (better error estimate, but takes longer)
        frMap[(ptBin1, ptBin2)].Print("v")
        simPdf.getParameters(combDataMap[(ptBin1, ptBin2)]).writeToFile(buildNamePars("parameters_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".txt")
        hEff.SetBinContent(ptBin1+1, ptBin2+1, efficiency.getVal())
        hEff.SetBinError(ptBin1+1, ptBin2+1, efficiency.getError())
        hEff1D.SetBinContent(ptBin1+1, efficiency.getVal())
        hEff1D.SetBinError(ptBin1+1, efficiency.getError())
        plotResults(ptBin1, ptBin2, combDataMap[(ptBin1,ptBin2)], canvas2, canvas3)

canvas2.Print("RooFitCanvas.pdf")
canvas3.Print("RooFitCanvasPass.pdf")

canvas4 = ROOT.TCanvas("efficiency", "efficiency", 800, 800)
hEff.Draw("COLZTEXTE")
canvas4.Print("Efficiency.pdf")
canvas5 = ROOT.TCanvas("efficiency1D", "efficiency1D", 600, 600)
hEff1D.Draw()
canvas5.Print("Efficiency1D.pdf")








