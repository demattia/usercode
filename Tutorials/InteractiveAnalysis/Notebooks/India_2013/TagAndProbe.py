# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Tag and Probe

# <markdowncell>

# This exercise will show how to do an efficiency measurement using the tag and probe method.
# 
# The data has been processed and saved in tree for easier use. The code of the simple tree producer is available:
# 
# https://github.com/demattia/usercode/tree/master/DisplacedVertex/Analysis/TagAndProbe
# 
# The tree contains globalMuons (called "muons"), generalTracks (called "tracks"), triggerNames and a map of triggerFilterNames <-> triggerObjects.

# <codecell>

import ROOT
# from utils import *
# import rootnotes
from ROOT import RooRealVar, RooFormulaVar, RooVoigtian, RooChebychev, RooArgList, RooArgSet, \
    RooAddPdf, RooDataSet, RooCategory, RooSimultaneous, RooGenericPdf, RooGaussian, RooWorkspace
import bisect

# <markdowncell>

# This loads the definition of the classes for the objects stored in the tree.

# <codecell>

ROOT.gROOT.LoadMacro("Loader.C+")

# <markdowncell>

# Some utility functions and a class to pass parameters to functions

# <codecell>

class Properties:
    """Stores the min and max mass and the binning in pt"""
    def __init__(self, minMass, maxMass, ptBinsX, ptBinsY, triggerMatchDeltaR, NoBkgd, minDeltaR):
        self.minMass = minMass
        self.maxMass = maxMass
        self.ptBinsX = ptBinsX
        self.ptBinsY = ptBinsY
        self.triggerMatchDeltaR = triggerMatchDeltaR
        self.NoBkgd = NoBkgd
        self.minDeltaR = minDeltaR


muMass = 0.1057
        
def computeCosineAndMass(mu1, mu2):
    muon1 = ROOT.TLorentzVector()
    muon2 = ROOT.TLorentzVector()
    muon1.SetPtEtaPhiM(mu1.pt, mu1.eta, mu1.phi, muMass)
    muon2.SetPtEtaPhiM(mu2.pt, mu2.eta, mu2.phi, muMass)
    cosine = math.cos(muon1.Angle(muon2.Vect()))
    mass = (muon1+muon2).M()
    return (cosine, mass)

def passDileptonSelection(track1, track2, cosine):
    if track1.charge == track2.charge: return False
    if deltaR(track1.phi, track1.eta, track2.phi, track2.eta) < 0.2: return False
    if cosine <= -0.79: return False
    return True

def fillSingleCandidate(mass, p, track1, track2, histoMap, datasetMap):
    cosineAndMass = computeCosineAndMass(track1, track2)
    if not passDileptonSelection(track1, track2, cosineAndMass[0]): return False
    if cosineAndMass[1] < p.minMass or cosineAndMass[1] > p.maxMass: return False
    mass.setVal(cosineAndMass[1])
    bins = find_bins(track1.pt, track2.pt, p.ptBinsX, p.ptBinsY)
    if bins[0] == -1 or bins[1] == -1: return False
    if track1.charge == track2.charge: return False
    histoMap[bins].Fill(mass.getVal())
    #print bins, cosineAndMass[1], track1.pt, track2.pt
    datasetMap[bins].add(RooArgSet(mass))
    return True
        
def fillCandidates_tnp(mass, properties, matchedTags, probes, histoMap, datasetMap):
    for tag in matchedTags:
        for probe in probes:
            if deltaR(tag.phi,tag.eta, probe.phi, probe.eta) > properties.minDeltaR :
                fillSingleCandidate(mass, properties, tag, probe, histoMap, datasetMap)

def find_bins(pt1, pt2, ptBinsX, ptBinsY):
    return (bisect.bisect_right(ptBinsX, pt1)-1, bisect.bisect_right(ptBinsY, pt2)-1)

def find_position(bin1, bin2, ptBinsX):
    return bin1+len(ptBinsX)*bin2

def find_position_NoOverflow(bin1, bin2, ptBinsX):
    return bin1+(len(ptBinsX)-1)*bin2

def passSelectionGlobalMuon(track):
    if ( track.pt > 20 and abs(track.eta) < 2) and track.isolation/track.pt < 0.1 \
    and track.trackerLayersWithMeasurement >= 6  and track.dxy < 30. and track.dz < 30.:
        return True

def passSelection(track):
    if ( track.pt > 26 and abs(track.eta) < 2) and track.isolation/track.pt < 0.1 \
    and track.trackerLayersWithMeasurement >= 6 and track.dxy < 30. and track.dz < 30 and track.trackQuality:
        return True

def fillTriggerMatchedGlobalMuon(track, triggerObjects, matchedTracks, p):
    for triggerMuon in triggerObjects:
        if (deltaR(triggerMuon.phi, triggerMuon.eta, track.phi, track.eta) < p.triggerMatchDeltaR) and passSelectionGlobalMuon(track):
            matchedTracks.append(track)

def fillTriggerMatchedTrack(track, triggerObjects, matchedTracks, p):
    for triggerMuon in triggerObjects:
        if (deltaR(triggerMuon.phi, triggerMuon.eta, track.phi, track.eta) < p.triggerMatchDeltaR) and passSelection(track):
            matchedTracks.append(track)

def buildName(baseName, ptBin1, ptBin2, ptBinsX, ptBinsY):
    return baseName+str(ptBinsX[ptBin1])+"_"+str(ptBinsY[ptBin2])

def buildNamePars(baseName, ptBin1, ptBin12, ptBin2, ptBin22, ptBinsX, ptBinsY):
    return baseName+str(ptBinsX[ptBin1])+"_"+str(ptBinsX[ptBin12])+"_"+str(ptBinsY[ptBin2])+"_"+str(ptBinsY[ptBin22])

# <markdowncell>

# Function to build the workspace containing the pdfs to be used for the simultaneous fit of all probes and probes passing the trigger.

# <codecell>

def buildPdf(ws, p):

    mass = RooRealVar("mass", "mass", p.minMass, p.maxMass)
    getattr(ws,'import')(mass)

    # Construct signal pdf
    mean = RooRealVar("mean", "mean", 90, 85, 95)
    width = RooRealVar("width", "width", 2.4952, 1, 3)
    width.setConstant(ROOT.kTRUE)
    sigma = RooRealVar("sigma", "sigma", 1.2, 0.2, 10)
    signalAll = RooVoigtian("signalAll", "signalAll", mass, mean, width, sigma)

    turnOnAll = RooRealVar("turnOnAll","turnOnAll", 80., 40., 150.)
    widthAll_bkg = RooRealVar("widthAll","widthAll", 2., 0., 50.)
    decayAll_bkg = RooRealVar("decayAll","decayAll", 80., 20., 150.)
    meanB = RooRealVar("meanB", "meanB", 90, 60, 130)
    sigmaB = RooRealVar("sigmaB", "sigmaB", 10, 1, 20)
    bkg_a1 = RooRealVar("bkg_a1", "bkg_a1", 0., -2., 2.)
    bkg_a2 = RooRealVar("bkg_a2", "bkg_a2", 0., -2., 2.)
    backgroundAll = RooGaussian("backgroundAll", "backgroundAll", mass, meanB, sigmaB)

    # Construct composite pdf
    sigAll = RooRealVar("sigAll", "sigAll", 2000, 0, 100000)
    bkgAll = RooRealVar("bkgAll", "bkgAll", 100, 0, 10000)
    modelAll = RooAddPdf("modelAll", "modelAll", RooArgList(signalAll, backgroundAll), RooArgList(sigAll, bkgAll))
    if p.NoBkgd:
        modelAll = RooAddPdf("modelAll", "modelAll", RooArgList(signalAll), RooArgList(sigAll))
    # Define pdf for all probes

    # Construct signal pdf.
    # NOTE that sigma is shared with the signal sample model
    signalPass = RooVoigtian("signalPass","signalPass",mass,mean,width,sigma)
    # Construct the background pdf
    backgroundPass = RooGaussian("backgroundPass", "backgroundPass", mass, meanB, sigmaB)

    # Construct the composite model
    efficiency = RooRealVar("efficiency","efficiency",0.9,0.3,1.)
    sigPass = RooFormulaVar("sigPass", "@0*@1", RooArgList(sigAll, efficiency))
    bkgPass = RooRealVar("bkgPass", "bkgPass", 100, 0, 10000)
    modelPass = RooAddPdf("modelPass", "modelPass", RooArgList(signalPass, backgroundPass), RooArgList(sigPass, bkgPass))
    if p.NoBkgd:
        modelPass = RooAddPdf("modelPass", "modelPass", RooArgList(signalPass), RooArgList(sigPass))

    frac = RooRealVar("frac", "frac", 0.8, 0., 1.)

    # Define combined pdf for simultaneous fit

    # Define category to distinguish physics and control samples events
    sample = RooCategory("sample","sample")
    sample.defineType("all")
    sample.defineType("pass")

    simPdf = RooSimultaneous("simPdf","simultaneous pdf",sample)

    # Associate model with the physics state and model_ctl with the control state
    simPdf.addPdf(modelAll,"all")
    simPdf.addPdf(modelPass,"pass")
    # ws.import(simPdf)
    getattr(ws,'import')(simPdf)

# <markdowncell>

# Define the triggers to use and load the input tree

# <codecell>

# Triggers for tag and for probe
tagTrigger = "IsoMu24_v"
probeTrigger = "HLT_L2DoubleMu23_NoVertex_v"

# Define cuts and some useful variables
triggerMatchDeltaR = 0.1
minMass = 80
maxMass = 100
minDeltaR = 0.2

# --------#

MC = True
NoBkgd = True

# --------#

# <headingcell level=3>

# Event loop to fill histograms and datasets

# <codecell>

# %%rootprint

ROOT.gROOT.LoadMacro("Loader.C+")

# Load the input file
tree = ROOT.TChain("T")
tree.Add("Z_mumu_MC_tag_and_probe.root")

# tree.Add("TagAndProbe_ZMuMu.root")



# Define binning
ptBinsTag = [20, 10000]
ptBinsProbe = [26, 30, 35, 40, 45, 50, 60, 70]



p = Properties(minMass, maxMass, ptBinsTag, ptBinsProbe, triggerMatchDeltaR, NoBkgd, minDeltaR)

# Prepare the workspace
ws = RooWorkspace("ws", "workspace")
buildPdf(ws, p)

mass = ws.var("mass")
sample = ws.cat("sample")
simPdf = ws.pdf("simPdf")
efficiency = ws.var("efficiency")
meanB = ws.var("meanB")

# Prepare datasets and histograms
datasetAllMap = {}
datasetPassMap = {}
hAllMap = {}
hPassMap = {}
for ptBin1 in range(0, len(ptBinsTag)):
    for ptBin2 in range(0, len(ptBinsProbe)):
        datasetAllMap[(ptBin1, ptBin2)] = RooDataSet(buildName("datasetAll_", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                                     buildName("datasetAll_", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                                     RooArgSet(mass))
        datasetPassMap[(ptBin1, ptBin2)] = RooDataSet(buildName("datasetPass_", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                                      buildName("datasetPass_", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                                      RooArgSet(mass))
        hAllMap[(ptBin1,ptBin2)] = ROOT.TH1F(buildName("hAll_", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                             buildName("All probes", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                             100, minMass, maxMass)
        hPassMap[(ptBin1,ptBin2)] = ROOT.TH1F(buildName("hPass_", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                              buildName("Passing probes", ptBin1, ptBin2, ptBinsTag, ptBinsProbe),
                                              100, minMass, maxMass)



# Event loop
allCandidates = 0
passCandidates = 0

processedEvents = 0

totEvents = 100000
#totEvents = tree.GetEntries()
progress = 0


for event in tree:
    if processedEvents > totEvents:
        break
    processedEvents += 1

    # Get the trigger outcome for this event
    tagTriggerFired = False
    probeTriggerFired = False
    for name in event.triggerNames:
        if name.find(tagTrigger) != -1: tagTriggerFired = True
    # If none of the two triggers fired we can skip the event
    if not tagTriggerFired:
        continue

    # IsoMu24 filter name
    tagTriggerObjects = event.triggerFilterObjectsMap["hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15"]
    # DoubelMu23NoVertex filter name
    probeTriggerObjects = event.triggerFilterObjectsMap["hltL2DoubleMu23NoVertexL2PreFiltered"]
    # tagTriggerObjects = []
    # probeTriggerObjects = []

    matchedMuonsTagTrigger = []
    passingMuonsProbe = []
    allMuonsProbe = []
    
    # Find tag muons
    for tagMuon in event.muons:
        # Find a matching trigger object in DeltaR
        fillTriggerMatchedGlobalMuon(tagMuon, tagTriggerObjects, matchedMuonsTagTrigger, p)

    for probeMuon in event.tracks:
        # Find a matching trigger object in DeltaR
        if passSelection(probeMuon):
            fillTriggerMatchedTrack(probeMuon, probeTriggerObjects, passingMuonsProbe, p)
            allMuonsProbe.append(probeMuon)

    # Fill datasets and histograms
    fillCandidates_tnp(mass, p, matchedMuonsTagTrigger, passingMuonsProbe, hPassMap, datasetPassMap)
    fillCandidates_tnp(mass, p, matchedMuonsTagTrigger, allMuonsProbe, hAllMap, datasetAllMap)

# <headingcell level=3>

# Simultaneous fits to extract the efficiencies

# <markdowncell>

# First, some control plots before the fits

# <codecell>

print "Note: you can expand the figure dragging the border on the lower right\n"

# canvas = rootnotes.canvas("AllAndPassCanvas", (1000, 800))
# # canvas.Divide(len(ptBinsTag)-1, len(ptBinsProbe))
# canvas.Divide(4,2)
# for ptBin1 in range(0, len(ptBinsTag)-1):
#     for ptBin2 in range(0, len(ptBinsProbe)):
#         canvas.cd(find_position_NoOverflow(ptBin1, ptBin2, ptBinsTag)+1)
#         hAllMap[(ptBin1, ptBin2)].Draw()
#         hPassMap[(ptBin1, ptBin2)].Draw("same")
#         hPassMap[(ptBin1, ptBin2)].SetLineColor(2)

# canvas.Print("AllAndPassCanvas.pdf")
# canvas

# <markdowncell>

# Utility function to save the fit results to local files

# <codecell>

# Plot the fit results
def simultaneousFit(ptBin1, ptBin2, combData):
    #, canvas2, canvas3):
    frame1 = mass.frame(ROOT.RooFit.Bins(30),ROOT.RooFit.Title("All events"))
    # Plot all data tagged as physics sample
    combData.plotOn(frame1,ROOT.RooFit.Cut("sample==sample::all"))
    # Plot "physics" slice of simultaneous pdf.
    # NBL You _must_ project the sample index category with data using ProjWData
    # as a RooSimultaneous makes no prediction on the shape in the index category
    # and can thus not be integrated
    simPdf.plotOn(frame1, ROOT.RooFit.Slice(sample,"all"), ROOT.RooFit.ProjWData(RooArgSet(sample),combData))
    simPdf.plotOn(frame1, ROOT.RooFit.Slice(sample,"all"), ROOT.RooFit.Components("backgroundAll"),
                  ROOT.RooFit.ProjWData(RooArgSet(sample),combData), ROOT.RooFit.LineStyle(ROOT.kDashed))

    # The same plot for the control sample slice
    frame2 = mass.frame(ROOT.RooFit.Bins(30),ROOT.RooFit.Title("Passing events"))
    combData.plotOn(frame2,ROOT.RooFit.Cut("sample==sample::pass"))
    simPdf.plotOn(frame2,ROOT.RooFit.Slice(sample,"pass"),ROOT.RooFit.ProjWData(RooArgSet(sample),combData))
    simPdf.plotOn(frame2,ROOT.RooFit.Slice(sample,"pass"),ROOT.RooFit.Components("backgroundPass"),
                  ROOT.RooFit.ProjWData(RooArgSet(sample),combData),ROOT.RooFit.LineStyle(ROOT.kDashed))

    frame1.GetYaxis().SetTitleOffset(1.4)
    frame1.Draw()
    frame1.SetName(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe))
    frame1.SaveAs(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe)+".root")

    frame2.GetYaxis().SetTitleOffset(1.4)
    frame2.Draw()
    frame2.SetName(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe))
    frame2.SaveAs(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe)+".root")

# <markdowncell>

# Perform the simultaneous fit to all probes and passing probes to extract the efficiency and fill the plots using the function above

# <codecell>

# Build the combined dataset
combDataMap = {}
frMap = {}

from array import array
hEff = ROOT.TH2D("hEff", "hEff", len(ptBinsTag)-1, array('d',ptBinsTag), len(ptBinsProbe)-1, array('d',ptBinsProbe))
hEff1D = ROOT.TH1D("hEff1D", "hEff1D", len(ptBinsTag)-1, array('d',ptBinsTag))

# Construct combined dataset in (x,sample) and perform simultaneous fit
# Skip the last, overflow, bin from fitting to save time. It does not appear in the final plots.
# Note that the following code assumes an extra bin to build the name of the output file.
for ptBin1 in range(0, len(ptBinsTag)-1):
    for ptBin2 in range(0, len(ptBinsProbe)-1):
        combDataMap[(ptBin1, ptBin2)] = RooDataSet("combData"+"_"+str(ptBin1)+"_"+str(ptBin2),
                                                   "combined data "+str(ptBin1)+"_"+str(ptBin2),
                                                   RooArgSet(mass),ROOT.RooFit.Index(sample),
                                                   ROOT.RooFit.Import("all",datasetAllMap[(ptBin1, ptBin2)]),
                                                   ROOT.RooFit.Import("pass",datasetPassMap[(ptBin1, ptBin2)]))
        # meanB.setRange(ptBinsTag[ptBin1]+45, ptBinsTag[ptBin1+1]+45)
        # meanB.setVal((ptBinsTag[ptBin1]+ptBinsTag[ptBin1+1])/2+45)
        frMap[(ptBin1, ptBin2)] = simPdf.fitTo(combDataMap[(ptBin1,ptBin2)], ROOT.RooFit.Save(ROOT.kTRUE),
                                               ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.NumCPU(4))
        # Use this for minos (better error estimate, but takes longer)
        frMap[(ptBin1, ptBin2)].Print("v")
        simPdf.getParameters(combDataMap[(ptBin1, ptBin2)]).writeToFile(buildNamePars("parameters_", ptBin1, ptBin1+1,
                                                                                      ptBin2, ptBin2+1, ptBinsTag,
                                                                                      ptBinsProbe)+".txt")
        hEff.SetBinContent(ptBin1+1, ptBin2+1, efficiency.getVal())
        hEff.SetBinError(ptBin1+1, ptBin2+1, efficiency.getError())
        hEff1D.SetBinContent(ptBin1+1, efficiency.getVal())
        hEff1D.SetBinError(ptBin1+1, efficiency.getError())
        simultaneousFit(ptBin1, ptBin2, combDataMap[(ptBin1,ptBin2)])

# <headingcell level=3>

# Make all plots

# <codecell>

import ROOT
from utils import *

# Plot the fit results

# def plotResults(ptBin1, ptBin2, canvas2, canvas3):
#     canvas2.cd(find_position_NoOverflow(ptBin1, ptBin2, ptBinsTag)+1)
#     inputAll = ROOT.TFile.Open(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe)+".root")
#     frame1 = inputAll.FindObjectAny(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe))
#     frame1.Draw()
# 
#     canvas3.cd(find_position_NoOverflow(ptBin1, ptBin2, ptBinsTag)+1)
#     inputPass= ROOT.TFile.Open(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe)+".root")
#     frame2 = inputPass.FindObjectAny(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe))
#     frame2.Draw()
# 
# 
# from array import array
# hEff = ROOT.TH2D("hEff", "hEff", len(ptBinsTag)-1, array('d',ptBinsTag), len(ptBinsProbe)-1, array('d',ptBinsProbe))
# hEff1D = ROOT.TH1D("hEff1D", "hEff1D", len(ptBinsProbe)-1, array('d',ptBinsProbe))
# 
# # Canvases for fit results
# canvas2 = rootnotes.canvas("RooFitCanvas", (400, 2000))
# canvas2.Divide(len(ptBinsTag)-1,len(ptBinsProbe)-1)
# canvas3 = rootnotes.canvas("RooFitCanvasPass", (400, 2000))
# canvas3.Divide(len(ptBinsTag)-1,len(ptBinsProbe)-1)
# 
# # Construct combined dataset in (x,sample) and perform simultaneous fit
# # Skip the last, overflow, bin from fitting to save time. It does not appear in the final plots.
# # Note that the following code assumes an extra bin to build the name of the output file.
# for ptBin1 in range(0, len(ptBinsTag)-1):
#     for ptBin2 in range(0, len(ptBinsProbe)-1):
#         plotResults(ptBin1, ptBin2, canvas2, canvas3)
#         for line in open(buildNamePars("parameters_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsTag, ptBinsProbe)+".txt"):
#             if line.find("efficiency") != -1:
#                 eff = float(line.split()[2])
#                 effErr = float(line.split()[4])
#                 print "eff["+str(ptBinsTag[ptBin1])+"_"+str(ptBinsTag[ptBin1+1])+\
#                 ", "+str(ptBinsProbe[ptBin2])+"_"+str(ptBinsProbe[ptBin2+1])+"] = ", eff, "+/-", effErr
#                 hEff.SetBinContent(ptBin1+1, ptBin2+1, eff)
#                 hEff.SetBinError(ptBin1+1, ptBin2+1, effErr)
#                 hEff1D.SetBinContent(ptBin2+1, eff)
#                 hEff1D.SetBinError(ptBin2+1, effErr)
# 
# canvas4 = rootnotes.canvas("efficiency", (600, 600))
# hEff.Draw("COLZTEXTE")
# hEff.SetStats(ROOT.kFALSE)
# 
# canvas5 = rootnotes.canvas("efficiency1D", (600, 600))
# hEff1D.GetXaxis().SetTitle("probe muon p_{T} [GeV/c]")
# hEff1D.GetYaxis().SetTitle("passing/all")
# hEff1D.Draw()
# hEff1D.GetYaxis().SetRangeUser(0,1.1)
# hEff1D.SetStats(ROOT.kFALSE)
# 
# # <codecell>
# 
# canvas2
# 
# # <codecell>
# 
# canvas3
# 
# # <codecell>
# 
# canvas4
# 
# # <codecell>
# 
# canvas5

# <headingcell level=3>

# Questions

# <markdowncell>

# **Question 1**: Redo the fit adding the background. What is the effect on the results?

# <codecell>


# <markdowncell>

# **Question 2**: add binning in $\eta$ and compute the efficiency vs $p_T$ and $\eta$

# <codecell>


