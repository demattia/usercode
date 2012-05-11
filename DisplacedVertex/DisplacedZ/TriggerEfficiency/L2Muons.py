import math
import ROOT
import sys
from DataFormats.FWLite import Events, Handle

# events = Events ('reco_1.root')

events = Events(
    [
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_1.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_2.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_3.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_4.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_5.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_6.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_7.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_8.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_9.root',
    'castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/DisplacedZ/5_2_5/reco_10.root'
    ]
)

# create handle outside of loop
handle  = Handle ('std::vector<reco::Track>')
handleTracks  = Handle ('std::vector<reco::Track>')
handleTrigger = Handle ('edm::TriggerResults')

# a label is just a tuple of strings that is initialized just like an edm::InputTag
label = ("hltL2Muons")
labelTracks = ("generalTracks")
labelTrigger = ("TriggerResults")

# Disable drawing to the screen
# ROOT.gROOT.SetBatch()        # don't pop up canvases

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
# Note that the genParticle distances are in mm

muonPtLowHist = ROOT.TH1F ("MuonPtLow", "lowest pt of the two muons", 500, 0, 200)
muonPtLowHist.GetXaxis().SetTitle("p_{T} [GeV/c]")
muonPtHighHist = ROOT.TH1F ("MuonPtHigh", "highest pt of the two muons", 500, 0, 200)
muonPtHighHist.GetXaxis().SetTitle("p_{T} [GeV/c]")
muons3DAngleHist = ROOT.TH1F ("Muons3DAngle", "3D angle between the two muons", 100, -3.15, 3.15)
muons3DAngleHist.GetXaxis().SetTitle("3D angle")

passTrigger = 0
passOfflinePtCuts = 0
totalEvents = 0
# loop over events
for event in events:
    totalEvents += 1
    try:
        # Note that there are also triggerResults for SIM and RECO which are not needed here. Take the HLT results.
        # If not specified it will take the last one by default which is the RECO.
        event.getByLabel(labelTrigger, "", "HLT", handleTrigger)
        triggerResults = handleTrigger.product()
        # print "number of paths = ", triggerResults.size()
    except:
        print "trigger results error"
        continue
    try:
        event.getByLabel(label, handle)
        muons = handle.product()
        event.getByLabel(labelTracks, handleTracks)
        tracks = handleTracks.product()
    except:
        print "No hlt muon found"
        continue

    print "Analyzing event", totalEvents

    muonsArray = []
    for muon in muons:
        muonsArray.append(muon)

    tracksArray = []
    for track in tracks:
        tracksArray.append(track)
    tracksArray = sorted(tracksArray, key=lambda track: track.pt(), reverse=True)

    # if ptLow > 23 and ptHigh > 23 and angle < 2.5 and validChambersOne > 1 and validChambersTwo > 1:
    if triggerResults.accept(147):
        passTrigger += 1
        if len(tracksArray) > 1 and min(tracksArray[0].pt(), tracksArray[1].pt()) > 26:
            passOfflinePtCuts += 1


    if len(muonsArray) < 2:
        print "Less than 2 hlt muons found"
        continue
    ptLow = min(muonsArray[0].pt(), muonsArray[1].pt())
    ptHigh = max(muonsArray[0].pt(), muonsArray[1].pt())
    muonPtLowHist.Fill( ptLow )
    muonPtHighHist.Fill( ptHigh )
    angle = math.acos((muonsArray[0].px()*muonsArray[1].px() + muonsArray[0].py()*muonsArray[1].py() + muonsArray[0].pz()*muonsArray[1].pz())/(muonsArray[0].p()*muonsArray[1].p()));
    muons3DAngleHist.Fill(angle)

    # validChambersOne = muonsArray[0].hitPattern().dtStationsWithAnyHits() + muonsArray[0].hitPattern().cscStationsWithAnyHits()
    # validChambersTwo = muonsArray[1].hitPattern().dtStationsWithAnyHits() + muonsArray[1].hitPattern().cscStationsWithAnyHits()

triggerEfficiency = passTrigger/float(totalEvents)
print "trigger efficiency =", triggerEfficiency, "+/-", math.sqrt(triggerEfficiency*(1-triggerEfficiency)/totalEvents)
passOfflinePtCutsEfficiency = passOfflinePtCuts/float(totalEvents)
print "offline pt cuts efficiency =", passOfflinePtCutsEfficiency, "+/-", math.sqrt(passOfflinePtCutsEfficiency*(1-passOfflinePtCutsEfficiency)/totalEvents)

# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
c1.Divide(2,2)
c1.cd(1)
muonPtLowHist.Draw()

lineLow = ROOT.TLine(23,0,23,muonPtLowHist.GetMaximum());
lineLow.SetLineColor(2)
lineLow.SetLineWidth(2)
lineLow.Draw()

c1.cd(2)
muonPtHighHist.Draw()

lineHigh = ROOT.TLine(23,0,23,muonPtHighHist.GetMaximum());
lineHigh.SetLineColor(2)
lineHigh.SetLineWidth(2)
lineHigh.Draw()

c1.cd(3)
muons3DAngleHist.Draw()
lineAngle = ROOT.TLine(2.5,0,2.5,muons3DAngleHist.GetMaximum());
lineAngle.SetLineColor(2)
lineAngle.SetLineWidth(2)
lineAngle.Draw()

c1.Print("L2Muons.png")
