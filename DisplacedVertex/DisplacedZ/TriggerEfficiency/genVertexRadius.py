import math
import ROOT
import sys
from DataFormats.FWLite import Events, Handle

# events = Events ('reco_1.root')
# events = Events ('PromptZ/reco_1.root')

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
handle  = Handle ('std::vector<reco::GenParticle>')

# a label is just a tuple of strings that is initialized just like an edm::InputTag
label = ("genParticles")

# Disable drawing to the screen
# ROOT.gROOT.SetBatch()        # don't pop up canvases

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
# Note that the genParticle distances are in mm
vertexRadiusHist = ROOT.TH1F ("vertexRadiusHist", "Vertex radius", 500, -2, 2.5)
vertexRadiusHist.GetXaxis().SetTitle("Vertex radius (mm)");
ZPtHist = ROOT.TH1F ("ZPtHist", "Pt of the Z", 500, 0, 200)
ZPtHist.GetXaxis().SetTitle("p_{T} [GeV/c]");
# ZMassHist = ROOT.TH1F ("ZMass", "generated Z mass", 100, 70, 120)
# ZMassHist.GetXaxis().SetTitle("Mass [GeV/c^{2}]");

muonPtLowHist = ROOT.TH1F ("MuonPtLow", "lowest pt of the two muons", 500, 0, 200)
muonPtLowHist.GetXaxis().SetTitle("p_{T} [GeV/c]");
muonPtHighHist = ROOT.TH1F ("MuonPtHigh", "highest pt of the two muons", 500, 0, 200)
muonPtHighHist.GetXaxis().SetTitle("p_{T} [GeV/c]");

# loop over events
for event in events:
    # use getByLabel, just like in cmsRun
    event.getByLabel(label, handle)

    # get the product
    genParticles = handle.product()

    muonsPt = []

    for genParticle in genParticles:
        if genParticle.pdgId() == 23:
            ZPtHist.Fill( genParticle.pt() )
            # ZMassHist.Fill( genParticle.mass() )
        if abs(genParticle.pdgId()) == 11 or abs(genParticle.pdgId()) == 13 or abs(genParticle.pdgId()) == 15:
            muonsPt.append(genParticle.pt())
            if genParticle.mother().pdgId() == 23:
                vertexRadiusHist.Fill( math.sqrt(genParticle.vx()**2 + genParticle.vy()**2) )


    if len(muonsPt) < 2:
        print "Error, less than 2 stable muons found"
        sys.exit()
    muonPtLowHist.Fill( min(muonsPt[0], muonsPt[1]) )
    muonPtHighHist.Fill( max(muonsPt[0], muonsPt[1]) )


                # inner4v = ROOT.TLorentzVector (innerMuon.px(), innerMuon.py(),
                #                                innerMuon.pz(), innerMuon.energy())
                # outer4v = ROOT.TLorentzVector (outerMuon.px(), outerMuon.py(),
                #                                outerMuon.pz(), outerMuon.energy())

# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
c1.Divide(2,2)
c1.cd(1)
vertexRadiusHist.Draw()
c1.cd(2)
ZPtHist.Draw()
c1.cd(3)
# ZMassHist.Draw()
# c1.cd(4)
muonPtLowHist.Draw()
c1.cd(4)
muonPtHighHist.Draw()
c1.Print("Z.png")
