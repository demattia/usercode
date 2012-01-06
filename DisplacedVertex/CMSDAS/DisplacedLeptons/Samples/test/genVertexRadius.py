import math
import ROOT
from DataFormats.FWLite import Events, Handle

# Reminder: the generation is H0 -> A0 A0, with the A0 forced to decay to leptons. The H0 is 35 and the A0 is 36.
events = Events ('Pythia_H0_pyupda_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root')

# create handle outside of loop
handle  = Handle ('std::vector<reco::GenParticle>')

# a label is just a tuple of strings that is initialized just like an edm::InputTag
label = ("genParticles")

# Disable drawing to the screen
# ROOT.gROOT.SetBatch()        # don't pop up canvases

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
# Note that the genParticle distances are in mm
vertexRadiusHist = ROOT.TH1F ("vertexRadiusHist", "Vertex radius", 500, 0, 500)
vertexRadiusHist.GetXaxis().SetTitle("Vertex radius (mm)");

# loop over events
for event in events:
    # use getByLabel, just like in cmsRun
    event.getByLabel(label, handle)
    
    # get the product
    genParticles = handle.product()

    for genParticle in genParticles:
        if abs(genParticle.pdgId()) == 11 or abs(genParticle.pdgId()) == 13 or abs(genParticle.pdgId()) == 15:
            # if genParticle.mother().pdgId() == genParticle.pdgId() and (genParticle.status() == 1 or genParticle.pdgId() == 15) and genParticle.mother().status() == 3:
            # The second request after the "and" is not needed, it is only there as a reminder
            if genParticle.mother().mother().pdgId() == 36 and genParticle.mother().pdgId() == genParticle.pdgId():
                vertexRadiusHist.Fill( math.sqrt(genParticle.vx()**2 + genParticle.vy()**2) )
                # inner4v = ROOT.TLorentzVector (innerMuon.px(), innerMuon.py(),
                #                                innerMuon.pz(), innerMuon.energy())
                # outer4v = ROOT.TLorentzVector (outerMuon.px(), outerMuon.py(),
                #                                outerMuon.pz(), outerMuon.energy())
            
# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
vertexRadiusHist.Draw()
c1.Print ("vertexRadius.png")
            
