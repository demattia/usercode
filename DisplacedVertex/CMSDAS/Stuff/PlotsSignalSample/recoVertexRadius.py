import os
import sys
import math
import ROOT
from DataFormats.FWLite import Events, Handle

# Reminder: the generation is H0 -> A0 A0, with the A0 forced to decay to leptons. The H0 is 35 and the A0 is 36.
# create handle outside of loop
nfile = int(sys.argv[1])
print nfile
file = open("file.txt")
#print (list(file)[1])
#filename = list(file)[nfile]	
filename = 'root://xrootd.rcac.purdue.edu/%s' % (list(file)[nfile])
filename = filename.rstrip('\n')
print filename
 
sample = open("sample.txt")
samplename = list(sample)[nfile]

#events = Events('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_200_050F/pat/PATtuple_1_1_ZcD.root')
events = Events (filename)
handleMuon  = Handle ('std::vector<pat::Muon>')
handle  = Handle ('std::vector<reco::GenParticle>')
labelMuon = ("selectedPatMuons")
label = ("genParticles")

# a label is just a tuple of strings that is initialized just like an edm::InputTag
#labelTrack = ("generalTracks")
#labelStandAlone = ("standAloneMuons")
##labelMuon = ("muons")

# Disable drawing to the screen
# ROOT.gROOT.SetBatch()        # don't pop up canvases

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
# Note that the genParticle distances are in mm
transIPMuonHist = ROOT.TH1F ("transIPMuonHist", "Transverse Impact Parameter", 500, 0, 200)
transIPMuonHist.GetXaxis().SetTitle("Transverse IP(cm)");
ptMuHist = ROOT.TH1F ("ptMuHist", "Global Muon p_{T}", 500, 0, 500)
ptMuHist.GetXaxis().SetTitle("p_{T}(GeV/c)")
ptGenMuHist = ROOT.TH1F ("ptGenMuHist", "Gen Muon p_{T}", 500, 0, 500)
ptGenMuHist.GetXaxis().SetTitle("p_{T}(GeV/c)")
VertexGenMuHist = ROOT.TH1F ("VertexGenMuHist", "Vertex Radius", 500, 0, 200)
VertexGenMuHist.GetXaxis().SetTitle("Vertex Radius(cm)");
ptGenMuCutHist = ROOT.TH1F ("ptGenMuCutHist", "Gen Muon p_{T} > 25", 500, 0, 500)
ptGenMuCutHist.GetXaxis().SetTitle("p_{T}(GeV/c)")
ptGenMuFromH0 = ROOT.TH1F ("ptGenMuFromH0", "Gen Muon p_{T} from H0", 500, 0, 500)
ptGenMuFromH0.GetXaxis().SetTitle("p_{T}(GeV/c)")
VertexGenMuPtCutHist = ROOT.TH1F ("VertexGenMuPtCutHist", "Vertex Radius", 500, 0, 200)
VertexGenMuPtCutHist.GetXaxis().SetTitle("Vertex Radius(cm)");
VertexGenMuFromH0 = ROOT.TH1F ("VertexGenMuFromH0", "Vertex Radius", 500, 0, 200)
VertexGenMuFromH0.GetXaxis().SetTitle("Vertex Radius(cm)");
transIPMuonPtCutHist = ROOT.TH1F ("transIPMuonPtCutHist", "Transverse Impact Parameter", 500, 0, 200)
transIPMuonPtCutHist.GetXaxis().SetTitle("Transverse IP(cm)");
ptMuCutHist = ROOT.TH1F ("ptMuCutHist", "Global Muon p_{T} > 25", 500, 0, 500)
ptMuCutHist.GetXaxis().SetTitle("p_{T}(GeV/c)")

# loop over events
for event in events:
    # use getByLabel, just like in cmsRun
#    event.getByLabel(labelTrack, handleTrack)
#    event.getByLabel(labelStandAlone, handleStandAlone)
    event.getByLabel(labelMuon, handleMuon)
    
    # get the product
#    tracks = handleTrack.product()
#    standAloneMuons = handleStandAlone.product()
    muons = handleMuon.product()

    event.getByLabel(label, handle)

    # get the product
    genParticles = handle.product()

    for genParticle in genParticles:
	if abs(genParticle.pdgId()) == 13 and genParticle.status()==1 and genParticle.pt() > 3:
		ptGenMuHist.Fill(genParticle.pt())
		VertexGenMuHist.Fill( math.sqrt(genParticle.vx()**2 + genParticle.vy()**2) )
		if genParticle.pt()>25:
			ptGenMuCutHist.Fill(genParticle.pt())
			VertexGenMuPtCutHist.Fill( math.sqrt(genParticle.vx()**2 + genParticle.vy()**2) )
#		if genParticle.mother().mother().pdgId() == 36 and genParticle.mother().pdgId() == genParticle.pdgId():
#			ptGenMuFromH0.Fill(genParticle.pt())
#			VertexGenMuFromH0.Fill(math.sqrt(genParticle.vx()**2 + genParticle.vy()**2) )
#    for track in tracks:
#        vertexRadiusTrackHist.Fill( math.sqrt(track.vx()**2 + track.vy()**2) )
#    for standAloneMuon in standAloneMuons:
#        vertexRadiusStandAloneHist.Fill( math.sqrt(standAloneMuon.vx()**2 + standAloneMuon.vy()**2) )
    for muon in muons:
	if muon.isGlobalMuon()==1:
        	transIPMuonHist.Fill( math.sqrt(muon.vx()**2 + muon.vy()**2) )
        	ptMuHist.Fill( muon.pt())
		if muon.pt() > 25:
			transIPMuonPtCutHist.Fill( math.sqrt(muon.vx()**2 + muon.vy()**2) )
			ptMuCutHist.Fill( muon.pt())
# make a canvas, draw, and save it
samplename = samplename.rstrip('\n')
print samplename
c1 = ROOT.TCanvas()
ptGenMuHist.Draw()
c1.Print ('GenPtMuon_Signal_%s.png' % (samplename))
 
c2 = ROOT.TCanvas()
VertexGenMuHist.Draw()
c2.SetLogy()
c2.Print ('VertexGenMu_Signal_%s.png' % (samplename))

c3 = ROOT.TCanvas()
transIPMuonHist.Draw()
c3.SetLogy()
c3.Print ('transIPMuon_Signal_%s_glb.png' % (samplename))

c4 = ROOT.TCanvas()
ptMuHist.Draw()
c4.SetLogy()
c4.Print ('PtMuon_Signal_%s_glb.png' % (samplename))

c5 = ROOT.TCanvas()
ptGenMuCutHist.Draw()
c5.Print ('PtGenMuon_Pt25_Signal_%s.png' % (samplename))

c6 = ROOT.TCanvas()
VertexGenMuPtCutHist.Draw()
c6.SetLogy()
c6.Print ('VertexGenMu_Pt25_Signal_%s.png' % (samplename))

c7 = ROOT.TCanvas()
transIPMuonPtCutHist.Draw()
c7.SetLogy()
c7.Print ('transIPMuon_Pt25_Signal_%s_glb.png' % (samplename))

c8 = ROOT.TCanvas()
ptMuCutHist.Draw()
c8.Print ('PtMuon_Pt25_Signal_%s_glb.png' % (samplename))


