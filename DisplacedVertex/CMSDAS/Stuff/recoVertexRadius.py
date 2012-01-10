import math
import ROOT
from DataFormats.FWLite import Events, Handle

# Reminder: the generation is H0 -> A0 A0, with the A0 forced to decay to leptons. The H0 is 35 and the A0 is 36.
# create handle outside of loop
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_200_050F/pat/PATtuple_1_1_ZcD.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7//Signal_400_005L/pat/PATtuple_6_1_i1B.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_400_020F/pat/PATtuple_4_1_BWV.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_120_020F//pat/PATtuple_5_1_nlV.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_1000_150F//pat/PATtuple_1_2_2Fr.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_1000_020F//pat/PATtuple_9_1_bRA.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_1000_050F/pat/PATtuple_10_1_HXU.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_120_050F/pat/PATtuple_3_2_cUJ.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_1000_350F/pat/PATtuple_3_1_9mQ.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_200_020F/pat/PATtuple_2_1_T2M.root')
#events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_400_150F/pat/PATtuple_5_1_uR7.root')
events = Events ('root://xrootd.rcac.purdue.edu//store/user/demattia/longlived/CMSSW_4_2_7/Signal_400_050F/pat/PATtuple_9_1_nd6.root')
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
transIPTrackHist = ROOT.TH1F ("transIPTrackHist", "Transverse Impact Parameter", 500, 0, 200)
transIPTrackHist.GetXaxis().SetTitle("Transverse IP(cm)");
transIPStandAloneHist = ROOT.TH1F ("vertexRadiusStandAloneHist", "Transverse Impact Parameter", 500, 0, 200)
transIPStandAloneHist.GetXaxis().SetTitle("Vertex radius (cm)");
transIPMuonHist = ROOT.TH1F ("vertexRadiusMuonHist", "Transverse Impact Parameter", 500, 0, 200)
transIPMuonHist.GetXaxis().SetTitle("Transverse IP(cm)");
ptMuHist = ROOT.TH1F ("ptMuHist", "ptMuHist", 500, 0, 500)
ptMuHist.GetXaxis().SetTitle("p_{T}(GeV/c)")
ptGenMuHist = ROOT.TH1F ("ptGenMuHist", "ptGenMuHist", 500, 0, 500)
ptGenMuHist.GetXaxis().SetTitle("p_{T}(GeV/c)")
transIPGenMuHist = ROOT.TH1F ("transIPGenMuHist", "transIPGenMuHist", 500, 0, 200)
transIPGenMuHist.GetXaxis().SetTitle("Transverse IP(cm)");

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
	if abs(genParticle.pdgId()) == 13:
		ptGenMuHist.Fill(genParticle.pt())
		transIPGenMuHist.Fill( math.sqrt(genParticle.vx()**2 + genParticle.vy()**2) )
#    for track in tracks:
#        vertexRadiusTrackHist.Fill( math.sqrt(track.vx()**2 + track.vy()**2) )
#    for standAloneMuon in standAloneMuons:
#        vertexRadiusStandAloneHist.Fill( math.sqrt(standAloneMuon.vx()**2 + standAloneMuon.vy()**2) )
    for muon in muons:
	if muon.isGlobalMuon()==1:
        	transIPMuonHist.Fill( math.sqrt(muon.vx()**2 + muon.vy()**2) )
        	ptMuHist.Fill( muon.pt())
# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
ptGenMuHist.Draw()
c1.Print ("GenPtMuon_Signal_400_050F.png")
 
c2 = ROOT.TCanvas()
transIPGenMuHist.Draw()
c2.Print ("transIPGenMu_Signal_400_050F.png")

c3 = ROOT.TCanvas()
transIPMuonHist.Draw()
c3.Print ("transIPMuon_Signal_400_050F_glb.png")
#c3.SetLogy()
c4 = ROOT.TCanvas()
ptMuHist.Draw()
#c4.SetLogy()
c4.Print ("PtMuon_Signal_400_050F_glb.png")

