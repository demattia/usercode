#!/usr/bin/python

print "Importing ROOT"
from ROOT import *
print "Importing cmstools"
from PhysicsTools.PythonAnalysis.cmstools import *
#from cmstools import *

#from EventChain import *

### prepare the FWLite autoloading mechanism 
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()

# open the file and access the event tree 
# events = EventTree("/data2/demattia/Data/MuonGun/MuonGun_0.root")
eventsChain = TChain("Events")
eventsChain.Add("/data2/demattia/Data/MuonGun/MuonGun_0.root")
events = EventTree(eventsChain.GetTree())

# Create a TH2F to store the values
h = TH2F("ptVSeta","muon pt vs eta", 100, -3.2, 3.2, 100, 0, 120);
hProfile = TProfile("ptVSetaProfile","muon pt vs eta", 100, -3.2, 3.2, 0, 120);

# loop over the events
iEvent = 0
for event in events:
    muons = event.getProduct("muons")
    source = event.source
    iEvent += 1
    print "iEvent = ", iEvent
    # print "number of muons = ", muons.size()
    for muon in muons:
        h.Fill(muon.eta(), muon.pt())
        hProfile.Fill(muon.eta(), muon.pt())

c = TCanvas("name","title",1000,800)
c.cd()
h.Draw()
c.Print("ptVSeta.pdf")

cProfile = TCanvas("nameProfile","title",1000,800)
cProfile.cd()
hProfile.Draw()
cProfile.Print("ptVSetaProfile.pdf")
