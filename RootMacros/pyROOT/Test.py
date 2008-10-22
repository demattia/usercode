#!/usr/bin/python

print "Importing ROOT"
from ROOT import *
print "Importing cmstools"
from PhysicsTools.PythonAnalysis.cmstools import *
#from cmstools import *

### prepare the FWLite autoloading mechanism 
gSystem.Load("libFWCoreFWLite.so") 
AutoLibraryLoader.enable() 

# open the file and access the event tree 
events = EventTree("/data2/demattia/Data/Z/Filter_Z_11.root")

# loop over the events
iEvent = 0
for event in events:
    prod = event.getProduct("muons")
    source = event.source
    iEvent += 1
    print "iEvent = ", iEvent
    print "number of muons = ", prod.size()
