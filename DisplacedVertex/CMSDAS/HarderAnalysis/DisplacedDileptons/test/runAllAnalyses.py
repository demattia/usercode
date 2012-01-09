#!/bin/env python

"""
...
"""

import os,sys

data = [
    "Data_Mu_Run2011A1", 
    "Data_Mu_Run2011A2", 
   #"Data_Mu_Run2011A3", 
   #"Data_Mu_Run2011A4", 
   #"Data_Mu_Run2011B1", 
    "Data_Photon_Run2011A1", 
    "Data_Photon_Run2011A2", 
   #"Data_Photon_Run2011A3", 
   #"Data_Photon_Run2011A4", 
   #"Data_Photon_Run2011B1", 
    ] 

mc_signal = [
    "Signal_1000_020F", 
    "Signal_1000_050F", 
    "Signal_1000_150F", 
    "Signal_1000_350F", 
    "Signal_120_020F", 
    "Signal_120_050F", 
    "Signal_200_020F", 
    "Signal_200_050F", 
   #"Signal_400_005L", 
   #"Signal_400_020F", 
    "Signal_400_050F", 
    "Signal_400_150F", 
    ]


mc_background = [
    "TTbar", 
    "WW", 
    "WZ", 
    "ZZ", 
    "Zee10", 
    "ZeeJets120", 
    "ZeeJets170", 
    "ZeeJets20", 
    "ZeeJets230", 
    "ZeeJets300", 
    "ZeeJets30", 
    "ZeeJets50", 
    "ZeeJets80", 
    "Zee", 
    "Zmumu10", 
    "ZmumuJets120", 
    "ZmumuJets170", 
    "ZmumuJets20", 
    "ZmumuJets230", 
    "ZmumuJets300", 
    "ZmumuJets30", 
    "ZmumuJets50", 
    "ZmumuJets80", 
    "Zmumu", 
    "Ztautau",
   #"DisplacedE_50GeV_stdRECO", 
   #"DisplacedMu_50GeV_stdRECO", 
    "QCDem1000", 
    "QCDem1400", 
    "QCDem170", 
    "QCDem1800", 
   #"QCDem20", 
    "QCDem300", 
   #"QCDem30", 
    "QCDem470", 
    "QCDem600", 
    "QCDem800", 
    "QCDem80", 
    "QCDmu120", 
    "QCDmu150", 
    "QCDmu15", 
    "QCDmu20", 
    "QCDmu30", 
    "QCDmu50", 
    "QCDmu80",
    ]

sample_test = ["TTbar"]

sample_list = []

nn = "enter one of:\n\tdata \n\tsignal \n\tbackground \n\tall \n\ttest" 
if len(sys.argv)!=2:
    print nn
    sys.exit(1)

if sys.argv[1] == "all" :
    sample_list += data
    sample_list += mc_signal
    sample_list += mc_background
elif sys.argv[1] == "data" :
    sample_list += data
elif sys.argv[1] == "signal" :
    sample_list += mc_signal
elif sys.argv[1] == "background" :
    sample_list += mc_background
elif sys.argv[1] == "test" :
    sample_list += sample_test
else :
    print "invalid input >>", sys.argv[1], ">> try again!"
    print nn

print "processing:", sample_list

    
ldir = os.environ['LOCALRT'] + "/src/DisplacedLeptons/Samples/python/samples/"

for fn in sample_list :
    cmd = "./run_analysis.py " + ldir + fn + "_cff.py"
    print cmd
    os.system(cmd)


print "Done submitting."

