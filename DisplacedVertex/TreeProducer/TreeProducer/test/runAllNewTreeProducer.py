#!/bin/env python

"""
...
"""

import os,sys

data = [
   "Data_Mu_Run2012A1", 
   # "Data_Mu_Run2012B1", 
   "Data_Photon_Run2012A1", 
   "Data_Photon_Run2012B1", 
    ]

mc_signal = [
   # "Signal_1000_020F", 
   # "Signal_1000_050F", 
   # "Signal_1000_150F", 
   # "Signal_1000_350F", 
   # "Signal_120_020F", 
   # "Signal_120_050F", 
   # "Signal_200_020F", 
   # "Signal_200_050F", 
   # "Signal_400_005L", 
   # "Signal_400_020F", 
   # "Signal_400_050F", 
   # "Signal_400_150F",
   "HZZ"
    ]

mc_vv = [
   "TTJets", 
   "WW", 
   "WZ", 
   "ZZ", 
   "DYJets10",
   "DYJets50",
   #"DisplacedE_50GeV_stdRECO", 
   #"DisplacedMu_50GeV_stdRECO", 
    ]

mc_qcd = [
"QCDmu15",
"QCDmu30",
"QCDmu120",
"QCDmu170",
"QCDmu300",
"QCDmu470",
"QCDmu600",
"QCDmu800",
"QCDmu1000",
"QCDem20",
"QCDem30",
"QCDem80",
"QCDem170",
"QCDem250",
"QCDem350",
         # "QCDem20",
         # "QCDem30",
         # "QCDem80",
         # "QCDem170",
         # "QCDem250",
         # "QCDem350" ,
         # "QCDmu15", 
         # "QCDmu20", 
         # "QCDmu30", 
         # "QCDmu50", 
         # "QCDmu80",
         # "QCDmu120", 
         # "QCDmu170",
         # "QCDmu300", 
         # "QCDmu470",
         # "QCDmu600",
         # "QCDmu800",
         # "QCDmu1000",
    ]

mc_benchmark = [
   # "DisplacedMu_50GeV_stdRECO",
   # "DisplacedE_50GeV_stdRECO"
    ]

mc_background = []
mc_background += mc_vv
mc_background += mc_qcd

sample_test = [ "DisplacedMu_50GeV_stdRECO", "Zmumu", "Data_Mu_Run2011A1", "Data_Photon_Run2011A1", "Signal_1000_350F"]

sample_list = []

nn = "enter one of:\n\tdata \n\tsignal \n\tbackground \n\tqcd \n\tvectorboson \n\tall \n\ttest\n\n and specify runmode" 
if len(sys.argv)!=3:
    print nn
    sys.exit(1)

if sys.argv[1] == "all" :
    sample_list += data
    sample_list += mc_signal
    sample_list += mc_background
    sample_list += mc_benchmark
elif sys.argv[1] == "data" :
    sample_list += data
elif sys.argv[1] == "signal" :
    sample_list += mc_signal
elif sys.argv[1] == "background" :
    sample_list += mc_background
elif sys.argv[1] == "qcd" :
    sample_list += mc_qcd
elif sys.argv[1] == "vectorboson" :
    sample_list += mc_vv
elif sys.argv[1] == "test" :
    sample_list += sample_test
elif sys.argv[1] == "selected" :
    sample_list += data
    sample_list += mc_vv
else :
    print "invalid input >>", sys.argv[1], ">> try again!"
    print nn
    sys.exit(1)

runmode=sys.argv[2]

print "processing:", sample_list

    
ldir = os.environ['LOCALRT'] + "/src/SampleFiles/Samples/python/"
sdir = os.environ['LOCALRT'] + "/src/TreeProducer/TreeProducer/test/"

for fn in sample_list :
    cmd = sdir + "run_newTreeProducer.py " + ldir + fn + "_cff.py " + runmode
    print cmd
    os.system(cmd)


print "Done submitting."

