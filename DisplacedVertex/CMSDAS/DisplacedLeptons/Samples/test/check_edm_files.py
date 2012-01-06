#!/usr/bin/env python
from DisplacedLeptons.Samples.sampletools import *

# check command line argument
import os,sys
if len(sys.argv)<2:
    print "error: need to provide a valid sample cff file as argument"
    sys.exit(1)
filename=sys.argv[1]
if len(sys.argv)==3:
    print "doing quick check mode, not attempting to open root file"
    detailed_check=0
    if sys.argv[2].find(".py")>0: filename=sys.argv[2]
else:
    detailed_check=1

sample=AnalysisSample(filename)

# look at all categories of files
if len(sample.sampleSimFiles)>0:
    print "SIM FILES:"
    check_files(sample,sample.sampleSimFiles,detailed_check)
if len(sample.sampleRecoFiles)>0:
    print "RECO FILES:"
    check_files(sample,sample.sampleRecoFiles,detailed_check)
if len(sample.samplePatFiles)>0:
    print "PAT FILES:"
    check_files(sample,sample.samplePatFiles,detailed_check)
