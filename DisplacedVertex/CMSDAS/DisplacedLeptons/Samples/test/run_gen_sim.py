#!/bin/env python

#################################
# GENERAL CHECKS AND PREPARATIONS
#################################

# check command line argument
import os
import sys
if len(sys.argv) != 2:
    print "error: need to provide a valid sample cff file as argument"
    sys.exit(1)
from DisplacedLeptons.Samples.sampletools import *
sample = AnalysisSample(sys.argv[1])

if sample.sampleRelease != os.getenv("CMSSW_VERSION"):
    print "error: trying to generate", sample.sampleRelease, "MC in a", os.getenv("CMSSW_VERSION"),\
          "release directory"
    sys.exit(1)


#############################################
# CREATE cmsRun CONFIG FILE WITH cmsDriver.py
#############################################

# check generator configuration file specified in sample file
configfile = sample.cmssw_base + "/src/" + sample.packageName + "/" + sample.sampleGeneratorConfig
if not os.path.isfile(configfile):
    print "error: could not find sampleGeneratorConfig specified in", sample.cff
    sys.exit(1)
genconffile = open(configfile, "r")
content = ""
for line in genconffile.readlines():
    content += line
genconffile.close()
import FWCore.ParameterSet.Config as cms
code = compile(content, "<string>", "exec")
exec(code)
try:
    if generator.comEnergy.value() != sample.sampleCMSEnergy:
        print "error: different centre-of-mass energies in", sample.cff, "and", configfile
        sys.exit(1)
except:
    # generator.comEnergy not set. this usually means we are dealing with a particle gun.
    if sample.sampleCMSEnergy != 0:
        print "error: sampleCMSEnergy set, but nothing to compare it to in the generator config"
        sys.exit(1)

# make sure we are using unique run numbers (at least amongst the self-made samples)
sample.check_run_number()

# now run cmsDriver
sample.create_work_dir("gen")
sample.run_cmsDriver(configfile, "GEN,SIM,DIGI,L1,DIGI2RAW,HLT:GRun", "GEN-SIM-RAW", "RAWSIM")


##########
# RUN CRAB
##########

# generation step always to be done on the GRID. this gives us good mechanisms to create
# unique run numbers and random seeds.
sample.run_on = "GRID"

# try to make sure the output directory is ready and we are not going to overwrite files
sample.check_dir(sample.sampleSimFiles, "output")

# create crab config file
sample.run_crab(sample.sampleSimFiles, 250)

# monitor jobs
# job_watcher(sample.workdir)
