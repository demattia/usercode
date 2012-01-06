#!/bin/env python

#################################
# GENERAL CHECKS AND PREPARATIONS
#################################

# check command line argument
import os,sys
if len(sys.argv)!=2:
    print "error: need to provide a valid sample cff file as argument"
    sys.exit(1)
from DisplacedLeptons.Samples.sampletools import *
sample=AnalysisSample(sys.argv[1])

# check CMSSW version
if sample.sampleRelease>os.getenv("CMSSW_VERSION"):
    print "error: trying to analyse",sample.sampleRelease,"sample in a",\
          os.getenv("CMSSW_VERSION"),"release directory"
    sys.exit(1)

# do we have input?
if len(sample.samplePatFiles)==0:
    print "error: no samplePatFiles defined in",sample.cff
    sys.exit(1)
filelist=sample.check_dir(sample.samplePatFiles,"input",1)


##########################
# SET UP CMSSW CONFIG FILE
##########################

sample.create_work_dir("analysis")
infile=open(sample.cmssw_base+"/src/HarderAnalysis/DisplacedDileptons"\
            +"/test/main_cfg.py","r")
sample.driverconf=sample.workdir+"/main_cfg.py"
outfile=open(sample.driverconf,"w")
written=0
pythonId=sample.packageNamePython+".samples"
for line in infile.readlines():
    if line.find(pythonId)>0:
        # replace original sample include line
        if written==0:
            written=1
            outfile.write("from "+pythonId+"."+sample.id+"_cff import *\n")
    else:
        outfile.write(line)
infile.close()
outfile.close()


##############
# RUN THE JOBS
##############

# check whether we are logged on to the right system to run these jobs
if (sample.read_from=="RAL" or sample.read_from=="LOCAL"):
    sample.run_on="RAL"
elif sample.read_from=="FNAL":
    sample.run_on="FNAL"
else:
    print "error: reading from",sample.read_from,"does not work!"
    sys.exit(1)

# now submit jobs
sample.run_batch_jobs(filelist,[],sample.run_on)

# give jobs time to be registered
import time
time.sleep(5)

# monitor jobs
# job_watcher(sample.workdir)
