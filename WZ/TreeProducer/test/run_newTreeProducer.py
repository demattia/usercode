#!/bin/env python

#################################
# GENERAL CHECKS AND PREPARATIONS
#################################

# check command line argument
import os,sys
if len(sys.argv)!=3:
    print "ERROR: Need to provide a valid sample cff file and the run mode as arguments"
    print "Runmode: Choose from CERN, RAL.  Others added once supported..."
    sys.exit(1)

from SampleFiles.Tools.AnalysisSample import *
sample=AnalysisSample(sys.argv[1])

runmode=sys.argv[2]

# Check runmode is supported
if not ( runmode=="CERN" or runmode=="RAL" or runmode=="FNAL" ):
    print "ERROR: Can't run at this location/this way"
    print "Runmode: Choose from CERN, RAL or FNAL.  Others added once supported..."
    sys.exit(1)

# check CMSSW version
if sample.sampleRelease>os.getenv("CMSSW_VERSION"):
    print "error: trying to analyse",sample.sampleRelease,"sample in a",\
          os.getenv("CMSSW_VERSION"),"release directory"
    sys.exit(1)

# do we have input?
if len(sample.samplePatFiles)==0:
    print "error: no samplePatFiles defined in",sample.cff
    sys.exit(1)
    
#filelist=sample.check_dir(sample.samplePatFiles,"input",1)
filelist=sample.samplePatFiles

##########################
# SET UP CMSSW CONFIG FILE
##########################

sample.create_work_dir("analysis")
# sample.create_work_dir_temp("analysis", "/uscmst1b_scratch/lpc1/3DayLifetime/"+os.environ["USER"])
infile=open(sample.cmssw_base+"/src/TreeProducer/TreeProducer/test/newTreeProducer_cfg.py","r")
sample.driverconf=sample.workdir+"/newTreeProducer_cfg.py"
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

# now submit jobs
sample.run_batch_jobs(filelist,runmode)

