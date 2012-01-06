#!/bin/env python

def cmssw_version(release):
    # tool to make sure CMSSW releases can be used in string comparison
    # despite different numbers of digits
    newrelease=""
    for part in release.split("_"):
        try:
            num=int(part)
            newrelease+="%03i"%num
        except:
            newrelease+=part
    return newrelease


#################################
# GENERAL CHECKS AND PREPARATIONS
#################################

# check command line argument
import os,sys
if len(sys.argv)<2:
    print "error: need to provide a valid sample cff file as argument"
    sys.exit(1)
filename=sys.argv[1]
if len(sys.argv)==3:
    print "starting jobs only, will not run job watcher"
    run_job_watcher=0
    if sys.argv[2].find(".py")>0: filename=sys.argv[2]
else:
    run_job_watcher=1
from DisplacedLeptons.Samples.sampletools import *
sample=AnalysisSample(filename)

# check CMSSW version
if cmssw_version(sample.sampleRelease)>cmssw_version(os.getenv("CMSSW_VERSION")):
    print "error: trying to read",sample.sampleRelease,"MC in a",os.getenv("CMSSW_VERSION"),\
          "release directory"
    sys.exit(1)

# do we have exactly one input source?
if len(sample.sampleRecoFiles)==0 and sample.sampleDataSet=="":
    print "error: neither sampleRecoFiles nor sampleDataSet defined in",sample.cff
    sys.exit(1)
if len(sample.sampleRecoFiles)>0 and sample.sampleDataSet!="":
    print "error: both sampleRecoFiles and sampleDataSet defined in",sample.cff
    sys.exit(1)
if sample.sampleDataSet!="":
    sample.run_on="GRID"
else:
    sample.run_on="RAL"

# make sure input files exist and are accessible
if len(sample.sampleRecoFiles)>0: sample.check_dir(sample.sampleRecoFiles,"input")

# do we have output files defined?
if len(sample.samplePatFiles)==0:
    print "error: no PAT files defined as output in",sample.cff
    sys.exit(1)

# find out whether we should do re-reconstruction with non-standard settings
# we guess this from the sample cff name, reco and pat file names and directories
rereco=0
noreco=0
if sample.cff.lower().find("rereco")>0: rereco=1
if sample.cff.lower().find("modreco")>0: noreco=1
if sample.cff.lower().find("stdreco")>0: noreco=1
for entry in sample.sampleRecoFiles:
    if entry.lower().find("modreco")>0: noreco=1
for entry in sample.samplePatFiles:
    if entry.lower().find("rereco")>0: rereco=1
    if entry.lower().find("modreco")>0: noreco=1
    if entry.lower().find("stdreco")>0: noreco=1
# default, if nothing is mentioned, is not to do rereco
if noreco==0 and rereco==0: noreco=1
if noreco==rereco:
    print "error: cannot identify whether or not to do re-reco in",sample.cff
    sys.exit(1)


# make sure the output directory is ready
sample.check_dir(sample.samplePatFiles,"output")


##########################
# SET UP CMSSW CONFIG FILE
##########################

sample.create_work_dir("pat")
infile=open(sample.cmssw_base+"/src/"+sample.packageName\
            +"/test/makePATtuple_cfg.py","r")
sample.driverconf=sample.workdir+"/makePATtuple_cfg.py"
outfile=open(sample.driverconf,"w")
written=0
pythonId=sample.packageNamePython+".samples"
for line in infile.readlines():
    if line.find(pythonId)>0:
        # replace original sample include line
        if written==0:
            written=1
            outfile.write("from "+pythonId+"."+sample.id+"_cff import *\n")
    elif line.find("rereco=")>=0:
        outfile.write("rereco=%i\n"%rereco)
    else:
        outfile.write(line)
infile.close()
outfile.close()


##############
# RUN THE JOBS
##############

# are we running on the GRID or locally?
if sample.sampleDataSet!="":
    if sample.store_at=="LOCAL":
        print "error: GRID processing requested with output to local files"
        sys.exit(1)
    sample.run_on="GRID"
elif (sample.read_from=="RAL" or sample.read_from=="LOCAL")\
         and (sample.store_at=="RAL" or sample.store_at=="LOCAL"):
    sample.run_on="RAL"
else:
    print "error: reading from",sample.read_from,\
          "and store at",sample.store_to,"does not work!"
    sys.exit(1)

# now submit jobs
if sample.run_on=="GRID":
    events_per_job=50000 # appropriate for data with prefilter and QCD-like background MC
    if sample.sampleDataSet.find("Z")>=0 or sample.sampleDataSet.find("W")>=0\
       or sample.sampleDataSet.find("LongLived")>=0:
        events_per_job=10000 # much less prefilter reduction on these samples
    sample.run_crab(sample.samplePatFiles,events_per_job)
else:
    sample.run_batch_jobs(sample.sampleRecoFiles,sample.samplePatFiles,sample.run_on)

# monitor jobs (after a short delay to make sure jobs have been registered)
if run_job_watcher:
    import time
    time.sleep(5)
    job_watcher(sample.workdir)
