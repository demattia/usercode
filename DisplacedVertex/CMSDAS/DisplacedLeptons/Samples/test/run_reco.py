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
if sample.sampleRelease!=os.getenv("CMSSW_VERSION"):
    print "error: trying to reconstruct",sample.sampleRelease,"MC in a",\
          os.getenv("CMSSW_VERSION"),"release directory"
    sys.exit(1)

# do we have exactly one input source?
if len(sample.sampleSimFiles)==0 and sample.sampleDataSet=="":
    print "error: neither sampleSimFiles nor sampleDataSet defined in",sample.cff
    sys.exit(1)
if len(sample.sampleSimFiles)>0 and sample.sampleDataSet!="":
    print "error: both sampleSimFiles and sampleDataSet defined in",sample.cff
    sys.exit(1)
if sample.sampleDataSet!="":
    sample.run_on="GRID"
else:
    sample.run_on=sample.where_am_i

# make sure input files exist and are accessible
if len(sample.sampleSimFiles)>0: sample.check_dir(sample.sampleSimFiles,"input")

# do we have output files defined?
if len(sample.sampleRecoFiles)==0:
    print "error: no reco files defined as output in",sample.cff
    sys.exit(1)

# when running on local input files, do we have matching output?
if len(sample.sampleSimFiles)>0 and len(sample.sampleSimFiles)!=len(sample.sampleRecoFiles):
    print "error: different number of sim and reco files in",sample.cff
    sys.exit(1)

# we need a reference to a generator config file for cmsDriver
if sample.sampleGeneratorConfig=="" and sample.sampleGeneratorReference=="":
    print "error: no reference to generator config (sampleGeneratorConfig"\
          +" or sampleGeneratorReference) in",sample.cff
    sys.exit(1)
if sample.sampleGeneratorConfig=="": sample.sampleGeneratorConfig=sample.sampleGeneratorReference

# find out whether we should use standard or modified reconstruction
# we guess this from the sample cff name, reco and pat file names and directories
modreco=0
stdreco=0
if sample.cff.lower().find("modreco")>0: modreco=1
if sample.cff.lower().find("stdreco")>0: stdreco=1
for entry in sample.sampleRecoFiles+sample.samplePatFiles:
    if entry.lower().find("modreco")>0: modreco=1
    if entry.lower().find("stdreco")>0: stdreco=1
if modreco==1 and stdreco==1:
    print "error: mixture of stdreco and modreco in",sample.cff
    sys.exit(1)
if modreco==0 and stdreco==0:
    print "error: cannot identify whether to use modreco or stdreco in",sample.cff
    sys.exit(1)
recomode="stdreco"
if modreco: recomode="modreco"


#############################################
# CREATE cmsRun CONFIG FILE WITH cmsDriver.py
#############################################

# create workdir
sample.create_work_dir(recomode)

# check generator configuration file specified in sample file
if len(sample.sampleSimFiles)>0:
    # this is a locally produced sample. we should have a local generator config file
    configfile=sample.cmssw_base+"/src/"+sample.packageName+"/"+sample.sampleGeneratorConfig
else:
    # looks like we are running over a GEN-SIM-RAW sample from DBS.
    os.system("cd "+sample.workdir+"; cvs co "+sample.sampleGeneratorConfig)
    configfile=sample.workdir+"/"+sample.sampleGeneratorConfig
if not os.path.isfile(configfile):
    print "error: could not find sampleGeneratorConfig specified in",sample.cff
    sys.exit(1)
                        

# now run cmsDriver
sample.run_cmsDriver(configfile,"RAW2DIGI,RECO","GEN-SIM-RECO","RECOSIM")

# for modified reco, we need to include additional stuff in the config file
if recomode=="modreco":
    os.system("cp -p "+sample.driverconf+" "+sample.driverconf+".orig")
    configfile=open(sample.driverconf,"a")
    configfile.write("\nprocess.load(\""+sample.packageNamePython+".modreco.modReco_cff\")\n")
    configfile.close()

# make sure the output directory is ready
sample.check_dir(sample.sampleRecoFiles,"output")


##########
# RUN JOBS
##########

# are we running on the GRID or locally?
if sample.sampleDataSet!="":
    if sample.store_at=="LOCAL":
        print "error: GRID processing requested with output to local files"
        sys.exit(1)
    sample.run_on="GRID"
elif sample.where_am_i=="RAL" and (sample.read_from=="RAL" or sample.read_from=="LOCAL")\
         and (sample.store_at=="RAL" or sample.store_at=="LOCAL"):
    sample.run_on="RAL"
elif sample.where_am_i=="CERN" and (sample.read_from=="CERN" or sample.read_from=="LOCAL")\
         and (sample.store_at=="CERN" or sample.store_at=="LOCAL"):
    sample.run_on="CERN"
else:
    print "error: running at",sample.where_am_i,"to read from",sample.read_from,\
          "and store at",sample.store_to,"does not work!"
    sys.exit(1)

if sample.run_on=="GRID":
    sample.run_crab(sample.sampleRecoFiles,250)
else:
    sample.run_batch_jobs(sample.sampleSimFiles,sample.sampleRecoFiles,sample.run_on)


# monitor jobs (after a short delay to make sure jobs have been registered)
import time
time.sleep(5)
job_watcher(sample.workdir)
