#!/usr/bin/env python
import os,sys
from DisplacedLeptons.Samples.sampletools import *

workdir=os.getenv("CMSSW_BASE","NONE")
if workdir=="NONE":
    print "error: need to set up release area first"
    sys.exit(1)
workdir=workdir+"/src"

crabdirs=os.popen("cd "+workdir+"; ls -1d workdirs/*/workdir_*").readlines()
status_survey={}
for crabdir in crabdirs:
    jobdir=crabdir[:crabdir.rfind("/")]
    print "=====",jobdir
    job_watcher_grid(jobdir,1)

