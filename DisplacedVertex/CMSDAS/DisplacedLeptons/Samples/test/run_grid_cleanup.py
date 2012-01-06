#!/bin/env python

#################################
# GENERAL CHECKS AND PREPARATIONS
#################################

# check command line argument
import os,sys
if len(sys.argv)!=2:
    print "error: need to provide a valid working directory as argument"
    sys.exit(1)

from DisplacedLeptons.Samples.sampletools import *
grid_cleanup(sys.argv[1])
