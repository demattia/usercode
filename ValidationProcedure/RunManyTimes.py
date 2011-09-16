#!/usr/bin/python

import sys
import os
import math

if len(sys.argv) != 3:
    print "Usage: RunManyTimes CFG N"
    print "where CFG is the cfg to be run and N is the number of jobs to run"
    sys.exit()

vSizes = []
vSizeAverage = 0

cfg = sys.argv[1]
numTrials = int(sys.argv[2])

for i in range(0, numTrials):
    os.system("cmsRun "+cfg+" > output.log")
    file = open("output.log")
    for line in file:
        if line.find("MemoryReport") != -1:
            print "VSIZE =", line.split()[4]
            vSizes.append(float(line.split()[4]))
            vSizeAverage += vSizes[i]

vSizeAverage /= numTrials

rms = 0
for value in vSizes:
    rms += (vSizeAverage - value)**2
    rms = math.sqrt(rms)/numTrials

print "mean memory used +/- rms =", vSizeAverage, "+/-", rms

os.system("rm output.log")
