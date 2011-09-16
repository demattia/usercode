#!/usr/bin/python

import os
import math

vSizes = []
vSizeAverage = 0

numTrials = 50

for i in range(0, numTrials):
    # print "cmsRun SiStripQualityStatistics_cfg.py > output.log"
    # os.system("cmsRun SiStripQualityStatistics_cfg.py > output.log")
    os.system("cmsRun SimpleCheck.py > output.log")
    file = open("output.log")
    for line in file:
        if line.find("MemoryReport") != -1:
            print "VSIZE =", line.split()[4]
            vSizes.append(float(line.split()[4]))
            vSizeAverage += vSizes[i]

# print vSizes
vSizeAverage /= numTrials
# print vSizeAverage

rms = 0
for value in vSizes:
    rms += (vSizeAverage - value)**2
    rms = math.sqrt(rms)/numTrials

print "mean memory used +/- rms =", vSizeAverage, "+/-", rms
