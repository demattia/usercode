#!/usr/bin/python

'''
This script parses the CMSSW_N.stdout files in the current dir and
computes the efficiency of the filter selection. It needs to find
the lines "Total number of events read" and "Total number of events
written" number of events, as saved by MuScleFitFilter.
'''

import os

totalReadEvents = 0
totalWrittenEvents = 0

for i in os.popen('ls'):

    if( i.find("out") != -1 ):
        f = open(i.strip())
        for line in f:
            if( line.find("Total number of events read") != -1 ):
                readEvents = line.split()[6]
            if( line.find("Total number of events written") != -1 ):
                writtenEvents = line.split()[6]
        totalReadEvents += int(readEvents)
        totalWrittenEvents += int(writtenEvents)

print "totalReadEvents = ", totalReadEvents
print "totalWrittenEvents = ", totalWrittenEvents
print "Filter efficiency = ", float(totalWrittenEvents)/float(totalReadEvents)
