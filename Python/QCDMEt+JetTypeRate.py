#!/usr/bin/python

import string
import sys
from math import sqrt

type = sys.argv[1]

print
print "type = ",
print type
print

# Define the tuple containing the cross sections
# Including high luminosity factor
crossSections = [ 163000000.*0.01,
                  21600000.*0.01,
                  3080000.*0.01,
                  494000.*0.01,
                  101000.*0.01,
                  24500.*0.01,
                  6240.*0.01]

# Define the tuple with the names
qcdNames = [ "30-50", "50-80", "80-120", "120-170",
             "170-230", "230-300", "300-380" ]

singleBinMultijetRate =    [ 0, 0, 0, 0, 0, 0, 0 ]
singleBinMEtJetRate =      [ 0, 0, 0, 0, 0, 0, 0 ]
binEvents =                [ 0, 0, 0, 0, 0, 0, 0 ]
binEffPassingMultijet = [ 0, 0, 0, 0, 0, 0, 0 ]
binEffPassingMEtJet =   [ 0, 0, 0, 0, 0, 0, 0 ]

counter = 0
cumulativeMultijetRate = 0.
cumulativeMEtJetRate = 0.
# Loop on the cross sections
for x in crossSections:

    # Open a file in read mode
    name = "/data/demattia/Merge/TK3/Efficiency_QCD_" + qcdNames[counter] + "_tk3.txt"
    f=open(name, 'r')

    # ',' at the end needed to avoid newline
    for line in f:
        # No parenthesys needed, put them for clarity
        # if "Eff multijet" in line:

        if (line.find("Eff MEt + " + type + " jet") != -1):
            tupleLine = string.split(line)
            cumulativeMEtJetRate += float(tupleLine[6])*x
            singleBinMEtJetRate[counter] = float(tupleLine[6])*x
            binEffPassingMEtJet[counter] = float(tupleLine[6])

        #elif (line.find("Eff Multijet + forward jet") != -1) & (line.find("no-forward") == -1):
        #    """ Using string.split and .join to produce a
        #    tuple and emulate a simple awk"""
        #    tupleLine = string.split(line)
        #    cumulativeMultijetRate += float(tupleLine[3])*x
        #    singleBinMultijetRate[counter] = float(tupleLine[3])*x
        #    binEffPassingMultijet[counter] = float(tupleLine[3])

        elif "Total events" in line:
            tupleLine = string.split(line)
            binEvents[counter] = int(tupleLine[3])
    counter = counter+1

# Print the Multijet trigger rates in each qcd bin
#print "Multijet rates"
#print "--------------"
#counter = 0
#effTotalErr = 0.
#for rate in singleBinMultijetRate:
#    print "QCD_" + qcdNames[counter] + " Multijet rate = ",
#    print rate,
#    print " +- ",
#    effErr = sqrt(binEffPassingMultijet[counter]*binEvents[counter])/binEvents[counter]
#    effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
#    print crossSections[counter]*effErr
#    counter += 1

#print
#print "Multijet cumulative rate = ",
#print cumulativeMultijetRate,
#print " +- ",
#print sqrt(effTotalErr)
#print

# Print the MEt + Jet trigger rates in each qcd bin
print "MEt + " + type + " jet rates"
print "---------------"
counter = 0
effTotalErr = 0.
for rate in singleBinMEtJetRate:
    print "QCD_" + qcdNames[counter] + " MEt + " + type + " jet rate = ",
    print rate,
    print " +- ",
    effErr = sqrt(binEffPassingMEtJet[counter]*binEvents[counter])/binEvents[counter]
    effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
    print crossSections[counter]*effErr
    counter += 1

print
print "MEt + " + type + " jet cumulative rate = ",
print cumulativeMEtJetRate,
print " +- ",
print sqrt(effTotalErr)
print

# Print the total number of events in each bin
print "Number of events"
print "----------------"
counter = 0
for num in binEvents:
    print "Events in qcd bin " + qcdNames[counter] + " = ",
    print num
    counter += 1
