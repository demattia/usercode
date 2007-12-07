#!/usr/bin/python

import string
from math import sqrt

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
singleBinMultijetOfflineRate = [0] * 7
singleBinMEtJetRate =      [ 0, 0, 0, 0, 0, 0, 0 ]
singleBinMEtJetOfflineRate =      [0] * 7
binEvents =                [ 0, 0, 0, 0, 0, 0, 0 ]
binEffPassingMultijet = [ 0, 0, 0, 0, 0, 0, 0 ]
binEffPassingMultijetOffline = [0] * 7
binEffPassingMEtJet =   [ 0, 0, 0, 0, 0, 0, 0 ]
binEffPassingMEtJetOffline = [0] * 7

counter = 0
cumulativeMultijetRate = 0.
cumulativeMultijetOfflineRate = 0.
cumulativeMEtJetRate = 0.
cumulativeMEtJetOfflineRate = 0.
# Loop on the cross sections
for x in crossSections:

    # Open a file in read mode
    name = "/data/demattia/Merge/TK3/Efficiency_QCD_" + qcdNames[counter] + "_tk3.txt"
    f=open(name, 'r')

    # ',' at the end needed to avoid newline
    for line in f:
        # No parenthesys needed, put them for clarity
        # if "Eff multijet" in line:
        if (line.find("Eff multijet") != -1) & (line.find("no-forward") == -1):
            """ Using string.split and .join to produce a
            tuple and emulate a simple awk"""
            tupleLine = string.split(line)
            cumulativeMultijetRate += float(tupleLine[3])*x
            singleBinMultijetRate[counter] = float(tupleLine[3])*x
            binEffPassingMultijet[counter] = float(tupleLine[3])
        elif (line.find("offline efficiency after multijet") != -1):
            tupleLine = string.split(line)
            cumulativeMultijetOfflineRate += float(tupleLine[6])*x
            singleBinMultijetOfflineRate[counter] = float(tupleLine[6])*x
            binEffPassingMultijetOffline[counter] = float(tupleLine[6])
        elif (line.find("Eff MEt + Jet") != -1) & (line.find("no-forward") == -1):
            tupleLine = string.split(line)
            cumulativeMEtJetRate += float(tupleLine[5])*x
            singleBinMEtJetRate[counter] = float(tupleLine[5])*x
            binEffPassingMEtJet[counter] = float(tupleLine[5])
        elif (line.find("offline efficiency after MEt + Jet") != -1):
            tupleLine = string.split(line)
            cumulativeMEtJetOfflineRate += float(tupleLine[8])*x
            singleBinMEtJetOfflineRate[counter] = float(tupleLine[8])*x
            binEffPassingMEtJetOffline[counter] = float(tupleLine[8])
        elif "Total events" in line:
            tupleLine = string.split(line)
            binEvents[counter] = int(tupleLine[3])
    counter = counter+1

# Print the Multijet trigger rates in each qcd bin
print "Multijet rates"
print "--------------"
counter = 0
effTotalErr = 0.
for rate in singleBinMultijetRate:
#    print "events passing = ",
#    print binEffPassingMultijet[counter]*binEvents[counter]
    print "QCD_" + qcdNames[counter] + " Multijet rate = ",
    print rate,
    print " +- ",
    effErr = sqrt(binEffPassingMultijet[counter]*binEvents[counter])/binEvents[counter]
    effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
    print crossSections[counter]*effErr
    counter += 1

print
print "Offline after Multijet cumulative rate = ",
print cumulativeMultijetRate,
print " +- ",
print sqrt(effTotalErr)
print

# Print the offline rate after Multijet
print "Offline after Multijet rates"
print "----------------------------"
counter = 0
effTotalErr = 0.
for rate in singleBinMultijetOfflineRate:
#    print "events passing = ",
#    print binEffPassingMultijetOffline[counter]*binEvents[counter]
    print "QCD_" + qcdNames[counter] + " Offline after Multijet rate = ",
    print rate,
    print " +- ",
    effErr = sqrt(binEffPassingMultijetOffline[counter]*binEvents[counter])/binEvents[counter]
    effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
    print crossSections[counter]*effErr
    counter += 1

print
print "Offline after Multijet cumulative rate = ",
print cumulativeMultijetOfflineRate,
print " +- ",
print sqrt(effTotalErr)
print

# Print the MEt + Jet trigger rates in each qcd bin
print "MEt + Jet rates"
print "---------------"
counter = 0
effTotalErr = 0.
for rate in singleBinMEtJetRate:
    print "QCD_" + qcdNames[counter] + " MEt + Jet rate = ",
    print rate,
    print " +- ",
    effErr = sqrt(binEffPassingMEtJet[counter]*binEvents[counter])/binEvents[counter]
    effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
    print crossSections[counter]*effErr
    counter += 1

print
print "MEt + Jet cumulative rate = ",
print cumulativeMEtJetRate,
print " +- ",
print sqrt(effTotalErr)
print

# Print the Offline rates after MEt + Jet trigger in each qcd bin
print "Offline after MEt + Jet rates"
print "---------------"
counter = 0
effTotalErr = 0.
for rate in singleBinMEtJetOfflineRate:
    print "QCD_" + qcdNames[counter] + " MEt + Jet rate = ",
    print rate,
    print " +- ",
    effErr = sqrt(binEffPassingMEtJetOffline[counter]*binEvents[counter])/binEvents[counter]
    effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
    print crossSections[counter]*effErr
    counter += 1

print
print "MEt + Jet cumulative rate = ",
print cumulativeMEtJetOfflineRate,
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
