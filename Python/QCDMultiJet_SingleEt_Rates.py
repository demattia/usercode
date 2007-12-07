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

EtBins = ( "Et1", "Et2", "Et3", "Et4" )

for bin in EtBins:

    singleBinEtRate = [0] * 7
    # singleBinEt2Rate = [ 0, 0, 0, 0, 0, 0, 0 ]
    # singleBinEt3Rate = [ 0, 0, 0, 0, 0, 0, 0 ]
    # singleBinEt4Rate = [ 0, 0, 0, 0, 0, 0, 0 ]
    binEvents =       [0] * 7
    binEffPassingEt = [0] * 7
    # binEffPassingEt2 = [ 0, 0, 0, 0, 0, 0, 0 ]
    # binEffPassingEt3 = ( 0, 0, 0, 0, 0, 0, 0 )
    # binEffPassingEt4 = (0) * 7

    counter = 0
    cumulativeEtRate = 0.

    # Loop on the cross sections
    for x in crossSections:

        # Open a file in read mode
        name = "/data/demattia/Merge/TK3/Efficiency_QCD_" + qcdNames[counter] + "_tk3.txt"
        f=open(name, 'r')

        # ',' at the end needed to avoid newline
        for line in f:
            # No parenthesys needed, put them for clarity
            # if "Eff multijet" in line:
            if (line.find(bin) != -1) & (line.find("no-forward") == -1) & (line.find("central") == -1) & (line.find("tau") == -1) & (line.find("forward") == -1):
                tupleLine = string.split(line)
                # Carefull, this will not work if the power is more than 9
                tupleLine = string.split(tupleLine[3], "e-0")
                if len(tupleLine) == 2:
                    value = float(tupleLine[0])*10**(-float(tupleLine[1]))
                elif len(tupleLine) == 1:
                    value = float(tupleLine[0])
                # print "counter = ",
                # print counter,
                # print ", value = ",
                # print value
                cumulativeEtRate += float(value)*x
                singleBinEtRate[counter] = float(value)*x
                binEffPassingEt[counter] = float(value)

            elif "Total events" in line:
                tupleLine = string.split(line)
                # print tupleLine
                binEvents[counter] = int(tupleLine[3])
                counter = counter+1

    # Print the Multijet trigger rates in each qcd bin
    print bin +" rates"
    print "-----------"
    counter = 0
    effTotalErr = 0.
    for rate in singleBinEtRate:
        #    print "events passing = ",
        #    print binEffPassingMultijet[counter]*binEvents[counter]
        # print "counter = ",
        # print counter
        print "QCD_" + qcdNames[counter] + "_" + bin + " rate = ",
        print rate,
        print " +- ",
        effErr = sqrt(binEffPassingEt[counter]*binEvents[counter])/binEvents[counter]
        effTotalErr += effErr*effErr*crossSections[counter]*crossSections[counter]
        print crossSections[counter]*effErr
        counter += 1

    print
    print bin + " cumulative rate = ",
    print cumulativeEtRate,
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
