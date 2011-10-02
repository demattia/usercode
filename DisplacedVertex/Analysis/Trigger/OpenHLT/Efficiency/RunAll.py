#!/usr/bin/python

import os
import math

def computeEff(names):
    # Loop on the first file to initialize the map
    effs={}
    file = open("parallelCutEfficiency_"+names[0]+".txt")
    for line in file:
        effs[line.split()[0]] = float(line.split()[1])
        # print effs[line.split()[0]]

    # Loop on the rest to sum the values
    for j in range(1, len(names)):
        file = open("parallelCutEfficiency_"+names[j]+".txt")
        for line in file:
            effs[line.split()[0]] += float(line.split()[1])
            # print effs[line.split()[0]]

    # Take the average over the samples
    for key in sorted(effs.iterkeys()):
        # The last one is the background
        effs[key] = effs[key]/len(names)
        print "eff["+key+"] =", effs[key]

    return effs

# names=['MH120MFF50', 'MH200MFF20', 'MH200MFF50', 'MH400MFF05', 'MH400MFF20', 'MH400MFF50', 'MH400MFF150', 'MH1000MFF20', 'MH1000MFF50', 'MH1000MFF150', 'MH1000MFF350', 'part2']
# Some files still need to be copied, leave them out for now
namesSignal=['MH120MFF50', 'MH400MFF50', 'MH1000MFF20', 'MH1000MFF150']
namesBackground=['part2']

startingParallelDiffCut = 1.5


# names=['part2']

for name in namesSignal:
    print "running on signal sample", name
    os.system("./Run.sh \""+name+"\""+" "+str(startingParallelDiffCut))
    # print "./Run.sh \""+name+"\""+" "+str(startingParallelDiffCut)

effsSignal = computeEff(namesSignal)

for name in namesBackground:
    print "running on background sample", name
    os.system("./Run.sh \""+name+"\""+" \""+str(startingParallelDiffCut)+"\"")

effsBackground = computeEff(namesBackground)

# Draw the histogram
from ROOT import gROOT, TCanvas, TH1F, TProfile

gROOT.Reset()
c1 = TCanvas( 'c1', '<eff(signal)>/sqrt(eff(background))' )
 
#
# Create a one dimensional function and draw it
#

histo = TProfile('effSoverSqrtEffB', '<eff(signal)>/sqrt(eff(background))', int(len(effsSignal)), float(sorted(effsSignal.iterkeys())[0]), float(sorted(effsSignal.iterkeys())[len(effsSignal)-1]))
for key in sorted(effsSignal.iterkeys()):
    print "eff(signal)/sqrt(eff(background))["+key+"] =", effsSignal[key]/math.sqrt(effsBackground[key])
    histo.Fill(float(key), effsSignal[key]/math.sqrt(effsBackground[key]))
histo.Draw()
histo.SetXTitle("parallelism cut")
c1.SaveAs("canvas.pdf")
