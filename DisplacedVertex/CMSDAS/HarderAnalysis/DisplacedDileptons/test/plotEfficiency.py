#!/bin/env python

import os,sys,math
import ROOT
import time

from HarderAnalysis.DisplacedDileptons.mergeSamples import *
from HarderAnalysis.DisplacedDileptons.common import CMSPlotDecoration, set_ana_folder

myDir = os.environ["LOCALRT"]+"/src/workdirs/"
get_samples(myDir)

wdir = getWorkdirs(myDir)
workdirs_benchmark_e   =  wdir.workdirs_benchmark_e 
workdirs_benchmark_mu  =  wdir.workdirs_benchmark_mu 
workdirs_signal        =  wdir.workdirs_signal 

muAnalysis="muTrackAnalysis"
eAnalysis ="eTrackAnalysis"

fformat = ".gif"

folder = set_ana_folder()
benchmarkfolder = folder[0]
plotfolder      = folder[1]


#############################################
### DRAW EFFICIENCY HISTOGRAM
#############################################

def efficiencyPlot(workdir,histname1,histname2,filename,title,xlabel):

    # check whether normalization histogram is in prefilter file
    prefilter=0
    if histname2.find("Prefilter")>=0: prefilter=1

    # get histograms
    histfile=ROOT.TFile.Open(workdir+"/histograms.root")
    enum=histfile.Get(histname1)
    if prefilter:
        prefilterfile=ROOT.TFile.Open(workdir+"/prefilter.root")
        denom=prefilterfile.Get(histname2)
    else:
        denom=histfile.Get(histname2)
    try:
        nref=denom.Integral()
    except:
        if prefilter:
            print "ERROR: histogram",histname2,"not found in",\
                  workdir+"/prefilter.root"
        else:
            print "ERROR: histogram",histname2,"not found in",\
                  workdir+"/histograms.root"
            pass
        return -1
    try:
        nsel=enum.Integral()
    except:
        nsel=0
        pass
    print workdir.split("/")[-1],title,"average efficiency:",nsel/nref
    canv=ROOT.TCanvas()
    try:
        effi=enum.Clone()
        if enum.GetNbinsX()<denom.GetNbinsX():
            denom.Rebin(denom.GetNbinsX()/enum.GetNbinsX())
            pass
    except:
        effi=denom.Clone()
        effi.Reset()
        pass
    effi.Divide(denom)
    for i in range(effi.GetNbinsX()):
        n=denom.GetBinContent(i+1)
        if n>0:
            p=effi.GetBinContent(i+1)
            if p<0 or p>1:
                print "ERROR: efficiency",p,"in bin",i+1
            effi.SetBinError(i+1,math.sqrt(p*(1-p)/n))
        else:
            effi.SetBinError(i+1,0)
            pass
        pass
    effi.SetTitle(title)
    effi.GetXaxis().SetTitle(xlabel)
    effi.GetYaxis().SetTitle("efficiency")
    effi.SetMarkerStyle(20)
    effi.SetMinimum(0)
    effi.SetMaximum(1.1)
    effi.Draw()
    CMSPlotDecoration(filename)
    canv.Update()
    canv.Print(filename)
    canv.Clear()
    histfile.Close()
    if nref>0:
         return nsel/nref
    else:
        return -1
    pass


#############################################
### EFFICIENCY PLOTS FROM MC
#############################################


def makeEfficiencyPlots():

    print "-----------------------------------------------"
    print "--- making plots of uncorrected MC efficiencies"
    print "-----------------------------------------------"
    begintime=time.time()
    
    # tracking efficiency for electrons
    for workdir in workdirs_benchmark_e:
        sampleID=workdir.split("/")[-1].replace("_analysis","")
        efficiencyPlot(workdir,eAnalysis+"/leptons/trueLeptonRadWithTrack",
                       eAnalysis+"/gen/leptonProdVtxRadius2D",
                       benchmarkfolder+"/track_effi_electrons"+fformat,
                       "track reconstruction efficiency as function of impact parameter",
                       "d_{0} [cm]")
        pass

    # tracking efficiency for muons
    for workdir in workdirs_benchmark_mu:
        sampleID=workdir.split("/")[-1].replace("_analysis","")
        efficiencyPlot(workdir,muAnalysis+"/leptons/trueLeptonRadWithTrack",
                       muAnalysis+"/gen/leptonProdVtxRadius2D",
                       benchmarkfolder+"/track_effi_muons"+fformat,
                       "track reconstruction efficiency as function of impact parameter",
                       "d_{0} [cm]")
        pass

    # displaced lepton reconstruction efficiency
    for workdir in workdirs_benchmark_e:
        sampleID=workdir.split("/")[-1].replace("_analysis","")
        efficiencyPlot(workdir,
                       eAnalysis+"/electrons/trueProdVtxRad",
                       eAnalysis+"/gen/leptonProdVtxRadius2D",
                       benchmarkfolder+"/electron_effi_"+sampleID+fformat,
                       "electron reconstruction efficiency as function of impact parameter",
                       "d_{0} [cm]")
        pass
    for workdir in workdirs_benchmark_mu:
        sampleID=workdir.split("/")[-1].replace("_analysis","")
        efficiencyPlot(workdir,
                       muAnalysis+"/muons/trueProdVtxRad",
                       muAnalysis+"/gen/leptonProdVtxRadius2D",
                       benchmarkfolder+"/muon_effi_"+sampleID+fformat,
                       "muon reconstruction efficiency as function of impact parameter",
                       "d_{0} [cm]")
        pass

    # dilepton reconstruction efficiency plots as function of 2d decay length
    # these are uncorrected efficiencies
    # and taking the average efficiency from this plot will overestimate the
    # actual average efficiency somewhat because these plots cut off the low
    # efficiency region beyond 100cm radius.
    for workdir in workdirs_signal:
        sampleID=workdir.split("/")[-1].replace("_analysis","")
        efficiencyPlot(workdir,
                       muAnalysis+"/dileptons/trueDecayLength2D_1",
                       "isoTrackPrefilter/decayLength2D_1muon",
                       benchmarkfolder+"/dimuon1_effi_"+sampleID+fformat,
                       "dimuon reconstruction efficiency"+
                       " as function of decay length, 1 dimuon",
                       "L_{xy} [cm]")
        
        efficiencyPlot(workdir,
                       muAnalysis+"/dileptons/trueDecayLength2D_2",
                       "isoTrackPrefilter/decayLength2D_2muon",
                       benchmarkfolder+"/dimuon2_effi_"+sampleID+fformat,
                       "dimuon reconstruction efficiency"+
                       " as function of decay length, 2 dimuons",
                       "L_{xy} [cm]")
        
        efficiencyPlot(workdir,
                       eAnalysis+"/dileptons/trueDecayLength2D_1",
                       "isoTrackPrefilter/decayLength2D_1elec",
                       benchmarkfolder+"/dielectron1_effi_"+sampleID+fformat,
                       "dielectron reconstruction efficiency"+
                       " as function of decay length, 1 dielectron",
                       "L_{xy} [cm]")
    
        efficiencyPlot(workdir,
                       eAnalysis+"/dileptons/trueDecayLength2D_2",
                       "isoTrackPrefilter/decayLength2D_2elec",
                       benchmarkfolder+"/dielectron2_effi_"+sampleID+fformat,
                       "dielectron reconstruction efficiency"+
                       " as function of decay length, 2 dielectrons",
                       "L_{xy} [cm]")
        pass

    endtime=time.time()
    print "+++this took",endtime-begintime,"seconds"
    return



makeEfficiencyPlots()


