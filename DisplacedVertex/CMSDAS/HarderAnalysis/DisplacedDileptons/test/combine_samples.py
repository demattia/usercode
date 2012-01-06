#!/bin/env python

# script to create final plots and numbers for the displaced dilepton analysis
# Kristian Harder, RAL
# with contributions from Mark Baber, RAL/U Bristol

do_all_plots=0
draw_overflow=1
log_plots=1

import os,sys,math
import ROOT
import time
import datetime
import tempfile

from DisplacedLeptons.Samples.sampletools import srm_path

def mergeHistogramsFiles(workdirs):
    for dir in workdirs:
        if os.path.isfile(dir+"/histograms.root"):
            continue
        command = "cd "+dir+"; hadd histograms.root"
        listDir = os.listdir(dir)
        for fileName in listDir:
            # print fileName
            if fileName.startswith("histograms_") > 0:
                # print "found file: "+fileName
                command += " "+fileName
        print "merging histograms files in "+dir
        os.system(command)

#############################################
### DEFINE SAMPLES TO USE
#############################################



# define workdir area
#myDir = "/net/unixfsrv/ppd/ryd54680/CMSSW_4_2_4_patch1/src/workdirs/"
#myDir = "/home/ppd/pff62257/CMS/CMSSW_4_2_8_patch7/src/workdirs/"
myDir = "/uscms_data/d3/demattia/CMSDAS/CMSSW_4_2_7/src/workdirs/"


workdirs_data_mu = [
    myDir + "Data_Mu_Run2011A1_analysis_20120103",
    myDir + "Data_Mu_Run2011A2_analysis_20120103",
    # myDir + "Data_Mu_Run2011A3_analysis_20111215",
    # myDir + "Data_Mu_Run2011A4_analysis_20111215",
    # myDir + "Data_Mu_Run2011B1_analysis_20111215",
    ]
workdirs_data_e = [
    myDir + "Data_Photon_Run2011A1_analysis_20120103",
    # myDir + "Data_Photon_Run2011A2_analysis_20120103", # <- running
    # myDir + "Data_Photon_Run2011A3_analysis_20111215",
    # myDir + "Data_Photon_Run2011A4_analysis_20111215",
    # myDir + "Data_Photon_Run2011B1_analysis_20111215",
    ]

workdirs_background_mu = [
    myDir + "QCDmu15_analysis_20120104",
    myDir + "QCDmu20_analysis_20120104",
    myDir + "QCDmu30_analysis_20120104",
    myDir + "QCDmu50_analysis_20120104",
    # myDir + "QCDmu80_analysis_20120104", # <- running
    myDir + "QCDmu120_analysis_20120104",
    # myDir + "QCDmu150_analysis_20120103", # needs to be checked
    myDir + "TTbar_analysis_20120104_0",
    myDir + "Zmumu_analysis_20120104",
    myDir + "Zmumu10_analysis_20120104",
    myDir + "ZmumuJets20_analysis_20120104",
    myDir + "ZmumuJets30_analysis_20120104",
    myDir + "ZmumuJets50_analysis_20120104",
    myDir + "ZmumuJets80_analysis_20120104",
    myDir + "ZmumuJets120_analysis_20120104",
    myDir + "ZmumuJets170_analysis_20120104",
    myDir + "ZmumuJets230_analysis_20120104",
    myDir + "ZmumuJets300_analysis_20120104",
    myDir + "Ztautau_analysis_20120104",
    myDir + "WW_analysis_20120104",
    myDir + "WZ_analysis_20120104",
    myDir + "ZZ_analysis_20120104"
    ]

workdirs_background_e = [
    # myDir + "QCDem20_analysis_20120104",
    # myDir + "QCDem30_analysis_20120101",
    myDir + "QCDem80_analysis_20120104",
    # myDir + "QCDem170_analysis_20120104", # needs to be checked
    # myDir + "QCDem300_analysis_20120104",
    myDir + "QCDem470_analysis_20120104",
    myDir + "QCDem600_analysis_20120104",
    # myDir + "QCDem800_analysis_20120104",
    myDir + "QCDem1000_analysis_20120104",
    myDir + "QCDem1400_analysis_20120104",
    myDir + "QCDem1800_analysis_20120104",
    myDir + "TTbar_analysis_20120104_0",
    myDir + "Zee_analysis_20120104",
    myDir + "Zee10_analysis_20120104",
    myDir + "ZeeJets20_analysis_20120104",
    myDir + "ZeeJets30_analysis_20120104",
    myDir + "ZeeJets50_analysis_20120104",
    myDir + "ZeeJets80_analysis_20120104",
    myDir + "ZeeJets120_analysis_20120104",
    myDir + "ZeeJets170_analysis_20120104",
    myDir + "ZeeJets230_analysis_20120104",
    myDir + "ZeeJets300_analysis_20120104",
    myDir + "Ztautau_analysis_20120104",
    myDir + "WW_analysis_20120104",
    myDir + "WZ_analysis_20120104",
    myDir + "ZZ_analysis_20120104"
    ]

workdirs_benchmark_e = [
    # myDir.replace("4_2_7","4_2_5")+"DisplacedE_50GeV_stdRECO_analysis_20110808", # <- generating
    # myDir.replace("4_2_7","4_2_5")+"DisplacedE_50GeV_modRECO_analysis_20110808",
    #myDir.replace("4_2_7","4_2_5")+"DisplacedE_50GeV_reRECO_analysis_20110808"
    ]

workdirs_benchmark_mu = [
    # myDir.replace("4_2_7","4_2_5")+"DisplacedMu_50GeV_stdRECO_analysis_20110808", # <- generating
    # myDir.replace("4_2_7","4_2_5")+"DisplacedMu_50GeV_modRECO_analysis_20110808",
    #myDir.replace("4_2_7","4_2_5")+"DisplacedMu_50GeV_reRECO_analysis_20110808"
    ]

workdirs_signal = [
    myDir + "Signal_120_020F_analysis_20120104",
    myDir + "Signal_120_050F_analysis_20120104",
    myDir + "Signal_200_020F_analysis_20120104",
    myDir + "Signal_200_050F_analysis_20120104",
    # myDir + "Signal_400_005L_analysis_20120104",
    # myDir + "Signal_400_020F_analysis_20120104",
    myDir + "Signal_400_050F_analysis_20120104",
    myDir + "Signal_400_150F_analysis_20120104",
    myDir + "Signal_1000_020F_analysis_20120104",
    myDir + "Signal_1000_050F_analysis_20120104",
    myDir + "Signal_1000_150F_analysis_20120104",
    myDir + "Signal_1000_350F_analysis_20120104"
    ]

# Check if histograms.root was created for all the dirs
mergeHistogramsFiles(workdirs_data_mu)
mergeHistogramsFiles(workdirs_data_e)
mergeHistogramsFiles(workdirs_background_mu)
mergeHistogramsFiles(workdirs_background_e)
mergeHistogramsFiles(workdirs_benchmark_mu)
mergeHistogramsFiles(workdirs_benchmark_e)
mergeHistogramsFiles(workdirs_signal)

############################################################################
### DEFINE ANALYSIS CHANNELS
############################################################################

muAnalysis="muTrackAnalysis"
eAnalysis="eTrackAnalysis"


############################################################################
### DEFINE REPLACEMENT TRIGGERS TO BE USED WHEN ACTUAL TRIGGER MISSING IN MC
############################################################################

replacementTrigger={}
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v1"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v2"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v3"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v4"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v5"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v6"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu23_NoVertex_v7"]="HLT_L2DoubleMu23_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu30_NoVertex_v1"]="MC_L2DoubleMu30_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu30_NoVertex_v2"]="MC_L2DoubleMu30_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu30_NoVertex_v3"]="MC_L2DoubleMu30_NoVertex_v1"
replacementTrigger["HLT_L2DoubleMu30_NoVertex_v4"]="MC_L2DoubleMu30_NoVertex_v1"

replacementTrigger["HLT_DoublePhoton33_v1"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_v2"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_v3"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_v4"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_v5"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_HEVT_v1"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_HEVT_v2"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_HEVT_v3"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton33_HEVT_v4"]="HLT_DoublePhoton33_v2"
replacementTrigger["HLT_DoublePhoton38_HEVT_v1"]="MC_DoublePhoton38_v2"
replacementTrigger["HLT_DoublePhoton38_HEVT_v2"]="MC_DoublePhoton38_v2"
replacementTrigger["HLT_DoublePhoton38_HEVT_v3"]="MC_DoublePhoton38_v2"
replacementTrigger["HLT_DoublePhoton43_HEVT_v1"]="MC_DoublePhoton43_v2"


#############################################
### DEFINE HISTOGRAM PROPERTIES
#############################################

axisRange={}
boolrange=[-0.1,1.1]
axisRange["mass"]=[0,500]
axisRange["mass_corr"]=[0,500]
axisRange["mass_triggercorr"]=[0,500]
axisRange["mass_scalecorr"]=[0,500]
axisRange["leptonEtaL"]=[-5,5]
axisRange["leptonEtaH"]=[-5,5]
axisRange["leptonD0L"]=[0,5]
axisRange["leptonD0H"]=[0,5]
axisRange["leptonD0significanceL"]=[-10,10]
axisRange["leptonD0significanceH"]=[-10,10]
axisRange["leptonAbsD0significanceL"]=[0,10]
axisRange["leptonAbsD0significanceH"]=[0,10]
axisRange["decayLengthSignificance2D"]=[0,20]
axisRange["decayLength2D"]=[0,100]
axisRange["cosThetaStar"]=[-1.1,1.1]
axisRange["leptonQualityL"]=boolrange
axisRange["leptonQualityH"]=boolrange
axisRange["leptonPtL"]=[0,100]
axisRange["leptonPtH"]=[0,100]
axisRange["trackerIsolationL"]=[0,20.]
axisRange["trackerIsolationH"]=[0,20.]
axisRange["oppositeCharge"]=boolrange
axisRange["validTracks"]=boolrange
axisRange["validVertex"]=boolrange
axisRange["vetoBackToBack"]=[-1.1,1.1]
axisRange["dPhi"]=[0,1.0]
axisRange["dPhicorr"]=[0,1.0]
axisRange["dPhitriggerCorr"]=[0,1.0]
axisRange["maxHitsBeforeVertex"]=[-1,5]
axisRange["maxHitsMissedAfterVertex"]=[-1,5]
axisRange["vertexChi2"]=[0,20]
axisRange["eta_cand"]=[-5,5]
axisRange["phi_cand"]=[-6.5,6.5]
axisRange["numCaloMatches"]=[-0.2,2.2]
axisRange["numTrigMatches"]=[-0.2,2.2]
axisRange["numStandAloneMuons"]=[-0.2,2.2]
axisRange["numPrimaryVertices"]=[0,25]
axisRange["deltaRBetweenLeptons"]=[0,6.4]

histNBins={}
histNBins["mass"]=100
histNBins["mass_corr"]=100
histNBins["mass_triggercorr"]=100
histNBins["mass_scalecorr"]=100
histNBins["leptonD0L"]=25
histNBins["leptonD0H"]=25
histNBins["leptonD0significanceL"]=25
histNBins["leptonD0significanceH"]=25
histNBins["leptonAbsD0significanceL"]=25
histNBins["leptonAbsD0significanceH"]=25
histNBins["decayLengthSignificance2D"]=25
histNBins["decayLength2D"]=25
histNBins["cosThetaStar"]=25
histNBins["leptonQualityL"]=12
histNBins["leptonQualityH"]=12
histNBins["leptonPtL"]=100
histNBins["leptonPtH"]=100
histNBins["trackerIsolationL"]=50
histNBins["trackerIsolationH"]=50
histNBins["oppositeCharge"]=12
histNBins["validTracks"]=12
histNBins["validVertex"]=12
histNBins["vetoBackToBack"]=25
histNBins["dPhi"]=25
histNBins["dPhicorr"]=25
histNBins["dPhitriggerCorr"]=25
histNBins["maxHitsBeforeVertex"]=6
histNBins["maxHitsMissedAfterVertex"]=6
histNBins["vertexChi2"]=20
histNBins["leptonEtaL"]=50
histNBins["leptonEtaH"]=50
histNBins["eta_cand"]=50
histNBins["phi_cand"]=50
histNBins["numCaloMatches"]=24
histNBins["numTrigMatches"]=24
histNBins["numStandAloneMuons"]=24
histNBins["numPrimaryVertices"]=25
histNBins["deltaRBetweenLeptons"]=32

histXLegend={}
histXLegend["massPrompt"]=0.7
histXLegend["mass"]=0.33
histXLegend["mass_corr"]=0.33
histXLegend["mass_triggercorr"]=0.33
histXLegend["mass_scalecorr"]=0.33
histXLegend["leptonD0L"]=0.6
histXLegend["leptonD0H"]=0.6
histXLegend["leptonD0significanceL"]=0.6
histXLegend["leptonD0significanceH"]=0.6
histXLegend["leptonAbsD0significanceL"]=0.6
histXLegend["leptonAbsD0significanceH"]=0.6
histXLegend["decayLengthSignificance2D"]=0.6
histXLegend["decayLength2D"]=0.6
histXLegend["cosThetaStar"]=0.6
histXLegend["leptonQualityL"]=0.3
histXLegend["leptonQualityH"]=0.3
histXLegend["leptonPtL"]=0.6
histXLegend["leptonPtH"]=0.6
histXLegend["trackerIsolationL"]=0.6
histXLegend["trackerIsolationH"]=0.6
histXLegend["oppositeCharge"]=0.3
histXLegend["validTracks"]=0.3
histXLegend["validVertex"]=0.3
histXLegend["vetoBackToBack"]=0.6
histXLegend["dPhi"]=0.6
histXLegend["dPhicorr"]=0.6
histXLegend["dPhitriggerCorr"]=0.6
histXLegend["maxHitsBeforeVertex"]=0.6
histXLegend["maxHitsMissedAfterVertex"]=0.6
histXLegend["vertexChi2"]=0.6
histXLegend["leptonEtaL"]=0.6
histXLegend["leptonEtaH"]=0.6
histXLegend["eta_cand"]=0.6
histXLegend["phi_cand"]=0.6
histXLegend["numCaloMatches"]=0.3
histXLegend["numTrigMatches"]=0.3
histXLegend["numStandAloneMuons"]=0.3
histXLegend["numPrimaryVertices"]=0.6
histXLegend["numPV_central"]=0.6
histXLegend["numPV_low"]=0.6
histXLegend["numPV_high"]=0.6
histXLegend["numPV_veryhigh"]=0.6
histXLegend["deltaRBetweenLeptons"]=0.6
histXLegend["number_of_actual_leptons"]=0.2

histTitle={}
histTitle["mass"]="dilepton invariant mass"
histTitle["mass_corr"]="dilepton invariant mass"
histTitle["mass_triggercorr"]="dilepton invariant mass"
histTitle["mass_scalecorr"]="dilepton invariant mass"
histTitle["leptonD0L"]="lepton d_{xy}"
histTitle["leptonD0H"]="lepton d_{xy}"
histTitle["leptonD0significanceL"]="lepton d_{xy}/#sigma"
histTitle["leptonD0significanceH"]="lepton d_{xy}/#sigma"
histTitle["leptonAbsD0significanceL"]="lepton |d_{xy}/#sigma|"
histTitle["leptonAbsD0significanceH"]="lepton |d_{xy}/#sigma|"
histTitle["decayLengthSignificance2D"]="transverse decay length significance"
histTitle["decayLength2D"]="transverse decay length"
histTitle["cosThetaStar"]="cos(decay helicity angle)"
histTitle["leptonQualityL"]="lepton quality cuts, lower p_t"
histTitle["leptonQualityH"]="lepton quality cuts, higher p_t"
histTitle["leptonPtL"]="lower lepton p_t"
histTitle["leptonPtH"]="higher lepton p_t"
histTitle["trackerIsolationL"]="tracker isolation, lower p_t"
histTitle["trackerIsolationH"]="tracker isolation, higher p_t"
histTitle["oppositeCharge"]="opposite charge"
histTitle["validTracks"]="valid tracks"
histTitle["validVertex"]="valid vertex"
histTitle["vetoBackToBack"]="cos(angle between leptons)"
histTitle["dPhi"]="#Delta #varphi"
histTitle["dPhicorr"]="p_t corrected #Delta #varphi"
histTitle["dPhitriggerCorr"]="trigger corrected #Delta #varphi"
histTitle["maxHitsBeforeVertex"]="number of hits before vertex"
histTitle["maxHitsMissedAfterVertex"]="missed hits after vertex"
histTitle["vertexChi2"]="vertex #chi^{2}/NDF"
histTitle["leptonEtaL"]="lepton #eta"
histTitle["leptonEtaH"]="lepton #eta"
histTitle["eta_cand"]="dilepton #eta"
histTitle["phi_cand"]="dilepton #varphi"
histTitle["numCaloMatches"]="number of matched superclusters"
histTitle["numTrigMatches"]="number of matched trigger objects"
histTitle["numStandAloneMuons"]="number of standalone muons in candidate"
histTitle["numPrimaryVertices"]="number of primary vertices"
histTitle["numPV_central"]="number of primary vertices"
histTitle["numPV_low"]="number of primary vertices, MC low"
histTitle["numPV_high"]="number of primary vertices, MC high"
histTitle["numPV_veryhigh"]="number of primary vertices, MC very high"
histTitle["deltaRBetweenLeptons"]="#Delta R between leptons"
histTitle["number_of_actual_leptons"]="number of true leptons in candidate"

axisTitle={}
axisTitle["mass"]="mass [GeV/c^{2}]"
axisTitle["mass_corr"]="mass [GeV/c^{2}]"
axisTitle["mass_triggercorr"]="mass [GeV/c^{2}]"
axisTitle["mass_scalecorr"]="mass [GeV/c^{2}]"
axisTitle["leptonD0L"]="lepton d_{xy} (cm)"
axisTitle["leptonD0H"]="lepton d_{xy} (cm)"
axisTitle["leptonD0significanceL"]="d_{xy}/#sigma"
axisTitle["leptonD0significanceH"]="d_{xy}/#sigma"
axisTitle["leptonAbsD0significanceL"]="|d_{xy}/#sigma|"
axisTitle["leptonAbsD0significanceH"]="|d_{xy}/#sigma|"
axisTitle["decayLengthSignificance2D"]="L_{xy}/#sigma"
axisTitle["decayLength2D"]="L_{xy} [cm]"
axisTitle["cosThetaStar"]="cos(#Theta*)"
axisTitle["leptonQualityL"]="result"
axisTitle["leptonQualityH"]="result"
axisTitle["leptonPtL"]="p_{t} [GeV/c]"
axisTitle["leptonPtH"]="p_{t} [GeV/c]"
axisTitle["trackerIsolationL"]="#Sigma p_{t} [GeV/c]"
axisTitle["trackerIsolationH"]="#Sigma p_{t} [GeV/c]"
axisTitle["oppositeCharge"]="result"
axisTitle["validTracks"]="result"
axisTitle["validVertex"]="result"
axisTitle["vetoBackToBack"]="cos(#alpha)"
axisTitle["dPhi"]="#Delta #varphi"
axisTitle["dPhicorr"]="#Delta #varphi"
axisTitle["dPhitriggerCorr"]="#Delta #varphi"
axisTitle["maxHitsBeforeVertex"]="number of hits"
axisTitle["maxHitsMissedAfterVertex"]="number of hits"
axisTitle["vertexChi2"]="#chi^{2}/NDF"
axisTitle["leptonEtaL"]="#eta"
axisTitle["leptonEtaH"]="#eta"
axisTitle["eta_cand"]="#eta"
axisTitle["phi_cand"]="#varphi"
axisTitle["numCaloMatches"]="number of matched superclusters"
axisTitle["numTrigMatches"]="number of matched trigger objects"
axisTitle["numStandAloneMuons"]="number of standalone muons"
axisTitle["numPrimaryVertices"]="number of primary vertices"
axisTitle["numPV_central"]="number of primary vertices"
axisTitle["numPV_low"]="number of primary vertices"
axisTitle["numPV_high"]="number of primary vertices"
axisTitle["numPV_veryhigh"]="number of primary vertices"
axisTitle["deltaRBetweenLeptons"]="#Delta R"
axisTitle["number_of_actual_leptons"]="number of leptons"

legendEntry={}
legendEntry["Zee_"]="Z/#gamma* #rightarrow ee"
legendEntry["Zmumu_"]="Z/#gamma* #rightarrow #mu#mu"
legendEntry["Ztautau_"]="Z/#gamma* #rightarrow #tau#tau"
legendEntry["ZeeJets_"]="Z/#gamma* #rightarrow ee + jet"
legendEntry["ZmumuJets_"]="Z/#gamma* #rightarrow #mu#mu + jet"
legendEntry["QCDmu_"]="QCD"
legendEntry["QCDem_"]="QCD"
legendEntry["TTbar_"]="tt"
legendEntry["WW_"]="WW"
legendEntry["WZ_"]="WZ"
legendEntry["ZZ_"]="ZZ"
legendEntry["Signal_120_020F_"] = "H(120) #rightarrow_XX(20), 1 pb"
legendEntry["Signal_120_050F_"] = "H(120) #rightarrow XX(50), 1 pb"
legendEntry["Signal_200_020F_"] = "H(200) #rightarrow XX(20), 1 pb"
legendEntry["Signal_200_050F_"] = "H(200) #rightarrow XX(50), 1 pb"
legendEntry["Signal_400_005L_"] = "H(400) #rightarrow XX(5), 1 pb"
legendEntry["Signal_400_020F_"] = "H(400) #rightarrow XX(20), 1 pb"
legendEntry["Signal_400_050F_"] = "H(400) #rightarrow XX(50), 1 pb"
legendEntry["Signal_400_150F_"] = "H(400) #rightarrow XX(150), 1 pb"
legendEntry["Signal_1000_020F_"] = "H(1000) #rightarrow XX(20), 1 pb"
legendEntry["Signal_1000_050F_"] = "H(1000) #rightarrow XX(50), 1 pb"
legendEntry["Signal_1000_150F_"] = "H(1000) #rightarrow XX(150), 1 pb"
legendEntry["Signal_1000_350F_"] = "H(1000) #rightarrow XX(350), 1 pb"


histColour={}
histColour["Zee_"]=       ROOT.kRed
histColour["Zmumu_"]=     ROOT.kRed
histColour["Ztautau_"]=   ROOT.kOrange  +7
histColour["ZeeJets_"]=   ROOT.kYellow
histColour["ZmumuJets_"]= ROOT.kYellow
histColour["QCDmu_"]=     ROOT.kBlue
histColour["QCDem_"]=     ROOT.kBlue
histColour["TTbar_"]=     ROOT.kSpring  -9
histColour["WW_"]=        ROOT.kMagenta +1
histColour["WZ_"]=        ROOT.kCyan    -9
histColour["ZZ_"]=        ROOT.kGreen   +3


#############################################
### DRAW EFFICIENCY HISTOGRAM
#############################################

def efficiencyPlot(workdir,histname1,histname2,filename,title,xlabel):

    # check whether normalization histogram is in prefilter file
    prefilter=0
    if histname2.find("Prefilter")>=0: prefilter=1

    # check whether files exist
    if not os.path.exists(workdir+"/histograms.root"):
        print "ERROR: missing file",workdir+"/histograms.root"
        return -1
    if prefilter and not os.path.exists(workdir+"/prefilter.root"):
        print "ERROR: missing file",workdir+"/prefilter.root"
        return -1

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
    CMSPrelim(filename,1)
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
### PLOT DECORATION
#############################################

def CMSPrelim(channel,MConly=0):
    info=ROOT.TLatex()
    info.SetNDC()
    if channel.find("mu")>=0 and not MConly:
        info.DrawLatex(0.18,0.95,\
                       "CMS Preliminary #sqrt{s}=7 TeV L=1.2 fb^{-1}")
    elif not MConly:
        info.DrawLatex(0.18,0.95,\
                       "CMS Preliminary #sqrt{s}=7 TeV L=1.1 fb^{-1}")
    else:
        info.DrawLatex(0.18,0.95,\
                       "CMS Simulation")
        pass
    if channel.find("mu")>=0:
        info.DrawLatex(0.18,0.83,"#mu^{+}#mu^{-}")
    else:
        info.DrawLatex(0.18,0.83,"e^{+}e^{-}")
        pass
    return


#############################################
### DRAW OVERFLOW BIN IN HISTOGRAM
#############################################

def drawhist(h,opt):

    if draw_overflow:
        # get overflow bin content
        of=h.GetBinContent(h.GetNbinsX()+1)
        ofe=h.GetBinError(h.GetNbinsX()+1)
        # get last bin content
        lb=h.GetBinContent(h.GetNbinsX())
        lbe=h.GetBinError(h.GetNbinsX())
        # add them both up in last bin
        h.SetBinContent(h.GetNbinsX(),lb+of)
        h.SetBinError(h.GetNbinsX(),math.sqrt(lbe*lbe+ofe*ofe))
        pass
    h.Draw(opt)
    return



#############################################
### DRAW DATA/MC OVERLAY PLOT FROM HISTOGRAMS
#############################################

def overlayPlot(SampleTriggerCombination,backgrounddirs,weights,\
                histname,filename,use_color=1):
    if len(backgrounddirs)!=len(weights):
        print "overlayPlot: list of weights does not match list of workdirs"
        sys.exit(1)
        pass

    # find list of triggers for MC
    mctriggers=[]
    for [dataworkdir,trigger,lumi] in SampleTriggerCombination:
        mctrig=replacementTrigger[trigger]
        if not mctrig in mctriggers: mctriggers.append(mctrig)
    if len(mctriggers)==0:
        mctriggers.append("anyTrigger")
    if len(mctriggers)!=1:
        print "ERROR: more than one MC trigger requested."
        print "       This means using the same events several times."
        print "       And most likely the normalization is done per trigger."
        sys.exit(1)
    # construct corresponding sample/trigger combinations for MC
    MCSampleTriggerCombination=[]
    for i in range(len(backgrounddirs)):
        for trigger in mctriggers:
            MCSampleTriggerCombination.append([backgrounddirs[i],trigger,weights[i]])
 
    # fill histogram for data samples
    datahisto=0
    for [workdir,trigger,lumi] in SampleTriggerCombination:
        if not datahisto:
            histfile=ROOT.TFile.Open(workdir+"/histograms.root","R")
            try:
                datahisto=histfile.Get(histname.replace("TRIGGER",trigger)).Clone()
            except:
                print "WARNING: histogram",histname,"not found in",workdir
                return
            datahisto.SetDirectory(0)
            histfile.Close()
        else:
            histfile=ROOT.TFile.Open(workdir+"/histograms.root","R")
            datahisto.Add(histfile.Get(histname.replace("TRIGGER",trigger)))
            histfile.Close()
            pass
        pass

    # background MC plots: certain samples are split in pthat bins.
    # these should show up only as a single contribution in the plots.
    # thus define groups of samples to be plotted as one:
    groupList = []
    
    # Initialise a dictionary of histograms for each group
    valHistList = {}

    # fill histogram for background samples
    for [workdir,trigger,weight] in MCSampleTriggerCombination:
        if weight<=0: continue
        histfile=ROOT.TFile.Open(workdir+"/histograms.root","R")
        try:
            valhist=histfile.Get(histname.replace("TRIGGER",trigger)).Clone()
        except:
            valhist=0
            pass
        if not valhist: continue
        valhist.SetDirectory(0)
        histfile.Close()
        # Extract the file name from the workdir name
        fileName = workdir.split("/")[-1].split("analysis")[0]
        # Sort the file into a group.
        # if group does not exist add its name as a new group and add a
        # histogram to the histogram lists
        groupName=""
        for group in groupList:
            if fileName.find(group)>=0:
                groupName = group
                pass
            pass
        if groupName=="":
            # File does not fit in any defined group

            # background MC plots: certain samples are split in pthat bins.
            # these should show up only as a single contribution in the plots.
            # thus define groups of samples to be plotted as one:
            if fileName.find("QCDem")>=0:
                groupName = "QCDem_"
            elif fileName.find("QCDmu")>=0:
                groupName = "QCDmu_"
            elif fileName.find("ZeeJets")>=0:
                groupName = "ZeeJets_"
            elif fileName.find("ZmumuJets")>=0:
                groupName = "ZmumuJets_"
            else:
                groupName = fileName
                pass

                
            groupList.append(groupName)

           
            # Add a histogram to each of the histogram lists
            valHistList[groupName]=ROOT.TH1F()
            pass
        
        # Add the values to the respective group histogram
        if valHistList[groupName].GetEntries()==0:
            valHistList[groupName]=valhist
            valHistList[groupName].Scale(weight)
        else:
            valHistList[groupName].Add(valhist,weight)
            pass
        pass


    #Move the signal MC to the back of list to ensure it is plotted on top of the other distributions
    for group in groupList:
        if group.find("Signal") >= 0:
            groupList.append(groupList.pop(groupList.index(group))) #Remove element and append to end of the list
            pass
        pass



    # Make the histograms cumulative, with the first histogram containing
    # just its own values and the last histogram containing total sum of
    # all the histograms
    for i in range(1,len(groupList)):

###########################################################
        valHistList[groupList[i]]+=valHistList[groupList[i-1]]
        
        if groupList[i].find("Signal") >= 0: 
             valHistList[groupList[i]].SetLineColor(2)  # Make the line colour unique
        else:
            valHistList[group].SetFillColor(histColour[group])
            pass
        pass
    # Reverse order of the list so that when plotted the first element is on top
    groupList.reverse()

    canv=ROOT.TCanvas("tmpcanv")

    # define legend
    if use_color:
        valLegend = ROOT.TLegend(histXLegend[histname.split("/")[-1]],0.4,histXLegend[histname.split("/")[-1]]+0.3,0.7)
    else:
        valLegend = ROOT.TLegend(histXLegend[histname.split("/")[-1]],0.55,histXLegend[histname.split("/")[-1]]+0.3,0.7)
        pass
    valLegend.SetFillColor(ROOT.kWhite)
    
    # draw background for each of the groups on the same histogram
    scalefactor=1.0
    if not use_color:
        valLegend.AddEntry(valHistList[groupList[0]],"MC","f")
    for group in groupList:
        if use_color:
            if group.find("Signal")>=0:
                valLegend.AddEntry(valHistList[group],legendEntry[group],"l")
                pass
            else:
                valLegend.AddEntry(valHistList[group],legendEntry[group],"f")
                pass
            pass
        valHistList[group].SetMarkerStyle(0)
        if (groupList.index(group) == 0):
            # first histogram to be drawn
            try:
                scalefactor=datahisto.Integral()/valHistList[group].Integral()
            except:
                pass
            print "WARNING: MC scale factor",scalefactor
            valHistList[group].Scale(scalefactor)
            xmin=valHistList[group].GetXaxis().GetXmin()
            xmax=valHistList[group].GetXaxis().GetXmax()
            try:
                ymax=max(datahisto.GetMaximum(),valHistList[group].GetMaximum())*1.3
            except:
                ymax=valHistList[group].GetMaximum()*1.3
                pass
            valHistList[group].SetMaximum(ymax)
            valHistList[group].GetXaxis().SetTitle(axisTitle[histname.split("/")[-1]])
            valHistList[group].GetYaxis().SetTitle("Entries")
            if ymax>1e5: valHistList[group].GetYaxis().SetTitleOffset(1.5)
            valHistList[group].SetTitle(histTitle[histname.split("/")[-1]])
            drawhist(valHistList[group],"h")
        elif use_color:
            # overlay
            valHistList[group].Scale(scalefactor)
            drawhist(valHistList[group],"hsame")
            pass
        pass

    if datahisto:
        datahisto.SetMarkerStyle(21)
        drawhist(datahisto,"esame")
        valLegend.AddEntry(datahisto,"data","p")
        pass
    valLegend.Draw()
    CMSPrelim(histname.split("/")[0])
    canv.Update()
    canv.Print(filename)

    valHistList.clear()
    if datahisto: datahisto.Delete()
    canv.Close()
    return


#############################################
### GET HISTOGRAM OF A TREE QUANTITY
#############################################

# note that this is a particularly complicated way of doing things.
# it would be a lot more elegant to read the tree in this python script
# and fill the histogram right here. however, reading the tree event by
# event in python is way too slow, and using functions like TTree::Project
# or TTree::Draw seems to be associated with massive memory leaks, leading
# to crashes of this script when running into the 4G memory limit per process...

def get_histogram(workdir,treedir,varname,nbins,minval,maxval,minmass,maxmass):
    # get temporary root file name
    tempfilename="/tmp/temproot_"+os.getenv("USER","default")+".root"
    
    command="TFile* infile = new TFile(\""+workdir+"/histograms.root\");\n"\
             +"TTree* tree = infile->Get(\""+treedir+"/bigTree\");\n"\
             +"TFile* outfile = new TFile(\""+tempfilename+"\",\"RECREATE\");\n"\
             +"TH1F* valhist = new TH1F(\"valhist\",\"valhist\","\
             +str(nbins)+","+str(minval)+","+str(maxval)+");\n"\
             +"TH1F* vhloose = new TH1F(\"vhloose\",\"vhloose\","\
             +str(nbins)+","+str(minval)+","+str(maxval)+");\n"\
             +"TH1F* masspasshist = new TH1F(\"masspasshist\",\"masspasshist\","\
             +str(100)+","+str(minmass)+","+str(maxmass)+");\n"\
             +"TH1F* massfailhist = new TH1F(\"massfailhist\",\"massfailhist\","\
             +str(100)+","+str(minmass)+","+str(maxmass)+");\n"\
             +"tree->Project(\"valhist\",\""+varname+"\",\"_weight*(passesAllOtherCuts)\");\n"\
             +"tree->Project(\"vhloose\",\""+varname+"\",\"_weight*(passesAllOtherCutsIgnoreLifetime)\");\n"\
             +"tree->Project(\"masspasshist\",\"mass\",\"_weight*("+varname+"_pass && passesAllOtherCuts)\");\n"
    if maxmass<10:
        command+="tree->Project(\"massfailhist\",\"mass\",\"_weight*(passesAllOtherCuts)\");\n"
    else:
        command+="tree->Project(\"massfailhist\",\"mass\",\"_weight*(!"+varname+"_pass && passesAllOtherCutsIgnoreLifetime)\");\n"
        pass
    command+="valhist->Write();\n"\
             +"vhloose->Write();\n"\
             +"masspasshist->Write();\n"\
             +"massfailhist->Write();\n"\
             +"outfile.Close();\n"\
             +"infile.Close();\n"\
             +"gApplication->Terminate();\n"
    os.system("root -b -l << EOF\n"+command+"\nEOF\n")
    histfile=ROOT.TFile.Open(tempfilename)
    valhist=histfile.Get("valhist")
    clonedvalhist=valhist.Clone("tmpval")
    clonedvalhist.SetDirectory(0)
    vhloose=histfile.Get("vhloose")
    clonedvhloose=vhloose.Clone("tmpvalloose")
    clonedvhloose.SetDirectory(0)
    masspasshist=histfile.Get("masspasshist")
    clonedmasspasshist=masspasshist.Clone("tmpmasspass")
    clonedmasspasshist.SetDirectory(0)
    massfailhist=histfile.Get("massfailhist")
    clonedmassfailhist=massfailhist.Clone("tmpmassfail")
    clonedmassfailhist.SetDirectory(0)
    histfile.Close()
    os.system("rm "+tempfilename)
    return [clonedvalhist,clonedvhloose,clonedmasspasshist,clonedmassfailhist]


def get_masses(workdir,treedir,weight=1.0):
    command="TFile* infile = new TFile(\""+workdir+"/histograms.root\");\n"\
             +"TTree* tree = infile->Get(\""+treedir+"/bigTree\");\n"\
             +"tree->SetScanField(0);\n"\
             +"tree->Scan(\"mass:_weight\","\
             +"\"passesAllCuts\"\n);"\
             +"infile.Close();\n"\
             +"gApplication->Terminate();\n"
    result=os.popen("root -b -l << EOF\n"+command+"\nEOF\n").readlines()
    masslist=[]
    for entry in result:
        fields=entry.replace("*","").strip().split()
        if len(fields)!=3: continue
        try:
            newmass=float(fields[1].strip())
            newweight=float(fields[2].strip())*weight
            masslist.append([newmass,newweight])
        except:
            pass
        pass
    return masslist


#############################################
### DRAW DATA/MC OVERLAY PLOT FROM ROOT TREE
#############################################


def treeOverlayPlot(SampleTriggerCombination,backgrounddirs,signaldirs,weights,\
                    signalWeights,lepton_name,treedir,varname,filename,
                    efficiencyListData=[],
                    efficiencyListMC=[],efficiencyListSignal=[]):
    if len(backgrounddirs)!=len(weights):
        print "treeOverlayPlot: list of weights does not match list of workdirs"
        sys.exit(1)
        pass
    if len(signaldirs)!=len(signalWeights):
        print "treeOverlayPlot: list of signal weights"+\
              "does not match list of signal workdirs"
        sys.exit(1)
        pass

    # Keep only the signal MC we are interested in
    for i in range(len(signaldirs)):
        if signaldirs[i].find("Signal_1000_350F") >= 0:

            tempDir = signaldirs[i]
            tempWeight = signalWeights[i]
            signaldirs = [tempDir]
            signalWeights =[tempWeight]
            break


    
    #####################################################################

    # find list of triggers for MC
    mctriggers=[]
    for [dataworkdir,trigger,lumi] in SampleTriggerCombination:
        mctrig=replacementTrigger[trigger]
        if not mctrig in mctriggers: mctriggers.append(mctrig)
    if len(mctriggers)!=1:
        print "ERROR: more than one MC trigger requested."
        print "       This means using the same events several times."
        print "       And most likely the normalization is done per trigger."
        sys.exit(1)
    # construct corresponding sample/trigger combinations for MC
    MCSampleTriggerCombination=[]
    for i in range(len(backgrounddirs)):
        for trigger in mctriggers:
            MCSampleTriggerCombination.append([backgrounddirs[i],trigger,weights[i]])

    ########################################################################
    MCSignalSampleTriggerCombination=[]
    for i in range(len(signaldirs)):
        for trigger in mctriggers:
            MCSignalSampleTriggerCombination.append([signaldirs[i],trigger,signalWeights[i]])
            
        
    # now loop over all trees to get minimum and maximum values (for histogram boundaries)
    found_tree=0
    files_with_data_trigger=[]
    files_with_no_trigger_whatsoever=[]
    for [workdir,trigger,lumi] in SampleTriggerCombination+MCSampleTriggerCombination+MCSignalSampleTriggerCombination:
        histfile=ROOT.TFile.Open(workdir+"/histograms.root")
        try:
            roottree=histfile.Get(treedir.replace("TRIGGER",trigger)+"/bigTree")
            numentries=roottree.GetEntries()
            files_with_data_trigger.append(workdir)
        except:
            numentries=-1
            files_with_no_trigger_whatsoever.append(workdir)
            pass
        if numentries>=0:
            newMinValue=roottree.GetMinimum(varname)
            newMaxValue=roottree.GetMaximum(varname)
            try:
                if newMinValue<minValue: minValue=newMinValue
                if newMaxValue<maxValue: maxValue=newMaxValue
                pass
            except:
                minValue=newMinValue
                maxValue=newMaxValue
                pass
            newMinMass=roottree.GetMinimum("mass")
            newMaxMass=roottree.GetMaximum("mass")
            try:
                if newMinMass<minMass: minMass=newMinMass
                if newMaxMass<maxMass: maxMass=newMaxMass
                pass
            except:
                minMass=newMinMass
                maxMass=newMaxMass
                pass
            pass
        if histfile: histfile.Close()
        pass


    if len(files_with_data_trigger)==0:
        print "ERROR: variable",varname,"not found in data samples"
        return [[],[]]
    if len(files_with_no_trigger_whatsoever)>0:
        print "WARNING: the following samples have no suitable trigger information whatsoever:"
        for sample in files_with_no_trigger_whatsoever:
            print "        ",sample

    # create histogram
    histname=varname.split("/")[-1]
    nbins=histNBins[histname]
    histMinMass=axisRange["mass"][0]
    histMaxMass=axisRange["mass"][1]
    if (minMass>=0 and minMass<histMinMass) or maxMass>histMaxMass:
        print "WARNING: candidate masses outside of histogram boundaries:",histname,
        print "histo range",histMinMass,"-",histMaxMass,", actual range",minMass,"-",maxMass 
        pass
    if axisRange.has_key(histname):
        histMinValue=axisRange[histname][0]
        histMaxValue=axisRange[histname][1]
    else:
        # round extreme values to two significant digits
        histMinValue=float("%.2e"%minValue)
        histMaxValue=float("%.2e"%maxValue)
        print "WARNING: no user-defined histogram boundaries:",histname
        pass

    if varname=="deltaRBetweenLeptons": histMaxMass=5 # J/psi

    datahisto=ROOT.TH1F(histname,varname,nbins,histMinValue,histMaxValue)
    datahisto_loose=ROOT.TH1F(histname+"_loose",varname+" loose",nbins,histMinValue,histMaxValue)
    mchisto=ROOT.TH1F(histname+"_mc",varname,nbins,histMinValue,histMaxValue)
    mchisto_loose=ROOT.TH1F(histname+"_mcloose",varname+" loose",nbins,histMinValue,histMaxValue)
    nbins_mass=100
    datahistoMassPass=ROOT.TH1F(histname+"_MassPass",varname,nbins_mass,histMinMass,histMaxMass)
    mchistoMassPass=ROOT.TH1F(histname+"_MassPass_mc",varname,nbins_mass,histMinMass,histMaxMass)
    datahistoMassFail=ROOT.TH1F(histname+"_MassFail",varname,nbins_mass,histMinMass,histMaxMass)
    mchistoMassFail=ROOT.TH1F(histname+"_MassFail_mc",varname,nbins_mass,histMinMass,histMaxMass)

    # determine whether we want to fill mass plots as well
    massplots=0
    if varname=="decayLengthSignificance2D": massplots=1
    if varname=="deltaRBetweenLeptons": massplots=1
    datamasslist=[]
    bkgmasslist=[]

    # Initialise the efficiency variables
    passesThisAndAllOtherCuts_data = 0
    passesAllOtherCuts_data        = 0

    # fill histogram for data samples
    for [workdir,trigger,lumi] in SampleTriggerCombination:
            
        [valhist,valhist_loose,masspasshist,massfailhist]=get_histogram(workdir,treedir.replace("TRIGGER",trigger),varname,nbins,histMinValue,histMaxValue,histMinMass,histMaxMass)
        datahisto.Add(valhist)
        datahisto_loose.Add(valhist_loose)
        datahistoMassPass.Add(masspasshist)
        datahistoMassFail.Add(massfailhist)

        # Store variables to calculate the cut efficiencies
        passesThisAndAllOtherCuts_data += valhist_loose.Integral()\
                                          -massfailhist.Integral()
        passesAllOtherCuts_data        += valhist_loose.Integral()

        valhist.Delete()
        valhist_loose.Delete()
        masspasshist.Delete()
        massfailhist.Delete()
        if massplots:
            datamasslist+=get_masses(workdir,treedir.replace("TRIGGER",trigger))
            pass
        pass


    groupList = []
    
    # Initialise the efficiency variables
    passesThisAndAllOtherCuts_bkg = 0
    passesAllOtherCuts_bkg        = 0
    passesThisAndAllOtherCuts_sig = 0
    passesAllOtherCuts_sig        = 0

    # Initialise a dictionary of histograms for each group
    valHistList = {}
    valHistList_loose = {}
    massPassHistList = {}
    massFailHistList = {}

    # fill histogram for background samples
    for [workdir,trigger,weight] in MCSampleTriggerCombination + MCSignalSampleTriggerCombination:
        if weight<=0: continue



        #Generate list of types to iterate through 
        if workdir.find("Signal") >= 0:
            #Sample is signal
            
            #Account for fact that 'leptons' does not contain a 'partial' tree
            if varname.find("dileptons") >= 0:
                #Tree is dileptons
                typeList = ["_background","_signal","_partial"]
            else:
                #Tree is leptons
                typeList = ["_background","_signal"]
                pass
            pass
           
        else:
            #Sample is background
            typeList = ["_background"]
            pass
        pass
    
    
        for type in typeList:
                  
             #print "\n", varname.replace("_background",type).replace("TRIGGER",trigger), "\n"

             [tempvalhist,tempvalhist_loose,tempmasspasshist,tempmassfailhist]=get_histogram(workdir,treedir.replace("_background",type).replace("TRIGGER",trigger),varname,nbins,histMinValue,histMaxValue,histMinMass,histMaxMass)

             
             if type == "_background":
                 #Initialise histograms
                 valhist = tempvalhist.Clone()
                 valhist_loose = tempvalhist_loose.Clone()
                 masspasshist = tempmasspasshist.Clone()
                 massfailhist = tempmassfailhist.Clone()
             else:
                 valhist.Add(tempvalhist)
                 valhist_loose.Add(tempvalhist_loose)
                 masspasshist.Add(tempmasspasshist)
                 massfailhist.Add(tempmassfailhist)
                 pass
             if type=="_signal":
                 passesThisAndAllOtherCuts_sig += tempvalhist_loose.Integral()\
                                                  -tempmassfailhist.Integral()
                 passesAllOtherCuts_sig        += tempvalhist_loose.Integral()

             tempvalhist.Delete()
             tempvalhist_loose.Delete()
             tempmasspasshist.Delete()
             tempmassfailhist.Delete()


     
       
        # Extract the file name from the workdir name
        fileName = workdir.split("/")[-1].split("analysis")[0]
        # Sort the file into a group.
        # if group does not exist add its name as a new group and add a
        # histogram to the histogram lists
        for group in groupList:
            if fileName.find(group.strip("_"))>=0:
                groupName = group
                break
            pass
        else:
            # File does not fit in any defined group

            # background MC plots: certain samples are split in pthat bins.
            # these should show up only as a single contribution in the plots.
            # thus define groups of samples to be plotted as one:

            if fileName.find("QCDem")>=0:
                groupName = "QCDem_"
            elif fileName.find("QCDmu")>=0:
                groupName = "QCDmu_"
            elif fileName.find("ZeeJets")>=0:
                groupName = "ZeeJets_"
            elif fileName.find("ZmumuJets")>=0:
                groupName = "ZmumuJets_"
            else:
                groupName = fileName
                
            groupList.append(groupName)

           
            # Add a histogram to each of the histogram lists
            valHistList[groupName]= ROOT.TH1F(histname+"_mc_"+groupName,varname,
                                              nbins,histMinValue,histMaxValue)
            valHistList_loose[groupName]= ROOT.TH1F(histname+"_mcloose_"+groupName,
                                                    varname+" loose",
                                                    nbins,histMinValue,histMaxValue)
            massPassHistList[groupName]= ROOT.TH1F(histname+"_MassPass_mc_"+
                                                   groupName,varname,100,
                                                   histMinMass,histMaxMass)
            massFailHistList[groupName]= ROOT.TH1F(histname+"_MassFail_mc_"+
                                                   groupName,varname,100,
                                                   histMinMass,histMaxMass)
            pass
        
        # Add the values to the respective group histogram
        valHistList[groupName].Add(valhist,weight)
        valHistList_loose[groupName].Add(valhist_loose,weight)
        massPassHistList[groupName].Add(masspasshist,weight)
        massFailHistList[groupName].Add(massfailhist,weight)

        # Store variables to calculate the efficiency
        if workdir.find("Signal")<0:
            passesThisAndAllOtherCuts_bkg += valhist_loose.Integral()\
                                             -massfailhist.Integral()
            passesAllOtherCuts_bkg        += valhist_loose.Integral()
            pass
        
        valhist.Delete()
        valhist_loose.Delete()
        masspasshist.Delete()
        massfailhist.Delete()
        if massplots:
            if workdir.find("Signal")>=0:
                #The sample is signal do not add it to bkgmasslist
                pass
            else:
                bkgmasslist+=get_masses(workdir,treedir.replace("TRIGGER",trigger),weight)
                pass
            
            pass
        pass




    #---------------------
    # cut efficiencies
    #---------------------
   
    cutName = varname.split("/")[-1]
 
    #Calculate the cut efficiency for data
    if (passesAllOtherCuts_data > 0):
        efficiency_data = 100.*float(passesThisAndAllOtherCuts_data)/float(passesAllOtherCuts_data) 
    else:
        efficiency_data = 0
          
    #Store data in a dictionary for easy passing later 
    efficiencyData = {'cutName': cutName,'passesThisAndAllOtherCuts': passesThisAndAllOtherCuts_data,
                      'passesAllOtherCuts': passesAllOtherCuts_data,'efficiency': efficiency_data}

    efficiencyListData.append(efficiencyData)

    #Calculate the cut efficiency for background MC
    if (passesAllOtherCuts_bkg > 0):
        efficiency_bkg = 100.*float(passesThisAndAllOtherCuts_bkg)/float(passesAllOtherCuts_bkg) 
    else:
        efficiency_bkg = 0
          
    #Store data in a dictionary for easy passing later 
    efficiencyMC = {'cutName': cutName,'passesThisAndAllOtherCuts': passesThisAndAllOtherCuts_bkg,
                      'passesAllOtherCuts': passesAllOtherCuts_bkg,'efficiency': efficiency_bkg}

    efficiencyListMC.append(efficiencyMC)

    
    #Calculate the cut efficiency for signal MC
    if (passesAllOtherCuts_sig > 0):
        efficiency_sig = 100.*float(passesThisAndAllOtherCuts_sig)/float(passesAllOtherCuts_sig) 
    else:
        efficiency_sig = 0
          
    #Store data in a dictionary for easy passing later 
    efficiencySignal = {'cutName': cutName,'passesThisAndAllOtherCuts': passesThisAndAllOtherCuts_sig,
                      'passesAllOtherCuts': passesAllOtherCuts_sig,'efficiency': efficiency_sig}

    efficiencyListSignal.append(efficiencySignal)

    
    # Make the histograms cumulative, with the first histogram containing
    # just its own values and the last histogram containing total sum of
    # all the histograms





############################################################################################################
#
# ADDING SIGNAL MC CONTRIBUTION
#
# This requires that a Signal sample is added to SampleTriggerCombination, this can be achieved by
# either: adding Signal as another function argument or by adding a Signal workdir to the list of
# background workdirs. Would this lead to problems elsewhere where signal is not required?
############################################################################################################

    #Move the signal MC to the back of list to ensure it is plotted on top of the other distributions
    tempGroupList = []
    for group in groupList:
        if group.find("Signal") >= 0:
            tempGroupList.append(group)
        else:
            tempGroupList.insert(0,group)
            pass
        pass

    groupList = tempGroupList
   
    tempHist = {}
    for group in groupList:
        for histList in [valHistList,valHistList_loose,massPassHistList,massFailHistList]:

            histIndex = [valHistList,valHistList_loose,massPassHistList,massFailHistList].index(histList)
            groupIndex = groupList.index(group)
       
            if (groupList.index(group) == 0):
                # first histogram index
                tempHist[histIndex] = histList[group].Clone()
                histList[group].SetFillColor(histColour[group])
            else:    
                tempHist[histIndex].Add(histList[group])  
                histList[group]=tempHist[histIndex].Clone()

                if group.find("Signal") >= 0:        # Sample is signal MC
                    histList[group].SetLineColor(2)  # Make the line colour unique
                    #histList[group].SetFillColor(0)  # Make the fill white
                else:                    
                    # Pick a color relative to the position in the list
                    #histList[group].SetFillColor(groupIndex+2)


                    histList[group].SetFillColor(histColour[group])
                    pass
                pass
            pass
        pass
    tempHist.clear()
    # Reverse order of the list so that when plotted the first element is on top
    groupList.reverse()

    canv=ROOT.TCanvas("tmpcanv")
    if log_plots: canv.SetLogy()

    # determine histogram boundaries
    ymax=max(datahisto.GetMaximum(),valHistList[groupList[0]].GetMaximum())
    ymax_loose=max(datahisto_loose.GetMaximum(),valHistList_loose[groupList[0]].GetMaximum())
    if log_plots:
        ymax*=10
        ymax_loose*=10
    else:
        ymax*=1.2
        ymax_loose*=1.2
        pass

    # define legend
    valLegend = ROOT.TLegend(histXLegend[cutName],0.62,histXLegend[cutName]+0.3,0.92)
    valLegend.SetFillColor(ROOT.kWhite)

    # draw background for each of the groups on the same histogram
    for group in groupList:

        ############################################################
        
        if group.find("Signal")>=0:
            valLegend.AddEntry(valHistList[group],legendEntry[group],"l")
            pass
        else:
            valLegend.AddEntry(valHistList[group],legendEntry[group],"f")
            pass
        
        ############################################################

        
        valHistList[group].SetMaximum(ymax)
        valHistList[group].SetMarkerStyle(0)
        if (groupList.index(group) == 0):
            # first histogram to be drawn
            valHistList[group].GetXaxis().SetTitle(axisTitle[cutName])
            valHistList[group].GetYaxis().SetTitle("Entries")
            valHistList[group].SetTitle(histTitle[histname]+", "+lepton_name+" channel")
            drawhist(valHistList[group],"h")
        else:
            # overlay
            drawhist(valHistList[group],"hsame")
            pass
        pass
    
    datahisto.SetMarkerStyle(21)
    drawhist(datahisto,"esame")
    valLegend.AddEntry(datahisto,"data","p")
    valLegend.Draw()
    CMSPrelim(lepton_name)
    canv.Update()
    canv.Print(filename)

    # also dump these histograms into a root file
    plotfile=ROOT.TFile.Open(filename.replace(".pdf",".root"),"RECREATE")
    datahisto.Write()
    for group in groupList:
        if group.find("Signal")>=0:
            signalHisto=valHistList[group].Clone(datahisto.GetName()+"_SignalMC")
            continue
        backgroundHisto=valHistList[group].Clone(datahisto.GetName()+"_BackgroundMC")
        break
    for i in range(signalHisto.GetNbinsX()):
        sumval=signalHisto.GetBinContent(i+1)
        sumerr=signalHisto.GetBinError(i+1)
        bkgval=backgroundHisto.GetBinContent(i+1)
        bkgerr=backgroundHisto.GetBinError(i+1)
        signalHisto.SetBinContent(i+1,sumval-bkgval)
        signalHisto.SetBinError(i+1,math.sqrt(sumerr*sumerr-bkgerr*bkgerr))
        pass
    signalHisto.Write()
    backgroundHisto.Write()
    plotfile.Close()

    # draw background for each of the groups on the same histogram
    for group in groupList:

        valHistList_loose[group].SetMaximum(ymax_loose)
        valHistList_loose[group].SetMarkerStyle(0)
        if (groupList.index(group) == 0):
            # first histogram to be drawn
            valHistList_loose[group].GetXaxis().SetTitle(axisTitle[cutName])
            valHistList_loose[group].GetYaxis().SetTitle("Entries")
            valHistList_loose[group].SetTitle(histTitle[histname]+", "+lepton_name+" channel, prompt")
            drawhist(valHistList_loose[group],"h")
        else:
            # overlay
            drawhist(valHistList_loose[group],"hsame")
            pass
        pass
    
    datahisto_loose.SetMarkerStyle(21)
    drawhist(datahisto_loose,"esame")
    valLegend.Draw()
    CMSPrelim(lepton_name)
    canv.Update()
    canv.Print(filename.replace(".","_prompt."))

    # also dump these histograms into a root file
    plotfile=ROOT.TFile.Open(filename.replace(".","_prompt.").replace(".pdf",".root"),"RECREATE")
    datahisto_loose.Write()
    for group in groupList:
        if group.find("Signal")>=0:
            signalHisto=valHistList_loose[group].Clone(datahisto.GetName()+"_SignalMC")
            continue
        backgroundHisto=valHistList_loose[group].Clone(datahisto.GetName()+"_BackgroundMC")
        break
    for i in range(signalHisto.GetNbinsX()):
        sumval=signalHisto.GetBinContent(i+1)
        sumerr=signalHisto.GetBinError(i+1)
        bkgval=backgroundHisto.GetBinContent(i+1)
        bkgerr=backgroundHisto.GetBinError(i+1)
        signalHisto.SetBinContent(i+1,sumval-bkgval)
        signalHisto.SetBinError(i+1,math.sqrt(sumerr*sumerr-bkgerr*bkgerr))
        pass
    signalHisto.Write()
    backgroundHisto.Write()
    plotfile.Close()

    if massplots:
        # draw mass distribution for candidates passing all cuts
        ymax=max(datahistoMassPass.GetMaximum(),massPassHistList[groupList[0]].GetMaximum())
        if log_plots:
            ymax*=10
        else:
            ymax*=1.2
            pass
        
        # legend
        massPassLegend = ROOT.TLegend(histXLegend["mass"],0.62,histXLegend["mass"]+0.3,0.92)
        massPassLegend.SetFillColor(ROOT.kWhite)
                          
        # draw background for each of the groups on the same histogram
        for group in groupList:
                      
############################################################
            if group.find("Signal")>=0:
                massPassLegend.AddEntry(massPassHistList[group],legendEntry[group],"l")
                pass
            else:
                massPassLegend.AddEntry(massPassHistList[group],legendEntry[group],"f")
                pass
              ############################################################


            
            massPassHistList[group].SetMaximum(ymax)
            massPassHistList[group].SetMarkerStyle(0)
            if (groupList.index(group) == 0):
                massPassHistList[group].SetTitle("di"+lepton_name+" invariant mass, passing all cuts")
                massPassHistList[group].GetXaxis().SetTitle("mass [GeV/c^{2}]")
                massPassHistList[group].GetXaxis().SetNdivisions(1005)
                massPassHistList[group].GetYaxis().SetTitle("Entries")
                drawhist(massPassHistList[group],"h")
            else:
                drawhist(massPassHistList[group],"hsame")
                pass
            pass
        # draw data
        datahistoMassPass.SetMarkerStyle(21)
        drawhist(datahistoMassPass,"esame")
        massPassLegend.AddEntry(datahistoMassPass,"data","p")
        massPassLegend.Draw()
        CMSPrelim(lepton_name)
        canv.Update()
        if filename.find("deltaRBetweenLeptons")<0:
            canv.Print(filename.replace("decayLengthSignificance2D","massDisplaced"))
            pass
        
        # draw mass distribution for candidates passing all other cuts but failing this one
        ymax=max(datahistoMassFail.GetMaximum(),massFailHistList[groupList[0]].GetMaximum())
        if log_plots:
            ymax*=10
        else:
            ymax*=1.2
            pass
        
        # legend
        massFailLegend = ROOT.TLegend(histXLegend["massPrompt"],0.62,histXLegend["massPrompt"]+0.3,0.92)
        massFailLegend.SetFillColor(ROOT.kWhite)

        # Draw background for each of the groups on the same histogram
        for group in groupList:

            ############################################################
        
            if group.find("Signal")>=0:
                massFailLegend.AddEntry(massFailHistList[group],legendEntry[group],"l")
                pass
            else:
                massFailLegend.AddEntry(massFailHistList[group],legendEntry[group],"f")
                pass
        
        ############################################################
            
            massFailHistList[group].SetMaximum(ymax)
            massFailHistList[group].SetMarkerStyle(0)
            if (groupList.index(group) == 0):
                massFailHistList[group].SetTitle("di"+lepton_name+" mass, prompt candidates")
                massFailHistList[group].GetXaxis().SetTitle("mass [GeV/c^{2}]")
                massFailHistList[group].GetXaxis().SetNdivisions(1005)
                massFailHistList[group].GetYaxis().SetTitle("Entries")
                drawhist(massFailHistList[group],"h")
            else:
                drawhist(massFailHistList[group],"hsame")
                pass
            pass
        # draw data
        datahistoMassFail.SetMarkerStyle(21)
        drawhist(datahistoMassFail,"esame")
        massFailLegend.AddEntry(datahistoMassFail,"data","p")
        massFailLegend.Draw()
        CMSPrelim(lepton_name)
        canv.Update()
        canv.Print(filename.replace("decayLengthSignificance2D","massPrompt").replace("deltaRBetweenLeptons","massDisplacedNoDeltaR"))
        pass
        
    valHistList.clear()
    massPassHistList.clear()
    massFailHistList.clear()
    datahisto.Delete()
    datahistoMassPass.Delete()
    datahistoMassFail.Delete()
    canv.Close()
    return [datamasslist,bkgmasslist]


###################################################
### GET WEIGHTED NUMBER OF CANDIDATES PASSING CUTS
###################################################

def number_of_surviving_candidates(SampleTriggerCombination,candtype,analysis_directory,cuts):
    print "KHDEBUG: disabled"
    return 1
    if candtype!="background" and candtype!="signal":
        print "ERROR: illegal candidate category"
        sys.exit(1)
        pass
    overall_integral=0.0
    for [workdir,trigger,weight] in SampleTriggerCombination:
        # if this is real data then the weight is actually luminosity
        if workdir.find("Data")>=0: weight=1.0
        # the big tree containing cut results for all variables
        treename=analysis_directory+"/dileptons_"+candtype+"_"+trigger\
                  +"/bigTree"
        # make ROOT print the number of entries passing a cut.
        # in order to do this, we plot an arbitrary quantity of which
        # we definitely know the range of values (e.g. oppositeCharge)
        # into a dummy histogram with just one bin
        command="TFile* infile = new TFile(\""+workdir+"/histograms.root\");\n"\
                 +"TTree* tree = infile->Get(\""+treename+"\");\n"\
                 +"TH1F* valhist = new TH1F(\"valhist\",\"valhist\","\
                 +"1,-10.,10.);\n"\
                 +"tree->Project(\"valhist\",\"oppositeCharge\",\"_weight*("+cuts+")\");\n"\
                 +"cout << \"RES\" << \"ULT \" << valhist->Integral() << endl;\n"\
                 +"infile.Close();\n"\
                 +"gApplication->Terminate();\n"
        result=os.popen("root -b -l << EOF\n"+command+"\nEOF\n").readlines()
        integral=-999;
        for line in result:
            if line.find("RESULT")>=0:
                integral=float(line.split()[1])
                pass
            pass
        if (integral>=0):
            overall_integral+=integral*weight
        else:
            print "ERROR for",workdir,trigger
            pass
        pass
    return overall_integral
   
    

#############################################
### SCAN SAMPLE DESCRIPTION FILE
#############################################

def sample_cff_code(workdir):
    # find sample description file
    if not os.path.isdir(workdir):
        print "ERROR: directory",workdir,"does not exist"
        sys.exit(1)
    sample_cff = ""
    for filename in os.listdir(workdir):
        if filename.find("_cff.py")>0:
            sample_cff=filename
            pass
        pass
    if sample_cff=="":
        print "ERROR: sample description file not found in",workdir
        sys.exit(1)
        pass
    samplefile=open(workdir+"/"+sample_cff,"r")
    content=""
    for line in samplefile.readlines():
        content+=line
        pass
    samplefile.close()
    code=compile(content,"<string>","exec")
    return code


#############################################
### GET MC SAMPLE SIZE AND CROSS-SECTION
#############################################

def get_sample_weight(mcworkdir):
    
    # get cross-section and sample size from sample description file
    sampleXSec=-1
    sampleNumEvents=0
    sampleBaseDir=""
    samplePatFiles=[]
    exec(sample_cff_code(mcworkdir))
    if sampleXSec==-1:
        print "ERROR: sample description file does not contain sampleXSec"
        sys.exit(1)
        pass
    if sampleBaseDir=="":
        print "ERROR: sample description file does not contain sampleBaseDir"
        sys.exit(1)
        pass
    if sampleNumEvents==0:
        print "ERROR: sample description file does not contain sampleNumEvents"
        sys.exit(1)
        pass
    if len(samplePatFiles)==0:
        print "ERROR: sample description file does not contain samplePatFiles"
        sys.exit(1)
        pass

    # check number of processed events (after prefilter, if any)
    histfile=ROOT.TFile.Open(mcworkdir+"/histograms.root")
    num_events_processed=0
    try:
        num_events_processed=histfile.Get("prefilterPassTwo/numSignalPass").GetEntries()
    except:
        pass

    # how many of them ended up being analyzed? (= not rejected by goodcoll)
    num_events_analyzed=0
    found=0
    try:
        num_events_analyzed=histfile.Get(eAnalysis+"/trig/eventsVsRun").GetEntries()
        found=1
    except:
        pass
    try:
        num_events_analyzed=histfile.Get(muAnalysis+"/trig/eventsVsRun").GetEntries()
        found=1
    except:
        pass
    if not found:
        print "ERROR: unable to retrieve histogram with number of events from",mcworkdir
        return 0.0

    # check number of events before prefilter (if any)
    numevents_before_prefilter=-num_events_processed
    numevents_after_prefilter=-num_events_processed
    if os.path.exists(mcworkdir+"/prefilter.root"):
        histfile=ROOT.TFile.Open(mcworkdir+"/prefilter.root","R")
        try:
            numevents_before_prefilter=histfile.Get("isoTrackPrefilter/numSignal").GetEntries()
            numevents_after_prefilter=histfile.Get("isoTrackPrefilter/numSignalPass").GetEntries()
        except:
            print "ERROR: prefilter file in",mcworkdir,"seems to be corrupt"
            pass
        pass

    # this weight will have to be multiplied with integrated lumi to get the
    # expected number of events from this physics process
    sampleweight=sampleXSec/(num_events_processed*1.0\
                             *numevents_before_prefilter\
                             /numevents_after_prefilter)

    print "MC sample",mcworkdir.split("/")[-1],":"
    print "   "\
          +str(int(.5+abs(numevents_before_prefilter)/sampleNumEvents*100.))\
          +"% ("+str(abs(int(numevents_before_prefilter)))\
          +"/"+str(sampleNumEvents)\
          +") of total sample size available in PAT"
    print "   resulting sample weight (to be multiplied by integrated lumi):",\
          sampleweight
    if numevents_before_prefilter>0:
        print "   fraction of events passing the prefilter: ",\
              float(numevents_after_prefilter)/numevents_before_prefilter
    else:
        print "   no prefilter applied"
        pass
    print "   fraction of events passing the goodcollision filter:",\
          float(num_events_analyzed)/num_events_processed

    # cross-check that we didn't lose any events in analysis stage
    # (through exceptions or unreadable files). the problem is that
    # we do lose a small number of events from the primary vertex filter.
    # to be addressed.
    if num_events_processed!=abs(numevents_after_prefilter):
        print "ERROR:",mcworkdir,"has",numevents_after_prefilter,\
              "passing the prefilter, but",\
              num_events_processed,"events processed in analysis"

    return sampleweight
    

#############################################
### PROCESS INDIVIDUAL DATASET
#############################################

def process_dataset(dataworkdir,lepton_name,analysis_directory):

    # get dCache directory of this sample from sample description file
    sampleBaseDir=""
    exec(sample_cff_code(dataworkdir))
    if sampleBaseDir=="":
        print "ERROR: sample description file does not contain sampleBaseDir"
        sys.exit(1)
        pass

    # now get analysis config file (for the list of triggers that was looked at)
    if not os.path.exists(dataworkdir+"/main_cfg.py"):
        print "ERROR: could not find main_cfg.py in workdir"
        sys.exit(1)
        pass
    ehltPaths=[]
    muhltPaths=[]
    configfile=open(dataworkdir+"/main_cfg.py","r")
    content=""
    for line in configfile.readlines():
        if line.find("HLT")>0 and line.find("pathNames")<0:
            content+=line.replace("process.","").replace("Analysis.","").replace("cms.vstring(","[").replace(")","]")
            pass
        pass
    configfile.close()
    code=compile(content,"<string>","exec")
    exec(code)
    if len(ehltPaths)==0 or len(muhltPaths)==0:
        print "ERROR: process.eTrackAnalysis.hltPaths or process.muTrackAnalysis.hltPaths not found in main_cfg.py"
        sys.exit(1)
        pass
    if lepton_name=="muon":
        selectedHLTPaths=muhltPaths
        unbiasedHLTPaths=ehltPaths
    else:
        selectedHLTPaths=ehltPaths
        unbiasedHLTPaths=muhltPaths
        pass

    # get maximum run range from histogram file
    histfile=ROOT.TFile.Open(dataworkdir+"/histograms.root")
    runrangehist=histfile.Get(analysis_directory+"/trig/eventsVsRun")
    try:
        runMin=int(runrangehist.GetXaxis().GetXmin())
        runMax=int(runrangehist.GetXaxis().GetXmax())
    except:
        # looks like this histogram file is corrupt or otherwise unusable
        print "ERROR: histogram.root unusable."
        sys.exit(1)
    histfile.Close()

    # get luminosity information
    lumiOverview=dataworkdir+"/lumiOverview.csv"
    lumiDetails=dataworkdir+"/lumiResult.csv"
    os.system("rm -f "+lumiOverview+" "+lumiDetails)
    print "srmcp -2 \""+sampleBaseDir.replace("root://xrootd.rcac.purdue.edu", "srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=")+"/lumiOverview.csv\" \"file:///"+lumiOverview+"\""
    os.system("srmcp -2 \""+sampleBaseDir.replace("root://xrootd.rcac.purdue.edu", "srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=")+"/lumiOverview.csv\" \"file:///"+lumiOverview+"\"")
    # os.system("dccp "+sampleBaseDir+"/lumiOverview.csv "+lumiOverview+" &> /dev/null")
    if not os.path.exists(lumiOverview):
        # ok, for some reason we cannot get these files from dcache.
        # check whether the make_pat working directory still exists
        # where those files were probably created in
        patdirs=os.popen("ls -1rtd "+dataworkdir[:dataworkdir.find("_analysis")]\
                         +"_pat*").readlines()
        if len(patdirs)>0:
            print "WARNING: using lumiOverview.csv from",patdirs[-1].strip()
            os.system("cp "+patdirs[-1].strip()+"/lumiOverview.csv "\
                      +lumiOverview+" &> /dev/null")
            pass
        pass
    os.system("srmcp -2 \""+sampleBaseDir.replace("root://xrootd.rcac.purdue.edu", "srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=")+"/lumiResult.csv\" \"file:///"+lumiDetails+"\"")
    # os.system("dccp "+sampleBaseDir+"/lumiResult.csv "+lumiDetails+" &> /dev/null")
    if not os.path.exists(lumiDetails):
        # same procedure as above
        patdirs=os.popen("ls -1rtd "+dataworkdir[:dataworkdir.find("_analysis")]\
                         +"_pat*").readlines()
        if len(patdirs)>0:
            print "WARNING: using lumiResult.csv from",patdirs[-1].strip()
            os.system("cp "+patdirs[-1].strip()+"/lumiResult.csv "\
                      +lumiDetails+" &> /dev/null")
            pass
        pass
    if not os.path.exists(lumiOverview) or not os.path.exists(lumiDetails):
        # ok, now give up
        print "ERROR: unable to retrieve lumi files from "+sampleBaseDir
        sys.exit(1)
        pass

    # read resulting csv file for delivered lumi and overall recorded lumi
    lumiCSV=open(lumiOverview)
    lumi_delivered=0.0
    lumi_recorded=0.0
    recLumiPerRun={}
    for line in lumiCSV.readlines():
        if len(line.split(","))==4:
            # old lumiCalc format
            (run,dlumi,rlumi,hlumi)=line.split(",")
        else:
            # lumiCalc2.py format
            entries=line.split(",")
            run=entries[0]
            dlumi=entries[2]
            rlumi=entries[-1]
            pass
        try:
            run=int(run)
            dlumi=float(dlumi)/1e6 # convert from /microbarn to /picobarn
            rlumi=float(rlumi)/1e6 # convert from /microbarn to /picobarn
        except:
            continue
        if run<runMin or run>runMax:
            print "ERROR: json contains run",run,"outside expected range [",\
                  runMin,",",runMax,"]"
            #sys.exit(1)
            pass
        lumi_delivered+=dlumi
        lumi_recorded+=rlumi
        recLumiPerRun[run]=rlumi
        pass
    lumiCSV.close()
    print "  delivered lumi for this dataset:",lumi_delivered,"/pb"
    print "  recorded lumi for this dataset:",lumi_recorded,"/pb"

    # read and analyze resulting csv file for recorded lumi
    lumiCSV=open(lumiDetails)
    lumiPerRun={}
    lumiPerTrigger={}
    firstRun={}
    lastRun={}
    maxLumi={}
    previousBestTrigger="None"
    triggerWithMaxLumi={}
    for line in lumiCSV.readlines():
        if line.startswith("Run"):
            continue
        # read line from csv
        if len(line.split(","))==3:
            # old style lumiCalc result
            (run,hltpath,lumi)=line.split(",")
        else:
            # new style lumiCalc2.py
            # print "line = "+line
            modifiedString = ""
            foundQuote = False
            for char in line:
                if char == "\"":
                    foundQuote = not foundQuote
                if char == "," and foundQuote:
                    modifiedString += "###"
                else:
                    modifiedString += char
            # modifiedString = modifiedString.split(",")
            # print "modifiedString = "+modifiedString
            # for elem in modifiedString:
            #     print elem.replace("###", ",")
        
            
            # beforeLS = line.split("[")[0]
            # afterLS = line.split("]")[-1]
            # print "beforeLS = "+beforeLS
            # print "afterLS = "+afterLS
            # run = beforeLS.split(",")[0]
            # print afterLS.split(',')
            (run, lumiSL, recorded, hltpath, l1bit, lumi) = modifiedString.split(",")
            run = run.replace("###", ",")
            lumiSL = lumiSL.replace("###", ",")
            recorded = recorded.replace("###", ",")
            hltpath = hltpath.replace("###", ",")
            l1bit = l1bit.replace("###", ",")
            lumi = lumi.replace("###", ",")
            # print run
            # print hltpath
            # print l1bit
            # print lumi
            # (run,selectedLS,recorded,hltpath,l1bit,lumi)=line.split(",")
            pass
        if not hltpath in selectedHLTPaths: continue
        try:
            run=int(run)
            lumi=float(lumi)/1e6 # convert from /microbarn to /picobarn
        except:
            continue
        # accumulate luminosity for each trigger
        if not lumiPerTrigger.has_key(hltpath):
            lumiPerTrigger[hltpath]=0.0
            lumiPerRun[hltpath]=ROOT.TH1F("lumiPerRun_"+hltpath,
                                          "recorded lumi per run for "+hltpath,
                                          runMax-runMin+1,runMin,runMax)
            firstRun[hltpath]=runMax+1
            lastRun[hltpath]=0
            pass
        lumiPerTrigger[hltpath]+=lumi
        lumiPerRun[hltpath].Fill(run,lumi)
        if lumi>0:
            if run<firstRun[hltpath]: firstRun[hltpath]=run
            if run>lastRun[hltpath]: lastRun[hltpath]=run
            pass
        # find best trigger for each run
        if not maxLumi.has_key(run):
            maxLumi[run]=0.0
            triggerWithMaxLumi[run]="None"
            pass
        if lumi>maxLumi[run]:
            if maxLumi[run]>0:
                # there already was one active trigger from our selection here.
                # this is potentially dangerous because we might double-count
                # data from this run when looking at selections with different
                # triggers. (Though we will double-count the luminosity too,
                # and thus things might be ok except for statistical correlation)
                print "WARNING: more than one analysis trigger in run",run
            maxLumi[run]=lumi
            triggerWithMaxLumi[run]=hltpath
            pass
        if lumi==maxLumi[run] and hltpath==previousBestTrigger:
            triggerWithMaxLumi[run]=hltpath
        previousBestTrigger=triggerWithMaxLumi[run]
        pass
    lumiCSV.close()

    # summarize results by trigger
    for trigger in selectedHLTPaths:
        if not lumiPerTrigger.has_key(trigger):
            lumiPerTrigger[trigger]=0.0
            firstRun[trigger]=0
            lastRun[trigger]=0
            pass
        print "  recorded lumi for %40s : %10.8f in run range %i-%i"\
              %(trigger,lumiPerTrigger[trigger],firstRun[trigger],lastRun[trigger])
        pass
    # add dummy trigger
    lumiPerTrigger["anyTrigger"]=lumi_recorded
    firstRun["anyTrigger"]=runMin
    lastRun["anyTrigger"]=runMax
    
    # summarize results by run
    runlist=[]
    lumisum=0.0
    best_triggers=[]
    for run in maxLumi.keys():
        if recLumiPerRun[run]>0:
            if not triggerWithMaxLumi[run] in best_triggers:
                best_triggers.append(triggerWithMaxLumi[run])
                pass
            runlist.append("  run "+str(run)+": best trigger "+triggerWithMaxLumi[run]+" live %4.1f%%"%(100*maxLumi[run]/recLumiPerRun[run]))
        else:
            runlist.append("  run "+str(run)+": no data recorded")
        lumisum+=maxLumi[run]
        pass
    runlist.sort()
    for line in runlist:
        print line
        pass
    
    print "  => total lumi with best triggers (",best_triggers,") is",lumisum,\
          "/pb out of",lumi_recorded,"/pb\n"
    if (lumisum-lumi_recorded)>1e-5:
        print "ERROR: sum of trigger lumis (",lumisum,") exceeds total recorded (",\
              lumi_recorded,")"
        sys.exit(1)

    best_trigger_lumi=[]
    for trigger in best_triggers:
        best_trigger_lumi.append(lumiPerTrigger[trigger])
    return [lumisum,best_triggers,best_trigger_lumi]


#############################################
### PREPARE DATA/MC OVERLAY PLOTS
#############################################


def makePlots(SampleTriggerCombination,workdirs_background,workdirs_signal,
              mcweights,mcSignalWeights,
              lepton_name,analysis_directory,plotfolder,limitfolder):


    # get total luminosity for this set of samples
    lumisum=0.0
    for [dataworkdir,trigger,lumi] in SampleTriggerCombination:
        lumisum+=lumi
        sampleID=dataworkdir.split("/")[-1].replace("_analysis","")+"_"+trigger
        pass

    # corresponding normalization for background samples
    mcfactors=[]
    for weight in mcweights:
        mcfactors.append(weight*lumisum)
        pass

#################################################################
    # corresponding normalization for signal samples
    mcfactorsSignal=[]
    for weight in mcSignalWeights:
        mcfactorsSignal.append(weight*lumisum)
        pass
#################################################################

    # Make directory to store plots
    folder = plotfolder + "/" + sampleID + "/"
    if not os.path.isdir(folder):
        os.mkdir(folder)
        pass

    print "preparing plots for the following samples and triggers:"
    for [dataworkdir,trigger,lumi] in SampleTriggerCombination:
        print "   ",trigger,"in",dataworkdir
    print "-> using directory",folder
        
    #---------------------
    # PLOTS
    #---------------------

    # primary vertex distributions
    for plot in ["numPV_central","numPV_low","numPV_high","numPV_veryhigh"]:
        overlayPlot(SampleTriggerCombination,workdirs_background,mcfactors,
                    analysis_directory+"/weights/"+plot,
                    folder+lepton_name+"_"+plot+".pdf",0)
        pass

    # fraction of actual true leptons in candidate (MC only)
    overlayPlot([],workdirs_background,mcfactors,
                analysis_directory+"/dileptons/number_of_actual_leptons",
                folder+"number_of_actual_"+lepton_name+"s.pdf",0)
    
    # Lists to store all the required data about the cut efficiencies
    efficiencyListData = []
    efficiencyListMC = []
    efficiencyListSignal = []
           
    # decay length significance. we treat this cut in a special way
    # because its inversion gives us the control distribution of prompt
    # candidates, and thus we want additional plots/diagnostics here
    [datamasslist,bkgmasslist]=treeOverlayPlot(SampleTriggerCombination,
                                               workdirs_background,workdirs_signal,mcfactors,
                                               mcfactorsSignal,
                                               lepton_name,
                                               analysis_directory+\
                                               "/dileptons_background_TRIGGER/",
                                               "decayLengthSignificance2D",
                                               folder+"di"+lepton_name+\
                                               "_decayLengthSignificance2D.pdf",
                                               efficiencyListData,efficiencyListMC,
                                               efficiencyListSignal)
    
    # store mass lists for limit calculation
    minmass=999
    maxmass=0
    datalistfile=open(limitfolder+"/masses_data_"+lepton_name+".txt","w")
    for entry in datamasslist:
        datalistfile.write("%10.4f %10.4f\n"%(entry[0],entry[1]))
        if entry[0]>maxmass: maxmass=entry[0]
        if entry[0]<minmass: minmass=entry[0]
        pass
    datalistfile.close()
    bkglistfile=open(limitfolder+"/masses_backgroundMC_"+lepton_name+".txt","w")
    for entry in bkgmasslist:
        bkglistfile.write("%10.4f %10.4f\n"%(entry[0],entry[1]))
        if entry[0]>maxmass: maxmass=entry[0]
        if entry[0]<minmass: minmass=entry[0]
        pass
    bkglistfile.close()
    # now also create histograms
    datahist=ROOT.TH1F("datamass","mass distribution of candidates in data",
                       100,minmass*0.9,maxmass*1.1)
    for entry in datamasslist: datahist.Fill(entry[0],entry[1])
    outfile=ROOT.TFile.Open(limitfolder+"/masses_data_"+lepton_name+".root",
                            "RECREATE")
    datahist.Write()
    outfile.Close()
    bkghist=ROOT.TH1F("backgroundmass","mass distribution of candidates in MC",
                       100,minmass*0.9,maxmass*1.1)
    for entry in bkgmasslist: bkghist.Fill(entry[0],entry[1])
    outfile=ROOT.TFile.Open(limitfolder+"/masses_backgroundMC_"+lepton_name+".root",
                            "RECREATE")
    bkghist.Write()
    outfile.Close()

    # get list of all quantities that are available to plot
    cutNames=[]
    for [workdir,trigger,lumi] in SampleTriggerCombination:
        histfile=ROOT.TFile.Open(workdir+"/histograms.root")
        for analysisDir in [eAnalysis,muAnalysis]:
            if histfile.Get("/"+analysisDir+"/dileptons_background_"+trigger):
                histfile.cd("/"+analysisDir+"/dileptons_background_"+trigger)
                keylist=ROOT.gDirectory.GetListOfKeys()
                for key in keylist:
                    entry=key.GetName()
                    if not key.ReadObj().IsA().InheritsFrom("TTree"):
                        continue
                    # exclude certain entries we definitely do not want
                    if entry in ["decayLengthSignificance2D", # already done
                                 "bigTree", # different format
                                 "dPhitriggerCorr",
                                 "differentTrigObjects",
                                 "leptonQualityL",
                                 "leptonQualityH",
                                 "mass_calocorr", # masses taken from elsewhere
                                 "mass_corr",
                                 "mass_scalecorr",
                                 "mass_triggercorr",
                                 "validTracks",
                                 "validVertex"]: continue
                    if not entry in cutNames: cutNames.append(entry)
                    pass
                pass
            pass
        histfile.Close()
        pass


    if not do_all_plots:
        # only do the most important plots (for the PAS)
        cutNames=["trackerIsolationL","trackerIsolationH","dPhicorr"]
        pass
    
    for cutName in cutNames:
        treeOverlayPlot(SampleTriggerCombination,workdirs_background,workdirs_signal,
                        mcfactors,mcfactorsSignal,lepton_name,
                        analysis_directory+"/dileptons_background_TRIGGER/",cutName,
                        folder+"di"+lepton_name+"_"+cutName+".pdf",efficiencyListData,
                        efficiencyListMC,efficiencyListSignal)
        pass

    listOfCuts=[
        "eta_cand_pass",
        "numStandAloneMuons_pass",
        "leptonPtL_pass && leptonPtH_pass",
        "leptonEtaL_pass && leptonEtaH_pass",
        "leptonAbsD0significanceL_pass && leptonAbsD0significanceH_pass",
        "trackerIsolationL_pass && trackerIsolationH_pass",
        "oppositeCharge_pass",
        "vetoBackToBack_pass",
        "deltaRBetweenLeptons_pass",
        "numCaloMatches_pass",
        "numTrigMatches_pass",
        "differentTrigObjects_pass",
        "validTracks_pass && validVertex_pass && vertexChi2_pass",
        "dPhicorr_pass",
        "maxHitsBeforeVertex_pass",
        "decayLengthSignificance2D_pass",
        "passesAllCuts"
        ]
    if analysis_directory.find("mu")>=0:
        listOfCuts.remove("numCaloMatches_pass")
        pass

    print "==========efficiency output for data,",analysis_directory
    appliedCuts="eta_cand_pass"
    for additionalCut in listOfCuts:
        appliedCuts+=" && "+additionalCut
        sum=number_of_surviving_candidates(SampleTriggerCombination,"background",
                                           analysis_directory,appliedCuts)
        print "%-50s"%additionalCut,sum
        pass
    print "======================================"

    # now measure efficiency to find a matching supercluster
    # for candidates passing all cuts except trigger match
    if analysis_directory[0]=="e":
        print "==========efficiency to match supercluster in data"
        appliedCuts="eta_cand_pass"
        for additionalCut in listOfCuts:
            if additionalCut=="numTrigMatches_pass": continue
            if additionalCut=="numCaloMatches_pass": continue
            if additionalCut=="differentTrigObjects_pass": continue
            if additionalCut=="decayLengthSignificance2D_pass": continue
            if additionalCut.find("leptonAbsD0significance")>=0: continue
            if additionalCut=="dPhicorr_pass": continue
            if additionalCut=="passesAllCuts": continue
            appliedCuts+=" && "+additionalCut
            pass
        sum1=number_of_surviving_candidates(SampleTriggerCombination,
                                            "background",
                                            analysis_directory,appliedCuts)
        appliedCuts+=" && numCaloMatches_pass"
        sum2=number_of_surviving_candidates(SampleTriggerCombination,
                                            "background",
                                            analysis_directory,appliedCuts)
        print "fraction of candidates with matching supercluster:",\
              sum2,"/",sum1
        print "======================================"
        pass

    # find list of triggers for MC
    mctriggers=[]
    for [dataworkdir,trigger,lumi] in SampleTriggerCombination:
        mctrig=replacementTrigger[trigger]
        if not mctrig in mctriggers: mctriggers.append(mctrig)
    if len(mctriggers)==0:
        mctriggers.append("anyTrigger")
    if len(mctriggers)!=1:
        print "ERROR: more than one MC trigger requested."
        print "       This means using the same events several times."
        print "       And most likely the normalization is done per trigger."
        sys.exit(1)
    # construct corresponding sample/trigger combinations for background MC
    MCSampleTriggerCombination=[]
    for i in range(len(workdirs_background)):
        for trigger in mctriggers:
            MCSampleTriggerCombination.append([workdirs_background[i],trigger,mcfactors[i]])
 
    print "==========efficiency output for background MC,",analysis_directory
    appliedCuts="eta_cand_pass"
    for additionalCut in listOfCuts:
        appliedCuts+=" && "+additionalCut
        sum=number_of_surviving_candidates(MCSampleTriggerCombination,
                                           "background",
                                           analysis_directory,appliedCuts)
        print "%-50s"%additionalCut,sum
        pass
    print "======================================"

    # construct corresponding sample/trigger combinations for signal MC
    for i in range(len(workdirs_signal)):
        MCSampleTriggerCombination=[]
        for trigger in mctriggers:
            MCSampleTriggerCombination.append([workdirs_signal[i],trigger,1.0])
            pass
        print "==========efficiency output for "+workdirs_signal[i].split("/")[-1],",",analysis_directory
        appliedCuts="eta_cand_pass"
        for additionalCut in listOfCuts:
            appliedCuts+=" && "+additionalCut
            sum=number_of_surviving_candidates(MCSampleTriggerCombination,
                                               "signal",
                                               analysis_directory,appliedCuts)
            print "%-50s"%additionalCut,sum
            pass
        print "======================================"
        pass
    
    
    #---------------------
    # Efficiency Output
    #---------------------

    #Widths of each column
    col1Width = 26
    col2Width = 25
    col3Width = 30
    col4Width = 21

    #Print header for efficiency table
    print "\n"
    print "CUT EFFICIENCIES: DATA",lepton_name
    print " "*(col1Width - len("CutName"))                        +"CutName",
    print " "*(10 + col2Width - len("PassesThisAndAllOtherCuts"))+"PassesThisAndAllOtherCuts",
    print " "*(-6 + col3Width - len("PassesAllOtherCuts"))       +"PassesAllOtherCuts",
    print " "*(-4 + col4Width - len("Efficiency"))               +"Efficiency"
    print "#"*104

    #Print the efficiency data
    for data in efficiencyListData:

        #The string below contains a single row of the table with all elements alligned
        formattedStr = " "*(col1Width - len(data['cutName']))+str(data['cutName'])
        formattedStr+= " "*(col2Width - len(str(data['passesThisAndAllOtherCuts'])))
        formattedStr+= str(data['passesThisAndAllOtherCuts'])
        formattedStr+= " "*(col3Width - len(str(data['passesAllOtherCuts'])))+str(data['passesAllOtherCuts'])
        formattedStr+= " "*(col4Width - len(str(float("%.2e"%data['efficiency']))))
        formattedStr+= str(float("%.2e"%data['efficiency']))+"%"    
        print formattedStr
    
    #Print header for efficiency table
    print "\n"
    print "CUT EFFICIENCIES: BACKGROUND MC", lepton_name
    print " "*(col1Width - len("CutName"))                        +"CutName",
    print " "*(10 + col2Width - len("PassesThisAndAllOtherCuts"))+"PassesThisAndAllOtherCuts",
    print " "*(-6 + col3Width - len("PassesAllOtherCuts"))       +"PassesAllOtherCuts",
    print " "*(-4 + col4Width - len("Efficiency"))               +"Efficiency"
    print "#"*104

    #Print the efficiency data
    for data in efficiencyListMC:

        #The string below contains a single row of the table with all elements alligned
        formattedStr = " "*(col1Width - len(data['cutName']))+str(data['cutName'])
        formattedStr+= " "*(col2Width - len(str(data['passesThisAndAllOtherCuts'])))
        formattedStr+= str(data['passesThisAndAllOtherCuts'])
        formattedStr+= " "*(col3Width - len(str(data['passesAllOtherCuts'])))+str(data['passesAllOtherCuts'])
        formattedStr+= " "*(col4Width - len(str(float("%.2e"%data['efficiency']))))
        formattedStr+= str(float("%.2e"%data['efficiency']))+"%"    
        print formattedStr
    
    #Print header for efficiency table
    print "\n"
    print "CUT EFFICIENCIES: SIGNAL MC", lepton_name
    print " "*(col1Width - len("CutName"))                        +"CutName",
    print " "*(10 + col2Width - len("PassesThisAndAllOtherCuts"))+"PassesThisAndAllOtherCuts",
    print " "*(-6 + col3Width - len("PassesAllOtherCuts"))       +"PassesAllOtherCuts",
    print " "*(-4 + col4Width - len("Efficiency"))               +"Efficiency"
    print "#"*104

    #Print the efficiency data
    for data in efficiencyListSignal:

        #The string below contains a single row of the table with all elements alligned
        formattedStr = " "*(col1Width - len(data['cutName']))+str(data['cutName'])
        formattedStr+= " "*(col2Width - len(str(data['passesThisAndAllOtherCuts'])))
        formattedStr+= str(data['passesThisAndAllOtherCuts'])
        formattedStr+= " "*(col3Width - len(str(data['passesAllOtherCuts'])))+str(data['passesAllOtherCuts'])
        formattedStr+= " "*(col4Width - len(str(float("%.2e"%data['efficiency']))))
        formattedStr+= str(float("%.2e"%data['efficiency']))+"%"    
        print formattedStr
    
    return



# trigger efficiency measurements on unbiased sample
def triggerEfficiency(SampleTriggerCombination,treeDir):
    
    # Initialise the count of the number of events for the different triggers
    sumTrig = 0
    anyTrig = 0

    curDir = ""
    
    for [workdir,trigger,lumi] in SampleTriggerCombination:
        if curDir != workdir: #Check the file is not already open
            f = ROOT.TFile.Open(workdir + "/histograms.root")
            anyTrig += f.Get(treeDir+"/dileptons_background_anyTrigger/"
                             +"numTrigMatches").GetEntries("passesAllOtherCutsIgnoreLifetime")
            curDir = workdir
            pass
        
        sumTrig += f.Get(treeDir+"/dileptons_background_"+trigger+"/"
                         +"numTrigMatches").GetEntries("passesAllOtherCutsIgnoreLifetime")
        pass

    # Calculate the trigger efficiency
    if (anyTrig > 0):
        efficiency = 100.*float(sumTrig)/float(anyTrig)
        pass
    else:
        efficiency = 0
        pass
    
    # Print filename
    print "\nTrigger efficiency results\n"
   
    # Print header for efficiency table
    print " "*(8  - len("sumTrig"))    + "sumTrig",
    print " "*(6  - len("anyTig"))     + "anyTrig",
    print " "*(11 - len("Efficiency")) + "Efficiency"
    print "#"*28

    # Print the efficiency data
    outputStr =  " "*(8  - len(str(sumTrig)))+ str(sumTrig)
    outputStr+=  " "*(8  - len(str(anyTrig)))+ str(anyTrig)
    outputStr+=  " "*(11 - len(str(float("%.2e"%efficiency))))+ str(float("%.2e"%efficiency))+"%\n"
    print outputStr

    f.Close()
    return
    

#=======================================================================
#=======================================================================
#==== MAIN
#=======================================================================
#=======================================================================


    
#############################################
### PREPARATIONS
#############################################

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.Macro( 'HarderAnalysis/DisplacedDileptons/test/tdrstyle.C' )
#ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadLeftMargin(0.15);
ROOT.gStyle.SetPadRightMargin(0.05);
#ROOT.gStyle.SetTitleYOffset(1.5)


#############################################
### OBTAIN BACKGROUND MC SCALE FACTORS
#############################################

print "====================================="
print "EVALUATING MC SAMPLE WEIGHTS"
print "====================================="
mcweights_e=[]
for workdir in workdirs_background_e:
    mcweights_e.append(get_sample_weight(workdir))
mcweights_mu=[]
for workdir in workdirs_background_mu:
    mcweights_mu.append(get_sample_weight(workdir))
mcweights_sig=[]
for workdir in workdirs_signal:
    # signal cross-section unknown, but we define some value
    # for making signal drawn on top of background MC look good
    mcweights_sig.append(get_sample_weight(workdir))
print "\n"


##############################################
### LUMI+TRIGGER ANALYSIS FOR ALL DATA SAMPLES
##############################################

lumisum=0.0
eSampleTriggerCombination1=[]
eSampleTriggerCombination2=[]
eTriggerInMuSampleCombination1=[]
for workdir in workdirs_data_e:
    print "====================================="
    print "DATA:",workdir
    print "====================================="
    [lumi,best_triggers,best_trigger_lumi]=process_dataset(workdir,"electron",eAnalysis)
    # now put the sample/trigger combinations into categories that are being treated together
    for i in range(len(best_triggers)):
        if best_triggers[i].find("DoublePhoton33")>=0:
            eSampleTriggerCombination1.append([workdir,best_triggers[i],best_trigger_lumi[i]])
            eTriggerInMuSampleCombination1.append([workdirs_data_mu[workdirs_data_e.index(workdir)],best_triggers[i],best_trigger_lumi[i]])
        else:
            eSampleTriggerCombination2.append([workdir,best_triggers[i],best_trigger_lumi[i]])
    lumisum+=lumi
    print "\n"
    pass
print "------ total lumi processed in electron channel:",lumisum,"/pb --------\n\n"

lumisum=0.0
muSampleTriggerCombination1=[]
muSampleTriggerCombination2=[]
for workdir in workdirs_data_mu:
    print "====================================="
    print "DATA:",workdir
    print "====================================="
    [lumi,best_triggers,best_trigger_lumi]=process_dataset(workdir,"muon",muAnalysis)
    # now put the sample/trigger combinations into categories that are being treated together
    for i in range(len(best_triggers)):
        if best_triggers[i].find("HLT_L2DoubleMu23_NoVertex")>=0:
            muSampleTriggerCombination1.append([workdir,best_triggers[i],best_trigger_lumi[i]])
        else:
            muSampleTriggerCombination2.append([workdir,best_triggers[i],best_trigger_lumi[i]])
    lumisum+=lumi
    print "\n"
    pass
print "------ total lumi processed in muon channel:",lumisum,"/pb --------\n\n"


##############################################
### OUTPUT DIRECTORY STRUCTURE
##############################################

# Get current date to use for a folder name
t = datetime.datetime.now()
curDate = str(t.year) + "_" + str(t.month) + "_" + str(t.day)
plotfolder="Analysis/"+curDate
if os.path.exists(plotfolder):
    num=0
    while os.path.exists(plotfolder+"_rev"+str(num)): num+=1
    plotfolder+="_rev"+str(num)
    pass
# Check the analysis folder exists
if not os.path.isdir("Analysis"):
    os.mkdir("Analysis")
    pass
# Make a folder with the current date if it does not exist
if not os.path.isdir(plotfolder):
    os.mkdir(plotfolder)
    pass
# subfolder for all the stuff needed for limit calculation
limitfolder=plotfolder+"/limits"
if not os.path.isdir(limitfolder):
    os.mkdir(limitfolder)
# subfolder for various efficiency and benchmark plots
benchmarkfolder=plotfolder+"/benchmarks"
if not os.path.isdir(benchmarkfolder):
    os.mkdir(benchmarkfolder)


##############################################
### DATA/MC PLOTS FOR SELECTED DATA/TRIGGERS
##############################################

if len(eSampleTriggerCombination1)>0:
    makePlots(eSampleTriggerCombination1,workdirs_background_e,
              workdirs_signal,mcweights_e,mcweights_sig,
              "electron",eAnalysis,plotfolder,limitfolder)
    pass
if len(eSampleTriggerCombination2)>0:
    makePlots(eSampleTriggerCombination2,workdirs_background_e,
              workdirs_signal,mcweights_e,mcweights_sig,
              "electron",eAnalysis,plotfolder,limitfolder)
    pass
if len(muSampleTriggerCombination1)>0:
    makePlots(muSampleTriggerCombination1,workdirs_background_mu,
              workdirs_signal,mcweights_mu,mcweights_sig,
              "muon",muAnalysis,plotfolder,limitfolder)
    pass
if len(muSampleTriggerCombination2)>0:
    makePlots(muSampleTriggerCombination2,workdirs_background_mu,
              workdirs_signal,mcweights_mu,mcweights_sig,
              "muon",muAnalysis,plotfolder,limitfolder)
    pass


#############################################
### EFFICIENCY PLOTS FROM MC
#############################################

# tracking efficiency for electrons
for workdir in workdirs_benchmark_e:
    sampleID=workdir.split("/")[-1].replace("_analysis","")
    efficiencyPlot(workdir,eAnalysis+"/leptons/trueLeptonRadWithTrack",
                   eAnalysis+"/gen/leptonProdVtxRadius2D",
                   benchmarkfolder+"/track_effi_electrons.pdf",
                   "track reconstruction efficiency as function of impact parameter",
                   "d_{0} [cm]")
    pass

# tracking efficiency for muons
for workdir in workdirs_benchmark_mu:
    sampleID=workdir.split("/")[-1].replace("_analysis","")
    efficiencyPlot(workdir,muAnalysis+"/leptons/trueLeptonRadWithTrack",
                   muAnalysis+"/gen/leptonProdVtxRadius2D",
                   benchmarkfolder+"/track_effi_muons.pdf",
                   "track reconstruction efficiency as function of impact parameter",
                   "d_{0} [cm]")
    pass

# displaced lepton efficiency before and after rereco
for workdir in workdirs_benchmark_e:
    sampleID=workdir.split("/")[-1].replace("_analysis","")
    efficiencyPlot(workdir,
                   eAnalysis+"/electrons/trueProdVtxRad",
                   eAnalysis+"/gen/leptonProdVtxRadius2D",
                   benchmarkfolder+"/electron_effi_"+sampleID+".pdf",
                   "electron reconstruction efficiency as function of impact parameter",
                   "d_{0} [cm]")
    pass
for workdir in workdirs_benchmark_mu:
    sampleID=workdir.split("/")[-1].replace("_analysis","")
    efficiencyPlot(workdir,
                   muAnalysis+"/muons/trueProdVtxRad",
                   muAnalysis+"/gen/leptonProdVtxRadius2D",
                   benchmarkfolder+"/muon_effi_"+sampleID+".pdf",
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
                   benchmarkfolder+"/dimuon1_effi_"+sampleID+".pdf",
                   "dimuon reconstruction efficiency"+
                   " as function of decay length, 1 dimuon",
                   "L_{xy} [cm]")

    efficiencyPlot(workdir,
                   muAnalysis+"/dileptons/trueDecayLength2D_2",
                   "isoTrackPrefilter/decayLength2D_2muon",
                   benchmarkfolder+"/dimuon2_effi_"+sampleID+".pdf",
                   "dimuon reconstruction efficiency"+
                   " as function of decay length, 2 dimuons",
                   "L_{xy} [cm]")
    
    efficiencyPlot(workdir,
                   eAnalysis+"/dileptons/trueDecayLength2D_1",
                   "isoTrackPrefilter/decayLength2D_1elec",
                   benchmarkfolder+"/dielectron1_effi_"+sampleID+".pdf",
                   "dielectron reconstruction efficiency"+
                   " as function of decay length, 1 dielectron",
                   "L_{xy} [cm]")
    
    efficiencyPlot(workdir,
                   eAnalysis+"/dileptons/trueDecayLength2D_2",
                   "isoTrackPrefilter/decayLength2D_2elec",
                   benchmarkfolder+"/dielectron2_effi_"+sampleID+".pdf",
                   "dielectron reconstruction efficiency"+
                   " as function of decay length, 2 dielectrons",
                   "L_{xy} [cm]")
    pass

# trigger efficiency in MC
mctriggers=[]
for datatrigger,mctrigger in replacementTrigger.iteritems():
    if not mctrigger in mctriggers: mctriggers.append(mctrigger)
for workdir in workdirs_signal:
    [edenom,edenom_loose,dm1,dm2]=get_histogram(workdir,eAnalysis\
                                   +"/dileptons_signal_anyTrigger/numTrigMatches",
                                   1,-10,10,0,3000)
    [mudenom,mudenom_loose,dm1,dm2]=get_histogram(workdir,muAnalysis\
                                    +"/dileptons_signal_anyTrigger/numTrigMatches",
                                    1,-10,10,0,3000)
    for mctrigger in mctriggers:
        if mctrigger.find("Photon")>=0 and edenom.Integral()>0:
            [enum,enum_loose,em1,em2]=get_histogram(workdir,eAnalysis\
                                         +"/dileptons_signal_%s/numTrigMatches"%mctrigger,
                                         1,-10,10,0,3000)
            print "MC trigger efficiency of",mctrigger,"in",eAnalysis,\
                  "channel of",workdir.split("/")[-1],\
                  ":",enum.Integral()/edenom.Integral()
        elif mctrigger.find("Mu")>=0 and mudenom.Integral()>0:
            [enum,enum_loose,em1,em2]=get_histogram(workdir,muAnalysis\
                                         +"/dileptons_signal_%s/numTrigMatches"%mctrigger,
                                         1,-10,10,0,3000)
            print "MC trigger efficiency of",mctrigger,"in",muAnalysis,\
                  "channel of",workdir.split("/")[-1],\
                  ":",enum.Integral()/mudenom.Integral()
            pass
        pass
    pass

# ctau factors for lifetime reweighting
# here we use k=3 and k=1./3
ctfact_list=[1./30.,1./20.,1./10.,1./3.,0.5,2.0,3.,10.,20.,30.]


# corrected efficiencies, systematics
mutrack_effi1=[]
mutrack_effi2=[]
etrack_effi1=[]
etrack_effi2=[]
mutrack_effi1_pileupuncertainty=[]
mutrack_effi2_pileupuncertainty=[]
etrack_effi1_pileupuncertainty=[]
etrack_effi2_pileupuncertainty=[]
mutrack_effi1_statisticaluncertainty=[]
mutrack_effi2_statisticaluncertainty=[]
etrack_effi1_statisticaluncertainty=[]
etrack_effi2_statisticaluncertainty=[]
mutrack_effi1_uncorrected=[]
mutrack_effi2_uncorrected=[]
etrack_effi1_uncorrected=[]
etrack_effi2_uncorrected=[]
mutrack_effi1_lifetime=[]
mutrack_effi2_lifetime=[]
etrack_effi1_lifetime=[]
etrack_effi2_lifetime=[]
mutrack_effi1_lifetime_statuncert=[]
mutrack_effi2_lifetime_statuncert=[]
etrack_effi1_lifetime_statuncert=[]
etrack_effi2_lifetime_statuncert=[]
for entry in ctfact_list:
    mutrack_effi1_lifetime.append([])
    mutrack_effi2_lifetime.append([])
    etrack_effi1_lifetime.append([])
    etrack_effi2_lifetime.append([])
    mutrack_effi1_lifetime_statuncert.append([])
    mutrack_effi2_lifetime_statuncert.append([])
    etrack_effi1_lifetime_statuncert.append([])
    etrack_effi2_lifetime_statuncert.append([])
    pass

for treename in\
        eAnalysis+"/dileptons_signal_HLT_DoublePhoton33_v2/phi_cand",\
        muAnalysis+"/dileptons_signal_HLT_L2DoubleMu23_NoVertex_v1/phi_cand":
    for numDecays in [1,2]:
        for workdir in workdirs_signal:
            # extract ctau from DBS name
            ctau=0
            dbsname=os.popen("grep sampleDataSet "+workdir+"/*_cff.py").readlines()
            for line in dbsname:
                for field in line.split("_"):
                    if field.find("CTau-")==0:
                        ctau=int(field.split("-")[1])
                        pass
                    pass
                pass
            ctau/=10. # convert from Pythia units (mm) to cm
            # extract efficiency
            sampleID=workdir.split("/")[-1].replace("_analysis","")
            enumfile=ROOT.TFile.Open(workdir+"/histograms.root")
            enumtree=enumfile.Get(treename)
            # uncorrected, and central value with pile-up correction
            enumtree.Project("valhist","value",\
                             #"weight*(passesAllCuts"+\
                             "(0.5*(weight_up+weight_down))*(passesAllCuts"+\
                             " && numDecays==%i)"%numDecays)
            enum_uncorrected=enumfile.Get("valhist").GetEntries()
            enum_central=enumfile.Get("valhist").Integral()
            try:
                enum_central_relerr=1./math.sqrt(enumfile.Get("valhist").GetEffectiveEntries())
            except:
                enum_central_relerr=0
                pass
            # pile-up variation up and down
            enumtree.Project("valhist","value",
                             "weight_up*(passesAllCuts"+\
                             " && numDecays==%i)"%numDecays)
            enum_up=enumfile.Get("valhist").Integral()
            enumtree.Project("valhist","value","weight_down*(passesAllCuts"+\
                             " && numDecays==%i)"%numDecays)
            enum_down=enumfile.Get("valhist").Integral()
            # lifetime variation up and down
            # weight (from Ian): (1/newCtau)*exp(-ctau/newCtau)/(1/oldCtau)/exp(-ctau/oldCtau)
            # if we reweight by a specific factor k, this gives
            # (1./k)*exp(-ctau/k/oldCtau)/exp(-ctau/oldCtau)
            # =(1./k)*exp(ctau*(1./oldCtau-1./k/oldCtau))
            # =(1./k)*exp(ctau*(1./k)*(k-1.)/oldCtau)
            enum_list=[]
            enum_relerr_list=[]
            for ctfact in ctfact_list:
                project="(0.5*(weight_up+weight_down)"+\
                         "*1./%f*exp(-ctau1/%f/%f)"%(ctfact,ctfact,ctau)+\
                         "/exp(-ctau1/%f)"%ctau+\
                         "*1./%f*exp(-ctau2/%f/%f)"%(ctfact,ctfact,ctau)+\
                         "/exp(-ctau2/%f))"%ctau+\
                         "*(passesAllCuts"+\
                         " && numDecays==%i)"%numDecays
                enumtree.Project("valhist","value",project)
                enum_list.append(enumfile.Get("valhist").Integral())
                try:
                    enum_relerr_list.append(1./math.sqrt(enumfile.Get("valhist").GetEffectiveEntries()))
                except:
                    enum_relerr_list.append(0.0)
                    pass
                pass
            enumfile.Close()
            
            # now open prefilter file to get normalization
            # (numSignal in proper decay channel)
            denomfile=ROOT.TFile.Open(workdir+"/prefilter.root")
            if treename.find("eTrack")>=0:
                denomhist=denomfile.Get("isoTrackPrefilter/numSignalE")
            elif treename.find("muTrack")>=0:
                denomhist=denomfile.Get("isoTrackPrefilter/numSignalMu")
                pass
            denom_uncorrected=denomhist.GetBinContent(numDecays+1)*numDecays

            # and get overall normalization correction due to lumi weights
            # from prefilter.root
            pileuphist=denomfile.Get("isoTrackPrefilter/pileup_3bx")

            pileup_weights_MC=[ 0.104109,
                                0.0703573,
                                0.0698445,
                                0.0698254,
                                0.0697054,
                                0.0697907,
                                0.0696751,
                                0.0694486,
                                0.0680332,
                                0.0651044,
                                0.0598036,
                                0.0527395,
                                0.0439513,
                                0.0352202,
                                0.0266714,
                                0.019411,
                                0.0133974,
                                0.00898536,
                                0.0057516,
                                0.00351493,
                                0.00212087,
                                0.00122891,
                                0.00070592,
                                0.000384744,
                                0.000219377 ]
            pileup_data = [ 0.019091,
                            0.0293974,
                            0.0667931,
                            0.108859,
                            0.139533,
                            0.149342,
                            0.138629,
                            0.114582,
                            0.0859364,
                            0.059324,
                            0.0381123,
                            0.0229881,
                            0.0131129,
                            0.00711764,
                            0.00369635,
                            0.00184543,
                            0.000889604,
                            0.000415683,
                            0.000188921,
                            0.000146288,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0 ]
            # basis of calculation for pile-up shifting weights
            p1_minus = [
                -0.677786,
                -0.619614,
                -0.49465,
                -0.357963,
                -0.238359,
                -0.110002,
                0.0348629,
                0.191263,
                0.347648,
                0.516615,
                0.679646,
                0.836673,
                0.97764,
                1.135,
                1.29922,
                1.42467,
                1.55901,
                1.61762,
                1.67275,
                1.96008
                ]
            p2_minus = [
                0.526164,
                0.251816,
                0.11049,
                0.026917,
                -0.0464692,
                -0.087022,
                -0.0931581,
                -0.0714295,
                -0.0331772,
                0.0347473,
                0.108658,
                0.193048,
                0.272314,
                0.376357,
                0.4964,
                0.58854,
                0.684959,
                0.731063,
                0.760044,
                1.02386
                ]
            p1_plus = [
                -0.739059,
                -0.594445,
                -0.477276,
                -0.359707,
                -0.233573,
                -0.103458,
                0.0373401,
                0.176571,
                0.337617,
                0.499074,
                0.675126,
                0.840522,
                1.00917,
                1.15847,
                1.23816,
                1.44271,
                1.52982,
                1.46385,
                1.5802,
                0.988689
                ]
            p2_plus = [
                0.208068,
                0.130033,
                0.0850356,
                0.0448344,
                0.000749832,
                -0.0331347,
                -0.0653281,
                -0.0746009,
                -0.0800667,
                -0.0527636,
                -0.00402649,
                0.103338,
                0.261261,
                0.491084,
                0.857966,
                1.19495,
                1.75071,
                2.65559,
                3.35433,
                5.48835
                ]
            shift=0.6
            pweight_up=[]
            pweight_down=[]
            for i in range(len(p1_plus)):
                pweight_up.append(1+p1_plus[i]*shift+p2_plus[i]*shift*shift)
                pweight_down.append(1-p1_minus[i]*shift+p2_minus[i]*shift*shift)
                pass
            for i in range(5):
                # not quite correct, but these bins are hardly populated anyway
                pweight_up.append(1.0)
                pweight_down.append(1.0)
                pass

            # note: denominator for lifetime reweighting is simply whatever
            # we get with central lumi weights. ctau weights should preserve
            # normalization.
            
            # now calculate the weighted number of events
            sum=0.0
            sumw=0.0
            sumw2=0.0
            sumw_up=0.0
            sumw_down=0.0
            for i in range(pileuphist.GetNbinsX()):
                if i<len(pileup_data) and i<len(pileup_weights_MC):
                    w=pileup_data[i]/pileup_weights_MC[i]
                    n=pileuphist.GetBinContent(i+1)
                    sumw+=w*n
                    sumw2+=w*w*n
                    sumw_up+=pweight_up[i]*w*n
                    sumw_down+=pweight_down[i]*w*n
                    sum+=n
                    pass
                pass
            denom_central=denom_uncorrected*sumw/sum
            denom_up=denom_uncorrected*sumw_up/sum
            denom_down=denom_uncorrected*sumw_down/sum
            effi_uncorrected="%-70s %10.4f"%(sampleID,enum_uncorrected/denom_uncorrected)
            effi=enum_central/denom_central
            effi_corrected="%-70s %10.4f"%(sampleID,enum_central/denom_central)
            effi_e1=(enum_up/denom_up)-(enum_central/denom_central)
            effi_e2=(enum_down/denom_down)-(enum_central/denom_central)
            effi_pileuperror="%-50s %10.4f %10.4f"%(sampleID,
                                                    abs(max(effi_e1,effi_e2)),
                                                    -abs(min(effi_e1,effi_e2)))
            effi_statisticalerror="%-50s %10.4f %10.4f"%(sampleID,
                                                         effi*enum_central_relerr,
                                                         -effi*enum_central_relerr)

            effi_lifetime=[]
            effi_lifetime_staterror=[]
            for i in range(len(enum_list)):
                enum=enum_list[i]
                enum_relerr=enum_relerr_list[i]
                effi=enum/denom_central
                effi_lifetime.append("%-70s %10.4f"%(sampleID,effi))
                effi_lifetime_staterror.append("%-50s %10.4f %10.4f"\
                                               %(sampleID,
                                                 effi*enum_relerr,
                                                 -effi*enum_relerr))
                pass
            
            if treename[0]=="e" and numDecays==1:
                etrack_effi1.append(effi_corrected)
                etrack_effi1_uncorrected.append(effi_uncorrected)
                for i in range(len(effi_lifetime)):
                    etrack_effi1_lifetime[i].append(effi_lifetime[i])
                    etrack_effi1_lifetime_statuncert[i].append(effi_lifetime_staterror[i])
                    pass
                etrack_effi1_pileupuncertainty.append(effi_pileuperror)
                etrack_effi1_statisticaluncertainty.append(effi_statisticalerror)
            elif treename[0]=="e" and numDecays==2:
                etrack_effi2.append(effi_corrected)
                etrack_effi2_uncorrected.append(effi_uncorrected)
                for i in range(len(effi_lifetime)):
                    etrack_effi2_lifetime[i].append(effi_lifetime[i])
                    etrack_effi2_lifetime_statuncert[i].append(effi_lifetime_staterror[i])
                    pass
                etrack_effi2_pileupuncertainty.append(effi_pileuperror)
                etrack_effi2_statisticaluncertainty.append(effi_statisticalerror)
            elif treename[0]=="m" and numDecays==1:
                mutrack_effi1.append(effi_corrected)
                mutrack_effi1_uncorrected.append(effi_uncorrected)
                for i in range(len(effi_lifetime)):
                    mutrack_effi1_lifetime[i].append(effi_lifetime[i])
                    mutrack_effi1_lifetime_statuncert[i].append(effi_lifetime_staterror[i])
                    pass
                mutrack_effi1_pileupuncertainty.append(effi_pileuperror)
                mutrack_effi1_statisticaluncertainty.append(effi_statisticalerror)
            elif treename[0]=="m" and numDecays==2:
                mutrack_effi2.append(effi_corrected)
                mutrack_effi2_uncorrected.append(effi_uncorrected)
                for i in range(len(effi_lifetime)):
                    mutrack_effi2_lifetime[i].append(effi_lifetime[i])
                    mutrack_effi2_lifetime_statuncert[i].append(effi_lifetime_staterror[i])
                    pass
                mutrack_effi2_pileupuncertainty.append(effi_pileuperror)
                mutrack_effi2_statisticaluncertainty.append(effi_statisticalerror)
                pass
            pass
        pass
    pass

mu1file=open(limitfolder+"/efficiencies_dimuon1.txt","w")
for entry in mutrack_effi1: mu1file.write("%s\n"%entry)
mu1file.close()
mu2file=open(limitfolder+"/efficiencies_dimuon2.txt","w")
for entry in mutrack_effi2: mu2file.write("%s\n"%entry)
mu2file.close()
e1file=open(limitfolder+"/efficiencies_dielectron1.txt","w")
for entry in etrack_effi1: e1file.write("%s\n"%entry)
e1file.close()
e2file=open(limitfolder+"/efficiencies_dielectron2.txt","w")
for entry in etrack_effi2: e2file.write("%s\n"%entry)
e2file.close()

for i in range(len(ctfact_list)):
    filespec="ctaufact_%.3f"%ctfact_list[i]
    mu1file=open(limitfolder+"/efficiencies_dimuon1_"+filespec+".txt","w")
    for entry in mutrack_effi1_lifetime[i]: mu1file.write("%s\n"%entry)
    mu1file.close()
    mu2file=open(limitfolder+"/efficiencies_dimuon2_"+filespec+".txt","w")
    for entry in mutrack_effi2_lifetime[i]: mu2file.write("%s\n"%entry)
    mu2file.close()
    e1file=open(limitfolder+"/efficiencies_dielectron1_"+filespec+".txt","w")
    for entry in etrack_effi1_lifetime[i]: e1file.write("%s\n"%entry)
    e1file.close()
    e2file=open(limitfolder+"/efficiencies_dielectron2_"+filespec+".txt","w")
    for entry in etrack_effi2_lifetime[i]: e2file.write("%s\n"%entry)
    e2file.close()

    mu1file=open(limitfolder+"/efficiencies_dimuon1_"+filespec+"_statistical_uncertainty.txt","w")
    for entry in mutrack_effi1_lifetime_statuncert[i]: mu1file.write("%s\n"%entry)
    mu1file.close()
    mu2file=open(limitfolder+"/efficiencies_dimuon2_"+filespec+"_statistical_uncertainty.txt","w")
    for entry in mutrack_effi2_lifetime_statuncert[i]: mu2file.write("%s\n"%entry)
    mu2file.close()
    e1file=open(limitfolder+"/efficiencies_dielectron1_"+filespec+"_statistical_uncertainty.txt","w")
    for entry in etrack_effi1_lifetime_statuncert[i]: e1file.write("%s\n"%entry)
    e1file.close()
    e2file=open(limitfolder+"/efficiencies_dielectron2_"+filespec+"_statistical_uncertainty.txt","w")
    for entry in etrack_effi2_lifetime_statuncert[i]: e2file.write("%s\n"%entry)
    e2file.close()
    pass

mu1file=open(limitfolder+"/efficiencies_dimuon1_uncorrected.txt","w")
for entry in mutrack_effi1_uncorrected: mu1file.write("%s\n"%entry)
mu1file.close()
mu2file=open(limitfolder+"/efficiencies_dimuon2_uncorrected.txt","w")
for entry in mutrack_effi2_uncorrected: mu2file.write("%s\n"%entry)
mu2file.close()
e1file=open(limitfolder+"/efficiencies_dielectron1_uncorrected.txt","w")
for entry in etrack_effi1_uncorrected: e1file.write("%s\n"%entry)
e1file.close()
e2file=open(limitfolder+"/efficiencies_dielectron2_uncorrected.txt","w")
for entry in etrack_effi2_uncorrected: e2file.write("%s\n"%entry)
e2file.close()

mu1file=open(limitfolder+"/efficiencies_dimuon1_pileup_uncertainty.txt","w")
for entry in mutrack_effi1_pileupuncertainty: mu1file.write("%s\n"%entry)
mu1file.close()
mu2file=open(limitfolder+"/efficiencies_dimuon2_pileup_uncertainty.txt","w")
for entry in mutrack_effi2_pileupuncertainty: mu2file.write("%s\n"%entry)
mu2file.close()
e1file=open(limitfolder+"/efficiencies_dielectron1_pileup_uncertainty.txt","w")
for entry in etrack_effi1_pileupuncertainty: e1file.write("%s\n"%entry)
e1file.close()
e2file=open(limitfolder+"/efficiencies_dielectron2_pileup_uncertainty.txt","w")
for entry in etrack_effi2_pileupuncertainty: e2file.write("%s\n"%entry)
e2file.close()

mu1file=open(limitfolder+"/efficiencies_dimuon1_statistical_uncertainty.txt","w")
for entry in mutrack_effi1_statisticaluncertainty: mu1file.write("%s\n"%entry)
mu1file.close()
mu2file=open(limitfolder+"/efficiencies_dimuon2_statistical_uncertainty.txt","w")
for entry in mutrack_effi2_statisticaluncertainty: mu2file.write("%s\n"%entry)
mu2file.close()
e1file=open(limitfolder+"/efficiencies_dielectron1_statistical_uncertainty.txt","w")
for entry in etrack_effi1_statisticaluncertainty: e1file.write("%s\n"%entry)
e1file.close()
e2file=open(limitfolder+"/efficiencies_dielectron2_statistical_uncertainty.txt","w")
for entry in etrack_effi2_statisticaluncertainty: e2file.write("%s\n"%entry)
e2file.close()

    
# signal shape for limit setting code
for workdir in workdirs_signal:
    sampleID=workdir.split("/")[-1].replace("_analysis","")
    masspeak=float(sampleID.split("_")[2].replace("F","").replace("L",""))
    lowerbound=0.5*masspeak
    if lowerbound<50: lowerbound=0
    upperbound=1.5*masspeak
    if upperbound<100: upperbound=2*masspeak
    histfile=ROOT.TFile.Open(workdir+"/histograms.root")
    [valhist,valhist_loose,masshist,rejecthist]=get_histogram(workdir,
                                                muAnalysis+"/dileptons_signal_HLT_L2DoubleMu23_NoVertex_v1/mass_corr",100,0,100,lowerbound,upperbound)
    masshist.SetTitle("dimuon mass distribution, signal")
    masshist.SetName("dileptonmass_signal")
    outfile=ROOT.TFile.Open(limitfolder+"/signal_shape_muon_"+sampleID+".root","RECREATE")
    masshist.Write()
    outfile.Close()
    [valhist,valhist_loose,masshist,rejecthist]=get_histogram(workdir,
                                                eAnalysis+"/dileptons_signal_HLT_DoublePhoton33_v2/mass_corr",100,0,100,lowerbound,upperbound)
    masshist.SetTitle("dielectron mass distribution, signal")
    masshist.SetName("dileptonmass_signal")
    outfile=ROOT.TFile.Open(limitfolder+"/signal_shape_electron_"+sampleID+".root","RECREATE")
    masshist.Write()
    outfile.Close()


    histfile.Close()


    
#############################################
### SPECIFIC STUDIES
#############################################

# deltaR between muons. potential problem region for trigger at small deltaR
for workdir in workdirs_signal:
    histfile=ROOT.TFile.Open(workdir+"/histograms.root")
    deltaRtree=histfile.Get(muAnalysis+"/dileptons_signal_anyTrigger/deltaRBetweenLeptons")
    canv=ROOT.TCanvas()
    deltaRtree.Draw("value","passesAllOtherCutsIgnoreLifetime")
    hist=0
    for item in canv.GetListOfPrimitives():
        if item.IsA().InheritsFrom("TH1"): hist=item
        pass
    if hist:
        hist.GetXaxis().SetTitle("deltaR")
        title="#Delta R btw leptons, "+workdir.split("/")[-1]
        title=title[:title.find("_analysis")]
        hist.SetTitle(title)
        pass
    CMSPrelim(muAnalysis,1)
    canv.Update()
    canv.Print(benchmarkfolder+"/deltaR_"+workdir.split("/")[-1]+".pdf")
    histfile.Close()
    pass
    

# track isolation distributions
for workdir in workdirs_signal:
    histfile=ROOT.TFile.Open(workdir+"/histograms.root")
    for analysisDir in [eAnalysis,muAnalysis]:
        isolationTree=histfile.Get(analysisDir+"/dileptons_signal_anyTrigger/trackerIsolationL")
        if isolationTree:
            canv=ROOT.TCanvas()
            canv.SetLogy()
            isolationTree.Draw("value","passesAllOtherCuts && value<100","h")
            hist=0
            for item in canv.GetListOfPrimitives():
                if item.IsA().InheritsFrom("TH1"): hist=item
                pass
            if hist:
                hist.GetXaxis().SetTitle("#Sigma p_{t} [GeV/c]")
                title="track #Sigma p_{t} within a #Delta R < 0.3 cone, "+workdir.split("/")[-1]
                title=title[:title.find("_analysis")]
                if analysisDir==eAnalysis:
                    title+=", electron channel"
                else:
                    title+=", muon channel"
                    pass
                hist.SetTitle(title)
                pass
            CMSPrelim(analysisDir,1)
            canv.Update()
            canv.Print(benchmarkfolder+"/isolation_"+analysisDir+"_"+workdir.split("/")[-1]+".pdf")
            pass
        pass
    histfile.Close()
    pass
    
# track isolation dependence on number of pile-up events
for workdir in workdirs_signal:
    histfile=ROOT.TFile.Open(workdir+"/histograms.root")
    for analysisDir in [eAnalysis,muAnalysis]:
        isoHist=histfile.Get(analysisDir+"/dileptons/isolation_vs_pileup_signal")
        if not isoHist: continue
        canv=ROOT.TCanvas()
        hist=isoHist.ProfileY()
        hist.SetMarkerStyle(20)
        hist.SetMinimum(0)
        hist.GetXaxis().SetTitle("pile-up vertices")
        title="average isolation #Sigma p_{t} vs pile-up, "
        if analysisDir==eAnalysis:
            title+="electron channel"
        else:
            title+="muon channel"
            pass
        hist.SetTitle(title)
        hist.Draw()
        canv.Update()
        CMSPrelim(analysisDir,1)
        canv.Print(benchmarkfolder+"/isolationVsPileup_"+analysisDir+"_"+workdir.split("/")[-1]+".pdf")
        pass
    histfile.Close()
    pass

#############################################
### PHOTON TRIGGER EFFICIENCY MEASUREMENT
#############################################

#triggerEfficiency(eTriggerInMuSampleCombination1,"/eTrackAnalysis")

#############################################
### RESOLUTION+BIAS PLOTS FROM MC
#############################################

# currently not running muAnalysis!
sys.exit(0)

# muon momentum bias
for workdir in workdirs_benchmark_mu:
    histfile=ROOT.TFile.Open(workdir+"/histograms.root")
    # note these plots are not available in the mutrack histograms
    pt_vs_radius_standalone=histfile.Get("muAnalysis/leptons/pt_vs_radius_standalone")
    pt_vs_radius_matched=histfile.Get("muAnalysis/leptons/pt_vs_radius_allothers")
    canv=ROOT.TCanvas()
    pt_vs_radius_standalone.SetMarkerColor(ROOT.kBlue)
    if workdir.find("stdRECO")>0:
        pt_vs_radius_standalone.SetTitle("reconstructed p_t vs impact parameter, standard RECO")
    else:
        pt_vs_radius_standalone.SetTitle("reconstructed p_t vs impact parameter, no IP constraint")
        pass
    pt_vs_radius_standalone.GetXaxis().SetTitle("d_{xy} [cm]")
    pt_vs_radius_standalone.GetYaxis().SetTitle("p_{t} [GeV/c]")
    pt_vs_radius_standalone.Draw()
    line=ROOT.TLine(0,50,100,50)
    line.SetLineWidth(3)
    line.Draw("same")
    pt_vs_radius_matched.SetMarkerColor(ROOT.kRed)
    pt_vs_radius_matched.Draw("same")
    CMSPrelim(muAnalysis,1)
    canv.Update()
    canv.Print(benchmarkfolder+"/muon_momentum_vs_radius_"+workdir.split("/")[-1]+".pdf")
    canv.Close()
    pass
