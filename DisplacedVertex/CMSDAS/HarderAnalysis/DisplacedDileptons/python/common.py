#!/bin/env python

import ROOT
import time
import os,sys,math

##############################################
### OUTPUT DIRECTORY STRUCTURE
##############################################

def set_ana_folder():

    # Get current date to use for a folder name
    plotfolder="Analysis/"+time.strftime("%Y%m%d")
    if os.path.exists(plotfolder):
        num=0
        while os.path.exists(plotfolder+"_rev%03i"%num): num+=1
        plotfolder+="_rev%03i"%num
        pass
    # Check the analysis folder exists
    if not os.path.isdir("Analysis"):
        os.mkdir("Analysis")
        pass
    # Make a folder with the current date if it does not exist
    if not os.path.isdir(plotfolder):
        os.mkdir(plotfolder)
        pass
    # subfolder for various efficiency and benchmark plots
    benchmarkfolder=plotfolder+"/benchmarks"
    if not os.path.isdir(benchmarkfolder):
        os.mkdir(benchmarkfolder)
    output = [benchmarkfolder, plotfolder]
    return output


def CMSPlotDecoration(channel,lumisum=0):
    info=ROOT.TLatex()
    info.SetNDC()
    if lumisum>0:
        info.DrawLatex(0.18,0.95,\
                       "CMS Preliminary #sqrt{s}=7 TeV L=%3.1f fb^{-1}"%(lumisum/1000.))
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
### DEFINE HISTOGRAM PROPERTIES
#############################################

    
def setAxisRange():
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
    return axisRange
    
def setHistNBins():
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
    return histNBins


def setHistXLegend():
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
    return histXLegend

    
def setHistTitle():
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
    return histTitle


def setAxisTitle():    
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
    return axisTitle
    
def setLegendEntry():
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
    return legendEntry
        
def setHistColour():
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
    return histColour


############################################################################
### DEFINE REPLACEMENT TRIGGERS TO BE USED WHEN ACTUAL TRIGGER MISSING IN MC
############################################################################

def replace_trig_name():
    replacementTrigger = {}
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

    return replacementTrigger


# tools: draw overflow bin in histogram

def drawhist(h,opt):
    draw_overflow=1
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
