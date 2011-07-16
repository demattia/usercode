#!/bin/bash

SIMDIR='/home/demattia/TrackingEfficiencyFromCosmics/Simulation/MoreStatistics'
DATADIR='/home/demattia/TrackingEfficiencyFromCosmics/Data'

# Produce efficiency plots
cp EfficiencyAnalyzer_cfg.py $SIMDIR
cd $SIMDIR
cmsRun EfficiencyAnalyzer_cfg.py
cd -
root -l -b -q 'SetStyle.C("'$SIMDIR'/", "Sim")'

cp EfficiencyAnalyzer_cfg.py $DATADIR
cd $DATADIR
cmsRun EfficiencyAnalyzer_cfg.py
cd -
root -l -b -q 'SetStyle.C("'$DATADIR'/", "Data")'

# Restyle control plots
mkdir -p plots_Sim
root -l -b -q ProduceStyledPlots.C'("'${SIMDIR}'/TrackingEfficiencyFromCosmics.root", "Sim")'
mkdir -p plots_Data
root -l -b -q ProduceStyledPlots.C'("'${DATADIR}'/TrackingEfficiencyFromCosmics.root", "Data")'

# Produce comparison plots (with the proper style)
mkdir -p comparedPlots
root -l -b -q 'ComparePlots.C("'${DATADIR}'/TrackingEfficiencyFromCosmics.root", "'${SIMDIR}'/TrackingEfficiencyFromCosmics.root")'
root -l -b -q 'CompareDxy.C("'${DATADIR}'/TrackingEfficiencyFromCosmics.root", "'${SIMDIR}'/TrackingEfficiencyFromCosmics.root")'

# Efficiency ratios
cat EfficiencyRatioProducer_cfg.py | sed s@SIMDIR@${SIMDIR}@ | sed s@DATADIR@${DATADIR}@ > EffRatioProd_cfg.py
cmsRun EffRatioProd_cfg.py
# Restyle efficiency ratio plot
root -l -b -q 'SetEffRatioStyle.C("", "")'


# Move all the plots in a subdir
mkdir -p plotsToCopy
cp styled_EffVsDxyRatio.p*                            plotsToCopy
cp styled_EffVsDxy_Sim.p*                             plotsToCopy
cp styled_EffVsDxy_Data.p*                            plotsToCopy

cp plots_Sim/standAloneToGenDeltaDxy.p*               plotsToCopy
cp plots_Sim/cleanedStandAloneToGenDeltaDxy.p*        plotsToCopy

cp comparedPlots/standAloneMuons_pt.p*                plotsToCopy
cp comparedPlots/standAloneMuons_eta.p*               plotsToCopy
cp comparedPlots/standAloneMuons_phi.p*               plotsToCopy
cp comparedPlots/standAloneMuons_dxy.p*               plotsToCopy
cp comparedPlots/standAloneMuons_dz.p*                plotsToCopy
cp comparedPlots/standAloneMuons_chi2.p*              plotsToCopy
cp comparedPlots/standAloneMuons_NValidHits.p*        plotsToCopy

cp comparedPlots/standAloneDelta_DeltaPhi.pdf         plotsToCopy
cp comparedPlots/standAloneDelta_DeltaDxy.pdf         plotsToCopy
cp comparedPlots/standAloneDelta_DeltaDz.pdf          plotsToCopy

cp comparedPlots/cleanedStandAloneMuons_pt.p*         plotsToCopy
cp comparedPlots/cleanedStandAloneMuons_eta.p*        plotsToCopy
cp comparedPlots/cleanedStandAloneMuons_phi.p*        plotsToCopy
cp comparedPlots/cleanedStandAloneMuons_dxy.p*        plotsToCopy
cp comparedPlots/cleanedStandAloneMuons_dz.p*         plotsToCopy
cp comparedPlots/cleanedStandAloneMuons_chi2.p*       plotsToCopy
cp comparedPlots/cleanedStandAloneMuons_NValidHits.p* plotsToCopy

cp comparedPlots/generalTracks_pt.p*                  plotsToCopy
cp comparedPlots/generalTracks_eta.p*                 plotsToCopy
cp comparedPlots/generalTracks_phi.p*                 plotsToCopy
cp comparedPlots/generalTracks_dxy.p*                 plotsToCopy
cp comparedPlots/generalTracks_dz.p*                  plotsToCopy
cp comparedPlots/generalTracks_chi2.p*                plotsToCopy
cp comparedPlots/generalTracks_NValidHits.p*          plotsToCopy

cp comparedPlots/minDeltaR.p*                         plotsToCopy
cp compareDxy.p*                                      plotsToCopy
cp styled_EffVsDxyRatio.p*                            plotsToCopy
