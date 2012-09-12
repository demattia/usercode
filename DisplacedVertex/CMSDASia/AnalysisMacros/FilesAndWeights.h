#ifndef FILESANDWEIGHTS_H
#define FILESANDWEIGHTS_H

#include <map>
#include <TString.h>

std::map<TString, double> filesAndWeightsMap( bool electrons )
{
  // The weights are the cross sections in pb. The number of processed events is read from a histograms in the tree root file.
   
  // Files and weights map
  std::map<TString, double> fw;
  //  bool elecrton_ = false;
  // Signal
  //  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs//Signal_200_050F_analysis_20120831"] = 1.;

  // Vector bosons
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/ZZ_analysis_20120831"] = 8.3;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/WZ_analysis_20120831"] = 22;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/WW_analysis_20120831"] = 57.1;
  // fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/WJetsToLNu_analysis_20120831"] = 362570.2;

  // TTbar
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/TTJets_analysis_20120831"] = 225.2;

  // DY JETS
  // DYJets10 cross section too small...
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/DYJets10_analysis_20120831"] = 915;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/DYJets50_analysis_20120831"] = 3503.7;

//

if (!electrons)
{
  // QCD mu
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu15_analysis_20120831"] = 2738580.0;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu20_analysis_20120831"] = 1865500.0;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu30_analysis_20120831"] = 806298.0;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu50_analysis_20120831"] = 176187.6;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu80_analysis_20120831"] = 40448.0;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu120_analysis_20120831"] = 7463.94;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu170_analysis_20120831"] = 2299.752;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu300_analysis_20120831"] = 151.8048;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu470_analysis_20120831"] = 11.79648;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu600_analysis_20120831"] = 2.690196;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu800_analysis_20120831"] = 0.368781;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDmu1000_analysis_20120831"] = 0.0849078;

  // Data
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Mu_Run2012A1_analysis_20120831"] = 1;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Mu_Run2012B1_analysis_20120831"] = 1.;
  // fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Mu_Run2012C1_analysis_20120831"] = 1.;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Mu_Run2012C2_analysis_20120903"] = 1.;

 }
else
 {
  // QCD em
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDem20_analysis_20120831"] = 2914860.0;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDem30_analysis_20120831"] = 4615893.0;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDem80_analysis_20120831"] = 183294.9;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDem170_analysis_20120831"] = 4586.52;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDem250_analysis_20120831"] = 556.75;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/QCDem350_analysis_20120831"] = 89.1;

  // Data
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Photon_Run2012A1_analysis_20120831"] = 1.;
  fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Photon_Run2012B1_analysis_20120831"] = 1.;
  // fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Photon_Run2012C1_analysis_20120831"] = 1.;
  // fw["/afs/cern.ch/work/d/demattia/CMSDAS/CMSSW_5_3_3/src/workdirs/Data_Photon_Run2012C2_analysis_20120831"] = 1.;
 }
  return fw;
}

#endif
