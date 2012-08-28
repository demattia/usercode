#ifndef FILESANDWEIGHTS_H
#define FILESANDWEIGHTS_H

#include <map>
#include <TString.h>

std::map<TString, double> filesAndWeightsMap( bool electrons )
{
  // The weights are the cross sections in pb. The number of processed events is read from a histograms in the tree root file.
   
  // Files and weights map
  std::map<TString, double> fw;

  // Signal
  // HZZ_analysis_20120722
  // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/Signal"] = 10000.;
  
  // Vector bosons
  // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/ZZ_analysis_20120814"] = 8.3; // This is a guess at the moment sigma(7 TeV) * 1.2
  fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/ZZ4l_analysis_20120817"] = 9.4*3*0.033658*3*0.033658; // This is the SM expectation. ATLAS measured ~ 9.4 but consistent within the uncertanties
  fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/WW_analysis_20120814"] = 57.1;
  fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/WZ_analysis_20120816"] = 32.*(0.03658*3)*(0.108*3);

  // TTbar
  fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/TTJets_analysis_20120817"] = 225.197; // MCFM, LO:136.3;
  
  // DY+jets
  fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/DYJets10_analysis_20120814"] = 915.;
  fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/DYJets50_analysis_20120814"] = 3503.7; // LO:2950.0;
  
  if (!electrons) {
    // QCD mu
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu15_analysis_20120814"] = 7.022E8 * 0.0039;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu20_analysis_20120722"] = 2.87E8 * 0.0065;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu30_analysis_20120814"] = 6.609E7 * 0.0122;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu50_analysis_20120722"] = 8081400.0 * 0.0218;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu80_analysis_20120722"] = 1024000.0 * 0.0395;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu120_analysis_20120814"] = 157800 * 0.0473;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu170_analysis_20120814"] = 34020 * 0.0676;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu300_analysis_20120814"] = 1757 * 0.0864;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu470_analysis_20120814"] = 115.2 * 0.1024;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu600_analysis_20120814"] = 27.01 * 0.0996;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu800_analysis_20120814"] = 3.57 * 0.1033;
    // fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDmu1000_analysis_20120814"] = 0.774 * 0.1097;
  
    // Data
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/Data_Mu_Run2012A1_analysis_20120814"] = 1.;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/Data_Mu_Run2012B1_analysis_20120816"] = 1.;
  }
  else {
    // QCD em
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDem20_analysis_20120814"] = 2.886E8 * 0.0101;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDem30_analysis_20120814"] = 7.433E7 * 0.0621;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDem80_analysis_20120814"] = 1191000.0 * 0.1539;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDem170_analysis_20120814"] = 30990 * 0.148;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDem250_analysis_20120814"] = 4250.0 * 0.131;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/QCDem350_analysis_20120814"] = 810.0 * 0.11;
  
    // Data
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/Data_Photon_Run2012A1_analysis_20120814"] = 1.;
    fw["/uscmst1b_scratch/lpc1/3DayLifetime/demattia/NewAnalysis/CMSSW_5_2_6/src/workdirs/Data_Photon_Run2012B1_analysis_20120814"] = 1.;
  }
  return fw;
}

#endif
