#ifndef HIVARIABLES_H
#define HIVARIABLES_H
/**
 * The HiVariables class evaluates some higher level variables and 
 * fills the corresponding histograms.
 * The variables are: SumEt; Mtot (total mass of the jets);
 * Centrality (SumEt/Mtot); Aplanarity; Sphericity; Thrust;
 *
 * The Jet struct is a container of the jets variables ot be passed
 * to HiVariables by the method Fill(), which requires: a vector of
 * Jet objects, the fCurrent sample number and the TOTAL number of
 * events. It uses the HistoSampler to rescale the histograms for
 * the different cross sections.
 *
 * Author M. De Mattia - 9/10/2007
 */

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <cstring>
// For the Matrix
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include "AnalysisExamples/AnalysisClasses/interface/SimpleJet.h"

#include <iostream>

using namespace std;

class HiVariables {

  SimpleJet * _tempJet;
  float _temp_pt;
  float _temp_pz;
  float _temp_e;
  TLorentzVector _momJ;
  TVector3 _bst;
  float _thrust[3];
  float _SumEt;
  float _Centrality;
  float _Aplanarity;
  float _Sphericity;
  float _Y;
  float _Thrust;
  float _HjetMass;
  float _LjetMass;
  float _MjetBroadening;
  float _mjetBroadening;
  char _histo_name[11];

 public:

  TH1F * SumEt;
  TH1F * Centrality;
  TH1F * Aplanarity;
  TH1F * Sphericity;
  TH2F * Y_S;
  TH1F * Thrust;
  TH1F * HeavyJetMass;
  TH1F * LightJetMass;
  TH1F * JetMassDifference;
  TH1F * MajorBroadening;
  TH1F * MinorBroadening;
  TH1F * BroadeningDifference;

  HiVariables(const char * Histo_Name);

  ~HiVariables() {
//     delete SumEt;
//     delete Centrality;
//     delete Aplanarity;
//     delete Sphericity;
//     delete Y_S;
//     delete Thrust;
//     delete HeavyJetMass;
//     delete LightJetMass;
//     delete JetMassDifference;
//     delete MajorBroadening;
//     delete MinorBroadening;
//     delete BroadeningDifference;
  }
  void Fill(vector<SimpleJet> &);
  TH1F * SumEtPlot() {
    return(SumEt);
  }
  TH1F * CentralityPlot() {
    return(Centrality);
  }
  TH1F * AplanarityPlot() {
    return(Aplanarity);
  }
  TH1F * SphericityPlot() {
    return(Sphericity);
  }
  TH2F * Y_SPlot() {
    return(Y_S);
  }
  TH1F * ThrustPlot() {
    return(Thrust);
  }
  TH1F * HeavyJetMassPlot() {
    return(HeavyJetMass);
  }
  TH1F * LightJetMassPlot() {
    return(LightJetMass);
  }
  TH1F * JetMassDifferencePlot() {
    return(JetMassDifference);
  }
  TH1F * MajorBroadeningPlot() {
    return(MajorBroadening);
  }
  TH1F * MinorBroadeningPlot() {
    return(MinorBroadening);
  }
  TH1F * BroadeningDifferencePlot() {
    return(BroadeningDifference);
  }
  // To return the values
  float Get_SumEt() {
    return(_SumEt);
  }
  float Get_Centrality() {
    return(_Centrality);
  }
  float Get_Aplanarity() {
    return(_Aplanarity);
  }
  float Get_Sphericity() {
    return(_Sphericity);
  }
  float Get_Y() {
    return(_Y);
  }
  float Get_Thrust() {
    return(_Thrust);
  }
  float Get_HeavyJetMass() {
    return(_HjetMass);
  }
  float Get_LightJetMass() {
    return(_LjetMass);
  }
  float Get_MajorBroadening() {
    return(_MjetBroadening);
  }
  float Get_MinorBroadening() {
    return(_mjetBroadening);
  }
  /// To write a canvas with all the plots
  void Plot();

  /** This functions uses _momJB, it must be
   * called after _momJB has been filled
   */
  void EvalThrust(float **, int);

  // Member functions
  float Psi(float * a, float * b) {
    float scalar = sqrt((a[0]+b[0])+(a[1]+b[1])+(a[2]+b[2]));
    return (scalar>0) ? 1. : -1.;
  }
};

#endif // HIVARIABLES_H
