#include "TString.h"
#include <sstream>
#include <iostream>

TString Selection(const bool endcaps, const bool data, const bool cut_based = false, const bool blinding = false, const int splitting = 0, const TString & maxRun = "") {
  TString cuts = "";

  // TString trigger("((mu1_HLT_DoubleMu2BarrelBsL3 && mu2_HLT_DoubleMu2BarrelBsL3) || (mu1_HLT_DoubleMu2BsL3 && mu2_HLT_DoubleMu2BsL3) || (mu1_HLT_DoubleMu2Dimuon6BsL3 && mu2_HLT_DoubleMu2Dimuon6BsL3) || (mu1_HLT_DoubleMu3BsL3 && mu2_HLT_DoubleMu3BsL3) || (mu1_HLT_VertexmumuFilterBs345 && mu2_HLT_VertexmumuFilterBs345) || (mu1_HLT_VertexmumuFilterBs3p545 && mu2_HLT_VertexmumuFilterBs3p545) || (mu1_HLT_VertexmumuFilterBs4 && mu2_HLT_VertexmumuFilterBs4) || (mu1_HLT_VertexmumuFilterBs47 && mu2_HLT_VertexmumuFilterBs47) || (mu1_HLT_VertexmumuFilterBs6 && mu2_HLT_VertexmumuFilterBs6))");
  TString trigger("((mu1_HLT_VertexmumuFilterBs345 && mu2_HLT_VertexmumuFilterBs345) || (mu1_HLT_VertexmumuFilterBs3p545 && mu2_HLT_VertexmumuFilterBs3p545) || (mu1_HLT_VertexmumuFilterBs4 && mu2_HLT_VertexmumuFilterBs4) || (mu1_HLT_VertexmumuFilterBs47 && mu2_HLT_VertexmumuFilterBs47))");
  cuts += trigger;

  TString split("");
  std::ostringstream convert;
  convert << splitting;
  TString splittingString(convert.str());
  std::cout << "splittingString = " << splittingString << std::endl;
  if( splitting != -1 ) split = " && (event%3 == "+splittingString+")";
  cuts += split;
  if( maxRun != "" && data ) cuts += " && (run <= "+maxRun+")";

  // Muon-id: GlobalMuon prompt tight
  // TString muId = "(mu1_GMPT && mu2_GMPT)";
  // No muon-id, it will be applied later.
  // TString muId = "";
  TString muId = "mu1_GM && (mu1_globalChi2 < 10.) && (mu1_numberOfValidMuonHits > 0.) && (mu1_numMatchedStations > 1.) && (mu1_numberOfValidPixelHits > 0.) && (mu1_trackerLayersWithMeasurement > 5.) && mu2_GM && (mu2_globalChi2 < 10.) && (mu2_numberOfValidMuonHits > 0.) && (mu2_numMatchedStations > 1.) && (mu2_numberOfValidPixelHits > 0.) && (mu2_trackerLayersWithMeasurement > 5.)";
  if( muId != "" ) cuts += " && " + muId;

  // Preselection cuts
  TString preselection("(kinMass > 4.9 && kinMass < 5.9 && pt > 5. && pt < 9999. && mu1_pt > 4. && mu1_pt < 999. && mu2_pt > 4. && mu2_pt < 999. && l3d < 2. && l3dsig > 0. && l3dsig < 200. && NChi2 < 20. && delta3d < 0.1 && delta3d/delta3dErr < 5. && dca < 0.1 && acos(cosAlpha3D) < 1.0 && ntrk20 < 21 && minDca < 2.5 && isolation20 > 0. && (mu1_charge*mu2_charge == -1) && (lxysig > 3) && (abs(pvlip) < 1.0) && (abs(pvlip)/abs(pvlipErr) < 5.0) && kinMassError < 20)");
  cuts += " && " + preselection;
  //for new tree: mass -> kinMass, me<20 ie kinMassError < 20; isolation -> isolation20 ;  
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HeavyFlavorAnalysis/Bs2MuMu/macros2/preselection.cc?revision=1.13&view=markup
  //main analysis pre-selection: gmuid && (5.00<pt)&&(pt<9999.00) && (4.00<m1pt)&&(m1pt<999.00) &&(4.00<m2pt)&&(m2pt<999.00) && (flsxy>2.00) && (fl3d<2.00) && (pvip<0.10) && !(TMath::IsNaN(pvips)) && (pvips>0) && (pvips<5.00) && abs(pvlip) < 1.00 && abs(pvlips) < 5.00 &&  (closetrk<21) && (fls3d>0.00) && (fls3d<200.00) && (docatrk<2.50) && (maxdoca<0.10) && (chi2dof<20.00) && (iso>0.00) && (alpha<1.00) && (me<0.20)

  TString barrelCuts("(mu1_eta < 1.4 && mu1_eta > -1.4 && mu2_eta < 1.4 && mu2_eta > -1.4)");
  if( endcaps ) {
    TString endcapsCuts("!"+barrelCuts+" && mu1_eta < 2.4 && mu2_eta < 2.4 && mu1_eta > -2.4 && mu2_eta > -2.4");
    cuts += " && " + endcapsCuts;
  }
  else {
    cuts += " && " + barrelCuts;
  }

  // Cut-based analysis
  if( cut_based ) {
    TString lifetimeCuts("(delta3d < 0.008 && delta3d/delta3dErr < 2.000)");
    TString trackCuts("(isolation20 > 0.8 && dca > 0.015 && ntrk20 < 2 && pvw8 > 0.6)"); 
    cuts += " && " + lifetimeCuts + " && " + trackCuts;
    if( endcaps ) {
      TString endcapsPtCuts("((mu1_pt > mu2_pt && mu1_pt > 4.5 && mu2_pt > 4.2)||(mu1_pt < mu2_pt && mu1_pt > 4.2 && mu2_pt > 4.5))");
      TString candidateEndcapsCuts("(pt > 8.5 && cosAlpha3D > 0.99955 && l3dsig > 15.0 && NChi2 < 1.8)");
      cuts += " && " + endcapsPtCuts + " && " + candidateEndcapsCuts;
    } else {
      TString barrelPtCuts("((mu1_pt > mu2_pt && mu1_pt > 4.5 && mu2_pt > 4.0)||(mu1_pt < mu2_pt && mu1_pt > 4.0 && mu2_pt > 4.5))");
      TString candidateBarrelCuts("(pt > 6.5 && cosAlpha3D > 0.99875 && l3dsig > 13.0 && NChi2 < 2.2)");
      cuts += " && " + barrelPtCuts + " && " + candidateBarrelCuts;
    }
  }

  TString blindingCuts("&& ((kinMass > 4.9 && kinMass < 5.2) || (kinMass > 5.45 && kinMass < 5.9))");
  if( !blinding || !data ) blindingCuts = "";
  cuts += blindingCuts;

  return cuts;
}

// mvaMuonID selection cut
#define mvaMuonIDCut 0.3769
bool mvaMuonIDSelection(const double & mva1, const double & mva2 )
{
  return true;
  return( (mva1 > mvaMuonIDCut) && (mva2 > mvaMuonIDCut) );
}


// TString Selection2(const bool endcaps, const bool data, const bool cut_based = false, const bool blinding = false)
TString Selection2(const bool endcaps, const bool cut_based = false)
{
  TString cuts("");
  if( cut_based ) {
    TString lifetimeCuts("pvip < 0.008 && pvips < 2.000");
    TString trackCuts("pvw8 > 0.6 && m1q == -m2q && iso > 0.8 && docatrk > 0.015 && closetrk < 2");
    cuts += lifetimeCuts + " && " + trackCuts;
    if( endcaps ) {
      TString endcapsCuts("!(m1eta < 1.4 && m1eta > -1.4 && m2eta < 1.4 && m2eta > -1.4) && ((m1pt > m2pt && m1pt > 4.5 && m2pt > 4.2)||(m1pt < m2pt && m1pt > 4.2 && m2pt > 4.5))");
      TString candidateEndcapsCuts("pt > 8.5 && alpha < 0.030 && fls3d > 15.0 && NChi2 < 1.8");
      cuts += " && " + endcapsCuts + " && " + candidateEndcapsCuts;
    }
    else {
      TString barrelCuts("(m1eta < 1.4 && m1eta > -1.4 && m2eta < 1.4 && m2eta > -1.4) && ((m1pt > m2pt && m1pt > 4.5 && m2pt > 4.0)||(m1pt < m2pt && m1pt > 4.0 && m2pt > 4.5))");
      TString candidateBarrelCuts("pt > 6.5 && alpha < 0.050 && fls3d > 13.0 && NChi2 < 2.2");
      cuts += " && " + barrelCuts + " && " + candidateBarrelCuts;
    }
  }
  return( cuts );
}
