#define MuonAnalyzer_cxx
#include "MuonAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TMath.h"
#include "Math/VectorUtil.h" 

#include <cmath>

// #define M_PI 3.1415926535897932385

inline double deltaBarePhi(double phi1, double phi2)
{
  double dphi = phi2-phi1; 
  if ( dphi > M_PI ) {
    dphi -= 2.0*M_PI;
  } else if ( dphi <= -M_PI ) {
    dphi += 2.0*M_PI;
  }
  return dphi;
}

inline double deltaPhi(float phi1, float phi2)
{
  using ROOT::Math::VectorUtil::Phi_mpi_pi;
  return deltaBarePhi(Phi_mpi_pi(phi2),Phi_mpi_pi(phi1));
}

inline double deltaR(const float & eta1, const float & phi1, const float & eta2, const float & phi2)
{
  return std::sqrt(std::pow(deltaPhi(phi1, phi2),2) + std::pow(eta1 - eta2, 2));
}

void savePlot(const TString & name, TH1F * histo, const TString & xTitle, const bool logY = false)
{
   TCanvas * canvas = new TCanvas(name, name);
   canvas->Draw();
   histo->Draw();
   histo->GetXaxis()->SetTitle(xTitle);
   canvas->SetLogy(logY);
   canvas->SaveAs(name+".png");
   canvas->SaveAs(name+".root");
}
void savePlot(const TString & name, TH2F * histo, const TString & xTitle, const TString & yTitle)
{
   TCanvas * canvas = new TCanvas(name, name);
   canvas->Draw();
   histo->Draw();
   histo->GetXaxis()->SetTitle(xTitle);
   histo->GetYaxis()->SetTitle(yTitle);
   histo->SetMarkerStyle(1);
   canvas->SaveAs(name+".png");
}

void MuonAnalyzer::Loop()
{

  gStyle->SetOptStat(0);

//   In a ROOT session, you can do:
//      Root > .L MuonAnalyzer.C
//      Root > MuonAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  TH2F * deltaRvsVxy = new TH2F("deltaRvsVxy", "deltaRvsVxy", 100, 0, 1, 1000, 0, 100);
  TH2F * deltaRvsGenEta = new TH2F("deltaRvsGenEta", "deltaRvsGenEta", 100, 0, 1, 100, -2.4, 2.4);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for( int i=0; i<tracks_; ++i ) {
      double deltaRValue = deltaR(tracks_eta[i], tracks_phi[i], tracks_genEta[i], tracks_genPhi[i]);
      deltaRvsVxy->Fill( deltaRValue, sqrt(pow(tracks_genVx[i], 2) + pow(tracks_genVy[i], 2)) );
      deltaRvsGenEta->Fill( deltaRValue, tracks_genEta[i] );
    }
  }

  savePlot("deltaRvsDxy", deltaRvsVxy, "#Delta R(reco, gen)", "vertex radius (cm)");
  savePlot("deltaRvsGenEta", deltaRvsGenEta, "#Delta R(reco, gen)", "gen #eta");
}
