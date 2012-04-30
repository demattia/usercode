#define AnalyzeTreeClass_cxx
#include "AnalyzeTreeClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

bool AnalyzeTreeClass::acceptanceCuts()
{
  if( fabs(leptonEtaH) < 2. &&
      fabs(leptonEtaL) < 2. &&
      leptonPtL > 33. &&
      leptonPtH > 33. ) return true;
  return false;
}

bool AnalyzeTreeClass::trackSelectionCuts()
{
  if( leptonQualityH_pass &&
      leptonQualityL_pass &&

      leptonAbsD0significanceH > 2. && // for muons
      leptonAbsD0significanceL > 2. && // for muons

      // leptonAbsD0significanceH > 3. && // for electrons
      // leptonAbsD0significanceL > 3. && // for electrons

      trackerIsolationH_pass &&
      trackerIsolationL_pass
      // Did not see the number of hits in the bigTree, is it applied before?
      ) return true;
  return false;
}

bool AnalyzeTreeClass::dileptonSelectionCuts()
{
  if( oppositeCharge_pass &&
      vertexChi2 < 5 &&
      maxHitsBeforeVertex <= 1 &&
      deltaRBetweenLeptons > 0.2 && // muons only
      vetoBackToBack > -0.95 && // muons only

      dPhicorr < 0.2 && // muons
      decayLengthSignificance2D > 5 && // muons

      // dPhicorr < 0.8 && // electrons
      // decayLengthSignificance2D > 8 && // electrons

      numStandAloneMuons == 0 && // we are not using standAloneMuons      

      _mass > 15.
      ) return true;
  return false;
}

bool AnalyzeTreeClass::analysisCuts()
{
  if( acceptanceCuts() &&
      trackSelectionCuts() &&
      dileptonSelectionCuts()
      ) return true;
  return false;
}

void AnalyzeTreeClass::Loop()
{
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  TFile * outputFile = new TFile("output.root", "RECREATE");
  outputFile->cd();
  TH1F * histoDecayLength2D = new TH1F("decayLength2D", "L_{xy}", 100000, 0, 100);
  TH1F * histoMass = new TH1F("mass", "mass", 6000, 0, 600);
  TH1F * controlAccept = new TH1F("controlAccept", "controlAccept", 100000, -3.2, 3.2);
  TH1F * controlReject = new TH1F("controlReject", "controlReject", 100000, -3.2, 3.2);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( analysisCuts() ) {
      histoDecayLength2D->Fill(decayLength2D);
      histoMass->Fill(_mass);

      if( !passesAllCuts ) {
      	std::cout << "we have a problem" << std::endl;
      }
    }

    // // Control variables
    // if( dPhicorr_pass ) {
    //   controlAccept->Fill(dPhicorr);
    // }
    // else {
    //   controlReject->Fill(dPhicorr);
    // }


    // To be fast
    // if( ientry > 1000 ) break;
  }

  outputFile->Write();
  outputFile->Close();

}
