#define computeEfficiency_cxx
#include "computeEfficiency.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <locale>

const double muMass = 105.65837;

TLorentzVector computeEfficiency::fromPtEtaPhiToPxPyPz( const double & pt, const double & eta, const double & phi )
{
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  double tmp = 2*atan(exp(-eta));
  double pz = pt*cos(tmp)/sin(tmp);
  double E = sqrt(px*px+py*py+pz*pz+muMass*muMass);

  return TLorentzVector(px,py,pz,E);
}

void computeEfficiency::Loop()
{
  gROOT->SetBatch(true);
  // Apply the default trigger cuts
  defaultTriggerCuts_ = true;
  bool validChambersCut = true;
  double maxParallelCutDouble = 3.2;

  // One more bin because it will use both extremes
  Int_t ptBins = ptMax_ - ptMin_ + 1;
  Int_t parallelBins = (3.2 - startingParallelCut_)*10 + 1;

  TFile * outputFile = new TFile(outputFileName_, "RECREATE");
  outputFile->cd();
  entriesHisto_ = new TH1F("entriesHisto", "entries histogram", 1, 0, 1);
  effMap_ = new TH2F("effMap", "efficiency map", ptBins, ptMin_, ptMax_+1, parallelBins, startingParallelCut_, maxParallelCutDouble+0.1);

  // gROOT->Reset();
  gStyle->SetOptStat(0);
//   In a ROOT session, you can do:
//      Root > .L computeEfficiency.C
//      Root > computeEfficiency t
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

   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("NohMuL2NoVtx",1);
   fChain->SetBranchStatus("ohMuL2NoVtxPt",1);
   fChain->SetBranchStatus("ohMuL2NoVtxPhi",1);
   fChain->SetBranchStatus("ohMuL2NoVtxEta",1);
   fChain->SetBranchStatus("ohMuL2NoVtxChg",1);
   fChain->SetBranchStatus("ohMuL2NoVtxPtErr",1);
   fChain->SetBranchStatus("ohMuL2NoVtxDr",1);
   fChain->SetBranchStatus("ohMuL2NoVtxDz",1);
   fChain->SetBranchStatus("ohMuL2NoVtxL1idx",1);
   fChain->SetBranchStatus("ohMuL2NoVtxNhits",1);
   fChain->SetBranchStatus("ohMuL2NoVtxNchambers",1);

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      if( jentry%100 == 0 ) std::cout << "Analyzing entry number " << jentry << std::endl;
      // std::cout << "Number of L2 NoVtx muons = " << NohMuL2NoVtx << std::endl;

      parallelDiff_ = -99.;
      if( NohMuL2NoVtx > 1 &&
          ((ohMuL2NoVtxPt[0] > 23. && ohMuL2NoVtxPt[1] > 23. && ohMuL2NoVtxNhits[0] && ohMuL2NoVtxNhits[1]) || !defaultTriggerCuts_) &&
          ((ohMuL2NoVtxNchambers[0] > 1) || !validChambersCut)) {

        // Fill the bin 0 with 1 count
        entriesHisto_->Fill(0);

        TLorentzVector firstMuon = fromPtEtaPhiToPxPyPz(ohMuL2NoVtxPt[0], ohMuL2NoVtxEta[0], ohMuL2NoVtxPhi[0]);
        TLorentzVector secondMuon = fromPtEtaPhiToPxPyPz(ohMuL2NoVtxPt[1], ohMuL2NoVtxEta[1], ohMuL2NoVtxPhi[1]);
        double px1 = firstMuon.Px();
        double py1 = firstMuon.Py();
        double pz1 = firstMuon.Pz();
        double px2 = secondMuon.Px();
        double py2 = secondMuon.Py();
        double pz2 = secondMuon.Pz();
        parallelDiff_ = acos((px1*px2 + py1*py2 + pz1*pz2)/sqrt(px1*px1 + py1*py1 + pz1*pz1)/sqrt(px2*px2 + py2*py2 + pz2*pz2));

        double ptCut = 0.;
        for( Int_t iPt = 0; iPt < ptBins; ++iPt ) {
          ptCut = iPt + ptMin_;
          // std::cout << "ptCut = " << ptCut << std::endl;
          for( Int_t iParallelCut = 0; iParallelCut < parallelBins; ++iParallelCut ) {
            parallelDiffCut_ = iParallelCut/10. + startingParallelCut_;
            // std::cout << "parallelDiffCut = " << parallelDiffCut_ << std::endl;
            // Initialize for the selections
            if( ohMuL2NoVtxPt[0] > ptCut && ohMuL2NoVtxPt[1] > ptCut && parallelDiff_ < parallelDiffCut_ ) {
              // effMap_->Fill(ptCut, parallelDiffCut_);
              // Workaround since it is not filling the bins properly.
              effMap_->SetBinContent(iPt+1, iParallelCut+1, effMap_->GetBinContent(iPt+1, iParallelCut+1) + 1);
            }
          }
        }
      }
      // if (Cut(ientry) < 0) continue;
   }
   // effMap_->SaveAs("effMap.root");
   effMap_->Write();
   entriesHisto_->Write();
   outputFile->Write();
}
