#define checkOpenHLT_cxx
#include "checkOpenHLT.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <locale>

const double muMass = 105.65837;

TLorentzVector checkOpenHLT::fromPtEtaPhiToPxPyPz( const double & pt, const double & eta, const double & phi )
{
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  double tmp = 2*atan(exp(-eta));
  double pz = pt*cos(tmp)/sin(tmp);
  double E = sqrt(px*px+py*py+pz*pz+muMass*muMass);

  return TLorentzVector(px,py,pz,E);
}

void checkOpenHLT::fillNameArray(std::string * nameArray)
{
  nameArray[0] = "first";
  nameArray[1] = "second";
  nameArray[2] = "third";
  nameArray[3] = "fourth";
}

void checkOpenHLT::prepareHistograms(const TString & name, const int bins, const double & min, const double & max, const TString & axisTitle)
{
  std::locale loc;
  std::string nameArray[4];
  fillNameArray(nameArray);

  for( int i=0; i<4; ++i ) {
    std::string nameCopy(nameArray[i]);
    TString namePart(name + std::toupper(nameArray[i][0], loc) + nameCopy.erase(0,1));
    histoMap_.insert(std::make_pair(namePart, new TH1F(namePart, name + "of the " + nameArray[i] + " muon", bins, min, max)));
    histoMap_[namePart]->GetXaxis()->SetTitle(axisTitle);
    for( int j=i+1; j<4; ++j ) {
      std::stringstream ss;
      ss << i << "_" << j;
      TString correlationName(name+"Correlation_"+ss.str());
      histoMap_.insert(std::make_pair(correlationName, new TH2F(correlationName, name + "correlation of the " + nameArray[i] +
                                                                " and " + nameArray[j] + " muons",
                                                                bins, min, max,
                                                                bins, min, max)));
      ss.str("");
      ss << " (" << i << ")";
      histoMap_[correlationName]->GetXaxis()->SetTitle(axisTitle + ss.str());
      ss.str("");
      ss << " (" << j << ")";
      histoMap_[correlationName]->GetYaxis()->SetTitle(axisTitle + ss.str());
    }
  }
}

void checkOpenHLT::prepareAllHistograms(const TString & name, TFile * outputFile)
{
  outputFile->cd();
  outputFile->mkdir(name);
  outputFile->cd(name);
  prepareHistograms("pt"+name, 100, 0., 300., "p_{T} [GeV/c]");
  prepareHistograms("eta"+name, 100, -3., 3., "#eta");
  prepareHistograms("phi"+name, 100, -3.2, 3.2, "#phi");
  // prepareHistograms("chg"+name, 2, -1, 1);
  prepareHistograms("nhits"+name, 30, 0, 30, "number of valid hits");
  prepareHistograms("nchambers"+name, 10, 0, 10, "number of valid chambers");
  histoMap_.insert(std::make_pair("parallelism"+name, new TH1F("Parallelism"+name, "parallelism"+name, 100, 0., 3.2)));
}

void checkOpenHLT::applyCuts(const int arraySize, const bool selectOnChambers, const double & parallelDiff,
                             const bool selectOnParallelism, bool * selectionArray)
{
  // Initialize for the selections
  selectionArray[0] = false;
  selectionArray[1] = false;
  selectionArray[2] = false;
  selectionArray[3] = false;
  for( int i=0; i<arraySize; ++i ) {
    if( ohMuL2NoVtxNhits[i] &&
        ((ohMuL2NoVtxNchambers[i] > 1) || !selectOnChambers) &&
        ((parallelDiff < 2.) || !selectOnParallelism) ) {
      selectionArray[i] = true;
    }
  }
}

void checkOpenHLT::fillAllHistograms(const TString & name, const int arraySize, const bool * selectionArray)
{
  fillHistograms("pt"+name, ohMuL2NoVtxPt, arraySize, selectionArray);
  fillHistograms("eta"+name, ohMuL2NoVtxEta, arraySize, selectionArray);
  fillHistograms("phi"+name, ohMuL2NoVtxPhi, arraySize, selectionArray);
  // fillHistograms("chg"+name, ohMuL2NoVtxChg, arraySize, selectionArray);
  fillHistograms("nhits"+name, ohMuL2NoVtxNhits, arraySize, selectionArray);
  fillHistograms("nchambers"+name, ohMuL2NoVtxNchambers, arraySize, selectionArray);
  if( arraySize > 1 && selectionArray[0] && selectionArray[1] ) histoMap_["parallelism"+name]->Fill(parallelDiff_);
}

void checkOpenHLT::saveHistograms(const TString & name)
{
  std::locale loc;
  std::string nameArray[4];
  fillNameArray(nameArray);

  for( int i=0; i<4; ++i ) {
    std::string nameCopy(nameArray[i]);
    TString namePart(name + std::toupper(nameArray[i][0], loc) + nameCopy.erase(0,1));
    TCanvas * canvas = new TCanvas;
    canvas->cd();
    histoMap_[namePart]->Draw();
    canvas->Print(dir_+namePart+".pdf");
    canvas->Print(dir_+namePart+".gif");
    for( int j=i+1; j<4; ++j ) {
      std::stringstream ss;
      ss << i << "_" << j;
      TString correlationName(name+"Correlation_"+ss.str());
      canvas->cd();
      histoMap_[correlationName]->Draw();
      histoMap_[correlationName]->SetMarkerStyle(1);
      // canvas->Print(dir_+correlationName+".pdf");
      // canvas->Print(dir_+correlationName+".gif");
    }
  }
}

void checkOpenHLT::saveAllHistograms(const TString & name)
{
  saveHistograms("pt"+name);
  saveHistograms("eta"+name);
  saveHistograms("phi"+name);
  saveHistograms("nhits"+name);
  saveHistograms("nchambers"+name);
  saveHistogram((TH1F*)histoMap_["parallelism"+name]);
}

void checkOpenHLT::saveHistogram(TH1F * histo)
{
  TCanvas * canvas = new TCanvas;
  canvas->cd();
  histo->Draw();
  // canvas->Print(dir_+TString(histo->GetName())+".pdf");
  // canvas->Print(dir_+TString(histo->GetName())+".gif");
}

void checkOpenHLT::Loop()
{
  // Apply the default trigger cuts
  defaultTriggerCuts_ = true;

  // gROOT->Reset();
  gStyle->SetOptStat(0);
//   In a ROOT session, you can do:
//      Root > .L checkOpenHLT.C
//      Root > checkOpenHLT t
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

   // Setup all histograms
   TFile * outputFile = new TFile("CheckOpenHLT.root", "RECREATE");
   TH1F * numMuons = new TH1F("NumMuons", "Number of muons", 5, 0, 4);

   TString noCutsName("_NoCuts_");
   TString oneValidHitName("_OneValidHit_");
   TString oneValidChamberName("_OneValidChamber_");
   TString parallelismCutName("_ParallelismCut_");

   prepareAllHistograms(noCutsName, outputFile);
   prepareAllHistograms(oneValidHitName, outputFile);
   prepareAllHistograms(oneValidChamberName, outputFile);
   prepareAllHistograms(parallelismCutName, outputFile);

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
      if( NohMuL2NoVtx > 1 ) {
        TLorentzVector firstMuon = fromPtEtaPhiToPxPyPz(ohMuL2NoVtxPt[0], ohMuL2NoVtxEta[0], ohMuL2NoVtxPhi[0]);
        TLorentzVector secondMuon = fromPtEtaPhiToPxPyPz(ohMuL2NoVtxPt[1], ohMuL2NoVtxEta[1], ohMuL2NoVtxPhi[1]);
        double px1 = firstMuon.Px();
        double py1 = firstMuon.Py();
        double pz1 = firstMuon.Pz();
        double px2 = secondMuon.Px();
        double py2 = secondMuon.Py();
        double pz2 = secondMuon.Pz();
        parallelDiff_ = acos((px1*px2 + py1*py2 + pz1*pz2)/sqrt(px1*px1 + py1*py1 + pz1*pz1)/sqrt(px2*px2 + py2*py2 + pz2*pz2));
      }

      numMuons->Fill(NohMuL2NoVtx);

      if( NohMuL2NoVtx > 0 ) {

        // Skip if need to apply the default trigger cuts and they do not pass the pt cut
        if( defaultTriggerCuts_ && !(NohMuL2NoVtx > 1 && ohMuL2NoVtxPt[0] > 23 && ohMuL2NoVtxPt[1] > 23) ) continue;

        int arraySize = std::min(NohMuL2NoVtx, 4);

        bool selectOnChambers = false;
        bool selectOnParallelism = false;

        bool selectionArray[4];

        // No cuts
        selectionArray[0] = true;
        selectionArray[1] = true;
        selectionArray[2] = true;
        selectionArray[3] = true;
        fillAllHistograms(noCutsName, arraySize, selectionArray);

        // Fill histograms for the > 0 valid hit cut
        applyCuts(arraySize, selectOnChambers, parallelDiff_, selectOnParallelism, selectionArray);
        fillAllHistograms(oneValidHitName, arraySize, selectionArray);

        // One valid chamber cut
        selectOnChambers = true;
        applyCuts(arraySize, selectOnChambers, parallelDiff_, selectOnParallelism, selectionArray);
        fillAllHistograms(oneValidChamberName, arraySize, selectionArray);

        // Anti-parallel cut
        selectOnParallelism = true;
        applyCuts(arraySize, selectOnChambers, parallelDiff_, selectOnParallelism, selectionArray);
        fillAllHistograms(parallelismCutName, arraySize, selectionArray);
      }
      // if (Cut(ientry) < 0) continue;
   }
   saveAllHistograms(noCutsName);
   saveAllHistograms(oneValidHitName);
   saveAllHistograms(oneValidChamberName);
   saveAllHistograms(parallelismCutName);

   outputFile->Write();
}
