#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "tmvaMuonID/TMVAClassification_BDT.class.C"
#include "Common/Selection.h"

// This function adds the mvaMuonId variables to the trees. It also applies the cut.

void AddMuonID(const TString & file_name = "Barrel") {

  TFile * inputFile = new TFile(file_name);
  TTree *nt1 = (TTree*)inputFile->Get("probe_tree");
  Float_t mu1_validFrac, mu1_globalChi2, mu1_chi2LocPos, mu1_trkEHitsOut, mu1_segComp, mu1_glbTrackProb, mu1_chi2LocMom, mu1_trkVHits, mu1_pt, mu1_eta;
  Float_t mu2_validFrac, mu2_globalChi2, mu2_chi2LocPos, mu2_trkEHitsOut, mu2_segComp, mu2_glbTrackProb, mu2_chi2LocMom, mu2_trkVHits, mu2_pt, mu2_eta;
  Float_t mu1_MVAMuonID, mu2_MVAMuonID;
  Float_t cosa, alpha;
  nt1->SetBranchAddress("mu1_validFrac", &mu1_validFrac);
  nt1->SetBranchAddress("mu1_globalChi2", &mu1_globalChi2);
  nt1->SetBranchAddress("mu1_chi2LocPos", &mu1_chi2LocPos);
  nt1->SetBranchAddress("mu1_trkEHitsOut", &mu1_trkEHitsOut);
  nt1->SetBranchAddress("mu1_segComp", &mu1_segComp);
  nt1->SetBranchAddress("mu1_glbTrackProb", &mu1_glbTrackProb);
  nt1->SetBranchAddress("mu1_chi2LocMom", &mu1_chi2LocMom);
  nt1->SetBranchAddress("mu1_trkVHits", &mu1_trkVHits);
  nt1->SetBranchAddress("m1pt", &mu1_pt);
  nt1->SetBranchAddress("m1eta", &mu1_eta);
  nt1->SetBranchAddress("mu2_validFrac", &mu2_validFrac);
  nt1->SetBranchAddress("mu2_globalChi2", &mu2_globalChi2);
  nt1->SetBranchAddress("mu2_chi2LocPos", &mu2_chi2LocPos);
  nt1->SetBranchAddress("mu2_trkEHitsOut", &mu2_trkEHitsOut);
  nt1->SetBranchAddress("mu2_segComp", &mu2_segComp);
  nt1->SetBranchAddress("mu2_glbTrackProb", &mu2_glbTrackProb);
  nt1->SetBranchAddress("mu2_chi2LocMom", &mu2_chi2LocMom);
  nt1->SetBranchAddress("mu2_trkVHits", &mu2_trkVHits);
  nt1->SetBranchAddress("m2pt", &mu2_pt);
  nt1->SetBranchAddress("m2eta", &mu2_eta);
  nt1->SetBranchAddress("cosa", &cosa);

  TFile *g = new TFile(file_name+"_muonID.root", "RECREATE");
  TTree *newtree = nt1->CloneTree(0); // Do no copy the data yet

  // DO NOT delete the old tree and close the old file
  // add the branch to the new tree and try to fill it
  newtree->Branch("mu1_MVAMuonID", &mu1_MVAMuonID);
  newtree->Branch("mu2_MVAMuonID", &mu2_MVAMuonID);
  newtree->Branch("alpha", &alpha);

  std::vector<std::string> theInputVariables;
  theInputVariables.push_back("trkValidFract");
  theInputVariables.push_back("glbNChi2");
  theInputVariables.push_back("pt");
  theInputVariables.push_back("eta");
  theInputVariables.push_back("segComp");
  theInputVariables.push_back("chi2LocMom");
  theInputVariables.push_back("chi2LocPos");
  theInputVariables.push_back("glbTrackProb");
  theInputVariables.push_back("NTrkVHits");
  theInputVariables.push_back("NTrkEHitsOut");

  ReadBDT muonID(theInputVariables);

  std::vector<double> inputValues;
  inputValues.resize(10, 0.);


  int nentries=nt1->GetEntries();
  cout<<nentries<<endl;

  for( int i=0; i < nentries; i++){
    if( i%1000 == 0 ) std::cout << "processed: " << i << " events" << std::endl;
    nt1->GetEntry(i);

    inputValues[0] = mu1_validFrac;
    inputValues[1] = mu1_globalChi2;
    inputValues[2] = mu1_pt;
    inputValues[3] = mu1_eta;
    inputValues[4] = mu1_segComp;
    inputValues[5] = mu1_chi2LocMom;
    inputValues[6] = mu1_chi2LocPos;
    inputValues[7] = mu1_glbTrackProb;
    inputValues[8] = mu1_trkVHits;
    inputValues[9] = mu1_trkEHitsOut;
    mu1_MVAMuonID = muonID.GetMvaValue( inputValues );

    inputValues[0] = mu2_validFrac;
    inputValues[1] = mu2_globalChi2;
    inputValues[2] = mu2_pt;
    inputValues[3] = mu2_eta;
    inputValues[4] = mu2_segComp;
    inputValues[5] = mu2_chi2LocMom;
    inputValues[6] = mu2_chi2LocPos;
    inputValues[7] = mu2_glbTrackProb;
    inputValues[8] = mu2_trkVHits;
    inputValues[9] = mu2_trkEHitsOut;
    mu2_MVAMuonID = muonID.GetMvaValue( inputValues );

    alpha = acos(cosa);
    if( mvaMuonIDSelection(mu1_MVAMuonID, mu2_MVAMuonID) ) newtree->Fill();
  }

  g->Write();

  std::cout << "AddMuonID(" << file_name << ")... done." << std::endl;

}
