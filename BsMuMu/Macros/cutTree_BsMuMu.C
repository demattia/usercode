#include "Selection.h"

void cutTree_BsMuMu( TString inFile="selection_Run2012A.root",
		     TString outFile="data_preselection.root",
		     const bool endcaps, const bool data,
		     const bool cut_based, const bool blinding,
		     const int splitting )
{
  TFile infile(inFile);
  gDirectory->Cd("detailedDimuonTree");
  TTree* inTree = (TTree*)gROOT->FindObject("probe_tree");
  if(!inTree){
    cout<<"Could not access yield tree!"<<endl;
    return;
  }

  TFile outfile(outFile,"recreate");

  TString selection(Selection(endcaps, data, cut_based, blinding, splitting));

  // std::cout << "SELECTION = " << selection << std::endl;

  outTree = (TTree*)inTree->CopyTree(selection);

  //close everything
  outTree->Write();
  outfile.Close();
  infile.Close();
}
