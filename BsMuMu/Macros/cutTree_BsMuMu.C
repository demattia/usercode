#include "Common/Selection.h"
#include "setdirs.h"

#define False false
#define True true

void cutTree_BsMuMu( TString inFile,
		     TString outFile,
		     const bool endcaps, const bool data,
		     const bool cut_based, const bool blinding,
		     const int splitting, const TString & maxRun ) {

  TFile infile(inFile);
  gDirectory->Cd("detailedDimuonTree");
  TTree* inTree = (TTree*)gROOT->FindObject("probe_tree");
  if(!inTree) {
    cout<<"Could not access yield tree!"<<endl;
    return;
  }

  TFile outfile(outFile,"recreate");
  
  TString selection(Selection(endcaps, data, cut_based, blinding, splitting, maxRun));
  //std::cout << "SELECTION = " << selection << std::endl;

  TTree* outTree = (TTree*)inTree->CopyTree(selection);

  //close everything
  outTree->Write();
  outfile.Close();
  infile.Close();
}
