#include "../Selection.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <iostream>

// Apply the analysis cuts or the preselection cuts and filter the tree.
// This macro is specific to the tree from the main analysis.
void cutTree_BsMuMu( const bool data, const TString & eb, const TString & index = "0", const TString & inFile="tmva-trees-3000.root" )
{
  TString outFile;
  TFile infile(inFile);
  if( data ) {
    outFile += "data_afterCuts_";
    gDirectory->Cd("sidebandChan"+eb+"Events"+index);
  }
  else {
    outFile += "MC_afterCuts_";
    gDirectory->Cd("signalChan"+eb+"Events"+index);
  }
  outFile += eb+"_"+index+".root";

  TTree* inTree = (TTree*)gROOT->FindObject("events");
  if(!inTree){
    cout<<"Could not access yield tree!"<<endl;
    return;
  }

  std::cout << "file" << std::endl;

  TFile outfile(outFile,"recreate");

  std::cout << "file2" << std::endl;

  bool endcaps = 1;
  if( eb == "0" ) endcaps = 0;
  TString cuts(Selection2(endcaps, false));

  TTree* outTree = (TTree*)inTree->CopyTree(cuts);

  //close everything
  outTree->Write();
  outfile.Close();
  infile.Close();
}
