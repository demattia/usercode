//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 19 13:21:45 2012 by ROOT version 5.32/00
// from TTree outputTree/outputTree
// found on file: histograms.root
//////////////////////////////////////////////////////////

#ifndef treeAnalyzer_h
#define treeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TreeProducer/TreeProducer/interface/Candidates.h"
#include "TreeProducer/TreeProducer/interface/TreeCandidate.h"
#include "TreeProducer/TreeProducer/interface/TreeLepton.h"
#include "TreeProducer/TreeProducer/src/LinkDef.h"
#include <TH1F.h>
#include <vector>
#include <string>
// #ifndef __CINT__
// #endif

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class treeAnalyzer {
public :
  int totalProcessedEvents_;
  int eventsPassingTrigger_;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Candidates      *candidates;
   std::vector<std::string> *triggers;

   // List of branches
   TBranch        *b_candidates;   //!
   TBranch        *b_triggers;   //!

   // treeAnalyzer(TTree *tree=0);
   treeAnalyzer(TString fileName = "", const double & weight = 1., const bool electrons = false);
   virtual ~treeAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   double weight_;
   bool electrons_;
   TString dirName_;
};

#endif

#ifdef treeAnalyzer_cxx
// treeAnalyzer::treeAnalyzer(TTree *tree) : fChain(0) 
// {
// // if parameter tree is not specified (or zero), connect the file
// // used to generate this class and read the Tree.
//    if (tree == 0) {
//       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histograms.root");
//       if (!f || !f->IsOpen()) {
//          f = new TFile("histograms.root");
//       }
//       TDirectory * dir = (TDirectory*)f->Get("histograms.root:/muTrackAnalysis");
//       dir->GetObject("outputTree",tree);
//    }
//    Init(tree);
// }

treeAnalyzer::treeAnalyzer(TString fileName, const double & weight, const bool electrons) : fChain(0), weight_(weight), electrons_(electrons), dirName_("")
{
  TString dirName("muTrackAnalysis");
  if( electrons_ ) dirName = "eTrackAnalysis";

  TChain * f = new TChain(dirName+"/outputTree");


  if( fileName == "" ) {
    // fileName = "/uscms_data/d3/demattia/CMSDAS/TreeProduction/CMSSW_4_2_7/src/workdirs/Signal_200_050F_analysis_20120616/histograms.root";
    fileName = "/uscms_data/d3/demattia/CMSDAS/TreeProduction/CMSSW_4_2_7/src/workdirs/Zmumu_analysis_20120616/histograms.root";
  }
  TFile tempFile(fileName, "READ");
  TDirectory * dir = (TDirectory*)tempFile.Get(dirName);
   totalProcessedEvents_ = ((TH1F*)dir->Get("totalProcessedEvents"))->GetBinContent(1);
//  eventsPassingTrigger_ = ((TH1F*)dir->Get("eventsPassingTrigger"))->GetBinContent(1);

  if( weight_ != 1. ) weight_ /= totalProcessedEvents_;

  // Get the name of the directory where the sample is. The output file will be written in the same directory.
  if( std::string(fileName).find_last_of("/") != std::string::npos ) {
    // dirName_ = std::string(fileName).substr(0, std::string(fileName).find_last_of("/"));
    std::string tempDir(std::string(fileName).substr(0, std::string(fileName).find_last_of("/")));
    dirName_ = tempDir.substr(std::string(tempDir).find_last_of("/")+1);
  }

  f->Add(fileName);
  Init(f);
}

treeAnalyzer::~treeAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   candidates = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("candidates", &candidates, &b_candidates);
   fChain->SetBranchAddress("triggers", &triggers, &b_triggers);
   Notify();
}

Bool_t treeAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeAnalyzer_cxx
