//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun  1 10:44:55 2012 by ROOT version 5.32/00
// from TTree treeAnalyzer/treeAnalyzer
// found on file: histograms.root
//////////////////////////////////////////////////////////

#ifndef treeAnalyzer_h
#define treeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>

// Candidates
#include "Candidates.h"
#include "TreeCandidate.h"
#include "TreeLepton.h"
#include "LinkDef.h"

#include "commonTools.C"

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <map>

// For pileup weights
#include "PileupReweighting.C"


// Fixed size dimensions of array or collections stored in the TTree if any.

class treeAnalyzer {
  public :
    int totalProcessedEvents_;
    int eventsPassingTrigger_;

   void     Loop();

   // Analysis cuts
   bool acceptanceCuts( std::vector<TreeCandidate>::const_iterator cand );
   bool trackSelectionCuts( std::vector<TreeCandidate>::const_iterator cand, const bool removeIsolationCut = false );
   bool lifetimeRelatedCuts( std::vector<TreeCandidate>::const_iterator candconst, bool removeLifetimeRelatedCuts = false, const bool decayLengthSigniInverted = false);
   bool dileptonSelectionCuts( std::vector<TreeCandidate>::const_iterator cand, const std::vector<TreeLepton> & leptons );
   bool analysisCuts(const std::vector<TreeCandidate>::const_iterator & cand, const std::vector<TreeLepton> & leptons, const bool removeIsolationCut, const bool removeLifetimeRelatedCuts, const bool decayLengthSigniInverted );

   bool triggerMatching( std::vector<TreeCandidate>::const_iterator cand, const std::vector<TreeLepton> & leptons );
   void initializeCuts();
   double ptCut_;
   double d0SignificanceCut_;
   double dPhiCorrCut_;
   double decayLengthSignificance2DCut_;

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
    // virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    Cut();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    double weight_;

    PileupReweighting puweights_;

    // Electron channel if true, muon if false
    bool electrons_;

    TString dirName_;
};


// treeAnalyzer::treeAnalyzer(TTree *tree) : fChain(0)
treeAnalyzer::treeAnalyzer(TString fileName, const double & weight, const bool electrons) : fChain(0), weight_(weight), electrons_(electrons), dirName_("")
{
  TString dirName("muTrackAnalysis");
  if( electrons_ ) dirName = "eTrackAnalysis";

  TChain * f = new TChain(dirName+"/outputTree");

  if( fileName == "" ) {
    // fileName = "/uscms_data/d3/demattia/CMSDAS/TreeProduction/CMSSW_4_2_7/src/workdirs/Signal_200_050F_analysis_20120616/histograms.root";
    fileName = "/afs/cern.ch/work/e/ejclemen/DisplacedLeptons2012/CMSSW_5_2_5/src/TreeProducer/TreeProducer/test/histograms.root";
  }
  TFile tempFile(fileName, "READ");
  TDirectory * dir = (TDirectory*)tempFile.Get(dirName);
  totalProcessedEvents_ = ((TH1F*)dir->Get("totalProcessedEvents"))->GetBinContent(1);
  //  eventsPassingTrigger_ = ((TH1F*)dir->Get("eventsPassingTrigger"))->GetBinContent(1);

  std::cout << "File : " << fileName.Data() << std::endl;

  if( weight_ != 1 ) weight_ /= totalProcessedEvents_;
  std::cout << "Weight : " << weight_ << std::endl;

  // Get the name of the directory where the sample is. The output file will be written in the same directory.
  if( std::string(fileName).find_last_of("/") != std::string::npos ) {
    // dirName_ = std::string(fileName).substr(0, std::string(fileName).find_last_of("/"));
    std::string tempDir(std::string(fileName).substr(0, std::string(fileName).find_last_of("/")));
    dirName_ = tempDir.substr(std::string(tempDir).find_last_of("/")+1);
  }

  f->Add(fileName);
  Init(f);

  // Initialise PU weights
  TString puFile = "PileupFiles/pileup_muon.root";
  if ( electrons ) puFile = "PileupFiles/pileup_electron.root";
  puweights_ = PileupReweighting(puFile);
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
  triggers = 0;
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
// Int_t treeAnalyzer::Cut(Long64_t entry)
Int_t treeAnalyzer::Cut()
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef treeAnalyzer_h
