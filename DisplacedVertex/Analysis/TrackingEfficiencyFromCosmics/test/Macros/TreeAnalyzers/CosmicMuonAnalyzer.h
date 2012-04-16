//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 30 14:32:13 2012 by ROOT version 5.27/06b
// from TTree T/MuonsTree
// found on file: cosmicMuons1Leg.root
//////////////////////////////////////////////////////////

#ifndef CosmicMuonAnalyzer_h
#define CosmicMuonAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

#include "../../../../RootTreeProducers/interface/Track.h"
#include "../../../../RootTreeProducers/interface/GenParticle.h"
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class Track+;
#pragma link C++ class std::vector<Track>+;
#pragma link C++ class GenParticle+;
#pragma link C++ class std::vector<GenParticle>+;
#endif

bool MC = true;

class CosmicMuonAnalyzer {
public :
  // TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          run;
   std::vector<Track>   *tracks;
   std::vector<Track>   *muons;
   std::vector<GenParticle> *genParticles;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_tracks;   //!
   TBranch        *b_muons;   //!
   TBranch        *b_genParticles;   //!

   CosmicMuonAnalyzer(TTree *tree=0);
   virtual ~CosmicMuonAnalyzer();
   // virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CosmicMuonAnalyzer_cxx
CosmicMuonAnalyzer::CosmicMuonAnalyzer(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("");
      if (!f) {
         f = new TFile("");
      }
      tree = (TTree*)gDirectory->Get("T");

   }
   Init(tree);
}

CosmicMuonAnalyzer::~CosmicMuonAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CosmicMuonAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CosmicMuonAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CosmicMuonAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tracks = 0;
   muons = 0;
   genParticles = 0;
   // Set branch addresses and branch pointers
   // if (!tree) return;
   // fChain = tree;

   // Read all the trees
   fChain = new TChain("T");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg.root");
   // fChain->Add("../../cosmicMuons1Leg.root");

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("tracks", &tracks, &b_tracks);
   fChain->SetBranchAddress("muons", &muons, &b_muons);
   if( MC ) {
     fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
   }
   Notify();
}

Bool_t CosmicMuonAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CosmicMuonAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
// Int_t CosmicMuonAnalyzer::Cut(Long64_t entry)
// {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
//   return 1;
// }
#endif // #ifdef CosmicMuonAnalyzer_cxx
