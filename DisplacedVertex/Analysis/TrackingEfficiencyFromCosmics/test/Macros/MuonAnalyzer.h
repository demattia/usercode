//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 28 16:01:46 2011 by ROOT version 5.27/06b
// from TTree T/Muons
// found on file: globalMuons_reco.root
//////////////////////////////////////////////////////////

#ifndef MuonAnalyzer_h
#define MuonAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
   const Int_t kMaxtracks = 4;

class MuonAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          run;
   Int_t           tracks_;
   UInt_t          tracks_fUniqueID[kMaxtracks];   //[tracks_]
   UInt_t          tracks_fBits[kMaxtracks];   //[tracks_]
   Double_t        tracks_pt[kMaxtracks];   //[tracks_]
   Double_t        tracks_ptError[kMaxtracks];   //[tracks_]
   Double_t        tracks_eta[kMaxtracks];   //[tracks_]
   Double_t        tracks_etaError[kMaxtracks];   //[tracks_]
   Double_t        tracks_phi[kMaxtracks];   //[tracks_]
   Double_t        tracks_phiError[kMaxtracks];   //[tracks_]
   Int_t           tracks_charge[kMaxtracks];   //[tracks_]
   Double_t        tracks_dxy[kMaxtracks];   //[tracks_]
   Double_t        tracks_dxyError[kMaxtracks];   //[tracks_]
   Double_t        tracks_dz[kMaxtracks];   //[tracks_]
   Double_t        tracks_dzError[kMaxtracks];   //[tracks_]
   Double_t        tracks_vx[kMaxtracks];   //[tracks_]
   Double_t        tracks_vy[kMaxtracks];   //[tracks_]
   Double_t        tracks_vz[kMaxtracks];   //[tracks_]
   Double_t        tracks_chi2[kMaxtracks];   //[tracks_]
   Double_t        tracks_normalizedChi2[kMaxtracks];   //[tracks_]
   Double_t        tracks_referencePointRadius[kMaxtracks];   //[tracks_]
   Double_t        tracks_referencePointZ[kMaxtracks];   //[tracks_]
   Int_t           tracks_nHits[kMaxtracks];   //[tracks_]
   Int_t           tracks_nValidHits[kMaxtracks];   //[tracks_]
   Int_t           tracks_nValidPlusInvalidHits[kMaxtracks];   //[tracks_]
   Double_t        tracks_innermostHitRadius[kMaxtracks];   //[tracks_]
   Double_t        tracks_innermostHitZ[kMaxtracks];   //[tracks_]
   Double_t        tracks_genPt[kMaxtracks];   //[tracks_]
   Double_t        tracks_genEta[kMaxtracks];   //[tracks_]
   Double_t        tracks_genPhi[kMaxtracks];   //[tracks_]
   Int_t           tracks_genCharge[kMaxtracks];   //[tracks_]
   Double_t        tracks_genDxy[kMaxtracks];   //[tracks_]
   Double_t        tracks_genDxyError[kMaxtracks];   //[tracks_]
   Double_t        tracks_genDz[kMaxtracks];   //[tracks_]
   Double_t        tracks_genDzError[kMaxtracks];   //[tracks_]
   Double_t        tracks_genVx[kMaxtracks];   //[tracks_]
   Double_t        tracks_genVy[kMaxtracks];   //[tracks_]
   Double_t        tracks_genVz[kMaxtracks];   //[tracks_]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_tracks_;   //!
   TBranch        *b_tracks_fUniqueID;   //!
   TBranch        *b_tracks_fBits;   //!
   TBranch        *b_tracks_pt;   //!
   TBranch        *b_tracks_ptError;   //!
   TBranch        *b_tracks_eta;   //!
   TBranch        *b_tracks_etaError;   //!
   TBranch        *b_tracks_phi;   //!
   TBranch        *b_tracks_phiError;   //!
   TBranch        *b_tracks_charge;   //!
   TBranch        *b_tracks_dxy;   //!
   TBranch        *b_tracks_dxyError;   //!
   TBranch        *b_tracks_dz;   //!
   TBranch        *b_tracks_dzError;   //!
   TBranch        *b_tracks_vx;   //!
   TBranch        *b_tracks_vy;   //!
   TBranch        *b_tracks_vz;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_normalizedChi2;   //!
   TBranch        *b_tracks_referencePointRadius;   //!
   TBranch        *b_tracks_referencePointZ;   //!
   TBranch        *b_tracks_nHits;   //!
   TBranch        *b_tracks_nValidHits;   //!
   TBranch        *b_tracks_nValidPlusInvalidHits;   //!
   TBranch        *b_tracks_innermostHitRadius;   //!
   TBranch        *b_tracks_innermostHitZ;   //!
   TBranch        *b_tracks_genPt;   //!
   TBranch        *b_tracks_genEta;   //!
   TBranch        *b_tracks_genPhi;   //!
   TBranch        *b_tracks_genCharge;   //!
   TBranch        *b_tracks_genDxy;   //!
   TBranch        *b_tracks_genDxyError;   //!
   TBranch        *b_tracks_genDz;   //!
   TBranch        *b_tracks_genDzError;   //!
   TBranch        *b_tracks_genVx;   //!
   TBranch        *b_tracks_genVy;   //!
   TBranch        *b_tracks_genVz;   //!

   // MuonAnalyzer(TTree *tree=0);
   MuonAnalyzer(const TString & treeName);
   virtual ~MuonAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MuonAnalyzer_cxx
// MuonAnalyzer::MuonAnalyzer(TTree *tree)
MuonAnalyzer::MuonAnalyzer(const TString & treeName)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  // if (tree == 0) {
  TTree * tree;
  // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("globalMuons_reco.root");
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(treeName);
  if (!f) {
    // f = new TFile("globalMuons_reco.root");
    f = new TFile(treeName);
  }
  tree = (TTree*)gDirectory->Get("T");
  // }
  Init(tree);
}

MuonAnalyzer::~MuonAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonAnalyzer::LoadTree(Long64_t entry)
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

void MuonAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("tracks", &tracks_, &b_tracks_);
   fChain->SetBranchAddress("tracks.fUniqueID", tracks_fUniqueID, &b_tracks_fUniqueID);
   fChain->SetBranchAddress("tracks.fBits", tracks_fBits, &b_tracks_fBits);
   fChain->SetBranchAddress("tracks.pt", tracks_pt, &b_tracks_pt);
   fChain->SetBranchAddress("tracks.ptError", tracks_ptError, &b_tracks_ptError);
   fChain->SetBranchAddress("tracks.eta", tracks_eta, &b_tracks_eta);
   fChain->SetBranchAddress("tracks.etaError", tracks_etaError, &b_tracks_etaError);
   fChain->SetBranchAddress("tracks.phi", tracks_phi, &b_tracks_phi);
   fChain->SetBranchAddress("tracks.phiError", tracks_phiError, &b_tracks_phiError);
   fChain->SetBranchAddress("tracks.charge", tracks_charge, &b_tracks_charge);
   fChain->SetBranchAddress("tracks.dxy", tracks_dxy, &b_tracks_dxy);
   fChain->SetBranchAddress("tracks.dxyError", tracks_dxyError, &b_tracks_dxyError);
   fChain->SetBranchAddress("tracks.dz", tracks_dz, &b_tracks_dz);
   fChain->SetBranchAddress("tracks.dzError", tracks_dzError, &b_tracks_dzError);
   fChain->SetBranchAddress("tracks.vx", tracks_vx, &b_tracks_vx);
   fChain->SetBranchAddress("tracks.vy", tracks_vy, &b_tracks_vy);
   fChain->SetBranchAddress("tracks.vz", tracks_vz, &b_tracks_vz);
   fChain->SetBranchAddress("tracks.chi2", tracks_chi2, &b_tracks_chi2);
   fChain->SetBranchAddress("tracks.normalizedChi2", tracks_normalizedChi2, &b_tracks_normalizedChi2);
   fChain->SetBranchAddress("tracks.referencePointRadius", tracks_referencePointRadius, &b_tracks_referencePointRadius);
   fChain->SetBranchAddress("tracks.referencePointZ", tracks_referencePointZ, &b_tracks_referencePointZ);
   fChain->SetBranchAddress("tracks.nHits", tracks_nHits, &b_tracks_nHits);
   fChain->SetBranchAddress("tracks.nValidHits", tracks_nValidHits, &b_tracks_nValidHits);
   fChain->SetBranchAddress("tracks.nValidPlusInvalidHits", tracks_nValidPlusInvalidHits, &b_tracks_nValidPlusInvalidHits);
   fChain->SetBranchAddress("tracks.innermostHitRadius", tracks_innermostHitRadius, &b_tracks_innermostHitRadius);
   fChain->SetBranchAddress("tracks.innermostHitZ", tracks_innermostHitZ, &b_tracks_innermostHitZ);
   fChain->SetBranchAddress("tracks.genPt", tracks_genPt, &b_tracks_genPt);
   fChain->SetBranchAddress("tracks.genEta", tracks_genEta, &b_tracks_genEta);
   fChain->SetBranchAddress("tracks.genPhi", tracks_genPhi, &b_tracks_genPhi);
   fChain->SetBranchAddress("tracks.genCharge", tracks_genCharge, &b_tracks_genCharge);
   fChain->SetBranchAddress("tracks.genDxy", tracks_genDxy, &b_tracks_genDxy);
   fChain->SetBranchAddress("tracks.genDxyError", tracks_genDxyError, &b_tracks_genDxyError);
   fChain->SetBranchAddress("tracks.genDz", tracks_genDz, &b_tracks_genDz);
   fChain->SetBranchAddress("tracks.genDzError", tracks_genDzError, &b_tracks_genDzError);
   fChain->SetBranchAddress("tracks.genVx", tracks_genVx, &b_tracks_genVx);
   fChain->SetBranchAddress("tracks.genVy", tracks_genVy, &b_tracks_genVy);
   fChain->SetBranchAddress("tracks.genVz", tracks_genVz, &b_tracks_genVz);
   Notify();
}

Bool_t MuonAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuonAnalyzer_cxx
