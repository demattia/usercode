//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 28 18:52:07 2012 by ROOT version 5.27/06b
// from TTree bigTree/bigTree
// found on file: Data_Mu_Run2011A4_analysis_20120219/histograms.root
//////////////////////////////////////////////////////////

#ifndef AnalyzeTreeClass_h
#define AnalyzeTreeClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AnalyzeTreeClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         cosThetaStar;
   Float_t         cosThetaStar_pass;
   Float_t         dPhi;
   Float_t         dPhi_pass;
   Float_t         dPhicorr;
   Float_t         dPhicorr_pass;
   Float_t         dPhitriggerCorr;
   Float_t         dPhitriggerCorr_pass;
   Float_t         decayLength2D;
   Float_t         decayLength2D_pass;
   Float_t         decayLengthSignificance2D;
   Float_t         decayLengthSignificance2D_pass;
   Float_t         deltaRBetweenLeptons;
   Float_t         deltaRBetweenLeptons_pass;
   Float_t         differentTrigObjects;
   Float_t         differentTrigObjects_pass;
   Float_t         eta_cand;
   Float_t         eta_cand_pass;
   Float_t         leptonAbsD0significanceH;
   Float_t         leptonAbsD0significanceH_pass;
   Float_t         leptonAbsD0significanceL;
   Float_t         leptonAbsD0significanceL_pass;
   Float_t         leptonD0H;
   Float_t         leptonD0H_pass;
   Float_t         leptonD0L;
   Float_t         leptonD0L_pass;
   Float_t         leptonD0significanceH;
   Float_t         leptonD0significanceH_pass;
   Float_t         leptonD0significanceL;
   Float_t         leptonD0significanceL_pass;
   Float_t         leptonEtaH;
   Float_t         leptonEtaH_pass;
   Float_t         leptonEtaL;
   Float_t         leptonEtaL_pass;
   Float_t         leptonPtH;
   Float_t         leptonPtH_pass;
   Float_t         leptonPtL;
   Float_t         leptonPtL_pass;
   Float_t         leptonQualityH;
   Float_t         leptonQualityH_pass;
   Float_t         leptonQualityL;
   Float_t         leptonQualityL_pass;
   Float_t         mass_calocorr;
   Float_t         mass_calocorr_pass;
   Float_t         mass_corr;
   Float_t         mass_corr_pass;
   Float_t         mass_scalecorr;
   Float_t         mass_scalecorr_pass;
   Float_t         mass_triggercorr;
   Float_t         mass_triggercorr_pass;
   Float_t         maxHitsBeforeVertex;
   Float_t         maxHitsBeforeVertex_pass;
   Float_t         maxHitsMissedAfterVertex;
   Float_t         maxHitsMissedAfterVertex_pass;
   Float_t         numPrimaryVertices;
   Float_t         numPrimaryVertices_pass;
   Float_t         numStandAloneMuons;
   Float_t         numStandAloneMuons_pass;
   Float_t         numTrigMatches;
   Float_t         numTrigMatches_pass;
   Float_t         oppositeCharge;
   Float_t         oppositeCharge_pass;
   Float_t         phi_cand;
   Float_t         phi_cand_pass;
   Float_t         trackerIsolationH;
   Float_t         trackerIsolationH_pass;
   Float_t         trackerIsolationL;
   Float_t         trackerIsolationL_pass;
   Float_t         validTracks;
   Float_t         validTracks_pass;
   Float_t         validVertex;
   Float_t         validVertex_pass;
   Float_t         vertexChi2;
   Float_t         vertexChi2_pass;
   Float_t         vetoBackToBack;
   Float_t         vetoBackToBack_pass;
   Float_t         _mass;
   Float_t         _pileup1BX;
   Float_t         _numDecays;
   Float_t         _ctau1;
   Float_t         _ctau2;
   Float_t         _leptonD01;
   Float_t         _leptonD02;
   Bool_t          passesAllCuts;
   Bool_t          passesAllCutsIgnoreLifetime;

   // List of branches
   TBranch        *b_cosThetaStar;   //!
   TBranch        *b_cosThetaStar_pass;   //!
   TBranch        *b_dPhi;   //!
   TBranch        *b_dPhi_pass;   //!
   TBranch        *b_dPhicorr;   //!
   TBranch        *b_dPhicorr_pass;   //!
   TBranch        *b_dPhitriggerCorr;   //!
   TBranch        *b_dPhitriggerCorr_pass;   //!
   TBranch        *b_decayLength2D;   //!
   TBranch        *b_decayLength2D_pass;   //!
   TBranch        *b_decayLengthSignificance2D;   //!
   TBranch        *b_decayLengthSignificance2D_pass;   //!
   TBranch        *b_deltaRBetweenLeptons;   //!
   TBranch        *b_deltaRBetweenLeptons_pass;   //!
   TBranch        *b_differentTrigObjects;   //!
   TBranch        *b_differentTrigObjects_pass;   //!
   TBranch        *b_eta_cand;   //!
   TBranch        *b_eta_cand_pass;   //!
   TBranch        *b_leptonAbsD0significanceH;   //!
   TBranch        *b_leptonAbsD0significanceH_pass;   //!
   TBranch        *b_leptonAbsD0significanceL;   //!
   TBranch        *b_leptonAbsD0significanceL_pass;   //!
   TBranch        *b_leptonD0H;   //!
   TBranch        *b_leptonD0H_pass;   //!
   TBranch        *b_leptonD0L;   //!
   TBranch        *b_leptonD0L_pass;   //!
   TBranch        *b_leptonD0significanceH;   //!
   TBranch        *b_leptonD0significanceH_pass;   //!
   TBranch        *b_leptonD0significanceL;   //!
   TBranch        *b_leptonD0significanceL_pass;   //!
   TBranch        *b_leptonEtaH;   //!
   TBranch        *b_leptonEtaH_pass;   //!
   TBranch        *b_leptonEtaL;   //!
   TBranch        *b_leptonEtaL_pass;   //!
   TBranch        *b_leptonPtH;   //!
   TBranch        *b_leptonPtH_pass;   //!
   TBranch        *b_leptonPtL;   //!
   TBranch        *b_leptonPtL_pass;   //!
   TBranch        *b_leptonQualityH;   //!
   TBranch        *b_leptonQualityH_pass;   //!
   TBranch        *b_leptonQualityL;   //!
   TBranch        *b_leptonQualityL_pass;   //!
   TBranch        *b_mass_calocorr;   //!
   TBranch        *b_mass_calocorr_pass;   //!
   TBranch        *b_mass_corr;   //!
   TBranch        *b_mass_corr_pass;   //!
   TBranch        *b_mass_scalecorr;   //!
   TBranch        *b_mass_scalecorr_pass;   //!
   TBranch        *b_mass_triggercorr;   //!
   TBranch        *b_mass_triggercorr_pass;   //!
   TBranch        *b_maxHitsBeforeVertex;   //!
   TBranch        *b_maxHitsBeforeVertex_pass;   //!
   TBranch        *b_maxHitsMissedAfterVertex;   //!
   TBranch        *b_maxHitsMissedAfterVertex_pass;   //!
   TBranch        *b_numPrimaryVertices;   //!
   TBranch        *b_numPrimaryVertices_pass;   //!
   TBranch        *b_numStandAloneMuons;   //!
   TBranch        *b_numStandAloneMuons_pass;   //!
   TBranch        *b_numTrigMatches;   //!
   TBranch        *b_numTrigMatches_pass;   //!
   TBranch        *b_oppositeCharge;   //!
   TBranch        *b_oppositeCharge_pass;   //!
   TBranch        *b_phi_cand;   //!
   TBranch        *b_phi_cand_pass;   //!
   TBranch        *b_trackerIsolationH;   //!
   TBranch        *b_trackerIsolationH_pass;   //!
   TBranch        *b_trackerIsolationL;   //!
   TBranch        *b_trackerIsolationL_pass;   //!
   TBranch        *b_validTracks;   //!
   TBranch        *b_validTracks_pass;   //!
   TBranch        *b_validVertex;   //!
   TBranch        *b_validVertex_pass;   //!
   TBranch        *b_vertexChi2;   //!
   TBranch        *b_vertexChi2_pass;   //!
   TBranch        *b_vetoBackToBack;   //!
   TBranch        *b_vetoBackToBack_pass;   //!
   TBranch        *b__mass;   //!
   TBranch        *b__pileup1BX;   //!
   TBranch        *b__numDecays;   //!
   TBranch        *b__ctau1;   //!
   TBranch        *b__ctau2;   //!
   TBranch        *b__leptonD01;   //!
   TBranch        *b__leptonD02;   //!
   TBranch        *b_passesAllCuts;   //!
   TBranch        *b_passesAllCutsIgnoreLifetime;   //!

   AnalyzeTreeClass(TTree *tree=0);
   virtual ~AnalyzeTreeClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual bool acceptanceCuts();
   virtual bool analysisCuts();
   virtual bool trackSelectionCuts();
   virtual bool dileptonSelectionCuts();
};

#endif

#ifdef AnalyzeTreeClass_cxx
AnalyzeTreeClass::AnalyzeTreeClass(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Data_Mu_Run2011A4_analysis_20120219/histograms.root");
    if (!f) {
      f = new TFile("Data_Mu_Run2011A4_analysis_20120219/histograms.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("muTrackAnalysis");
    TDirectory * innerDir = (TDirectory*)dir->Get("dileptons_background_HLT_L2DoubleMu30_NoVertex_v3");
    tree = (TTree*)innerDir->Get("bigTree");
    // tree = (TTree*)gDirectory->Get("bigTree");
  }
  Init(tree);
}

AnalyzeTreeClass::~AnalyzeTreeClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeTreeClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeTreeClass::LoadTree(Long64_t entry)
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

void AnalyzeTreeClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("cosThetaStar", &cosThetaStar, &b_cosThetaStar);
   fChain->SetBranchAddress("cosThetaStar_pass", &cosThetaStar_pass, &b_cosThetaStar_pass);
   fChain->SetBranchAddress("dPhi", &dPhi, &b_dPhi);
   fChain->SetBranchAddress("dPhi_pass", &dPhi_pass, &b_dPhi_pass);
   fChain->SetBranchAddress("dPhicorr", &dPhicorr, &b_dPhicorr);
   fChain->SetBranchAddress("dPhicorr_pass", &dPhicorr_pass, &b_dPhicorr_pass);
   fChain->SetBranchAddress("dPhitriggerCorr", &dPhitriggerCorr, &b_dPhitriggerCorr);
   fChain->SetBranchAddress("dPhitriggerCorr_pass", &dPhitriggerCorr_pass, &b_dPhitriggerCorr_pass);
   fChain->SetBranchAddress("decayLength2D", &decayLength2D, &b_decayLength2D);
   fChain->SetBranchAddress("decayLength2D_pass", &decayLength2D_pass, &b_decayLength2D_pass);
   fChain->SetBranchAddress("decayLengthSignificance2D", &decayLengthSignificance2D, &b_decayLengthSignificance2D);
   fChain->SetBranchAddress("decayLengthSignificance2D_pass", &decayLengthSignificance2D_pass, &b_decayLengthSignificance2D_pass);
   fChain->SetBranchAddress("deltaRBetweenLeptons", &deltaRBetweenLeptons, &b_deltaRBetweenLeptons);
   fChain->SetBranchAddress("deltaRBetweenLeptons_pass", &deltaRBetweenLeptons_pass, &b_deltaRBetweenLeptons_pass);
   fChain->SetBranchAddress("differentTrigObjects", &differentTrigObjects, &b_differentTrigObjects);
   fChain->SetBranchAddress("differentTrigObjects_pass", &differentTrigObjects_pass, &b_differentTrigObjects_pass);
   fChain->SetBranchAddress("eta_cand", &eta_cand, &b_eta_cand);
   fChain->SetBranchAddress("eta_cand_pass", &eta_cand_pass, &b_eta_cand_pass);
   fChain->SetBranchAddress("leptonAbsD0significanceH", &leptonAbsD0significanceH, &b_leptonAbsD0significanceH);
   fChain->SetBranchAddress("leptonAbsD0significanceH_pass", &leptonAbsD0significanceH_pass, &b_leptonAbsD0significanceH_pass);
   fChain->SetBranchAddress("leptonAbsD0significanceL", &leptonAbsD0significanceL, &b_leptonAbsD0significanceL);
   fChain->SetBranchAddress("leptonAbsD0significanceL_pass", &leptonAbsD0significanceL_pass, &b_leptonAbsD0significanceL_pass);
   fChain->SetBranchAddress("leptonD0H", &leptonD0H, &b_leptonD0H);
   fChain->SetBranchAddress("leptonD0H_pass", &leptonD0H_pass, &b_leptonD0H_pass);
   fChain->SetBranchAddress("leptonD0L", &leptonD0L, &b_leptonD0L);
   fChain->SetBranchAddress("leptonD0L_pass", &leptonD0L_pass, &b_leptonD0L_pass);
   fChain->SetBranchAddress("leptonD0significanceH", &leptonD0significanceH, &b_leptonD0significanceH);
   fChain->SetBranchAddress("leptonD0significanceH_pass", &leptonD0significanceH_pass, &b_leptonD0significanceH_pass);
   fChain->SetBranchAddress("leptonD0significanceL", &leptonD0significanceL, &b_leptonD0significanceL);
   fChain->SetBranchAddress("leptonD0significanceL_pass", &leptonD0significanceL_pass, &b_leptonD0significanceL_pass);
   fChain->SetBranchAddress("leptonEtaH", &leptonEtaH, &b_leptonEtaH);
   fChain->SetBranchAddress("leptonEtaH_pass", &leptonEtaH_pass, &b_leptonEtaH_pass);
   fChain->SetBranchAddress("leptonEtaL", &leptonEtaL, &b_leptonEtaL);
   fChain->SetBranchAddress("leptonEtaL_pass", &leptonEtaL_pass, &b_leptonEtaL_pass);
   fChain->SetBranchAddress("leptonPtH", &leptonPtH, &b_leptonPtH);
   fChain->SetBranchAddress("leptonPtH_pass", &leptonPtH_pass, &b_leptonPtH_pass);
   fChain->SetBranchAddress("leptonPtL", &leptonPtL, &b_leptonPtL);
   fChain->SetBranchAddress("leptonPtL_pass", &leptonPtL_pass, &b_leptonPtL_pass);
   fChain->SetBranchAddress("leptonQualityH", &leptonQualityH, &b_leptonQualityH);
   fChain->SetBranchAddress("leptonQualityH_pass", &leptonQualityH_pass, &b_leptonQualityH_pass);
   fChain->SetBranchAddress("leptonQualityL", &leptonQualityL, &b_leptonQualityL);
   fChain->SetBranchAddress("leptonQualityL_pass", &leptonQualityL_pass, &b_leptonQualityL_pass);
   fChain->SetBranchAddress("mass_calocorr", &mass_calocorr, &b_mass_calocorr);
   fChain->SetBranchAddress("mass_calocorr_pass", &mass_calocorr_pass, &b_mass_calocorr_pass);
   fChain->SetBranchAddress("mass_corr", &mass_corr, &b_mass_corr);
   fChain->SetBranchAddress("mass_corr_pass", &mass_corr_pass, &b_mass_corr_pass);
   fChain->SetBranchAddress("mass_scalecorr", &mass_scalecorr, &b_mass_scalecorr);
   fChain->SetBranchAddress("mass_scalecorr_pass", &mass_scalecorr_pass, &b_mass_scalecorr_pass);
   fChain->SetBranchAddress("mass_triggercorr", &mass_triggercorr, &b_mass_triggercorr);
   fChain->SetBranchAddress("mass_triggercorr_pass", &mass_triggercorr_pass, &b_mass_triggercorr_pass);
   fChain->SetBranchAddress("maxHitsBeforeVertex", &maxHitsBeforeVertex, &b_maxHitsBeforeVertex);
   fChain->SetBranchAddress("maxHitsBeforeVertex_pass", &maxHitsBeforeVertex_pass, &b_maxHitsBeforeVertex_pass);
   fChain->SetBranchAddress("maxHitsMissedAfterVertex", &maxHitsMissedAfterVertex, &b_maxHitsMissedAfterVertex);
   fChain->SetBranchAddress("maxHitsMissedAfterVertex_pass", &maxHitsMissedAfterVertex_pass, &b_maxHitsMissedAfterVertex_pass);
   fChain->SetBranchAddress("numPrimaryVertices", &numPrimaryVertices, &b_numPrimaryVertices);
   fChain->SetBranchAddress("numPrimaryVertices_pass", &numPrimaryVertices_pass, &b_numPrimaryVertices_pass);
   fChain->SetBranchAddress("numStandAloneMuons", &numStandAloneMuons, &b_numStandAloneMuons);
   fChain->SetBranchAddress("numStandAloneMuons_pass", &numStandAloneMuons_pass, &b_numStandAloneMuons_pass);
   fChain->SetBranchAddress("numTrigMatches", &numTrigMatches, &b_numTrigMatches);
   fChain->SetBranchAddress("numTrigMatches_pass", &numTrigMatches_pass, &b_numTrigMatches_pass);
   fChain->SetBranchAddress("oppositeCharge", &oppositeCharge, &b_oppositeCharge);
   fChain->SetBranchAddress("oppositeCharge_pass", &oppositeCharge_pass, &b_oppositeCharge_pass);
   fChain->SetBranchAddress("phi_cand", &phi_cand, &b_phi_cand);
   fChain->SetBranchAddress("phi_cand_pass", &phi_cand_pass, &b_phi_cand_pass);
   fChain->SetBranchAddress("trackerIsolationH", &trackerIsolationH, &b_trackerIsolationH);
   fChain->SetBranchAddress("trackerIsolationH_pass", &trackerIsolationH_pass, &b_trackerIsolationH_pass);
   fChain->SetBranchAddress("trackerIsolationL", &trackerIsolationL, &b_trackerIsolationL);
   fChain->SetBranchAddress("trackerIsolationL_pass", &trackerIsolationL_pass, &b_trackerIsolationL_pass);
   fChain->SetBranchAddress("validTracks", &validTracks, &b_validTracks);
   fChain->SetBranchAddress("validTracks_pass", &validTracks_pass, &b_validTracks_pass);
   fChain->SetBranchAddress("validVertex", &validVertex, &b_validVertex);
   fChain->SetBranchAddress("validVertex_pass", &validVertex_pass, &b_validVertex_pass);
   fChain->SetBranchAddress("vertexChi2", &vertexChi2, &b_vertexChi2);
   fChain->SetBranchAddress("vertexChi2_pass", &vertexChi2_pass, &b_vertexChi2_pass);
   fChain->SetBranchAddress("vetoBackToBack", &vetoBackToBack, &b_vetoBackToBack);
   fChain->SetBranchAddress("vetoBackToBack_pass", &vetoBackToBack_pass, &b_vetoBackToBack_pass);
   fChain->SetBranchAddress("_mass", &_mass, &b__mass);
   fChain->SetBranchAddress("_pileup1BX", &_pileup1BX, &b__pileup1BX);
   fChain->SetBranchAddress("_numDecays", &_numDecays, &b__numDecays);
   fChain->SetBranchAddress("_ctau1", &_ctau1, &b__ctau1);
   fChain->SetBranchAddress("_ctau2", &_ctau2, &b__ctau2);
   fChain->SetBranchAddress("_leptonD01", &_leptonD01, &b__leptonD01);
   fChain->SetBranchAddress("_leptonD02", &_leptonD02, &b__leptonD02);
   fChain->SetBranchAddress("passesAllCuts", &passesAllCuts, &b_passesAllCuts);
   fChain->SetBranchAddress("passesAllCutsIgnoreLifetime", &passesAllCutsIgnoreLifetime, &b_passesAllCutsIgnoreLifetime);
   Notify();
}

Bool_t AnalyzeTreeClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyzeTreeClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeTreeClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyzeTreeClass_cxx
