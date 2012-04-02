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

class CosmicMuonAnalyzer {
public :
  // TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          run;
   vector<Track>   *tracks;
   vector<Track>   *muons;
   vector<GenParticle> *genParticles;

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../cosmicMuons1Leg.root");
      if (!f) {
         f = new TFile("../../cosmicMuons1Leg.root");
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
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_1.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_10.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_100.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_101.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_102.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_103.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_104.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_105.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_106.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_107.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_108.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_109.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_11.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_110.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_111.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_112.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_113.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_114.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_115.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_116.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_117.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_118.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_119.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_12.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_120.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_121.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_122.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_123.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_124.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_125.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_126.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_127.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_128.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_129.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_13.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_130.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_131.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_132.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_133.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_134.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_135.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_136.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_137.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_138.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_139.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_14.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_140.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_141.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_142.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_143.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_144.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_147.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_148.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_149.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_15.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_150.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_151.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_152.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_153.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_154.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_155.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_156.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_157.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_158.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_159.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_16.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_160.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_161.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_162.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_163.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_164.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_166.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_167.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_168.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_169.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_17.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_170.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_171.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_172.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_173.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_174.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_175.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_176.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_177.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_178.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_179.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_18.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_180.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_181.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_182.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_183.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_184.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_185.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_186.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_187.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_188.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_189.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_19.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_190.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_191.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_192.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_194.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_195.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_196.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_197.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_198.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_199.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_2.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_20.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_200.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_201.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_202.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_205.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_206.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_207.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_208.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_209.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_21.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_210.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_211.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_212.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_213.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_214.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_216.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_217.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_218.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_219.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_22.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_220.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_221.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_222.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_223.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_224.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_226.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_227.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_228.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_229.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_23.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_230.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_231.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_232.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_233.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_234.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_235.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_236.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_237.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_238.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_239.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_24.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_240.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_241.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_242.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_243.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_244.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_245.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_246.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_247.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_248.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_249.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_25.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_250.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_251.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_252.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_253.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_254.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_255.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_256.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_257.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_258.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_259.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_26.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_260.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_261.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_262.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_263.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_264.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_265.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_266.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_267.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_268.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_269.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_27.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_270.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_271.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_272.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_273.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_274.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_275.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_276.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_277.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_278.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_279.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_28.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_280.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_281.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_282.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_283.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_284.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_285.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_286.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_287.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_288.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_289.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_29.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_290.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_291.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_292.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_293.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_294.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_295.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_296.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_297.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_298.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_299.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_3.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_30.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_300.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_301.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_302.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_303.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_304.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_305.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_306.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_307.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_308.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_309.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_31.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_310.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_311.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_312.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_313.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_314.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_315.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_316.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_317.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_318.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_319.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_32.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_320.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_321.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_322.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_323.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_325.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_326.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_327.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_328.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_329.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_33.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_330.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_331.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_333.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_334.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_335.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_336.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_339.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_34.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_340.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_341.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_342.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_343.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_344.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_347.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_348.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_349.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_35.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_350.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_351.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_352.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_353.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_354.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_355.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_356.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_357.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_358.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_359.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_36.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_360.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_361.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_362.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_363.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_364.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_365.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_366.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_367.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_368.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_369.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_37.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_370.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_371.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_373.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_374.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_375.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_376.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_377.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_378.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_379.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_38.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_380.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_381.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_382.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_383.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_384.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_385.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_387.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_388.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_389.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_39.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_390.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_391.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_392.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_393.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_394.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_395.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_396.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_397.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_399.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_4.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_40.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_400.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_401.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_402.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_403.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_404.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_405.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_406.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_407.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_408.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_409.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_41.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_410.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_411.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_412.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_413.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_414.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_415.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_416.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_417.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_418.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_419.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_42.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_420.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_421.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_422.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_424.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_425.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_426.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_427.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_428.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_429.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_43.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_430.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_431.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_432.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_433.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_434.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_435.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_436.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_437.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_438.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_439.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_44.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_440.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_441.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_442.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_443.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_444.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_445.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_446.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_447.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_448.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_449.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_45.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_450.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_451.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_452.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_454.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_455.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_456.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_457.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_458.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_459.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_46.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_460.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_461.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_462.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_463.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_465.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_466.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_468.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_47.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_470.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_472.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_473.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_474.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_475.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_476.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_477.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_478.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_479.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_48.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_480.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_481.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_482.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_483.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_484.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_485.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_486.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_487.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_488.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_489.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_49.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_490.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_491.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_492.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_493.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_494.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_495.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_496.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_497.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_498.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_499.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_5.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_50.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_500.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_501.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_502.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_503.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_504.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_505.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_506.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_507.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_508.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_509.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_51.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_510.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_511.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_512.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_513.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_514.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_516.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_517.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_518.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_519.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_52.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_521.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_522.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_523.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_524.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_525.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_528.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_529.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_53.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_531.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_532.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_533.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_534.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_535.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_536.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_537.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_538.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_539.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_54.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_540.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_541.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_545.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_546.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_547.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_548.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_549.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_55.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_551.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_552.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_553.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_554.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_555.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_556.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_557.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_558.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_559.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_56.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_560.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_561.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_562.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_563.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_564.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_565.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_566.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_567.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_568.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_569.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_57.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_570.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_571.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_572.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_573.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_574.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_575.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_576.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_577.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_578.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_579.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_580.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_581.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_582.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_583.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_584.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_585.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_586.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_587.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_588.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_589.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_59.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_590.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_591.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_592.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_593.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_594.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_595.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_596.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_597.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_598.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_599.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_6.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_60.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_600.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_601.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_602.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_603.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_604.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_605.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_606.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_607.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_608.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_609.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_61.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_610.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_611.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_613.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_614.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_615.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_616.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_617.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_618.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_619.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_62.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_620.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_621.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_622.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_623.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_624.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_625.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_626.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_627.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_628.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_629.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_63.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_630.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_631.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_633.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_634.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_635.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_636.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_637.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_638.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_639.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_64.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_640.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_641.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_642.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_643.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_644.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_645.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_646.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_647.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_648.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_649.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_65.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_650.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_651.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_652.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_653.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_654.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_655.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_656.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_657.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_658.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_659.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_66.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_660.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_661.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_662.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_663.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_664.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_665.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_666.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_667.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_668.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_669.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_67.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_670.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_671.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_672.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_673.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_674.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_675.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_676.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_677.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_678.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_679.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_68.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_680.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_681.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_682.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_683.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_684.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_685.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_686.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_687.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_688.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_69.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_690.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_691.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_692.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_693.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_694.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_695.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_697.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_698.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_699.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_7.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_70.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_701.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_702.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_703.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_705.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_706.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_707.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_708.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_71.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_710.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_712.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_713.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_714.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_715.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_716.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_717.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_718.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_719.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_720.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_721.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_722.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_723.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_724.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_725.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_726.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_727.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_728.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_729.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_730.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_731.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_732.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_733.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_734.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_735.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_736.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_737.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_738.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_739.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_74.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_740.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_741.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_742.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_743.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_744.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_745.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_746.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_747.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_748.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_749.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_75.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_750.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_751.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_752.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_753.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_754.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_755.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_756.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_757.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_758.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_759.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_76.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_760.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_761.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_762.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_763.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_764.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_765.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_766.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_767.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_768.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_769.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_77.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_78.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_79.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_8.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_80.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_81.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_82.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_83.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_84.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_85.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_87.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_88.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_89.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_9.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_90.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_91.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_92.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_93.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_94.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_95.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_96.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_97.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_98.root");
   fChain->Add("/home/demattia/RootTrees/Cosmics/MC/cosmicMuons1Leg_99.root");




   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("tracks", &tracks, &b_tracks);
   fChain->SetBranchAddress("muons", &muons, &b_muons);
   fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
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
