#include "YieldTree.C"

void BookHistos();
void WriteHistos(Bool_t normalize, Double_t intLumi, Char_t *fNameOut);
//==========================================
void runYield(Char_t *fNameOut = "ppMuMu.root",//histos_upsilon1SMC_7TeV.root", //fileName where the histos will be stored
		  Bool_t matchMC = kTRUE,//doesn't matter now
		  Bool_t removeQQ = kTRUE,
		  Bool_t printGoodEvents = kFALSE
		  ){

  Bool_t normalize = kFALSE; //set to true *only* if you want to normalise the histos; 
                             //in standard mode: should be kFALSE
  Double_t intLumi = 16.675; //only relevant if "normalise"=kTRUE; 
                             //2.36 TeV: 13557/813=16.675; 900 GeV: 6724./404=16.644

  YieldTree tree;

  BookHistos();

  tree.Loop(removeQQ, matchMC, printGoodEvents);
  WriteHistos(normalize, intLumi, fNameOut);

}
//==========================================
void BookHistos(){
  t = new TTree("upsilonYield","upsilonYield");
  t->Branch("invariantMass",&invariantMass,"invariantMass/F");
  t->Branch("upsPt",&upsPt,"upsPt/F");
  t->Branch("upsRapidity",&upsRapidity,"upsRapidity/F");
  t->Branch("genUpsPt",&genUpsPt,"genUpsPt/F");
  t->Branch("genUpsRap",&genUpsRap,"genUpsRap/F");
  t->Branch("genMinMuPt",&genMinMuPt,"genMinMuPt/F");
  t->Branch("genMinMuEta",&genMinMuEta,"genMinMuEta/F");
  t->Branch("genPosMuPt",&genPosMuPt,"genPosMuPt/F");
  t->Branch("genPosMuEta",&genPosMuEta,"genPosMuEta/F");
  t->Branch("upsRapidity",&upsRapidity,"upsRapidity/F");
  t->Branch("muPlusPt",&muPlusPt,"muPlusPt/F");
  t->Branch("muPlusEta",&muPlusEta,"muPlusEta/F");
  t->Branch("muMinusPt",&muMinusPt,"muMinusPt/F");
  t->Branch("muMinusEta",&muMinusEta,"muMinusEta/F");
  t->Branch("McInvMass",&McInvMass,"McInvMass/F");
  t->Branch("GenRef",&GenRef,"GenRef/I");
  t->Branch("SampleFlag",&SampleFlag,"SampleFlag/I");

}

//==========================================
void WriteHistos(Bool_t normalize, Double_t intLumi, Char_t *fNameOut){
  TFile *fOut = new TFile(fNameOut, "RECREATE");
  t->Write();   
  fOut->Close();
}


