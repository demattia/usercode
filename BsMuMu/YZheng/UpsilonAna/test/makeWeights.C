Float_t getLow(TH2F* h, Float_t pt, Float_t eta){
  int i = h->FindBin(pt,eta);
  return h->GetBinContent(i) - h->GetBinError(i);
}

Float_t get(TH2F* h, Float_t pt, Float_t eta){
  int i = h->FindBin(pt,eta);
  return h->GetBinContent(i);
}

Float_t getHigh(TH2F* h, Float_t pt, Float_t eta){
  int i = h->FindBin(pt,eta);
  return h->GetBinContent(i) + h->GetBinError(i);
}

void makeWeights(TString inFile="upsilonYield.root", TString accFile="acceptance.root", TString idFile="fit_result_MuFromTkUpsTkM.root", TString trigFile="fit_result_HltFromUpsTm.root", TString outFile="upsilonYieldWeighted.root"){
  // open file with the yield
  TFile infile(inFile);
  gDirectory->Cd("upsilonYield");
  TTree* inTree = (TTree*)gROOT->FindObject("upsilonYield");
  if(!inTree){
    cout<<"Could not access yield tree!"<<endl;
    return;
  }
  Float_t invariantMass, upsPt, upsRapidity, muPlusPt, muPlusEta, muMinusPt, muMinusEta;
  inTree->SetBranchAddress("invariantMass",&invariantMass);
  inTree->SetBranchAddress("upsPt"        ,&upsPt);
  inTree->SetBranchAddress("upsRapidity"  ,&upsRapidity);
  inTree->SetBranchAddress("muPlusPt"     ,&muPlusPt);
  inTree->SetBranchAddress("muPlusEta"    ,&muPlusEta);
  inTree->SetBranchAddress("muMinusPt"    ,&muMinusPt);
  inTree->SetBranchAddress("muMinusEta"   ,&muMinusEta);

  // open file with acceptance
  TFile accfile(accFile);
  TH2F* acc = (TH2F*)gROOT->FindObject("acceptance");
  if(!acc){
    cout<<"Could not access acceptance histogram!"<<endl;
    return;
  }
  // open file with id efficiency
  TFile idfile(idFile);
  TH2F* id = (TH2F*)gROOT->FindObject("fit_eff_pt_eta");
  if(!acc){
    cout<<"Could not access id efficiency histogram!"<<endl;
    return;
  }

  // open file with trigger efficiency
  TFile trigfile(trigFile);
  TH2F* trig = (TH2F*)gROOT->FindObject("fit_eff_pt_eta");
  if(!acc){
    cout<<"Could not access trigger efficiency histogram!"<<endl;
    return;
  }

  // create file for the yield and weights
  TFile outfile(outFile,"recreate");
  gDirectory->mkdir("upsilonYield")->cd();
  TTree *outTree = inTree->CloneTree(0);
  Float_t weight;
  outTree->Branch("weight", &weight, "weight/F");

  //loop through the candidates and calculate weights
  for (Int_t i=0;i<inTree->GetEntries(); i++) {
    if(i%1000==0) std::cout<<i<<std::endl;
    inTree->GetEntry(i);
    weight = get(acc,upsPt,upsRapidity) * get(id,muPlusPt,muPlusEta) * get(id,muMinusPt,muMinusEta) * 
             ( 1. - (1.-get(trig,muPlusPt,muPlusEta)) * (1.-get(trig,muMinusPt,muMinusEta)) );
    if(weight>0){
      weight = 1./weight;
      outTree->Fill();
    }else{
      std::cout<<"skipping event: "<<invariantMass<<" "<<upsPt<<" "<<upsRapidity<<" "<<muPlusPt<<" "<<muPlusEta<<" "<<muMinusPt<<" "<<muMinusEta<<" "<<std::endl;
    }
  }

  //close everything
  outTree->Write();
  outfile.Close();
  infile.Close();
  idfile.Close();
  trigfile.Close();
}
