/**
 * This macro shows how to access FW objects in an uncompiled,
 * nameless FWLite macro. It uses directly the TChain.
 * It accesses both reco::Muons and MC information.
 */
{
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");

  gROOT->SetBatch(true);

  TH2F * ptVSetaGen = new TH2F("ptVSetaGen","genMuon pt vs eta", 100, -3.2, 3.2, 100, 0, 120);
  TProfile * ptVSetaGenProfile = new TProfile("ptVSetaGenProfile", "genMuon pt vs eta", 100, -3.2, 3.2, 0, 120);

  TH2F * ptVSeta = new TH2F("ptVSeta","muon pt vs eta", 100, -3.2, 3.2, 100, 0, 120);
  TProfile * ptVSetaProfile = new TProfile("ptVSetaProfile", "muon pt vs eta", 100, -3.2, 3.2, 0, 120);

  TChain events("Events");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_0.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_1.root");
  events->GetEntry();

  // connection of products and branches
  reco::MuonCollection muons;
  TBranch *muonsBranch;
  events->SetBranchAddress(events->GetAlias("muons"),&muons,&muonsBranch); 
  cout << "muonsBranch = " << muonsBranch <<endl;
  edm::HepMCProduct evtMC;
  TBranch *sourceBranch;
  events->SetBranchAddress(events->GetAlias("source"),&evtMC,&sourceBranch); 
  cout << "sourceBranch = " << sourceBranch <<endl;

  HepMC::GenEvent * Evt;
  reco::MuonCollection *recoMuon;
  for( unsigned int index = 0; index < events->GetEntries(); index++) { 
    sourceBranch->SetAddress(&evtMC);
    sourceBranch->GetEntry(index);
    muonsBranch->SetAddress(&muons);
    muonsBranch->GetEntry(index);
    events->GetEntry(index,0);
    Evt = evtMC.GetEvent();

    // DO NOT USE THE particle_const_iterator OR IT WILL
    // JUST EXIT NOT WORKING AND NOT GIVING ANY ERROR MESSAGE
    HepMC::GenEvent::particle_iterator part=Evt->particles_begin();
    for( ; part!=Evt->particles_end(); ++part ) {
      double pt = sqrt(pow((*part)->momentum().px(),2)+pow((*part)->momentum().px(),2));
      double pz = (*part)->momentum().pz();
      double theta = atan2(pt,pz);
      double eta = -log(tan(theta/2.));
      ptVSetaGen->Fill( eta, pt );
      ptVSetaGenProfile->Fill( eta, pt );
    }
    for(unsigned int muonId=0; muonId<muons.size(); ++muonId) {
      //       cout << "muons["<<muonId<<"] = " << muons[muonId].pt()
      //            << ", " << muons[muonId].eta() << endl;
      ptVSeta->Fill(muons[muonId].eta(), muons[muonId].pt());
      ptVSetaProfile->Fill(muons[muonId].eta(), muons[muonId].pt());
    }
  }

  TCanvas genMuonCanvas("ptVSetaGen","gen muon pt VS eta", 1000, 800);
  genMuonCanvas.cd();
  ptVSetaGen->Draw();
  genMuonCanvas->Print("ptVSetaGen.pdf");
  TCanvas genMuonCanvasProfile("ptVSetaGen","gen muon pt VS eta", 1000, 800);
  genMuonCanvasProfile.cd();
  ptVSetaGenProfile->Draw();
  genMuonCanvasProfile.Print("ptVSetaGenProfile.pdf");

  TCanvas recoMuonCanvas("ptVSeta","reco muon pt VS eta", 1000, 800);
  recoMuonCanvas.cd();
  ptVSeta->Draw();
  recoMuonCanvas->Print("ptVSetaReco.pdf");
  TCanvas recoMuonCanvasProfile("ptVSeta","reco muon pt VS eta", 1000, 800);
  recoMuonCanvasProfile.cd();
  ptVSetaProfile->Draw();
  recoMuonCanvasProfile.Print("ptVSetaRecoProfile.pdf");
}
