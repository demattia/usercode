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
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_2.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_3.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_4.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_5.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_6.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_7.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_8.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_9.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_10.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_11.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_12.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_13.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_14.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_15.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_16.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_17.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_18.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_19.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_20.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_21.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_22.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_23.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_24.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_25.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_26.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_27.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_28.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_29.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_30.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_31.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_32.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_33.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_34.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_35.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_36.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_37.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_38.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_39.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_40.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_41.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_42.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_43.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_44.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_45.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_46.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_47.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_48.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_49.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_50.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_51.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_52.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_53.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_54.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_55.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_56.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_57.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_58.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_59.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_60.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_61.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_62.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_63.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_64.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_65.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_66.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_67.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_68.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_69.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_70.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_71.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_72.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_73.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_74.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_75.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_76.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_77.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_78.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_79.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_80.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_81.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_82.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_83.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_84.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_85.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_86.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_87.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_88.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_89.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_90.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_91.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_92.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_93.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_94.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_95.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_96.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_97.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_98.root");
  events.Add("/data2/demattia/Data/MuonGun/MuonGun_99.root");

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
    cout << "Event number = " << index << endl;
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
