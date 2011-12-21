divide(const TString & ptscale = "", const TString & writeMode = "recreate",
       const char * file="acceptance0_50_newAccCuts_basehistos_1.root", const char * outname="acceptance0_50_newAccCuts.root") {

  if( ptscale != "" && ptscale != "_ptup" && ptscale != "_ptdown" ) {
    std::cout << "Error: ptscale = "+ptscale+". Possible ptscale values are \"\" (no scale correction), \"_ptup\" and \"ptdown\"" << std::endl;
    exit(1);
  }
  std::cout << "pt scale set to " << ptscale << std::endl;


  TFile* f1 = new TFile(file);

  TH2F* genPtRap1 = (TH2F*) f1->Get("genUps1");
  TH2F* genPtRap2 = (TH2F*) f1->Get("genUps2");
  TH2F* genPtRap3 = (TH2F*) f1->Get("genUps3");
  TH2F* genPtRap4 = (TH2F*) f1->Get("genUps4");
  TH2F* genPtRap5 = (TH2F*) f1->Get("genUps5");
  TH1F* genPt1 = (TH1F*) f1->Get("genPt1");
  TH1F* genPt2 = (TH1F*) f1->Get("genPt2");
  TH1F* genPt3 = (TH1F*) f1->Get("genPt3");
  TH1F* genPt4 = (TH1F*) f1->Get("genPt4");
  TH1F* genPt5 = (TH1F*) f1->Get("genPt5");

  cout << "entries = " << genPtRap1->GetEntries() << endl;

  TH2F* recoPtRap1 = (TH2F*) f1->Get("recoUps1");
  TH2F* recoPtRap2 = (TH2F*) f1->Get("recoUps2");
  TH2F* recoPtRap3 = (TH2F*) f1->Get("recoUps3");
  TH2F* recoPtRap4 = (TH2F*) f1->Get("recoUps4");
  TH2F* recoPtRap5 = (TH2F*) f1->Get("recoUps5");
  TH1F* recoPt1 = (TH1F*) f1->Get("recoUpsPt1");
  TH1F* recoPt2 = (TH1F*) f1->Get("recoUpsPt2");
  TH1F* recoPt3 = (TH1F*) f1->Get("recoUpsPt3");
  TH1F* recoPt4 = (TH1F*) f1->Get("recoUpsPt4");
  TH1F* recoPt5 = (TH1F*) f1->Get("recoUpsPt5");

  cout<<"test"<<recoPt2->GetEntries()<<endl;
  TFile out(outname,writeMode);
  TH2F* acc1 = (TH2F*)genPtRap1->Clone("acc_pt35_eta16_pt25_eta24_trk_unpol"+ptscale);
  acc1->Sumw2();
  acc1->Divide(recoPtRap1,genPtRap1,1,1,"B");

  TH2F* acc2 = (TH2F*)genPtRap2->Clone("acc_pt35_eta16_pt25_eta24_trk_helT"+ptscale);
  acc2->Sumw2();
  acc2->Divide(recoPtRap2,genPtRap2,1,1,"B");

  TH2F* acc3 = (TH2F*)genPtRap3->Clone("acc_pt35_eta16_pt25_eta24_trk_helL"+ptscale);
  acc3->Sumw2();
  acc3->Divide(recoPtRap3,genPtRap3,1,1,"B");

  TH2F* acc4 = (TH2F*)genPtRap4->Clone("acc_pt35_eta16_pt25_eta24_trk_csT"+ptscale);
  acc4->Sumw2();
  acc4->Divide(recoPtRap4,genPtRap4,1,1,"B");

  TH2F* acc5 = (TH2F*)genPtRap5->Clone("acc_pt35_eta16_pt25_eta24_trk_csL"+ptscale);
  acc5->Sumw2();
  acc5->Divide(recoPtRap5,genPtRap5,1,1,"B");

  TH1F* accPt1 = (TH1F*)recoPt1->Clone("acc_pt35_eta16_pt25_eta24_trk_unpol_pt"+ptscale);
  accPt1->Sumw2();
  accPt1->Divide(recoPt1,genPt1,1,1,"B");

  TH1F* accPt2 = (TH1F*)recoPt2->Clone("acc_pt35_eta16_pt25_eta24_trk_helT_pt"+ptscale);
  accPt2->Sumw2();
  accPt2->Divide(recoPt2,genPt2,1,1,"B");

  TH1F* accPt3 = (TH1F*)recoPt3->Clone("acc_pt35_eta16_pt25_eta24_trk_helL_pt"+ptscale);
  accPt3->Sumw2();
  accPt3->Divide(recoPt3,genPt3,1,1,"B");

  TH1F* accPt4 = (TH1F*)recoPt4->Clone("acc_pt35_eta16_pt25_eta24_trk_csT_pt"+ptscale);
  accPt4->Sumw2();
  accPt4->Divide(recoPt4,genPt4,1,1,"B");

  TH1F* accPt5 = (TH1F*)recoPt5->Clone("acc_pt35_eta16_pt25_eta24_trk_csL_pt"+ptscale);
  accPt5->Sumw2();
  accPt5->Divide(recoPt5,genPt5,1,1,"B");

  acc1->Write();
  acc2->Write();
  acc3->Write();
  acc4->Write();
  acc5->Write();
  accPt1->Write();
  accPt2->Write();
  accPt3->Write();
  accPt4->Write();
  accPt5->Write();

  out.Close();
}
