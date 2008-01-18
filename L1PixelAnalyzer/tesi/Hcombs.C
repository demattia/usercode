{
  g=new TFile("./root/TDAna_TTH_120.root");
  g->cd();
  
  TH1D * H1 =  dynamic_cast<TH1D*>(g->Get("Hnotpt"));
  TH1D * H2 =  dynamic_cast<TH1D*>(g->Get("Hpt"));
  TH1D * H3 =  dynamic_cast<TH1D*>(g->Get("Hnoteta"));
  TH1D * H4 =  dynamic_cast<TH1D*>(g->Get("Heta"));
  TH1D * H5 =  dynamic_cast<TH1D*>(g->Get("Hnotdr"));
  TH1D * H6 =  dynamic_cast<TH1D*>(g->Get("Hdr"));
  TH1D * H7 =  dynamic_cast<TH1D*>(g->Get("MHnot"));
  TH1D * H8 =  dynamic_cast<TH1D*>(g->Get("MHbest"));
  H1->Sumw2();
  H2->Sumw2();
  H3->Sumw2();
  H4->Sumw2();
  H5->Sumw2();
  H6->Sumw2();
  H7->Sumw2();
  H8->Sumw2();
  double prob = H1->KolmogorovTest(H2);
  cout << "KS test on H Pt: P=" << prob << endl;
  prob = H3->KolmogorovTest(H4);
  cout <<  "KS test on H Eta: P=" << prob << endl;
  prob = H5->KolmogorovTest(H6);
  cout <<  "KS test on H dr: P=" << prob << endl;
  prob = H7->KolmogorovTest(H8);
  cout <<  "KS test on H mass: P=" << prob << endl;

  TCanvas * b = new TCanvas ( "b", "Higgs reconstruction kinematics", 500, 500 );
  b->Divide(2,2);
  
  b->cd(1);
  H2->DrawNormalized("HISTO");
  H1->SetLineColor(kRed);
  H1->SetMarkerColor(kRed);
  H1->SetMarkerSize(0.3);
  H1->DrawNormalized("PESAME");
  b->cd(2);
  H4->DrawNormalized("HISTO");
  H3->SetLineColor(kRed);
  H3->SetMarkerColor(kRed);
  H3->SetMarkerSize(0.3);
  H3->DrawNormalized("PESAME");
  b->cd(3);
  H6->DrawNormalized("HISTO");
  H5->SetLineColor(kRed);
  H5->SetMarkerColor(kRed);
  H5->SetMarkerSize(0.3);
  H5->DrawNormalized("PESAME");
  b->cd(4);
  H8->DrawNormalized("HISTO");
  H7->SetLineColor(kRed);
  H7->SetMarkerColor(kRed);
  H7->SetMarkerSize(0.3);
  H7->DrawNormalized("PESAME");
  b->Print("./ps/Hcombs.ps");
}
