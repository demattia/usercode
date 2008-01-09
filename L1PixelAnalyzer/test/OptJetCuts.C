{
  TFile * G = new TFile("./root/TDAna_ttH_120_tk3.root");

  G->cd();

  TCanvas * g = new TCanvas ("g", "TTh, Mh=120 GeV", 500, 700 );
  g->Divide(2,4);
  g->cd(1);
  N04->SetLineColor(kBlue);
  N04->Draw("LEGO");
  g->cd(2);
  N04->Draw("CONT4");
  g->cd(3);
  Perf04->SetLineColor(kBlue);
  Perf04->Draw("LEGO");
  g->cd(4);
  Perf04->Draw("CONT4");
  g->cd(5);
  Drmedall->SetLineColor(kBlue);
  Drmedall->Draw("LEGO");
  g->cd(6);
  Drmedall->Draw("CONT4");
  g->cd(7);
  Nlo->SetLineColor(kBlue);
  Nlo->Draw("LEGO");
  g->cd(8);
  Nlo->Draw("CONT4");
  g->Print("OptJetCuts1.ps");

  TCanvas * g2 = new TCanvas ("g2", "TTh, Mh=120 GeV", 500, 500 );
  g2->Divide(2,2);
  g2->cd(1);
  Hrecfrac->SetLineColor(kBlue);
  Hrecfrac->Draw("LEGO");
  g2->cd(2);
  Hrecfrac->Draw("CONT4");
  g2->cd(3);
  Trecfrac->SetLineColor(kBlue);
  Trecfrac->Draw("LEGO");
  g2->cd(4);
  Trecfrac->Draw("CONT4");
  g2->Print("OptJetCuts2.ps");

  gStyle->SetOptFit(111111);
  TCanvas * g3 = new TCanvas ("g3", "TTh, Mh=120 GeV", 500, 500 );
  g3->Divide(1,3);
  g3->cd(1);
  MHbest->SetLineColor(kBlue);
  MHbest->Fit("gaus","","",80., 160.);
  g3->cd(2);
  MTbest->SetLineColor(kBlue);
  MTbest->Fit("gaus","","",150., 250.);
  g3->cd(3);
  MWbest->SetLineColor(kBlue);
  MWbest->Fit("gaus","","",50., 130.);
  g3->Print("OptJetCuts3.ps");

}
