{
  TFile * G = new TFile("./root/TDAna_ttH_120_tk3.root");

  TH1F * R = new TH1F ( "R", "Fraction of events with 4 HEM tags in first N jets", 
			20, -0.5, 19.5 );

  G->cd();
  E4NJSSS->Sumw2();
  N4NJSSS->Sumw2();
  for ( int i=1; i<21; i++ ) {
    double t = E4NJSSS->GetBinContent(i);
    double a = N4NJSSS->GetBinContent(i);
    double f = 0;
    if ( a>0 ) f=t/a;
    if ( a>0 ) sf = sqrt(f*(1-f)/a);
    R->SetBinContent(i,f);
    R->SetBinError(i,sf);
  }

  TCanvas * g = new TCanvas ("g", "TTh Mh=120 GeV", 500, 500 );
  g->Divide(1,2);
  g->cd(1);
  R->SetLineColor(kBlue);
  R->SetMarkerColor(kBlue);
  R->SetMarkerStyle(24);
  E4NJSSS->SetLineColor(kRed);
  E4NJSSS->SetLineWidth(2);
  N4NJSSS->SetLineWidth(2);
  N4NJSSS->Draw("HISTO");
  E4NJSSS->Draw("SAMEHISTO");
  g->cd(2);
  R->Draw("PE");
  g->Print("Njcut_4t.ps");

}
