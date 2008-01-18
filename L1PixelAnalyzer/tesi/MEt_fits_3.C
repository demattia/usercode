{

  gStyle->SetOptFit(111111);

  // This macro extracts information on the MET resolution by fitting the combination
  // of several TProfiles
  // --------------------------------------------------------------------------------
  TFile * QCD[8]; 
  QCD[0] = new TFile("./root/TDAna_QCD30-50_tk3.root");
  QCD[1] = new TFile("./root/TDAna_QCD50-80_tk3.root");
  QCD[2] = new TFile("./root/TDAna_QCD80-120_tk3.root");
  QCD[3] = new TFile("./root/TDAna_QCD120-170_tk3.root");
  QCD[4] = new TFile("./root/TDAna_QCD170-230_tk3.root");
  QCD[5] = new TFile("./root/TDAna_QCD230-300_tk3.root");
  QCD[6] = new TFile("./root/TDAna_QCD300-380_tk3.root");
  QCD[7] = new TFile("./root/TDAna_QCD380incl_tk3.root");
  double QCDxs[8] = { 155929000., 20938850., 2949713., 499656., 100995.,  23855., 6391., 2821.};
  double N[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};

  //double QCDxs[8] = { 1., 1., 1., 1., 1., 1., 1., 1. };
  //double N[8] = { 1., 1., 1., 1., 1., 1., 1., 1. };

  TH1F * Res = new TH1F ( "Res", "MEt resolution vs SumEt", 20, 0., 4000. );

  TH1F * Slices[20];
  char nameh[30];
  for ( int i=0; i<20; i++ ) {
    sprintf (nameh,"Slices%d", i );
    Slices[i]= new TH1F ( nameh, nameh, 100, -200., 200. );
  }

  double tot[20][100]={0.};
  double s2_tot[20][100]={0.};
  for ( int i=0; i<8; i++ ) {
    QCD[i]->cd();
    for ( int ibin=1; ibin<=100; ibin++ ) {
      int ibin2=(ibin-1)/5;
      for ( int jbin=1; jbin<=100; jbin++ ) {
	double t=MEx_SumEt->GetBinContent(ibin,jbin);
	tot[ibin2][jbin]+=t*QCDxs[i]/N[i];
	s2_tot[ibin2][jbin]+=t*pow(QCDxs[i]/N[i],2);
      }
    }
  }

  for ( int ibin2=0; ibin2<20; ibin2++ ) {
    for ( int jbin=1; jbin<=100; jbin++ ) {
      Slices[ibin2]->SetBinContent(jbin,tot[ibin2][jbin]);
      Slices[ibin2]->SetBinError(jbin,sqrt(s2_tot[ibin2][jbin]));
    }
  }

  TCanvas * c = new TCanvas ("c", "MEt fits", 600, 600 );
  c->Divide(4,4);
  TF1 * G = new TF1 ( "G", "[0]*exp(-0.5*pow(x/[1],2))", -200., 200. );
  G->SetLineWidth(1);
  G->SetLineColor(kRed);
  for ( int ibin2=5; ibin2<20; ibin2++ ) {
    c->cd(ibin2-4);
    G->SetParameters(10000., 20.+2.4*(double)ibin2 );
    Slices[ibin2]->Fit("G","V","",-200.,200.);
    Res->SetBinContent(ibin2,G->GetParameter(1));
    Res->SetBinError(ibin2,G->GetParError(1));
  }

  TF1 * Q = new TF1 ( "Q", "[0]*pow(x,[1])", 0., 4000. );
  Q->SetParameters(0.05,0.8);
  Q->SetLineWidth(1);
  Q->SetLineColor(kRed);

  TCanvas * b = new TCanvas ("b", "missing Et comparisons", 700, 700 );
  b->cd();
  Res->Fit("Q","V","",800,4000);
  Res->Draw("PE");
  b->Print("MetResFits3_g.ps");

  c->cd(16);
  Res->Draw("PE");
  c->Print("MetResFits3_s.ps");

}
