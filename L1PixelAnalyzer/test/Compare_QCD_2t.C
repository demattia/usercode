Compare_QCD ( TString pippo ) 
{

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
  double NQCD[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};

  TFile * QCDOLD[8];     
  QCDOLD[0] = new TFile("./rootold/TDAna_QCD30-50_tk3.root");
  QCDOLD[1] = new TFile("./rootold/TDAna_QCD50-80_tk3.root");
  QCDOLD[2] = new TFile("./rootold/TDAna_QCD80-120_tk3.root");
  QCDOLD[3] = new TFile("./rootold/TDAna_QCD120-170_tk3.root");
  QCDOLD[4] = new TFile("./rootold/TDAna_QCD170-230_tk3.root");
  QCDOLD[5] = new TFile("./rootold/TDAna_QCD230-300_tk3.root");
  QCDOLD[6] = new TFile("./rootold/TDAna_QCD300-380_tk3.root");
  QCDOLD[7] = new TFile("./rootold/TDAna_QCD380incl_tk3.root");

  double Lumfactor = 100000.;

  TH1D * H = dynamic_cast<TH1D*>(QCD[0]->Get(pippo));
  double minx=H->GetBinLowEdge(1);
  double maxx=50.*H->GetBinWidth(1);
  TH1D * Histo_QCD = new TH1D ( pippo+"_QCD", "", 50, minx, maxx );
  TH1D * R_QCD = new TH1D ( pippo+"_QCD", "", 50, minx, maxx );
  TH1F * Histo_QCDOLD = new TH1F ( pippo+"_QCDOLD", "", 50, minx,maxx );

  // Extract sum histograms with the right normalization and errors
  // --------------------------------------------------------------
  double totQCD[50]={0.};
  double s2_totQCD[50]={0.};
  for ( int i=0; i<8; i++ ) {
    cout << "Processing QCD file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(QCD[i]->Get(pippo));
    TH1D * HistoW = dynamic_cast<TH1D*>(QCD[i]->Get(pippo+"W"));
    for ( int ibin=1; ibin<=50; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      double s2t=HistoW->GetBinContent(ibin);
      totQCD[ibin-1]+=t*QCDxs[i]/NQCD[i]*Lumfactor;
      s2_totQCD[ibin-1]+=s2t*pow(QCDxs[i]/NQCD[i]*Lumfactor,2);
    }
  }
  double totQCDOLD[50]={0.};
  double totNQCDOLD[50]={0.};
  double s2_totQCDOLD[50]={0.};
  for ( int i=0; i<8; i++ ) {
    cout << "Processing QCD OLD file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(QCDOLD[i]->Get(pippo));
    for ( int ibin=1; ibin<=50; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totQCDOLD[ibin-1]+=t*QCDxs[i]/NQCD[i]*Lumfactor;
      totNQCDOLD[ibin-1]+=t;
      s2_totQCDOLD[ibin-1]+=t*pow(QCDxs[i]/NQCD[i]*Lumfactor,2);
    }
  }
  // Once grandtotals of weights are computed for each bin, we can
  // add to the total s2 the Poisson contribution 1/sqrt(N) * T
  // -------------------------------------------------------------
  for ( int ibin=1; ibin<=50; ibin++ ) {
    if ( totNQCDOLD[ibin-1]==0 ) totNQCDOLD[ibin-1]=1;
    s2_totQCD[ibin-1]+=pow(totQCD[ibin-1],2)/totNQCDOLD[ibin-1];
  }

  // OK now fill total histograms
  // ----------------------------
  double nQCD=0.;
  double s2_NQCD=0.;
  double nQCDOLD=0.;
  double s2_NQCDOLD=0.;
  for ( int ibin=1; ibin<=50; ibin++ ) {
    nQCD+=totQCD[ibin-1];
    s2_NQCD+=s2_totQCD[ibin-1];
    nQCDOLD+=totQCDOLD[ibin-1];
    s2_NQCDOLD+=s2_totQCDOLD[ibin-1];
    Histo_QCD->SetBinContent(ibin,totQCD[ibin-1]);
    Histo_QCDOLD->SetBinContent(ibin,totQCDOLD[ibin-1]);
    Histo_QCD->SetBinError(ibin,sqrt(s2_totQCD[ibin-1]));
    Histo_QCDOLD->SetBinError(ibin,sqrt(s2_totQCDOLD[ibin-1]));
    double R=1.;
    double s_R;
    if ( totQCDOLD[ibin-1]>0 && totQCD[ibin-1] ) {
      R = totQCDOLD[ibin-1]/totQCD[ibin-1];
      s_R = R*sqrt(s2_totQCD[ibin-1]/pow(totQCD[ibin-1],2)+
		   s2_totQCDOLD[ibin-1]/pow(totQCDOLD[ibin-1],2));
      cout << ibin-1 << " " << totQCD[ibin-1] << "+-" << sqrt(s2_totQCD[ibin-1])/totQCD[ibin-1] 
	   << " / " <<  totQCDOLD[ibin-1] <<  "+-" <<sqrt(s2_totQCDOLD[ibin-1])/totQCDOLD[ibin-1]  
	   << " = " << R << "+-" << s_R << endl;
    }
    R_QCD->SetBinContent(ibin,R);
    R_QCD->SetBinError(ibin,s_R);
  }
  cout << "Totals: N(seen) = " << nQCDOLD << "+-" << sqrt(s2_NQCDOLD) << endl;
  cout << "        N(pred) = " << nQCD    << "+-" << sqrt(s2_NQCD) << endl;

  TCanvas * b = new TCanvas ("b", "Kinematics comparison", 700, 700 );
  b->Divide(1,2);

  b->cd(1);
  //b->GetPad(1)->SetLogy();
  Histo_QCD->SetLineColor(kRed);
  Histo_QCD->Draw("PE");
  Histo_QCDOLD->SetLineColor(kBlue);
  Histo_QCDOLD->Draw("PESAME");
  b->cd(2);
  R_QCD->SetMinimum(0.);
  R_QCD->SetMaximum(4.);
  R_QCD->Draw("PE");

  b->Print(pippo+".ps");

}
