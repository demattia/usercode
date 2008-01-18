Compare_W ( TString pippo ) 
{

  TFile * W[11];
  W[0] = new TFile ("./root/TDAna_W0w_tk3.root");
  W[1] = new TFile ("./root/TDAna_W10w_tk3.root");
  W[2] = new TFile ("./root/TDAna_W11w_tk3.root");
  W[3] = new TFile ("./root/TDAna_W20w_tk3.root");
  W[4] = new TFile ("./root/TDAna_W21w_tk3.root");
  W[5] = new TFile ("./root/TDAna_W30w_tk3.root");
  W[6] = new TFile ("./root/TDAna_W31w_tk3.root");
  W[7] = new TFile ("./root/TDAna_W40w_tk3.root");
  W[8] = new TFile ("./root/TDAna_W41w_tk3.root");
  W[9] = new TFile ("./root/TDAna_W50w_tk3.root");
  W[10] = new TFile ("./root/TDAna_W51w_tk3.root");
  double Wxs[11] = { 45000., 9200., 250., 2500., 225., 590., 100., 125., 40., 85., 40. };
  double NW[11] = { 88000., 40000., 100530., 99523., 105255., 79000., 
		    88258., 83038., 30796., 59022., 41865. };

  TFile * WOLD[11];
  WOLD[0] = new TFile ("./root/TDAna_W0_tk3.root");
  WOLD[1] = new TFile ("./root/TDAna_W10_tk3.root");
  WOLD[2] = new TFile ("./root/TDAna_W11_tk3.root");
  WOLD[3] = new TFile ("./root/TDAna_W20_tk3.root");
  WOLD[4] = new TFile ("./root/TDAna_W21_tk3.root");
  WOLD[5] = new TFile ("./root/TDAna_W30_tk3.root");
  WOLD[6] = new TFile ("./root/TDAna_W31_tk3.root");
  WOLD[7] = new TFile ("./root/TDAna_W40_tk3.root");
  WOLD[8] = new TFile ("./root/TDAna_W41_tk3.root");
  WOLD[9] = new TFile ("./root/TDAna_W50_tk3.root");
  WOLD[10] = new TFile ("./root/TDAna_W51_tk3.root");

  double Lumfactor = 100000.;

  TH1D * H = dynamic_cast<TH1D*>(W[0]->Get(pippo));
  double minx=H->GetBinLowEdge(1);
  double maxx=50.*H->GetBinWidth(1);
  TH1D * Histo_W = new TH1D ( pippo+"_W", "", 50, minx, maxx );
  TH1D * R_W = new TH1D ( pippo+"_W", "", 50, minx, maxx );
  TH1F * Histo_WOLD = new TH1F ( pippo+"_WOLD", "", 50, minx,maxx );

  // Extract sum histograms with the right normalization and errors
  // --------------------------------------------------------------
  double totW[50]={0.};
  double s2_totW[50]={0.};
  for ( int i=0; i<11; i++ ) {
    if ( i!=4 ) {
    cout << "Processing W file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(W[i]->Get(pippo));
    TH1D * HistoW = dynamic_cast<TH1D*>(W[i]->Get(pippo+"W"));
    for ( int ibin=1; ibin<=50; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      double s2t=HistoW->GetBinContent(ibin);
      totW[ibin-1]+=t*Wxs[i]/NW[i]*Lumfactor;
      s2_totW[ibin-1]+=s2t*pow(Wxs[i]/NW[i]*Lumfactor,2);
    }
    }
  }
  double totWOLD[50]={0.};
  double totNWOLD[50]={0.};
  double s2_totWOLD[50]={0.};
  for ( int i=0; i<11; i++ ) {
    if ( i!=4 ) {
    cout << "Processing W OLD file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(WOLD[i]->Get(pippo));
    for ( int ibin=1; ibin<=50; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totWOLD[ibin-1]+=t*Wxs[i]/NW[i]*Lumfactor;
      totNWOLD[ibin-1]+=t;
      s2_totWOLD[ibin-1]+=t*pow(Wxs[i]/NW[i]*Lumfactor,2);
    }
    }
  }
  // Once grandtotals of weights are computed for each bin, we can
  // add to the total s2 the Poisson contribution 1/sqrt(N) * T
  // -------------------------------------------------------------
  for ( int ibin=1; ibin<=50; ibin++ ) {
    if ( totNWOLD[ibin-1]==0 ) totNWOLD[ibin-1]=1;
    s2_totW[ibin-1]+=pow(totW[ibin-1],2)/totNWOLD[ibin-1];
  }

  // OK now fill total histograms
  // ----------------------------
  double nW=0.;
  double s2_NW=0.;
  double nWOLD=0.;
  double s2_NWOLD=0.;
  for ( int ibin=1; ibin<=50; ibin++ ) {
    nW+=totW[ibin-1];
    s2_NW+=s2_totW[ibin-1];
    nWOLD+=totWOLD[ibin-1];
    s2_NWOLD+=s2_totWOLD[ibin-1];
    Histo_W->SetBinContent(ibin,totW[ibin-1]);
    Histo_WOLD->SetBinContent(ibin,totWOLD[ibin-1]);
    Histo_W->SetBinError(ibin,sqrt(s2_totW[ibin-1]));
    Histo_WOLD->SetBinError(ibin,sqrt(s2_totWOLD[ibin-1]));
    double R=1.;
    double s_R;
    if ( totWOLD[ibin-1]>0 && totW[ibin-1] ) {
      R = totWOLD[ibin-1]/totW[ibin-1];
      s_R = R*sqrt(s2_totW[ibin-1]/pow(totW[ibin-1],2)+
		   s2_totWOLD[ibin-1]/pow(totWOLD[ibin-1],2));
      cout << ibin-1 << " " << totW[ibin-1] << "+-" << sqrt(s2_totW[ibin-1])/totW[ibin-1] 
	   << " / " <<  totWOLD[ibin-1] <<  "+-" <<sqrt(s2_totWOLD[ibin-1])/totWOLD[ibin-1]  
	   << " = " << R << "+-" << s_R << endl;
    }
    R_W->SetBinContent(ibin,R);
    R_W->SetBinError(ibin,s_R);
  }
  cout << "Totals: N(seen) = " << nWOLD << "+-" << sqrt(s2_NWOLD) << endl;
  cout << "        N(pred) = " << nW    << "+-" << sqrt(s2_NW) << endl;

  TCanvas * b = new TCanvas ("b", "Kinematics comparison", 700, 700 );
  b->Divide(1,2);

  b->cd(1);
  //b->GetPad(1)->SetLogy();
  Histo_W->SetLineColor(kRed);
  Histo_W->Draw("PE");
  Histo_WOLD->SetLineColor(kBlue);
  Histo_WOLD->Draw("PESAME");
  b->cd(2);
  R_W->SetMinimum(0.);
  R_W->SetMaximum(4.);
  R_W->Draw("PE");

  b->Print(pippo+".ps");

}
