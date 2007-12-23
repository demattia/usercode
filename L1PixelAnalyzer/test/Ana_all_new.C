void Ana_all (TString var, TString sel ) 
{

  const int nbins=50;
  TString pippo = var+sel;

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

  TFile * TTH = new TFile("./root/TDAna_ttH_120_tk3.root");
  double TTHxs = 0.667 ;
  double NTTH = 62000.;

  TFile * TT[5];
  TT[0] = new TFile("./root/TDAna_TT0_tk3.root");
  TT[1] = new TFile("./root/TDAna_TT1_tk3.root");
  TT[2] = new TFile("./root/TDAna_TT2_tk3.root");
  TT[3] = new TFile("./root/TDAna_TT3_tk3.root");
  TT[4] = new TFile("./root/TDAna_TT4_tk3.root");
  // double TTxs[5] = { 619., 176., 34.,  6., 1.5 };  // from web
  double TTxs[5] = { 434., 162., 43., 10., 1.9 };     // from note
  double NTT[5] = { 57900., 66000., 98159., 14768., 5304. };

  double Lumfactor = 100000; // 100/fb of luminosity assumed in histograms
  
  TH1D * H = dynamic_cast<TH1D*>(TTH->Get(pippo));
  double minx=H->GetBinLowEdge(1);
  double maxx=nbins*H->GetBinWidth(1)+minx;
  TH1D * Histo_QCD = new TH1D ( pippo+"_QCD", "", nbins, minx, maxx );
  TH1D * Histo_TTH = new TH1D ( pippo+"_TTH", "", nbins, minx, maxx );
  TH1D * Histo_TT = new TH1D ( pippo+"_TT", "", nbins, minx, maxx );
  TH1D * Histo_W = new TH1D ( pippo+"_W", "", nbins, minx, maxx );

  TH1D * Histo_TOT = new TH1D ( pippo+"_TOT", "", nbins, minx, maxx );

  // Extract sum histograms with the right normalization and errors
  // --------------------------------------------------------------
  double totWW[11][nbins]={0.};
  double totW[nbins]={0.};
  double s2_totW[nbins]={0.};
  double totNW[11][nbins]={0.};
  for ( int i=0; i<11; i++ ) {
    cout << "Processing W file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(W[i]->Get(pippo));
    TH1D * HistoW = dynamic_cast<TH1D*>(W[i]->Get(pippo+"W"));
    // For W, we need also total entries in histograms to add a
    // Poisson fluke contribution to total errors from matrix:
    // ----------------------------------------------------------
    TH1D * HistoN = dynamic_cast<TH1D*>(W[i]->Get(pippo+"N"));    
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      double s2t=HistoW->GetBinContent(ibin);
      double n=HistoN->GetBinContent(ibin);
      totWW[i][ibin-1]+=t*Wxs[i]/NW[i]*Lumfactor;
      s2_totW[ibin-1]+=s2t*pow(Wxs[i]/NW[i]*Lumfactor,2);
      totNW[i][ibin-1]+=n;
    }
  }
  // Once grandtotals of weights are computed for each bin, we can
  // add to the total s2 the Poisson contribution 1/sqrt(N) * T
  // -------------------------------------------------------------
  for ( int i=0; i<11; i++ ) {
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      totW[ibin-1]+=totWW[i][ibin-1];
      if ( totNW[i][ibin-1]>0 ) {
	s2_totW[ibin-1]+=pow(totWW[i][ibin-1],2)/totNW[i][ibin-1];
      }
    }
  }

  double totWQCD[8][nbins]={0.};
  double totQCD[nbins]={0.};
  double s2_totQCD[nbins]={0.};
  double totNQCD[8][nbins]={0.};
  for ( int i=0; i<8; i++ ) {
    cout << "Processing QCD file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(QCD[i]->Get(pippo));
    TH1D * HistoW = dynamic_cast<TH1D*>(QCD[i]->Get(pippo+"W"));
    // For QCD, we need also total entries in histograms to add a
    // Poisson fluke contribution to total errors from matrix:
    // ----------------------------------------------------------
    TH1D * HistoN = dynamic_cast<TH1D*>(QCD[i]->Get(pippo+"N"));    
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      double s2t=HistoW->GetBinContent(ibin);
      double n=HistoN->GetBinContent(ibin);
      totWQCD[i][ibin-1]+=t*QCDxs[i]/NQCD[i]*Lumfactor;
      s2_totQCD[ibin-1]+=s2t*pow(QCDxs[i]/NQCD[i]*Lumfactor,2);
      totNQCD[i][ibin-1]+=n;
    }
  }
  // Once grandtotals of weights are computed for each bin, we can
  // add to the total s2 the Poisson contribution 1/sqrt(N) * T
  // -------------------------------------------------------------
  for ( int i=0; i<8; i++ ) {
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      totQCD[ibin-1]+=totWQCD[i][ibin-1];
      if ( totNQCD[i][ibin-1]>0 ) {
	s2_totQCD[ibin-1]+=pow(totWQCD[i][ibin-1],2)/totNQCD[i][ibin-1];
      }
    }
  }

  double totTT[nbins]={0.};
  double s2_totTT[nbins]={0.};
  for ( int i=0; i<5; i++ ) {
    cout << "Processing TT file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(TT[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totTT[ibin-1]+=t*TTxs[i]/NTT[i]*Lumfactor;
      s2_totTT[ibin-1]+=t*pow(TTxs[i]/NTT[i]*Lumfactor,2);
    }
  }
  double totTTH[nbins]={0.};
  double s2_totTTH[nbins]={0.};
  cout << "Processing TTH file " << " ..." << endl;
  TH1D * Histo = dynamic_cast<TH1D*>(TTH->Get(pippo));
  for ( int ibin=1; ibin<=nbins; ibin++ ) {
    double t=Histo->GetBinContent(ibin);
    totTTH[ibin-1]+=t*TTHxs/NTTH*Lumfactor;
    s2_totTTH[ibin-1]+=t*pow(TTHxs/NTTH*Lumfactor,2);
  }

  // OK now fill total histograms
  // ----------------------------
  for ( int ibin=1; ibin<=nbins; ibin++ ) {
    Histo_QCD->SetBinContent(ibin,totQCD[ibin-1]);
    Histo_TTH->SetBinContent(ibin,totTTH[ibin-1]);
    Histo_TT->SetBinContent(ibin,totTT[ibin-1]);
    Histo_W->SetBinContent(ibin,totW[ibin-1]);
    Histo_QCD->SetBinError(ibin,sqrt(s2_totQCD[ibin-1]));
    Histo_TTH->SetBinError(ibin,sqrt(s2_totTTH[ibin-1]));
    Histo_TT->SetBinError(ibin,sqrt(s2_totTT[ibin-1]));
    Histo_W->SetBinError(ibin,sqrt(s2_totW[ibin-1]));
    double grandtot = totQCD[ibin-1]+totTTH[ibin-1]+totTT[ibin-1]+totW[ibin-1];
    double grandtote= sqrt(s2_totQCD[ibin-1]+s2_totTTH[ibin-1]+s2_totTT[ibin-1]+s2_totW[ibin-1]);
    Histo_TOT->SetBinContent(ibin,grandtot);
    Histo_TOT->SetBinError(ibin,grandtote);
  }

  TCanvas * b = new TCanvas ("b", "Kinematics comparison", 700, 700 );
  b->Divide(1,2);

  b->cd(1);
  b->GetPad(1)->SetLogy();
  Histo_TOT->SetMarkerStyle(20);
  Histo_TOT->SetMarkerSize(0.4);
  Histo_TOT->SetMinimum(1.);
  Histo_TOT->DrawCopy("PE");
  Histo_QCD->SetMarkerStyle(21);
  Histo_QCD->SetMarkerSize(0.4);
  Histo_QCD->SetMarkerColor(kRed);
  Histo_QCD->SetLineColor(kRed);
  Histo_QCD->DrawCopy("PESAME");
  Histo_TTH->SetMarkerStyle(24);
  Histo_TTH->SetMarkerSize(0.4);
  Histo_TTH->SetMarkerColor(kBlue);
  Histo_TTH->SetLineColor(kBlue);
  Histo_TTH->DrawCopy("PESAME");
  Histo_TT->SetMarkerStyle(25);
  Histo_TT->SetMarkerSize(0.4);
  Histo_TT->SetMarkerColor(kGreen);
  Histo_TT->SetLineColor(kGreen);
  Histo_TT->DrawCopy("PESAME");
  Histo_W->SetMarkerStyle(26);
  Histo_W->SetMarkerSize(0.4);
  Histo_W->SetMarkerColor(kCyan);
  Histo_W->SetLineColor(kCyan);
  Histo_W->DrawCopy("PESAME");
  Histo_TOT->DrawCopy("PESAME");

  b->cd(2);
  Histo_TOT->SetNormFactor(1);
  Histo_TOT->Draw("PE");
  Histo_QCD->SetNormFactor(1);
  Histo_QCD->Draw("PESAME");
  Histo_TTH->SetNormFactor(1);
  Histo_TTH->Draw("PESAME");
  Histo_TT->SetNormFactor(1);
  Histo_TT->Draw("PESAME");
  Histo_W->SetNormFactor(1);
  Histo_W->Draw("PESAME");
  Histo_TOT->Draw("PESAME");
  b->Print(pippo+".ps");

  // Close files
  // -----------
  //   for ( int i=0; i<8; i++ ) {
  //     QCD[i]->Close();
  //     QCDOLD[i]->Close();
  //   }
  //   for ( int i=0; i<10; i++ ) {
  //     W[i]->Close();
  //   }
  //   for ( int i=0; i<5; i++ ) {
  //     TT[i]->Close();
  //   }
  //   TTH->Close();

}
