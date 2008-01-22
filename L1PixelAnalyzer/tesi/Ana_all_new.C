void Ana_all_new (TString var, TString sel, double mincut=-10000, double maxcut=10000 ) 
{

  // - Controllare tutti i numeri eventi
  // - Sezione d'urto single top ? Modificare!

  const int nbins=50;
  TString pippo = var+sel;

  TFile * QCD[16];     
  QCD[0] = new TFile("./root/TDAna_QCD_30-50.root");
  QCD[1] = new TFile("./root/TDAna_QCD_50-80.root");
  QCD[2] = new TFile("./root/TDAna_QCD_80-120.root");
  QCD[3] = new TFile("./root/TDAna_QCD_120-170.root");
  QCD[4] = new TFile("./root/TDAna_QCD_170-230.root");
  QCD[5] = new TFile("./root/TDAna_QCD_230-300.root");
  QCD[6] = new TFile("./root/TDAna_QCD_300-380.root");
  QCD[7] = new TFile("./root/TDAna_QCD_380-incl.root");
  QCD[8] = new TFile("./root/TDAna_EXTRA_QCD_30-50.root");
  QCD[9] = new TFile("./root/TDAna_EXTRA_QCD_50-80.root");
  QCD[10] = new TFile("./root/TDAna_EXTRA_QCD_80-120.root");
  QCD[11] = new TFile("./root/TDAna_EXTRA_QCD_120-170.root");
  QCD[12] = new TFile("./root/TDAna_EXTRA_QCD_170-230.root");
  QCD[13] = new TFile("./root/TDAna_EXTRA_QCD_230-300.root");
  QCD[14] = new TFile("./root/TDAna_EXTRA_QCD_300-380.root");
  QCD[15] = new TFile("./root/TDAna_EXTRA_QCD_380-incl.root");
  double QCDxs[16] = { 155929000., 20938850., 2949713., 499656., 
		       100995.,  23855., 6391., 2821., 155929000., 
		       20938850., 2949713., 499656., 100995.,  
		       23855., 6391., 2821. };
  double NQCD[16] = { 170000., 82000., 104000., 96000., 100000., 
		      102000., 112000., 102000., 130000., 282000.,
		      282000., 312000., 436000., 346000., 
		      418000., 406000. }; // last 4 to be checked
  // Sum QCD numbers so that we do not need to change the code below:
  // by adding numbers for each xs we can just sum all samples as always
  // -------------------------------------------------------------------
  for ( int i=0; i<16; i++ ) {
    if ( i<8 ) {
      NQCD[i] = NQCD[i]+NQCD[i+8];
    } else {
      NQCD[i] = NQCD[i-8]; // they have already changed, so this is correct
    }
  }

  TFile * QCDMULT[8];     
  QCDMULT[0] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_30-50.root");
  QCDMULT[1] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_50-80.root");
  QCDMULT[2] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_80-120.root");
  QCDMULT[3] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_120-170.root");
  QCDMULT[4] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_170-230.root");
  QCDMULT[5] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_230-300.root");
  QCDMULT[6] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_300-380.root");
  QCDMULT[7] = new TFile("./root/TDAna_MULTI_EXTRA_QCD_380-incl.root");
  double QCDMxs[8] =  { 155929000., 20938850., 2949713., 499656., 
			100995.,  23855., 6391., 2821. };
  double NQCDM[8] = { 11229443., 26347818., 24971508.,  29514603., 
		      40608575., 32268604., 37943909., 33232000. };

  TFile * W[11];
  W[0] = new TFile ("./root/TDAna_W_0JETS.root");
  W[1] = new TFile ("./root/TDAna_W_1JETS_0ptw100.root");
  W[2] = new TFile ("./root/TDAna_W_1JETS_100ptw300.root");
  W[3] = new TFile ("./root/TDAna_W_2JETS_0ptw100.root");
  W[4] = new TFile ("./root/TDAna_W_2JETS_100ptw300.root");
  W[5] = new TFile ("./root/TDAna_W_3JETS_0ptw100.root"); 
  W[6] = new TFile ("./root/TDAna_W_3JETS_100ptw300.root");
  W[7] = new TFile ("./root/TDAna_W_4JETS_0ptw100.root");
  W[8] = new TFile ("./root/TDAna_W_4JETS_100ptw300.root");
  W[9] = new TFile ("./root/TDAna_W_5JETS_0ptw100.root");
  W[10] = new TFile ("./root/TDAna_W_5JETS_100ptw300.root");
  double Wxs[11] = { 45000., 9200., 250., 2500., 225., 590., 
		     100., 125., 40., 85., 40. };
  double NW[11] = { 88000., 40000., 100530., 99523., 105255., 79000., 
		    88258., 83038., 30796., 59022., 41865. };

  TFile * TTH = new TFile("./root/TDAna_TTH_120.root");
  double TTHxs = 0.667 ;
  double NTTH = 1634000.; // 1652000.; // 96000.;

  TFile * TT[5];
  TT[0] = new TFile("./root/TDAna_TT_0JETS.root");
  TT[1] = new TFile("./root/TDAna_TT_1JETS.root");
  TT[2] = new TFile("./root/TDAna_TT_2JETS.root");
  TT[3] = new TFile("./root/TDAna_TT_3JETS.root");
  TT[4] = new TFile("./root/TDAna_TT_4JETS.root");
  // double TTxs[5] = { 619., 176., 34.,  6., 1.5 };  // from web
  double TTxs[5] = { 434., 162., 43., 10., 1.9 };     // from note
  double NTT[5] = { 57900., 66000., 98159., 14768., 5352. };

  TFile * TTPYT;
  TTPYT = new TFile("./root/TDAna_TTBAR.root");
  double TTPYTxs = 650.9;
  double NTTPYT = 972000.;

  TFile * ST[3];
  ST[0] = new TFile("./root/TDAna_SINGLETOP_TQ_TQB_LHC_E.root");
  ST[1] = new TFile("./root/TDAna_SINGLETOP_TQ_TQB_LHC_MU.root");
  ST[2] = new TFile("./root/TDAna_SINGLETOP_TQ_TQB_LHC_TAU.root");
  double STxs[3] = { 27.43, 26.97, 28.71 };
  double NST[3] = { 92000, 94000, 94000 }; 

  TFile * ZNN[2];
  ZNN[0] = new TFile("./root/TDAna_ZNUNUJETS_120-170.root");
  ZNN[1] = new TFile("./root/TDAna_ZNUNUJETS_170-230.root");
  double ZNNxs[2] = { 51.47, 15.52 };
  double NZNN[2] = { 29897., 25600. };

  ///////////////////////////////////////////////////////////////////

  double Lumfactor = 100000; // 100/fb of luminosity assumed in histograms
  
  TH1D * H = dynamic_cast<TH1D*>(TTH->Get(pippo));
  double minx=H->GetBinLowEdge(1);
  double maxx=nbins*H->GetBinWidth(1)+minx;

  int imin=(int)((mincut-minx)/(maxx-minx)*50);
  int imax=(int)((maxcut-minx)/(maxx-minx)*50);
  if ( imin<1 ) imin=1;
  if ( imax>nbins ) imax=nbins;

  TH1D * Histo_QCD = new TH1D ( pippo+"_QCD", "", nbins, minx, maxx );
  TH1D * Histo_QCDM = new TH1D ( pippo+"_QCDM", "", nbins, minx, maxx );
  TH1D * Histo_TTH = new TH1D ( pippo+"_TTH", "", nbins, minx, maxx );
  TH1D * Histo_TT = new TH1D ( pippo+"_TT", "", nbins, minx, maxx );
  TH1D * Histo_TTPYT = new TH1D ( pippo+"_TTPYT", "", nbins, minx, maxx );
  TH1D * Histo_W = new TH1D ( pippo+"_W", "", nbins, minx, maxx );
  TH1D * Histo_ST = new TH1D ( pippo+"_ST", "", nbins, minx, maxx );
  TH1D * Histo_ZNN = new TH1D ( pippo+"_ZNN", "", nbins, minx, maxx );

  TH1D * Histo_TOT = new TH1D ( pippo+"_TOT", "", nbins, minx, maxx );

  // Extract sum histograms with the right normalization and errors
  // --------------------------------------------------------------

  // QCD plus EXTRA QCD:
  // -------------------
  double totWQCD[16][nbins]={0.};
  double totQCD[nbins]={0.};
  double s2_totQCD[nbins]={0.};
  double totNQCD[16][nbins]={0.};
  for ( int i=0; i<16; i++ ) {
    cout << "Processing QCD file #" << i << " ..." << endl;
    if ( i==14 ) i++;
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
  for ( int i=0; i<16; i++ ) {
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      totQCD[ibin-1]+=totWQCD[i][ibin-1];
      if ( totNQCD[i][ibin-1]>0 ) {
	s2_totQCD[ibin-1]+=pow(totWQCD[i][ibin-1],2)/totNQCD[i][ibin-1];
      }
    }
  }
  // IF WE WANT TO RENORMALIZE QCD MATRIX-WEIGHTED HISTOGRAMS TO NUMBER EXPECTATION:
  // -------------------------------------------------------------------------------
//   if ( renormalize_toNqcd ) {
//     double tt=0.;
//     for ( int ibin=1; ibin<=nbins; ibin++ ) {
//       tw += totQCD[ibin-1];
//       tn += totNQCD[i]
//     }
//     for ( int ibin=1; ibin<=nbins; ibin++ ) {
//       totQCD[ibin-1];
//     }
//   }
  
  // QCD Multiplied:
  // ---------------
  double totQCDM[nbins]={0.};
  double s2_totQCDM[nbins]={0.};
  for ( int i=0; i<8; i++ ) { 
    cout << "Processing QCD MULT file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(QCDMULT[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totQCDM[ibin-1]+=t*QCDxs[i]/NQCDM[i]*Lumfactor;
      s2_totQCDM[ibin-1]+=t*pow(QCDMxs[i]/NQCDM[i]*Lumfactor,2);
    }
  }

  // TTH:
  // ----
  double totTTH[nbins]={0.};
  double s2_totTTH[nbins]={0.};
  cout << "Processing TTH file " << " ..." << endl;
  TH1D * Histo = dynamic_cast<TH1D*>(TTH->Get(pippo));
  for ( int ibin=1; ibin<=nbins; ibin++ ) {
    double t=Histo->GetBinContent(ibin);
    totTTH[ibin-1]+=t*TTHxs/NTTH*Lumfactor;
    s2_totTTH[ibin-1]+=t*pow(TTHxs/NTTH*Lumfactor,2);
  }

  // TT AlpGen:
  // ----------
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

  // TTPYT:
  // ------
  double totTTPYT[nbins]={0.};
  double s2_totTTPYT[nbins]={0.};
  cout << "Processing TTPYT file " << " ..." << endl;
  TH1D * Histo = dynamic_cast<TH1D*>(TTPYT->Get(pippo));
  for ( int ibin=1; ibin<=nbins; ibin++ ) {
    double t=Histo->GetBinContent(ibin);
    totTTPYT[ibin-1]+=t*TTPYTxs/NTTPYT*Lumfactor;
    s2_totTTPYT[ibin-1]+=t*pow(TTPYTxs/NTTPYT*Lumfactor,2);
  }

  // W files:
  // --------
  double totWW[11][nbins]={0.};
  double totW[nbins]={0.};
  double s2_totW[nbins]={0.};
  double totNW[11][nbins]={0.};
  for ( int i=0; i<11; i++ ) {
    cout << "Processing W file #" << i << " ..." << endl;
    if ( i==5 ) i++;
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

//   // ST :
//   // ----
   double totST[nbins]={0.};
   double s2_totST[nbins]={0.};
  for ( int i=0; i<3; i++ ) {
    cout << "Processing ST file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(ST[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totST[ibin-1]+=t*STxs[i]/NST[i]*Lumfactor;
      s2_totST[ibin-1]+=t*pow(STxs[i]/NST[i]*Lumfactor,2);
    }
  }

  // ZNN :
  // -----
  double totZNN[nbins]={0.};
  double s2_totZNN[nbins]={0.};
  for ( int i=0; i<2; i++ ) {
    cout << "Processing ZNN file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(ZNN[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totZNN[ibin-1]+=t*ZNNxs[i]/NZNN[i]*Lumfactor;
      s2_totZNN[ibin-1]+=t*pow(ZNNxs[i]/NZNN[i]*Lumfactor,2);
    }
  }

  // OK now fill total histograms
  // ----------------------------
  for ( int ibin=1; ibin<=nbins; ibin++ ) {
    Histo_QCD->SetBinContent(ibin,totQCD[ibin-1]);
    Histo_QCDM->SetBinContent(ibin,totQCDM[ibin-1]); // not in total
    Histo_TTH->SetBinContent(ibin,totTTH[ibin-1]);
    Histo_TT->SetBinContent(ibin,totTT[ibin-1]);
    Histo_TTPYT->SetBinContent(ibin,totTTPYT[ibin-1]); // not in total
    Histo_W->SetBinContent(ibin,totW[ibin-1]);
    Histo_ST->SetBinContent(ibin,totST[ibin-1]);
    Histo_ZNN->SetBinContent(ibin,totZNN[ibin-1]);
    Histo_QCD->SetBinError(ibin,sqrt(s2_totQCD[ibin-1]));
    Histo_QCDM->SetBinError(ibin,sqrt(s2_totQCDM[ibin-1])); // not in total
    Histo_TTH->SetBinError(ibin,sqrt(s2_totTTH[ibin-1]));
    Histo_TT->SetBinError(ibin,sqrt(s2_totTT[ibin-1]));
    Histo_TTPYT->SetBinError(ibin,sqrt(s2_totTTPYT[ibin-1])); // not in total
    Histo_W->SetBinError(ibin,sqrt(s2_totW[ibin-1]));
    Histo_ST->SetBinError(ibin,sqrt(s2_totST[ibin-1]));
    Histo_ZNN->SetBinError(ibin,sqrt(s2_totZNN[ibin-1]));
    double grandtot = (totQCD[ibin-1]+totTTH[ibin-1]+totTT[ibin-1]+
		       totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    double grandtote= sqrt(s2_totQCD[ibin-1]+s2_totTTH[ibin-1]+s2_totTT[ibin-1]+
			   s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT->SetBinContent(ibin,grandtot);
    Histo_TOT->SetBinError(ibin,grandtote);
  }

  double B_int = Histo_TOT->Integral(imin,imax);
  double S_int = Histo_TTH->Integral(imin,imax);
  cout << "Selected events: " << B_int << "B, " << S_int << "S. Significance = " << S_int/sqrt(B_int) << endl;
  cout << "QCD comparison: P= " << Histo_QCD->KolmogorovTest(Histo_QCDM) << endl;
  cout << "TOP comparison: P= " << Histo_TT->KolmogorovTest(Histo_TTPYT) << endl;
  
  TCanvas * QCDcomp = new TCanvas ( "QCDcomp", "Comparison QCD / QCD mult", 600, 600 );
  QCDcomp->Divide(1,2);
  QCDcomp->cd(1);
  QCDcomp->GetPad(1)->SetLogy();
  Histo_QCD->SetMarkerStyle(21);
  Histo_QCD->SetMarkerSize(0.4);
  Histo_QCD->SetMarkerColor(kRed);
  Histo_QCD->SetLineColor(kRed);
  Histo_QCD->DrawCopy("PE");
  Histo_QCDM->SetMarkerStyle(21);
  Histo_QCDM->SetMarkerSize(0.4);
  Histo_QCDM->SetMarkerColor(kBlue);
  Histo_QCDM->SetLineColor(kBlue);
  Histo_QCDM->DrawCopy("PESAME");
  QCDcomp->cd(2);
  Histo_QCD->SetMarkerStyle(21);
  Histo_QCD->SetMarkerSize(0.4);
  Histo_QCD->SetMarkerColor(kRed);
  Histo_QCD->SetLineColor(kRed);
  Histo_QCD->DrawCopy("PE");
  Histo_QCDM->SetMarkerStyle(21);
  Histo_QCDM->SetMarkerSize(0.4);
  Histo_QCDM->SetMarkerColor(kBlue);
  Histo_QCDM->SetLineColor(kBlue);
  Histo_QCDM->DrawCopy("PESAME");
  QCDcomp->Print("./ps/"+pippo+"_QCDcomp.ps");

  TCanvas * TTcomp = new TCanvas ( "TTcomp", "Comparison TT AlpGen / TT Pythia", 600, 600 );
  TTcomp->Divide(1,2);
  TTcomp->cd(1);
  TTcomp->GetPad(1)->SetLogy();
  Histo_TT->SetMarkerStyle(21);
  Histo_TT->SetMarkerSize(0.4);
  Histo_TT->SetMarkerColor(kRed);
  Histo_TT->SetLineColor(kRed);
  Histo_TT->DrawCopy("PE");
  Histo_TTPYT->SetMarkerStyle(21);
  Histo_TTPYT->SetMarkerSize(0.4);
  Histo_TTPYT->SetMarkerColor(kBlue);
  Histo_TTPYT->SetLineColor(kBlue);
  Histo_TTPYT->DrawCopy("PESAME");
  TTcomp->cd(2);
  Histo_TT->SetMarkerStyle(21);
  Histo_TT->SetMarkerSize(0.4);
  Histo_TT->SetMarkerColor(kRed);
  Histo_TT->SetLineColor(kRed);
  Histo_TT->DrawCopy("PE");
  Histo_TTPYT->SetMarkerStyle(21);
  Histo_TTPYT->SetMarkerSize(0.4);
  Histo_TTPYT->SetMarkerColor(kBlue);
  Histo_TTPYT->SetLineColor(kBlue);
  Histo_TTPYT->DrawCopy("PESAME");
  TTcomp->Print("./ps/"+pippo+"_TTcomp.ps");

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
  Histo_ST->SetMarkerStyle(27);
  Histo_ST->SetMarkerSize(0.4);
  Histo_ST->SetMarkerColor(kBlack);
  Histo_ST->SetLineColor(kBlack);
  Histo_ST->DrawCopy("SAMEHISTO");
  Histo_ZNN->SetMarkerStyle(27);
  Histo_ZNN->SetMarkerSize(0.4);
  Histo_ZNN->SetMarkerColor(kGreen);
  Histo_ZNN->SetLineColor(kGreen);
  Histo_ZNN->DrawCopy("SAMEHISTO");
  Histo_TOT->DrawCopy("PESAME");
  b->cd(2);
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
  Histo_ST->SetMarkerStyle(27);
  Histo_ST->SetMarkerSize(0.4);
  Histo_ST->SetMarkerColor(kBlack);
  Histo_ST->SetLineColor(kBlack);
  Histo_ST->DrawCopy("SAMEHISTO");
  Histo_ZNN->SetMarkerStyle(27);
  Histo_ZNN->SetMarkerSize(0.4);
  Histo_ZNN->SetMarkerColor(kGreen);
  Histo_ZNN->SetLineColor(kGreen);
  Histo_ZNN->DrawCopy("SAMEHISTO");
  Histo_TOT->DrawCopy("PESAME");
  b->Print("./ps/"+pippo+".ps");

  // Histograms to fit:
  // ------------------
  TFile * xmia = new TFile("xmia.root","RECREATE");
  xmia->cd();
  Histo_TOT->Write();
  Histo_TTH->Write();
  xmia->Close();


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
  
  ////////////////////////////////////////////////////////
  
//   // Now compute max quality factor for this distribution
//   // -----------------------------------------------------
//   TH1F * S_ = new TH1F ( "S", "S", 50, minx, maxx );
//   TH1F * B_ = new TH1F ( "B", "B", 50, minx, maxx );
//   TH1F * ES_ = new TH1F ( "ES", "ES", 50, minx, maxx );
//   TH1F * EB_ = new TH1F ( "EB", "EB", 50, minx, maxx);  
//   TH1F * R_ = new TH1F ( "R", "R", 50, minx, maxx );
//   TH1F * Q_ = new TH1F ( "Q", "Q", 50, minx, maxx );
  
//   double IS = Histo_TTH->Integral();
//   double IB = Histo_TOT->Integral();
//   double es[nbins];
//   double s_es[nbins];
//   double eb[nbins];
//   double s_eb[nbins];
//   double maxS=0;
//   double maxB=0;
//   for ( int ibin=1; ibin<=nbins; ibin++ ) {
//     S_->SetBinContent(ibin,Histo_TTH->GetBinContent(ibin)/IS);
//     S_->SetBinError(ibin,Histo_TTH->GetBinError(ibin)/IS);
//     if ( maxS<Histo_TTH->GetBinContent(ibin)/IS ) maxS=Histo_TTH->GetBinContent(ibin)/IS;
//     B_->SetBinContent(ibin,Histo_TOT->GetBinContent(ibin)/IB);
//     B_->SetBinError(ibin,Histo_TOT->GetBinError(ibin)/IB);
//     if ( maxB<Histo_TOT->GetBinContent(ibin)/IB ) maxB=Histo_TOT->GetBinContent(ibin)/IB;
//     double fs=0;    // signal failing cut
//     double s2_fs=0;
//     double ps=0;    // signal passing cut
//     double s2_ps=0;
//     double fb=0;    // background failing cut
//     double s2_fb=0;
//     double pb=0;    // background passing cut
//     double s2_pb=0;
//     for ( int jbin=1; jbin<=nbins; jbin++ ) {
//       if ( jbin<ibin ) { 
// 	fs    += Histo_TTH->GetBinContent(jbin);
// 	s2_fs += Histo_TTH->GetBinError(jbin);
//       }
//       if ( jbin>=ibin ) {
// 	ps    += Histo_TTH->GetBinContent(jbin);
// 	s2_ps += Histo_TTH->GetBinError(jbin);
//       }
//     }
//     for ( int jbin=1; jbin<=nbins; jbin++ ) {
//       if ( jbin<ibin ) { 
// 	fb    += Histo_TOT->GetBinContent(jbin);
// 	s2_fb += Histo_TOT->GetBinError(jbin);
//       }
//       if ( jbin>=ibin ) {
// 	pb    += Histo_TOT->GetBinContent(jbin);
// 	s2_pb += Histo_TOT->GetBinError(jbin);
//       }
//     }
//     if ( fs+ps>0 ) {
//       es[ibin-1]=ps/(fs+ps);
//       s_es[ibin-1]=sqrt(ps*ps*s2_fs+fs*fs*s2_ps)/pow(fs+ps,2);
//     } else {
//       es[ibin-1]=0;
//       s_es[ibin-1]=0;
//     }
//     ES->SetBinContent(ibin,es[ibin-1]);
//     ES->SetBinError(ibin,s_es[ibin-1]);
//     if ( fb+pb>0 ) {
//       eb[ibin-1]=pb/(fb+pb);
//       s_eb[ibin-1]=sqrt(pb*pb*s2_fb+fb*fb*s2_pb)/pow(fb+pb,2);
//     } else {
//       eb[ibin-1]=0;
//       s_eb[ibin-1]=0;
//     }
//     EB->SetBinContent(ibin,eb[ibin-1]);
//     EB->SetBinError(ibin,s_eb[ibin-1]);
//   }
  
//   double R;
//   double s_R;
//   double Q;
//   double s_Q;
//   double maxQm1s=0;
//   double s_maxQm1s=0;
//   double x_maxQ=0;
//   for ( int ibin=1; ibin<=nbins; ibin++ ) {
//     if ( eb[ibin-1]>0 ) {
//       R = es[ibin-1]/eb[ibin-1];
//       s_R = sqrt(s_es[ibin-1]*s_es[ibin-1]/eb[ibin-1]/eb[ibin-1]+
// 		 es[ibin-1]*es[ibin-1]*s_eb[ibin-1]*s_eb[ibin-1]/pow(eb[ibin-1],4));
//       Q = es[ibin-1]*es[ibin-1]/eb[ibin-1];
//       s_Q = sqrt(pow(2*es[ibin-1]/eb[ibin-1]*s_es[ibin-1],2)+pow(es[ibin-1]/eb[ibin-1],4)*pow(s_eb[ibin-1],2));
//     } else {
//       R = es[ibin-1]/0.1;
//       s_R = sqrt(pow(s_es[ibin-1]/0.1,2)+
// 		 pow(es[ibin-1]*s_eb[ibin-1],2)/pow(0.1,4));
//       Q = es[ibin-1]*es[ibin-1]/0.1;
//       s_Q = sqrt(pow(2*es[ibin-1]/0.1*s_es[ibin-1],2)+pow(es[ibin-1]/0.1,4)*pow(s_eb[ibin-1],2));
//     }
//     R_->SetBinContent(ibin,R);
//     R_->SetBinError(ibin,s_R);
//     Q_->SetBinContent(ibin,Q);
//     Q_->SetBinError(ibin,s_Q);
//     if ( maxQm1s<Q-s_Q ) {
//       maxQm1s=Q;
//       s_maxQm1s=s_Q;
//       x_maxQ=(double)ibin*(maxx-minx)/(double)nbins;
//     } 
//   }
//   double maxX=maxS;
//   if (maxB>maxX) maxX=maxB;
//   cout << endl;
//   cout << "Maximum quality factor: " << maxQm1s << "+-" << s_maxQm1s << " at x > " << x_maxQ << endl; 
//   cout << "--------------------------------------------------------------------" << endl;
  
//   TCanvas * b = new TCanvas ("b", "Cut optimization", 500, 500 );
//   b->Divide(2,2);
//   b->cd(1);
//   S_->SetMaximum(1.2*maxX);
//   S_->SetLineColor(kRed);
//   S_->SetMarkerColor(kRed);
//   S_->Draw("PE");
//   B_->SetLineColor(kBlue);
//   B_->SetMarkerColor(kBlue);
//   B_->Draw("PESAME");
//   b->cd(2);
//   ES_->SetLineColor(kRed);
//   ES_->SetMarkerColor(kRed);
//   ES_->Draw("PE");
//   EB_->SetLineColor(kBlue);
//   EB_->SetMarkerColor(kBlue);
//   EB_->Draw("PESAME");
//   b->cd(3);
//   R_->Draw("PE");
//   b->cd(4);
//   Q_->Draw("PE");
//   b->Print(pippo+"_opt.ps");


}
