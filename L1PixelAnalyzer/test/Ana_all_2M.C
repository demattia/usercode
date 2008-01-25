void Ana_all_new (TString var, TString sel, bool dops, bool renormalize ) 
{

  double mincut=-10000.;
  double maxcut=+10000.;

  // - Controllare tutti i numeri eventi
  // - Sezione d'urto single top ? Modificare!

  const int nbins=50;
  TString pippo1[9];
  TString S[9] = { "0","1","2","3","4","5","6","7","8"};
  for ( int i = 0; i < 9; ++i ) {
    pippo1[i] = var+sel+S[i];
  }
  TString pippo = var+sel;

  TFile * QCD[16];     
  QCD[0] = new TFile("TDAna/TDAna_QCD_30-50.root");
  QCD[1] = new TFile("TDAna/TDAna_QCD_50-80.root");
  QCD[2] = new TFile("TDAna/TDAna_QCD_80-120.root");
  QCD[3] = new TFile("TDAna/TDAna_QCD_120-170.root");
  QCD[4] = new TFile("TDAna/TDAna_QCD_170-230.root");
  QCD[5] = new TFile("TDAna/TDAna_QCD_230-300.root");
  QCD[6] = new TFile("TDAna/TDAna_QCD_300-380.root");
  QCD[7] = new TFile("TDAna/TDAna_QCD_380-incl.root");
  QCD[8] = new TFile("TDAna/TDAna_EXTRA_QCD_30-50.root");
  QCD[9] = new TFile("TDAna/TDAna_EXTRA_QCD_50-80.root");
  QCD[10] = new TFile("TDAna/TDAna_EXTRA_QCD_80-120.root");
  QCD[11] = new TFile("TDAna/TDAna_EXTRA_QCD_120-170.root");
  QCD[12] = new TFile("TDAna/TDAna_EXTRA_QCD_170-230.root");
  QCD[13] = new TFile("TDAna/TDAna_EXTRA_QCD_230-300.root");
  QCD[14] = new TFile("TDAna/TDAna_EXTRA_QCD_300-380.root");
  QCD[15] = new TFile("TDAna/TDAna_EXTRA_QCD_380-incl.root");
//   double QCDxs[16] = { 155929000., 20938850., 2949713., 499656., 
// 		       100995.,  23855., 6391., 2821., 155929000., 
// 		       20938850., 2949713., 499656., 100995.,  
// 		       23855., 6391., 2821. };
  double QCDxs[16] = { 163000000., 21600000., 3080000., 494000., 101000., 24500., 6240., 2821., 
		       163000000., 21600000., 3080000., 494000., 101000., 24500., 6240., 2821. };
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
  // Numbers from ascii files "decay".
  // ---------------------------------
  // MultiJet trigger
  double QCDentries[16] = {0.,0.,0.,0.,1.,13.,21.,28.,
			   0.,0.,0.,1.,5.,33.,66.,149.};
  // PixelJet trigger
//  double QCDentries[16] = {0.,0.,0.,0.,2.,8.,18.,27.,
// 			     0.,0.,0.,1.,5.,30.,67.,141.};
  double sigmaQCD[16];
  for ( int i=0; i<16; i++ ) { 
    if (QCDentries[i]<=10) sigmaQCD[i] = (1 + TMath::Sqrt(QCDentries[i] + 0.75));
    else sigmaQCD[i] = 1./sqrt(QCDentries[i]);
  };
  

  TFile * QCDMULT[8];     
  QCDMULT[0] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_30-50.root");
  QCDMULT[1] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_50-80.root");
  QCDMULT[2] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_80-120.root");
  QCDMULT[3] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_120-170.root");
  QCDMULT[4] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_170-230.root");
  QCDMULT[5] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_230-300.root");
  QCDMULT[6] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_300-380.root");
  QCDMULT[7] = new TFile("TDAna/TDAna_MULTI_EXTRA_QCD_380-incl.root");
  double QCDMxs[8] =  { 163000000., 21600000., 3080000., 494000., 
			101000., 245000., 6240., 2821. };
  //                    155929000., 20938850., 2949713., 499656., 
  //			100995.,  23855., 6391., 2821. };
  double NQCDM[8] = { 11229443., 26347818., 24971508.,  29514603., 
		      40608575., 32268604., 37943909., 33232000. };
  double QCDMentries[8]={0.};

  TFile * W[11];
  W[0] = new TFile ("TDAna/TDAna_W_0JETS.root");
  W[1] = new TFile ("TDAna/TDAna_W_1JETS_0ptw100.root");
  W[2] = new TFile ("TDAna/TDAna_W_1JETS_100ptw300.root");
  W[3] = new TFile ("TDAna/TDAna_W_2JETS_0ptw100.root");
  W[4] = new TFile ("TDAna/TDAna_W_2JETS_100ptw300.root");
  W[5] = new TFile ("TDAna/TDAna_W_3JETS_0ptw100.root"); 
  W[6] = new TFile ("TDAna/TDAna_W_3JETS_100ptw300.root");
  W[7] = new TFile ("TDAna/TDAna_W_4JETS_0ptw100.root");
  W[8] = new TFile ("TDAna/TDAna_W_4JETS_100ptw300.root");
  W[9] = new TFile ("TDAna/TDAna_W_5JETS_0ptw100.root");
  W[10] = new TFile ("TDAna/TDAna_W_5JETS_100ptw300.root");
  double Wxs[11] = { 45000., 9200., 250., 2500., 225., 590., 
		     100., 125., 40., 85., 40. };
  double NW[11] = { 88000., 40000., 100530., 99523., 105255., 79000., 
		    88258., 83038., 30796., 59022., 41865. };
  // Numbers from ascii files "decay".
  // ---------------------------------
  // MultiJet trigger
  double Wentries[11] = {0.,0.,0.,0.,5.,100.,10.,1.,4.,0.,5.};
//  // PixelJet trigger
//  double Wentries[11] = {0.,0.,0.,0.,7.,100.,9.,2.,5.,0.,6.};

  TFile * TTH = new TFile("TDAna/TDAna_TTH_120.root");
  double TTHxs = 0.667 ;
  double NTTH = 1634000.; // 1652000.; // 96000.;

  TFile * TT[5];
  TT[0] = new TFile("TDAna/TDAna_TT_0JETS.root");
  TT[1] = new TFile("TDAna/TDAna_TT_1JETS.root");
  TT[2] = new TFile("TDAna/TDAna_TT_2JETS.root");
  TT[3] = new TFile("TDAna/TDAna_TT_3JETS.root");
  TT[4] = new TFile("TDAna/TDAna_TT_4JETS.root");
  // double TTxs[5] = { 619., 176., 34.,  6., 1.5 };  // from web
  double TTxs[5] = { 434., 162., 43., 10., 1.9 };     // from note
  double NTT[5] = { 57900., 66000., 98159., 14768., 5352. };
  // Numbers which should be equal to those in ascii file "decay" 
  // (but we take them from histos to make sure we do not screw up)
  // --------------------------------------------------------------
  double TTentries[5]={0.};

  TFile * TTPYT;
  TTPYT = new TFile("TDAna/TDAna_TTBAR.root");
  double TTPYTxs = 650.9;
  double NTTPYT = 972000.;

  TFile * ST[3];
  ST[0] = new TFile("TDAna/TDAna_SINGLETOP_TQ_TQB_LHC_E.root");
  ST[1] = new TFile("TDAna/TDAna_SINGLETOP_TQ_TQB_LHC_MU.root");
  ST[2] = new TFile("TDAna/TDAna_SINGLETOP_TQ_TQB_LHC_TAU.root");
  double STxs[3] = { 27.43, 26.97, 28.71 };
  double NST[3] = { 92000, 94000, 94000 }; 

  TFile * ZNN[2];
  ZNN[0] = new TFile("TDAna/TDAna_ZNUNUJETS_120-170.root");
  ZNN[1] = new TFile("TDAna/TDAna_ZNUNUJETS_170-230.root");
  double ZNNxs[2] = { 51.47, 15.52 };
  double NZNN[2] = { 29897., 25600. };
  // Get Z entries from ascii files:
  // ---------------------------------
  // MultiJet trigger
  double Zentries[2] = { 0., 0. };
  //  // PixelJet trigger
  //  double Zentries[2] = { 0., 0. };

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

  TH1D * Histo_TOT[9];
  for ( int i=0; i<9; i++ ) {
    Histo_TOT[i] = new TH1D ( pippo1[i], pippo1[i], nbins, minx, maxx );
  }

  // Extract sum histograms with the right normalization and errors
  // --------------------------------------------------------------

  // QCD plus EXTRA QCD:
  // -------------------
  double totWQCD[16][nbins]={0.};
  double totQCD[nbins]={0.};
  double s2_totQCD[nbins]={0.};
  double totNQCD[16][nbins]={0.};
  double normQCD=0.;
  double normQCDnth[16]={0.};
  for ( int i=0; i<16; i++ ) {
    cout << "Processing QCD file #" << i << " ..." << endl;
    if ( i==2 ) i++; /////////////////////////// kludge
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
    normQCDnth[i]=QCDentries[i]*QCDxs[i]/NQCD[i]*Lumfactor;
    normQCD+=normQCDnth[i];
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
  if ( renormalize ) {
    double tw=0.;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      tw += totQCD[ibin-1];
    }
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      totQCD[ibin-1]=totQCD[ibin-1]*(normQCD/tw);
      s2_totQCD[ibin-1]=s2_totQCD[ibin-1]*pow(normQCD/tw,2);
    }
  }
  
  // QCD Multiplied:
  // ---------------
  double totQCDM[nbins]={0.};
  double s2_totQCDM[nbins]={0.};
  double normQCDM = 0.;
  double normQCDMnth[8] = {0.};
  for ( int i=0; i<8; i++ ) { 
    cout << "Processing QCD MULT file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(QCDMULT[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totQCDM[ibin-1]+=t*QCDxs[i]/NQCDM[i]*Lumfactor;
      s2_totQCDM[ibin-1]+=t*pow(QCDMxs[i]/NQCDM[i]*Lumfactor,2);
      QCDMentries[i]+=t;
    }
    normQCDMnth[i]=QCDMentries[i]*QCDMxs[i]/NQCDM[i]*Lumfactor;
    normQCDM+=normQCDMnth[i];
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
  double normTT=0.;
  double normTTnth[5]={0.};
  for ( int i=0; i<5; i++ ) { 
    cout << "Processing TT file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(TT[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totTT[ibin-1]+=t*TTxs[i]/NTT[i]*Lumfactor;
      s2_totTT[ibin-1]+=t*pow(TTxs[i]/NTT[i]*Lumfactor,2);
      TTentries[i]+=t;
    }
    normTTnth[i]=TTentries[i]*TTxs[i]/NTT[i]*Lumfactor;
    normTT+=normTTnth[i];
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
  double normW=0.;
  for ( int i=0; i<11; i++ ) {
    cout << "Processing W file #" << i << " ..." << endl;
    if ( i==6 ) i++;
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
    normW += Wentries[i]*Wxs[i]/NW[i]*Lumfactor;
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
  // IF WE WANT TO RENORMALIZE W+JETS MATRIX-WEIGHTED HISTOGRAMS TO NUMBER EXPECTATION:
  // ----------------------------------------------------------------------------------
  if ( renormalize ) {
    double tw=0.;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      tw += totW[ibin-1];
    }
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      totW[ibin-1]=totW[ibin-1]*(normW/tw);
      s2_totW[ibin-1]=s2_totW[ibin-1]*pow(normW/tw,2);
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
  double normZ=0.;
  for ( int i=0; i<2; i++ ) {
    cout << "Processing ZNN file #" << i << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(ZNN[i]->Get(pippo));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totZNN[ibin-1]+=t*ZNNxs[i]/NZNN[i]*Lumfactor;
      s2_totZNN[ibin-1]+=t*pow(ZNNxs[i]/NZNN[i]*Lumfactor,2);
    }
    normZ+=Zentries[i]*ZNNxs[i]/NZNN[i]*Lumfactor;
  }
  // IF WE WANT TO RENORMALIZE Z MATRIX-WEIGHTED HISTOGRAMS TO NUMBER EXPECTATION:
  // -----------------------------------------------------------------------------
  if ( renormalize ) {
    double tw=0.;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      tw += totZNN[ibin-1];
    }
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      totZNN[ibin-1]=totZNN[ibin-1]*(normZ/tw);
      s2_totZNN[ibin-1]=s2_totZNN[ibin-1]*pow(normZ/tw,2);
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
  }

  cout << "Filling nine bgrs" << QCDMentries[0] << " " << QCDMentries[1]  
       << " " << QCDMentries[2] << " " << QCDMentries[3] << " " << QCDMentries[4]  
       << " " << QCDMentries[5] << " " << QCDMentries[6] << " " << QCDMentries[7]  
    << endl;

  // Fill nine different backgrounds based on +- 1-sigma variations in QCD and/or TOP normalization
  // ----------------------------------------------------------------------------------------------
  // QCD: take integral, renormalize it to integral+-sqrt(integral);
  // TOP: take integral, renormalize it to integral+-sqrt(integral);
  double sigmaTTint=0.;
  double sigmaQCDMint=0.;
  for ( int isample=0; isample<8; isample++ ) {
    if ( isample<5 ) { 
      // TTentries[] is already the square of the error on the number of 
      // entries in the sample
      if ( TTentries[isample]<=10 ) sigmaTTint+=(1 + TMath::Sqrt(TTentries[isample] + 0.75));
      else sigmaTTint+=pow(normTTnth[isample],2)/TTentries[isample];
      
    }
   //    sigmaQCDint+=pow(normQCDnth[isample]*sigmaQCD[isample],2);
    if (QCDMentries[i]<=10) sigmaQCDMint+=(1 + TMath::Sqrt(QCDMentries[isample] + 0.75));
    else sigmaQCDMint+=pow(normQCDMnth[isample],2)/QCDMentries[isample];
  }
  sigmaTTint=sqrt(sigmaTTint);
  sigmaQCDMint=sqrt(sigmaQCDMint);

  cout << "Filling nine bgrs 2 " << endl;

  // Nine cases of fluctuations in QCD and/or TT:
  // --------------------------------------------
  double grandtot=0.;;
  double grandtote=0.;;
  double qcdtmp=0.;
  double s2_qcdtmp=0.;
  double tttmp=0.;
  double s2_tttmp=0.;
  for ( int ibin=1; ibin<=nbins; ibin++ ) {
    cout << "ibin" << ibin << endl;
    // case 0    
    grandtot = (totQCD[ibin-1]+totTTH[ibin-1]+totTT[ibin-1]+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_totQCD[ibin-1]+s2_totTTH[ibin-1]+s2_totTT[ibin-1]+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[0]->SetBinContent(ibin,grandtot);
    Histo_TOT[0]->SetBinError(ibin,grandtote);
    // example: grandtot[ifluctuation]=(totQCD[]*normQCDM+sqrt(normQCDM))
    cout << "ibin" << ibin << endl;
    // case 1: QCD+1 sigma
    // -------------------
    qcdtmp=totQCD[ibin-1]*(normQCDM+sigmaQCDMint)/normQCDM;
    s2_qcdtmp=s2_totQCD[ibin-1]*(normQCDM+sigmaQCDMint)/normQCDM;
    grandtot = (qcdtmp+totTTH[ibin-1]+totTT[ibin-1]+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_qcdtmp+s2_totTTH[ibin-1]+s2_totTT[ibin-1]+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[1]->SetBinContent(ibin,grandtot);
    Histo_TOT[1]->SetBinError(ibin,grandtote);
     cout << "ibin" << ibin << endl;
   // case 2: TT+1 sigma
    // ------------------
    tttmp=totTT[ibin-1]*(normTT+sigmaTTint)/normTT;
    s2_tttmp=s2_totTT[ibin-1]*(normTT+sigmaTTint)/normTT;
    grandtot = (totQCD[ibin-1]+totTTH[ibin-1]+tttmp+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_totQCD[ibin-1]+s2_totTTH[ibin-1]+s2_tttmp+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[2]->SetBinContent(ibin,grandtot);
    Histo_TOT[2]->SetBinError(ibin,grandtote);
     cout << "ibin" << ibin << endl;
   // case 3: QCD-1 sigma
    // -------------------
    qcdtmp=totQCD[ibin-1]*(normQCDM-sigmaQCDMint)/normQCDM;
    s2_qcdtmp=s2_totQCD[ibin-1]*(normQCDM-sigmaQCDMint)/normQCDM;
    grandtot = (qcdtmp+totTTH[ibin-1]+totTT[ibin-1]+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_qcdtmp+s2_totTTH[ibin-1]+s2_totTT[ibin-1]+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[3]->SetBinContent(ibin,grandtot);
    Histo_TOT[3]->SetBinError(ibin,grandtote);

    cout << "ibin" << ibin << endl;
    // case 4: TT-1 sigma
    // ------------------
    tttmp=totTT[ibin-1]*(normTT-sigmaTTint)/normTT;
    s2_tttmp=s2_totTT[ibin-1]*(normTT-sigmaTTint)/normTT;
    grandtot = (totQCD[ibin-1]+totTTH[ibin-1]+tttmp+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_totQCD[ibin-1]+s2_totTTH[ibin-1]+s2_tttmp+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[4]->SetBinContent(ibin,grandtot);
    Histo_TOT[4]->SetBinError(ibin,grandtote);

     cout << "ibin" << ibin << endl;
   // case 5: QCD,TT -1 sigma
    // ------------------
    qcdtmp=totQCD[ibin-1]*(normQCDM-sigmaQCDMint)/normQCDM;
    s2_qcdtmp=s2_totQCD[ibin-1]*(normQCDM-sigmaQCDMint)/normQCDM;
    tttmp=totTT[ibin-1]*(normTT-sigmaTTint)/normTT;
    s2_tttmp=s2_totTT[ibin-1]*(normTT-sigmaTTint)/normTT;
    grandtot = (qcdtmp+totTTH[ibin-1]+tttmp+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_qcdtmp+s2_totTTH[ibin-1]+s2_tttmp+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[5]->SetBinContent(ibin,grandtot);
    Histo_TOT[5]->SetBinError(ibin,grandtote);

    cout << "ibin" << ibin << endl;
    // case 6: QCD,TT +1 sigma
    // ------------------
    qcdtmp=totQCD[ibin-1]*(normQCDM+sigmaQCDMint)/normQCDM;
    s2_qcdtmp=s2_totQCD[ibin-1]*(normQCDM+sigmaQCDMint)/normQCDM;
    tttmp=totTT[ibin-1]*(normTT+sigmaTTint)/normTT;
    s2_tttmp=s2_totTT[ibin-1]*(normTT+sigmaTTint)/normTT;
    grandtot = (qcdtmp+totTTH[ibin-1]+tttmp+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_qcdtmp+s2_totTTH[ibin-1]+s2_tttmp+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[6]->SetBinContent(ibin,grandtot);
    Histo_TOT[6]->SetBinError(ibin,grandtote);

    cout << "ibin" << ibin << endl;
    // case 7: QCD+1,TT-1 sigma
    // ------------------
    qcdtmp=totQCD[ibin-1]*(normQCDM+sigmaQCDMint)/normQCDM;
    s2_qcdtmp=s2_totQCD[ibin-1]*(normQCDM+sigmaQCDMint)/normQCDM;
    tttmp=totTT[ibin-1]*(normTT-sigmaTTint)/normTT;
    s2_tttmp=s2_totTT[ibin-1]*(normTT-sigmaTTint)/normTT;
    grandtot = (qcdtmp+totTTH[ibin-1]+tttmp+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_qcdtmp+s2_totTTH[ibin-1]+s2_tttmp+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[7]->SetBinContent(ibin,grandtot);
    Histo_TOT[7]->SetBinError(ibin,grandtote);

    cout << "ibin" << ibin << endl;
    // case 8: QCD-1,TT+1 sigma
    // ------------------
    qcdtmp=totQCD[ibin-1]*(normQCDM-sigmaQCDMint)/normQCDM;
    s2_qcdtmp=s2_totQCD[ibin-1]*(normQCDM-sigmaQCDMint)/normQCDM;
    tttmp=totTT[ibin-1]*(normTT+sigmaTTint)/normTT;
    s2_tttmp=s2_totTT[ibin-1]*(normTT+sigmaTTint)/normTT;
    grandtot = (qcdtmp+totTTH[ibin-1]+tttmp+
		totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
    grandtote= sqrt(s2_qcdtmp+s2_totTTH[ibin-1]+s2_tttmp+
		    s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
    Histo_TOT[8]->SetBinContent(ibin,grandtot);
    Histo_TOT[8]->SetBinError(ibin,grandtote);
  }

  double S_int = Histo_TTH->Integral(imin,imax);
  double B_int[9] = {0.};
  for (int i = 0; i < 9; ++i ) {
    B_int[i] = Histo_TOT[i]->Integral(imin,imax);
    cout << "Selected events: " << B_int[i] << "B, " << S_int << "S. Significance = " << S_int/sqrt(B_int[i]) << endl;
  }
  cout << "TOP comparison: P= " << Histo_TT->KolmogorovTest(Histo_TTPYT) << endl;
  cout << "QCD comparison: P= " << Histo_QCD->KolmogorovTest(Histo_QCDM) << endl;
  
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
  if ( dops ) QCDcomp->Print("./ps/"+pippo+"_QCDcomp.ps");

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
  if ( dops ) TTcomp->Print("./ps/"+pippo+"_TTcomp.ps");

  TCanvas * b = new TCanvas ("b", "Kinematics comparison", 700, 700 );
  b->Divide(1,2);  
  b->cd(1);
  b->GetPad(1)->SetLogy();
  Histo_TOT[0]->SetMarkerStyle(20);
  Histo_TOT[0]->SetMarkerSize(0.4);
  Histo_TOT[0]->SetMinimum(1.);
  Histo_TOT[0]->DrawCopy("PE");
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
  Histo_TOT[0]->DrawCopy("PESAME");
  b->cd(2);
  Histo_TOT[0]->SetMarkerStyle(20);
  Histo_TOT[0]->SetMarkerSize(0.4);
  Histo_TOT[0]->SetMinimum(1.);
  Histo_TOT[0]->DrawCopy("PE");
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
  Histo_TOT[0]->DrawCopy("PESAME");
  if ( dops ) b->Print("./ps/"+pippo+".ps");

   // Histograms to fit:
   // ------------------
   TFile * xmia = new TFile("xmia.root","RECREATE");
   xmia->cd();
   for ( int i = 0; i < 9; ++i ) Histo_TOT[i]->Write();
   Histo_TTH->Write();
   xmia->Close();



}
