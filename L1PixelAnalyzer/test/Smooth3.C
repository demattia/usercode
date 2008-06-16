void Smooth (TString sel, double frac_true_err, bool dops, bool doroot ) {

  // -C8 : maximum quality factor 1.1314+-0.00503415 at LR > 0.4
  // -M8 : maximum quality factor 1.03056+-0.00499927 at LR > 1.4
  // -C6 : maximum quality factor 1.18157+-0.00673922 at LR > 0.8
  // -M6 : maximum quality factor 1.01414+-0.00240137 at LR > 1.2
  // MEt : maximum quality factor 1.00301+-0.00088719 at LR > 0.8
  // -MEtSig : maximum quality factor 1.01873+-0.00367832 at LR > 1.4
  // -CorrSumEt : maximum quality factor 1.02956+-0.00487776 at LR > 1.4
  // GoodHt : maximum quality factor 1.05092+-0.00244139 at LR > 0.4
  // M45bestall : maximum quality factor 1.03104+-0.00370672 at LR > 1
  // -Chi2mass : maximum quality factor 1.13193+-0.00467282 at LR > 0.4
  // Chi2extall : maximum quality factor 1.07473+-0.00716133 at LR > 1.4
  // Mbbnoh : maximum quality factor 2.17059+-0.0321631 at LR > 0.8
  // DPbbnoh : maximum quality factor 1.30261+-0.0105512 at LR > 0.8
  // SumHED4 : maximum quality factor 1.08131+-0.00809236 at LR > 1.4
  // SumHED6 : maximum quality factor 1.09023+-0.0160105 at LR > 1.8
  // DP12 : maximum quality factor 2.29376+-0.0270215 at LR > 0.4
  // MEtDPM : maximum quality factor 1+-0.1 at LR > 0
  // -MEtDP1 : maximum quality factor 1.06697+-0.00357291 at LR > 0.4
  // MEtDP2 : maximum quality factor 1+-0.1 at LR > 0
  // -M45best : maximum quality factor 1.02626+-0.00341181 at LR > 1.2
  // -M_others : maximum quality factor 1.11418+-0.00859052 at LR > 1.4
  // -Et6 : maximum quality factor 1.02596+-0.00268856 at LR > 1
  
  const int nbins = 50;
  
  const int nvars=25;
  TString var[nvars] = { "C8", "M8", "C6", "M6", "MEt", "MEtSig", "CorrSumEt", "GoodHt", 
			 "Hbestcomb", "Chi2mass", "Mbbnoh", "DPbbnoh", 
			 "SumHED4", "SumHED6", "MEtDPM", "MEtDP1", "MEtDP2",
			 "M_others", "Et6", "Scprod", "Thdeta", "M5", "M3best", "Mwbest",
                         "TTMS1"};

  TString pippo[nvars];
  TString pippotot[nvars];
  TString pippotth[nvars];
  TString pippototS[nvars];
  TString pippotthS[nvars];
  for ( int i=0; i<nvars; i++ ) { 
    pippo[i] = var[i]+sel; 
    pippotot[i]=var[i]+sel+"_bgr";
    pippotth[i]=var[i]+sel+"_sig";
    pippototS[i]=var[i]+sel+"_bgrS";
    pippotthS[i]=var[i]+sel+"_sigS";
  }
  
  TFile * QCD[16];     
  QCD[0] = new TFile("./rootnew/TDAna_QCD_30-50.root");
  QCD[1] = new TFile("./rootnew/TDAna_QCD_50-80.root");
  QCD[2] = new TFile("./rootnew/TDAna_QCD_80-120.root");
  QCD[3] = new TFile("./rootnew/TDAna_QCD_120-170.root");
  QCD[4] = new TFile("./rootnew/TDAna_QCD_170-230.root");
  QCD[5] = new TFile("./rootnew/TDAna_QCD_230-300.root");
  QCD[6] = new TFile("./rootnew/TDAna_QCD_300-380.root");
  QCD[7] = new TFile("./rootnew/TDAna_QCD_380-incl.root");
  QCD[8] = new TFile("./rootnew/TDAna_EXTRA_QCD_30-50.root");
  QCD[9] = new TFile("./rootnew/TDAna_EXTRA_QCD_50-80.root");
  QCD[10] = new TFile("./rootnew/TDAna_EXTRA_QCD_80-120.root");
  QCD[11] = new TFile("./rootnew/TDAna_EXTRA_QCD_120-170.root");
  QCD[12] = new TFile("./rootnew/TDAna_EXTRA_QCD_170-230.root");
  QCD[13] = new TFile("./rootnew/TDAna_EXTRA_QCD_230-300.root");
  QCD[14] = new TFile("./rootnew/TDAna_EXTRA_QCD_300-380.root");
  QCD[15] = new TFile("./rootnew/TDAna_EXTRA_QCD_380-incl.root");
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
  QCDMULT[0] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_30-50.root");
  QCDMULT[1] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_50-80.root");
  QCDMULT[2] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_80-120.root");
  QCDMULT[3] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_120-170.root");
  QCDMULT[4] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_170-230.root");
  QCDMULT[5] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_230-300.root");
  QCDMULT[6] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_300-380.root");
  //QCDMULT[7] = new TFile("./rootnew/TDAna_MULTI_EXTRA_QCD_380-incl.root");
  double QCDMxs[8] =  { 155929000., 20938850., 2949713., 499656., 
			100995.,  23855., 6391., 2821. };
  double NQCDM[8] = { 11229443., 26347818., 24971508.,  29514603., 
		      40608575., 32268604., 37943909., 33232000. };

  TFile * W[11];
  W[0] = new TFile ("./rootnew/TDAna_W_0JETS.root");
  W[1] = new TFile ("./rootnew/TDAna_W_1JETS_0ptw100.root");
  W[2] = new TFile ("./rootnew/TDAna_W_1JETS_100ptw300.root");
  W[3] = new TFile ("./rootnew/TDAna_W_2JETS_0ptw100.root");
  W[4] = new TFile ("./rootnew/TDAna_W_2JETS_100ptw300.root");
  W[5] = new TFile ("./rootnew/TDAna_W_3JETS_0ptw100.root"); 
  W[6] = new TFile ("./rootnew/TDAna_W_3JETS_100ptw300.root");
  W[7] = new TFile ("./rootnew/TDAna_W_4JETS_0ptw100.root");
  W[8] = new TFile ("./rootnew/TDAna_W_4JETS_100ptw300.root");
  W[9] = new TFile ("./rootnew/TDAna_W_5JETS_0ptw100.root");
  W[10] = new TFile ("./rootnew/TDAna_W_5JETS_100ptw300.root");
  double Wxs[11] = { 45000., 9200., 250., 2500., 225., 590., 
		     100., 125., 40., 85., 40. };
  double NW[11] = { 88000., 40000., 100530., 99523., 105255., 79000., 
		    88258., 83038., 30796., 59022., 41865. };

  TFile * TTH = new TFile("./rootnew/TDAna_TTH_120.root");
  double TTHxs = 0.667 ;
  double NTTH = 1634000.; // 1652000.; // 96000.;

  TFile * TT[5];
  TT[0] = new TFile("./rootnew/TDAna_TT_0JETS.root");
  TT[1] = new TFile("./rootnew/TDAna_TT_1JETS.root");
  TT[2] = new TFile("./rootnew/TDAna_TT_2JETS.root");
  TT[3] = new TFile("./rootnew/TDAna_TT_3JETS.root");
  TT[4] = new TFile("./rootnew/TDAna_TT_4JETS.root");
  // double TTxs[5] = { 619., 176., 34.,  6., 1.5 };  // from web
  double TTxs[5] = { 434., 162., 43., 10., 1.9 };     // from note
  double NTT[5] = { 57900., 66000., 98159., 14768., 5352. };

  TFile * TTPYT;
  TTPYT = new TFile("./rootnew/TDAna_TTBAR.root");
  double TTPYTxs = 650.9;
  double NTTPYT = 972000.;

  TFile * ST[3];
  ST[0] = new TFile("./rootnew/TDAna_SINGLETOP_TQ_TQB_LHC_E.root");
  //ST[1] = new TFile("./rootnew/TDAna_SINGLETOP_TQ_TQB_LHC_MU.root");
  ST[2] = new TFile("./rootnew/TDAna_SINGLETOP_TQ_TQB_LHC_TAU.root");
  double STxs[3] = { 27.43, 26.97, 28.71 };
  double NST[3] = { 92000., 94000., 94000. };

  TFile * ZNN[2];
  ZNN[0] = new TFile("./rootnew/TDAna_ZNUNUJETS_120-170.root");
  ZNN[1] = new TFile("./rootnew/TDAna_ZNUNUJETS_170-230.root");
  double ZNNxs[2] = { 51.47, 15.52 };
  double NZNN[2] = { 29897., 25600. };

  double Lumfactor = 100000; // 100/fb of luminosity assumed in histograms

  //   TH1D * Histo_QCD;
  //   TH1D * Histo_QCDM;
  //   TH1D * Histo_TT;
  //   TH1D * Histo_TTPYT;
  //   TH1D * Histo_W;
  //   TH1D * Histo_ST;
  //   TH1D * Histo_ZNN;

  TH1D * Histo_TOT[nvars];
  TH1D * Histo_TTH[nvars];
  TH1D * Histo_TOTS[nvars];
  TH1D * Histo_TTHS[nvars];
  for ( int i=0; i<nvars; i++ ) {
    cout << i << endl;
    TH1D * H = dynamic_cast<TH1D*>(TTH->Get(pippo[i]));
    double minx=H->GetBinLowEdge(1);
    double maxx=nbins*H->GetBinWidth(1);
    Histo_TOT[i] = new TH1D ( pippotot[i],pippotot[i], nbins, minx, maxx );
    Histo_TTH[i] = new TH1D ( pippotth[i],pippotth[i], nbins, minx, maxx );
    Histo_TOTS[i] = new TH1D ( pippototS[i],pippototS[i], nbins, minx, maxx );
    Histo_TTHS[i] = new TH1D ( pippotthS[i],pippotthS[i], nbins, minx, maxx );
    //     Histo_QCD = new TH1D ( pippo+"_QCD", "", nbins, minx, maxx );
    //     Histo_QCDM = new TH1D ( pippo+"_QCDM", "", nbins, minx, maxx );
    //     Histo_TT = new TH1D ( pippo+"_TT", "", nbins, minx, maxx );
    //     Histo_TTPYT = new TH1D ( pippo+"_TTPYT", "", nbins, minx, maxx );
    //     Histo_W = new TH1D ( pippo+"_W", "", nbins, minx, maxx );
    //     Histo_ST = new TH1D ( pippo+"_ST", "", nbins, minx, maxx );
    //     Histo_ZNN = new TH1D ( pippo+"_ZNN", "", nbins, minx, maxx );
  }
  
  cout << "Starting loop on variables needing smoothing" << endl;

  // Loop on variables 
  // -----------------
  for ( int ivar=0; ivar<nvars; ivar++ ) {

    //if ( ivar==2 ) { // kludge start
    
    // Extract sum histograms with the right normalization and errors
    // --------------------------------------------------------------

    // QCD plus EXTRA QCD:
    // -------------------
    double totWQCD[16][nbins]={0.};
    double totQCD[nbins]={0.};
    double s2_totQCD[nbins]={0.};
    double totNQCD[16][nbins]={0.};
    for ( int i=0; i<16; i++ ) {
      //if ( i==2 ) i++;                             //      <<<<<<<<<<<<--- NB (AND BELOW AGAIN!)
      // cout << "Processing QCD file #" << i << " ..." << endl;
      TH1D * Histo = dynamic_cast<TH1D*>(QCD[i]->Get(pippo[ivar]));
      TH1D * HistoW = dynamic_cast<TH1D*>(QCD[i]->Get(pippo[ivar]+"W"));
      // For QCD, we need also total entries in histograms to add a
      // Poisson fluke contribution to total errors from matrix:
      // ----------------------------------------------------------
      TH1D * HistoN = dynamic_cast<TH1D*>(QCD[i]->Get(pippo[ivar]+"N"));    
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
      if ( i==2 ) i++;        //   <<<<<<<<<<<<<< ----------------------- NB
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
        totQCD[ibin-1]+=totWQCD[i][ibin-1];
        if ( totNQCD[i][ibin-1]>0 ) {
          s2_totQCD[ibin-1]+=pow(totWQCD[i][ibin-1],2)/totNQCD[i][ibin-1];
        }
      }
    }

    // QCD Multiplied:
    // ---------------
    double totQCDM[nbins]={0.};
    double s2_totQCDM[nbins]={0.};
    for ( int i=0; i<7; i++ ) { ////////////// NB  
      // cout << "Processing QCD MULT file #" << i << " ..." << endl;
      TH1D * Histo = dynamic_cast<TH1D*>(QCDMULT[i]->Get(pippo[ivar]));
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
    // cout << "Processing TTH file " << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(TTH->Get(pippo[ivar]));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totTTH[ibin-1]+=t*TTHxs/NTTH*Lumfactor;
      s2_totTTH[ibin-1]+=t*pow(TTHxs/NTTH*Lumfactor,2);
    }

    // TT AlpGen:
    // ----------
    double totTT[nbins]={0.};
    double s2_totTT[nbins]={0.};
    for ( int i=0; i<5; i++ ) { //////////////////////// NNBB
      // cout << "Processing TT file #" << i << " ..." << endl;
      TH1D * Histo = dynamic_cast<TH1D*>(TT[i]->Get(pippo[ivar]));
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
    // cout << "Processing TTPYT file " << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(TTPYT->Get(pippo[ivar]));
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
      // if ( i==5 ) i++; //////////////// NNBB
      // cout << "Processing W file #" << i << " ..." << endl;
      TH1D * Histo = dynamic_cast<TH1D*>(W[i]->Get(pippo[ivar]));
      TH1D * HistoW = dynamic_cast<TH1D*>(W[i]->Get(pippo[ivar]+"W"));
      // For W, we need also total entries in histograms to add a
      // Poisson fluke contribution to total errors from matrix:
      // ----------------------------------------------------------
      TH1D * HistoN = dynamic_cast<TH1D*>(W[i]->Get(pippo[ivar]+"N"));    
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

    // ST :
    // ----
    double totST[nbins]={0.};
    double s2_totST[nbins]={0.};
    for ( int i=0; i<3; i++ ) {
      cout << "Processing ST file #" << i << " ..." << endl;
      if ( i==1 ) i++;
      TH1D * Histo = dynamic_cast<TH1D*>(ST[i]->Get(pippo[ivar]));
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
      // cout << "Processing ZNN file #" << i << " ..." << endl;
      TH1D * Histo = dynamic_cast<TH1D*>(ZNN[i]->Get(pippo[ivar]));
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
        double t=Histo->GetBinContent(ibin);
        totZNN[ibin-1]+=t*ZNNxs[i]/NZNN[i]*Lumfactor;
        s2_totZNN[ibin-1]+=t*pow(ZNNxs[i]/NZNN[i]*Lumfactor,2);
      }
    }

    // OK now fill total histograms
    // ----------------------------
    double total_sig=0.;
    double total_bgr=0.;
    double grandtot[nbins]={0.};
    double grandtote[nbins]={0.};
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      Histo_TTH[ivar]->SetBinContent(ibin,totTTH[ibin-1]);
      Histo_TTHS[ivar]->SetBinContent(ibin,totTTH[ibin-1]);
      Histo_TTH[ivar]->SetBinError(ibin,sqrt(s2_totTTH[ibin-1]));
      Histo_TTHS[ivar]->SetBinError(ibin,sqrt(s2_totTTH[ibin-1]));
      grandtot[ibin-1] = (totQCD[ibin-1]+totTTH[ibin-1]+totTT[ibin-1]+
                          totW[ibin-1]+totST[ibin-1]+totZNN[ibin-1]);
      grandtote[ibin-1]= sqrt(s2_totQCD[ibin-1]+s2_totTTH[ibin-1]+s2_totTT[ibin-1]+
                              s2_totW[ibin-1]+s2_totST[ibin-1]+s2_totZNN[ibin-1]);
      total_sig+=totTTH[ibin-1];
      total_bgr+=grandtot[ibin-1];
      Histo_TOT[ivar]->SetBinContent(ibin,grandtot[ibin-1]);
      Histo_TOT[ivar]->SetBinError(ibin,grandtote[ibin-1]);
      //       Histo_TOTS[ivar]->SetBinContent(ibin,grandtot);
      //       Histo_TOTS[ivar]->SetBinError(ibin,grandtote);
    }


    //////////////////////////////
    
    // Preliminary: we set to zero bins with less than 0.01% of the integral
    // ---------------------------------------------------------------------
    for ( int ibin=0; ibin<nbins; ibin++ ) {
      if ( grandtot[ibin]/total_bgr<0.0001 ) grandtot[ibin]=0.;
    }
    // Pre-pre-smoothing: kill bins with abnormally large errors, 
    // replace with average of neighbors
    // ----------------------------------------------------------
    double averror=0.;
    double nav=0.;
    for ( int ibin=0; ibin<nbins; ibin++ ) {
      if ( grandtot[ibin]>0 && grandtote[ibin]>0 ) {
	averror+=grandtote[ibin]/grandtot[ibin];
	nav++;
      }
    }
    if ( nav>0 ) {
      averror=averror/nav;
      double err;
      double nlterr[100]={0.};
      double max67pc;
      int i67=-1;
      for ( int i=0; i<100 && i67==-1; i++ ) {
	err=2*averror/100.*(double)i;
	for ( int ibin=0; ibin<nbins; ibin++ ) {
	  if ( grandtot[ibin]>0 && grandtote[ibin]>0 ) {
	    if ( grandtote[ibin]/grandtot[ibin]<err ) nlterr[i]++;
	  }
	}
	nlterr[i]=nlterr[i]/nav;
	if ( nlterr[i]>frac_true_err ) {
	  max67pc=err;
	  i67=i;
	}
      }
    }

    // try three-point line fits: y=a+bx
    // ---------------------------------
    for ( int ibin=0; ibin<nbins-2; ibin++ ) { // consider triplets 0/1/2 ... n-3/n-2/n-1
      if ( grandtote[ibin]*grandtote[ibin+1]*grandtote[ibin+2]>0 ) {
	double sumx=0.;
	double sumx2=0.;
	double sumy=0.;
	double sumxy=0.;
	double sums=0.;
	for ( int in=0; in<3; in++ ) {
	  double s2 = pow(grandtote[ibin+in],2);
	  double x = (in-1)*(maxx-minx)/nbins; // x=0 for central bin
	  sumx += x/s2;
	  sumx2+= pow(x,2)/s2;
	  sumy += grandtot[ibin+in]/s2;
	  sumxy+= x*grandtot[ibin+in]/s2;
	  sums += 1/s2;
	}
	double delta = sums*sumx2-sumx*sumx;
	double a = (sumx2*sumy-sumx*sumxy)/delta;
	double b = (sums*sumxy-sumx*sumy)/delta;
	double s2_a = sumx2/delta;
	double s2_b = sums/delta;
	// compare difference between central bin and fit with fit error
	if ( a - grandtot[ibin+1] > 2*sqrt(s2_a) || grandtote[ibin+1]>max67pc*grandtot[ibin+1] ) {
	  cout << var[ivar] << ": smoothing bin #" << ibin << " since fit = " 
	       << a << "+-" << sqrt(s2_a) << " while content = " << grandtot[ibin+1] 
	       << "+-" << grandtote[ibin+1] << endl;
	  cout << a - grandtot[ibin+1] << " " << 2*sqrt(s2_a)  << " " <<  grandtote[ibin+1] 
	       << " " << max67pc*grandtot[ibin+1] << endl;  
	  cout << "contents: " 
	       << "  i=-1: " << grandtot[ibin] << "+-" << grandtote[ibin] 
	       << "  i=0 : " << grandtot[ibin+1] << "+-" << grandtote[ibin+1] 
	       << "  i=1 : " << grandtot[ibin+2] << "+-" << grandtote[ibin+2] 
	       << "  a= " << a << "+-" << sqrt(s2_a) 
	       << "  b= " << b << "+-" << sqrt(s2_b) << endl;
	  grandtot[ibin+1]=a;
	  grandtote[ibin+1]=sqrt(s2_a);
	}
      }
    }
    //////////////////////////////
    

    //     // Preliminary: we set to zero bins with less than 0.01% of the integral
    //     // ---------------------------------------------------------------------
    //     for ( int ibin=0; ibin<nbins; ibin++ ) {
    //       if ( grandtot[ibin]/total_bgr<0.0001 ) grandtot[ibin]=0.;
    //     }
    //     // Pre-pre-smoothing: kill bins with abnormally large errors, 
    //     // replace with average of neighbors
    //     // ----------------------------------------------------------
    //     double averror=0.;
    //     double nav=0.;
    //     for ( int ibin=0; ibin<nbins; ibin++ ) {
    //       if ( grandtot[ibin]>0 && grandtote[ibin]>0 ) {
    // 	averror+=grandtote[ibin]/grandtot[ibin];
    // 	nav++;
    //       }
    //     }
    //     if ( nav>0 ) {
    //       averror=averror/nav;
    //       double err;
    //       double nlterr[100]={0.};
    //       double max67pc;
    //       int i67=-1;
    //       for ( int i=0; i<100 && i67==-1; i++ ) {
    // 	err=2*averror/100.*(double)i;
    // 	for ( int ibin=0; ibin<nbins; ibin++ ) {
    // 	  if ( grandtot[ibin]>0 && grandtote[ibin]>0 ) {
    // 	    if ( grandtote[ibin]/grandtot[ibin]<err ) nlterr[i]++;
    // 	  }
    // 	}
    // 	nlterr[i]=nlterr[i]/nav;
    // 	if ( nlterr[i]>frac_true_err ) {
    // 	  max67pc=err;
    // 	  i67=i;
    // 	}
    //       }
    //       cout << "var#" << ivar<< " averror=" << averror << " nav=" 
    // 	   << nav << " max67pc=" << max67pc << endl; 

    //       // Ok now we know which bins to disregard. Let us correct them:
    //       // ------------------------------------------------------------
    //       int inext;
    //       int iprev;
    //       for ( int ibin=0; ibin<nbins; ibin++ ) {

    // 	cout << "Bin #" << ibin << ":" << grandtot[ibin] << "+-" << grandtote[ibin] << endl;

    // 	if ( grandtot[ibin]>0. && grandtote[ibin]>0. ) {
    // 	  if ( grandtote[ibin]/grandtot[ibin]>2*max67pc ) { // consider this bin for smoothing
    // 	    cout << "      will be considered for smoothing." << endl;
    // 	    // Starting bin:
    // 	    // -------------
    // 	    if ( ibin==0 ) {
    // 	      inext=0;
    // 	      // look at the three next ones
    // 	      // ---------------------------
    // 	      for ( int jbin=ibin+1; jbin<ibin+4 && inext==0; jbin++ ) {
    // 		if ( grandtot[jbin]>0. && grandtote[jbin]>0. ) {
    // 		  if ( grandtote[jbin]/grandtot[jbin]<2*max67pc ) inext=jbin;
    // 		}
    // 	      }
    // 	      if ( inext==0 ) inext=1;
    // 	      iprev=inext; // This will set bin zero value at bin inext value	 
    // 	      // Other bins:
    // 	      // -----------
    // 	    } else if ( ibin<nbins-1 ) {
    // 	      iprev=ibin-1;
    // 	      inext=0;
    // 	      // look at the next ones
    // 	      // --------------------
    // 	      for ( int jbin=ibin+1; jbin<ibin+4 && jbin<nbins && inext==0; jbin++ ) {
    // 		if ( grandtot[jbin]>0. && grandtote[jbin]>0. ) {
    // 		  if ( grandtote[jbin]/grandtot[jbin]<2*max67pc ) inext=jbin;
    // 		}
    // 	      }
    // 	      if ( inext==0 ) inext=iprev;
    // 	    } else if ( ibin==nbins-1 ) { 
    // 	      iprev=ibin-1;
    // 	      inext=ibin-1;
    // 	    } 
    // 	    // Compute new value of bin, taking into account the prescription that
    // 	    // when averaging, an average of x and zero is given value zero.
    // 	    // Also, do not replace bins which are consistent with their neighbors 
    // 	    // within the errors in the average of the latter
    // 	    // ------------------------------------------------------------------
    // 	    if ( grandtot[iprev]>0. && grandtot[inext]>0. ) {
    // 	      //	      double av = (grandtot[iprev]+grandtot[inext])/2;
    // 	      //	      double s_av = sqrt((pow(grandtote[iprev],2)+pow(grandtote[inext],2))/4.);
    // 	      double av = (grandtot[iprev]/pow(grandtote[iprev],2)+
    // 			   grandtot[inext]/pow(grandtote[inext],2))/
    // 		(1/pow(grandtote[iprev],2)+1/pow(grandtote[inext],2));
    // 	      double s_av = sqrt(1/(1/pow(grandtote[iprev],2)+1/pow(grandtote[inext],2)));
    // 	      // NB we insert a kludge to avoid considering downward flukes 
    // 	      // to bins too low to have a small relative error and which, for
    // 	      // our construction of a mixed template, are mostly real
    // 	      if ( fabs(grandtot[ibin]-av)>3*s_av && grandtot[ibin]>grandtot[iprev]) {//////////nb
    // 		//	      if ( fabs(grandtot[ibin]-av)>3*s_av && 
    // 		// (grandtot[ibin]-grandtot[iprev])/(grandtot[inext]-grandtot[ibin])<0 ) {
    // 		grandtot[ibin]=av;
    // 		grandtote[ibin]=s_av; 
    // 	      }
    // 	    } else {
    // 	      grandtot[ibin]=0.;
    // 	      grandtote[ibin]=sqrt((pow(grandtote[iprev],2)+pow(grandtote[inext],2))/4.);
    // 	    }
    // 	  }
    // 	} 
    //       }
    //     }


    
    double total_bgr_s=0.;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      total_bgr_s+=grandtot[ibin-1];
      // cout << "filling Histo_TOTS " << ibin << endl;
      Histo_TOTS[ivar]->SetBinContent(ibin,grandtot[ibin-1]);
      Histo_TOTS[ivar]->SetBinError(ibin,grandtote[ibin-1]);
    } 
    
    // Now do the regular smoothing
    // ----------------------------
    //
    // Histo_TOTS[ivar]->Smooth(1);
    // Histo_TTHS[ivar]->Smooth(1);
    double sumtot=0.;
    double sumtth=0.;
    double sumtotS=0.;
    double sumtthS=0.; 
//     if ( ivar==11 ) Histo_TOTS[ivar]->SetBinContent(1,0.); //  Kludge to avoid the bin at zero
//     if ( ivar==11 ) Histo_TTHS[ivar]->SetBinContent(1,0.); //  Kludge to avoid the bin at zero
//     if ( ivar==11 ) Histo_TOTS[ivar]->SetBinError(1,0.);
//     if ( ivar==11 ) Histo_TTHS[ivar]->SetBinError(1,0.);
//     if ( ivar==11 ) Histo_TOT[ivar]->SetBinContent(1,0.); //  Kludge to avoid the bin at zero
//     if ( ivar==11 ) Histo_TTH[ivar]->SetBinContent(1,0.); //  Kludge to avoid the bin at zero
//     if ( ivar==11 ) Histo_TOT[ivar]->SetBinError(1,0.);
//     if ( ivar==11 ) Histo_TTH[ivar]->SetBinError(1,0.);
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      sumtotS+=Histo_TOTS[ivar]->GetBinContent(ibin);
      sumtthS+=Histo_TTHS[ivar]->GetBinContent(ibin);
      sumtot+=Histo_TOT[ivar]->GetBinContent(ibin);
      sumtth+=Histo_TTH[ivar]->GetBinContent(ibin);
    }
    if ( sumtotS>0 ) Histo_TOTS[ivar]->Scale(1./sumtotS);
    if ( sumtthS>0 ) Histo_TTHS[ivar]->Scale(1./sumtthS);
    if ( sumtot>0 ) Histo_TOT[ivar]->Scale(1./sumtot);
    if ( sumtth>0 ) Histo_TTH[ivar]->Scale(1./sumtth);

    //} // kludge end
    
  } // end of ivar loop

  
  cout << "Done, now plotting and writing histos." << endl;
  
  // File where smoothed histos are stored
  // -------------------------------------
  TString fname;
  fname="functionfile" + sel + ".root";
  if ( doroot ) TFile * Smoothed = new TFile(fname,"RECREATE");
  if ( doroot ) Smoothed->cd();
  
  TCanvas * b1 = new TCanvas ("b1", "Kinematics comparison", 600, 600 );
  b1->Divide(2,4);
  for ( int ivar=0; ivar<8; ivar++ ) {
    b1->cd(ivar+1);
    double a = Histo_TOTS[ivar]->GetMaximum(1000000000.);
    if (a<Histo_TTHS[ivar]->GetMaximum(100000000.) ) a=Histo_TTHS[ivar]->GetMaximum(100000000.);
    Histo_TOTS[ivar]->SetMaximum(1.2*a);
    Histo_TOTS[ivar]->SetMinimum(0.);
    Histo_TOTS[ivar]->Draw();
    Histo_TTHS[ivar]->SetLineColor(kBlue);
    Histo_TTHS[ivar]->Draw("PESAME");
  }
  if ( dops ) b1->Print("./ps/Smooth_svsb_"+sel+"1.ps");
  
  TCanvas * b2 = new TCanvas ("b2", "Kinematics comparison", 600, 600 );
  b2->Divide(2,4);
  for ( int ivar=8; ivar<16; ivar++ ) {
    b2->cd(ivar-7);
    double a = Histo_TOTS[ivar]->GetMaximum(1000000000.);
    if (a<Histo_TTHS[ivar]->GetMaximum(100000000.) ) a=Histo_TTHS[ivar]->GetMaximum(100000000.);
    Histo_TOTS[ivar]->SetMaximum(1.2*a);
    Histo_TOTS[ivar]->SetMinimum(0.);
    Histo_TOTS[ivar]->Draw();
    Histo_TTHS[ivar]->SetLineColor(kBlue);
    Histo_TTHS[ivar]->Draw("PESAME");
  }
  if ( dops ) b2->Print("./ps/Smooth_svsb_"+sel+"2.ps");
  
  TCanvas * b3 = new TCanvas ("b3", "Kinematics comparison", 600, 600 );
  b3->Divide(2,4);
  for ( int ivar=16; ivar<24; ivar++ ) {
    b3->cd(ivar-15);
    double a = Histo_TOTS[ivar]->GetMaximum(1000000000.);
    if ( a<Histo_TTHS[ivar]->GetMaximum(100000000.) ) a=Histo_TTHS[ivar]->GetMaximum(100000000.);
    Histo_TOTS[ivar]->SetMaximum(1.2*a);
    Histo_TOTS[ivar]->SetMinimum(0.);
    Histo_TOTS[ivar]->Draw();
    Histo_TTHS[ivar]->SetLineColor(kBlue);
    Histo_TTHS[ivar]->Draw("PESAME");
  }
  if ( dops) b3->Print("./ps/Smooth_svsb_"+sel+"3.ps");
  
  TCanvas * b4 = new TCanvas ("b4", "Kinematics comparison", 600, 600 );
  b4->Divide(2,4);
  for ( int ivar=24; ivar<nvars; ivar++ ) {
    b4->cd(ivar-23);
    double a = Histo_TOTS[ivar]->GetMaximum(1000000000.);
    if (a<Histo_TTHS[ivar]->GetMaximum(100000000.) ) a=Histo_TTHS[ivar]->GetMaximum(100000000.);
    Histo_TOTS[ivar]->SetMaximum(1.2*a);
    Histo_TOTS[ivar]->SetMinimum(0.);
    Histo_TOTS[ivar]->Draw();
    Histo_TTHS[ivar]->SetLineColor(kBlue);
    Histo_TTHS[ivar]->Draw("PESAME");
  }
  if ( dops) b4->Print("./ps/Smooth_svsb_"+sel+"4.ps");

  TCanvas * c1 = new TCanvas ("c1", "Kinematics comparison", 600, 600 );
  c1->Divide(2,4);
  for ( int ivar=0; ivar<8; ivar++ ) {
    c1->cd(ivar+1);
    double a = Histo_TOT[ivar]->GetMaximum(1000000000.);
    if ( a<Histo_TOTS[ivar]->GetMaximum(100000000.) ) a=Histo_TOTS[ivar]->GetMaximum(100000000.);
    Histo_TOT[ivar]->SetMaximum(1.2*a);
    Histo_TOT[ivar]->SetMinimum(0.);
    Histo_TOT[ivar]->SetLineColor(kRed);
    Histo_TOT[ivar]->Draw("PE");    
    Histo_TOTS[ivar]->Draw("PESAME");
    if ( doroot ) {
      Histo_TOT[ivar]->Write();
      Histo_TTH[ivar]->Write();
      Histo_TOTS[ivar]->Write();
      Histo_TTHS[ivar]->Write();
    }
  }
  if ( dops )  c1->Print("./ps/Smooth_check_"+sel+"1.ps");
  
  TCanvas * c2 = new TCanvas ("c2", "Kinematics comparison", 600, 600 );
  c2->Divide(2,4);
  for ( int ivar=8; ivar<16; ivar++ ) {
    c2->cd(ivar-7);
    double a = Histo_TOT[ivar]->GetMaximum(1000000000.);
    if ( a<Histo_TOTS[ivar]->GetMaximum(100000000.) ) a=Histo_TOTS[ivar]->GetMaximum(100000000.);
    Histo_TOT[ivar]->SetMaximum(1.2*a);
    Histo_TOT[ivar]->SetMinimum(0.);
    Histo_TOT[ivar]->SetLineColor(kRed);
    Histo_TOT[ivar]->Draw("PE");    
    Histo_TOTS[ivar]->Draw("PESAME");
    if ( doroot ) {
      Histo_TOT[ivar]->Write();
      Histo_TTH[ivar]->Write();
      Histo_TOTS[ivar]->Write();
      Histo_TTHS[ivar]->Write();
    }
    if ( dops ) c2->Print("./ps/Smooth_check_"+sel+"2.ps");
  }
  TCanvas * c3 = new TCanvas ("c3", "Kinematics comparison", 600, 600 );
  c3->Divide(2,4);
  for ( int ivar=16; ivar<24; ivar++ ) {
    c3->cd(ivar-15);
    double a = Histo_TOT[ivar]->GetMaximum(1000000000.);
    if ( a<Histo_TOTS[ivar]->GetMaximum(100000000.) ) a=Histo_TOTS[ivar]->GetMaximum(100000000.);
    Histo_TOT[ivar]->SetMaximum(1.2*a);
    Histo_TOT[ivar]->SetMinimum(0.);
    Histo_TOT[ivar]->SetLineColor(kRed);
    Histo_TOT[ivar]->Draw("PE");    
    Histo_TOTS[ivar]->Draw("PESAME");
    if ( doroot ) {
      Histo_TOT[ivar]->Write();
      Histo_TTH[ivar]->Write();
      Histo_TOTS[ivar]->Write();
      Histo_TTHS[ivar]->Write();
    }
  }
  if ( dops ) c3->Print("./ps/Smooth_check_"+sel+"3.ps");
  
  TCanvas * c4 = new TCanvas ("c4", "Kinematics comparison", 600, 600 );
  c4->Divide(2,4);
  for ( int ivar=24; ivar<nvars; ivar++ ) {
    c4->cd(ivar-23);
    double a = Histo_TOT[ivar]->GetMaximum(1000000000.);
    if (a<Histo_TOTS[ivar]->GetMaximum(100000000.) ) a=Histo_TOTS[ivar]->GetMaximum(100000000.);
    Histo_TOT[ivar]->SetMaximum(1.2*a);
    Histo_TOT[ivar]->SetMinimum(0.);
    Histo_TOT[ivar]->SetLineColor(kRed);
    Histo_TOT[ivar]->Draw("PE");    
    Histo_TOTS[ivar]->Draw("PESAME");
    if ( doroot ) {
      Histo_TOT[ivar]->Write();
      Histo_TTH[ivar]->Write();
      Histo_TOTS[ivar]->Write();
      Histo_TTHS[ivar]->Write();
    }
    }
  if ( dops ) c4->Print("./ps/Smooth_check_"+sel+"4.ps");
  
  if ( doroot ) Smoothed->Close();
  
}
