void Qvalues_lik (TString sel ) 
{

  // Al livello SS si ha:
  // C6 : maximum quality factor 1.32383+-0.0147123 at LR > 1.6
  // Thdeta : maximum quality factor 1.14013+-0.00640858 at LR > 1.2
  // Et6 : maximum quality factor 1.10808+-0.0165823 at LR > 1.8
  // Mbbnoh : maximum quality factor 1.10478+-0.096053 at LR > 3.6
  // GoodHt : maximum quality factor 1.04881+-0.00237432 at LR > 0.4
  // MEtSig : maximum quality factor 1.03554+-0.00867608 at LR > 1.6
  // M8 : maximum quality factor 1.02226+-0.00460643 at LR > 1.4
  // MEtDP2 : maximum quality factor 1.0201+-0.00191511 at LR > 0.4
  
  // C8 : maximum quality factor 1.12895+-0.00503226 at LR > 0.4
  // M6 : maximum quality factor 1.00924+-0.00221231 at LR > 1.2
  // MEt : maximum quality factor 1.00703+-0.00241349 at LR > 1
  // CorrSumEt : maximum quality factor 1.02076+-0.00401397 at LR > 1.4
  // Hbestcomb : maximum quality factor 1.00574+-0.00318069 at LR > 1.4
  // Chi2mass : maximum quality factor 1+-0.1 at LR > 0
  // DPbbnoh : maximum quality factor 1.01703+-0.00562878 at LR > 1.4
  // SumHED4 : maximum quality factor 1.1178+-0.0165147 at LR > 1.8
  // SumHED6 : maximum quality factor 1.03873+-0.00649612 at LR > 1.4
  // MEtDPM : maximum quality factor 1+-0.1 at LR > 0
  // MEtDP1 : maximum quality factor 1.01565+-0.00412023 at LR > 1.4
  // M_others : maximum quality factor 1+-0.1 at LR > 0
  // Scprod : maximum quality factor 1.01074+-0.00200097 at LR > 1.2
  // M5 : maximum quality factor 1.01335+-0.00219015 at LR > 0.8
  
  // In the attempt at defining a yardstick with which to decide which variables
  // are the most discriminating for separating signal from sum of backgrounds,
  // we could consider:
  // - the max KS difference between tthS and totS histograms
  // - the Q-value of a cut "x<x_cut" or "x>x_cut" selection, with varying x_cut
  // - the Q-value of a cut "x_min<x<x_max" with varying x_min, x_max (this would
  //   be a good idea for variables such as higgs mass, which peaks more in a particular
  //   region for the signal)
  // - the KS between two distributions: the LR distr. obtained from background-distributed
  //   events and the LR distr. obtained from signal-distributed events.
  //   
  // Let's first try this latter approach which seems the best motivated (the KS of the LR
  // distributions is the correct measure to quantify how much a 2-comp fit will separate the
  // two).

  const int nvars=24;
  TString var[nvars] = { "C8", "M8", "C6", "M6", "MEt", "MEtSig", "CorrSumEt", "GoodHt", 
			 "Hbestcomb", "Chi2mass", "Mbbnoh", "DPbbnoh", 
			 "SumHED4", "SumHED6", "MEtDPM", "MEtDP1", "MEtDP2",
			 "M_others", "Et6", "Scprod", "Thdeta", "M5", "M3best", "Mwbest" }

  const int nbins=50;
  const int nbinsl=20;
  TString pippo[nvars];
  TString pippotot[nvars];
  TString pippotth[nvars];
  TString pippototS[nvars];
  TString pippotthS[nvars];
  for ( int ivar=0; ivar<nvars; ivar++ ) { 
    pippo[ivar] = var[ivar]+sel; 
    pippotot[ivar]=var[ivar]+sel+"_bgr";
    pippotth[ivar]=var[ivar]+sel+"_sig";
    pippototS[ivar]=var[ivar]+sel+"_bgrS";
    pippotthS[ivar]=var[ivar]+sel+"_sigS";
  }

  TString fname;
  fname="functionfile"+sel+".root";
  TFile * Smoothed = new TFile(fname);
  Smoothed->cd();

  TCanvas * b = new TCanvas ("b", "Cut optimization", 500, 500 );
  b->Divide(2,2);
  for ( int ivar=0; ivar<nvars; ivar++ ) {

    TH1D * Histo_TTHS = dynamic_cast<TH1D*>(Smoothed->Get(pippotthS[ivar]));
    TH1D * Histo_TOTS = dynamic_cast<TH1D*>(Smoothed->Get(pippototS[ivar]));
    double minx=Histo_TTHS->GetBinLowEdge(1);
    double maxx=nbins*Histo_TTHS->GetBinWidth(1)+minx;
    
    ////////////////////////////////////////////////////////
    // Check of Likelihood ratio per each variable
    // -------------------------------------------
    TH1D * LS = new TH1D ( "LS", "LS", nbinsl, -2., 2. );
    TH1D * LB = new TH1D ( "LB", "LB", nbinsl, -2., 2. );
    LS->Sumw2();
    LB->Sumw2();
    for ( int i=0; i<10000; i++ ) {
      // Fish at random from bgr shape
      double x = Histo_TOTS->GetRandom();
      int bx;
      double r;
      bx= (int)((x-minx)/(maxx-minx)*nbins)+1;
      if ( bx>=1 && bx<=nbins ) {
	double fs = Histo_TTHS->GetBinContent(bx);
	double fb = Histo_TOTS->GetBinContent(bx);
	if (fb>0 && fs>0 ) {
	  r = fs/fb;
	} else {
	  if ( fb==0. ) cout << "problem with getrandom" << endl;
	  if ( fs==0. ) r = 0.01; 
	}
      }
      double rel_lik;
      if ( r>=0 ) {
	rel_lik = log(r);
      } else {
	rel_lik=-1.99;	    
      }
      if ( rel_lik<-1.99 ) rel_lik=-1.99;
      if ( rel_lik>1.99 ) rel_lik=1.99;
      LB->Fill(rel_lik);
    }
    for ( int i=0; i<10000; i++ ) {
      // Fish at random from sig shape
      double x = Histo_TTHS->GetRandom();
      int bx;
      double r;
      bx= (int)((x-minx)/(maxx-minx)*nbins)+1;
      if ( bx>=1 && bx<=nbins ) {
	double fs = Histo_TTHS->GetBinContent(bx);
	double fb = Histo_TOTS->GetBinContent(bx);
	if (fb>0 && fs>0 ) {
	  r = fs/fb;
	} else {
	  if ( fs==0 ) cout << "problem with getrandom" << endl;
	  if ( fb==0 ) r = 100.; 
	}
      }
      double rel_lik;
      if ( r>=0 ) {
	rel_lik = log(r);
      } else {
	rel_lik=-1.99;	    
      }
      if ( rel_lik<-1.99 ) rel_lik=-1.99;
      if ( rel_lik>1.99 ) rel_lik=1.99;
      LS->Fill(rel_lik);
    }

    ////////////////////////////////////////////////////////

    // Now compute max quality factor for this distribution
    // -----------------------------------------------------
    TH1D * S_  = new TH1D ( "S",  "S "+var[ivar],  nbinsl, -2., 2. );     // Signal distribution
    TH1D * B_  = new TH1D ( "B",  "B "+var[ivar],  nbinsl, -2., 2. );     // Total background distribution
    TH1D * ES_ = new TH1D ( "ES", "ES "+var[ivar], nbinsl, -2., 2. );  // Efficiency of signal 
    TH1D * EB_ = new TH1D ( "EB", "EB "+var[ivar], nbinsl, -2., 2. );   // Efficiency of background
    TH1D * R_  = new TH1D ( "R",  "R "+var[ivar],  nbinsl, -2., 2. );     // Esig/Ebgr distribution
    TH1D * Q_  = new TH1D ( "Q",  "Q "+var[ivar],  nbinsl, -2., 2. );     // Q-value
    
    // First of all, fill normalized distributions
    // -------------------------------------------
    double IS = LS->Integral();
    double IB = LB->Integral();
    double maxS=0;
    double maxB=0;
    for ( int ibin=1; ibin<=nbinsl; ibin++ ) {
      double cont_tth = LS->GetBinContent(ibin);
      double erro_tth = LS->GetBinError(ibin);
      S_->SetBinContent(ibin,cont_tth/IS);
      S_->SetBinError(ibin,erro_tth/IS);
      if ( maxS<cont_tth/IS ) maxS=cont_tth/IS;
      double cont_tot = LB->GetBinContent(ibin);
      double erro_tot = LB->GetBinError(ibin);
      B_->SetBinContent(ibin,cont_tot/IB);
      B_->SetBinError(ibin,erro_tot/IB);
      if ( maxB<cont_tot/IB ) maxB=cont_tot/IB;
    }

    ////////////////////////////////////////////
    // Efficiency and gains for x>xcut 
    // -------------------------------
    double es[nbinsl];
    double s_es[nbinsl];
    double eb[nbinsl];
    double s_eb[nbinsl];
    for ( int ibin=1; ibin<=nbinsl; ibin++ ) {

      // Get signal efficiency curve
      // ---------------------------
      double fs=0.;    // signal failing cut
      double s2_fs=0.;
      double ps=0.;    // signal passing cut
      double s2_ps=0.;
      for ( int jbin=1; jbin<=nbinsl; jbin++ ) {
	if ( jbin<ibin ) { 
	  fs    += S_->GetBinContent(jbin);
	  s2_fs += pow(S_->GetBinError(jbin),2);
	}
	if ( jbin>=ibin ) {
	  ps    += S_->GetBinContent(jbin);
	  s2_ps += pow(S_->GetBinError(jbin),2);
	}
      }
      if ( fs+ps>0. ) {
	es[ibin-1]=ps/(fs+ps);
	s_es[ibin-1]=sqrt(ps*ps*s2_fs+fs*fs*s2_ps)/pow(fs+ps,2);
      } else {
	es[ibin-1]=0;
	s_es[ibin-1]=0;
      }
      ES_->SetBinContent(ibin,es[ibin-1]);
      ES_->SetBinError(ibin,s_es[ibin-1]);

      // Get background efficiency curve
      // -------------------------------
      double fb=0.;    // background failing cut
      double s2_fb=0.;
      double pb=0.;    // background passing cut
      double s2_pb=0.;
      for ( int jbin=1; jbin<=nbinsl; jbin++ ) {
	if ( jbin<ibin ) { 
	  fb    += B_->GetBinContent(jbin);
	  s2_fb += pow(B_->GetBinError(jbin),2);
	}
	if ( jbin>=ibin ) {
	  pb    += B_->GetBinContent(jbin);
	  s2_pb += pow(B_->GetBinError(jbin),2);
	}
      }
      if ( fb+pb>0 ) {
	eb[ibin-1]=pb/(fb+pb);
	s_eb[ibin-1]=sqrt(pb*pb*s2_fb+fb*fb*s2_pb)/pow(fb+pb,2);
      } else {
	eb[ibin-1]=0.;
	s_eb[ibin-1]=0.;
      }
      EB_->SetBinContent(ibin,eb[ibin-1]);
      EB_->SetBinError(ibin,s_eb[ibin-1]);
    }

    // Now compute efficiency ratio and quality factor
    // -----------------------------------------------
    double maxQm1s=1.;
    double s_maxQm1s=0.1;
    double x_maxQ=0;
    for ( int ibin=1; ibin<=nbinsl; ibin++ ) {
      double R;
      double s_R;
      double Q=1.; 
      double s_Q=0.1;
      if ( eb[ibin-1]>0. ) {
	R = es[ibin-1]/eb[ibin-1];
	s_R = sqrt(s_es[ibin-1]*s_es[ibin-1]/eb[ibin-1]/eb[ibin-1]+
		   es[ibin-1]*es[ibin-1]*s_eb[ibin-1]*s_eb[ibin-1]/pow(eb[ibin-1],4));
	Q = es[ibin-1]*es[ibin-1]/eb[ibin-1];
	s_Q = sqrt(pow(2*es[ibin-1]/eb[ibin-1]*s_es[ibin-1],2)+
		   pow(es[ibin-1]/eb[ibin-1],4)*pow(s_eb[ibin-1],2));
	R_->SetBinContent(ibin,R);
	R_->SetBinError(ibin,s_R);
	Q_->SetBinContent(ibin,Q);
	Q_->SetBinError(ibin,s_Q);
	if ( maxQm1s<Q-s_Q && (maxQm1s-1.)/s_maxQm1s<(Q-1.)/s_Q ) {
	  maxQm1s=Q;
	  s_maxQm1s=s_Q;
	  x_maxQ=(double)ibin*4/(double)nbinsl;
	} 
      } 
    }
    double maxX=maxS;
    if (maxB>maxX) maxX=maxB;
    cout << var[ivar] << " : maximum quality factor " << maxQm1s << "+-" << s_maxQm1s 
	 << " at LR > " << x_maxQ << endl; 
    
    b->cd(1);
    S_->SetMaximum(1.2*maxX);
    S_->SetLineColor(kRed);
    S_->SetMarkerColor(kRed);
    S_->DrawCopy("PE");
    B_->SetLineColor(kBlue);
    B_->SetMarkerColor(kBlue);
    B_->DrawCopy("PESAME");
    b->cd(2);
    ES_->SetLineColor(kRed);
    ES_->SetMarkerColor(kRed);
    ES_->DrawCopy("PE");
    EB_->SetLineColor(kBlue);
    EB_->SetMarkerColor(kBlue);
    EB_->DrawCopy("PESAME");
    b->cd(3);
    R_->DrawCopy("PE");
    b->cd(4);
    Q_->DrawCopy("PE");
    b->Print("./ps/"+pippo[ivar]+"_opt_L.ps");

    delete S_;
    delete B_;
    delete ES_;
    delete EB_;
    delete R_;
    delete Q_;

    delete Histo_TTHS;
    delete Histo_TOTS;
    delete LS;
    delete LB;

  } // end ivar loop

}
