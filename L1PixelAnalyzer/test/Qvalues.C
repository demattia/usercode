void Qvalues (TString sel ) 
{

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

  const int nvars=17;
  TString var[nvars] = { "C8", "M8", "C6", "M6", "MEt", "MEtSig", "CorrSumEt", "GoodHt", 
		      "M45bestall", "Chi2mass", "Chi2extall", "Mbbnoh", "DPbbnoh", 
		      "SumHED4", "SumHED6", "DP12", "MEtDP2" };
  const int nbins=50;
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

    ////////////////////////////////////////////////////////

    // Now compute max quality factor for this distribution
    // -----------------------------------------------------
    TH1D * S_  = new TH1D ( "S",  "S "+var[ivar],  nbins, minx, maxx );     // Signal distribution
    TH1D * B_  = new TH1D ( "B",  "B "+var[ivar],  nbins, minx, maxx );     // Total background distribution
    TH1D * ES_ = new TH1D ( "ES", "ES "+var[ivar], nbins, minx, maxx );  // Efficiency of signal 
    TH1D * EB_ = new TH1D ( "EB", "EB "+var[ivar], nbins, minx, maxx);   // Efficiency of background
    TH1D * R_  = new TH1D ( "R",  "R "+var[ivar],  nbins, minx, maxx );     // Esig/Ebgr distribution
    TH1D * Q_  = new TH1D ( "Q",  "Q "+var[ivar],  nbins, minx, maxx );     // Q-value
    
    // First of all, fill normalized distributions
    // -------------------------------------------
    double IS = Histo_TTHS->Integral();
    double IB = Histo_TOTS->Integral();
    double maxS=0;
    double maxB=0;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double cont_tth = Histo_TTHS->GetBinContent(ibin);
      double erro_tth = Histo_TTHS->GetBinError(ibin);
      S_->SetBinContent(ibin,cont_tth/IS);
      S_->SetBinError(ibin,erro_tth/IS);
      if ( maxS<cont_tth/IS ) maxS=cont_tth/IS;
      double cont_tot = Histo_TOTS->GetBinContent(ibin);
      double erro_tot = Histo_TOTS->GetBinError(ibin);
      B_->SetBinContent(ibin,cont_tot/IB);
      B_->SetBinError(ibin,erro_tot/IB);
      if ( maxB<cont_tot/IB ) maxB=cont_tot/IB;
    }

    ////////////////////////////////////////////
    // Efficiency and gains for x>xcut 
    // -------------------------------
    double es[nbins];
    double s_es[nbins];
    double eb[nbins];
    double s_eb[nbins];
    for ( int ibin=1; ibin<=nbins; ibin++ ) {

      // Get signal efficiency curve
      // ---------------------------
      double fs=0.;    // signal failing cut
      double s2_fs=0.;
      double ps=0.;    // signal passing cut
      double s2_ps=0.;
      for ( int jbin=1; jbin<=nbins; jbin++ ) {
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
      for ( int jbin=1; jbin<=nbins; jbin++ ) {
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
    double maxQm1s=0;
    double s_maxQm1s=0;
    double x_maxQ=0;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double R;
      double s_R;
      double Q;
      double s_Q;
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
	if ( maxQm1s<Q-s_Q ) {
	  maxQm1s=Q;
	  s_maxQm1s=s_Q;
	  x_maxQ=(double)ibin*(maxx-minx)/(double)nbins;
	} 
      } 
    }
    double maxX=maxS;
    if (maxB>maxX) maxX=maxB;
    cout << var[ivar] << " : maximum quality factor " << maxQm1s << "+-" << s_maxQm1s 
	 << " at " << var[ivar] << " > " << x_maxQ << endl; 
    
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
    b->Print("./ps/"+pippo[ivar]+"_opt_lt.ps");

    ////////////////////////////////////////////
    // Efficiency and gains for x<xcut 
    // -------------------------------
    for ( int ibin=1; ibin<=nbins; ibin++ ) {

      // Get signal efficiency curve
      // ---------------------------
      double fs=0.;    // signal failing cut
      double s2_fs=0.;
      double ps=0.;    // signal passing cut
      double s2_ps=0.;
      for ( int jbin=1; jbin<=nbins; jbin++ ) {
	if ( jbin>=ibin ) { 
	  fs    += S_->GetBinContent(jbin);
	  s2_fs += pow(S_->GetBinError(jbin),2);
	}
	if ( jbin<ibin ) {
	  ps    += S_->GetBinContent(jbin);
	  s2_ps += pow(S_->GetBinError(jbin),2);
	}
      }
      if ( fs+ps>0. ) {
	es[ibin-1]=ps/(fs+ps);
	s_es[ibin-1]=sqrt(ps*ps*s2_fs+fs*fs*s2_ps)/pow(fs+ps,2);
      } else {
	es[ibin-1]=0.;
	s_es[ibin-1]=0.;
      }
      ES_->SetBinContent(ibin,es[ibin-1]);
      ES_->SetBinError(ibin,s_es[ibin-1]);

      // Get background efficiency curve
      // -------------------------------
      double fb=0.;    // background failing cut
      double s2_fb=0.;
      double pb=0.;    // background passing cut
      double s2_pb=0.;
      for ( int jbin=1; jbin<=nbins; jbin++ ) {
	if ( jbin>=ibin ) { 
	  fb    += B_->GetBinContent(jbin);
	  s2_fb += pow(B_->GetBinError(jbin),2);
	}
	if ( jbin<ibin ) {
	  pb    += B_->GetBinContent(jbin);
	  s2_pb += pow(B_->GetBinError(jbin),2);
	}
      }
      if ( fb+pb>0. ) {
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
    maxQm1s=0.;
    s_maxQm1s=0.;
    x_maxQ=0.;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double R;
      double s_R;
      double Q;
      double s_Q;
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
	if ( maxQm1s<Q-s_Q ) {
	  maxQm1s=Q;
	  s_maxQm1s=s_Q;
	  x_maxQ=(double)ibin*(maxx-minx)/(double)nbins;
	} 
      } 
//       else {
// 	R = es[ibin-1]/0.1;
// 	s_R = sqrt(pow(s_es[ibin-1]/0.1,2)+
// 		   pow(es[ibin-1]*s_eb[ibin-1],2)/pow(0.1,4));
// 	Q = es[ibin-1]*es[ibin-1]/0.1;
// 	s_Q = sqrt(pow(2*es[ibin-1]/0.1*s_es[ibin-1],2)+pow(es[ibin-1]/0.1,4)*pow(s_eb[ibin-1],2));
//       }
    }
    double maxX=maxS;
    if (maxB>maxX) maxX=maxB;
    cout << var[ivar] << " : maximum quality factor " << maxQm1s << "+-" << s_maxQm1s 
	 << " at " << var[ivar] << " <= " << x_maxQ << endl; 
    
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
    b->Print("./ps/"+pippo[ivar]+"_opt_gt.ps");

    delete S_;
    delete B_;
    delete ES_;
    delete EB_;
    delete R_;
    delete Q_;

    delete Histo_TTHS;
    delete Histo_TOTS;

  } // end ivar loop

}
