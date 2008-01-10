void Smooth_1 (TString v, TString sel) 
{

  const int nbins = 50;

  const int nvars=1;
  TString var[nvars] = v;
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
  
  const int nqcdsamples=8;
  TFile * QCD[nqcdsamples];     
  QCD[0] = new TFile("./root/TDAna_QCD30-50_tk3.root");
  QCD[1] = new TFile("./root/TDAna_QCD50-80_tk3.root");
  QCD[2] = new TFile("./root/TDAna_QCD80-120_tk3.root");
  QCD[3] = new TFile("./root/TDAna_QCD120-170_tk3.root");
  QCD[4] = new TFile("./root/TDAna_QCD170-230_tk3.root");
  QCD[5] = new TFile("./root/TDAna_QCD230-300_tk3.root");
  QCD[6] = new TFile("./root/TDAna_QCD300-380_tk3.root");
  QCD[7] = new TFile("./root/TDAna_QCD380incl_tk3.root");
  double QCDxs[nqcdsamples] = { 155929000., 20938850., 2949713., 499656., 100995.,  23855., 6391., 2821.};
  double NQCD[nqcdsamples] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};


  const int nwsamples=11;
  TFile * W[nwsamples];
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
  double Wxs[nwsamples] = { 45000., 9200., 250., 2500., 225., 590., 100., 125., 40., 85., 40. };
  double NW[nwsamples] = { 88000., 40000., 100530., 99523., 105255., 79000., 
			   88258., 83038., 30796., 59022., 41865. };

  TFile * TTH = new TFile("./root/TDAna_ttH_120_tk3.root");
  double TTHxs = 0.667 ;
  double NTTH = 62000.; // 1652000.; // 62000.;

  const int nttsamples=5;
  TFile * TT[nttsamples];
  TT[0] = new TFile("./root/TDAna_TT0_tk3.root");
  TT[1] = new TFile("./root/TDAna_TT1_tk3.root");
  TT[2] = new TFile("./root/TDAna_TT2_tk3.root");
  TT[3] = new TFile("./root/TDAna_TT3_tk3.root");
  TT[4] = new TFile("./root/TDAna_TT4_tk3.root");
  // double TTxs[5] = { 619., 176., 34.,  6., 1.5 };  // from web
  double TTxs[nttsamples] = { 434., 162., 43., 10., 1.9 };     // from note
  double NTT[nttsamples] = { 57900., 66000., 98159., 14768., 5304. };

  double Lumfactor = 100000; // 100/fb of luminosity assumed in histograms

  TH1D * Histo_TOT[nvars];
  TH1D * Histo_TTH[nvars];
  TH1D * Histo_TOTS[nvars];
  TH1D * Histo_TTHS[nvars];
  for ( int i=0; i<nvars; i++ ) {
    cout << i << endl;
    TH1D * H = dynamic_cast<TH1D*>(TTH->Get(pippo[i]));
    double minx=H->GetBinLowEdge(1);
    double maxx=nbins*H->GetBinWidth(1);
    Histo_TOT[i] = new TH1D ( pippotot[i]," ", nbins, minx, maxx );
    Histo_TTH[i] = new TH1D ( pippotth[i]," ", nbins, minx, maxx );
    Histo_TOTS[i] = new TH1D ( pippototS[i]," ", nbins, minx, maxx );
    Histo_TTHS[i] = new TH1D ( pippotthS[i]," ", nbins, minx, maxx );
  }

  cout << "Starting loop on variables needing smoothing" << endl;

  // Loop on variables 
  // -----------------
  for ( int ivar=0; ivar<nvars; ivar++ ) {

    // Extract sum histograms with the right normalization and errors
    // --------------------------------------------------------------
    double totWW[nwsamples][nbins]={0.};
    double totW[nbins]={0.};
    double s2_totW[nbins]={0.};
    double totNW[nwsamples][nbins]={0.};
    for ( int i=0; i<nwsamples; i++ ) {
      cout << "Processing W file #" << i << " ..." << endl;
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
    for ( int i=0; i<nwsamples; i++ ) {
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
	totW[ibin-1]+=totWW[i][ibin-1];
	if ( totNW[i][ibin-1]>0 ) {
	  s2_totW[ibin-1]+=pow(totWW[i][ibin-1],2)/totNW[i][ibin-1];
	}
      }
    }
    
    double totWQCD[nqcdsamples][nbins]={0.};
    double totQCD[nbins]={0.};
    double s2_totQCD[nbins]={0.};
    double totNQCD[nqcdsamples][nbins]={0.};
    for ( int i=0; i<nqcdsamples; i++ ) {
      cout << "Processing QCD file #" << i << " ..." << endl;
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
    for ( int i=0; i<nqcdsamples; i++ ) {
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
	totQCD[ibin-1]+=totWQCD[i][ibin-1];
	if ( totNQCD[i][ibin-1]>0 ) {
	  s2_totQCD[ibin-1]+=pow(totWQCD[i][ibin-1],2)/totNQCD[i][ibin-1];
	}
      }
    }
    
    double totTT[nbins]={0.};
    double s2_totTT[nbins]={0.};
    for ( int i=0; i<nttsamples; i++ ) {
      cout << "Processing TT file #" << i << " ..." << endl;
      TH1D * Histo = dynamic_cast<TH1D*>(TT[i]->Get(pippo[ivar]));
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
	double t=Histo->GetBinContent(ibin);
	totTT[ibin-1]+=t*TTxs[i]/NTT[i]*Lumfactor;
	s2_totTT[ibin-1]+=t*pow(TTxs[i]/NTT[i]*Lumfactor,2);
      }
    }
    double totTTH[nbins]={0.};
    double s2_totTTH[nbins]={0.};
    cout << "Processing TTH file " << " ..." << endl;
    TH1D * Histo = dynamic_cast<TH1D*>(TTH->Get(pippo[ivar]));
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      double t=Histo->GetBinContent(ibin);
      totTTH[ibin-1]+=t*TTHxs/NTTH*Lumfactor;
      s2_totTTH[ibin-1]+=t*pow(TTHxs/NTTH*Lumfactor,2);
    }
    
    // OK now fill total histograms
    // ----------------------------
    double total_sig=0.;
    double total_bgr=0.;
    double grandtot[nbins]={0.};
    double grandtote[nbins]={0.};
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      Histo_TTH[ivar]->SetBinContent(ibin,totTTH[ibin-1]);
      Histo_TTH[ivar]->SetBinError(ibin,sqrt(s2_totTTH[ibin-1]));
      Histo_TTHS[ivar]->SetBinContent(ibin,totTTH[ibin-1]);
      Histo_TTHS[ivar]->SetBinError(ibin,sqrt(s2_totTTH[ibin-1]));
      grandtot[ibin-1] = totQCD[ibin-1]+totTT[ibin-1]+totW[ibin-1];
      grandtote[ibin-1]= sqrt(s2_totQCD[ibin-1]+s2_totTT[ibin-1]+s2_totW[ibin-1]);
      Histo_TOT[ivar]->SetBinContent(ibin,grandtot[ibin-1]);
      Histo_TOT[ivar]->SetBinError(ibin,grandtote[ibin-1]);
      Histo_TOTS[ivar]->SetBinContent(ibin,grandtot[ibin-1]);
      Histo_TOTS[ivar]->SetBinError(ibin,grandtote[ibin-1]);
      total_sig+=totTTH[ibin-1];
      total_bgr+=grandtot[ibin-1];
    }

    // Pre-smoothing algorithm that levels down spikes
    // -----------------------------------------------
    double mean;
    double errmean;
    double delta_bgr=0.;
    for ( int ibin=0; ibin<nbins-2; ibin++ ) {
      double delta=fabs(grandtot[ibin]-grandtot[ibin+1]);
      if ( delta>3*grandtote[ibin] && grandtot[ibin]>0 ) {  // the signal is a bin way off from the previous
	// Avoid averaging with bins with large error
	if ( delta>3*grandtote[ibin+2] ) {
	  mean=(grandtot[ibin]+grandtot[ibin+2])/2.;
	  errmean=sqrt(pow(grandtote[ibin],2)+pow(grandtote[ibin+2],2))/2.;
	  cout << ivar << " " << ibin << " " << mean <<"+-" << errmean << " " << grandtot[ibin+1];
	  if ( mean>0 && fabs(grandtot[ibin+1]-mean)>3.*errmean ) {
	    delta_bgr += grandtot[ibin+1]-mean;
	    grandtot[ibin+1]=mean;
	    grandtote[ibin+1]=errmean;
	    cout << ": doing it " << mean << " " << errmean << " " << delta_bgr << endl;
	    cout << grandtot[ibin] << " " << grandtot[ibin+1] << " " << grandtot[ibin+2] << endl;
	  } else cout << ": not doing it." << endl;
	}
      }
    }
    // fix bin 0 if needed
    // -------------------
    mean=(grandtot[1]+grandtot[2])/2;
    errmean=sqrt(pow(grandtote[1],2)+pow(grandtote[2],2))/2.;
    if ( mean>0 & fabs(grandtot[0]-mean)>3.*errmean ) {
      delta_bgr += grandtot[0]-mean;
      grandtot[0]=mean;
      grandtote[0]=errmean;
    } 
    // fix last bin if needed
    // ----------------------
    mean=(grandtot[nbins-3]+grandtot[nbins-2])/2;
    errmean=sqrt(pow(grandtote[nbins-3],2)+pow(grandtote[nbins-2],2))/2;
    if ( mean>0 & fabs(grandtot[nbins-1]-mean)>3.*errmean ) {
      delta_bgr += grandtot[nbins-1]-mean;
      grandtot[nbins-1]=mean;
      grandtote[nbins-1]=errmean;
    } 
    // Renormalize bgr histogram
    // -------------------------
    cout << ivar << " " << ": delta=" << delta_bgr << endl;
    if ( delta_bgr>0 && total_bgr>delta_bgr ) {
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
	Histo_TOTS[ivar]->SetBinContent(ibin,grandtot[ibin-1]*total_bgr/(total_bgr-delta_bgr));
	Histo_TOTS[ivar]->SetBinError(ibin,grandtote[ibin-1]*total_bgr/(total_bgr-delta_bgr));
      } 
    } 
    
    Histo_TOTS[ivar]->Smooth(1);
    Histo_TTHS[ivar]->Smooth(1);
    double sumtot=0.;
    double sumtth=0.;
    double sumtotS=0.;
    double sumtthS=0.;
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

  } // end of ivar loop

  TString fname;
  fname="./root/"+v+sel+".root";

  TFile * Smoothed = new TFile(fname,"RECREATE");
  Smoothed->cd();

  TCanvas * b1 = new TCanvas ("b1", "Kinematics comparison", 600, 600 );
  for ( int ivar=0; ivar<nvars; ivar++ ) {
    b1->cd(ivar+1);
    Histo_TOTS[ivar]->SetMinimum(0.);
    Histo_TOTS[ivar]->Draw();
    Histo_TTHS[ivar]->SetLineColor(kBlue);
    Histo_TTHS[ivar]->Draw("PESAME");
  }
  b1->Print("./ps/"+v+"_smooth_svsb_1.ps");

  TCanvas * b3 = new TCanvas ("b3", "Kinematics comparison", 600, 600 );
  for ( int ivar=0; ivar<nvars; ivar++ ) {
    b3->cd(ivar+1);
    Histo_TOT[ivar]->SetMinimum(0.);
    Histo_TOT[ivar]->SetLineColor(kRed);
    Histo_TOT[ivar]->Draw("PE");    
    Histo_TOTS[ivar]->Draw("PESAME");
    Histo_TOT[ivar]->Write();
    Histo_TTH[ivar]->Write();
    Histo_TOTS[ivar]->Write();
    Histo_TTHS[ivar]->Write();
  }
  b3->Print("./ps/"+v+"_smooth_check_1.ps");

  Smoothed->Close();

}
