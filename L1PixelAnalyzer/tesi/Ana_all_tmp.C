void Ana_all (TString var, TString sel ) 
{

  const int nbins=50;
  TString pippo = var+sel;

  TFile * QCD[8];     
  QCD[0] = new TFile("./TDAna_QCD30-50_tk3.root");
  QCD[1] = new TFile("./TDAna_QCD50-80_tk3.root");
  QCD[2] = new TFile("./TDAna_QCD80-120_tk3.root");
  QCD[3] = new TFile("./TDAna_QCD120-170_tk3.root");
  QCD[4] = new TFile("./TDAna_QCD170-230_tk3.root");
  QCD[5] = new TFile("./TDAna_QCD230-300_tk3.root");
  //QCD[6] = new TFile("./TDAna_QCD300-380_tk3.root");
  //QCD[7] = new TFile("./TDAna_QCD380incl_tk3.root");
  QCD[6] = new TFile("./smallrun/TDAna_QCD300-380_tk3.root");
  QCD[7] = new TFile("./smallrun/TDAna_QCD380incl_tk3.root");
  double QCDxs[8] = { 155929000., 20938850., 2949713., 499656., 100995.,  23855., 6391., 2821.};
  //double NQCD[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};
  double NQCD[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 20000., 20000.};

  TFile * W[11];
  W[0] = new TFile ("./TDAna_W0w_tk3.root");
  W[1] = new TFile ("./TDAna_W10w_tk3.root");
  W[2] = new TFile ("./TDAna_W11w_tk3.root");
  W[3] = new TFile ("./TDAna_W20w_tk3.root");
  W[4] = new TFile ("./TDAna_W21w_tk3.root");
  W[5] = new TFile ("./TDAna_W30w_tk3.root");
  W[6] = new TFile ("./TDAna_W31w_tk3.root");
  W[7] = new TFile ("./TDAna_W40w_tk3.root");
  W[8] = new TFile ("./TDAna_W41w_tk3.root");
  W[9] = new TFile ("./TDAna_W50w_tk3.root");
  W[10] = new TFile ("./TDAna_W51w_tk3.root");
  double Wxs[11] = { 45000., 9200., 250., 2500., 225., 590., 100., 125., 40., 85., 40. };
  double NW[11] = { 88000., 40000., 100530., 99523., 105255., 79000., 
		    88258., 83038., 30796., 59022., 41865. };

  TFile * TTH = new TFile("./TDAna_ttH_120_tk3.root");
  double TTHxs = 0.667 ;
  double NTTH = 1000000.; // 1652000.; // 96000.;

  TFile * TT[5];
  TT[0] = new TFile("./smallrun/TDAna_TT0_tk3.root");
  TT[1] = new TFile("./smallrun/TDAna_TT1_tk3.root");
  TT[2] = new TFile("./smallrun/TDAna_TT2_tk3.root");
  TT[3] = new TFile("./smallrun/TDAna_TT3_tk3.root");
  TT[4] = new TFile("./smallrun/TDAna_TT4_tk3.root");
  // double TTxs[5] = { 619., 176., 34.,  6., 1.5 };  // from web
  double TTxs[5] = { 434., 162., 43., 10., 1.9 };     // from note
  double NTT[5] = { 20000., 20000., 20000., 14768., 5304. };

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
