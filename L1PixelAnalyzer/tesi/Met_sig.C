{

  // This macro extracts information on the MET resolution by fitting the combination
  // of several TProfiles
  // --------------------------------------------------------------------------------
  TFile * TTH = new TFile("TDAna_ttH_120_tk3.root");
  //TFile * TTH = new TFile("TT_all.root");

  TFile * QCD[8]; 
  QCD[0] = new TFile("TDAna_QCD30-50_tk3.root");
  QCD[1] = new TFile("TDAna_QCD50-80_tk3.root");
  QCD[2] = new TFile("TDAna_QCD80-120_tk3.root");
  QCD[3] = new TFile("TDAna_QCD120-170_tk3.root");
  QCD[4] = new TFile("TDAna_QCD170-230_tk3.root");
  QCD[5] = new TFile("TDAna_QCD230-300_tk3.root");
  QCD[6] = new TFile("TDAna_QCD300-380_tk3.root");
  QCD[7] = new TFile("TDAna_QCD380incl_tk3.root");
  double QCDxs[8] = { 155929000., 20938850., 2949713., 499656., 100995.,  23855., 6391., 2821.};
  double N[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};
    
  TH1F * S_ = new TH1F ( "S", "S", 100, 0, 20 );
  TH1F * B_ = new TH1F ( "B", "B", 100, 0, 20 );
  TH1F * ES_ = new TH1F ( "ES", "ES", 100, 0, 20 );
  TH1F * EB_ = new TH1F ( "EB", "EB", 100, 0, 20 );  
  TH1F * R_ = new TH1F ( "R", "R", 100, 0, 20 );
  TH1F * Q_ = new TH1F ( "Q", "Q", 100, 0, 20 );

  TTH->cd();
  double es[100]={0.};
  double s_es[100]={0.};
  double IS=MEtSigS->Integral();
  double IScurr=IS;
  double maxS=0;
  for ( int ibin=1; ibin<=100; ibin++ ) {
    es[ibin-1] = IScurr/IS;
    s_es[ibin-1] = sqrt(es[ibin-1]*(1-es[ibin-1])/IS);
    IScurr -= MEtSigS->GetBinContent(ibin);
    ES_->SetBinContent(ibin,es[ibin-1]);
    ES_->SetBinError(ibin,s_es[ibin-1]);
    S_->SetBinContent(ibin,MEtSigS->GetBinContent(ibin)/IS);
    S_->SetBinError(ibin,MEtSigS->GetBinError(ibin)/IS);
    if ( MEtSigS->GetBinContent(ibin)/IS>maxS ) maxS=MEtSigS->GetBinContent(ibin)/IS;
  }

  double B[100]={0.};
  double s2_B[100]={0.};
  double IB = 0;
  double s2_IB = 0;
  double IBcurr = 0;
  double s2_IBcurr = 0;
  double eb[100];
  double s_eb[100];
  double maxB=0;
  for ( int i=0; i<8; i++ ) {
    QCD[i]->cd();
    for ( int ibin=1; ibin<=100; ibin++ ) {
      B[ibin-1]+=QCDxs[i]/N[i]*MEtSigS->GetBinContent(ibin);
      s2_B[ibin-1]+=(QCDxs[i]/N[i])*(QCDxs[i]/N[i])*MEtSigS->GetBinContent(ibin);
      IB += QCDxs[i]/N[i]*MEtSigS->GetBinContent(ibin);
      s2_IB += (QCDxs[i]/N[i])*(QCDxs[i]/N[i])*MEtSigS->GetBinContent(ibin);
    }
  }

  for ( int ibin=100; ibin>=1; ibin-- ) {
    IBcurr += B[ibin-1];
    s2_IBcurr += s2_B[ibin-1];
    eb[ibin-1] = IBcurr/IB;
    s_eb[ibin-1] = eb[ibin-1]*sqrt( (s2_IBcurr/IB/IB) + (s2_IB*IBcurr*IBcurr/IB/IB/IB/IB) ); 
    EB_->SetBinContent(ibin,eb[ibin-1]);
    EB_->SetBinError(ibin,s_eb[ibin-1]);
    B_->SetBinContent(ibin,B[ibin-1]/IB);
    B_->SetBinError(ibin,sqrt(s2_B[ibin-1])/IB);
    if ( B[ibin-1]/IB>maxB ) maxB=B[ibin-1]/IB;
  }

  double R;
  double s_R;
  double Q;
  double s_Q;
  double maxQm1s=0;
  double s_maxQm1s=0;
  double x_maxQ=0;
  for ( int ibin=1; ibin<=100; ibin++ ) {
    if ( eb[ibin-1]>0 ) {
      R = es[ibin-1]/eb[ibin-1];
      s_R = sqrt(s_es[ibin-1]*s_es[ibin-1]/eb[ibin-1]/eb[ibin-1]+
		 es[ibin-1]*es[ibin-1]*s_eb[ibin-1]*s_eb[ibin-1]/pow(eb[ibin-1],4));
      Q = es[ibin-1]*es[ibin-1]/eb[ibin-1];
      s_Q = sqrt(pow(2*es[ibin-1]/eb[ibin-1]*s_es[ibin-1],2)+pow(es[ibin-1]/eb[ibin-1],4)*pow(s_eb[ibin-1],2));
    } else {
      R = es[ibin-1]/0.1;
      s_R = sqrt(pow(s_es[ibin-1]/0.1,2)+
		 pow(es[ibin-1]*s_eb[ibin-1],2)/pow(0.1,4));
      Q = es[ibin-1]*es[ibin-1]/0.1;
      s_Q = sqrt(pow(2*es[ibin-1]/0.1*s_es[ibin-1],2)+pow(es[ibin-1]/0.1,4)*pow(s_eb[ibin-1],2));
    }
    R_->SetBinContent(ibin,R);
    R_->SetBinError(ibin,s_R);
    Q_->SetBinContent(ibin,Q);
    Q_->SetBinError(ibin,s_Q);
    if ( maxQm1s<Q-s_Q ) {
      maxQm1s=Q;
      s_maxQm1s=s_Q;
      x_maxQ=(double)ibin*20./100.;
    } 
  }
  double maxX=maxS;
  if (maxB>maxX) maxX=maxB;
  cout << endl;
  cout << "Maximum quality factor: " << maxQm1s << "+-" << s_maxQm1s << " at x > " << x_maxQ << endl; 
  cout << "--------------------------------------------------------------------" << endl;

  TCanvas * b = new TCanvas ("b", "missing Et comparisons", 500, 500 );
  b->Divide(2,2);
  b->cd(1);
  S_->SetMaximum(1.2*maxX);
  S_->SetLineColor(kRed);
  S_->SetMarkerColor(kRed);
  S_->Draw("PE");
  B_->SetLineColor(kBlue);
  B_->SetMarkerColor(kBlue);
  B_->Draw("PESAME");
  b->cd(2);
  ES_->SetLineColor(kRed);
  ES_->SetMarkerColor(kRed);
  ES_->Draw("PE");
  EB_->SetLineColor(kBlue);
  EB_->SetMarkerColor(kBlue);
  EB_->Draw("PESAME");
  b->cd(3);
  R_->Draw("PE");
  b->cd(4);
  Q_->Draw("PE");
  b->Print("MEtSigS_opt.ps");
}
