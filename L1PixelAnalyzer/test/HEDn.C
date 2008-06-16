void HED ( ) {
  
  const int nbins=200;
  
  TFile * QCD[8];     
  QCD[0] = new TFile("./root/TDAna_QCD30-50_tk3now.root");
  QCD[1] = new TFile("./root/TDAna_QCD50-80_tk3now.root");
  QCD[2] = new TFile("./root/TDAna_QCD80-120_tk3now.root");
  QCD[3] = new TFile("./root/TDAna_QCD120-170_tk3now.root");
  QCD[4] = new TFile("./root/TDAna_QCD170-230_tk3now.root");
  QCD[5] = new TFile("./root/TDAna_QCD230-300_tk3now.root");
  QCD[6] = new TFile("./root/TDAna_QCD300-380_tk3now.root");
  QCD[7] = new TFile("./root/TDAna_QCD380incl_tk3now.root");
  double QCDxs[8] = { 155929000., 20938850., 2949713., 499656., 100995.,  23855., 6391., 2821.};
  double NQCD[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};
  
  double Lumfactor = 100000; // 100/fb of luminosity assumed in histograms
  
  TH1D * H = dynamic_cast<TH1D*> (QCD[0]->Get("HED1"));
  
  cout << "All ok" << endl;
  double minx=H->GetBinLowEdge(1);
  double maxx=nbins*H->GetBinWidth(1)+minx;
  TH1D * HED[8]; 
  TH1D * HPD[8]; 
  for ( int i=0; i<8; i++ ) {
    char Ename[20];
    char Pname[20];
    sprintf ( Ename, "HED_%d", i );
    sprintf ( Pname, "HPD_%d", i );
    HED[i] = new TH1D ( Ename, Ename, nbins, minx, maxx );
    HPD[i] = new TH1D ( Pname, Pname, nbins, minx, maxx );
  }
  
  for ( int iN=0; iN<8; iN++ ) {
    char nameE[20];
    char nameP[20];
    sprintf ( nameE, "HED%d", iN+1 );
    sprintf ( nameP, "HPD%d", iN+1 );
    double tot1[nbins]={0.};
    double s2_tot1[nbins]={0.};
    double tot2[nbins]={0.};
    double s2_tot2[nbins]={0.};
    for ( int i=0; i<8; i++ ) {
      HEtmp = dynamic_cast<TH1D*>(QCD[i]->Get(nameE));
      HPtmp = dynamic_cast<TH1D*>(QCD[i]->Get(nameP));
      for ( int ibin=1; ibin<=nbins; ibin++ ) {
	double t1=HEtmp->GetBinContent(ibin);
	tot1[ibin-1]+=t1*QCDxs[i]/NQCD[i]*Lumfactor;
	s2_tot1[ibin-1]+=t1*pow(QCDxs[i]/NQCD[i]*Lumfactor,2);
	double t2=HPtmp->GetBinContent(ibin);
	tot2[ibin-1]+=t2*QCDxs[i]/NQCD[i]*Lumfactor;
	s2_tot2[ibin-1]+=t2*pow(QCDxs[i]/NQCD[i]*Lumfactor,2);
      }
    }
    
    // Now renormalize histograms
    // --------------------------
    double i1=0.;
    double i2=0.;
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      i1+=tot1[ibin-1];
      i2+=tot2[ibin-1];
    }
    for ( int ibin=1; ibin<=nbins; ibin++ ) {
      HED[iN]->SetBinContent(ibin,tot1[ibin-1]/i1);
      HED[iN]->SetBinError(ibin,sqrt(s2_tot1[ibin-1])/i1);
      HPD[iN]->SetBinContent(ibin,tot2[ibin-1]/i2);
      HPD[iN]->SetBinError(ibin,sqrt(s2_tot2[ibin-1])/i2);
    }
  } // iN
  
  TCanvas * b = new TCanvas ("b", "Kinematics comparison", 700, 700 );
  b->Divide(2,2);
  
  b->cd(1);
  b->GetPad(1)->SetLogy();
  HED[0]->SetMarkerStyle(20);
  HED[0]->SetMarkerSize(0.4);
  HED[0]->SetMarkerColor(kBlue);
  HED[0]->SetLineColor(kBlue);
  HED[0]->DrawCopy("PE");
  HED[1]->SetMarkerStyle(21);
  HED[1]->SetMarkerSize(0.4);
  HED[1]->SetMarkerColor(kRed);
  HED[1]->SetLineColor(kRed);
  HED[1]->DrawCopy("PESAME");
  HED[2]->SetMarkerStyle(24);
  HED[2]->SetMarkerSize(0.4);
  HED[2]->SetMarkerColor(kBlack);
  HED[2]->SetLineColor(kBlack);
  HED[2]->DrawCopy("PESAME");
  HED[3]->SetMarkerStyle(25);
  HED[3]->SetMarkerSize(0.4);
  HED[3]->SetMarkerColor(kGreen);
  HED[3]->SetLineColor(kGreen);
  HED[3]->DrawCopy("PESAME");
  b->cd(2);
  b->GetPad(2)->SetLogy();
  HED[4]->SetMarkerStyle(20);
  HED[4]->SetMarkerSize(0.4);
  HED[4]->SetMarkerColor(kBlue);
  HED[4]->SetLineColor(kBlue);
  HED[4]->DrawCopy("PE");
  HED[5]->SetMarkerStyle(21);
  HED[5]->SetMarkerSize(0.4);
  HED[5]->SetMarkerColor(kRed);
  HED[5]->SetLineColor(kRed);
  HED[5]->DrawCopy("PESAME");
  HED[6]->SetMarkerStyle(24);
  HED[6]->SetMarkerSize(0.4);
  HED[6]->SetMarkerColor(kBlack);
  HED[6]->SetLineColor(kBlack);
  HED[6]->DrawCopy("PESAME");
  HED[7]->SetMarkerStyle(25);
  HED[7]->SetMarkerSize(0.4);
  HED[7]->SetMarkerColor(kGreen);
  HED[7]->SetLineColor(kGreen);
  HED[7]->DrawCopy("PESAME");
  b->cd(3);
  b->GetPad(3)->SetLogy();
  HPD[0]->SetMinimum(0.00000000001);
  HPD[0]->SetMarkerStyle(20);
  HPD[0]->SetMarkerSize(0.4);
  HPD[0]->SetMarkerColor(kBlue);
  HPD[0]->SetLineColor(kBlue);
  HPD[0]->DrawCopy("PE");
  HPD[1]->SetMarkerStyle(21);
  HPD[1]->SetMarkerSize(0.4);
  HPD[1]->SetMarkerColor(kRed);
  HPD[1]->SetLineColor(kRed);
  HPD[1]->DrawCopy("PESAME");
  HPD[2]->SetMarkerStyle(24);
  HPD[2]->SetMarkerSize(0.4);
  HPD[2]->SetMarkerColor(kBlack);
  HPD[2]->SetLineColor(kBlack);
  HPD[2]->DrawCopy("PESAME");
  HPD[3]->SetMarkerStyle(25);
  HPD[3]->SetMarkerSize(0.4);
  HPD[3]->SetMarkerColor(kGreen);
  HPD[3]->SetLineColor(kGreen);
  HPD[3]->DrawCopy("PESAME");
  b->cd(4);
  b->GetPad(4)->SetLogy();
  HPD[4]->SetMarkerStyle(20);
  HPD[4]->SetMarkerSize(0.4);
  HPD[4]->SetMarkerColor(kBlue);
  HPD[4]->SetLineColor(kBlue);
  HPD[4]->DrawCopy("PE");
  HPD[5]->SetMarkerStyle(21);
  HPD[5]->SetMarkerSize(0.4);
  HPD[5]->SetMarkerColor(kRed);
  HPD[5]->SetLineColor(kRed);
  HPD[5]->DrawCopy("PESAME");
  HPD[6]->SetMarkerStyle(24);
  HPD[6]->SetMarkerSize(0.4);
  HPD[6]->SetMarkerColor(kBlack);
  HPD[6]->SetLineColor(kBlack);
  HPD[6]->DrawCopy("PESAME");
  HPD[7]->SetMarkerStyle(25);
  HPD[7]->SetMarkerSize(0.4);
  HPD[7]->SetMarkerColor(kGreen);
  HPD[7]->SetLineColor(kGreen);
  HPD[7]->DrawCopy("PESAME");
  b->Print("./ps/HEDn.ps");

  // Close files
  // -----------
  TFile * File = new TFile ("HEDn.root","RECREATE");
  File->cd();
  for ( int i=0; i<8; i++ ) {
    HED[i]->Write();
    HPD[i]->Write();
  }
  File->Close();

}
