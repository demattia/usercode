{

  gStyle->SetOptFit(111111);

  // This macro extracts information on the MET resolution by fitting the combination
  // of several TProfiles
  // --------------------------------------------------------------------------------
  TFile * QCD[8]; 
  QCD[0] = new TFile("./root/TDAna_QCD50-80_tk3.root");
  QCD[1] = new TFile("./root/TDAna_QCD50-80_tk3.root");
  QCD[2] = new TFile("./root/TDAna_QCD80-120_tk3.root");
  QCD[3] = new TFile("./root/TDAna_QCD120-170_tk3.root");
  QCD[4] = new TFile("./root/TDAna_QCD170-230_tk3.root");
  QCD[5] = new TFile("./root/TDAna_QCD230-300_tk3.root");
  QCD[6] = new TFile("./root/TDAna_QCD300-380_tk3.root");
  QCD[7] = new TFile("./root/TDAna_QCD380incl_tk3.root");
  double QCDxs[8] = { 155929000., 20938850., 2949713., 499656., 100995.,  23855., 6391., 2821.};
  double N[8] = { 86000., 78000., 104000., 96000., 100000., 102000., 112000., 102000.};

  //double QCDxs[8] = { 1., 1., 1., 1., 1., 1., 1., 1. };
  //double N[8] = { 1., 1., 1., 1., 1., 1., 1., 1. };


  TH1F * MEtRes = new TH1F("MEtRes", "MEtRes", 100, 0, 4000 );
  TH1F * UncorrMEtRes = new TH1F("UncorrMEtRes", "MEtRes", 100, 0, 4000 );
  TH1F * CorrMEtRes = new TH1F("CorrMEtRes", "MEtRes", 100, 0, 4000 );
  TH1F * MEtResC = new TH1F("MEtResC", "MEtRes", 100, 0, 4000 );
  TH1F * UncorrMEtResC = new TH1F("UncorrMEtResC", "MEtRes", 100, 0, 4000 );
  TH1F * CorrMEtResC = new TH1F("CorrMEtResC", "MEtRes", 100, 0, 4000 );
  TH1F * MEtResJ = new TH1F("MEtResJ", "MEtRes", 100, 0, 4000 );
  TH1F * UncorrMEtResJ = new TH1F("UncorrMEtResJ", "MEtRes", 100, 0, 4000 );
  TH1F * CorrMEtResJ = new TH1F("CorrMEtResJ", "MEtRes", 100, 0, 4000 );
   
  double varU[100]={0.};
  double totU[100]={0.};
  double varC[100]={0.};
  double totC[100]={0.};
  double varM[100]={0.};
  double totM[100]={0.};
  double varUC[100]={0.};
  double totUC[100]={0.};
  double varCC[100]={0.};
  double totCC[100]={0.};
  double varMC[100]={0.};
  double totMC[100]={0.};
  double varUJ[100]={0.};
  double totUJ[100]={0.};
  double varCJ[100]={0.};
  double totCJ[100]={0.};
  double varMJ[100]={0.};
  double totMJ[100]={0.};
  double N2;
  for ( int i=0; i<8; i++ ) {
    QCD[i]->cd();
    for ( int ibin=1; ibin<=100; ibin++ ) {
      for ( int jbin=1; jbin<=100; jbin++ ) {
	N2=UncorrMEt_SumEt->GetBinContent(ibin,jbin);
	varU[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totU[ibin]+=N2*QCDxs[i]/N[i];
	N2=CorrMEt_SumEt->GetBinContent(ibin,jbin);
	varC[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totC[ibin]+=N2*QCDxs[i]/N[i];
	N2=MEt_SumEt->GetBinContent(ibin,jbin);
	varM[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totM[ibin]+=N2*QCDxs[i]/N[i];
	N2=UncorrMEt_SumEtC->GetBinContent(ibin,jbin);
	varUC[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totUC[ibin]+=N2*QCDxs[i]/N[i];
	N2=CorrMEt_SumEtC->GetBinContent(ibin,jbin);
	varCC[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totCC[ibin]+=N2*QCDxs[i]/N[i];
	N2=MEt_SumEtC->GetBinContent(ibin,jbin);
	varMC[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totMC[ibin]+=N2*QCDxs[i]/N[i];
	N2=UncorrMEt_SumEtJ->GetBinContent(ibin,jbin);
	varUJ[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totUJ[ibin]+=N2*QCDxs[i]/N[i];
	N2=CorrMEt_SumEtJ->GetBinContent(ibin,jbin);
	varCJ[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totCJ[ibin]+=N2*QCDxs[i]/N[i];
	N2=MEt_SumEtJ->GetBinContent(ibin,jbin);
	varMJ[ibin]+=N2*QCDxs[i]/N[i]*pow((double)jbin*10,2);
	totMJ[ibin]+=N2*QCDxs[i]/N[i];
      }
    }
  }
  for ( int ibin=1; ibin<=100; ibin++ ) {
    double content;
    double error;
    if ( totU[ibin]>0 ) {
      content = sqrt(varU[ibin]/totU[ibin]);
      error = sqrt(varU[ibin])/totU[ibin];
    }
    UncorrMEtRes->SetBinContent(ibin,content);
    UncorrMEtRes->SetBinError(ibin,error);
    if ( totC[ibin]>0 ) {
      content = sqrt(varC[ibin]/totC[ibin]);
      error = sqrt(varC[ibin])/totC[ibin];
    }
    CorrMEtRes->SetBinContent(ibin,content);
    CorrMEtRes->SetBinError(ibin,error);
    if ( totM[ibin]>0 ) {
      content = sqrt(varM[ibin]/totM[ibin]);
      error = sqrt(varM[ibin])/totM[ibin];
    }
    MEtRes->SetBinContent(ibin,content);
    MEtRes->SetBinError(ibin,error);
    if ( totUC[ibin]>0 ) {
      content = sqrt(varUC[ibin]/totUC[ibin]);
      error = sqrt(varUC[ibin])/totUC[ibin];
    }
    UncorrMEtResC->SetBinContent(ibin,content);
    UncorrMEtResC->SetBinError(ibin,error);
    if ( totCC[ibin]>0 ) {
      content = sqrt(varCC[ibin]/totCC[ibin]);
      error = sqrt(varCC[ibin])/totCC[ibin];
    }
    CorrMEtResC->SetBinContent(ibin,content);
    CorrMEtResC->SetBinError(ibin,error);
    if ( totMC[ibin]>0 ) {
      content = sqrt(varMC[ibin]/totMC[ibin]);
      error = sqrt(varMC[ibin])/totMC[ibin];
    }
    MEtResC->SetBinContent(ibin,content);
    MEtResC->SetBinError(ibin,error);
    if ( totUJ[ibin]>0 ) {
      content = sqrt(varUJ[ibin]/totUJ[ibin]);
      error = sqrt(varUJ[ibin])/totUJ[ibin];
    }
    UncorrMEtResJ->SetBinContent(ibin,content);
    UncorrMEtResJ->SetBinError(ibin,error);
    if ( totCJ[ibin]>0 ) {
      content = sqrt(varCJ[ibin]/totCJ[ibin]);
      error = sqrt(varCJ[ibin])/totCJ[ibin];
    }
    CorrMEtResJ->SetBinContent(ibin,content);
    CorrMEtResJ->SetBinError(ibin,error);
    if ( totMJ[ibin]>0 ) {
      content = sqrt(varMJ[ibin]/totMJ[ibin]);
      error = sqrt(varMJ[ibin])/totMJ[ibin];
    }
    MEtResJ->SetBinContent(ibin,content);
    MEtResJ->SetBinError(ibin,error);
  }

  TCanvas * b = new TCanvas ("b", "missing Et comparisons", 700, 700 );
  b->Divide(3,3);
  b->cd(1);
  //  UncorrMEtRes->Fit("Line","V","",500.,4000.);
  UncorrMEtRes->SetMaximum(200);
  UncorrMEtRes->SetMinimum(0);
  UncorrMEtRes->Draw();
  b->cd(2);
//   Line->SetParameters(15,0.03);
//   CorrMEtRes->Fit("Line","V","",500.,4000.);
  CorrMEtRes->SetMaximum(200);
  CorrMEtRes->SetMinimum(0);
  CorrMEtRes->Draw();
  b->cd(3);
//   Line->SetParameters(15,0.03);
//   MEtRes->Fit("Line","V","",500.,4000.);
  MEtRes->SetMaximum(200);
  MEtRes->SetMinimum(0);
  MEtRes->Draw();
  b->cd(4);
//   Line->SetParameters(10,0.03);
//   UncorrMEtResC->Fit("Line","V","",500.,4000.);
  UncorrMEtResC->SetMaximum(200);
  UncorrMEtResC->SetMinimum(0);
  UncorrMEtResC->Draw();
  b->cd(5);
//   Line->SetParameters(15,0.03);
//   CorrMEtResC->Fit("Line","V","",500.,4000.);
  CorrMEtResC->SetMaximum(200);
  CorrMEtResC->SetMinimum(0);
  CorrMEtResC->Draw();
  b->cd(6);
//   Line->SetParameters(15,0.03);
//   MEtResC->Fit("Line","V","",500.,4000.);
  MEtResC->SetMaximum(200);
  MEtResC->SetMinimum(0);
  MEtResC->Draw();
  b->cd(7);
//   Line->SetParameters(20,0.04);
//   UncorrMEtResJ->Fit("Line","V","",500.,2000.);
  UncorrMEtResJ->SetMaximum(200);
  UncorrMEtResJ->SetMinimum(0);
  UncorrMEtResJ->Draw();
  b->cd(8);
//   Line->SetParameters(40,0.04);
//   CorrMEtResJ->Fit("Line","V","",500.,2000.);
  CorrMEtResJ->SetMaximum(200);
  CorrMEtResJ->SetMinimum(0);
  CorrMEtResJ->Draw();
  b->cd(9);
//   Line->SetParameters(20,0.04);
//   MEtResJ->Fit("Line","V","",500.,2000.);
  MEtResJ->SetMaximum(200);
  MEtResJ->SetMinimum(0);
  MEtResJ->Draw();
  b->Print("MetResFits2_all.ps");


  TF1 * Line = new TF1 ( "Line", "[0]+[1]*x", 0, 4000 );
  Line->SetLineWidth(1);
  Line->SetLineColor(kRed);
  Line->SetParameters(10,0.03);

  TCanvas * b2 = new TCanvas ("b2", "missing Et comparisons", 700, 700 );
  b2->cd();
  MEtRes->Fit("Line", "V", "", 500., 2000.);
  MEtRes->Draw("PE");

  b2->Print("MetResFits2_line.ps");


}
