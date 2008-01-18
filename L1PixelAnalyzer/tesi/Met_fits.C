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

  TProfile * MEtRes = new TProfile("MEtRes", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * UncorrMEtRes = new TProfile("UncorrMEtRes", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * CorrMEtRes = new TProfile("CorrMEtRes", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * MEtResC = new TProfile("MEtResC", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * UncorrMEtResC = new TProfile("UncorrMEtResC", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * CorrMEtResC = new TProfile("CorrMEtResC", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * MEtResJ = new TProfile("MEtResJ", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * UncorrMEtResJ = new TProfile("UncorrMEtResJ", "MEtRes", 100, 0, 4000, 0, 1000 );
  TProfile * CorrMEtResJ = new TProfile("CorrMEtResJ", "MEtRes", 100, 0, 4000, 0, 1000 );
   
  for ( int i=0; i<8; i++ ) {
    QCD[i]->cd();
    UncorrMEtRes->Add(UncorrMEt_SumEt->ProfileX(),QCDxs[i]/N[i]);
    CorrMEtRes->Add(CorrMEt_SumEt->ProfileX(),QCDxs[i]/N[i]);
    MEtRes->Add(MEt_SumEt->ProfileX(),QCDxs[i]/N[i]);
    UncorrMEtResC->Add(UncorrMEt_SumEtC->ProfileX(),QCDxs[i]/N[i]);
    CorrMEtResC->Add(CorrMEt_SumEtC->ProfileX(),QCDxs[i]/N[i]);
    MEtResC->Add(MEt_SumEtC->ProfileX(),QCDxs[i]/N[i]);
    UncorrMEtResJ->Add(UncorrMEt_SumEtJ->ProfileX(),QCDxs[i]/N[i]);
    CorrMEtResJ->Add(CorrMEt_SumEtJ->ProfileX(),QCDxs[i]/N[i]);
    MEtResJ->Add(MEt_SumEtJ->ProfileX(),QCDxs[i]/N[i]);
  }

  for ( int i=1; i<=100; i++ ) {
    if ( UncorrMEtRes->GetBinError(i)==0 ) UncorrMEtRes->SetBinError(i,UncorrMEtRes->GetBinError(i-1));
    if ( CorrMEtRes->GetBinError(i)==0 ) CorrMEtRes->SetBinError(i,CorrMEtRes->GetBinError(i-1));
    if ( MEtRes->GetBinError(i)==0 ) MEtRes->SetBinError(i,CorrMEtRes->GetBinError(i-1));
    if ( UncorrMEtResC->GetBinError(i)==0 ) UncorrMEtResC->SetBinError(i,UncorrMEtResC->GetBinError(i-1));
    if ( CorrMEtResC->GetBinError(i)==0 ) CorrMEtResC->SetBinError(i,CorrMEtResC->GetBinError(i-1));
    if ( MEtResC->GetBinError(i)==0 ) MEtResC->SetBinError(i,CorrMEtResC->GetBinError(i-1));
    if ( UncorrMEtResJ->GetBinError(i)==0 ) UncorrMEtResJ->SetBinError(i,UncorrMEtResJ->GetBinError(i-1));
    if ( CorrMEtResJ->GetBinError(i)==0 ) CorrMEtResJ->SetBinError(i,CorrMEtResJ->GetBinError(i-1));
    if ( MEtResJ->GetBinError(i)==0 ) MEtResJ->SetBinError(i,CorrMEtResJ->GetBinError(i-1));
  }

  TF1 * Line = new TF1 ( "Line", "[0]+[1]*x", 0, 4000 );
  Line->SetLineWidth(1);
  Line->SetLineColor(kRed);

  TCanvas * b = new TCanvas ("b", "missing Et comparisons", 700, 700 );
  b->Divide(3,3);
  b->cd(1);
  Line->SetParameters(10,0.03);
  UncorrMEtRes->Fit("Line","V","",500.,4000.);
  UncorrMEtRes->SetMaximum(200);
  UncorrMEtRes->Draw();
  b->cd(2);
  Line->SetParameters(15,0.03);
  CorrMEtRes->Fit("Line","V","",500.,4000.);
  CorrMEtRes->SetMaximum(200);
  CorrMEtRes->Draw();
  b->cd(3);
  Line->SetParameters(15,0.03);
  MEtRes->Fit("Line","V","",500.,4000.);
  MEtRes->SetMaximum(200);
  MEtRes->Draw();
  b->cd(4);
  Line->SetParameters(10,0.03);
  UncorrMEtResC->Fit("Line","V","",500.,4000.);
  UncorrMEtResC->SetMaximum(200);
  UncorrMEtResC->Draw();
  b->cd(5);
  Line->SetParameters(15,0.03);
  CorrMEtResC->Fit("Line","V","",500.,4000.);
  CorrMEtResC->SetMaximum(200);
  CorrMEtResC->Draw();
  b->cd(6);
  Line->SetParameters(15,0.03);
  MEtResC->Fit("Line","V","",500.,4000.);
  MEtResC->SetMaximum(200);
  MEtResC->Draw();
  b->cd(7);
  Line->SetParameters(20,0.04);
  UncorrMEtResJ->Fit("Line","V","",500.,2000.);
  UncorrMEtResJ->SetMaximum(200);
  UncorrMEtResJ->Draw();
  b->cd(8);
  Line->SetParameters(40,0.04);
  CorrMEtResJ->Fit("Line","V","",500.,2000.);
  CorrMEtResJ->SetMaximum(200);
  CorrMEtResJ->Draw();
  b->cd(9);
  Line->SetParameters(20,0.04);
  MEtResJ->Fit("Line","V","",500.,2000.);
  MEtResJ->SetMaximum(200);
  MEtResJ->Draw();
  b->Print("MetResFits_line.ps");

//   TF1 * Power = new TF1 ( "Power", "[0]+[1]*pow(x,[2])+[3]*pow(x,[4])", 0, 4000 );
//   Power->SetLineWidth(1);
//   Power->SetLineColor(kRed);

//   TCanvas * b2 = new TCanvas ("b2", "missing Et comparisons", 700, 700 );
//   b2->Divide(3,3);
//   b2->cd(1);
//   Power->SetParameters(10.,0.05,1.0,0.001,0.5);
//   UncorrMEtRes->Fit("Power","V","",0.,4000.);
//   UncorrMEtRes->SetMaximum(200);
//   UncorrMEtRes->Draw();
//   b2->cd(2);
//   Power->SetParameters(10.,0.04,0.9,0.001,0.5);
//   CorrMEtRes->Fit("Power","V","",0.,4000.);
//   CorrMEtRes->SetMaximum(200);
//   CorrMEtRes->Draw();
//   b2->cd(3);
//   Power->SetParameters(10.,0.05,1.0,0.001,0.5);
//   MEtRes->Fit("Power","V","",0.,4000.);
//   MEtRes->SetMaximum(200);
//   MEtRes->Draw();
//   b2->cd(4);
//   Power->SetParameters(10.,0.05,1.0,0.001,0.5);
//   UncorrMEtResC->Fit("Power","V","",0.,4000.);
//   UncorrMEtResC->SetMaximum(200);
//   UncorrMEtResC->Draw();
//   b2->cd(5);
//   Power->SetParameters(10.,0.04,1.0,0.001,0.5);
//   CorrMEtResC->Fit("Power","V","",0.,4000.);
//   CorrMEtResC->SetMaximum(200);
//   CorrMEtResC->Draw();
//   b2->cd(6);
//   Power->SetParameters(10.,0.04,1.0,0.001,0.5);
//   MEtResC->Fit("Power","V","",0.,4000.);
//   MEtResC->SetMaximum(200);
//   MEtResC->Draw();
//   b2->cd(7);
//   Power->SetParameters(30.,0.045,1.0,0.001,0.5);
//   UncorrMEtResJ->Fit("Power","V","",0.,4000.);
//   UncorrMEtResJ->SetMaximum(200);
//   UncorrMEtResJ->Draw();
//   b2->cd(8);
//   Power->SetParameters(40.,0.045,1.0,0.001,0.5);
//   CorrMEtResJ->Fit("Power","V","",0.,4000.);
//   CorrMEtResJ->SetMaximum(200);
//   CorrMEtResJ->Draw();
//   b2->cd(9);
//   Power->SetParameters(20.,0.07,1.0,0.001,0.5);
//   MEtResJ->Fit("Power","V","",0.,4000.);
//   MEtResJ->SetMaximum(200);
//   MEtResJ->Draw();
//   b2->Print("MetResFits_power.ps");

}
