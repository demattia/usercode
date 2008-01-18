Fit_DEt() {

  TFile * G = new TFile ("./root/TDAna_TTH_120.root");
  G->cd();

  double xmax=600;
  gStyle->SetOptFit(111111);
  TCanvas * h = new TCanvas ("h", "Energy correction functions", 600, 500 );
  h->Divide(2,2);
  TF1 * P2 = new TF1 ("P2", "[0]*x*x+[1]*x+[2]",0,1000);
  P2->SetParameters(1,1,0);
  P2->SetLineColor(kRed);
  P2->SetLineWidth(1);
  h->cd(1);
  DEtcb_prof->Draw("PE");
  h->cd(2);
  DEtb_prof->SetMarkerStyle(20);
  DEtb_prof->SetMarkerSize(0.5);
  DEtb_prof->Fit("P2","","",0,xmax);
  h->cd(3);
  DEtcb_prof->Draw("PE");
  h->cd(4);
  DEtq_prof->SetMarkerStyle(24);
  DEtq_prof->SetMarkerSize(0.5);
  DEtq_prof->Fit("P2","","",0,xmax);
  h->Print("./ps/Fit_DEt.ps");

}
