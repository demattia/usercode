{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TFile f("acceptance.root");

  TH2F* acc = gROOT->FindObject("acceptance");
  acc->GetZaxis()->SetRangeUser(-0.001,1.001);
  TCanvas c;
  acc->Draw("colz");
  c.SaveAs("acceptance.png");

  TH2F* acc1D = gROOT->FindObject("acceptance1D");
  acc1D->GetYaxis()->SetRangeUser(-0.001,1.001);
  TCanvas c1D;
  acc1D->Draw();
  c1D.SaveAs("acceptance1D.png");
}

