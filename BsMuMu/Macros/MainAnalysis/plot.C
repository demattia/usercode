// Plots the invariant mass distribution in the range used by the analysis.

void plot(const TString & index = "0")
{
  TFile *_file0 = TFile::Open("data_afterCuts_"+index+".root");
  tree = (TTree*)gDirectory->Get("events");
  tree->Draw("m>>mass(40, 4.9, 5.9)");
}
