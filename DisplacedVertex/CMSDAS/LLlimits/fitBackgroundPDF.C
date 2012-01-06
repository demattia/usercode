{
  TFile *_file0 = TFile::Open("limits/masses_backgroundMC_electron.root");
  TH1F *h = (TH1F*)_file0->Get("backgroundmass");
  TF1 *f = new TF1("f", "gaus(0) + [3]", 50, 350); // or change as necessary
  f->SetParameter(0, 50);
  f->SetParameter(1, 90);
  f->SetParameter(2, 2);
  f->SetParameter(3, 1);
  h->Fit(f, "LE"); // can also use chi2 by removing the LE
}
