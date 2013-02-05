using namespace RooFit;

void getEff(const char* inputFile, const char* inputDir, const char* inputEff, const char* output, const char* outputLo, const char* outputHi){
  TFile f(inputFile);
  f.cd(inputDir);
  RooDataSet* eff = (RooDataSet*)gROOT->FindObject(inputEff);
  const RooArgSet& set = *(eff->get());
  RooRealVar& pt = (RooRealVar&)set["pt"];
  RooRealVar& eta = (RooRealVar&)set["eta"];
  RooRealVar& e = (RooRealVar&)set["efficiency"];
  //TH2F* h = eff->createHistogram(pt, eta); it works on my laptop but on lxplus???
  TH2F *h = new TH2F("eff", "eff", pt.getBinning().numBins(), pt.getBinning().array(), eta.getBinning().numBins(), eta.getBinning().array());
  TH2F* hlo = h->Clone("eff_lo");
  TH2F* hhi = h->Clone("eff_hi");
  for(int i=0; i<eff->numEntries(); i++){
    eff->get(i);
    int b = h->FindBin(pt.getVal(), eta.getVal());
    h->SetBinContent(b, e.getVal());
    hlo->SetBinContent(b, e.getVal()+e.getErrorLo());
    hhi->SetBinContent(b, e.getVal()+e.getErrorHi());
  }
  TFile out(output,"recreate");
  h->Write("efficiency");
  out.Close();
  TFile outlo(outputLo,"recreate");
  hlo->Write("efficiency");
  outlo.Close();
  TFile outhi(outputHi,"recreate");
  hhi->Write("efficiency");
  outhi.Close();
}

