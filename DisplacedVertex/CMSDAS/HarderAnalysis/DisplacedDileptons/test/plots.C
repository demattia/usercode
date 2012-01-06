TH1F* add_hist(TH1F* h1, TH1F* h2) {
  if (h1==0) {
    return h2;
  } else if (h2==0) {
    return h1;
  } else {
    TH1F* h3 = h1->Clone();
    h3->Add(h2);
    return h3;
  }
}

void plots() {
  TCanvas* canv = new TCanvas();

  string channel="mu";
  string mode="mod";

  TFile* efile = new TFile((channel+"Histograms_"+mode+"reco3.root").c_str());
  efile->cd();
  if (channel=="mu") {
    efile->cd("analyzeDisplacedDileptons/dimuons");
  } else {
    efile->cd("analyzeDisplacedDileptons/dielectrons");
  }
  gDirectory->ls();
  TH1F* hist_all = (TH1F*)(gDirectory->Get("mass"));
  TH1F* hist_partial1   = (TH1F*)(gDirectory->Get("mass_signal_mixed"));
  TH1F* hist_partial2   = (TH1F*)(gDirectory->Get("mass_onesignal_oneother"));
  TH1F* hist_partial3   = (TH1F*)(gDirectory->Get("mass_other_resonance"));
  TH1F* hist_background = (TH1F*)(gDirectory->Get("mass_background"));

  gStyle->SetOptStat(0);
  canv->SetLogy();

  hist_all->SetFillStyle(1001);
  hist_all->SetFillColor(kYellow);
  if (channel=="mu") {
    hist_all->GetXaxis()->SetTitle("dimuon invariant mass (GeV/c^{2})");
  } else {
    hist_all->GetXaxis()->SetTitle("dielectron invariant mass (GeV/c^{2})");
  }
  hist_all->Draw("h0");

  if (hist_background) {
    hist_background->SetFillStyle(1001);
    hist_background->SetFillColor(kRed);
    hist_background->Draw("h0same");
  }

  TH1F* p1 = add_hist(hist_background,hist_partial1);
  TH1F* p2 = add_hist(p1,hist_partial2);
  TH1F* p3 = add_hist(p2,hist_partial3);

  if (p3) {
    p3->SetFillStyle(1001);
    p3->SetFillColor(kCyan);
    p3->Draw("h0same");
  }

  if (p2) {
    p2->SetFillStyle(1001);
    p2->SetFillColor(kBlue);
    p2->Draw("h0same");
  }

  if (p1) {
    p1->SetFillStyle(1001);
    p1->SetFillColor(kGreen);
    p1->Draw("h0same");
  }

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(hist_all,"fully reconstructed signal");
  if (hist_partial1) leg->AddEntry(p1,"swapped leptons");
  if (hist_partial2) leg->AddEntry(p2,"partial signal");
  if (hist_partial3) leg->AddEntry(p3,"other resonances");
  if (hist_background) leg->AddEntry(hist_background,"background");
  leg->Draw();

  canv->Update();
  canv->Print((channel+channel+"_mass_"+mode+".eps").c_str());
  canv->SetLogy(0);


  TH1F* effi = (TH1F*)gDirectory->Get("dilepton_reco_efficiency_decayLength2D");
  if (channel=="mu") {
    effi->SetTitle(("offline reconstruction efficiency of dimuons, "+mode+" reco").c_str());
  } else {
    effi->SetTitle(("offline reconstruction efficiency of dielectrons, "+mode+" reco").c_str());
  }
  effi->GetXaxis()->SetTitle("generated 2D decay length (cm)");
  effi->GetYaxis()->SetTitle("efficiency");
  effi->Scale(3); // work around normalisation bug: normalise to mu channel only
  effi->Draw();
  effi->SetMinimum(0);
  effi->SetMaximum(1.1);
  canv->Update();
  canv->Print((channel+channel+"_recoeffi_"+mode+".eps").c_str());

  TH1F* matcheffi = (TH1F*)gDirectory->Get("dilepton_match_efficiency_decayLength2D");
  matcheffi->SetTitle(("efficiency to have at least one trigger matched "+channel).c_str());
  matcheffi->GetXaxis()->SetTitle("generated 2D decay length (cm)");
  matcheffi->GetYaxis()->SetTitle("efficiency");
  matcheffi->Draw();
  matcheffi->SetMinimum(0);
  matcheffi->SetMaximum(1.1);
  canv->Update();
  canv->Print((channel+channel+"_matcheffi_"+mode+".eps").c_str());


  TH1F* costs = (TH1F*)gDirectory->Get("cosThetaStar");
  costs->SetTitle(("cos(#Theta*) btw positive "+channel+" and di-"+channel).c_str());
  costs->Rebin(5);
  costs->Fit("pol0");
  gStyle->SetOptFit(1);
  costs->GetXaxis()->SetTitle("cos(#Theta*)");
  costs->Draw();
  canv->Update();
  canv->Print((channel+channel+"_costs_"+mode+".eps").c_str());

  gApplication->Terminate();
}
