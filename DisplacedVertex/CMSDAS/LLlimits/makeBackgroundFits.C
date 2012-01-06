// Make the background fits.

const char *eleFile = "/afs/cern.ch/user/p/plujan/FromKristian/plots_with_root_files/Data_Photon_Run2011A2_20110813_HLT_DoublePhoton33_HEVT_v3/dielectron_decayLengthSignificance2D.root";
const char *muonFile = "/afs/cern.ch/user/p/plujan/FromKristian/plots_with_root_files/Data_Mu_Run2011A2_20110813_HLT_L2DoubleMu23_NoVertex_v4/dimuon_decayLengthSignificance2D.root";
  
void plotSingleBackground(TFile *f, float xfit, float xcut, bool isMu) {
  TH1F *h = f->Get("decayLengthSignificance2D_BackgroundMC");

  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(kFullSquare);

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);
  gPad->SetRightMargin(0.04);

  h->GetXaxis()->SetTitle("L_{xy}/#sigma");
  h->GetYaxis()->SetTitle("Number of background events");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->Draw();
  gPad->SetLogy();
  
  float xmax = h->GetXaxis()->GetXmax();
  float ymax = h->GetMaximum();

  TF1 *f1 = new TF1("f1", "[0]*exp(-x/[1]) + [2]*exp(-x/[3])");
  if (isMu) {
    f1->SetParameter(0, 750);
    f1->SetParameter(1, 0.25);
    f1->SetParameter(2, 90);
    f1->SetParameter(3, 0.6);
  } else {
    f1->SetParameter(0, 35);
    f1->SetParameter(1, 0.4);
    f1->SetParameter(2, 2);
    f1->SetParameter(3, 3);
  }

  float fitmin = 0;
  //if (!isMu) fitmin = 3;
  h->Fit(f1, "EM", "", fitmin, xfit);

  f1->SetLineColor(kRed);

  TF1 *f2 = new TF1("f2", "[0]*exp(-x/[1]) + [2]*exp(-x/[3])", xfit, xmax);
  for (int i=0; i<=3; i++)
    f2->SetParameter(i, f1->GetParameter(i));
  f2->SetLineColor(kRed);
  f2->SetLineStyle(kDashed);
  f2->Draw("same");
  
  if (!isMu) {
    TLine *l1 = new TLine(xfit, h->GetMinimum(), xfit, ymax);
    l1->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l1->SetLineStyle(kDashed);
    l1->Draw();
  }

  TLine *l2 = new TLine(xcut, h->GetMinimum(), xcut, ymax);
  l2->SetLineColor(kBlue);
  l2->SetLineWidth(2);
  l2->SetLineStyle(kDashed);
  l2->Draw();

  TArrow *a2 = new TArrow(xcut, ymax/50, xcut+2, ymax/50, 0.02, "|>");
  a2->SetLineWidth(2);
  a2->SetLineColor(kBlue);
  a2->SetFillColor(kBlue);
  a2->Draw();
  TText *at2 = new TText(xcut+2.5, ymax/50, "Signal region");
  at2->SetTextSize(0.03);
  at2->SetTextColor(kBlue);
  at2->SetTextAlign(12);
  at2->Draw();

  TLatex *t1;
  if (isMu) {
    h->SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.2 fb^{-1}"); 
    t1 = new TLatex(xmax-4, ymax/5, "#mu^{+}#mu^{-}");
  } else {
    h->SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.1 fb^{-1}"); 
    t1 = new TLatex(xmax-4, ymax/5, "e^{+}e^{-}");
  }
  t1->SetTextFont(42);
  t1->Draw();

  float intbkgnd = f2->Integral(xcut, xmax);

  TF1 *f3 = f2->Clone("f3");
  float interr2 = 0;
  float errplus2 = 0;
  float errminus2 = 0;
  for (int i=0; i<4; i++) {
    f3->SetParameter(i, f1->GetParameter(i)+f1->GetParError(i));
    float intplus = f3->Integral(xcut, xmax);
    f3->SetParameter(i, f1->GetParameter(i)-f1->GetParError(i));
    float intminus = f3->Integral(xcut, xmax);
    f3->SetParameter(i, f1->GetParameter(i));
    float interr = (intplus - intminus)/2;
    interr2 += interr*interr;
    errplus2 += (intplus - intbkgnd)*(intplus - intbkgnd);
    errminus2 += (intbkgnd - intminus)*(intbkgnd - intminus);
  }  
  //std::cout << intplus << " " << intminus << std::endl;
  // float interr = (intplus - intminus)/2;
  float interr = sqrt(interr2);
  float errplus = sqrt(errplus2);
  float errminus = sqrt(errminus2);

  std::cout << "Integral from " << xcut << " to " << xmax << " = " << intbkgnd
	    << " + " << errplus << " - " << errminus
	    << " (+/- " << interr << ")" << std::endl;

}

void makeBackgroundFits(void) {
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  // Some fixes to the TDR style
		   
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(62);

  TFile *fEle = new TFile(eleFile);
  TFile *fMu = new TFile(muonFile);

  c1 = new TCanvas("ele", "Electron channel background fit", 600, 600);
  plotSingleBackground(fEle, 8, 8, false);
  c1->Print("BackgroundFitElectrons.png");
  c1->Print("BackgroundFitElectrons.pdf");
  
  c2 = new TCanvas("mu", "Muon channel background fit", 600, 600);
  plotSingleBackground(fMu, 5, 5, true);
  c2->Print("BackgroundFitMuons.png");
  c2->Print("BackgroundFitMuons.pdf");
}
