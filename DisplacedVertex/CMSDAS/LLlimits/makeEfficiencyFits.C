// make the efficiency fits
// note: this only uses the statistical+pileup uncertainty,
// the tracking uncertainty is added later

// recommended usage: run parseEfficiencyFiles.pl
// and put the output that it produces in the appropriate place

// this then gives you some code that you can either paste into the
// dataFileXXX files, or in PlotLimitsOld.py (not recommended any more)

void makeEfficiencyFits(void) {
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);

  float eff_1000_x[4] = {20, 50, 150, 350};
  float ex_1000[4] = {0, 0, 0, 0};
  float eff_400_x[4] = {5, 20, 50, 150};
  float ex_400[4] = {0, 0, 0, 0};
  float eff_200_x[2] = {20, 50};
  float ex_200[2] = {0, 0};

  float effe1_1000_y[4] = {0.1549, 0.2274, 0.1983, 0.1063};
  float effe1_1000_ey[4] = {0.0023, 0.0029, 0.0026, 0.0021};

  float effe1_400_y[4] = {0.0000, 0.0976, 0.1030, 0.0543};
  float effe1_400_ey[4] = {0.0000, 0.0018, 0.0019, 0.0015};

  float effe1_200_y[2] = {0.0190, 0.0095};
  float effe1_200_ey[2] = {0.0009, 0.0006};

  float effe1down_1000_y[4] = {0.2146, 0.3161, 0.3188, 0.2204};
  float effe1down_1000_ey[4] = {0.0026, 0.0032, 0.0030, 0.0026};

  float effe1down_400_y[4] = {0.0000, 0.1396, 0.1554, 0.1131};
  float effe1down_400_ey[4] = {0.0000, 0.0021, 0.0022, 0.0021};

  float effe1down_200_y[2] = {0.0326, 0.0178};
  float effe1down_200_ey[2] = {0.0011, 0.0008};

  float effe1up_1000_y[4] = {0.0704, 0.1058, 0.0991, 0.0379};
  float effe1up_1000_ey[4] = {0.0016, 0.0023, 0.0020, 0.0015};

  float effe1up_400_y[4] = {0.0000, 0.0405, 0.0466, 0.0212};
  float effe1up_400_ey[4] = {0.0000, 0.0012, 0.0014, 0.0011};

  float effe1up_200_y[2] = {0.0069, 0.0230};
  float effe1up_200_ey[2] = {0.0005, 0.0009};

  float effe2_1000_y[4] = {0.1698, 0.2307, 0.1999, 0.1057};
  float effe2_1000_ey[4] = {0.0024, 0.0032, 0.0025, 0.0020};

  float effe2_400_y[4] = {0.0000, 0.1037, 0.1075, 0.0530};
  float effe2_400_ey[4] = {0.0000, 0.0020, 0.0019, 0.0015};

  float effe2_200_y[2] = {0.0193, 0.0106};
  float effe2_200_ey[2] = {0.0009, 0.0006};

  float effe2down_1000_y[4] = {0.2396, 0.3122, 0.3441, 0.2176};
  float effe2down_1000_ey[4] = {0.0027, 0.0035, 0.0030, 0.0026};

  float effe2down_400_y[4] = {0.0000, 0.1551, 0.1584, 0.1086};
  float effe2down_400_ey[4] = {0.0000, 0.0023, 0.0023, 0.0020};

  float effe2down_200_y[2] = {0.0290, 0.0220};
  float effe2down_200_ey[2] = {0.0011, 0.0009};

  float effe2up_1000_y[4] = {0.0761, 0.1067, 0.0850, 0.0445};
  float effe2up_1000_ey[4] = {0.0018, 0.0027, 0.0017, 0.0014};

  float effe2up_400_y[4] = {0.0000, 0.0436, 0.0520, 0.0257};
  float effe2up_400_ey[4] = {0.0000, 0.0015, 0.0014, 0.0011};

  float effe2up_200_y[2] = {0.0103, 0.0032};
  float effe2up_200_ey[2] = {0.0007, 0.0004};

  float effmu1_1000_y[4] = {0.0215, 0.3129, 0.4141, 0.2942};
  float effmu1_1000_ey[4] = {0.0009, 0.0030, 0.0041, 0.0029};

  float effmu1_400_y[4] = {0.0000, 0.2023, 0.3256, 0.2326};
  float effmu1_400_ey[4] = {0.0000, 0.0025, 0.0035, 0.0029};

  float effmu1_200_y[2] = {0.1434, 0.1165};
  float effmu1_200_ey[2] = {0.0022, 0.0021};

  float effmu1down_1000_y[4] = {0.0229, 0.4414, 0.5844, 0.5062};
  float effmu1down_1000_ey[4] = {0.0009, 0.0032, 0.0041, 0.0031};

  float effmu1down_400_y[4] = {0.0000, 0.2747, 0.4527, 0.4043};
  float effmu1down_400_ey[4] = {0.0000, 0.0028, 0.0036, 0.0033};

  float effmu1down_200_y[2] = {0.1980, 0.1893};
  float effmu1down_200_ey[2] = {0.0025, 0.0025};

  float effmu1up_1000_y[4] = {0.0171, 0.1691, 0.1932, 0.1198};
  float effmu1up_1000_ey[4] = {0.0008, 0.0024, 0.0035, 0.0021};

  float effmu1up_400_y[4] = {0.0000, 0.1014, 0.1660, 0.0941};
  float effmu1up_400_ey[4] = {0.0000, 0.0020, 0.0031, 0.0022};

  float effmu1up_200_y[2] = {0.0658, 0.0538};
  float effmu1up_200_ey[2] = {0.0016, 0.0015};

  float effmu2_1000_y[4] = {0.0332, 0.3976, 0.4347, 0.3049};
  float effmu2_1000_ey[4] = {0.0011, 0.0031, 0.0042, 0.0029};

  float effmu2_400_y[4] = {0.0000, 0.2610, 0.3469, 0.2503};
  float effmu2_400_ey[4] = {0.0000, 0.0029, 0.0045, 0.0028};

  float effmu2_200_y[2] = {0.1643, 0.1304};
  float effmu2_200_ey[2] = {0.0024, 0.0023};

  float effmu2down_1000_y[4] = {0.0370, 0.5599, 0.6078, 0.5263};
  float effmu2down_1000_ey[4] = {0.0012, 0.0032, 0.0042, 0.0031};

  float effmu2down_400_y[4] = {0.0000, 0.3694, 0.4844, 0.4608};
  float effmu2down_400_ey[4] = {0.0000, 0.0031, 0.0046, 0.0032};

  float effmu2down_200_y[2] = {0.2293, 0.2104};
  float effmu2down_200_ey[2] = {0.0027, 0.0028};

  float effmu2up_1000_y[4] = {0.0217, 0.2033, 0.2031, 0.1229};
  float effmu2up_1000_ey[4] = {0.0009, 0.0026, 0.0038, 0.0021};

  float effmu2up_400_y[4] = {0.0000, 0.1259, 0.1554, 0.1018};
  float effmu2up_400_ey[4] = {0.0000, 0.0022, 0.0041, 0.0020};

  float effmu2up_200_y[2] = {0.0744, 0.0504};
  float effmu2up_200_ey[2] = {0.0017, 0.0017};

  TCanvas *c1 = new TCanvas("c1", "Electron efficiency1", 800, 800);
  TGraphErrors *effe1_1000 = new TGraphErrors(4, eff_1000_x, effe1_1000_y, ex_1000, effe1_1000_ey);
  effe1_1000->Draw("ALP");
  effe1_1000->SetMinimum(0);
  effe1_1000->Fit("pol2","S");

  TGraphErrors *effe1_400 = new TGraphErrors(4, eff_400_x, effe1_400_y, ex_400, effe1_400_ey);
  effe1_400->Draw("LP same");
  effe1_400->SetLineColor(kBlue);
  effe1_400->Fit("pol2");

  TGraphErrors *effe1_200 = new TGraphErrors(2, eff_200_x, effe1_200_y, ex_200, effe1_200_ey);
  effe1_200->Draw("LP same");
  effe1_200->SetLineColor(kRed);
  effe1_200->Fit("pol1");

  TCanvas *c2 = new TCanvas("c2", "Electron efficiency2", 800, 800);
  TGraphErrors *effe2_1000 = new TGraphErrors(4, eff_1000_x, effe2_1000_y, ex_1000, effe2_1000_ey);
  effe2_1000->Draw("ALP");
  effe2_1000->SetMinimum(0);
  effe2_1000->Fit("pol2");

  TGraphErrors *effe2_400 = new TGraphErrors(4, eff_400_x, effe2_400_y, ex_400, effe2_400_ey);
  effe2_400->Draw("LP same");
  effe2_400->SetLineColor(kBlue);
  effe2_400->Fit("pol2");

  TGraphErrors *effe2_200 = new TGraphErrors(2, eff_200_x, effe2_200_y, ex_200, effe2_200_ey);
  effe2_200->Draw("LP same");
  effe2_200->SetLineColor(kBlue);
  effe2_200->Fit("pol1");

  TCanvas *c3 = new TCanvas("c3", "Muon efficiency1", 800, 800);
  TGraphErrors *effmu1_1000 = new TGraphErrors(4, eff_1000_x, effmu1_1000_y, ex_1000, effmu1_1000_ey);
  effmu1_1000->Draw("ALP");
  effmu1_1000->SetMinimum(0);
  effmu1_1000->Fit("pol2");

  TGraphErrors *effmu1_400 = new TGraphErrors(4, eff_400_x, effmu1_400_y, ex_400, effmu1_400_ey);
  effmu1_400->Draw("LP same");
  effmu1_400->SetLineColor(kBlue);
  effmu1_400->Fit("pol2");

  TGraphErrors *effmu1_200 = new TGraphErrors(2, eff_200_x, effmu1_200_y, ex_200, effmu1_200_ey);
  effmu1_200->Draw("LP same");
  effmu1_200->SetLineColor(kRed);
  effmu1_200->Fit("pol1");

  TCanvas *c4 = new TCanvas("c4", "Muon efficiency2", 800, 800);
  TGraphErrors *effmu2_1000 = new TGraphErrors(4, eff_1000_x, effmu2_1000_y, ex_1000, effmu2_1000_ey);
  effmu2_1000->Draw("ALP");
  effmu2_1000->SetMinimum(0);
  effmu2_1000->Fit("pol2");

  TGraphErrors *effmu2_400 = new TGraphErrors(4, eff_400_x, effmu2_400_y, ex_400, effmu2_400_ey);
  effmu2_400->Draw("LP same");
  effmu2_400->SetLineColor(kBlue);
  effmu2_400->Fit("pol2");

  TGraphErrors *effmu2_200 = new TGraphErrors(2, eff_200_x, effmu2_200_y, ex_200, effmu2_200_ey);
  effmu2_200->Draw("LP same");
  effmu2_200->SetLineColor(kBlue);
  effmu2_200->Fit("pol1");

  TCanvas *c5 = new TCanvas("c5", "Electron efficiency down1", 800, 800);
  TGraphErrors *effe1down_1000 = new TGraphErrors(4, eff_1000_x, effe1down_1000_y, ex_1000, effe1down_1000_ey);
  effe1down_1000->Draw("ALP");
  effe1down_1000->SetMinimum(0);
  effe1down_1000->Fit("pol2","S");

  TGraphErrors *effe1down_400 = new TGraphErrors(4, eff_400_x, effe1down_400_y, ex_400, effe1down_400_ey);
  effe1down_400->Draw("LP same");
  effe1down_400->SetLineColor(kBlue);
  effe1down_400->Fit("pol2");

  TGraphErrors *effe1down_200 = new TGraphErrors(2, eff_200_x, effe1down_200_y, ex_200, effe1down_200_ey);
  effe1down_200->Draw("LP same");
  effe1down_200->SetLineColor(kRed);
  effe1down_200->Fit("pol1");

  TCanvas *c6 = new TCanvas("c6", "Electron efficiency down2", 800, 800);
  TGraphErrors *effe2down_1000 = new TGraphErrors(4, eff_1000_x, effe2down_1000_y, ex_1000, effe2down_1000_ey);
  effe2down_1000->Draw("ALP");
  effe2down_1000->SetMinimum(0);
  effe2down_1000->Fit("pol2");

  TGraphErrors *effe2down_400 = new TGraphErrors(4, eff_400_x, effe2down_400_y, ex_400, effe2down_400_ey);
  effe2down_400->Draw("LP same");
  effe2down_400->SetLineColor(kBlue);
  effe2down_400->Fit("pol2");

  TGraphErrors *effe2down_200 = new TGraphErrors(2, eff_200_x, effe2down_200_y, ex_200, effe2down_200_ey);
  effe2down_200->Draw("LP same");
  effe2down_200->SetLineColor(kBlue);
  effe2down_200->Fit("pol1");

  TCanvas *c7 = new TCanvas("c7", "Muon efficiency down1", 800, 800);
  TGraphErrors *effmu1down_1000 = new TGraphErrors(4, eff_1000_x, effmu1down_1000_y, ex_1000, effmu1down_1000_ey);
  effmu1down_1000->Draw("ALP");
  effmu1down_1000->SetMinimum(0);
  effmu1down_1000->Fit("pol2");

  TGraphErrors *effmu1down_400 = new TGraphErrors(4, eff_400_x, effmu1down_400_y, ex_400, effmu1down_400_ey);
  effmu1down_400->Draw("LP same");
  effmu1down_400->SetLineColor(kBlue);
  effmu1down_400->Fit("pol2");

  TGraphErrors *effmu1down_200 = new TGraphErrors(2, eff_200_x, effmu1down_200_y, ex_200, effmu1down_200_ey);
  effmu1down_200->Draw("LP same");
  effmu1down_200->SetLineColor(kRed);
  effmu1down_200->Fit("pol1");

  TCanvas *c8 = new TCanvas("c8", "Muon efficiency down2", 800, 800);
  TGraphErrors *effmu2down_1000 = new TGraphErrors(4, eff_1000_x, effmu2down_1000_y, ex_1000, effmu2down_1000_ey);
  effmu2down_1000->Draw("ALP");
  effmu2down_1000->SetMinimum(0);
  effmu2down_1000->Fit("pol2");

  TGraphErrors *effmu2down_400 = new TGraphErrors(4, eff_400_x, effmu2down_400_y, ex_400, effmu2down_400_ey);
  effmu2down_400->Draw("LP same");
  effmu2down_400->SetLineColor(kBlue);
  effmu2down_400->Fit("pol2");

  TGraphErrors *effmu2down_200 = new TGraphErrors(2, eff_200_x, effmu2down_200_y, ex_200, effmu2down_200_ey);
  effmu2down_200->Draw("LP same");
  effmu2down_200->SetLineColor(kBlue);
  effmu2down_200->Fit("pol1");

  TCanvas *c9 = new TCanvas("c9", "Electron efficiency up1", 800, 800);
  TGraphErrors *effe1up_1000 = new TGraphErrors(4, eff_1000_x, effe1up_1000_y, ex_1000, effe1up_1000_ey);
  effe1up_1000->Draw("ALP");
  effe1up_1000->SetMinimum(0);
  effe1up_1000->Fit("pol2","S");

  TGraphErrors *effe1up_400 = new TGraphErrors(4, eff_400_x, effe1up_400_y, ex_400, effe1up_400_ey);
  effe1up_400->Draw("LP same");
  effe1up_400->SetLineColor(kBlue);
  effe1up_400->Fit("pol2");

  TGraphErrors *effe1up_200 = new TGraphErrors(2, eff_200_x, effe1up_200_y, ex_200, effe1up_200_ey);
  effe1up_200->Draw("LP same");
  effe1up_200->SetLineColor(kRed);
  effe1up_200->Fit("pol1");

  TCanvas *c10 = new TCanvas("c10", "Electron efficiency up2", 800, 800);
  TGraphErrors *effe2up_1000 = new TGraphErrors(4, eff_1000_x, effe2up_1000_y, ex_1000, effe2up_1000_ey);
  effe2up_1000->Draw("ALP");
  effe2up_1000->SetMinimum(0);
  effe2up_1000->Fit("pol2");

  TGraphErrors *effe2up_400 = new TGraphErrors(4, eff_400_x, effe2up_400_y, ex_400, effe2up_400_ey);
  effe2up_400->Draw("LP same");
  effe2up_400->SetLineColor(kBlue);
  effe2up_400->Fit("pol2");

  TGraphErrors *effe2up_200 = new TGraphErrors(2, eff_200_x, effe2up_200_y, ex_200, effe2up_200_ey);
  effe2up_200->Draw("LP same");
  effe2up_200->SetLineColor(kBlue);
  effe2up_200->Fit("pol1");

  TCanvas *c11 = new TCanvas("c11", "Muon efficiency up1", 800, 800);
  TGraphErrors *effmu1up_1000 = new TGraphErrors(4, eff_1000_x, effmu1up_1000_y, ex_1000, effmu1up_1000_ey);
  effmu1up_1000->Draw("ALP");
  effmu1up_1000->SetMinimum(0);
  effmu1up_1000->Fit("pol2");

  TGraphErrors *effmu1up_400 = new TGraphErrors(4, eff_400_x, effmu1up_400_y, ex_400, effmu1up_400_ey);
  effmu1up_400->Draw("LP same");
  effmu1up_400->SetLineColor(kBlue);
  effmu1up_400->Fit("pol2");

  TGraphErrors *effmu1up_200 = new TGraphErrors(2, eff_200_x, effmu1up_200_y, ex_200, effmu1up_200_ey);
  effmu1up_200->Draw("LP same");
  effmu1up_200->SetLineColor(kRed);
  effmu1up_200->Fit("pol1");

  TCanvas *c12 = new TCanvas("c12", "Muon efficiency up2", 800, 800);
  TGraphErrors *effmu2up_1000 = new TGraphErrors(4, eff_1000_x, effmu2up_1000_y, ex_1000, effmu2up_1000_ey);
  effmu2up_1000->Draw("ALP");
  effmu2up_1000->SetMinimum(0);
  effmu2up_1000->Fit("pol2");

  TGraphErrors *effmu2up_400 = new TGraphErrors(4, eff_400_x, effmu2up_400_y, ex_400, effmu2up_400_ey);
  effmu2up_400->Draw("LP same");
  effmu2up_400->SetLineColor(kBlue);
  effmu2up_400->Fit("pol2");

  TGraphErrors *effmu2up_200 = new TGraphErrors(2, eff_200_x, effmu2up_200_y, ex_200, effmu2up_200_ey);
  effmu2up_200->Draw("LP same");
  effmu2up_200->SetLineColor(kBlue);
  effmu2up_200->Fit("pol1");

  std::cout << "electron 1000" << std::endl;
  TF1 *e1_1000_fit = effe1_1000->GetFunction("pol2");
  std::cout << "SigEff1 " << e1_1000_fit->GetParameter(0) << " "
	    << e1_1000_fit->GetParameter(1) << " " << e1_1000_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff1E " << e1_1000_fit->GetParError(0) << " "
	    << e1_1000_fit->GetParError(1) << " " << e1_1000_fit->GetParError(2) << std::endl;
  TF1 *e2_1000_fit = effe2_1000->GetFunction("pol2");
  std::cout << "SigEff2 " << e2_1000_fit->GetParameter(0) << " "
	    << e2_1000_fit->GetParameter(1) << " " << e2_1000_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff2E " << e2_1000_fit->GetParError(0) << " "
	    << e2_1000_fit->GetParError(1) << " " << e2_1000_fit->GetParError(2) << std::endl;

  std::cout << "electron 400" << std::endl;
  TF1 *e1_400_fit = effe1_400->GetFunction("pol2");
  std::cout << "SigEff1 " << e1_400_fit->GetParameter(0) << " "
	    << e1_400_fit->GetParameter(1) << " " << e1_400_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff1E " << e1_400_fit->GetParError(0) << " "
	    << e1_400_fit->GetParError(1) << " " << e1_400_fit->GetParError(2) << std::endl;
  TF1 *e2_400_fit = effe2_400->GetFunction("pol2");
  std::cout << "SigEff2 " << e2_400_fit->GetParameter(0) << " "
	    << e2_400_fit->GetParameter(1) << " " << e2_400_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff2E " << e2_400_fit->GetParError(0) << " "
	    << e2_400_fit->GetParError(1) << " " << e2_400_fit->GetParError(2) << std::endl;

  std::cout << "electron 200" << std::endl;
  TF1 *e1_200_fit = effe1_200->GetFunction("pol1");
  std::cout << "SigEff1 " << e1_200_fit->GetParameter(0) << " "
	    << e1_200_fit->GetParameter(1) << std::endl;
  std::cout << "SigEff1E " << e1_200_fit->GetParError(0) << " "
	    << e1_200_fit->GetParError(1) << std::endl;
  TF1 *e2_200_fit = effe2_200->GetFunction("pol1");
  std::cout << "SigEff2 " << e2_200_fit->GetParameter(0) << " "
	    << e2_200_fit->GetParameter(1) << std::endl;
  std::cout << "SigEff2E " << e2_200_fit->GetParError(0) << " "
	    << e2_200_fit->GetParError(1) << std::endl;

  std::cout << "muon 1000" << std::endl;
  TF1 *mu1_1000_fit = effmu1_1000->GetFunction("pol2");
  std::cout << "SigEff1 " << mu1_1000_fit->GetParameter(0) << " "
	    << mu1_1000_fit->GetParameter(1) << " " << mu1_1000_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff1E " << mu1_1000_fit->GetParError(0) << " "
	    << mu1_1000_fit->GetParError(1) << " " << mu1_1000_fit->GetParError(2) << std::endl;
  TF1 *mu2_1000_fit = effmu2_1000->GetFunction("pol2");
  std::cout << "SigEff2 " << mu2_1000_fit->GetParameter(0) << " "
	    << mu2_1000_fit->GetParameter(1) << " " << mu2_1000_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff2E " << mu2_1000_fit->GetParError(0) << " "
	    << mu2_1000_fit->GetParError(1) << " " << mu2_1000_fit->GetParError(2) << std::endl;

  std::cout << "muon 400" << std::endl;
  TF1 *mu1_400_fit = effmu1_400->GetFunction("pol2");
  std::cout << "SigEff1 " << mu1_400_fit->GetParameter(0) << " "
	    << mu1_400_fit->GetParameter(1) << " " << mu1_400_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff1E " << mu1_400_fit->GetParError(0) << " "
	    << mu1_400_fit->GetParError(1) << " " << mu1_400_fit->GetParError(2) << std::endl;
  TF1 *mu2_400_fit = effmu2_400->GetFunction("pol2");
  std::cout << "SigEff2 " << mu2_400_fit->GetParameter(0) << " "
	    << mu2_400_fit->GetParameter(1) << " " << mu2_400_fit->GetParameter(2) << std::endl;
  std::cout << "SigEff2E " << mu2_400_fit->GetParError(0) << " "
	    << mu2_400_fit->GetParError(1) << " " << mu2_400_fit->GetParError(2) << std::endl;

  std::cout << "muon 200" << std::endl;
  TF1 *mu1_200_fit = effmu1_200->GetFunction("pol1");
  std::cout << "SigEff1 " << mu1_200_fit->GetParameter(0) << " "
	    << mu1_200_fit->GetParameter(1) << std::endl;
  std::cout << "SigEff1E " << mu1_200_fit->GetParError(0) << " "
	    << mu1_200_fit->GetParError(1) << std::endl;
  TF1 *mu2_200_fit = effmu2_200->GetFunction("pol1");
  std::cout << "SigEff2 " << mu2_200_fit->GetParameter(0) << " "
	    << mu2_200_fit->GetParameter(1) << std::endl;
  std::cout << "SigEff2E " << mu2_200_fit->GetParError(0) << " "
	    << mu2_200_fit->GetParError(1) << std::endl;

  TF1 *mu1down_1000_fit = effmu1down_1000->GetFunction("pol2");
  TF1 *mu2down_1000_fit = effmu2down_1000->GetFunction("pol2");
  TF1 *mu1down_400_fit = effmu1down_400->GetFunction("pol2");
  TF1 *mu2down_400_fit = effmu2down_400->GetFunction("pol2");
  TF1 *mu1down_200_fit = effmu1down_200->GetFunction("pol1");
  TF1 *mu2down_200_fit = effmu2down_200->GetFunction("pol1");
  TF1 *e1down_1000_fit = effe1down_1000->GetFunction("pol2");
  TF1 *e2down_1000_fit = effe2down_1000->GetFunction("pol2");
  TF1 *e1down_400_fit = effe1down_400->GetFunction("pol2");
  TF1 *e2down_400_fit = effe2down_400->GetFunction("pol2");
  TF1 *e1down_200_fit = effe1down_200->GetFunction("pol1");
  TF1 *e2down_200_fit = effe2down_200->GetFunction("pol1");
  TF1 *mu1up_1000_fit = effmu1up_1000->GetFunction("pol2");
  TF1 *mu2up_1000_fit = effmu2up_1000->GetFunction("pol2");
  TF1 *mu1up_400_fit = effmu1up_400->GetFunction("pol2");
  TF1 *mu2up_400_fit = effmu2up_400->GetFunction("pol2");
  TF1 *mu1up_200_fit = effmu1up_200->GetFunction("pol1");
  TF1 *mu2up_200_fit = effmu2up_200->GetFunction("pol1");
  TF1 *e1up_1000_fit = effe1up_1000->GetFunction("pol2");
  TF1 *e2up_1000_fit = effe2up_1000->GetFunction("pol2");
  TF1 *e1up_400_fit = effe1up_400->GetFunction("pol2");
  TF1 *e2up_400_fit = effe2up_400->GetFunction("pol2");
  TF1 *e1up_200_fit = effe1up_200->GetFunction("pol1");
  TF1 *e2up_200_fit = effe2up_200->GetFunction("pol1");

  std::cout << "if type == 1:" << std::endl;
  std::cout << "eff1 = " << mu1_1000_fit->GetParameter(0) << " + " << mu1_1000_fit->GetParameter(1) << "*mass "
	    << mu1_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << mu2_1000_fit->GetParameter(0) << " + " << mu2_1000_fit->GetParameter(1) << "*mass "
	    << mu2_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << mu1_1000_fit->GetParError(0) << " + " << mu1_1000_fit->GetParError(1) << "*mass + "
	    << mu1_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << mu2_1000_fit->GetParError(0) << " + " << mu2_1000_fit->GetParError(1) << "*mass + "
	    << mu2_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 2:" << std::endl;
  std::cout << "eff1 = " << mu1_400_fit->GetParameter(0) << " + " << mu1_400_fit->GetParameter(1) << "*mass "
	    << mu1_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << mu2_400_fit->GetParameter(0) << " + " << mu2_400_fit->GetParameter(1) << "*mass "
	    << mu2_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << mu1_400_fit->GetParError(0) << " + " << mu1_400_fit->GetParError(1) << "*mass + "
	    << mu1_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << mu2_400_fit->GetParError(0) << " + " << mu2_400_fit->GetParError(1) << "*mass + "
	    << mu2_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 3:" << std::endl;
  std::cout << "eff1 = " << mu1_200_fit->GetParameter(0) << " " << mu1_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff2 = " << mu2_200_fit->GetParameter(0) << " " << mu2_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff1e = " << mu1_200_fit->GetParError(0) << " + " << mu1_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "eff2e = " << mu2_200_fit->GetParError(0) << " + " << mu2_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "elif type == 5:" << std::endl;
  std::cout << "eff1 = " << e1_1000_fit->GetParameter(0) << " + " << e1_1000_fit->GetParameter(1) << "*mass "
	    << e1_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << e2_1000_fit->GetParameter(0) << " + " << e2_1000_fit->GetParameter(1) << "*mass "
	    << e2_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << e1_1000_fit->GetParError(0) << " + " << e1_1000_fit->GetParError(1) << "*mass + "
	    << e1_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << e2_1000_fit->GetParError(0) << " + " << e2_1000_fit->GetParError(1) << "*mass + "
	    << e2_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 6:" << std::endl;
  std::cout << "eff1 = " << e1_400_fit->GetParameter(0) << " + " << e1_400_fit->GetParameter(1) << "*mass "
	    << e1_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << e2_400_fit->GetParameter(0) << " + " << e2_400_fit->GetParameter(1) << "*mass "
	    << e2_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << e1_400_fit->GetParError(0) << " + " << e1_400_fit->GetParError(1) << "*mass + "
	    << e1_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << e2_400_fit->GetParError(0) << " + " << e2_400_fit->GetParError(1) << "*mass + "
	    << e2_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 7:" << std::endl;
  std::cout << "eff1 = " << e1_200_fit->GetParameter(0) << " " << e1_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff2 = " << e2_200_fit->GetParameter(0) << " " << e2_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff1e = " << e1_200_fit->GetParError(0) << " + " << e1_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "eff2e = " << e2_200_fit->GetParError(0) << " + " << e2_200_fit->GetParError(1) << "*mass"
	    << std::endl;
    
  std::cout << "elif type == 11:" << std::endl;
  std::cout << "eff1 = " << mu1down_1000_fit->GetParameter(0) << " + " << mu1down_1000_fit->GetParameter(1) << "*mass "
	    << mu1down_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << mu2down_1000_fit->GetParameter(0) << " + " << mu2down_1000_fit->GetParameter(1) << "*mass "
	    << mu2down_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << mu1down_1000_fit->GetParError(0) << " + " << mu1down_1000_fit->GetParError(1) << "*mass + "
	    << mu1down_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << mu2down_1000_fit->GetParError(0) << " + " << mu2down_1000_fit->GetParError(1) << "*mass + "
	    << mu2down_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 12:" << std::endl;
  std::cout << "eff1 = " << mu1down_400_fit->GetParameter(0) << " + " << mu1down_400_fit->GetParameter(1) << "*mass "
	    << mu1down_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << mu2down_400_fit->GetParameter(0) << " + " << mu2down_400_fit->GetParameter(1) << "*mass "
	    << mu2down_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << mu1down_400_fit->GetParError(0) << " + " << mu1down_400_fit->GetParError(1) << "*mass + "
	    << mu1down_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << mu2down_400_fit->GetParError(0) << " + " << mu2down_400_fit->GetParError(1) << "*mass + "
	    << mu2down_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 13:" << std::endl;
  std::cout << "eff1 = " << mu1down_200_fit->GetParameter(0) << " " << mu1down_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff2 = " << mu2down_200_fit->GetParameter(0) << " " << mu2down_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff1e = " << mu1down_200_fit->GetParError(0) << " + " << mu1down_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "eff2e = " << mu2down_200_fit->GetParError(0) << " + " << mu2down_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "elif type == 15:" << std::endl;
  std::cout << "eff1 = " << e1down_1000_fit->GetParameter(0) << " + " << e1down_1000_fit->GetParameter(1) << "*mass "
	    << e1down_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << e2down_1000_fit->GetParameter(0) << " + " << e2down_1000_fit->GetParameter(1) << "*mass "
	    << e2down_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << e1down_1000_fit->GetParError(0) << " + " << e1down_1000_fit->GetParError(1) << "*mass + "
	    << e1down_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << e2down_1000_fit->GetParError(0) << " + " << e2down_1000_fit->GetParError(1) << "*mass + "
	    << e2down_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 16:" << std::endl;
  std::cout << "eff1 = " << e1down_400_fit->GetParameter(0) << " + " << e1down_400_fit->GetParameter(1) << "*mass "
	    << e1down_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << e2down_400_fit->GetParameter(0) << " + " << e2down_400_fit->GetParameter(1) << "*mass "
	    << e2down_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << e1down_400_fit->GetParError(0) << " + " << e1down_400_fit->GetParError(1) << "*mass + "
	    << e1down_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << e2down_400_fit->GetParError(0) << " + " << e2down_400_fit->GetParError(1) << "*mass + "
	    << e2down_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 17:" << std::endl;
  std::cout << "eff1 = " << e1down_200_fit->GetParameter(0) << " " << e1down_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff2 = " << e2down_200_fit->GetParameter(0) << " " << e2down_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff1e = " << e1down_200_fit->GetParError(0) << " + " << e1down_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "eff2e = " << e2down_200_fit->GetParError(0) << " + " << e2down_200_fit->GetParError(1) << "*mass"
	    << std::endl;
    
  std::cout << "elif type == 21:" << std::endl;
  std::cout << "eff1 = " << mu1up_1000_fit->GetParameter(0) << " + " << mu1up_1000_fit->GetParameter(1) << "*mass "
	    << mu1up_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << mu2up_1000_fit->GetParameter(0) << " + " << mu2up_1000_fit->GetParameter(1) << "*mass "
	    << mu2up_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << mu1up_1000_fit->GetParError(0) << " + " << mu1up_1000_fit->GetParError(1) << "*mass + "
	    << mu1up_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << mu2up_1000_fit->GetParError(0) << " + " << mu2up_1000_fit->GetParError(1) << "*mass + "
	    << mu2up_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 22:" << std::endl;
  std::cout << "eff1 = " << mu1up_400_fit->GetParameter(0) << " + " << mu1up_400_fit->GetParameter(1) << "*mass "
	    << mu1up_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << mu2up_400_fit->GetParameter(0) << " + " << mu2up_400_fit->GetParameter(1) << "*mass "
	    << mu2up_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << mu1up_400_fit->GetParError(0) << " + " << mu1up_400_fit->GetParError(1) << "*mass + "
	    << mu1up_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << mu2up_400_fit->GetParError(0) << " + " << mu2up_400_fit->GetParError(1) << "*mass + "
	    << mu2up_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 23:" << std::endl;
  std::cout << "eff1 = " << mu1up_200_fit->GetParameter(0) << " " << mu1up_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff2 = " << mu2up_200_fit->GetParameter(0) << " " << mu2up_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff1e = " << mu1up_200_fit->GetParError(0) << " + " << mu1up_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "eff2e = " << mu2up_200_fit->GetParError(0) << " + " << mu2up_200_fit->GetParError(1) << "*mass"
	    << std::endl;

  std::cout << "elif type == 25:" << std::endl;
  std::cout << "eff1 = " << e1up_1000_fit->GetParameter(0) << " + " << e1up_1000_fit->GetParameter(1) << "*mass "
	    << e1up_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << e2up_1000_fit->GetParameter(0) << " + " << e2up_1000_fit->GetParameter(1) << "*mass "
	    << e2up_1000_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << e1up_1000_fit->GetParError(0) << " + " << e1up_1000_fit->GetParError(1) << "*mass + "
	    << e1up_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << e2up_1000_fit->GetParError(0) << " + " << e2up_1000_fit->GetParError(1) << "*mass + "
	    << e2up_1000_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 26:" << std::endl;
  std::cout << "eff1 = " << e1up_400_fit->GetParameter(0) << " + " << e1up_400_fit->GetParameter(1) << "*mass "
	    << e1up_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff2 = " << e2up_400_fit->GetParameter(0) << " + " << e2up_400_fit->GetParameter(1) << "*mass "
	    << e2up_400_fit->GetParameter(2) << "*mass*mass" << std::endl;
  std::cout << "eff1e = " << e1up_400_fit->GetParError(0) << " + " << e1up_400_fit->GetParError(1) << "*mass + "
	    << e1up_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "eff2e = " << e2up_400_fit->GetParError(0) << " + " << e2up_400_fit->GetParError(1) << "*mass + "
	    << e2up_400_fit->GetParError(2) << "*mass*mass" << std::endl;
  std::cout << "elif type == 27:" << std::endl;
  std::cout << "eff1 = " << e1up_200_fit->GetParameter(0) << " " << e1up_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff2 = " << e2up_200_fit->GetParameter(0) << " " << e2up_200_fit->GetParameter(1) << "*mass"
	    << std::endl;
  std::cout << "eff1e = " << e1up_200_fit->GetParError(0) << " + " << e1up_200_fit->GetParError(1) << "*mass"
	    << std::endl;
  std::cout << "eff2e = " << e2up_200_fit->GetParError(0) << " + " << e2up_200_fit->GetParError(1) << "*mass"
	    << std::endl;
    
}
