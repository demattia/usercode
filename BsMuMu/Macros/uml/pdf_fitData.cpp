#include "pdf_fitData.h"

pdf_fitData::pdf_fitData(bool print, string input_estimates, string range, int BF, bool SM, bool bd_constr, int simul, int simulbdt, int simulall, bool pee_, bool bdt_fit, bool final, string ch_s, int sig, bool asimov, bool syste, bool randomsyste, bool rare_constr, int nexp, bool bd, string years): pdf_analysis(print, ch_s, range, BF, SM, bd_constr, simul, simulbdt, simulall, pee_, bdt_fit, final) {
  cout << "fitData constructor" << endl;
  input_estimates_ = input_estimates;
  estimate_bs.resize(channels, -1);
  estimate_bd.resize(channels, -1);
  estimate_peak.resize(channels, -1);
  estimate_semi.resize(channels, -1);
  estimate_comb.resize(channels, -1);

  estimate2D_channel.resize(channels, vector<double> (channels_bdt, -1));
  estimate2D_bs.resize(channels, vector<double> (channels_bdt, -1));
  estimate2D_bd.resize(channels, vector<double> (channels_bdt, -1));
  estimate2D_peak.resize(channels, vector<double> (channels_bdt, -1));
  estimate2D_semi.resize(channels, vector<double> (channels_bdt, -1));
  estimate2D_comb.resize(channels, vector<double> (channels_bdt, -1));

  lumi = -1;

  random = false;

  RooRandom::randomGenerator()->SetSeed(0);

  sign = sig;

  syst = syste;
  randomsyst = randomsyste;
  rare_constr_ = rare_constr;

  /// the pdf
  pdfname = "pdf_ext_total";
  if (simul_ && !syst) pdfname = "pdf_ext_simul_noconstr";
  if (simul_ && syst) pdfname = "pdf_ext_simul";
  cout << red_color_bold << ">>>>>>>>>>>>>>> the name of the fitting pdf is " << pdfname << " <<<<<<<<<<<<<<<<<<<" << default_console_color << endl;

  NExp = nexp;
  Bd = bd;
  SMIsNull = false;
  years_ = years;

  asimov_ = asimov;

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);

  berns_ = false;
  freeze = false;
  stat_only_ = false;
  doubleNull = false;


}

pdf_fitData::~pdf_fitData() {
  cout << "pdf_fitData destructor" << endl;
}

void pdf_fitData::set_estimate() {
	if (simul_ && !simul_bdt_ && !simul_all_) {
		for (int i = 0; i < channels; i++) {
			estimate_bs[i] = ws_->var(pdf_analysis::name("N_bs", i))->getVal();
			estimate_bd[i] = ws_->var(pdf_analysis::name("N_bd", i))->getVal();
			estimate_peak[i] = ws_->var(pdf_analysis::name("N_peak", i))->getVal();
			estimate_semi[i] = ws_->var(pdf_analysis::name("N_semi", i))->getVal();
			estimate_comb[i] = ws_->var(pdf_analysis::name("N_comb", i))->getVal();
		}
	}
	else if (!simul_) {
		int i = atoi(ch_s_.c_str());
		estimate_bs[0] = ws_->var(pdf_analysis::name("N_bs", i))->getVal();
		estimate_bd[0] = ws_->var(pdf_analysis::name("N_bd", i))->getVal();
		estimate_peak[0] = ws_->var(pdf_analysis::name("N_peak", i))->getVal();
		estimate_semi[0] = ws_->var(pdf_analysis::name("N_semi", i))->getVal();
		estimate_comb[0] = ws_->var(pdf_analysis::name("N_comb", i))->getVal();
	}
	else {
		for (int i = 0; i < channels; i++) {
			for (int j = 0; j < bdt_index_max(i); j++) {
				estimate2D_bs[i][j] = ws_->var(pdf_analysis::name("N_bs", i, j))->getVal();
				estimate2D_bd[i][j] = ws_->var(pdf_analysis::name("N_bd", i, j))->getVal();
				estimate2D_peak[i][j] = ws_->var(pdf_analysis::name("N_peak", i, j))->getVal();
				estimate2D_semi[i][j] = ws_->var(pdf_analysis::name("N_semi", i, j))->getVal();
				estimate2D_comb[i][j] = ws_->var(pdf_analysis::name("N_comb", i, j))->getVal();
			}
		}
	}
}

void pdf_fitData::parse_estimate(string input) {
  char buffer[1024];
  char cutName[128];
  double cut;
  string name = input;
  cout << ">>>>>>>>>>>>>> " <<input << endl;
  if (name == "") name = input_estimates_;
  FILE *estimate_file = fopen(name.c_str(), "r");
  if (!estimate_file) {
  	cout << "no estimation from outside files" << endl;
  	lumi = 1;
  	return;
  }
  cout << "event estimates in " << name << " :" << endl;
  while (fgets(buffer, sizeof(buffer), estimate_file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '#') continue;
    sscanf(buffer, "%s %lf", cutName, &cut);
    if (!parse(cutName, cut)) {
      cout << "==> Error parsing variable " << cutName << endl;
      exit(EXIT_FAILURE);
    }
  }
  if (estimate_file) fclose(estimate_file);
  return;
}

bool pdf_fitData::parse(char *cutName, float cut) {
  if (lumi == -1) {
    char test_cut[128];
    sprintf(test_cut, "lumi");
    if (!strcmp(cutName, test_cut)) {
      lumi = (double)cut;
      cout << "lumi is: " << lumi << " times the true luminosity" << endl;
      return true;
    }
  }
  if (simul_ && !simul_bdt_ && !simul_all_) {
    for (int i = 0; i < channels; i++) {
      char test_cut[128];
      sprintf(test_cut, "bs_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_bs[i] = (double)cut;
        cout << "bs[" << i <<"]: " << estimate_bs[i] << endl;
        return true;
      }
      sprintf(test_cut, "bd_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_bd[i] = (double)cut;
        cout << "bd[" << i <<"]: " << estimate_bd[i] << endl;
        return true;
      }
      sprintf(test_cut, "peak_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_peak[i] = (double)cut;
        cout << "peak[" << i <<"]: " << estimate_peak[i] << endl;
        return true;
      }
      sprintf(test_cut, "semi_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_semi[i] = (double)cut;
        cout << "semi[" << i <<"]: " << estimate_semi[i] << endl;
        return true;
      }
      sprintf(test_cut, "comb_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_comb[i] = (double)cut;
        cout << "comb[" << i <<"]: " << estimate_comb[i] << endl;
        return true;
      }
    }
  }
  else if (!simul_) {
    int i = atoi(ch_s_.c_str());
    char test_cut[128];
    sprintf(test_cut, "bs_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_bs[0] = (double)cut;
      cout << "bs[" << i <<"]: " << estimate_bs[0] << endl;
      return true;
    }
    sprintf(test_cut, "bd_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_bd[0] = (double)cut;
      cout << "bd[" << i <<"]: " << estimate_bd[0] << endl;
      return true;
    }
    sprintf(test_cut, "peak_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_peak[0] = (double)cut;
      cout << "peak[" << i <<"]: " << estimate_peak[0] << endl;
      return true;
    }
    sprintf(test_cut, "semi_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_semi[0] = (double)cut;
      cout << "semi[" << i <<"]: " << estimate_semi[0] << endl;
      return true;
    }
    sprintf(test_cut, "comb_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_comb[0] = (double)cut;
      cout << "comb[" << i <<"]: " << estimate_comb[0] << endl;
      return true;
    }
    return true;
  }
  else {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
        char test_cut[128];
        sprintf(test_cut, "bs_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_bs[i][j] = (double)cut;
          cout << "bs[" << i <<"][" << j << "]: " << estimate2D_bs[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "bd_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_bd[i][j] = (double)cut;
          cout << "bd[" << i <<"][" << j << "]: " << estimate2D_bd[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "peak_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_peak[i][j] = (double)cut;
          cout << "peak[" << i <<"][" << j << "]: " << estimate2D_peak[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "semi_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_semi[i][j] = (double)cut;
          cout << "semi[" << i <<"][" << j << "]: " << estimate2D_semi[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "comb_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_comb[i][j] = (double)cut;
          cout << "comb[" << i <<"][" << j << "]: " << estimate2D_comb[i][j] << endl;
          return true;
        }
      }
    }
  }
  return true;
}

void pdf_fitData::print_estimate() {
	if (bdt_fit_) {
		for (int i = 0; i < channels; i++) {cout << name("bdt", i) << endl;
			cout << "minimum " << name("bdt", i) << " = " << ws_->var(name("bdt", i))->getMin() << endl;
		}
	}
	cout << "used estimates" << endl;
	cout << "canale bs bd peak semi comb" << endl;
	if (simul_ && !simul_bdt_ && !simul_all_) {
		for (int i = 0; i < channels; i++) {
			cout << i << " " ;
			if (BF_ == 0) cout << ws_->var(pdf_analysis::name("N_bs", i))->getVal() << " " << ws_->var(pdf_analysis::name("N_bd", i))->getVal();
			else if (BF_ == 1) {
				cout << " " << ws_->function(pdf_analysis::name("N_bs_formula", i))->getVal();
				cout << " " << ws_->var(pdf_analysis::name("N_bd", i))->getVal();
			}
			else {
				cout << " " << ws_->function(pdf_analysis::name("N_bs_formula", i))->getVal();
				cout << " " << ws_->function(pdf_analysis::name("N_bd_formula", i))->getVal();
			}
			cout << " " << ws_->var(pdf_analysis::name("N_peak", i))->getVal() << " " << ws_->var(pdf_analysis::name("N_semi", i))->getVal()  << " " << ws_->var(pdf_analysis::name("N_comb", i))->getVal() << endl;
			cout << " " << ws_->var(pdf_analysis::name("N_peak", i))->getError() << " " << ws_->var(pdf_analysis::name("N_semi", i))->getError()  << " " << ws_->var(pdf_analysis::name("N_comb", i))->getError() << endl << endl;
		}
	}
	else if (!simul_) {
		int i = atoi(ch_s_.c_str());
		cout << i << " " ;
		if (BF_ == 0) cout << ws_->var(pdf_analysis::name("N_bs", i))->getVal() << " " << ws_->var(pdf_analysis::name("N_bd", i))->getVal();
		else if (BF_ == 1) {
			cout << " " << ws_->function(pdf_analysis::name("N_bs_formula", i))->getVal();
			cout << " " << ws_->var(pdf_analysis::name("N_bd", i))->getVal();
		}
		else {
			cout << " " << ws_->function(pdf_analysis::name("N_bs_formula", i))->getVal();
			cout << " " << ws_->function(pdf_analysis::name("N_bd_formula", i))->getVal();
		}
		cout << " " << ws_->var(pdf_analysis::name("N_peak", i))->getVal() << " " << ws_->var(pdf_analysis::name("N_semi", i))->getVal()  << " " << ws_->var(pdf_analysis::name("N_comb", i))->getVal() << endl;
	}
	else {
		for (int i = 0; i < channels; i++) {
			for (int j = 0; j < bdt_index_max(i); j++) {
				cout << i << " " << j << " ";
				if (BF_ == 0) cout << ws_->var(pdf_analysis::name("N_bs", i, j))->getVal() << " " << ws_->var(pdf_analysis::name("N_bd", i, j))->getVal();
				else if (BF_ == 1) {
					cout << " " << ws_->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
					cout << " " << ws_->var(pdf_analysis::name("N_bd", i, j))->getVal();
				}
				else {
					cout << " " << ws_->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
					cout << " " << ws_->function(pdf_analysis::name("N_bd_formula", i, j))->getVal();
				}
				cout << " " << ws_->var(pdf_analysis::name("N_peak", i, j))->getVal() << " " << ws_->var(pdf_analysis::name("N_semi", i, j))->getVal()  << " " << ws_->var(pdf_analysis::name("N_comb", i, j))->getVal() << endl;
			}
		}
	}
	cout << "Bs BF estimated = " << ws_->var("BF_bs")->getVal() << endl;
}

RooFitResult* pdf_fitData::fit_pdf(bool do_not_import, string pdf_name, int strategy, bool extconstr) {
	if (pdf_name == "") pdf_name = pdfname;
  if (!simul_) {
  	pdf_name = "pdf_ext_total";
    RooAbsData* subdata = global_data->reduce(Form("etacat==etacat::etacat_%d", channel));
    global_data = (RooDataSet*)subdata;
  }
  //range_ = "sb_lo,sb_hi";
  cout << red_color_bold << ">>>>>>>>>>>>>>>>> fitting " << global_data->GetName() << " in range " << range_ << " with " << pdf_name << default_console_color << endl;
  global_data->Print();
  ws_->pdf(pdf_name.c_str())->Print();
  RooArgSet ext_constr( *ws_->pdf("fs_over_fu_gau"), *ws_->pdf("one_over_BRBR_gau"));
  RFR = ws_->pdf(pdf_name.c_str())->fitTo(*global_data, extconstr ? ExternalConstraints(ext_constr) : RooCmdArg::none() , Strategy(strategy) ,Extended(), Save(1), minos ? Minos(RooArgSet(*ws_->var("BF_bs"), *ws_->var("BF_bd"))) : Minos(0), Hesse(1));
  if (!do_not_import) ws_->import(*global_data);
  if (verbosity > 0) RFR->Print();

  if (stat_only_) {
  	cout << green_color_bold << "fit with statistical uncertainties only" << default_console_color << endl;
  	stat_error(true);
  	RooFitResult * rfr_stat = ws_->pdf(pdf_name.c_str())->fitTo(*global_data, Strategy(strategy), Extended(), Save(1), minos ? Minos(RooArgSet(*ws_->var("BF_bs"), *ws_->var("BF_bd"))) : Minos(0), Hesse(1));
  	rfr_stat->Print();
  }

  return RFR;
}

void pdf_fitData::print() {
  cout << "printing" << endl;
  RooPlot *rp = ws_->var("Mass")->frame();
  global_data->plotOn(rp, Binning(40));

  if (!pee) {
    ws_->pdf("pdf_ext_total")->plotOn(rp, FillColor(kYellow), Range(range_.c_str()), LineWidth(3), VisualizeError(*RFR), MoveToBack());
    ws_->pdf("pdf_ext_total")->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
  }
  else {
//    ws_->pdf("pdf_ext_total")->plotOn(rp, FillColor(kYellow), Range(range_.c_str()), LineWidth(3), VisualizeError(*RFR), MoveToBack(), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE));
    ws_->pdf("pdf_ext_total")->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE));
  }
  ws_->pdf("pdf_ext_total")->paramOn(rp, Layout(0.50, 0.9, 0.9));
  // components
  RooArgSet * set = ws_->pdf("pdf_ext_total")->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (!pee) {
        size_t found;
        found = name.find("pdf_bs");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_bd");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_comb");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kBlue - 5),   LineStyle(2), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_semi");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kGreen - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_peak");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kCyan - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
      }
      else {
        size_t found;
        found = name.find("pdf_bs");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_bd");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_comb");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE), LineColor(kBlue - 5),   LineStyle(2), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_semi");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE), LineColor(kGreen - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_peak");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *global_data, kFALSE), LineColor(kCyan - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
      }
    }
  }
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  c->Print((get_address("data_fitData_", "pdf_ext_total") + ".gif").c_str());
  c->Print((get_address("data_fitData_", "pdf_ext_total") + ".pdf").c_str());
  delete rp;
  delete c;
}

void pdf_fitData::print_each_channel(string var, string output, RooWorkspace* ws, RooDataSet* rds_) {
  if (ws == 0) ws = (RooWorkspace*)ws_->Clone();
  if (rds_ == 0) rds_ = (RooDataSet*)global_data->Clone();
  cout << red_color_bold << "printing" << default_console_color << endl;

  merge_mass_to_hist();

  std::vector<double> ssbweights;
  std::vector<double> bsnorms;

  for (int i = 0; i < channels; i++) {
  	if (var.find("bdt") != string::npos) var = name("bdt", i);
  	if (years_ == "0" && i > 1) break;
  	if (years_ == "1" && i < 2) continue;
    for (int j = 0; j < bdt_index_max(i); j++) {
      /// all texts
      string cut;
      string title;
      RooArgSet slice_set;
      RooArgSet projw_set;

      ws->cat("etacat")->setIndex(i);
      ws->cat("bdtcat")->setIndex(j);
      ws->cat("allcat")->setIndex(super_index(i, j));

      if (!simul_bdt_ && !simul_all_) {
        slice_set.add(*ws->cat("etacat"));
        projw_set.add(*ws->cat("etacat"));
        cut = Form("etacat==etacat::etacat_%d", i);
        title = get_title(i);
      }
      else if (simul_bdt_ && !simul_all_) {
        slice_set.add(*ws->cat("etacat"));
        projw_set.add(*ws->cat("etacat"));
        cut = Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", i, j);
        title = Form("%s - BDT bin %d", get_title(i).c_str(), j);
//        title = get_title(i);
        slice_set.add(*ws->cat("bdtcat"));
        projw_set.add(*ws->cat("bdtcat"));
      }
      else if (!simul_bdt_ && simul_all_) {
        slice_set.add(*ws->cat("allcat"));
        projw_set.add(*ws->cat("allcat"));
        cut = Form("allcat==allcat::allcat_%d", super_index(i, j));
        title = Form("%s - BDT bin %d", get_title(i).c_str(), j);
//        title = get_title(i);
      }
      int bins = 25;
      if (pee) projw_set.add(*ws->var("ReducedMassRes"));

      ws->var(var.c_str())->SetTitle("m_{#mu#mu}");
      ws->var(var.c_str())->setUnit("GeV");
      RooPlot* final_p = ws->var(var.c_str())->frame(Bins(bins), Title(title.c_str()));
      final_p->SetTitle("");
      rds_->plotOn(final_p, Cut(cut.c_str()), Invisible());
//      ws->pdf(pdfname.c_str())->plotOn(final_p, Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_), LineColor(kBlue), LineWidth(3));
//      if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_semi", i, j)), DrawOption("F"), FillColor(kGreen - 3), FillStyle(3004), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//      ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bs", i, j)), DrawOption("F"), FillColor(kRed),        FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//      ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bd", i, j)), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//      if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_peak", i, j)), DrawOption("F"), FillColor(kBlack), FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//
//
//      if (!berns_) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_comb", i, j)), LineColor(kBlue - 1),   LineStyle(2), LineWidth(3), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//      else ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_1comb", i, j)), LineColor(kBlue - 1),   LineStyle(2), LineWidth(3), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//      if (output == "") {
////        ws->pdf(pdfname.c_str())->plotOn(final_p, VisualizeError(*RFR, 1, 1), FillColor(kYellow), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_), MoveToBack());
//        if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_semi", i, j)), LineColor(kGreen),  LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//        ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bs", i, j)), LineColor(kBlack),        LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//        ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bd", i, j)), LineColor(kBlack), LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//        if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_peak", i, j)), LineColor(kBlack),  LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
//      }

      bool plotpdfs = true;

      RooArgList pdf_list("pdflist");
      pdf_list.add(*ws_->pdf(name("pdf_h_bs", i, j)));
      pdf_list.add(*ws_->pdf(name("pdf_h_bd", i, j)));
      if (berns_) pdf_list.add(*ws_->pdf(name("pdf_h_1comb", i, j)));
      else pdf_list.add(*ws_->pdf(name("pdf_h_comb", i, j)));
      pdf_list.add(*ws_->pdf(name("pdf_h_semi", i, j)));
      pdf_list.add(*ws_->pdf(name("pdf_h_peak", i, j)));
      RooArgList N_list("varlist");
      if (BF_ > 0) N_list.add(*ws_->function(name("N_bs_formula", i, j)));
      else N_list.add(*ws_->var(name("N_bs", i, j)));
      if ((SM_ || bd_constr_) && BF_ < 2) N_list.add(*ws_->function(name("N_bd_constr", i, j)));
      else if (BF_ > 1) N_list.add(*ws_->function(name("N_bd_formula", i, j)));
      else N_list.add(*ws_->var(name("N_bd", i, j)));
      N_list.add(*ws_->var(name("N_comb", i, j)));
      N_list.add(*ws_->var(name("N_semi", i, j)));
      N_list.add(*ws_->var(name("N_peak", i, j)));
      RooAddPdf pdf_ext_sum(name("pdf_h_ext_sum", i, j), "pdf_ext_sum", pdf_list, N_list);

      RooAbsPdf *fullsumpdf = ws->pdf(name("pdf_ext_sum", i, j));
      double fullsumpdfnorm = fullsumpdf->expectedEvents(*ws->var(var.c_str()));
      pdf_ext_sum.plotOn(final_p, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3), NumCPU(16), Name("fullpdf"));
      cout << "full printed" << endl;
      pdf_ext_sum.plotOn(final_p, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3));
      RooAbsPdf *bspdf = ws_->pdf(name("pdf_h_bs", i, j));
      double bsnorm = ws->function(name("N_bs_formula",i,j))->getVal();
      bspdf->plotOn(final_p, Normalization(bsnorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kRed),           FillStyle(3365), NumCPU(16), Name("bs"));
      bspdf->plotOn(final_p, Normalization(bsnorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kRed), LineWidth(2), LineStyle(1), NumCPU(16));
      cout << "bs printed" << endl;
      RooAbsPdf *bdpdf = ws_->pdf(name("pdf_h_bd", i, j));
      double bdnorm = ws->function(name("N_bd_formula",i,j))->getVal();
      bdpdf->plotOn(final_p, Normalization(bdnorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kViolet - 4),    FillStyle(3344), NumCPU(16), Name("bd"));
      bdpdf->plotOn(final_p, Normalization(bdnorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kViolet - 4), LineWidth(2), LineStyle(1), NumCPU(16));
      cout << "bd printed" << endl;
      RooAbsPdf *peakpdf = ws_->pdf(name("pdf_h_peak", i, j));
      double peaknorm = ws->var(name("N_peak",i,j))->getVal();
//      peakpdf->plotOn(final_p, Normalization(peaknorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kBlack),     FillStyle(3001), NumCPU(16));
      peakpdf->plotOn(final_p, Normalization(peaknorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kBlack), LineWidth(3), LineStyle(5), NumCPU(16), Name("peak"));
      cout << "peak printed" << endl;
      RooAbsPdf *semipdf = ws_->pdf(name("pdf_h_semi", i, j));
      double seminorm = ws->var(name("N_semi",i,j))->getVal();
//      semipdf->plotOn(final_p, Normalization(seminorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kGreen - 3), FillStyle(3004), NumCPU(16));
      semipdf->plotOn(final_p, Normalization(seminorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kGreen -3), LineWidth(4), LineStyle(2), NumCPU(16), Name("semi"));
      cout << "semi printed" << endl;
      RooAbsPdf *combpdf = ws_->pdf(name(berns_ ? "pdf_h_1comb" : "pdf_h_comb", i, j));
      double combnorm = ws->var(name("N_comb",i,j))->getVal();
      combpdf->plotOn(final_p, Normalization(combnorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kBlue - 1), LineWidth(3),  LineStyle(2), NumCPU(16), Name("comb"));
      cout << "comb printed" << endl;

      rds_->plotOn(final_p, Cut(cut.c_str()), Invisible());
      TH1D* histo_data = (TH1D*)rds_->createHistogram(name("histo_data",i,j), *ws_->var(var.c_str()), Cut(cut.c_str()), Binning(bins, 4.9, 5.9));
      histo_data->Sumw2(false);
      histo_data->SetBinErrorOption(TH1::kPoisson);
      histo_data->SetMarkerStyle(20);
      histo_data->SetLineColor(kBlack);
      final_p->SetMinimum(0);
      bool legend = true;
      if (legend) {
      	if (simul_all_ &&  i == 1 && j == 0) final_p->SetMaximum(50);
      	if (simul_all_ &&  i == 2 && j == 0) final_p->SetMaximum(290);
      	if (simul_all_ &&  i == 2 && j == 1) final_p->SetMaximum(50);
      	if (simul_all_ &&  i == 2 && j == 2) final_p->SetMaximum(15);
      	if (simul_all_ &&  i == 3 && j == 0) final_p->SetMaximum(170);
      	if (simul_all_ &&  i == 3 && j == 1) final_p->SetMaximum(22);
      	if (!simul_all_ && !simul_bdt_ && !bdt_fit_) final_p->SetMaximum(10);
      }
      TCanvas* final_c = new TCanvas("final_c", "final_c", 600, 600);
      final_p->GetYaxis()->SetTitleOffset(1.);
      final_p->GetXaxis()->SetTitleOffset(1.);
      final_p->GetXaxis()->SetLabelOffset(0.01);
      final_p->GetXaxis()->SetTitleSize(0.043);
      final_p->GetYaxis()->SetTitleSize(0.043);
      final_p->Draw();
      histo_data->Draw("Esame");

      TLatex Tlatex_i;
      Tlatex_i.SetTextFont(42);
      Tlatex_i.SetTextSize(0.035);
      Tlatex_i.SetTextAlign(11);
      Tlatex_i.SetNDC();
      Tlatex_i.DrawLatex(0.1,0.91, title.c_str());
//      if (simul_all_ || simul_bdt_) Tlatex_i.DrawLatex(0.75,0.85, Form("BDT bin %d", j));
      if (legend) {
      	final_c->Update();
      	TLegend *leg1 = new TLegend(0.55,0.62,0.86,0.87);
      	leg1->SetFillColor(kWhite);
      	leg1->SetLineColor(kWhite);
      	TLegendEntry *data_entry = new TLegendEntry(histo_data, "data", "lep");
      	data_entry->SetMarkerStyle(20);
      	leg1->AddEntry(data_entry, "data", "lep");
      	leg1->AddEntry(final_p->findObject("fullpdf"),"full PDF","L");
      	TLegendEntry *bs_entry = new TLegendEntry(final_p->findObject("bs"), "B^{0}_{s}#rightarrow#mu^{+}#mu^{-}", "f");
      	bs_entry->SetLineColor(kRed);
      	bs_entry->SetFillColor(kRed);
      	bs_entry->SetFillStyle(3365);
      	leg1->AddEntry(bs_entry,"B^{0}_{s}#rightarrow#mu^{+}#mu^{-}","f");
      	TLegendEntry *bd_entry = new TLegendEntry(final_p->findObject("bd"), "B^{0}#rightarrow#mu^{+}#mu^{-}", "f");
      	bd_entry->SetLineColor(kViolet - 4);
      	bd_entry->SetFillColor(kViolet - 4);
      	bd_entry->SetFillStyle(3344);
      	leg1->AddEntry(bd_entry,"B^{0}#rightarrow#mu^{+}#mu^{-}","f");
      	leg1->AddEntry(final_p->findObject("comb"),"combinatorial bkg","L");
      	leg1->AddEntry(final_p->findObject("semi"),"semileptonic bkg","L");
      	leg1->AddEntry(final_p->findObject("peak"),"peaking bkg","L");
      	leg1->Draw();
      }


      //compute s/(s+b) weights
      bsnorms.push_back(bsnorm);

      //set mass variable to peak position to compute s/(s+b) in b_s peak
      double meanbsval = ws->var(name("Mean_bs",i,j))->getVal();

      //create projections over reduced mass error to compute average pdf values in the mass dimension
      RooAbsPdf *fullsumpdfproj = pdf_ext_sum.createProjection(*ws->var("ReducedMassRes"));
      RooAbsPdf *bspdfproj = bspdf->createProjection(*ws->var("ReducedMassRes"));

      static_cast<RooRealVar*>(fullsumpdfproj->getObservables(*ws->var(var.c_str()))->find(*ws->var(var.c_str())))->setVal(meanbsval);
      static_cast<RooRealVar*>(bspdfproj->getObservables(*ws->var(var.c_str()))->find(*ws->var(var.c_str())))->setVal(meanbsval);

      double ssbweight = bsnorm*bspdfproj->getVal(*ws->var(var.c_str()))/(fullsumpdfnorm*fullsumpdfproj->getVal(*ws->var(var.c_str())));
      ssbweights.push_back(ssbweight);


      if (output == "" && false) {
        // legend
        RooArgSet* vars =  ws->pdf(pdfname.c_str())->getVariables();
        RooRealVar* N_bs;
        if (BF_ == 0) N_bs = (RooRealVar*)vars->find(pdf_analysis::name("N_bs", i, j));
        else if (BF_ < 3) N_bs = (RooRealVar*)vars->find("BF_bs");
        else if (BF_ == 3) N_bs = (RooRealVar*)vars->find("BF_bsbd");
        RooRealVar* N_bd;
        if (!bd_constr_ && !SM_ && BF_ < 2) {
          N_bd = (RooRealVar*)vars->find(pdf_analysis::name("N_bd", i, j));
        }
        else if (bd_constr_) N_bd = (RooRealVar*)vars->find("Bd_over_Bs");
        else if (BF_ > 1) N_bd = (RooRealVar*)vars->find("BF_bd");
        RooRealVar* N_comb = (RooRealVar*)vars->find(pdf_analysis::name("N_comb", i, j));
        vector <string> fitresult_tex_vec;
        if (BF_ == 0) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << fixed << "N(B_{s}) = " << N_bs->getVal() << " ^{+" << getErrorHigh(N_bs) << "}_{" << getErrorLow(N_bs) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
        }
        else {
          ostringstream fitresult_tex;
          string title("BF(B^{0}_{s})");
          string name("BF_bs");
          if (BF_ == 3) {
            title = "BF(B^{0}_{s})/BF(B^{0})";
            name = "BF_bsbd";
          }
          fitresult_tex << setprecision(2) << scientific << title << " = " << N_bs->getVal() << " ^{+" << getErrorHigh(N_bs) << "}_{" << getErrorLow(N_bs) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
          ostringstream fitresult_tex2;
          fitresult_tex2 << "(N(B^{0}_{s}) = " << setprecision(2) << fixed << ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          Double_t BF_bs_val = ws->var(name.c_str())->getVal();
          Double_t N_bs_ = ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          ws->var(name.c_str())->setVal(BF_bs_val + getErrorHigh(N_bs));
          Double_t N_bs_up = ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          Double_t N_bs_error_up = N_bs_up - N_bs_;
          fitresult_tex2 << " ^{+" << N_bs_error_up;
          ws->var(name.c_str())->setVal(BF_bs_val);
          ws->var(name.c_str())->setVal(BF_bs_val + getErrorLow(N_bs));
          Double_t N_bs_down = ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          Double_t N_bs_error_down = N_bs_ - N_bs_down;
          fitresult_tex2 << "}_{-" << N_bs_error_down << "}" << ")";
          ws->var(name.c_str())->setVal(BF_bs_val);
          fitresult_tex_vec.push_back(fitresult_tex2.str());
        }
        if (!bd_constr_ && !SM_ && BF_ < 2) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << fixed << "N(B_{d}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
        }
        else if (bd_constr_) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << fixed << "N(B_{d}) / N(B_{s}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
        }
        else if (BF_ > 1) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << scientific << "BF(B^{0}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
          ostringstream fitresult_tex2;
          fitresult_tex2 << "(N(B^{0}) = " << setprecision(2) << fixed << ws->function(name("N_bd_formula", i, j))->getVal();
          Double_t BF_bd_val = ws->var("BF_bd")->getVal();
          Double_t N_bd_ = ws->function(name("N_bd_formula", i, j))->getVal();
          ws->var("BF_bd")->setVal(BF_bd_val + getErrorHigh(N_bd));
          Double_t N_bd_up = ws->function(name("N_bd_formula", i, j))->getVal();
          Double_t N_bd_error_up = N_bd_up - N_bd_;
          fitresult_tex2 << " ^{+" << N_bd_error_up;
          ws->var("BF_bd")->setVal(BF_bd_val);
          ws->var("BF_bd")->setVal(BF_bd_val + getErrorLow(N_bd));
          Double_t N_bd_down = ws->function(name("N_bd_formula", i, j))->getVal();
          Double_t N_bd_error_down = N_bd_ - N_bd_down;
          fitresult_tex2 << "}_{-" << N_bd_error_down << "}" << ")";
          ws->var("BF_bd")->setVal(BF_bd_val);
          fitresult_tex_vec.push_back(fitresult_tex2.str());
        }
        ostringstream fitresult_tex;
        fitresult_tex << setprecision(2) << fixed << "N(comb. bkg) = " << N_comb->getVal() << " ^{+" << getErrorHigh(N_comb) << "}_{" << getErrorLow(N_comb) << "}";
        fitresult_tex_vec.push_back(fitresult_tex.str());
        if (BF_ > 0) {
        	RooRealVar* N_semi = (RooRealVar*)vars->find(pdf_analysis::name("N_semi", i, j));
        	RooRealVar* N_peak = (RooRealVar*)vars->find(pdf_analysis::name("N_peak", i, j));
        	ostringstream fitresult_tex2;
        	fitresult_tex2 << setprecision(2) << fixed << "N(semi bkg) = " << N_semi->getVal() << " ^{+" << getErrorHigh(N_semi) << "}_{" << getErrorLow(N_semi) << "}";
        	fitresult_tex_vec.push_back(fitresult_tex2.str());
        	ostringstream fitresult_tex3;
        	fitresult_tex3 << setprecision(2) << fixed << "N(peak bkg) = " << N_peak->getVal() << " ^{+" << getErrorHigh(N_peak) << "}_{" << getErrorLow(N_peak) << "}";
        	fitresult_tex_vec.push_back(fitresult_tex3.str());
        }
        TPaveText* fitresults;
        	if ((simul_all_ || simul_bdt_) && (i == 1 && j == 0 )) fitresults = new TPaveText(0.57, 0.76, 0.89, 0.89, "NDCR");
        	else if ((simul_all_ || simul_bdt_) && (i == 2 && j == 0 )) fitresults = new TPaveText(0.52, 0.20, 0.89, 0.49, "NDCR");
        	else if ((simul_all_ || simul_bdt_) && (i == 2 && j == 1 )) fitresults = new TPaveText(0.57, 0.76, 0.89, 0.89, "NDCR");
        	else if ((simul_all_ || simul_bdt_) && (i == 3 && j == 0 )) fitresults = new TPaveText(0.57, 0.18, 0.89, 0.45, "NDCR");
        	else if ((simul_all_ || simul_bdt_) && (i == 3 && j == 1 )) fitresults = new TPaveText(0.57, 0.76, 0.89, 0.89, "NDCR");
        	else fitresults = new TPaveText(0.57, 0.54, 0.89, 0.89, "NDCR");
        for (unsigned int jj = 0; jj < fitresult_tex_vec.size(); jj++) {
          fitresults->AddText(fitresult_tex_vec[jj].c_str());
        }
        fitresults->SetFillColor(0);
        fitresults->SetShadowColor(0);
        fitresults->SetTextSize(0.03);
        if (simul_all_ || simul_bdt_) {
        	if (i == 1 && j == 0 ) fitresults->SetTextSize(0.018);
        	if (i == 2 && j == 1 ) fitresults->SetTextSize(0.023);
        	if (i == 3 && j == 0 ) fitresults->SetTextSize(0.028);
        	if (i == 3 && j == 1 ) fitresults->SetTextSize(0.022);
        }
        fitresults->SetTextAlign(11);
        fitresults->SetLineColor(0);
        fitresults->Draw();
      }
      channel = i;
      channel_bdt = j;
      string output_name;
      if (output == "") output_name = get_address("data_" + var, "", true);
      else output_name = get_address(pdfname, output, true);
      final_c->Print( (output_name + ".gif").c_str() );
      final_c->Print( (output_name + ".pdf").c_str() );
      delete final_p;
      delete final_c;

      if (var.find("bdt") != string::npos) {
      	TH1* mass_bdt_h = ws_->pdf(name("pdf_ext_total", channel))->createHistogram("fit", *ws_->var(var.c_str()), Binning(50), YVar(*ws_->var("Mass"), Binning(30))) ;
      	mass_bdt_h->SetLineColor(kBlue) ;
      	mass_bdt_h->GetXaxis()->SetTitleOffset(2.) ;
      	mass_bdt_h->GetYaxis()->SetTitleOffset(2.) ;
      	mass_bdt_h->GetZaxis()->SetTitleOffset(2.5) ;
      	TCanvas* c_surf = new TCanvas("c_surf", "c_surf", 600, 600);
      	c_surf->SetTheta(33.9);
      	c_surf->SetPhi(245.3);
      	mass_bdt_h->Draw("surf");
      	c_surf->Print( (get_address("MassBDT", pdf_name) + ".gif").c_str());
      	c_surf->Print( (get_address("MassBDT", pdf_name) + ".pdf").c_str());
//      	c_surf->Print( (get_address("MassBDT", pdf_name) + ".C").c_str());
//      	c_surf->Print( (get_address("MassBDT", pdf_name) + ".root").c_str());
      	delete c_surf;
      }
      cout << "channel " << i << " " << j << "  entries =" << rds_->sumEntries(cut.c_str()) << endl;
    }
  }


  TLatex Tlatex;
  Tlatex.SetTextFont(42);
  Tlatex.SetTextSize(0.031);
  Tlatex.SetTextAlign(11);
  Tlatex.SetNDC();
  string latex_s("CMS - L = 5 fb^{-1} #sqrt{s} = 7 TeV, L = 20 fb^{-1} #sqrt{s} = 8 TeV");

  //compute normalization factor for s/(s+b) weights to conserve fitted b_s normalization
  double nbstotal = 0;
  double nbstotalweighted = 0;
  for (unsigned int i=0; i<ssbweights.size(); ++i) {
    nbstotal += bsnorms[i];
    nbstotalweighted += ssbweights[i]*bsnorms[i];
    printf("i = %i, ssbweight = %5f\n",i,ssbweights[i]);
  }
  double weightscale = nbstotal/nbstotalweighted;

  for (unsigned int i=0; i<ssbweights.size(); ++i) {
    printf("i = %i, ssbweightnorm = %5f\n",i,weightscale*ssbweights[i]);
  }

  RooArgList totalpdflist;
  RooArgList bspdflist;
  RooArgList bdpdflist;
  RooArgList peakpdflist;
  RooArgList semipdflist;
  RooArgList combpdflist;

  RooArgList totalcoeflist;
  RooArgList bscoeflist;
  RooArgList bdcoeflist;
  RooArgList peakcoeflist;
  RooArgList semicoeflist;
  RooArgList combcoeflist;

  RooRealVar weightvar("weightvar","",1.);

  RooArgSet dvars(*rds_->get());
  dvars.add(weightvar);

  RooDataSet wdata("wdata","",dvars,WeightVar(weightvar));

  //build lists of pdfs and normalizations for weighted sum, as well as weighted dataset
  int icat = 0;
  for (int i = 0; i < channels; i++) {
  	if (var.find("bdt") != string::npos) var = name("bdt", i);
  	if (years_ == "0" && i > 1) break;
  	if (years_ == "1" && i < 2) continue;
    for (int j = 0; j < bdt_index_max(i); j++) {

      string cut;
      if (!simul_bdt_ && !simul_all_) {
        cut = Form("etacat==etacat::etacat_%d", i);
      }
      else if (simul_bdt_ && !simul_all_) {
        cut = Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", i, j);
      }
      else if (!simul_bdt_ && simul_all_) {
        cut = Form("allcat==allcat::allcat_%d", super_index(i, j));
      }

      double weight = ssbweights[icat]*weightscale;

      //fill weighted data
      RooDataSet *dcat = (RooDataSet*)rds_->reduce(cut.c_str());
      for (int iev=0; iev<dcat->numEntries(); ++iev) {
      	const RooArgSet *dset = dcat->get(iev);
      	wdata.add(*dset,weight);
      }

      bspdflist.add(*ws_->pdf(name("pdf_h_bs", i, j)));
      bdpdflist.add(*ws_->pdf(name("pdf_h_bd", i, j)));
      if (berns_) combpdflist.add(*ws_->pdf(name("pdf_h_1comb", i, j)));
      else combpdflist.add(*ws_->pdf(name("pdf_h_comb", i, j)));
      semipdflist.add(*ws_->pdf(name("pdf_h_semi", i, j)));
      peakpdflist.add(*ws_->pdf(name("pdf_h_peak", i, j)));

      if (BF_ > 0) bscoeflist.add(RooConst(weight*ws_->function(name("N_bs_formula", i, j))->getVal()));
      else bscoeflist.add(RooConst(weight*ws_->var(name("N_bs", i, j))->getVal()));
      if ((SM_ || bd_constr_) && BF_ < 2) bdcoeflist.add(RooConst(weight*ws_->function(name("N_bd_constr", i, j))->getVal()));
      else if (BF_ > 1) bdcoeflist.add(RooConst(weight*ws_->function(name("N_bd_formula", i, j))->getVal()));
      else bdcoeflist.add(RooConst(weight*ws_->var(name("N_bd", i, j))->getVal()));
      combcoeflist.add(RooConst(weight*ws_->var(name("N_comb", i, j))->getVal()));
      semicoeflist.add(RooConst(weight*ws_->var(name("N_semi", i, j))->getVal()));
      peakcoeflist.add(RooConst(weight*ws_->var(name("N_peak", i, j))->getVal()));

      ++icat;
    }
  }

  totalpdflist.add(bspdflist);
  totalpdflist.add(bdpdflist);
  totalpdflist.add(peakpdflist);
  totalpdflist.add(semipdflist);
  totalpdflist.add(combpdflist);

  totalcoeflist.add(bscoeflist);
  totalcoeflist.add(bdcoeflist);
  totalcoeflist.add(peakcoeflist);
  totalcoeflist.add(semicoeflist);
  totalcoeflist.add(combcoeflist);

  RooAddPdf wpdfsum("wpdfsum", "", totalpdflist, totalcoeflist);
  RooAddPdf wpdfbs("wpdfbs", "", bspdflist, bscoeflist);
  RooAddPdf wpdfbd("wpdfbd", "", bdpdflist, bdcoeflist);
  RooAddPdf wpdfpeak("wpdfpeak", "", peakpdflist, peakcoeflist);
  RooAddPdf wpdfsemi("wpdfsemi", "", semipdflist, semicoeflist);
  RooAddPdf wpdfcomb("wpdfcomb", "", combpdflist, combcoeflist);


  //sum of b_s and b_d pdfs
  RooArgList totalsigpdflist;
  totalsigpdflist.add(bspdflist);
  totalsigpdflist.add(bdpdflist);

  RooArgList totalsigcoeflist;
  totalsigcoeflist.add(bscoeflist);
  totalsigcoeflist.add(bdcoeflist);

  RooAddPdf wsigsum("wsigsum","",totalsigpdflist,totalsigcoeflist);

  //sum of bkg pdfs
  RooArgList totalbkgpdflist;
  totalbkgpdflist.add(peakpdflist);
  totalbkgpdflist.add(semipdflist);
  totalbkgpdflist.add(combpdflist);

  RooArgList totalbkgcoeflist;
  totalbkgcoeflist.add(peakcoeflist);
  totalbkgcoeflist.add(semicoeflist);
  totalbkgcoeflist.add(combcoeflist);

  RooAddPdf wbkgsum("wbkgsum","",totalbkgpdflist,totalbkgcoeflist);
  double fullbkgnorm = wbkgsum.expectedEvents(*ws->var(var.c_str()));

  bool plotpdfs = true;

  //final weighted plot
  int bins = 25;

  //histogram of weigthted data
  TH1 *wdatahist = wdata.createHistogram("wdatahist",*ws->var(var.c_str()),Binning(bins,ws->var(var.c_str())->getMin(),ws->var(var.c_str())->getMax()));

  //histgram of sum of bkg pdfs
  TH1 *wbkghist = wbkgsum.createHistogram("wbkghist",*ws->var(var.c_str()),Binning(bins,ws->var(var.c_str())->getMin(),ws->var(var.c_str())->getMax()));
  wbkghist->Scale(fullbkgnorm/wbkghist->GetSumOfWeights());

//   TCanvas *csubtest = new TCanvas;
//   //wdatahist->Draw("E");
//   wbkghist->Draw("Hist");
//   csubtest->SaveAs("csubtest.eps");

 // return;

  //subtract background contribution without changing the error on each bin
  for (int ibin=1; ibin<(wdatahist->GetNbinsX()+1); ++ibin) {
    wdatahist->SetBinContent(ibin, wdatahist->GetBinContent(ibin) - wbkghist->GetBinContent(ibin));
  }

  RooDataHist wdatasub("wdatasub","",*ws->var(var.c_str()),wdatahist);

  RooPlot* final_p = ws->var(var.c_str())->frame(Bins(bins));

  wdata.plotOn(final_p, Invisible(), Name("data"));
  final_p->SetTitle("");
  double fullsumpdfnorm = wpdfsum.expectedEvents(*ws->var(var.c_str()));
  if (plotpdfs) wpdfsum.plotOn(final_p, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3), NumCPU(16), Name("fullpdf"));
  cout << "full printed" << endl;

  double bsnorm = wpdfbs.expectedEvents(*ws->var(var.c_str()));
  if (plotpdfs) wpdfbs.plotOn(final_p, Normalization(bsnorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kRed),           FillStyle(3001), NumCPU(16), Name("bs"));
  cout << "bs printed" << endl;

  double bdnorm = wpdfbd.expectedEvents(*ws->var(var.c_str()));
  if (plotpdfs) wpdfbd.plotOn(final_p, Normalization(bdnorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kViolet - 4),    FillStyle(3001), NumCPU(16), Name("bd"));
  cout << "bd printed" << endl;

  double peaknorm = wpdfpeak.expectedEvents(*ws->var(var.c_str()));
//  if (plotpdfs) wpdfpeak.plotOn(final_p, Normalization(peaknorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kBlack),     FillStyle(3001), NumCPU(16));
  if (plotpdfs) wpdfpeak.plotOn(final_p, Normalization(peaknorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kBlack), LineWidth(3), LineStyle(5), NumCPU(16), Name("peak"));
  cout << "peak printed" << endl;

  double seminorm = wpdfsemi.expectedEvents(*ws->var(var.c_str()));
//  if (plotpdfs) wpdfsemi.plotOn(final_p, Normalization(seminorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kGreen - 3), FillStyle(3004), NumCPU(16));
  if (plotpdfs) wpdfsemi.plotOn(final_p, Normalization(seminorm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kGreen -3), LineWidth(4), LineStyle(2), NumCPU(16), Name("semi"));
  cout << "semi printed" << endl;

  double combnorm = wpdfcomb.expectedEvents(*ws->var(var.c_str()));
  if (plotpdfs) wpdfcomb.plotOn(final_p, Normalization(combnorm, RooAbsReal::NumEvent), LineColor(kBlue - 1), LineWidth(3),  LineStyle(2), NumCPU(16), Name("comb"));
  cout << "comb printed" << endl;

  wdata.plotOn(final_p, DataError(RooAbsData::SumW2));

  final_p->GetYaxis()->SetTitle("S/(S+B) Weighted Events / ( 0.04 GeV)");
  final_p->GetYaxis()->SetTitleOffset(1.);
  final_p->GetXaxis()->SetTitleOffset(1.);
  final_p->GetXaxis()->SetLabelOffset(0.01);
  final_p->GetXaxis()->SetTitleSize(0.043);
  final_p->GetYaxis()->SetTitleSize(0.043);

  TCanvas* final_c = new TCanvas("final_c", "final_c", 600, 600);
  final_p->SetMinimum(0);
  final_p->Draw();
  Tlatex.DrawLatex(0.1,0.91, latex_s.c_str());
  final_c->Update();

  if (plotpdfs) {
  	TLegend *leg1 = new TLegend(0.60,0.63,0.88,0.87);
  	leg1->SetFillColor(kWhite);
  	leg1->SetLineColor(kWhite);
  	TLegendEntry *data_entry = new TLegendEntry(final_p->findObject("data"), "data", "lep");
  	data_entry->SetMarkerStyle(20);
  	leg1->AddEntry(data_entry, "data", "lep");
  	leg1->AddEntry(final_p->findObject("fullpdf"),"full PDF","L");
  	TLegendEntry *bs_entry = new TLegendEntry(final_p->findObject("bs"), "B_{s}#rightarrow#mu^{+}#mu^{-}", "f");
  	bs_entry->SetLineColor(kRed);
  	bs_entry->SetFillColor(kRed);
  	bs_entry->SetFillStyle(3001);
  	leg1->AddEntry(bs_entry,"B_{s}#rightarrow#mu^{+}#mu^{-}","f");
  	TLegendEntry *bd_entry = new TLegendEntry(final_p->findObject("bd"), "B_{d}#rightarrow#mu^{+}#mu^{-}", "f");
  	bd_entry->SetLineColor(kViolet - 4);
  	bd_entry->SetFillColor(kViolet - 4);
  	bd_entry->SetFillStyle(3001);
  	leg1->AddEntry(bd_entry,"B_{d}#rightarrow#mu^{+}#mu^{-}","f");
  	leg1->AddEntry(final_p->findObject("comb"),"combinatorial bkg","L");
  	leg1->AddEntry(final_p->findObject("semi"),"semileptonic bkg","L");
  	leg1->AddEntry(final_p->findObject("peak"),"peaking bkg","L");
  	leg1->Draw();
  }

  string output_name = "fig/data_weighted";
  final_c->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".gif").c_str() );
  final_c->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".pdf").c_str() );
  delete final_c;
  delete final_p;

  //background-subtracted plot

  RooPlot* final_psub = ws->var(var.c_str())->frame(Bins(bins));

  //wdatasub.plotOn(final_psub, Invisible(), Name("data"), DataError(RooAbsData::SumW2));
  final_psub->SetTitle("");
  double sigsumpdfnorm = wsigsum.expectedEvents(*ws->var(var.c_str()));
  if (plotpdfs) wsigsum.plotOn(final_psub, Normalization(sigsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3), NumCPU(16), Name("sigpdf"));
  cout << "full printed" << endl;

  if (plotpdfs) wpdfbs.plotOn(final_psub, Normalization(bsnorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kRed),           FillStyle(3001), NumCPU(16), Name("bs"));
  cout << "bs printed" << endl;

  if (plotpdfs) wpdfbd.plotOn(final_psub, Normalization(bdnorm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kViolet - 4),    FillStyle(3001), NumCPU(16), Name("bd"));
  cout << "bd printed" << endl;


  wdatasub.plotOn(final_psub, DataError(RooAbsData::SumW2));

  final_psub->GetYaxis()->SetTitle("S/(S+B) Weighted Events / ( 0.04 GeV)");
  final_psub->GetYaxis()->SetTitleOffset(1.2);
  final_psub->GetXaxis()->SetTitleOffset(1.2);
  final_psub->GetXaxis()->SetLabelOffset(0.01);

  TCanvas* final_csub = new TCanvas("final_csub", "final_csub", 600, 600);
  final_psub->SetMinimum(0);
  if (wdatahist->GetMinimum()) {
    final_psub->SetMinimum(1.5*wdatahist->GetMinimum());
  }
  final_psub->Draw();
  //wdatahist->Draw("ESAME");
  Tlatex.DrawLatex(0.1,0.91, latex_s.c_str());
  final_csub->Update();

  string output_name_sub = "fig/data_weighted_sub";
  final_csub->Print( (output_name_sub + (simul_all_? "_simulBdt" : "") + ".gif").c_str() );
  final_csub->Print( (output_name_sub + (simul_all_? "_simulBdt" : "") + ".pdf").c_str() );
  delete final_csub;
  delete final_psub;


  RooArgList siglist(bscoeflist);
  siglist.add(bdcoeflist);
  RooArgList bkglist(peakcoeflist);
  bkglist.add(semicoeflist);
  bkglist.add(combcoeflist);

  RooArgList sigpdflist(bspdflist);
  sigpdflist.add(bdpdflist);
  RooArgList bkgpdflist(peakpdflist);
  bkgpdflist.add(semipdflist);
  bkgpdflist.add(combpdflist);
  RooAddPdf wpdfsig("wpdfsig", "", sigpdflist, siglist);
  RooAddPdf wpdfbkg("wpdfbkg", "", bkgpdflist, bkglist);

  RooArgList totlist(siglist);
  totlist.add(bkglist);
  RooArgList totpdflist(sigpdflist);
  totpdflist.add(bkgpdflist);

  RooAddPdf wpdftot("wpdftot", "", totpdflist, totlist);

  RooPlot* final_pp = ws->var(var.c_str())->frame(Bins(bins));

  bool plotpdfs2 = true;
  wdata.plotOn(final_pp, Invisible(), Name("data"));
  final_pp->SetTitle("");
  if (plotpdfs2) wpdftot.plotOn(final_pp, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), FillColor(kYellow), DrawOption("F"), NumCPU(16), Name("sigpdf"));
  if (plotpdfs2) wpdftot.plotOn(final_pp, Components(bkgpdflist), FillColor(kWhite), DrawOption("F"), NumCPU(16));
  if (plotpdfs2) wpdftot.plotOn(final_pp, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3), NumCPU(16), Name("fullpdf"));
  if (plotpdfs2) wpdftot.plotOn(final_pp, Components(bkgpdflist), LineStyle(2), LineColor(kBlue), NumCPU(16), Name("bkgpdf"));
  cout << "full printed" << endl;
  wdata.plotOn(final_pp, DataError(RooAbsData::SumW2));

  final_pp->GetYaxis()->SetTitle("S/(S+B) Weighted Events / ( 0.04 GeV)");
  final_pp->GetYaxis()->SetTitleOffset(1.);
  final_pp->GetXaxis()->SetTitleOffset(1.);
  final_pp->GetXaxis()->SetLabelOffset(0.01);
  final_pp->GetXaxis()->SetTitleSize(0.043);
  final_pp->GetYaxis()->SetTitleSize(0.043);

  TCanvas* final_cc = new TCanvas("final_cc", "final_cc", 600, 600);
  final_pp->SetMinimum(0);
  final_pp->Draw();
  Tlatex.DrawLatex(0.1,0.91, latex_s.c_str());
  final_cc->Update();
  if (plotpdfs2) {
    TLegend *leg2 = new TLegend(0.60,0.63,0.88,0.87);
    leg2->SetFillColor(kWhite);
    leg2->SetLineColor(kWhite);
    TLegendEntry *data2_entry = new TLegendEntry(final_pp->findObject("data"), "data", "lep");
    data2_entry->SetMarkerStyle(20);
    leg2->AddEntry(data2_entry, "data", "lep");
    leg2->AddEntry(final_pp->findObject("fullpdf"),"full PDF","L");
    TLegendEntry *signal_entry = new TLegendEntry(final_pp->findObject("sigpdf"), "signals", "f");
    signal_entry->SetLineColor(kYellow);
    signal_entry->SetFillColor(kYellow);
    signal_entry->SetFillStyle(1001);
    leg2->AddEntry(signal_entry,"signals","f");
    TLegendEntry *bkg_entry = new TLegendEntry(final_pp->findObject("bkgpdf"), "bkgs", "l");
    bkg_entry->SetLineColor(kBlue);
    bkg_entry->SetLineStyle(2);
    bkg_entry->SetLineWidth(3);
    leg2->AddEntry(bkg_entry, "bkgs","l");
    leg2->Draw();
  }
  output_name = "fig/data_weighted2";
  final_cc->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".gif").c_str() );
  final_cc->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".pdf").c_str() );
  delete final_cc;
  delete final_pp;

  RooPlot* final_ppp = ws->var(var.c_str())->frame(Bins(bins));
  bool plotpdfs3 = true;
    wdata.plotOn(final_ppp, Invisible(), Name("data"));
    final_ppp->SetTitle("");
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), FillColor(kYellow), DrawOption("F"), NumCPU(16), Name("sigpdf"));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(bkgpdflist), FillColor(kWhite), DrawOption("F"), NumCPU(16));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlack), LineWidth(3), NumCPU(16), Name("fullpdf"));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(bkgpdflist), LineStyle(2), LineColor(kBlue), NumCPU(16));

    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(bkgpdflist), FillStyle(1001), DrawOption("F"), FillColor(kOrange), NumCPU(16), Name("peakpdf"));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(bkgpdflist), LineColor(kBlack), DrawOption("L"), FillColor(kBlack), LineWidth(2), LineStyle(2), NumCPU(16));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(RooArgSet(combpdflist, semipdflist)), FillStyle(1001), DrawOption("F"), FillColor(kOrange - 9), NumCPU(16), Name("semipdf"));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(RooArgSet(combpdflist, semipdflist)), LineColor(kBlack), LineStyle(2), DrawOption("L"), LineWidth(2), FillColor(kGreen), NumCPU(16));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(combpdflist), FillStyle(1001), DrawOption("F"), FillColor(kYellow - 10), NumCPU(16), Name("combpdf"));
    if (plotpdfs3) wpdftot.plotOn(final_ppp, Components(combpdflist), LineColor(kBlack), LineStyle(2), DrawOption("L"), FillColor(kBlack), LineWidth(2), NumCPU(16));
    cout << "full printed" << endl;
    wdata.plotOn(final_ppp, DataError(RooAbsData::SumW2));

    final_ppp->GetYaxis()->SetTitle("S/(S+B) Weighted Events / ( 0.04 GeV)");
    final_ppp->GetYaxis()->SetTitleOffset(1.);
    final_ppp->GetXaxis()->SetTitleOffset(1.);
    final_ppp->GetXaxis()->SetLabelOffset(0.01);
    final_ppp->GetXaxis()->SetTitleSize(0.043);
    final_ppp->GetYaxis()->SetTitleSize(0.043);

    TCanvas* final_ccc = new TCanvas("final_cc", "final_cc", 600, 600);
    final_ppp->SetMinimum(0);
    final_ppp->Draw();
    Tlatex.DrawLatex(0.1,0.91, latex_s.c_str());
    final_ccc->Update();
    TLegend *leg3 = new TLegend(0.60,0.63,0.88,0.87);
    leg3->SetFillColor(kWhite);
    leg3->SetLineColor(kWhite);
    TLegendEntry *data3_entry = new TLegendEntry(final_ppp->findObject("data"), "data", "lep");
    data3_entry->SetMarkerStyle(20);
    leg3->AddEntry(data3_entry, "data", "lep");
    leg3->AddEntry(final_ppp->findObject("fullpdf"),"full PDF","L");
    TLegendEntry *signal3_entry = new TLegendEntry(final_ppp->findObject("sigpdf"), "signals", "f");
    signal3_entry->SetLineColor(kBlack);
    signal3_entry->SetFillColor(kYellow);
    signal3_entry->SetFillStyle(1001);
    leg3->AddEntry(signal3_entry,"signals","f");
    TLegendEntry *peak_entry = new TLegendEntry(final_ppp->findObject("peakpdf"), "peak", "f");
    peak_entry->SetLineColor(kBlack);
    peak_entry->SetFillColor(kOrange);
    peak_entry->SetFillStyle(1001);
    leg3->AddEntry(peak_entry, "peaking bkg","f");
    TLegendEntry *semi_entry = new TLegendEntry(final_ppp->findObject("semipdf"), "semi", "f");
    semi_entry->SetLineColor(kBlack);
    semi_entry->SetFillColor(kOrange - 9);
    semi_entry->SetFillStyle(1001);
    leg3->AddEntry(semi_entry, "semileptonic bkg","f");
    TLegendEntry *comb_entry = new TLegendEntry(final_ppp->findObject("combpdf"), "comb", "f");
    comb_entry->SetLineColor(kBlack);
    comb_entry->SetFillColor(kYellow - 10);
    comb_entry->SetFillStyle(1001);
    leg3->AddEntry(comb_entry, "combinatorial bkg","f");
    leg3->Draw();

    output_name = "fig/data_weighted3";
    final_ccc->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".gif").c_str() );
    final_ccc->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".pdf").c_str() );
    delete final_ccc;
    delete final_ppp;

    ///4

    RooPlot* final_pppp = ws->var(var.c_str())->frame(Bins(bins));
    bool plotpdfs4 = true;
      wdata.plotOn(final_pppp, Invisible(), Name("data"));
      final_pppp->SetTitle("");
      if (plotpdfs4) wpdftot.plotOn(final_pppp, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), FillStyle(3001), FillColor(kRed), DrawOption("F"), NumCPU(16), Name("bspdf"));
      if (plotpdfs4) wpdftot.plotOn(final_pppp, Components(RooArgSet(bkgpdflist, bdpdflist)), FillStyle(3001), DrawOption("F"), FillColor(kViolet -4), NumCPU(16), Name("bdpdf"));
      if (plotpdfs4) wpdftot.plotOn(final_pppp, Components(RooArgSet(bkgpdflist, bdpdflist)), LineWidth(2), DrawOption("L"), FillColor(kViolet -4), NumCPU(16));
      if (plotpdfs4) wpdftot.plotOn(final_pppp, Components(bkgpdflist), FillColor(kWhite), DrawOption("F"), NumCPU(16));
      if (plotpdfs4) wpdftot.plotOn(final_pppp, Normalization(fullsumpdfnorm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3), NumCPU(16), Name("fullpdf"));
      if (plotpdfs4) wpdftot.plotOn(final_pppp, Components(bkgpdflist), LineStyle(2), LineColor(kBlue), NumCPU(16), Name("bkgpdf"));

      cout << "full printed" << endl;
      wdata.plotOn(final_pppp, DataError(RooAbsData::SumW2));

      final_pppp->GetYaxis()->SetTitle("S/(S+B) Weighted Events / ( 0.04 GeV)");
      final_pppp->GetYaxis()->SetTitleOffset(1.);
      final_pppp->GetXaxis()->SetTitleOffset(1.);
      final_pppp->GetXaxis()->SetLabelOffset(0.01);
      final_pppp->GetXaxis()->SetTitleSize(0.043);
      final_pppp->GetYaxis()->SetTitleSize(0.043);

      TCanvas* final_cccc = new TCanvas("final_cc", "final_cc", 600, 600);
      final_pppp->SetMinimum(0);
      final_pppp->Draw();
      Tlatex.DrawLatex(0.1,0.91, latex_s.c_str());
      final_cccc->Update();
      TLegend *leg4 = new TLegend(0.60,0.63,0.88,0.87);
      leg4->SetFillColor(kWhite);
      leg4->SetLineColor(kWhite);
      TLegendEntry *data4_entry = new TLegendEntry(final_pppp->findObject("data"), "data", "lep");
      data4_entry->SetMarkerStyle(20);
      leg4->AddEntry(data4_entry, "data", "lep");
      leg4->AddEntry(final_pppp->findObject("fullpdf"),"full PDF","L");
      TLegendEntry *bs4_entry = new TLegendEntry(final_pppp->findObject("bspdf"), "signals", "f");
      bs4_entry->SetLineColor(kRed);
      bs4_entry->SetFillColor(kRed);
      bs4_entry->SetFillStyle(3001);
      leg4->AddEntry(bs4_entry,"B_{s}#rightarrow#mu^{+}#mu^{-}","f");
      TLegendEntry *bd4_entry = new TLegendEntry(final_pppp->findObject("bdpdf"), "signals", "f");
      bd4_entry->SetLineColor(kViolet -4);
      bd4_entry->SetFillColor(kViolet -4);
      bd4_entry->SetFillStyle(3001);
      leg4->AddEntry(bd4_entry,"B_{s}#rightarrow#mu^{+}#mu^{-}","f");
      TLegendEntry *bkg_entry = new TLegendEntry(final_pppp->findObject("bkgpdf"), "bkgs", "f");
      bkg_entry->SetLineColor(kBlue);
      bkg_entry->SetLineStyle(2);
      bkg_entry->SetLineWidth(3);
      leg4->AddEntry(bkg_entry, "bkgs","l");
      leg4->Draw();

      output_name = "fig/data_weighted4";
      final_cccc->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".gif").c_str() );
      final_cccc->Print( (output_name + (simul_all_? "_simulBdt" : "") + ".pdf").c_str() );
      delete final_cccc;
      delete final_pppp;


}

void pdf_fitData::fill_dataset(RooDataSet* dataset, bool cut_b, vector <double> cut_, double bdt_cut, TTree* tree, int offset) {
  int events = 0;
  if (!strcmp(tree->GetName(), "SgData_bdt")) {
    TTree* reduced_tree = tree;
    Double_t m1eta_t, m2eta_t, m_t, eta_B_t, bdt_t, me_t, rme_t;
    Int_t evt_t;
    Bool_t muid_t;
    reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
    reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
    reduced_tree->SetBranchAddress("m", &m_t);
    reduced_tree->SetBranchAddress("me", &me_t);
    reduced_tree->SetBranchAddress("eta", &eta_B_t);
    reduced_tree->SetBranchAddress("bdt", &bdt_t);
    reduced_tree->SetBranchAddress("evt", &evt_t);
    reduced_tree->SetBranchAddress("muid",  &muid_t);
    for (int i = 0; i < reduced_tree->GetEntries(); i++) {
      reduced_tree->GetEntry(i);
      bool barrel = fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4;
      double mass = m_t;
      //mass scale: shifting data to MC pdfs
      if (barrel) mass = mass + 0.006;
      else mass = mass + 0.007;
      if (mass > 4.9 && mass < 5.9 && muid_t) {
      	if (me_t < 0.0 || me_t > 0.2) continue; //skip wrong mass scale
      	if (bdt_t < bdt_cut) continue;
        events++;
        Mass->setVal(mass);
        bdt->setVal(bdt_t);
        MassRes->setVal(me_t);
        Double_t rme_t = me_t/mass;
        ReducedMassRes->setVal(rme_t);
        Int_t eta_channel = -1;
        if (barrel) {
          eta_channel = 0 + offset*2;
          if (cut_b && bdt_t < cut_[eta_channel]) continue;
          channels_cat->setIndex(eta_channel);
        }
        else {
          eta_channel = 1 + offset*2;
          if (cut_b && bdt_t < cut_[eta_channel]) continue;
          channels_cat->setIndex(eta_channel);
        }
        int bdt_channel = bdt_index(eta_channel, bdt_t);
        if (simul_bdt_ || simul_all_) {
          if (bdt_channel == -1) continue; /// bdt < 0.1
          bdt_cat->setIndex(bdt_channel);
        }
        if (simul_all_) all_cat->setIndex(super_index(eta_channel, bdt_channel));
      	bdt_0->setVal(bdt_t);
      	bdt_1->setVal(bdt_t);
      	bdt_2->setVal(bdt_t);
      	bdt_3->setVal(bdt_t);
        RooArgSet varlist_tmp(*Mass, *ReducedMassRes, *channels_cat);
        if (bdt_fit_) {
        	varlist_tmp.add(*bdt);
        	varlist_tmp.add(*bdt_0);
        	varlist_tmp.add(*bdt_1);
        	varlist_tmp.add(*bdt_2);
        	varlist_tmp.add(*bdt_3);
        }
        if (simul_bdt_ || simul_all_) varlist_tmp.add(*bdt_cat);
        if (simul_all_) varlist_tmp.add(*all_cat);
        dataset->add(varlist_tmp);
      }
    }
  }
  else {
    cout << "tree name is not SgData_bdt" << endl;
    exit(EXIT_FAILURE);
  }
  cout << "total events = " << events << endl;
}

void pdf_fitData::define_dataset() {
  cout << red_color_bold << "defining dataset" << default_console_color << endl;
  RooArgSet varlist(*Mass, *ReducedMassRes, *channels_cat);
//	varlist.add(*bdt);
        	varlist.add(*bdt_0);
        	varlist.add(*bdt_1);
        	varlist.add(*bdt_2);
        	varlist.add(*bdt_3);
  if (simul_bdt_ || simul_all_) varlist.add(*bdt_cat);
  if (simul_all_) varlist.add(*all_cat);
  global_data = new RooDataSet("global_data", "global_data", varlist);
}

void pdf_fitData::make_dataset(bool cut_b, vector <double> cut_, double bdt_cut, TTree* tree, int offset) {
  cout << red_color_bold << "making dataset" << default_console_color << endl;

  if (!random) fill_dataset(global_data, cut_b, cut_, bdt_cut, tree, offset);

  else {
    if (syst && randomsyst) randomize_constraints(ws_);

    if (!simul_ || (simul_ && !simul_bdt_ && !simul_all_)) {
    	RooArgSet vars(*ws_->var("Mass"), *ws_->var("ReducedMassRes"), *ws_->cat("etacat"));
    	if (bdt_fit_) {
    		for (int i = 0; i < channels; i++) {
    			if (years_ == "0" && i > 1) break;
    			if (years_ == "1" && i < 2) continue;
    			vars.add(*ws_->var(name("bdt", i)));
    		}
    	}
      global_data = ws_->pdf(pdfname.c_str())->generate(vars, Extended());
    }
    else if (simul_ && (simul_bdt_ || simul_all_)) {
      RooArgSet set(*ws_->var("Mass"), *ws_->var("ReducedMassRes"), *ws_->cat("etacat"), *ws_->cat("bdtcat"));
      if (simul_all_) set.add(*ws_->cat("allcat"));
      global_data = new RooDataSet("global_data", "global_data", set);
      for (int i = 0; i < channels; i++) {
      	if (years_ == "0" && i > 1) break;
      	if (years_ == "1" && i < 2) continue;
        for (int j = 0; j < bdt_index_max(i); j++) {
          RooDataSet* data_i = ws_->pdf(name("pdf_ext_total", i, j))->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("ReducedMassRes")), Extended());
          channels_cat->setIndex(i);
          bdt_cat->setIndex(j);
          data_i->addColumn(*channels_cat);
          data_i->addColumn(*bdt_cat);
          if (simul_all_) {
            all_cat->setIndex(super_index(i, j));
            data_i->addColumn(*all_cat);
          }
          global_data->append(*data_i);
        }
      }
    }
  }
  global_data->SetName("global_data");
  cout << " entries = " <<  global_data->sumEntries() << endl;
}

void pdf_fitData::make_pdf_input(string root_s) {
  ws_file_input = new TFile(root_s.c_str());
  if (!ws_file_input) {cout << root_s.c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
  ws_input = (RooWorkspace*)ws_file_input->Get("ws");
  if (!ws_input) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
  cout << "ws file: " << root_s << endl;
  set_ws(ws_input);
}

void pdf_fitData::set_starting_N() {
  cout << red_color_bold << "making pdf" << default_console_color << endl;
  if (random || 1) {
    if (simul_ && !simul_bdt_ && !simul_all_) {
      for (int i = 0; i < channels; i++) {
      	ws_->var(name("N_bu", i))->setVal(N_bu_val[i][0]);
      	ws_->var(name("N_bu", i))->setError(N_bu_err[i][0]);
      	ws_->var(name("N_bs", i))->setVal(estimate_bs[i]);
        ws_->var(name("N_bd", i))->setVal(estimate_bd[i]);
        ws_->var(name("N_peak", i))->setVal(estimate_peak[i]);
        ws_->var(name("N_semi", i))->setVal(estimate_semi[i]);
        ws_->var(name("N_comb", i))->setVal(estimate_comb[i]);
      }
    }
    else if (!simul_) {
    	ws_->var("N_bu")->setVal(N_bu_val[ch_i_][0]);
    	ws_->var("N_bu")->setError(N_bu_err[ch_i_][0]);
      ws_->var("N_bs")->setVal(estimate_bs[ch_i_]);
      ws_->var("N_bd")->setVal(estimate_bd[ch_i_]);
      ws_->var("N_peak")->setVal(estimate_peak[ch_i_]);
      ws_->var("N_semi")->setVal(estimate_semi[ch_i_]);
      ws_->var("N_comb")->setVal(estimate_comb[ch_i_]);
    }
    else {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
        	ws_->var(name("N_bu", i, j))->setVal(N_bu_val[i][j]);
        	ws_->var(name("N_bu", i, j))->setError(N_bu_err[i][j]);
          ws_->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
          ws_->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
          ws_->var(name("N_peak", i, j))->setVal(estimate2D_peak[i][j]);
          ws_->var(name("N_semi", i, j))->setVal(estimate2D_semi[i][j] == 0 ? 0.1 : estimate2D_semi[i][j] );
          ws_->var(name("N_comb", i, j))->setVal( estimate2D_comb[i][j] == 0 ? 0.1 : estimate2D_comb[i][j] );
        }
      }
    }
    if (bd_constr_) ws_->var("Bd_over_Bs")->setVal(estimate_bd[0] / estimate_bs[0]);
  }
}

void pdf_fitData::set_final_pdf() {
	cout << red_color_bold << "setting final pdf" << default_console_color << endl;
	define_perchannel_pdf();
	if (simul_) define_simul();
	ws_->Print();
	cout << "done" << endl;
}

void pdf_fitData::define_perchannel_pdf () {
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (BF_ > 0) {
        define_constraints(i, j);
      }
      if (berns_) define_comb2(i, j);
      define_total_extended(i, j);
    }
  }
}

void pdf_fitData::define_constraints(int i, int j) {

  if (i == 0 && j == 0) {
  	RooGaussian fs_over_fu_gau("fs_over_fu_gau", "fs_over_fu_gau", *ws_->var("fs_over_fu"), RooConst(fs_over_fu_val), RooConst(fs_over_fu_err));
//  	RooLognormal fs_over_fu_gau("fs_over_fu_gau", "fs_over_fu_gau", *ws_->var("fs_over_fu"), RooConst(fs_over_fu_val), RooConst(1. + fs_over_fu_err/fs_over_fu_val));
    ws_->import(fs_over_fu_gau);

    RooGaussian one_over_BRBR_gau("one_over_BRBR_gau", "one_over_BRBR_gau", *ws_->var("one_over_BRBR"), RooConst(one_over_BRBR_val), RooConst(one_over_BRBR_err));
//    RooLognormal one_over_BRBR_gau("one_over_BRBR_gau", "one_over_BRBR_gau", *ws_->var("one_over_BRBR"), RooConst(one_over_BRBR_val), RooConst(1. + one_over_BRBR_err/one_over_BRBR_val));
    ws_->import(one_over_BRBR_gau);
  }

  RooGaussian N_bu_gau(name("N_bu_gau", i, j), "N_bu_gau", *ws_->var(name("N_bu", i, j)), RooConst(ws_->var(name("N_bu", i, j))->getVal()), RooConst(ws_->var(name("N_bu", i, j))->getError()));
//  RooLognormal N_bu_gau(name("N_bu_gau", i, j), "N_bu_gau", *ws_->var(name("N_bu", i, j)), RooConst(ws_->var(name("N_bu", i, j))->getVal()), RooConst(1. + ws_->var(name("N_bu", i, j))->getError() / ws_->var(name("N_bu", i, j))->getVal()));
  ws_->import(N_bu_gau);

  RooGaussian effratio_gau_bs(name("effratio_gau_bs", i, j), "effratio_gau_bs", *ws_->var(name("effratio_bs", i, j)), RooConst(effratio_bs_val[i][j]), RooConst(effratio_bs_err[i][j]));
  RooGaussian effratio_gau_bd(name("effratio_gau_bd", i, j), "effratio_gau_bd", *ws_->var(name("effratio_bd", i, j)), RooConst(effratio_bd_val[i][j]), RooConst(effratio_bd_err[i][j]));
//  RooLognormal effratio_gau_bs(name("effratio_gau_bs", i, j), "effratio_gau_bs", *ws_->var(name("effratio_bs", i, j)), RooConst(effratio_bs_val[i][j]), RooConst(1. + effratio_bs_err[i][j]/effratio_bs_val[i][j]));
//  RooLognormal effratio_gau_bd(name("effratio_gau_bd", i, j), "effratio_gau_bd", *ws_->var(name("effratio_bd", i, j)), RooConst(effratio_bd_val[i][j]), RooConst(1. + effratio_bd_err[i][j]/effratio_bd_val[i][j]));
  ws_->import(effratio_gau_bs);
  ws_->import(effratio_gau_bd);

  if (rare_constr_) {
  	if (stat_only_ && false) {
//  		RooRealVar N_semi_gau_sigma(name("N_semi_gau_sigma", i, j), "N_semi_gau_sigma", ws_->var(name("N_semi", i, j))->getError());
//  		RooGaussian N_semi_gau(name("N_semi_gau", i, j), "N_semi_gau", *ws_->var(name("N_semi", i, j)), RooConst(ws_->var(name("N_semi", i, j))->getVal()), N_semi_gau_sigma);
//  		ws_->import(N_semi_gau);
//
//  		RooRealVar N_peak_gau_sigma(name("N_peak_gau_sigma", i, j), "N_peak_gau_sigma", 1. + ws_->var(name("N_peak", i, j))->getError() / ws_->var(name("N_peak", i, j))->getVal());
//  		RooLognormal N_peak_gau(name("N_peak_gau", i, j), "N_peak_gau", *ws_->var(name("N_peak", i, j)), RooConst(ws_->var(name("N_peak", i, j))->getVal()), N_peak_gau_sigma);
//  		ws_->import(N_peak_gau);

  	}
  	else {
  		RooGaussian N_semi_gau(name("N_semi_gau", i, j), "N_semi_gau", *ws_->var(name("N_semi", i, j)), RooConst(ws_->var(name("N_semi", i, j))->getVal()), RooConst(ws_->var(name("N_semi", i, j))->getError()));
  	//  	RooLognormal N_semi_gau(name("N_semi_gau", i, j), "N_semi_gau", *ws_->var(name("N_semi", i, j)), RooConst(ws_->var(name("N_semi", i, j))->getVal()), RooConst(1. + ws_->var(name("N_semi", i, j))->getError() / ws_->var(name("N_semi", i, j))->getVal()));
  		ws_->import(N_semi_gau);

  		RooLognormal N_peak_gau(name("N_peak_gau", i, j), "N_peak_gau", *ws_->var(name("N_peak", i, j)), RooConst(ws_->var(name("N_peak", i, j))->getVal()), RooConst(1. + ws_->var(name("N_peak", i, j))->getError() / ws_->var(name("N_peak", i, j))->getVal()));
  		ws_->import(N_peak_gau);
  	}




  }

//  int eta = -1;
//  if (!simul_all_) eta = i / 2;
//  else {
//  	vector <int> indexes(get_EtaBdt_bins(i));
//  	eta = indexes[0] / 2;
//  }
//  RooGaussian Mean_gau_bs(name("Mean_gau_bs", i, j), "Mean_gau_bs", *ws_->var(name("Mean_bs", i, j)), RooConst(ws_->var(name("Mean_bs", i, j))->getVal()), RooConst(mass_scale_sys[eta]));
//  RooGaussian Mean_gau_bd(name("Mean_gau_bd", i, j), "Mean_gau_bd", *ws_->var(name("Mean_bd", i, j)), RooConst(ws_->var(name("Mean_bd", i, j))->getVal()), RooConst(mass_scale_sys[eta]));
//  ws_->import(Mean_gau_bs);
//  ws_->import(Mean_gau_bd);
//
//  RooGaussian Enne_gau_bs(name("Enne_gau_bs", i, j), "Enne_gau_bs", *ws_->var(name("Enne_bs", i, j)), RooConst(ws_->var(name("Enne_bs", i, j))->getVal()), RooConst(ws_->var(name("Enne_bs", i, j))->getError()));
//  RooGaussian Alpha_gau_bs(name("Alpha_gau_bs", i, j), "Alpha_gau_bs", *ws_->var(name("Alpha_bs", i, j)), RooConst(ws_->var(name("Alpha_bs", i, j))->getVal()), RooConst(ws_->var(name("Alpha_bs", i, j))->getError()));
//  RooGaussian Enne_gau_bd(name("Enne_gau_bd", i, j), "Enne_gau_bd", *ws_->var(name("Enne_bd", i, j)), RooConst(ws_->var(name("Enne_bd", i, j))->getVal()), RooConst(ws_->var(name("Enne_bd", i, j))->getError()));
//  RooGaussian Alpha_gau_bd(name("Alpha_gau_bd", i, j), "Alpha_gau_bd", *ws_->var(name("Alpha_bd", i, j)), RooConst(ws_->var(name("Alpha_bd", i, j))->getVal()), RooConst(ws_->var(name("Alpha_bd", i, j))->getError()));
//  ws_->import(Enne_gau_bs);
//  ws_->import(Alpha_gau_bs);
//  ws_->import(Enne_gau_bd);
//  ws_->import(Alpha_gau_bd);

}

void pdf_fitData::define_total_extended(int i, int j, bool h) {

  RooArgList pdf_list;
  if (!h) {
  	pdf_list.add(*ws_->pdf(name("pdf_bs", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_bd", i, j)));
  	if (berns_) pdf_list.add(*ws_->pdf(name("pdf_1comb", i, j)));
  	else pdf_list.add(*ws_->pdf(name("pdf_comb", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_semi", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_peak", i, j)));
  }
  else {
  	pdf_list.add(*ws_->pdf(name("pdf_h_bs", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_h_bd", i, j)));
  	if (berns_) pdf_list.add(*ws_->pdf(name("pdf_h_1comb", i, j)));
  	else pdf_list.add(*ws_->pdf(name("pdf_h_comb", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_h_semi", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_h_peak", i, j)));
  }


  RooArgList N_list("varlist");
  if (BF_ > 0) N_list.add(*ws_->function(name("N_bs_formula", i, j)));
  else N_list.add(*ws_->var(name("N_bs", i, j)));
  if ((SM_ || bd_constr_) && BF_ < 2) N_list.add(*ws_->function(name("N_bd_constr", i, j)));
  else if (BF_ > 1) N_list.add(*ws_->function(name("N_bd_formula", i, j)));
  else N_list.add(*ws_->var(name("N_bd", i, j)));

  RooFormulaVar N_comb2(name("N_comb2", i, j), "", "@0*@0", RooArgList(*ws_->var(name("N_comb", i, j))));
  ws_->import(N_comb2);
//  N_list.add(*ws_->function(name("N_comb2", i, j)));
  N_list.add(*ws_->var(name("N_comb", i, j)));
  N_list.add(*ws_->var(name("N_semi", i, j)));
  N_list.add(*ws_->var(name("N_peak", i, j)));


//  RooArgList constraints_list(*ws_->pdf(name("N_bu_gau", i, j)), *ws_->pdf("fs_over_fu_gau"), *ws_->pdf("one_over_BRBR_gau"), *ws_->pdf(name("effratio_gau_bs", i, j)));
  RooArgList constraints_list(*ws_->pdf(name("N_bu_gau", i, j)), *ws_->pdf(name("effratio_gau_bs", i, j)));
//  if (i == 0 && j == 0) {
//  	constraints_list.add(*ws_->pdf("fs_over_fu_gau"));
//  	constraints_list.add(*ws_->pdf("one_over_BRBR_gau"));
//  }

  if (rare_constr_) {
  	constraints_list.add(*ws_->pdf(name("N_peak_gau", i, j)));
  	constraints_list.add(*ws_->pdf(name("N_semi_gau", i, j)));
  }
  if (BF_ > 1) {
  	constraints_list.add(*ws_->pdf(name("effratio_gau_bd", i, j)));
  }


  if (!h) {

  	RooAddPdf pdf_ext_sum(name("pdf_ext_sum", i, j), "pdf_ext_sum", pdf_list, N_list);
  	RooProdPdf constraints_pdfs(name("pdf_constraints", i, j), "pdf_constraints", constraints_list);
  	RooProdPdf pdf_ext_total(name("pdf_ext_total", i, j), "pdf_ext_total", RooArgList(pdf_ext_sum, constraints_pdfs));
  	ws_->import(pdf_ext_total);
  }
  else {
  	RooAddPdf pdf_ext_sum(name("pdf_h_ext_sum", i, j), "pdf_ext_sum", pdf_list, N_list);
  	RooProdPdf constraints_pdfs(name("pdf_h_constraints", i, j), "pdf_constraints", constraints_list);
  	RooProdPdf pdf_ext_total(name("pdf_h_ext_total", i, j), "pdf_ext_total", RooArgList(pdf_ext_sum, constraints_pdfs));
  	ws_->import(pdf_ext_total);
  }
}

void pdf_fitData::define_simul() {
  if (!simul_bdt_ && !simul_all_) {
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *ws_->cat("etacat"));
//    RooSimultaneous pdf_sim("pdf_ext_simul0", "simultaneous pdf", *ws_->cat("etacat"));
    RooSimultaneous pdf_sim_noconstr("pdf_ext_simul_noconstr", "simultaneous pdf without constraints", *ws_->cat("etacat"));
    for (int i = 0; i < channels; i++) {
    	if (years_ == "0" && i > 1) break;
    	if (years_ == "1" && i < 2) continue;
      pdf_sim.addPdf(*ws_->pdf(name("pdf_ext_total", i)), name("etacat", i));
      pdf_sim_noconstr.addPdf(*ws_->pdf(name("pdf_ext_sum", i)), name("etacat", i));
    }

    RooArgList final_list(pdf_sim, *ws_->pdf("fs_over_fu_gau"), *ws_->pdf("one_over_BRBR_gau"));
    RooProdPdf final_pdf("pdf_ext_simul", "pdf_ext_simul", final_list);

//    ws_->import(final_pdf);
    ws_->import(pdf_sim);
    ws_->import(pdf_sim_noconstr);
    pdf_sim.graphVizTree("pdf_ext_simul.dot");
    pdf_sim_noconstr.graphVizTree("pdf_sim_noconstr.dot");
  }

  if (!simul_bdt_ && simul_all_) {
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *ws_->cat("allcat"));
//    RooSimultaneous pdf_sim("pdf_ext_simul0", "simultaneous pdf", *ws_->cat("allcat"));
    RooSimultaneous pdf_sim_noconstr("pdf_ext_simul_noconstr", "simultaneous pdf without constraints", *ws_->cat("allcat"));
    for (int i = 0; i < channels_all; i++) {
      vector <int> indexes(get_EtaBdt_bins(i));
      if (years_=="0" && indexes[0] > 1) continue;
      if (years_=="1" && indexes[0] < 2) continue;
//      if (indexes[1] == 0) continue;
      pdf_sim.addPdf(*ws_->pdf(name("pdf_ext_total", indexes[0], indexes[1])), Form("allcat_%d", i));
      pdf_sim_noconstr.addPdf(*ws_->pdf(name("pdf_ext_sum", indexes[0], indexes[1])), Form("allcat_%d", i));
    }

    RooArgList final_list(pdf_sim, *ws_->pdf("fs_over_fu_gau"), *ws_->pdf("one_over_BRBR_gau"));
    RooProdPdf final_pdf("pdf_ext_simul", "pdf_ext_simul", final_list);
    ws_->import(pdf_sim);
//    ws_->import(final_pdf);
    pdf_sim.graphVizTree("pdf_ext_simul_all.dot");
    ws_->import(pdf_sim_noconstr);
  }

  if (simul_bdt_ && !simul_all_) {
    RooSuperCategory* rsc = dynamic_cast<RooSuperCategory*> (ws_->obj("super_cat"));
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *rsc);
//    RooSimultaneous pdf_sim("pdf_ext_simul0", "simultaneous pdf", *rsc);
    RooSimultaneous pdf_sim_noconstr("pdf_ext_simul_noconstr", "simultaneous pdf without constraints", *rsc);
    for (int i = 0; i < channels; i++) {
    	if (years_=="0" && i > 1) break;
    	if (years_=="1" && i < 2) continue;
      for (int j = 0; j < bdt_index_max(i); j++) {
        RooArgSet icl = rsc->inputCatList();
        RooCategory* eta_c = (RooCategory*)icl.find("etacat");
        RooCategory* bdt_c = (RooCategory*)icl.find("bdtcat");
        eta_c->setIndex(i);
        bdt_c->setIndex(j);
        cout << rsc->getLabel() << " (" << rsc->getIndex() << ") " <<  i << " " << j << endl;
        pdf_sim.addPdf(*ws_->pdf(Form("pdf_ext_total_%d_%d", i, j)), rsc->getLabel());
        pdf_sim_noconstr.addPdf(*ws_->pdf(Form("pdf_ext_sum_%d_%d", i, j)), rsc->getLabel());
      }
    }

    RooArgList final_list(pdf_sim, *ws_->pdf("fs_over_fu_gau"), *ws_->pdf("one_over_BRBR_gau"));
    RooProdPdf final_pdf("pdf_ext_simul", "pdf_ext_simul", final_list);

    ws_->import(pdf_sim);
//    ws_->import(final_pdf);
    pdf_sim.graphVizTree("pdf_ext_simulBdt.dot");
    if (BF_ > 0) {
      ws_->import(pdf_sim_noconstr);
    }
  }

  if (simul_bdt_ && simul_all_) {
    cout << "simul_bdt_ can't be with simul_all_" << endl;
    exit(1);
  }
}

void pdf_fitData::set_syst() {

	ws_->var("fs_over_fu")->setConstant(!syst);
	ws_->var("one_over_BRBR")->setConstant(!syst);
	for (int i = 0; i < channels; i++) {
		for (int j = 0; j < bdt_index_max(i); j++) {
			ws_->var(name("N_bu", i, j))->setConstant(!syst);
			ws_->var(name("effratio_bs", i, j))->setConstant(!syst);
			if (rare_constr_) {
				ws_->var(name("N_peak", i, j))->setConstant(!syst);
				ws_->var(name("N_semi", i, j))->setConstant(!syst);
      }
			else {
				ws_->var(name("N_peak", i, j))->setConstant(true);
				ws_->var(name("N_semi", i, j))->setConstant(true);
			}
			if (BF_ > 1) ws_->var(name("effratio_bd", i, j))->setConstant(!syst);
//			ws_->var(name("Mean_bs", i, j))->setConstant(!syst);
//			ws_->var(name("Mean_bd", i, j))->setConstant(!syst);
//			ws_->var(name("Enne_bs", i, j))->setConstant(!syst);
//			ws_->var(name("Alpha_bs", i, j))->setConstant(!syst);
//			ws_->var(name("Enne_bd", i, j))->setConstant(!syst);
//			ws_->var(name("Alpha_bd", i, j))->setConstant(!syst);

//			if (bdt_fit_) ws_->var(name("beta_bd", i, j))->setConstant(true);

//			if (BF_ == 0) {
//				ws_->var(name("Mean_bs", i, j))->setConstant(false);
//				ws_->var(name("Mean_bd", i, j))->setConstant(false);
//			}
		}
	}

//	ws_->var("BF_bs")->setMin(-10);
//	ws_->var("BF_bd")->setMin(-10);
	for (int i = 0; i < channels; i++) {
			for (int j = 0; j < bdt_index_max(i); j++) {
//				ws_->var(name("N_bu", i, j))->setMin(-10);
//				ws_->var(name("N_peak", i, j))->setMin(-100);
//				ws_->var(name("N_semi", i, j))->setMin(-100);
//				ws_->var(name("N_comb", i, j))->setMin(-100);
			}
	}

	string constraints("");
	if (syst) {
		constraints += "fs_over_fu,one_over_BRBR";
		for (int i = 0; i < channels; i++) {
			if (years_=="0" && i > 1) break;
			if (years_=="1" && i < 2) continue;
			for (int j = 0; j < bdt_index_max(i); j++) {
				constraints += "," + (string)name("N_bu", i, j);
				constraints += "," + (string)name("effratio_bs", i, j);
				if (rare_constr_) {
					constraints += "," + (string)name("N_peak", i, j);
					constraints += "," + (string)name("N_semi", i, j);
				}
				if (BF_ > 1) constraints += "," + (string)name("effratio_bd", i, j);
//				constraints += "," + (string)name("Mean_bs", i, j);
//				constraints += "," + (string)name("Mean_bd", i, j);
//
//				constraints += "," + (string)name("Enne_bs", i, j);
//				constraints += "," + (string)name("Alpha_bs", i, j);
//				constraints += "," + (string)name("Enne_bd", i, j);
//				constraints += "," + (string)name("Alpha_bd", i, j);
			}
		}
	}
	ws_->defineSet("constr", constraints.c_str());
	cout << "constrain set: ";
	ws_->set("constr")->Print();
}

void pdf_fitData::save() {
  ostringstream output_ws;
  output_ws << "./output/ws_fitData_" << meth_;
  if (simul_)          output_ws << "_simul" << channels;
  else                 output_ws << "_" << ch_s_;
  if (SM_)             output_ws << "_SM";
  else if (bd_constr_) output_ws << "_BdConst";
  if (bdt_fit_)        output_ws << "_2D";
  if (pee)             output_ws << "_PEE";
  output_ws << ".root";
  ws_->SaveAs(output_ws.str().c_str());
}

void pdf_fitData::significance() {
	cout << red_color_bold << "significance" << default_console_color << endl;
	freeze_norm(freeze);

	if (sign < 0) return;

  if (sign == 0) sig_hand();
  else {
    make_models();
    if (sign == 1) sig_plhc();
    if (sign == 2) sig_plhts();
    if (sign == 3) sig_hybrid_plhts();
    if (sign == 4) sig_hybrid_roplhts();

    string ws_output("ws_/ws");
    if (simul_all_) ws_output += "_simulAll";
    if (simul_bdt_) ws_output += "_simulBdt";
    if (bdt_fit_) ws_output += "_2D";
    if (BF_>0) ws_output += "_BF";
    if (pee) ws_output += "_PEE";
    if (syst) ws_output += "_syst";
    ws_output += ".root";
    ws_->SaveAs(ws_output.c_str());
  }
}

Double_t pdf_fitData::sig_hand() {
  /// by hand
//  RooProdPdf* modelpdf = new RooProdPdf("model_pdf", "model_pdf", *ws_->pdf(pdfname.c_str()), *ws_->pdf("prior") );
//	ws_->import(*modelpdf);
	RooAbsPdf* modelpdf = ws_->pdf(pdfname.c_str());
  Double_t minNLL = fit_pdf(true, modelpdf->GetName(), 2, false)->minNll();
  RooRealVar *arg_var = 0;
  RooArgSet *all_vars = ws_->pdf(modelpdf->GetName())->getVariables();
  TIterator* vars_it = all_vars->createIterator();
  size_t found;
  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd && BF_ == 0) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";
  double null = 0.;
  if (SMIsNull && BF_ > 1) {
    if (Bd) null = Bd2MuMu_SM_BF_val;
    else null = Bs2MuMu_SM_BF_val;
  }
  else if (SMIsNull) {
    cout << "SMIsNull works only with BF = 2" << endl;
    exit(1);
  }
  if (!doubleNull) {
  	while ( (arg_var = (RooRealVar*)vars_it->Next())) {
  		string name(arg_var->GetName());
  		found = name.find(alt_name);
  		if (found != string::npos) {
  			arg_var->setVal(null);
  			arg_var->setConstant(1);
  		}
  		if (SM_ || bd_constr_){
  			found = name.find("Bd_over_Bs");
  			if (found!=string::npos) {
  				arg_var->setVal(0);
  				arg_var->setConstant(1);
  			}
  		}
  	}
  }
  else {
  	ws_->var("BF_bs")->setVal(0);
  	ws_->var("BF_bs")->setConstant(1);
  	ws_->var("BF_bd")->setVal(0);
  	ws_->var("BF_bd")->setConstant(1);
  }

  ///
  Double_t newNLL = fit_pdf(true, modelpdf->GetName(), 2, false)->minNll();
  ///

  if (!doubleNull) {
  	TIterator* vars_after = all_vars->createIterator();
  	while ( (arg_var = (RooRealVar*)vars_after->Next())) {
  		string name(arg_var->GetName());
  		found = name.find(alt_name);
  		if (found != string::npos) arg_var->setConstant(0);
  		if (SM_ || bd_constr_) {
  			found = name.find("Bd_over_Bs");
  			if (found!=string::npos) arg_var->setConstant(0);
  		}
  	}
  }
  else {
  	ws_->var("BF_bs")->setConstant(0);
  	ws_->var("BF_bd")->setConstant(0);
  }

  Double_t deltaLL = newNLL - minNLL;
  if (deltaLL > 0) deltaLL = -1*deltaLL;
  int dof = 1;
  if (doubleNull) dof = 2;
  Double_t signif;
  signif = RooStats::PValueToSignificance(0.5*TMath::Prob(-2*deltaLL, dof));
//  if (dof == 2) signif = RooStats::PValueToSignificance(0.5*TMath::Prob(-2*deltaLL, dof));
//  else signif = deltaLL>0 ? sqrt(2*deltaLL) : -sqrt(-2*deltaLL) ;
  if (verbosity > 0) {
    cout << "H1 minNLL = " << minNLL << endl;
    cout << "H0 minNLL = " << newNLL << endl;
    cout << "significance (by hand) = " << signif << endl << endl;
  }
  return signif;
}

void pdf_fitData::sig_plhc() {
  ModelConfig model;
  model.SetName("model");
  RooArgSet poi;
  RooArgSet CO;
//  if (pee) {
//    CO.add(*ws_->var("ReducedMassRes"));
//  }
  model.SetWorkspace(*ws_);
//  RooProdPdf* modelpdf = new RooProdPdf("model_pdf", "model_pdf", *ws_->pdf(pdfname.c_str()), *ws_->pdf("prior") );
  model.SetPdf(*ws_->pdf(pdfname.c_str()));
//  model.SetPdf(*modelpdf);

  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd && BF_ == 0) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";
  double null = 0.;
  if (SMIsNull && BF_ > 1) {
    if (Bd) null = Bd2MuMu_SM_BF_val;
    else null = Bs2MuMu_SM_BF_val;
  }
  else if (SMIsNull) {
    cout << "SMIsNull works only with BF = 2" << endl;
    exit(1);
  }
  if (!doubleNull) {
  	if (BF_ == 0) {
  		for (int i = 0; i < channels; i++) {
  			for (int j = 0; j < bdt_index_max(i); j++) {
  				poi.add(*ws_->var(name(alt_name.c_str(), i, j)));
  				poi.setRealValue(name(alt_name.c_str(), i, j), 0);
  			}
  		}
  	}
  	else {
  		poi.add(*ws_->var(alt_name.c_str()));
  		poi.setRealValue(alt_name.c_str(), null);
  	}
  	if (bd_constr_) {
  		poi.add(*ws_->var("Bd_over_Bs"));
  		poi.setRealValue("Bd_over_Bs", 0);
  	}
  }
  else {
  	poi.add(*ws_->var("BF_bs"));
  	poi.setRealValue("BF_bs", 0);
  	poi.add(*ws_->var("BF_bd"));
  	poi.setRealValue("BF_bd", 0);
  }



  ProfileLikelihoodCalculator plc;
  plc.SetData(*ws_->data("global_data"));
  plc.SetModel(model);
//  if (pee) plc.SetConditionalObservables(CO);
  plc.SetNullParameters(poi);
  HypoTestResult* htr = plc.GetHypoTest();
  cout << "ProfileLikelihoodCalculator: The p-value for the null is " << htr->NullPValue() << "; The significance for the null is " << htr->Significance() << endl;
  model.SetSnapshot(poi);
  ws_->import(model);
}

void pdf_fitData::make_models() {
  vector <vector <double> > N_bs(channels, vector <double> (channels_bdt) );
  vector <vector <double> > N_bd(channels, vector <double> (channels_bdt) );
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (BF_ == 0) N_bs[i][j] = ws_->var(name("N_bs", i, j))->getVal();
      else N_bs[i][j] = ws_->function(name("N_bs_formula", i, j))->getVal();
      if (!SM_ && !bd_constr_ && BF_ < 2) N_bd[i][j] = ws_->var(name("N_bd", i, j))->getVal();
      else if (BF_ > 1) N_bd[i][j] = ws_->function(name("N_bd_formula", i, j))->getVal();
    }
  }

  /// obs
  string observables("Mass");
  if (bdt_fit_) {
  	if (channels == 4) observables += ",bdt_0,bdt_1,bdt_2,bdt_3";
  	else observables += ",bdt_0,bdt_1";
//  	observables += ",bdt";
  }
  if (simul_) {
  	if (!simul_bdt_ && !simul_all_) observables += ",etacat";
  	else if (simul_bdt_) observables += ",bdtcat";
  	else if (simul_all_) observables += ",allcat";
  }
  ws_->defineSet("obs", observables.c_str());

  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd && BF_ == 0) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";

  /// poi
  ostringstream name_poi;
  if (BF_ == 0) {
    if (simul_) {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
          if (i != 0 || j != 0) name_poi << ",";
          name_poi << name(alt_name.c_str(), i, j);
        }
      }
    }
    else {
      name_poi << alt_name;
    }
  }
  else {
    name_poi << alt_name;
  }
  if (bd_constr_) name_poi << ",Bd_over_Bs";
  ws_->defineSet("poi", name_poi.str().c_str());

  /// nui
  RooArgSet nuisanceParams;
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      nuisanceParams.add(*ws_->var(name("N_comb", i, j)));
      if (rare_constr_) {
      	nuisanceParams.add(*ws_->var(name("N_peak", i, j)));
      	nuisanceParams.add(*ws_->var(name("N_semi", i, j)));
      }
      nuisanceParams.add(*ws_->var(name("exp_comb", i, j)));
      if (!SM_ && !bd_constr_ && BF_ < 2 && !Bd) nuisanceParams.add(*ws_->var(name("N_bd", i, j)));
      else if (Bd && BF_ < 1) nuisanceParams.add(*ws_->var(name("N_bs", i, j)));
      if (BF_ > 0 && syst) {
        nuisanceParams.add(*ws_->var(name("effratio_bs", i, j)));
        nuisanceParams.add(*ws_->var(name("N_bu", i, j)));
        if (BF_ > 1) {
          nuisanceParams.add(*ws_->var(name("effratio_bd", i, j)));
        }
      }
    }
  }
  if (BF_ > 0) {
    nuisanceParams.add(*ws_->var("fs_over_fu"));
    nuisanceParams.add(*ws_->var("one_over_BRBR"));
  }
  if (BF_ > 1 && !Bd) nuisanceParams.add(*ws_->var("BF_bd"));
  else if (BF_ > 1 && Bd) nuisanceParams.add(*ws_->var("BF_bs"));
  ws_->defineSet("nui", nuisanceParams);

  double null = 0.;
  if (SMIsNull && BF_ > 1) {
    if (Bd) null = Bd2MuMu_SM_BF_val;
    else null = Bs2MuMu_SM_BF_val;
  }
  else if (SMIsNull) {
    cout << "SMIsNull works only with BF = 2" << endl;
    exit(1);
  }

  RooArgSet CO;
  if (pee) {
    CO.add(*ws_->var("ReducedMassRes"));
    ws_->defineSet("CO", "ReducedMassRes");
  }

// RooProdPdf* modelpdf = new RooProdPdf("model_pdf", "model_pdf", *ws_->pdf(pdfname.c_str()), *ws_->pdf("prior") );
//  ws_->import(*modelpdf);

  ModelConfig* H1 = new ModelConfig("H1", "background + signal hypothesis", ws_);
//  if (pee) H1->SetConditionalObservables(*ws_->set("CO"));
  H1->SetPdf(*ws_->pdf(pdfname.c_str()));
  H1->SetParametersOfInterest(*ws_->set("poi"));
  H1->SetObservables(*ws_->set("obs"));
  H1->SetNuisanceParameters(*ws_->set("nui"));
  if (BF_ > 0 && syst) H1->SetConstraintParameters(*ws_->set("constr"));
  H1->SetSnapshot(*ws_->set("poi"));

  ModelConfig* H0 = new ModelConfig("H0", "null hypothesis", ws_);
//  if (pee) H0->SetConditionalObservables(*ws_->set("CO"));
  H0->SetPdf(*ws_->pdf(pdfname.c_str()));
  H0->SetParametersOfInterest(*ws_->set("poi"));
  H0->SetObservables(*ws_->set("obs"));
  H0->SetNuisanceParameters(*ws_->set("nui"));
  if (BF_ > 0 && syst) H0->SetConstraintParameters(*ws_->set("constr"));
  if (BF_ == 0) {
    if (simul_) {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
          ws_->var(name(alt_name.c_str(), i, j))->setVal(0.0);
        }
      }
    }
    else {
      ws_->var(alt_name.c_str())->setVal(0.0);
    }
  }
  else {
    ws_->var(alt_name.c_str())->setVal(null);
  }
  if (bd_constr_) {
    ws_->var("Bd_over_Bs")->setVal(0.0);
  }
  H0->SetSnapshot(*ws_->set("poi"));

  ws_->import(*H0);
  ws_->import(*H1);

  RooAbsPdf* nuisPdf = RooStats::MakeNuisancePdf(*H0,"nui_pdf");
  cout <<"nuisance pdf = " << endl;
  nuisPdf->Print();
  ws_->import(*nuisPdf);
}

void pdf_fitData::sig_plhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  ProfileLikelihoodTestStat pl_ts(*ws_->pdf(pdfname.c_str()));
//  pl_ts.SetPrintLevel(2);
  pl_ts.SetOneSidedDiscovery(true);
//  if (pee) pl_ts.SetConditionalObservables(*ws_->set("CO"));

  ProofConfig* pc = NULL;
  pc = new ProofConfig(*ws_, proof, Form("workers=%d", proof), kTRUE); // machine with 4 cores

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, NExp);
  if (pc && proof > 1) mcSampler_pl->SetProofConfig(pc);

  FrequentistCalculator frequCalc(*ws_->data("global_data"), *H1,*H0, mcSampler_pl); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
  frequCalc.SetToys(NExp, NExp/4.);
  HypoTestResult *htr_pl = frequCalc.GetHypoTest();
  plot_hypotest(htr_pl);
  cout << "ProfileLikelihoodTestStat + frequentist: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << " \\pm " << htr_pl->SignificanceError() << endl;
}

void pdf_fitData::sig_hybrid_plhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  ProfileLikelihoodTestStat pl_ts(*ws_->pdf(pdfname.c_str()));
//  pl_ts.SetPrintLevel(2);
  pl_ts.SetOneSidedDiscovery(true);
//  if (pee) pl_ts.SetConditionalObservables(*ws_->set("CO"));

  ProofConfig* pc = NULL;
  pc = new ProofConfig(*ws_, proof, Form("workers=%d", proof), kTRUE); // machine with "proof" cores

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, NExp);
  if (pc && proof > 1) mcSampler_pl->SetProofConfig(pc);

//  RooSimultaneous *sim  = dynamic_cast<RooSimultaneous *>(ws_->pdf("pdf_ext_simul"));
//  RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().clone(sim->indexCat().GetName());
//  for (int ic = 0, nc = cat->numBins((const char *)0); ic < nc; ++ic) {
//    cat->setBin(ic);
//    cout << cat->getLabel() << "";
//    RooAbsPdf* pdf_i = sim->getPdf(cat->getLabel());
//    cout << pdf_i->GetName() << endl;
//  }

  HybridCalculator hibrCalc(*ws_->data("global_data"), *H1, *H0, mcSampler_pl);
  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("nui_pdf"));
  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("nui_pdf"));
//  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("prior"));
//  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("prior"));
  hibrCalc.SetToys(NExp, NExp);

  HypoTestResult *htr_pl = hibrCalc.GetHypoTest();
  plot_hypotest(htr_pl);
  cout << "ProfileLikelihoodTestStat + hybrid: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << " \\pm " << htr_pl->SignificanceError() << endl;
}

void pdf_fitData::sig_hybrid_roplhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  RatioOfProfiledLikelihoodsTestStat ropl_ts(*H1->GetPdf(),*H0->GetPdf(), ws_->set("poi"));
  ropl_ts.SetSubtractMLE(false);
//  if (pee) ropl_ts.SetConditionalObservables(*ws_->set("CO"));

  ProofConfig* pc = NULL;
  pc = new ProofConfig(*ws_, proof, Form("workers=%d", proof), kTRUE); // machine with 4 cores

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(ropl_ts, NExp);
  if(pc && proof > 1) mcSampler_pl->SetProofConfig(pc);

  HybridCalculator hibrCalc(*ws_->data("global_data"), *H1, *H0, mcSampler_pl);
  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("nui_pdf"));
  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("nui_pdf"));
//  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("prior"));
//  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("prior"));
  hibrCalc.SetToys(NExp, NExp/4.);

  HypoTestResult *htr_pl = hibrCalc.GetHypoTest();
  plot_hypotest(htr_pl);
  cout << "RatioOfProfiledLikelihoodsTestStat + hybrid: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << " \\pm " << htr_pl->SignificanceError() << endl;
}

void pdf_fitData::plot_hypotest(HypoTestResult *hts) {
  hts->Print();
  TCanvas c_plot("c_plot", "c_plot", 600, 600);
  HypoTestPlot * plot = new HypoTestPlot(*hts, 100);
//  plot->Draw();
//  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".gif").c_str());
//  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".pdf").c_str());
//  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".root").c_str());
//  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".C").c_str());
//  if (NExp < 100) {
//  	delete plot;
//  	return;
//  }
  plot->SetLogYaxis(true);
  plot->Draw();
  c_plot.Print( (get_address(Form("hypotest_log_%d", sign), Bd ? "Bd" : "") + ".gif").c_str());
  c_plot.Print( (get_address(Form("hypotest_log_%d", sign), Bd ? "Bd" : "") + ".pdf").c_str());
  c_plot.Print( (get_address(Form("hypotest_log_%d", sign), Bd ? "Bd" : "") + ".root").c_str());
  c_plot.Print( (get_address(Form("hypotest_log_%d", sign), Bd ? "Bd" : "") + ".C").c_str());
  delete plot;
}

void pdf_fitData::make_prior() {
  cout << "making prior" << endl;
  vector <vector <RooGaussian*> > prior_bd(channels, vector <RooGaussian*> (channels_bdt));
  vector <vector <RooGaussian*> > prior_semi(channels, vector <RooGaussian*> (channels_bdt));
  vector <vector <RooGaussian*> > prior_comb(channels, vector <RooGaussian*> (channels_bdt));
  vector <vector <RooGaussian*> > prior_exp(channels, vector <RooGaussian*> (channels_bdt));

  vector <vector <RooGamma*> > prior_comb2(channels, vector <RooGamma*> (channels_bdt));

  RooArgList prior_list("prior_list");

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (!SM_ && !bd_constr_ && BF_ < 2 && !Bd) prior_bd[i][j] = new RooGaussian(name("prior_bd", i, j), name("prior_bd", i, j), *ws_->var(name("N_bd", i, j)), RooConst(ws_->var(name("N_bd", i, j))->getVal()), RooConst(ws_->var(name("N_bd", i, j))->getError()));
      else if (BF_ < 2 && Bd) prior_bd[i][j] = new RooGaussian(name("prior_bs", i, j), name("prior_bs", i, j), *ws_->var(name("N_bs", i, j)), RooConst(ws_->var(name("N_bs", i, j))->getVal()), RooConst(ws_->var(name("N_bs", i, j))->getError()));
      prior_semi[i][j] = new RooGaussian(name("prior_semi", i, j), name("prior_semi", i, j), *ws_->var(name("N_semi", i, j)), RooConst(ws_->var(name("N_semi", i, j))->getVal()), RooConst(ws_->var(name("N_semi", i, j))->getError()));
      prior_comb[i][j] = new RooGaussian(name("prior_comb", i, j), name("prior_comb", i, j), *ws_->var(name("N_comb", i, j)), RooConst(ws_->var(name("N_comb", i, j))->getVal()), RooConst(ws_->var(name("N_comb", i, j))->getError()));
      prior_exp[i][j] = new RooGaussian(name("prior_exp", i, j), name("prior_exp", i, j), *ws_->var(name("exp_comb", i, j)), RooConst(ws_->var(name("exp_comb", i, j))->getVal()), RooConst(ws_->var(name("exp_comb", i, j))->getError()));
      prior_list.add(*prior_bd[i][j]);
      if (!rare_constr_) {
      	prior_list.add(*prior_semi[i][j]);
      }
      prior_list.add(*prior_comb[i][j]);
      prior_list.add(*prior_exp[i][j]);

//      prior_comb2[i][j] = new RooGamma(name("prior_comb", i, j), name("prior_comb", i, j), *ws_->var(name("N_comb", i, j)), RooConst(ws_->var(name("N_comb", i, j))->getVal() + 1), RooConst(1.), RooConst(0.));
//      prior_list.add(*prior_comb2[i][j]);
    }
  }
  if (BF_ > 1 && !Bd) {
    RooGaussian* prior_bf_bd = new RooGaussian("prior_bf_bd", "prior_bf_bd", *ws_->var("BF_bd"), RooConst(ws_->var("BF_bd")->getVal()), RooConst(ws_->var("BF_bd")->getError()));
    prior_list.add(*prior_bf_bd);
  }
  else if (BF_ > 1 && Bd) {
    RooGaussian* prior_bf_bs = new RooGaussian("prior_bf_bs", "prior_bf_bs", *ws_->var("BF_bs"), RooConst(ws_->var("BF_bs")->getVal()), RooConst(ws_->var("BF_bs")->getError()));
    prior_list.add(*prior_bf_bs);
  }

  RooProdPdf prior("prior", "prior", prior_list);
  ws_->import(prior);
}

void pdf_fitData::setnewlumi() {
  if (BF_ > 0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
      	Double_t old_val = ws_->var(name("N_bu", i, j))->getVal();
      	Double_t old_err = ws_->var(name("N_bu", i, j))->getError();;
        Double_t new_val = old_val*lumi;
        Double_t new_err = old_err*sqrt(lumi);
        ws_->var(name("N_bu", i, j))->setMax(10*new_val);
        ws_->var(name("N_bu", i, j))->setVal(new_val);
        ws_->var(name("N_bu", i, j))->setError(new_err);
        cout << "channel " << i << " " << j << "; Bu expected " << new_val << " +- " << new_err << "; Bs expected: " << ws_->function(name("N_bs_formula", i, j))->getVal() << "; Bd expected: " << ws_->function(name("N_bd_formula", i, j))->getVal() << endl;
        ws_->var(name("N_peak", i, j))->setVal(ws_->var(name("N_peak", i, j))->getVal()*lumi);
        ws_->var(name("N_semi", i, j))->setVal(ws_->var(name("N_semi", i, j))->getVal()*lumi);
        ws_->var(name("N_comb", i, j))->setVal(ws_->var(name("N_comb", i, j))->getVal()*lumi);
        cout << "channel " << i << " " << j << "; peak expected: " << ws_->var(name("N_peak", i, j))->getVal() << "; semi expected: " << ws_->var(name("N_semi", i, j))->getVal()  << "; comb expected: " << ws_->var(name("N_comb", i, j))->getVal() << endl;
      }
    }
  }
}

void pdf_fitData::randomize_constraints(RooWorkspace* ws) {
   //generating random constrains
  if (BF_ == 0) {
    cout << "no BF!" << endl;
    return;
  }

//  RooDataSet* fs_over_fu_ds = ws->pdf("fs_over_fu_gau")->generate(RooArgSet(*ws->var("fs_over_fu")), 1);
//  ws->var("fs_over_fu")->setVal(fs_over_fu_ds->get(0)->getRealValue("fs_over_fu"));
//
//  RooDataSet* one_over_BRBR_ds = ws->pdf("one_over_BRBR_gau")->generate(RooArgSet(*ws->var("one_over_BRBR")), 1);
//  ws->var("one_over_BRBR")->setVal(one_over_BRBR_ds->get(0)->getRealValue("one_over_BRBR"));
//
//  for (int i = 0; i < channels; i++) {
//    for (int j = 0; j < bdt_index_max(i); j++) {
//      RooDataSet* N_bu_ds = ws->pdf(name("N_bu_gau", i, j))->generate(RooArgSet(*ws->var(name("N_bu", i, j))), 1);
//      RooDataSet* effratio_bs_ds = ws->pdf(name("effratio_gau_bs", i, j))->generate(RooArgSet(*ws->var(name("effratio_bs", i, j))), 1);
//      if (rare_constr_) {
//      	RooDataSet* N_peak_ds = ws->pdf(name("N_peak_gau", i, j))->generate(RooArgSet(*ws->var(name("N_peak", i, j))), 1);
//      	RooDataSet* N_semi_ds = ws->pdf(name("N_semi_gau", i, j))->generate(RooArgSet(*ws->var(name("N_semi", i, j))), 1);
//      	ws->var(name("N_peak", i, j))->setVal(N_peak_ds->get(0)->getRealValue(name("N_peak", i, j)));
//      	ws->var(name("N_semi", i, j))->setVal(N_semi_ds->get(0)->getRealValue(name("N_semi", i, j)));
//      }
//
//      ws->var(name("N_bu", i, j))->setVal(N_bu_ds->get(0)->getRealValue(name("N_bu", i, j)));
//      ws->var(name("effratio_bs", i, j))->setVal(effratio_bs_ds->get(0)->getRealValue(name("effratio_bs", i, j)));
//      if (BF_ > 1) {
//        RooDataSet* effratio_bd_ds = ws->pdf(name("effratio_gau_bd", i, j))->generate(RooArgSet(*ws->var(name("effratio_bd", i, j))), 1);
//        ws->var(name("effratio_bd", i, j))->setVal(effratio_bd_ds->get(0)->getRealValue(name("effratio_bd", i, j)));
//      }
//    }
//  }

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      RooGaussian exp_comb_gau(name("exp_comb_gau", i, j), "exp_comb_gau", *ws_->var(name("exp_comb", i, j)), RooConst(ws_->var(name("exp_comb", i, j))->getVal()), RooConst(ws_->var(name("exp_comb", i, j))->getVal()));
      RooDataSet* exp_comb_ds = exp_comb_gau.generate(RooArgSet(*ws->var(name("exp_comb", i, j))), 1);
      ws->var(name("exp_comb", i, j))->setVal(exp_comb_ds->get(0)->getRealValue(name("exp_comb", i, j)));
    }
  }
}

void pdf_fitData::extract_N_inRanges() {
  if (BF_ < 2) return;
  cout << red_color_bold << "extracting events in ranges..." << default_console_color << endl;
  string full_output = "output/yields.tex";
  FILE* file_out = fopen(full_output.c_str(), "w");

  string title_i[5] = {"$N_{B_s^0}$", "$N_{B^0}$", "$N_{\\mathrm{peak}}$", "$N_{\\mathrm{semi}}$", "$N_{\\mathrm{comb}}$"};
  string N_i[5] = {"N_bs_formula", "N_bd_formula", "N_peak", "N_semi", "N_comb"};
  string pdf_i[5] = {"pdf_bs", "pdf_bd", "pdf_peak", "pdf_semi", "pdf_comb"};

  fprintf(file_out, "\\begin{table}\n");
  fprintf(file_out, "\\centering\n");
  fprintf(file_out, "\\caption{Final invariant mass yields evaluated with the UML}\n");
  fprintf(file_out, "\\label{tab:UMLfinalyields}\n");
  fprintf(file_out, "\\begin{tabular}{|l|c|c|c|c|c|}\n");
  fprintf(file_out, "\\hline \n");

  for (int i = 0; i < channels; i++) {
    fprintf(file_out, "\\hline \n");
    fprintf(file_out, " \\multicolumn{6}{|c|}{Channel: %s %s} \\\\ \n", i%2==0? "barrel" : "endcap", i < 2 ? "2011" : "2012");
    fprintf(file_out, "\\hline \n");
    fprintf(file_out, "Variable  & low SB & $B^0$ window & $B_s^0$ window & high SB & \\textbf{all} \\\\ \n");
    fprintf(file_out, "\\hline \n");
    for (int j = 0; j < bdt_index_max(i); j++) {
      vector < Double_t> total(5, 0.);
      for (int l = 0; l < 5; l++) {
        fprintf(file_out, "%s ", title_i[l].c_str());
        for (unsigned int k = 0; k < massrange_names.size(); k++) {
          RooAbsReal* rar = ws_->pdf(name(pdf_i[l].c_str(), i, j))->createIntegral(RooArgSet(*ws_->var("Mass")), RooArgSet(*ws_->var("Mass")), massrange_names[k].c_str());
          Double_t N;
          if (l < 2) {
            RooAbsReal* rrv = ws_->function(name(N_i[l].c_str(), i, j));
            N = rrv->getVal() * rar->getVal();
          }
          else {
            RooRealVar* rrv = ws_->var(name(N_i[l].c_str(), i, j));
            N = rrv->getVal() * rar->getVal();
          }
          fprintf(file_out, "& %.2f ", N);
          total[k] += N;
        }
        Double_t N_all;
        if (l < 2) {
          N_all = ws_->function(name(N_i[l].c_str(), i, j))->getVal();
        }
        else {
          N_all = ws_->var(name(N_i[l].c_str(), i, j))->getVal();
        }
        fprintf(file_out, "& %.2f ", N_all);
        fprintf(file_out, " \\\\ \n");
        total[4] += N_all;
      }
      fprintf(file_out, "\\hline \n");
      fprintf(file_out, "$N_{\\mathrm{all}}$ ");
      for (unsigned int k = 0; k < massrange_names.size() + 1; k++) {
        fprintf(file_out, "& %.2f ", total[k]);
      }
      fprintf(file_out, " \\\\ \n");
    }
    fprintf(file_out, "\\hline \n");
  }
  fprintf(file_out, "\\end{tabular} \n");
  fprintf(file_out, "\\end{table} \n");
  fclose(file_out);
//  system(Form("cat %s", full_output.c_str()));
  cout << "tex file saved in " << full_output << endl;
}

void pdf_fitData::profile_NLL() {
  if (BF_ < 2) return;

//  RooNLLVar nll0("nll0","nll0",ws_->pdf(pdfname.c_str()),*global_data) ;

  // Construct unbinned likelihood
  RooArgSet ext_constr( *ws_->pdf("fs_over_fu_gau"), *ws_->pdf("one_over_BRBR_gau"));
  RooAbsReal* nll = ws_->pdf(pdfname.c_str())->createNLL(*global_data, NumCPU(16), Extended(), Strategy(2), Hesse(1)/*, ExternalConstraints(ext_constr)*/) ;
  // Minimize likelihood w.r.t all parameters before making plots
  RooMinuit(*nll).migrad() ;

  string var_alt("BF_bs");
  if (Bd) var_alt = "BF_bd";
  RooAbsReal* pll = nll->createProfile(*ws_->var(var_alt.c_str())) ;

  RooFormulaVar* double_pll = new RooFormulaVar("double_pll", "double_pll", "2*@0", RooArgList(*pll));
  // Plot the profile likelihood
  RooPlot* frame = ws_->var(var_alt.c_str())->frame(Bins(20), Range(0, Bd ? 1e-9 : 5e-9), Title(Form("profileLL in %s", var_alt.c_str()))) ;

//  double_pll->plotOn(frame, LineColor(kRed), ShiftToZero()) ;
//  nll0.plotOn(frame, ShiftToZero()) ;
//  TH1D * h = new TH1D("h", "h", 20, 0, 1e-9);
  frame->SetMinimum(0.);
  double_pll->plotOn(frame, LineColor(kRed), ShiftToZero(), Title("")) ;
  frame->SetTitle("");
  frame->SetMinimum(0.);

  TCanvas *c = new TCanvas("c","c",600, 600);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->SetYTitle("- 2 ln L");
  if (!Bd) frame->SetXTitle("BF (B^{0}_{s}#rightarrow #mu^{+}#mu^{-})");
  else frame->SetXTitle("BF (B^{0}#rightarrow #mu^{+}#mu^{-})");
  frame->Draw();
  TLatex Tlatex;
  Tlatex.SetTextFont(42);
  Tlatex.SetTextSize(0.031);
  Tlatex.SetTextAlign(11);
  Tlatex.SetNDC();
  string latex_s("CMS - L = 5 fb^{-1} #sqrt{s} = 7 TeV, L = 20 fb^{-1} #sqrt{s} = 8 TeV");
  Tlatex.DrawLatex(0.1,0.91, latex_s.c_str());
  c->Print((get_address("profileLL", var_alt, false) + ".gif").c_str());
  c->Print((get_address("profileLL", var_alt, false) + ".pdf").c_str());
  c->Print((get_address("profileLL", var_alt, false) + ".root").c_str());

  RooCurve * curve = (RooCurve*)c->GetPrimitive(Bd ? "double_pll_Norm[BF_bd]" : "double_pll_Norm[BF_bs]");
  Int_t n = curve->GetN();
  for (int i = 0; i < n; i++) {
  	Double_t x, y;
  	curve->GetPoint(i, x, y);
//  	if ( abs(y - 1) < 0.05) cout << i << " " << x << " " << y << endl;
  	cout << i << " " << x << " " << y << endl;
  }
  delete frame;
  delete c;

}

void pdf_fitData::doublescan() {
	// Construct unbinned likelihood
	RooAbsReal* nll = ws_->pdf(pdfname.c_str())->createNLL(*global_data, NumCPU(16), Extended()) ;
	// Minimize likelihood w.r.t all parameters before making plots
	RooMinuit(*nll).migrad() ;

   double BF[2] = {ws_->var("BF_bs")->getVal(),ws_->var("BF_bd")->getVal()};

   RooAbsReal* pll = nll->createProfile(RooArgSet(*ws_->var("BF_bs"), *ws_->var("BF_bd")));
   RooFormulaVar* double_pll = new RooFormulaVar("double_pll", "double_pll", "2*@0", RooArgList(*pll));
   TH2F *h2dpll = (TH2F*)double_pll->createHistogram("h2dpll", *ws_->var("BF_bs"), Binning(50,0.,10e-9), Scaling(kFALSE), YVar(*ws_->var("BF_bd"), Binning(50,0.,2.0e-9)));

   double contours[5] = {2.29574, 6.18007, 11.8291, 19.3339, 28.7437};
   h2dpll->SetContour(5, contours);

   TCanvas *c = new TCanvas("c","c",600, 600);

   h2dpll->GetYaxis()->SetTitle("B_{d}#rightarrow#mu^{+}#mu^{-} branching fraction");
   h2dpll->GetYaxis()->SetTitleOffset(1.);
   h2dpll->GetXaxis()->SetTitle("B_{s}#rightarrow#mu^{+}#mu^{-} branching fraction");
   h2dpll->GetXaxis()->SetTitleOffset(1.);
   h2dpll->GetXaxis()->SetLabelOffset(0.01);
   h2dpll->SetStats(0);
   h2dpll->SetTitle("");
   h2dpll->SetLineWidth(2);

   h2dpll->Draw("CONT1");

   TPolyMarker pm;
   pm.SetMarkerStyle(34);
   pm.SetMarkerSize(2.0);
   pm.DrawPolyMarker(1,&BF[0],&BF[1]);

   c->Print((get_address("doubleprofileLL", "", false) + ".gif").c_str());
   c->Print((get_address("doubleprofileLL", "", false) + ".pdf").c_str());
   c->Print((get_address("doubleprofileLL", "", false) + ".root").c_str());


}

void pdf_fitData::hack_ws(string frozen_ws_address) {
  cout << "hacking 2011 shape with 2012 shape from " << frozen_ws_address << endl;
  TFile * frozen_f = new TFile(frozen_ws_address.c_str());
  RooWorkspace * frozen_ws = (RooWorkspace*)frozen_f->Get("ws");
  for (int i = 2; i < 4; i++) {
    ws_->var(name("C0_semi", i-2))->setVal(frozen_ws->var(name("C0_semi", i))->getVal());
    ws_->var(name("C1_semi", i-2))->setVal(frozen_ws->var(name("C1_semi", i))->getVal());
    ws_->var(name("C2_semi", i-2))->setVal(frozen_ws->var(name("C2_semi", i))->getVal());
    ws_->var(name("C3_semi", i-2))->setVal(frozen_ws->var(name("C3_semi", i))->getVal());
    ws_->var(name("tau_semi", i-2))->setVal(frozen_ws->var(name("tau_semi", i))->getVal());
  }
}

void pdf_fitData::reset_minmax() {
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      ws_->var(name("N_peak", i, j))->setMax(ws_->var(name("N_peak", i, j))->getVal() + 10 * ws_->var(name("N_peak", i, j))->getError());
      ws_->var(name("N_semi", i, j))->setMax(ws_->var(name("N_semi", i, j))->getVal() + 10 * ws_->var(name("N_semi", i, j))->getError());
      ws_->var(name("N_comb", i, j))->setMax(ws_->var(name("N_comb", i, j))->getVal() + 10 * ws_->var(name("N_comb", i, j))->getError());
    }
  }
}

void pdf_fitData::print_gaussian_constraints() {
	cout << "gaussian constraints" << endl;
	RooArgSet * set2 = ws_->pdf(pdfname.c_str())->getComponents();
	TIterator* it2 = set2->createIterator();
	RooAbsPdf* var_Obj2 = 0;
	while((var_Obj2 = (RooAbsPdf*)it2->Next())){
		string name = var_Obj2->GetName();
		size_t found1 = name.find("gau");
		if (found1 != string::npos) var_Obj2->Print();
	}
}

void pdf_fitData::tweak_pdf(int free) {

	ws_->var("BF_bs")->setVal(Bs2MuMu_SM_BF_val);
	ws_->var("BF_bd")->setVal(Bd2MuMu_SM_BF_val);

	bool peak_const = true;
	bool semi_const = true;
	if (free == 1) semi_const = false;
	if (free == 2) peak_const = false;
	if (free == 3) {
		semi_const = false;
		peak_const = false;
	}

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
    	ws_->var(name("N_peak", i, j))->setConstant(peak_const);
    	ws_->var(name("N_semi", i, j))->setConstant(semi_const);

//    	ws_->var(name("N_comb", i, j))->setMin(-100);

//    	ws_->var(name("exp_comb", i, j))->setMax(0.0);
//    	ws_->var(name("exp_comb", i, j))->setVal(0.0);
//    	ws_->var(name("exp_comb", i, j))->setConstant(true);
    }
  }
//  ws_->var("BF_bs")->setVal(0.);
//  ws_->var("BF_bs")->setConstant(true);
//  ws_->var("BF_bd")->setVal(0.);
//  ws_->var("BF_bd")->setConstant(true);
//  ws_->var("fs_over_fu")->setVal(fs_over_fu_val);
//	ws_->var("fs_over_fu")->setConstant(1);
//	ws_->var("one_over_BRBR")->setConstant(1);
}

void pdf_fitData::freeze_norm(bool set) {
	if (set) cout << "freezing" << endl;
	ws_->var("fs_over_fu")->setConstant(set);
	ws_->var("one_over_BRBR")->setConstant(set);
	for (int i = 0; i < channels; i++) {
		for (int j = 0; j < bdt_index_max(i); j++) {
			ws_->var(name("N_bu", i, j))->setConstant(set);
		}
	}
}

void pdf_fitData::stat_error(bool stat_only) {

	ws_->var("fs_over_fu")->setConstant(stat_only);
	ws_->var("one_over_BRBR")->setConstant(stat_only);
	for (int i = 0; i < channels; i++) {
		for (int j = 0; j < bdt_index_max(i); j++) {
			ws_->var(name("N_bu", i, j))->setConstant(stat_only);
			ws_->var(name("effratio_bs", i, j))->setConstant(stat_only);
			ws_->var(name("effratio_bd", i, j))->setConstant(stat_only);

			ws_->var(name("N_semi", i, j))->setConstant(stat_only);
			ws_->var(name("N_peak", i, j))->setConstant(stat_only);

//			ws_->var(name("N_semi_gau_sigma", i, j))->setVal(ws_->var(name("N_semi", i, j))->getError());
//			ws_->var(name("N_peak_gau_sigma", i, j))->setVal(1. + ws_->var(name("N_peak", i, j))->getError() / ws_->var(name("N_peak", i, j))->getVal());
		}
	}
}

RooHistPdf * pdf_fitData::pdf_to_hist(RooAbsPdf* kpdf) {
	string name_temp(kpdf->GetName());
	name_temp += "_h";
	TH1 * h = kpdf->createHistogram(name_temp.c_str(), *ws_->var("ReducedMassRes"), Binning(100, 0.0009, 0.045));
	name_temp = kpdf->GetName();
	name_temp += "_rdh";
	RooDataHist * rdh = new RooDataHist(name_temp.c_str(), "histo", *ws_->var("ReducedMassRes"), h);
	name_temp = kpdf->GetName();
	name_temp += "_hp";
	RooHistPdf * hp = new RooHistPdf(name_temp.c_str(), "histo", RooArgList(*ws_->var("ReducedMassRes")), *rdh, 5);
	return hp;
}


void pdf_fitData::merge_mass_to_hist() {

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
    	for (int k = 0; k < 5; k++) {
    		string pdf_name = "ReducedMassRes_pdf_";
    		string final_name(pdf_name);
    		final_name += "hist_";

    		pdf_name += source[k];
    		final_name += source[k];

    		pdf_name = name(pdf_name, i, j);
    		final_name = name(final_name, i, j); cout << pdf_name << endl;
    		RooHistPdf * hp = pdf_to_hist(ws_->pdf(pdf_name.c_str()));

    		hp->SetName(final_name.c_str());
//    		cout << hp->GetName() << endl;
    		ws_->import(*hp);
    	}
      RooProdPdf pdf_bs(name("pdf_h_bs", i, j), "pdf_h_bs", *ws_->pdf(name("ReducedMassRes_pdf_hist_bs", i, j)), Conditional(*ws_->pdf(name("CB_bs", i, j)), *ws_->var("Mass")));
      RooProdPdf pdf_bd(name("pdf_h_bd", i, j), "pdf_h_bs", *ws_->pdf(name("ReducedMassRes_pdf_hist_bd", i, j)), Conditional(*ws_->pdf(name("CB_bd", i, j)), *ws_->var("Mass")));
      RooProdPdf pdf_peak(name("pdf_h_peak", i, j), "pdf_h_bs", *ws_->pdf(name("ReducedMassRes_pdf_hist_peak", i, j)), *ws_->pdf(name("mass_peak", i, j)));
      RooProdPdf pdf_semi(name("pdf_h_semi", i, j), "pdf_h_bs", *ws_->pdf(name("ReducedMassRes_pdf_hist_semi", i, j)), *ws_->pdf(name("keys_semi", i, j)));
      RooProdPdf pdf_comb(name("pdf_h_1comb", i, j), "pdf_h_bs", *ws_->pdf(name("ReducedMassRes_pdf_hist_comb", i, j)), *ws_->pdf(name("mass_1comb", i, j)));
      ws_->import(pdf_bs);
      ws_->import(pdf_bd);
      ws_->import(pdf_peak);
      ws_->import(pdf_semi);
      ws_->import(pdf_comb);
    }
  }

}

void pdf_fitData::save_for_cls() {

   int nbins = 0;
   for (int i = 0; i < channels; i++) nbins += bdt_index_max(i);

   FILE *fp = fopen("./output/datacard.txt","w");
   fprintf(fp,"imax %d number of bins\n",nbins);
   fprintf(fp,"jmax %d number of processes minus 1\n",nbins*5-1);
   fprintf(fp,"kmax * number of nuisance parameters\n");
   fprintf(fp,"----------------------------------------------------------------------------------------------------------------------------------\n");
   fprintf(fp,"shapes * * workspace.root ws:$PROCESS\n");

   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++)
           fprintf(fp,"shapes data_obs %s workspace.root ws:%s\n",name("ch",i,j),name("data_obs",i,j));

   fprintf(fp,"----------------------------------------------------------------------------------------------------------------------------------\n");

   fprintf(fp,"bin ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++)
           fprintf(fp,"%s ",name("ch",i,j));
   fprintf(fp,"\n");

   fprintf(fp,"observation ");
   for (int i = 0; i < channels; i++) {
		for (int j = 0; j < bdt_index_max(i); j++) {
           string cut;
           if (!simul_bdt_ && !simul_all_) {
               cut = Form("etacat==etacat::etacat_%d", i);
           }
           else if (simul_bdt_ && !simul_all_) {
               cut = Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", i, j);
           }
           else if (!simul_bdt_ && simul_all_) {
               cut = Form("allcat==allcat::allcat_%d", super_index(i, j));
           }

           fprintf(fp,"%g ",global_data->sumEntries(cut.c_str()));

           RooDataSet *sub = (RooDataSet*)global_data->reduce(cut.c_str());
           ws_->import(*sub,Rename(name("data_obs",i,j)));
       }
   }
   fprintf(fp,"\n");

   fprintf(fp,"----------------------------------------------------------------------------------------------------------------------------------\n");

   fprintf(fp,"bin ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++)
           for (int k=0;k<5;k++)
               fprintf(fp,"%s ",name("ch",i,j));
   fprintf(fp,"\n");

   fprintf(fp,"process ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           fprintf(fp,"%s ",name("pdf_bs", i, j));
           fprintf(fp,"%s ",name("pdf_bd", i, j));
           fprintf(fp,"%s ",name("pdf_1comb", i, j));
           fprintf(fp,"%s ",name("pdf_semi", i, j));
           fprintf(fp,"%s ",name("pdf_peak", i, j));

       }
   fprintf(fp,"\n");

   int counter_bkg = 1;
   int counter_sig = 0;
   fprintf(fp,"process ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {
               fprintf(fp,"%d ",counter_bkg);
               counter_bkg++;

               fprintf(fp,"%d ",counter_sig);
               counter_sig--;

               fprintf(fp,"%d ",counter_bkg);
               counter_bkg++;
               fprintf(fp,"%d ",counter_bkg);
               counter_bkg++;
               fprintf(fp,"%d ",counter_bkg);
               counter_bkg++;
       }
   fprintf(fp,"\n");

   fprintf(fp,"rate ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           if (BF_ > 0) fprintf(fp,"%g ",ws_->function(name("N_bs_formula", i, j))->getVal());
           else  fprintf(fp,"%g ",ws_->var(name("N_bs", i, j))->getVal());

           if ((SM_ || bd_constr_) && BF_ < 2) fprintf(fp,"%g ",ws_->function(name("N_bd_constr", i, j))->getVal());
           else if (BF_ > 1) fprintf(fp,"%g ",ws_->function(name("N_bd_formula", i, j))->getVal());
           else fprintf(fp,"%g ",ws_->var(name("N_bd", i, j))->getVal());

           fprintf(fp,"%g ",ws_->var(name("N_comb", i, j))->getVal());
           fprintf(fp,"%g ",ws_->var(name("N_semi", i, j))->getVal());
           fprintf(fp,"%g ",ws_->var(name("N_peak", i, j))->getVal());

       }
   fprintf(fp,"\n");

   fprintf(fp,"----------------------------------------------------------------------------------------------------------------------------------\n");

   fprintf(fp,"N_bu_gau lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooGaussian *gau = (RooGaussian*)ws_->pdf(name("N_bu_gau", i, j));
           RooConstVar *mean = (RooConstVar*)gau->findServer(0);
           RooConstVar *sigma = (RooConstVar*)gau->findServer(1);

           double k = 1. + sigma->getVal()/mean->getVal();

           fprintf(fp,"%g ",k);
           fprintf(fp,"%g ",k);
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
       }
   fprintf(fp,"\n");

   fprintf(fp,"fs_over_fu_gau lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooGaussian *gau = (RooGaussian*)ws_->pdf("fs_over_fu_gau");
           RooConstVar *mean = (RooConstVar*)gau->findServer(0);
           RooConstVar *sigma = (RooConstVar*)gau->findServer(1);

           double k = 1. + sigma->getVal()/mean->getVal();

           fprintf(fp,"%g ",k);
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
       }
   fprintf(fp,"\n");

   fprintf(fp,"one_over_BRBR_gau lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooGaussian *gau = (RooGaussian*)ws_->pdf("one_over_BRBR_gau");
           RooConstVar *mean = (RooConstVar*)gau->findServer(0);
           RooConstVar *sigma = (RooConstVar*)gau->findServer(1);

           double k = 1. + sigma->getVal()/mean->getVal();

           fprintf(fp,"%g ",k);
           fprintf(fp,"%g ",k);
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
       }
   fprintf(fp,"\n");

   fprintf(fp,"effratio_gau_bs lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooGaussian *gau = (RooGaussian*)ws_->pdf(name("effratio_gau_bs", i, j));
           RooConstVar *mean = (RooConstVar*)gau->findServer(0);
           RooConstVar *sigma = (RooConstVar*)gau->findServer(1);

           double k = 1. + sigma->getVal()/mean->getVal();

           fprintf(fp,"%g ",k);
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
       }
   fprintf(fp,"\n");

   fprintf(fp,"effratio_gau_bd lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooGaussian *gau = (RooGaussian*)ws_->pdf(name("effratio_gau_bd", i, j));
           RooConstVar *mean = (RooConstVar*)gau->findServer(0);
           RooConstVar *sigma = (RooConstVar*)gau->findServer(1);

           double k = 1. + sigma->getVal()/mean->getVal();

           fprintf(fp,"- ");
           fprintf(fp,"%g ",k);
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
       }
   fprintf(fp,"\n");

   fprintf(fp,"N_semi_gau lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooGaussian *gau = (RooGaussian*)ws_->pdf(name("N_semi_gau", i, j));
           RooConstVar *mean = (RooConstVar*)gau->findServer(0);
           RooConstVar *sigma = (RooConstVar*)gau->findServer(1);

           double k = 1. + sigma->getVal()/mean->getVal();

           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"%g ",k);
           fprintf(fp,"- ");
       }
   fprintf(fp,"\n");

   fprintf(fp,"N_peak_gau lnN ");
   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {

           RooLognormal *lnN = (RooLognormal*)ws_->pdf(name("N_peak_gau", i, j));
           RooConstVar *mu = (RooConstVar*)lnN->findServer(0);
           RooConstVar *kappa = (RooConstVar*)lnN->findServer(1);

           double k = kappa->getVal();

           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"- ");
           fprintf(fp,"%g ",k);
       }
   fprintf(fp,"\n");

   for (int i = 0; i < channels; i++)
		for (int j = 0; j < bdt_index_max(i); j++) {
           fprintf(fp,"%s param %g %g [0.,1.]\n",name("B_1comb", i, j),ws_->var(name("B_1comb", i, j))->getVal(),ws_->var(name("B_1comb", i, j))->getError());
       }
   fprintf(fp,"\n");

   fclose(fp);

   ws_->SaveAs("./output/workspace.root");

}
