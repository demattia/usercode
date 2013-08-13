#include "pdf_analysis.h"

pdf_analysis::pdf_analysis(bool print, string ch_s, string range, int BF, bool SM, bool bd_constr, int simul, int simulbdt, int simulall, bool pee_, bool bdt_fit, bool final) {
  cout << "analysis constructor" << endl;

  minos = false;
  print_ = print;
  meth_ = "bdt";
  ch_s_ = ch_s;
  ch_i_ = atoi(ch_s_.c_str());
  range_ = range;

  SM_ = SM;
  bd_constr_ = bd_constr;
  verbosity = 1;
  old_tree = false;

  no_legend = false;
  channel = atoi(ch_s.c_str());

  pee = pee_;

  channels = simul;
  channels_bdt = simulbdt;
  channels_all = simulall;
  supercatdim = 12;

  simul_ = simul > 1 ? true : false;
  simul_bdt_ = simulbdt > 1 ? true: false;
  simul_all_ = simulall > 1 ? true: false;

  bdt_fit_ = bdt_fit;
  BF_ = BF;
  newcomb_ = false;

  syst = false;
  randomsyst = false;

  if (simul_all_) channels_bdt = 4;
  default_console_color = "\033[0m";
  red_color_bold = "\033[1;31m";
  green_color_bold = "\033[1;32m";

  bdt_boundaries.resize(channels);

  bdtsplit = false;
  final_ = final;

  initialize();
//  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
//  RooMsgService::instance().getStream(0).removeTopic(Eval);
//  RooMsgService::instance().getStream(1).removeTopic(Eval);
//  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().Print() ;

  bdt_cuts.resize(channels);

}

pdf_analysis::~pdf_analysis() {
  cout << "pdf_analysis destructor" << endl;
}


void pdf_analysis::initialize () {
  cout << red_color_bold << "initialization" << default_console_color << endl;

  source.resize(5);
  source[0] = "bs";
  source[1] = "bd";
  source[2] = "comb";
  source[3] = "semi";
  source[4] = "peak";

  ws_ = new RooWorkspace("ws", "ws");
  Mass = new RooRealVar("Mass", "M_{#mu#mu}", 4.90, 5.90, "GeV");
  ws_->import(*Mass);
  massrange_names.resize(4);
  massrange_names[0] = "sb_lo";
  massrange_names[1] = "bd";
  massrange_names[2] = "bs";
  massrange_names[3] = "sb_hi";
  ws_->var("Mass")->setRange(massrange_names[0].c_str(), 4.90, 5.20);
  ws_->var("Mass")->setRange("blind", 5.20, 5.45);
  ws_->var("Mass")->setRange(massrange_names[1].c_str(), 5.20, 5.30);
  ws_->var("Mass")->setRange(massrange_names[2].c_str(), 5.30, 5.45);
  ws_->var("Mass")->setRange(massrange_names[3].c_str(), 5.45, 5.90);
  ws_->var("Mass")->setRange("overall", 4.90, 5.90);

  eta = new RooRealVar("eta", "Candidate pseudorapidity", -2.4, 2.4, "");
  ws_->import(*eta);
  ws_->var("eta")->setRange("eta_all", -2.4, 2.4);
  ws_->var("eta")->setRange("eta_barrel", -1.4, 1.4);
  ws_->var("eta")->setRange("eta_neg_endcap", -2.4, -1.4);
  ws_->var("eta")->setRange("eta_pos_endcap", 1.4, 2.4);

  m1eta = new RooRealVar("m1eta", "first muon pseudorapidity", -2.4, 2.4, "");
  ws_->import(*m1eta);
  m2eta = new RooRealVar("m2eta", "second muon pseudorapidity", -2.4, 2.4, "");
  ws_->import(*m2eta);
  weight = new RooRealVar("weight", "event weight", 0., 1000000.);
  ws_->import(*weight);

  bdt  = new RooRealVar("bdt", "bdt cut", -1., 1.);
  ws_->import(*bdt);
  bdt_0 = new RooRealVar("bdt_0", "bdt", -1., 1.);
  bdt_1 = new RooRealVar("bdt_1", "bdt", -1., 1.);
  bdt_2 = new RooRealVar("bdt_2", "bdt", -1., 1.);
  bdt_3 = new RooRealVar("bdt_3", "bdt", -1., 1.);
  ws_->import(*bdt_0);
  ws_->import(*bdt_1);
  ws_->import(*bdt_2);
  ws_->import(*bdt_3);

  fill_bdt_boundaries();
  channels_cat = new RooCategory("etacat", "eta channels");
  for (int i = 0; i < 4; i++) {
    channels_cat->defineType(Form("etacat_%d", i), i);
  }
  ws_->import(*channels_cat);

  bdt_cat = new RooCategory("bdtcat", "bdt channels");
  for (int i = -1; i < 4; i++) {
    bdt_cat->defineType(Form("bdtcat_%d", i), i);
  }
  ws_->import(*bdt_cat);

  all_cat = new RooCategory("allcat", "channels");
  for (int i = 0; i < super_index(channels, channels_bdt) + 1; i++) {
    all_cat->defineType(Form("allcat_%d", i), i);
  }
  ws_->import(*all_cat);

  RooArgSet cat_set(*ws_->cat("etacat"), *ws_->cat("bdtcat"));
  super_cat = new RooSuperCategory("super_cat", "super_cat", cat_set);
  ws_->import(*super_cat);

  MassRes = new RooRealVar("MassRes", "mass error", 0., 0.2, "GeV");
  ws_->import(*MassRes);

  ReducedMassRes = new RooRealVar("ReducedMassRes", "reduced mass error", 0.0009, 0.045);
  ws_->import(*ReducedMassRes);

//  obs = new RooArgSet(*ws_->var("Mass"), *ws_->var("bdt"), "obs");
  //ws_->import(*obs, RecycleConflictNodes());

  ////////////
  getBFnumbers("input/external_numbers.txt");
  ////////////
}

void pdf_analysis::define_N () {

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
    	RooRealVar N_bs(name("N_bs", i, j), "N_bs", 0, 1000);
    	ws_->import(N_bs);
      RooRealVar N_bd(name("N_bd", i, j), "N_bd", 0, 1000);
      ws_->import(N_bd);
      RooRealVar N_peak(name("N_peak", i, j), "N_peak", 0, 1000);
      ws_->import(N_peak);
      RooRealVar N_semi(name("N_semi", i, j), "N_semi", 0, 10000);
      ws_->import(N_semi);
      RooRealVar N_comb(name("N_comb", i, j), "N_comb", 0, 10000);
      ws_->import(N_comb);
    }
  }
}


void pdf_analysis::define_pdfs (bool semi) {

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      define_bs(i, j);
      define_bd(i, j);
      define_peak(i, j);
      define_semi(i, j, semi);
      define_comb(i, j);
      if (BF_ > 0) {
        define_bf(i, j);
      }
    }
  }
}

void pdf_analysis::fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error, bool hesse, bool setconstant) {

  pdf_name = "pdf_" + pdf;
  if (extended) pdf_name = "pdf_ext_" + pdf;
  
  RooAbsData* subdata;
  if (!simul_bdt_ && !simul_all_) subdata = data->reduce(Form("etacat==etacat::etacat_%d", channel));
  else subdata = data->reduce(Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", channel, channel_bdt));
//  rds_ = subdata;
  string rdh_name = subdata->GetName();
  cout << "**************************" << endl;
  cout << red_color_bold << "fitting " << rdh_name << endl;
  subdata->Print();
  cout << " with " << pdf_name << ":" << endl;
  ws_->pdf( pdf_name.c_str())->Print();
  cout << default_console_color << endl;
  if (!pee) {
    RFR = ws_->pdf( pdf_name.c_str())->fitTo(*subdata, Extended(extended), SumW2Error(sumw2error), NumCPU(2), Hesse(hesse), Save());
    if (print_) print(subdata);
  }
  else {
    RFR = ws_->pdf( pdf_name.c_str())->fitTo(*subdata, ConditionalObservables(*ws_->var("ReducedMassRes")), Extended(extended), SumW2Error(sumw2error), NumCPU(2), Hesse(hesse), Save());
    if (print_) print(subdata, pdf);
  }
  if (setconstant) set_pdf_constant(pdf_name);
}

void pdf_analysis::set_pdf_constant(string name) {
  RooArgSet * set = ws_->pdf(name.c_str())->getVariables();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (!(name == "Mass") && !(name == "Bd_over_Bs") && !(name == "MassRes") && !(name == "ReducedMassRes") && !(name == "bdt") && !(name == "etacat") && !(name == "bdtcat") && !(name == "allcat")
    		&& !(name == "bdt_0") && !(name == "bdt_1") && !(name == "bdt_2") && !(name == "bdt_3") ) {
    	size_t found;
      found = name.find("N_");
      if (found == string::npos) ws_->var(var_Obj->GetName())->setConstant(1);
    }
  }
}

void pdf_analysis::define_bs(int i, int j) {

  RooRealVar Mean_bs(name("Mean_bs", i, j), "Mean_bs", 5.35, 5.30, 5.40);
  RooRealVar Sigma_bs(name("Sigma_bs", i, j), "Sigma_bs", 0.02, 0.005, 0.2);
  RooRealVar Alpha_bs(name("Alpha_bs", i, j), "Alpha_bs", 2.8, 0.1, 3.0);
  RooRealVar Enne_bs(name("Enne_bs", i, j), "Enne_bs", 1., 0., 10.);

  if (!pee) {
    RooRealVar Sigma2_bs(name("Sigma2_bs", i, j), "Sigma2_bs", 0.04, 0.005, 0.2);
    RooRealVar CoeffGauss_bs(name("CoeffGauss_bs", i, j), "CoeffGauss_bs", 0.5, 0., 1.);

    RooGaussian Gau_bs(name("Gau_bs", i, j), "Gau_bs", *ws_->var("Mass"), Mean_bs, Sigma_bs);
    RooCBShape CB_bs(name("CB_bs", i, j), "CB_bs", *ws_->var("Mass"), Mean_bs, Sigma2_bs, Alpha_bs, Enne_bs);
    RooAddPdf pdf_bs_mass(name("pdf_bs_mass", i, j), "pdf_bs_mass", RooArgList(Gau_bs, CB_bs),  CoeffGauss_bs);
    if (!bdt_fit_) {
      RooAddPdf pdf_bs(name("pdf_bs", i, j), "pdf_bs", RooArgList(Gau_bs, CB_bs),  CoeffGauss_bs);
      ws_->import(pdf_bs);
    }
    else {
      RooProdPdf pdf_bs(name("pdf_bs", i, j), "pdf_bs", pdf_bs_mass, *ws_->pdf(name("bdt_pdf_bs", i, j)));
      ws_->import(pdf_bs);
    }
  }
  else {
    RooRealVar PeeK_bs(name("PeeK_bs", i, j), "PeeK_bs", 1., 0.1, 10.);
    RooFormulaVar SigmaRes_bs(name("SigmaRes_bs", i, j), "@0*@1*@2", RooArgList(*ws_->var("ReducedMassRes"), *ws_->var("Mass"), PeeK_bs));
    ws_->import(SigmaRes_bs);
    RooCBShape CB_bs(name("CB_bs", i, j), "CB_bs", *ws_->var("Mass"), Mean_bs, *ws_->function(name("SigmaRes_bs", i, j)), Alpha_bs, Enne_bs);
    if (!bdt_fit_) {
      RooProdPdf pdf_bs (name("pdf_bs", i, j), "pdf_bs", *ws_->pdf(name("ReducedMassRes_pdf_bs", i, j)), Conditional(CB_bs, *ws_->var("Mass")));
      ws_->import(pdf_bs);
    }
    else {
    	RooProdPdf pdf_bs_mass(name("pdf_bs_mass", i, j), "pdf_bs_mass", *ws_->pdf(name("ReducedMassRes_pdf_bs", i, j)), Conditional(CB_bs, *ws_->var("Mass")));
      RooProdPdf pdf_bs(name("pdf_bs", i, j), "pdf_bs", pdf_bs_mass, *ws_->pdf(name("bdt_pdf_bs", i, j)));
      ws_->import(pdf_bs);
    }
  }
}

void pdf_analysis::define_bd(int i, int j) {

  RooRealVar Mean_bd(name("Mean_bd", i, j), "Mean_bd", 5.25, 5.20, 5.30);
  RooRealVar Sigma_bd(name("Sigma_bd", i, j), "Sigma_bd", 0.02, 0.005, 0.2);
  RooRealVar Alpha_bd(name("Alpha_bd", i, j), "Alpha_bd", 2.8, 0.1, 3.0);
  RooRealVar Enne_bd(name("Enne_bd", i, j), "Enne_bd", 1., 0., 10.);

  if (!pee) {
    RooRealVar Sigma2_bd(name("Sigma2_bd", i, j), "Sigma2_bd", 0.04, 0.005, 0.2);
    RooRealVar CoeffGauss_bd(name("CoeffGauss_bd", i, j), "CoeffGauss_bd", 0.5, 0., 1.);

    RooGaussian Gau_bd(name("Gau_bd", i, j), "Gau_bd", *ws_->var("Mass"), Mean_bd, Sigma_bd);
    RooCBShape CB_bd(name("CB_bd", i, j), "CB_bd", *ws_->var("Mass"), Mean_bd, Sigma2_bd, Alpha_bd, Enne_bd);
    if (!bdt_fit_) {
      RooAddPdf pdf_bd(name("pdf_bd", i, j), "pdf_bd", RooArgList(Gau_bd, CB_bd),  CoeffGauss_bd);
      ws_->import(pdf_bd);
    }
    else {
      RooAddPdf pdf_bd_mass(name("pdf_bd_mass", i, j), "pdf_bd_mass", RooArgList(Gau_bd, CB_bd),  CoeffGauss_bd);
      RooProdPdf pdf_bd(name("pdf_bd", i, j), "pdf_bd", pdf_bd_mass, *ws_->pdf(name("bdt_pdf_bd", i, j)));
      ws_->import(pdf_bd);
    }
  }
  else {
    RooRealVar PeeK_bd(name("PeeK_bd", i, j), "PeeK_bd", 1., 0.1, 10.);
    RooFormulaVar SigmaRes_bd(name("SigmaRes_bd", i, j), "@0*@1*@2", RooArgList(*ws_->var("ReducedMassRes"), *ws_->var("Mass"), PeeK_bd));
    ws_->import(SigmaRes_bd);
    RooCBShape CB_bd(name("CB_bd", i, j), "CB_bd", *ws_->var("Mass"), Mean_bd, *ws_->function(name("SigmaRes_bd", i, j)), Alpha_bd, Enne_bd);
    if (!bdt_fit_) {
      RooProdPdf pdf_bd (name("pdf_bd", i, j), "pdf_bd", *ws_->pdf(name("ReducedMassRes_pdf_bd", i, j)), Conditional(CB_bd, *ws_->var("Mass")));
      ws_->import(pdf_bd);
    }
    else {
      RooProdPdf pdf_bd_mass(name("pdf_bd_mass", i, j), "pdf_bd_mass", *ws_->pdf(name("ReducedMassRes_pdf_bd", i, j)), Conditional(CB_bd, *ws_->var("Mass")));
      RooProdPdf pdf_bd(name("pdf_bd", i, j), "pdf_bd", pdf_bd_mass,*ws_->pdf(name("bdt_pdf_bd", i, j)));
      ws_->import(pdf_bd);
    }
  }

  if (SM_) {
    RooConstVar SM_Bd_over_Bs("SM_Bd_over_Bs", "SM_Bd_over_Bs", ratio_);
    ws_->import(SM_Bd_over_Bs);
    RooFormulaVar N_bd_constr(name("N_bd_constr", i, j), "N_bd_constr", "@0*@1", RooArgList(*ws_->var(name("N_bs", i, j)), *(RooConstVar*)ws_->obj("SM_Bd_over_Bs")));
    ws_->import(N_bd_constr);
  }
  else if (bd_constr_) {
    RooRealVar Bd_over_Bs("Bd_over_Bs", "Bd_over_Bs", 0., 10., "");
    ws_->import(Bd_over_Bs);
    RooFormulaVar N_bd_constr(name("N_bd_constr", i, j), "N_bd_constr", "@0*@1", RooArgList(*ws_->var(name("N_bs", i, j)), *ws_->var("Bd_over_Bs")));
    ws_->import(N_bd_constr);
  }
}

void pdf_analysis::define_peak(int i, int j) {

  if (old_tree) {
    RooRealVar Mean_peak(name("Mean_peak", i, j), "Mean_peak", 5.25, 5.20, 5.3);
    RooRealVar Sigma_peak(name("Sigma_peak", i, j), "Sigma_peak", 0.050, 0.01, 0.20);
    RooGaussian pdf_peak(name("pdf_peak", i, j), "pdf_peak", *ws_->var("Mass"), Mean_peak, Sigma_peak);
    ws_->import(pdf_peak);
  }
  else {
    RooRealVar Mean_peak(name("Mean_peak", i, j), "Mean_peak", 5.1, 4.9, 5.4);
    RooRealVar Sigma_peak(name("Sigma_peak", i, j), "Sigma_peak", 0.02, 0.005, 0.2);
    RooRealVar Sigma2_peak(name("Sigma2_peak", i, j), "Sigma2_peak", 0.04, 0.005, 0.2);
    RooRealVar Alpha_peak(name("Alpha_peak", i, j), "Alpha_peak", 2.8, 0., 100.0);
    RooRealVar Enne_peak(name("Enne_peak", i, j), "Enne_peak", 1., 0., 10.);
    RooRealVar CoeffGauss_peak(name("CoeffGauss_peak", i, j), "CoeffGauss_peak", 0.5, 0., 1.);
    RooGaussian Gau_peak(name("Gau_peak", i, j), "Gau_peak", *ws_->var("Mass"), Mean_peak, Sigma_peak);
    RooCBShape CB_peak(name("CB_peak", i, j), "CB_peak", *ws_->var("Mass"), Mean_peak, Sigma2_peak, Alpha_peak, Enne_peak);

    if (!pee) {
      if (!bdt_fit_) {
        RooAddPdf pdf_peak(name("pdf_peak", i, j), "pdf_peak", RooArgList(Gau_peak, CB_peak),  CoeffGauss_peak);
        ws_->import(pdf_peak);
      }
      else {
        RooAddPdf pdf_peak_mass(name("pdf_peak_mass", i, j), "pdf_peak_mass", RooArgList(Gau_peak, CB_peak),  CoeffGauss_peak);
        RooProdPdf pdf_peak(name("pdf_peak", i, j),"pdf_peak",pdf_peak_mass,*ws_->pdf(name("bdt_pdf_peak", i, j)));
        ws_->import(pdf_peak);
      }
    }
    else {
      RooAddPdf mass_peak(name("mass_peak", i, j), "mass_peak", RooArgList(Gau_peak, CB_peak),  CoeffGauss_peak);
      if (!bdt_fit_) {
        RooProdPdf pdf_peak (name("pdf_peak", i, j), "pdf_peak", *ws_->pdf(name("ReducedMassRes_pdf_peak", i, j)), Conditional(mass_peak, *ws_->var("Mass")));
        ws_->import(pdf_peak);
      }
      else {
        RooProdPdf pdf_peak_mass (name("pdf_peak_mass", i, j), "pdf_peak_mass", *ws_->pdf(name("ReducedMassRes_pdf_peak", i, j)), Conditional(mass_peak, *ws_->var("Mass")));
        RooProdPdf pdf_peak(name("pdf_peak", i, j),"pdf_peak",pdf_peak_mass,*ws_->pdf(name("bdt_pdf_peak", i, j)));
        ws_->import(pdf_peak);
      }
    }
  }
}

void pdf_analysis::define_semi(int i, int j, bool semi_keys) {

  if (old_tree) {
    RooRealVar m0_semi(name("m0_semi", i, j), "m0_semi", 5., 6.);
    RooRealVar c_semi(name("c_semi", i, j), "c_semi", 1., 0.1, 20);
    RooRealVar p_semi(name("p_semi", i, j), "p_semi", 0.5, 0.1, 5.);
    RooArgusBG ArgusBG(name("pdf_semi", i, j), "pdf_semi", *ws_->var("Mass"), m0_semi, c_semi, p_semi);
    if (!semi_keys) ws_->import(ArgusBG);
  }
  else {
    RooRealVar C0(name("C0_semi", i, j), "C0", -0.1, -5., 5., "");
    RooRealVar C1(name("C1_semi", i, j), "C1", 0.1, -5., 5., "");
    RooRealVar C2(name("C2_semi", i, j), "C2", -0.1, -5., 5., "");
    RooRealVar C3(name("C3_semi", i, j), "C3", 0.1, -5., 5., "");
    RooRealVar C4(name("C4_semi", i, j), "C4", -0.1, -5., 5., "");
    RooRealVar C5(name("C5_semi", i, j), "C5", 0.1, -5., 5., "");
    RooArgList poly_coeffs(C0, C1, C2/*, C3, C4, C5*/);
    RooChebychev poly(name("poly_semi", i, j), "poly", *ws_->var("Mass"), poly_coeffs);
    RooRealVar tau(name("tau_semi", i, j), "tau", -7,-20.,-0.1, "");
    RooExponential expo(name("expo_semi", i, j), "expo", *ws_->var("Mass"), tau);

    RooRealVar Mean_semi(name("Mean_semi", i, j), "Mean_semi", 5.1, 4.9, 5.5);
    RooRealVar Sigma_semi(name("Sigma_semi", i, j), "Sigma_semi", 0.050, 0.01, 0.50);
    RooGaussian gauss(name("gauss", i, j), "gauss", *ws_->var("Mass"), Mean_semi, Sigma_semi);

    if (!pee) {
      if (!bdt_fit_) {
        RooProdPdf pdf_semi(name("pdf_semi", i, j), "pdf_semi", expo, poly);
        if (!semi_keys) ws_->import(pdf_semi);
      }
      else {
      	if (!semi_keys) {
      		RooProdPdf pdf_semi_mass(name("pdf_semi_mass", i, j), "pdf_semi_mass", expo, poly);
      		RooProdPdf pdf_semi(name("pdf_semi", i, j),"pdf_semi", pdf_semi_mass, *ws_->pdf(name("bdt_pdf_semi", i, j)));
      		ws_->import(pdf_semi);
      	}
      	else {
//      		ws_->pdf(name("pdf_semi", i, j))->SetName(name("pdf_semi_mass", i, j));
      		RooProdPdf pdf_semi(name("pdf_semi", i, j),"pdf_semi", *ws_->pdf(name("keys_semi", i, j)), *ws_->pdf(name("bdt_pdf_semi", i, j)));
      		ws_->import(pdf_semi);
      	}
      }
    }
    else {
      RooProdPdf mass_semi(name("mass_semi", i, j), "pdf_semi", expo, poly);
      if (!bdt_fit_) {
      	if (!semi_keys) {
      		RooProdPdf pdf_semi (name("pdf_semi", i, j), "pdf_semi", *ws_->pdf(name("ReducedMassRes_pdf_semi", i, j)), Conditional(mass_semi, *ws_->var("Mass")));
      		ws_->import(pdf_semi);
      	}
      	else {
//      		ws_->pdf(name("pdf_semi", i, j))->SetName(name("pdf_semi_mass", i, j));
      		RooProdPdf pdf_semi (name("pdf_semi", i, j), "pdf_semi", *ws_->pdf(name("ReducedMassRes_pdf_semi", i, j)), Conditional(*ws_->pdf(name("keys_semi", i, j)), *ws_->var("Mass")));
      		ws_->import(pdf_semi);
      	}
      }
      else {
      	if (!semi_keys) {
      		RooProdPdf pdf_semi_mass (name("pdf_semi_mass", i, j), "pdf_semi_mass", *ws_->pdf(name("ReducedMassRes_pdf_semi", i, j)), Conditional(mass_semi, *ws_->var("Mass")));
      		RooProdPdf pdf_semi(name("pdf_semi", i, j),"pdf_semi",pdf_semi_mass,*ws_->pdf(name("bdt_pdf_semi", i, j)));
      		ws_->import(pdf_semi);
      	}
      	else {
//      		ws_->pdf(name("pdf_semi", i, j))->SetName(name("pdf_semi_masse", i, j));
      		RooProdPdf pdf_semi_mass (name("pdf_semi_mass", i, j), "pdf_semi_mass", *ws_->pdf(name("ReducedMassRes_pdf_semi", i, j)), Conditional(*ws_->pdf(name("keys_semi", i, j)), *ws_->var("Mass")));
      		RooProdPdf pdf_semi(name("pdf_semi", i, j),"pdf_semi", pdf_semi_mass,*ws_->pdf(name("bdt_pdf_semi", i, j)));
      		ws_->import(pdf_semi);
      	}
      }
    }
  }
}

void pdf_analysis::define_comb(int i, int j) {

  RooRealVar exp(name("exp_comb", i, j), "exp_comb", -1., -10., 10);

  if (!pee) {
    if (!bdt_fit_) {
      RooExponential pdf_comb(name("pdf_comb", i, j), "N_comb", *ws_->var("Mass"), exp);
      ws_->import(pdf_comb);
    }
    else {
      RooExponential pdf_comb_mass(name("pdf_comb_mass", i, j), "N_comb", *ws_->var("Mass"), exp);
      RooProdPdf pdf_comb(name("pdf_comb", i, j), "pdf_comb", pdf_comb_mass, *ws_->pdf(name("bdt_pdf_comb", i, j)));
      ws_->import(pdf_comb);
    }
  }
  else  {
    RooExponential mass_comb(name("mass_comb", i, j), "N_comb", *ws_->var("Mass"), exp);
    if (!bdt_fit_) {
      RooProdPdf pdf_comb (name("pdf_comb", i, j), "pdf_comb", *ws_->pdf(name("ReducedMassRes_pdf_comb", i, j)), Conditional(mass_comb, *ws_->var("Mass")));
      ws_->import(pdf_comb);
    }
    else {
      RooProdPdf pdf_comb_mass (name("pdf_comb_mass", i, j), "pdf_comb_mass", *ws_->pdf(name("ReducedMassRes_pdf_comb", i, j)), Conditional(mass_comb, *ws_->var("Mass")));
      RooProdPdf pdf_comb(name("pdf_comb", i, j),"pdf_comb", pdf_comb_mass, *ws_->pdf(name("bdt_pdf_comb", i, j)));
      ws_->import(pdf_comb);
    }
  }
}

void pdf_analysis::define_comb2(int i, int j) {

//  RooRealVar B1(name("B_1comb", i, j), "B_1comb", 1);
	RooRealVar B1(name("B_1comb", i, j), "B_1comb", 0.5, 0. , 1);
//	  RooRealVar B1(name("B_1comb", i, j), "B_1comb", 0., -3.1416, 3.1416);

//  RooRealVar B2(name("B_2comb", i, j), "B_2comb", 0, -100., 100.);
//  RooRealVar B2(name("B_2comb", i, j), "B_2comb", 0., -3.1416, 3.1416);
//  RooRealVar B2(name("B_2comb", i, j), "B_2comb", 0.5, 0., 1.);
  RooFormulaVar B2(name("B_2comb", i, j), "B_2comb", "1.-@0", RooArgList(B1));
//
//  RooFormulaVar bkgp2c(name("bkgp2c", i, j), "", "0. + 0.5*(1.-0.)*(TMath::Sin(@0)+1.0)", RooArgList(B2));
//	RooFormulaVar bkgp1c(name("bkgp1c", i, j), "", "0. + 0.5*(1.-0.)*(TMath::Sin(@0)+1.0)", RooArgList(B1));


//  RooFormulaVar B2(name("B_2comb", i, j), "B_2comb", "1.-@0", RooArgList(bkgp1c));

//  RooRealVar exp2(name("exp_2comb", i, j), "exp_2comb", -1., -10., 10);
//  RooRealVar exp_fract(name("exp_fract_2comb", i, j), "exp_fract_2comb", 0.5, 0., 1.);


  if (!pee) {
//    if (!bdt_fit_) {
//      RooExponential pdf_1comb(name("pdf_1comb", i, j), "pdf_1comb", *ws_->var("Mass"), exp1);
//      RooExponential pdf_2comb(name("pdf_2comb", i, j), "pdf_2comb", *ws_->var("Mass"), exp2);
//      RooAddPdf pdf_comb(name("pdf_12comb", i, j), ("pdf_comb"), RooArgList(pdf_1comb, pdf_2comb), RooArgList(exp_fract));
//      ws_->import(pdf_comb);
//    }
//    else {
//    	RooExponential pdf_1comb(name("pdf_1comb", i, j), "pdf_1comb", *ws_->var("Mass"), exp1);
//    	RooExponential pdf_2comb(name("pdf_2comb", i, j), "pdf_2comb", *ws_->var("Mass"), exp2);
//    	RooAddPdf pdf_comb_mass(name("pdf_comb_mass", i, j), ("pdf_comb"), RooArgList(pdf_1comb, pdf_2comb), RooArgList(exp_fract));
//      RooProdPdf pdf_comb(name("pdf_12comb", i, j), "pdf_comb", pdf_comb_mass, *ws_->pdf(name("bdt_pdf_comb", i, j)));
//      ws_->import(pdf_comb);
//    }
  }
  else  {
//  	RooArgList coeff(B1, bkgp2c);
  	RooArgList coeff(B1, B2);
//  	RooArgList coeff(bkgp1c, B2);

//  	RooPolynomial mass_comb(name("mass_1comb", i, j), "mass_comb", *ws_->var("Mass"), coeff);
//  	RooChebychev mass_comb(name("mass_1comb", i, j), "mass_comb", *ws_->var("Mass"), coeff);
//  	RooBernstein mass_comb("mass_1comb", "mass_comb",  *ws_->var("Mass") ,RooConst(1.), bkgp1positive);
  	RooBernstein mass_comb(name("mass_1comb", i, j), "mass_comb", *ws_->var("Mass"), coeff);
//  	RooAddPdf mass_comb(name("mass_comb", i, j), ("pdf_comb"), RooArgList(pdf_1comb, pdf_2comb), RooArgList(exp_fract));
    if (!bdt_fit_) {
      RooProdPdf pdf_comb (name("pdf_1comb", i, j), "pdf_comb", *ws_->pdf(name("ReducedMassRes_pdf_comb", i, j)), Conditional(mass_comb, *ws_->var("Mass")));
      ws_->import(pdf_comb);
    }
    else {
//      RooProdPdf pdf_comb_mass (name("pdf_comb_mass", i, j), "pdf_comb_mass", *ws_->pdf(name("ReducedMassRes_pdf_comb", i, j)), Conditional(mass_comb, *ws_->var("Mass")));
//      RooProdPdf pdf_comb(name("pdf_12comb", i, j),"pdf_comb", pdf_comb_mass, *ws_->pdf(name("bdt_pdf_comb", i, j)));
//      ws_->import(pdf_comb);
    }
  }
}


void pdf_analysis::define_signals(int i, int j) {

  RooRealVar N_signals(name("N_signals", i, j), "N_signals", 0., 1000);
  ws_->import(N_signals);

  RooRealVar bsfrac_signals(name("bsfrac_signals", i, j), "bsfrac_signals", 0.5, 0.0, 1.0);
  RooAddPdf pdf_signals(name("pdf_signals", i, j), "pdf_signals", RooArgSet(*ws_->pdf(name("pdf_bs", i, j)), *ws_->pdf(name("pdf_bd", i, j)) ), bsfrac_signals);
  ws_->import(pdf_signals);
}

void pdf_analysis::define_rare(int i, int j) {
  
  RooRealVar N_rare(name("N_rare", i, j), "N_rare", 0, 10000);
  ws_->import(N_rare);

  RooRealVar peakfrac_rare(name("peakfrac_rare", i, j), "peakfrac_rare", 0.5, 0.0, 1.0);
  RooAddPdf pdf_rare(name("pdf_rare", i, j), "pdf_rare", RooArgSet(*ws_->pdf(name("pdf_peak", i, j)), *ws_->pdf(name("pdf_semi", i, j)) ), peakfrac_rare);
  ws_->import(pdf_rare);
}

void pdf_analysis::define_rare3(int i, int j) {
  RooRealVar N_expo(name("N_expo3", i, j), "N_expo3", 0, 10000);
  ws_->import(N_expo);
  RooRealVar tau3(name("tau3", i, j), "tau3", -5, -20., -0.01);
  if (!pee) {
    RooExponential pdf_expo(name("pdf_expo3", i, j), "pdf_expo3", *ws_->var("Mass"), tau3);
    ws_->import(pdf_expo);
  }
  else {
    RooExponential expo3(name("expo3", i, j), "expo3", *ws_->var("Mass"), tau3);
    RooProdPdf pdf_expo(name("pdf_expo3", i, j), "pdf_expo3", *ws_->pdf(name("ReducedMassRes_pdf_rare", i, j)), Conditional(expo3, *ws_->var("Mass")));
    ws_->import(pdf_expo);
  }
  //ws_->factory("SUM::pdf_expo(expo1frac_rare[0.,1.]*Exponential::expo1(Mass,Alpha1[-10.,10.]),Exponential::expo2(Mass,Alpha2[-10.,10.]))");
}

void pdf_analysis::define_bf(int i, int j) {

  if (i == 0 && j == 0) {
    RooRealVar BF_bs("BF_bs", "Bs2MuMu branching fraction", Bs2MuMu_SM_BF_val, 0., 1e-8);
    RooRealVar BF_bd("BF_bd", "Bd2MuMu branching fraction", Bd2MuMu_SM_BF_val, 0., 1e-8);
    RooRealVar BF_bsbd("BF_bsbd", "Bs2MuMu / Bd2MuMu branching fraction ratio", Bd2MuMu_SM_BF_val/Bs2MuMu_SM_BF_val, 0., 20.);

    RooRealVar fs_over_fu("fs_over_fu", "fs_over_fu", fs_over_fu_val, max(fs_over_fu_val - 10*fs_over_fu_err, 0.), fs_over_fu_val + 10*fs_over_fu_err);
    RooRealVar one_over_BRBR("one_over_BRBR", "one_over_BRBR", one_over_BRBR_val, max(one_over_BRBR_val - 10*one_over_BRBR_err, 0.), one_over_BRBR_val + 10*one_over_BRBR_err);

    ws_->import(BF_bs);
    ws_->import(BF_bd);
    if (BF_ == 3) ws_->import(BF_bsbd);
    ws_->import(fs_over_fu);
    ws_->import(one_over_BRBR);
  }

  RooRealVar N_bu(name("N_bu", i, j), "N_bu", N_bu_val[i][j], max(N_bu_val[i][j] - 10*N_bu_err[i][j], 0.), N_bu_val[i][j] + 10*N_bu_err[i][j]);
  ws_->import(N_bu);

  RooRealVar effratio_bs(name("effratio_bs", i, j), "effratio_bs", effratio_bs_val[i][j], max(effratio_bs_val[i][j] - 10*effratio_bs_err[i][j], 0.), effratio_bs_val[i][j] + 10*effratio_bs_err[i][j]);
  RooRealVar effratio_bd(name("effratio_bd", i, j), "effratio_bd", effratio_bd_val[i][j], max(effratio_bd_val[i][j] - 10*effratio_bd_err[i][j], 0.), effratio_bd_val[i][j] + 10*effratio_bd_err[i][j]);

  RooFormulaVar N_bs_constr(name("N_bs_formula", i, j), "N_bs(i) = BF * K(i)",  "@0*@1*@2*@3*@4", RooArgList(*ws_->var("BF_bs"), *ws_->var(name("N_bu", i, j)), *ws_->var("fs_over_fu"), effratio_bs, *ws_->var("one_over_BRBR")));
  RooFormulaVar N_bsbd_constr(name("N_bs_formula", i, j), "N_bs(i) = BF * K(i)",  "@0*@1*@2*@3*@4*@5", RooArgList( *ws_->var("BF_bd"), *ws_->var("BF_bsbd"), *ws_->var(name("N_bu", i, j)), *ws_->var("fs_over_fu"), effratio_bs, *ws_->var("one_over_BRBR")));

  RooFormulaVar N_bd_constr(name("N_bd_formula", i, j), "N_bd(i) = BF * K(i)",  "@0*@1*@2*@3", RooArgList( *ws_->var("BF_bd"), *ws_->var(name("N_bu", i, j)), effratio_bd, *ws_->var("one_over_BRBR")));

  cout << "channel: " << i << ";" << j << " expected Bs = " << N_bs_constr.getVal() << endl;
  cout << "channel: " << i << ";" << j << " expected Bd = " << N_bd_constr.getVal() << endl;

  if (BF_ != 3) ws_->import(N_bs_constr);
  else ws_->import(N_bsbd_constr);
  ws_->import(N_bd_constr);
}

string pdf_analysis::define_pdf_sum(string name) {

  vector <string> pdfs;
  size_t found;
  found = name.find("bs");
  if (found != string::npos) pdfs.push_back("bs");
  found = name.find("bd");
  if (found != string::npos) pdfs.push_back("bd");
  found = name.find("rare");
  if (found != string::npos) pdfs.push_back("rare");
  found = name.find("comb");
  if (found != string::npos) pdfs.push_back("comb");
  found = name.find("semi");
  if (found != string::npos) pdfs.push_back("semi");
  found = name.find("hist");
  if (found != string::npos) pdfs.push_back("hist");
  found = name.find("expo3");
  if (found != string::npos) pdfs.push_back("expo3");
  found = name.find("peak");
  if (found != string::npos) pdfs.push_back("peak");

  string pdf_sum = "SUM::pdf_ext_";
  string title;
  for (unsigned int i = 0; i < pdfs.size(); i++) {
    pdf_sum += pdfs[i];
    title += pdfs[i];
  }
  pdf_sum += "(";
  for (unsigned int i = 0; i < pdfs.size(); i++) {
    pdf_sum += "N_";
    if (pdfs[i]=="hist" || pdfs[i]=="expo3") pdf_sum += "semi";
    else pdf_sum += pdfs[i];
    pdf_sum += "*pdf_";
    pdf_sum += pdfs[i];
    if (i != pdfs.size() -1) pdf_sum += ",";
  }
  pdf_sum += ")";
  cout << "formed pdf: " << pdf_sum << endl;
  ws_->factory(pdf_sum.c_str());
  return (title);
}

void pdf_analysis::print(RooAbsData* data, string output) {
  int colors[11] = {632, 400, 616, 432, 800, 416, 820, 840, 860, 880, 900};
  RooAbsData* subdata_res = data;
//  if (pee) subdata_res = (RooAbsData*)ws_->data(Form("MassRes_rdh_%s", output.c_str()))->Clone();
  RooPlot *rp = ws_->var("Mass")->frame();
  RooPlot *rp_bdt;
  if (bdt_fit_) rp_bdt = ws_->var(name("bdt", channel))->frame();
  data->plotOn(rp, Binning(40));
  if (bdt_fit_) data->plotOn(rp_bdt, Binning(100));
  if (!pee) {
    ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue));
    if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue));
  }
  else {
    ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *subdata_res, kFALSE));
    if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *subdata_res, kFALSE));

    TCanvas* cetad_mres = new TCanvas("cetad_mres", "cetad_mres", 600, 600);
    RooPlot *rp_res = ws_->var("ReducedMassRes")->frame();
    subdata_res->plotOn(rp_res);
    ws_->pdf(pdf_name.c_str())->plotOn(rp_res);
    rp_res->Draw();
    cetad_mres->Print( (get_address("MassEta", pdf_name) + ".gif").c_str());
    cetad_mres->Print( (get_address("MassEta", pdf_name) + ".pdf").c_str());
    delete rp_res;
    delete cetad_mres;
  }
  if(!no_legend) {
    ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
//    if (bdt_fit_) ws_->pdf(pdf_name.c_str())->paramOn(rp_bdt, Layout(0.0, 0.4, 0.9));
  }
  
  //components
  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (i > 11) i = 0;
      size_t found1 = pdf_name.find("total");
      //cout << name << endl;
      if (found1 == string::npos) {
        if (!pee) {
          ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
        }
        else {
          size_t found2 = pdf_name.find("SigmaRes");
          if (found2 == string::npos) {
            ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *subdata_res, kFALSE), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
          }
        }
      }
      else {
        if (!pee) {
          if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_comb" || name=="pdf_semi" || name=="pdf_peak") {
            ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
            if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
          }
        }
        else {
          size_t found2 = pdf_name.find("SigmaRes");
          if (found2 == string::npos) {
            if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_comb" || name=="pdf_semi" || name=="pdf_peak") {
              ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *subdata_res, kFALSE), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
              if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *subdata_res, kFALSE), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
            }
          }
        }
      }
      i++;
    }
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);  
  rp->Draw();
  c->Print((get_address(pdf_name, "") + ".pdf").c_str());
  c->Print((get_address(pdf_name, "") + ".gif").c_str());
  delete rp;
  delete c;

  if (bdt_fit_) {
    TCanvas* c_bdt = new TCanvas("c_bdt", "c_bdt", 600, 600);
    rp_bdt->Draw();
    c_bdt->Print((get_address("BDT", pdf_name) + ".pdf").c_str());
    c_bdt->Print((get_address("BDT", pdf_name) + ".gif").c_str());
    delete c_bdt;
  }
  if (bdt_fit_)	 delete rp_bdt;

  return;
}

void pdf_analysis::simsplit() {
  cout << "simsplitting" << endl;
  simul_ = true;
  ostringstream splitter;
  ostringstream splitted;
  splitter << "SIMCLONE::pdf_ext_simul(pdf_ext_total, $SplitParam({";
  RooArgSet* allvars = ws_->pdf("pdf_ext_total")->getVariables();
  TIterator* it = allvars->createIterator();
  RooRealVar* var_Obj = 0;
  int cont = 0;
  while ( (var_Obj = (RooRealVar*)it->Next()) ) {
    string name = var_Obj->GetName();
    if ( !(name == "Mass") && !(name == "Bd_over_Bs") && !(name == "SM_Bd_over_Bs") && !(name == "bdt") && !(name == "MassRes") && !(name == "ReducedMassRes")
    		&& !(name == "bdt_0") && !(name == "bdt_1")	&& !(name == "bdt_2")	&& !(name == "bdt_3")) {
      if (cont != 0) {
        splitter << ", ";
        splitted << ", ";
      }
      splitter << name;
      splitted << name;
      cont++;
    }
  }
  /// splitter << ",etapdf_bs"; /// warning
  splitter << "}, {etacat[";
  for (int i = 0; i < 2; i++) {
    if (i != 0) splitter << ",";
    splitter << "etacat_" << i ;
  }
  splitter << "],bdtcat[";
  for (int i = 0; i < 3; i++) {
    if (i != 0) splitter << ",";
    splitter << "bdtcat_" << i ;
  }
  splitter << "]}))";
  cout << splitter.str() << endl;
  ws_->factory(splitter.str().c_str());
  //RooSimWSTool sct(*ws_);

  ws_->Print();
  cout << splitted.str() << endl;
  //RooSimultaneous* model_sim2 = sct.build("pdf_ext_simul", "pdf_ext_total", SplitParam(splitted.str().c_str(), "etacat,bdtcat"));
  //model_sim2->Print();
}

double pdf_analysis::getErrorHigh(RooRealVar *var) {
  double value = var->getVal();
  if (value == var->getMax()) return 0;
  double error = var->getErrorHi();
  if (error == 0) return var->getError();
  else return error;
}

double pdf_analysis::getErrorLow(RooRealVar *var) {
  double value = var->getVal();
  if (value == var->getMin()) return 0;
  double error = var->getErrorLo();
  if (error == 0) {
    if (value - var->getError() < 0.0) return (-1 * value);
    else return (-1*var->getError());
  }
  else if (value + error < 0.0) {
    return -value;
  }
  return error;
}

TH1D* pdf_analysis::define_massRes_pdf(RooDataSet *rds, string name, bool rkeys) {
	cout << endl;
	ostringstream suffix;
	suffix << name;
	if (simul_) suffix << "_" << channel;
	if (simul_bdt_ || simul_all_) suffix << "_" << channel_bdt;

  RooArgList varlist(*ReducedMassRes, *weight);
  RooDataSet* subdata_bdt = new RooDataSet("subdata_bdt", "subdata_bdt", varlist, "weight");
  const RooArgSet* aRow;
  string name_h("ReducedMassRes_h_");
  name_h = name_h + suffix.str();
  TH1D histo(name_h.c_str(), rds->GetTitle(), 200, 0.0009, 0.045);
  TH1D * histo_fine = new TH1D(Form("%s_ReducedMassRes_%s_h", name.c_str(), suffix.str().c_str()), rds->GetTitle(), 1000, 0., 0.2); // no filled
  for (Int_t j = 0; j < rds->numEntries(); j++) {
    aRow = rds->get(j);
    RooRealVar* reducedmassres = (RooRealVar*)aRow->find("ReducedMassRes");
    double Weight = rds->weight();
    RooArgSet varlist_tmp_res(*reducedmassres);
    if (aRow->getCatIndex("etacat") == channel) {
      if ((!simul_bdt_ && !simul_all_) || aRow->getCatIndex("bdtcat") == channel_bdt || name == "comb") {
      	subdata_bdt->add(varlist_tmp_res, Weight);
        histo.Fill(reducedmassres->getVal(), Weight);
      }
    }
  }
  cout << name << " reduced resolution entries = " <<  histo.GetEntries() << endl;
  string name_rdh("ReducedMassRes_rdh_");
  name_rdh = name_rdh + suffix.str();
  RooDataHist *MassRes_rdh = new RooDataHist(name_rdh.c_str(), name_rdh.c_str(), *ws_->var("ReducedMassRes"), &histo);
  string name_pdf("ReducedMassRes_pdf_");
  name_pdf = name_pdf + suffix.str();
  RooHistPdf * MassRes_rhpdf = new RooHistPdf(name_pdf.c_str(), name_pdf.c_str(), RooArgList(*ws_->var("ReducedMassRes")), *MassRes_rdh);
  if (!rkeys) {
  	ws_->import(*MassRes_rhpdf);
  	print_pdf(MassRes_rhpdf, ws_->var("ReducedMassRes"));
  	return histo_fine;
  }

  if (subdata_bdt->numEntries() < 20000) {
  }
  else {
  	cout << "too many entries for RooKeysPdf, switching to random-generated distribution with 20000 events";
  	subdata_bdt = MassRes_rhpdf->generate(*ws_->var("ReducedMassRes"), 20000, Extended(false), AutoBinned(false));
  	cout << " ... generated" << endl;
  }
  RooKeysPdf *kest = new RooKeysPdf(name_pdf.c_str(), name_pdf.c_str(), *ws_->var("ReducedMassRes"), *subdata_bdt);
  ws_->import(*kest);
  print_pdf(kest, ws_->var("ReducedMassRes"));
  return histo_fine;
}

TH1D* pdf_analysis::define_semi_pdf(RooDataSet *rds, string name, bool rkeys) {
	cout << endl;
	ostringstream suffix;
	suffix << name;
	if (simul_) suffix << "_" << channel;
	if (simul_bdt_ || simul_all_) suffix << "_" << channel_bdt;

  RooArgList varlist(*Mass, *weight);
  RooDataSet* subdata_mass = new RooDataSet("subdata_mass", "subdata_mass", varlist, "weight");
  const RooArgSet* aRow;
  string name_h("Mass_h_");
  name_h = name_h + suffix.str();
  TH1D histo(name_h.c_str(), rds->GetTitle(), 200, 4.9, 5.9);
  TH1D * histo_fine = new TH1D(Form("%s_mass_%s_h", name.c_str(), suffix.str().c_str()), rds->GetTitle(), 1000, 4.9, 5.9); // no filled
  for (Int_t j = 0; j < rds->numEntries(); j++) {
    aRow = rds->get(j);
    RooRealVar* mass = (RooRealVar*)aRow->find("Mass");
    double Weight = rds->weight();
    RooArgSet varlist_tmp_mass(*mass);
    if (aRow->getCatIndex("etacat") == channel) {
      if ((!simul_bdt_ && !simul_all_) || aRow->getCatIndex("bdtcat") == channel_bdt) {
      	subdata_mass->add(varlist_tmp_mass, Weight);
        histo.Fill(mass->getVal(), Weight);
      }
    }
  }
  cout << name << " entries = " <<  histo.GetEntries() << endl;
  string name_rdh("Mass_rdh_");
  name_rdh = name_rdh + suffix.str();
  RooDataHist *Mass_rdh = new RooDataHist(name_rdh.c_str(), name_rdh.c_str(), *ws_->var("Mass"), &histo);
  string name_pdf("keys_");
  name_pdf = name_pdf + suffix.str();
  RooHistPdf * Mass_rhpdf = new RooHistPdf(name_pdf.c_str(), name_pdf.c_str(), RooArgList(*ws_->var("Mass")), *Mass_rdh);
  if (!rkeys) {
  	ws_->import(*Mass_rhpdf);
  	print_pdf(Mass_rhpdf, ws_->var("Mass"));
  	return histo_fine;
  }

  if (subdata_mass->numEntries() < 20000) {
  }
  else {
  	cout << "too many entries for RooKeysPdf, switching to random-generated distribution with 20000 events";
  	subdata_mass = Mass_rhpdf->generate(*ws_->var("Mass"), 20000, Extended(false), AutoBinned(false));
  	cout << " ... generated" << endl;
  }

  RooKeysPdf *kest = new RooKeysPdf(name_pdf.c_str(), name_pdf.c_str(), *ws_->var("Mass"), *subdata_mass, RooKeysPdf::MirrorBoth, 2.5) ;
  ws_->import(*kest);
  print_pdf(kest, ws_->var("Mass"));
  return histo_fine;
}

TH1D* pdf_analysis::define_bdt_pdf(RooDataSet *rds, string name, TFile* bdt_syst_f, bool rkeys, Double_t bdt_min) {
	cout << endl;
	ostringstream suffix;
	suffix << name;
	if (simul_) suffix << "_" << channel;
	if (simul_bdt_ || simul_all_) suffix << "_" << channel_bdt;
  RooArgList varlist(*ws_->var(pdf_analysis::name("bdt", channel)), *weight);
  RooDataSet* subdata_bdt = new RooDataSet("subdata_bdt", "subdata_bdt", varlist, "weight");
  const RooArgSet* aRow;
  TH1D histo(rds->GetTitle(), rds->GetTitle(), 200, ws_->var(pdf_analysis::name("bdt", channel))->getMin(), 1);
  TH1D * histo_fine = new TH1D(Form("bdt_%s_h", suffix.str().c_str()), rds->GetTitle(), 1000, ws_->var(pdf_analysis::name("bdt", channel))->getMin(), 1);
  for (Int_t j = 0; j < rds->numEntries(); j++) {
    aRow = rds->get(j);
    RooRealVar* BDT = (RooRealVar*)aRow->find(pdf_analysis::name("bdt", channel));
    double Weight = rds->weight();
    RooArgSet varlist_tmp_bdt(*BDT);
    if (aRow->getCatIndex("etacat") == channel) {
      if ((!simul_bdt_ && !simul_all_) || aRow->getCatIndex("bdtcat") == channel_bdt) {
      	subdata_bdt->add(varlist_tmp_bdt, Weight);
        histo.Fill(BDT->getVal(), Weight);
        histo_fine->Fill(BDT->getVal(), Weight);
      }
    }
  }
  cout << name  << " bdt entries = " <<  histo.GetEntries() << endl;

  string name_rdh("bdt_rdh_");
  name_rdh = name_rdh + suffix.str();
  RooDataHist *bdt_rdh = new RooDataHist(name_rdh.c_str(), name_rdh.c_str(), *ws_->var(pdf_analysis::name("bdt", channel)), &histo);
  string name_pdf("bdt_pdf_");
  name_pdf = name_pdf + suffix.str();
  RooHistPdf * bdt_rhpdf = new RooHistPdf(name_pdf.c_str(), name_pdf.c_str(), RooArgList(*ws_->var(pdf_analysis::name("bdt", channel))), *bdt_rdh);

  if (!rkeys) {
  	ws_->import(*bdt_rhpdf);
  	print_pdf(bdt_rhpdf, ws_->var(pdf_analysis::name("bdt", channel)));
  	return histo_fine;
  }
  else {
  	string name_rhpdf("bdt_rhpdf_");
  	name_rhpdf += suffix.str();
  	bdt_rhpdf->SetName(name_rhpdf.c_str());
  	if (subdata_bdt->numEntries() > 20000) {
  		cout << "too many entries for RooKeysPdf, switching to random-generated distribution with 20000 events";
  		subdata_bdt = bdt_rhpdf->generate(*ws_->var(pdf_analysis::name("bdt", channel)), 20000, Extended(false), AutoBinned(false));
  		cout << " ... generated" << endl;
  	}
  	string kest1_name("bdt_kest_");
  	kest1_name += suffix.str();
  	RooKeysPdf kest1(kest1_name.c_str(), kest1_name.c_str(), *ws_->var(pdf_analysis::name("bdt", channel)), *subdata_bdt, RooKeysPdf::MirrorBoth, 2);

  	if (bdt_syst_f) {
  		string name_rhpdf("bdt_rhpdf_");
  		name_rhpdf += suffix.str();
  		string bdt_syst_histo_name("bdt_div_jpsiK");
  		if (name.compare("comb") == 0) bdt_syst_histo_name = "bdt_div_5_1";
  		TH1D * bdt_syst_histo = (TH1D*)bdt_syst_f->Get(bdt_syst_histo_name.c_str());
  		TH1D bdt_biased_histo("bdt_biased_histo", "bdt_biased_histo", histo.GetNbinsX(), bdt_min, histo.GetXaxis()->GetXmax());
  		RooDataSet* subdata_biased_bdt = new RooDataSet("subdata_biased_bdt", "subdata_biased_bdt", varlist, "weight");
  		for (int i = 1; i <= bdt_syst_histo->GetNbinsX(); i++) {
  			double biased_value = (bdt_syst_histo->GetBinContent(i)) * (histo.GetBinContent(i));
  			bdt_biased_histo.SetBinContent(i, biased_value);
  			if (biased_value < 0) bdt_biased_histo.SetBinContent(i, 0);
  		}
  		if (subdata_bdt->numEntries() < 20000) {
  			for (int i = 0; i < subdata_bdt->numEntries(); i++) {
  				aRow = subdata_bdt->get(i);
  				RooRealVar* BDT = (RooRealVar*)aRow->find(pdf_analysis::name("bdt", channel));
  				double bdt_val = BDT->getVal();
  				int bin = bdt_syst_histo->FindBin(bdt_val);
  				float correction = bdt_syst_histo->GetBinContent(bin);
  				double Weight = subdata_bdt->weight() * correction;
  				RooArgSet varlist_tmp_bdt(*BDT);
  				subdata_biased_bdt->add(varlist_tmp_bdt, Weight);
  			}
  		}
  		else {
  			string name_biased_rdh("bdt_biased_rdh_");
  			name_biased_rdh += suffix.str();
  			RooDataHist *bdt_biased_rdh = new RooDataHist(name_biased_rdh.c_str(), name_biased_rdh.c_str(), *ws_->var(pdf_analysis::name("bdt", channel)), &bdt_biased_histo);
  			string name_biased_rhpdf("bdt_biased_rhpdf_");
  			name_biased_rhpdf += suffix.str();
  			RooHistPdf * bdt_biased_rhpdf = new RooHistPdf(name_biased_rhpdf.c_str(), name_biased_rhpdf.c_str(), RooArgList(*ws_->var(pdf_analysis::name("bdt", channel))), *bdt_biased_rdh);
  			subdata_biased_bdt = bdt_biased_rhpdf->generate(*ws_->var(pdf_analysis::name("bdt", channel)), 10000, Extended(false), AutoBinned(false));
  		}
  		string kest2_name("bdt_kest_biased_");
  		kest2_name += suffix.str();
  		RooKeysPdf kest2(kest2_name.c_str(), kest2_name.c_str(), *ws_->var(pdf_analysis::name("bdt", channel)), *subdata_biased_bdt, RooKeysPdf::MirrorLeft, 1.5) ;
  		string beta_name("beta_");
  		beta_name += suffix.str();
  		RooRealVar beta(beta_name.c_str(), beta_name.c_str(), 0., 0., 1.);
  		ws_->var(pdf_analysis::name("bdt", channel))->setBins(1000, "cache");
  		beta.setBins(100, "cache");
  		beta.setConstant(true);
  		RooIntegralMorph * bdt_rim = new RooIntegralMorph(name_pdf.c_str(), name_pdf.c_str(), kest1, kest2, *ws_->var(pdf_analysis::name("bdt", channel)), beta, true);
  		ws_->import(*bdt_rim);
  		print_pdf(bdt_rim, ws_->var(pdf_analysis::name("bdt", channel)));
  		return histo_fine;
  	}
  	else {
  		kest1.SetName(name_pdf.c_str());
  		ws_->import(kest1);
  		print_pdf(&kest1, ws_->var(pdf_analysis::name("bdt", channel)));
  		return histo_fine;
  	}
  }
}

void pdf_analysis::print_pdf(RooAbsPdf* pdf, RooRealVar * var) {
  RooPlot *rp = var->frame();
  pdf->plotOn(rp);
  pdf->paramOn(rp, Layout(0, 0.40, 0.90), ShowConstants(kTRUE));
  TCanvas canvas("canvas", "canvas", 600, 600);
//  canvas.SetLogy();
  rp->Draw();
  canvas.Print((get_address("distro", pdf->GetName(), kTRUE) + ".gif").c_str());
  canvas.Print((get_address("distro", pdf->GetName(), kTRUE) + ".pdf").c_str());
  delete rp;
}

const char* pdf_analysis::name(string name, int i, int j) {
//	if (name == "bdt") return name.c_str();
  if (!simul_) return name.c_str();
  else {
    if (!simul_bdt_ && !simul_all_) return Form("%s_%d", name.c_str(), i);
    else return Form("%s_%d_%d", name.c_str(), i, j);
  }
}

void pdf_analysis::set_bkg_normalization(string input) {
  FILE *estimate_file = fopen(input.c_str(), "r");
  if (!estimate_file) {cout << "file " << input << " does not exist" << endl; exit(1);}
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      char buffer[1024];
      char bkg_type[128];
      while (fgets(buffer, sizeof(buffer), estimate_file)) {
        if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
        if (buffer[0] == '#') continue;
        double value, error;
        string scan("%s\t%lf\t%lf");
//        if (simul_all_ || simul_bdt_) scan = "%s\t%lf\t%lf";
        sscanf(buffer, "%s\t%lf\t%lf", bkg_type, &value, &error);
        int index = i;
        if (!simul_) i = channel;


        if ( (!strcmp(bkg_type, name("N_peak", index, j))))  {
        	if ((simul_bdt_ || simul_all_ || bdt_fit_) && !final_) {
        		value *= peak_bdt_factor[i][j]; /// factor
        		error *= peak_bdt_factor[i][j]; /// factor
        	}
        	ws_->var(name("N_peak", i, j))->setVal(value);
        	ws_->var(name("N_peak", i, j))->setError(error);
        	cout << name("N_peak", i, j) << " set val to " << value << "; set error to " << error << endl;
        }
        if ( (!strcmp(bkg_type, name("N_semi", index, j))))  {
        	if ((simul_bdt_ || simul_all_ || bdt_fit_) && !final_) {
          	value *= semi_bdt_factor[i][j]; /// factor
          	error *= semi_bdt_factor[i][j]; /// factor
          }
          ws_->var(name("N_semi", i, j))->setVal(value);
          ws_->var(name("N_semi", i, j))->setError(error);
          cout << name("N_semi", i, j) << " set val to " << value << "; set error to " << error << endl;
        }
        if ( (!strcmp(bkg_type, name("N_comb", index, j))))  {
        	if ((simul_bdt_ || simul_all_ || bdt_fit_) && !final_) {
        		value *= comb_bdt_factor[i][j]; /// factor
        		error *= comb_bdt_factor[i][j]; /// factor
        	}
        	ws_->var(name("N_comb", i, j))->setVal(value);
        	ws_->var(name("N_comb", i, j))->setError(error);
        	cout << name("N_comb", i, j) << " set val to " << value << "; set error to " << error << endl;
        }
      }
      rewind(estimate_file);
    }
  }
  if (estimate_file) fclose(estimate_file);
}

void pdf_analysis::gen_and_fit(string pdfname) {
  pdf_name = pdfname;
  RooAbsPdf* pdf = ws_->pdf(pdf_name.c_str());
  pdf->Print();

  RooDataSet* protodata = (RooDataSet*)ws_->data(name("ReducedMassRes_rdh_bs", channel));
  protodata->Print();
  RooPlot* rp0 = ws_->var("ReducedMassRes")->frame();
  protodata->plotOn(rp0);
  TCanvas* canvas0 = new TCanvas("canvas0", "canvas0", 600, 600);
  rp0->Draw();
  canvas0->Print((get_address("sample0") + ".gif").c_str());
  canvas0->Print((get_address("sample0") + ".pdf").c_str());

  delete rp0;
  delete canvas0;

  RooDataSet* data = pdf->generate( RooArgSet(*ws_->var("Mass"), *ws_->var("ReducedMassRes"), *ws_->var("bdt_0")), 50);
  data->Print();
  RooPlot* rp_31 = ws_->var("Mass")->frame();
  RooPlot* rp_32 = ws_->var("ReducedMassRes")->frame();
  RooPlot* rp_33 = ws_->var("bdt_0")->frame();
  TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 1800, 600);
  canvas3->Divide(3,1);
  data->plotOn(rp_31);
  data->plotOn(rp_32);
  data->plotOn(rp_33);
  canvas3->cd(1);
  rp_31->Draw();
  canvas3->cd(2);
  rp_32->Draw();
  canvas3->cd(3);
  rp_33->Draw();
  canvas3->Print((get_address("dataa") + ".gif").c_str());
  canvas3->Print((get_address("dataa") + ".pdf").c_str());
  delete canvas3;
  delete rp_31;
  delete rp_32;
  delete rp_33;

  if (pee) pdf->fitTo(*data, ConditionalObservables(*ws_->var("ReducedMassRes")));
  else pdf->fitTo(*data);

  RooPlot* rp = ws_->var("Mass")->frame();
  data->plotOn(rp);
  pdf->paramOn(rp);
  if (pee) {
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *data, kFALSE), LineColor(kOrange));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *data, kFALSE), Components("pdf_bs"), LineColor(kRed));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *data, kFALSE), Components("pdf_bd"), LineColor(kBlue));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *data, kFALSE), Components("pdf_rare"), LineColor(kGreen));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("ReducedMassRes")), *data, kFALSE), Components("pdf_comb"), LineColor(kCyan));
  }
  else {
    pdf->plotOn(rp, LineColor(kOrange));
    pdf->plotOn(rp, Components("pdf_bs"), LineColor(kRed));
    pdf->plotOn(rp, Components("pdf_bd"), LineColor(kBlue));
    pdf->plotOn(rp, Components("pdf_comb"), LineColor(kCyan));
    pdf->plotOn(rp, Components("pdf_semi"), LineColor(kGreen));
    pdf->plotOn(rp, Components("pdf_peak"), LineColor(kGreen));
  }
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  rp->Draw();
  canvas->Print((get_address("sample") + ".gif").c_str());
  canvas->Print((get_address("sample") + ".pdf").c_str());
  delete rp;
  delete canvas;

  if (bdt_fit_) {
    RooPlot* rp_bdt = ws_->var("bdt_0")->frame();
    data->plotOn(rp_bdt);
    pdf->plotOn(rp_bdt, LineColor(kOrange));
    pdf->plotOn(rp_bdt, Components("pdf_bs"), LineColor(kRed));
    pdf->plotOn(rp_bdt, Components("pdf_bd"), LineColor(kBlue));
    pdf->plotOn(rp_bdt, Components("pdf_comb"), LineColor(kCyan));
    pdf->plotOn(rp_bdt, Components("pdf_semi"), LineColor(kGreen));
    pdf->plotOn(rp_bdt, Components("pdf_peak"), LineColor(kGreen));
    pdf->paramOn(rp_bdt);
    TCanvas* canvas_bdt = new TCanvas("canvas_bdt", "canvas_bdt", 600, 600);
    rp_bdt->Draw();
    canvas_bdt->Print((get_address("BDT_sample") + ".gif").c_str());
    canvas_bdt->Print((get_address("BDT_sample") + ".pdf").c_str());
    delete rp_bdt;
    delete canvas_bdt;
  }

  if (pee) {
    TH1* mass_eta_h = pdf->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var("MassRes"), Binning(30))) ;
    TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 600, 600);
    mass_eta_h->Draw("surf");
    canvas1->Print((get_address("sample1_signals") + ".gif").c_str());
    canvas1->Print((get_address("sample1_signals") + ".pdf").c_str());
    delete canvas1;
  }
}

string pdf_analysis::get_address(string name, string pdf, bool channeling) {
  ostringstream address;
  address <<  "fig/" << name;
  if (pdf == "") address << "_" << meth_;
  else address << "_" << pdf << "_" << meth_;
  if (simul_) {
    address  << "_simul_";
    if (channeling) address << channel;
  }
  if (!simul_) address  << "_" << ch_s_;
  if (simul_ && (simul_bdt_ || simul_all_)) {
    address << "_simulBdt_";
    if (channeling) address << channel_bdt;
  }
  if (BF_>0) address << "_BF" << BF_;
  if (SM_) address << "_SM";
  if (bd_constr_) address << "_BdConst";
  if (pee) address << "_PEE";
  if (bdt_fit_) address << "_2D";
  if (syst && !randomsyst) address << "_syst";
  else if (syst && randomsyst) address << "_randomsyst";
  return address.str();
}

void pdf_analysis::getBFnumbers(string numbers_filename) {
  eff_bd_val.resize(channels);
  eff_bs_val.resize(channels);
  eff_bu_val.resize(channels);
  N_bu_val.resize(channels);
  eff_bd_err.resize(channels);
  eff_bs_err.resize(channels);
  eff_bu_err.resize(channels);
  N_bu_err.resize(channels);
  eff_rel_err.resize(channels);

  if (simul_bdt_ || simul_all_ || bdt_fit_) {
  	bs_bdt_factor.resize(channels);
  	bd_bdt_factor.resize(channels);
  	bu_bdt_factor.resize(channels);
  	peak_bdt_factor.resize(channels);
  	semi_bdt_factor.resize(channels);
  	comb_bdt_factor.resize(channels);
  }
  for (int i = 0; i < channels; i++) {
  	int size = 1;
  	if (simul_bdt_ || simul_all_) size = bdt_boundaries[i].size() - 1;
  	eff_bd_val[i].resize(size);
  	eff_bs_val[i].resize(size);
  	eff_bu_val[i].resize(size);
  	N_bu_val[i].resize(size);
  	eff_bd_err[i].resize(size);
  	eff_bs_err[i].resize(size);
  	eff_bu_err[i].resize(size);
  	N_bu_err[i].resize(size);
  	eff_rel_err[i].resize(size);
  	if (simul_bdt_ || simul_all_ || bdt_fit_) {
  		bs_bdt_factor[i].resize(size);
  		bd_bdt_factor[i].resize(size);
  		bu_bdt_factor[i].resize(size);
  		peak_bdt_factor[i].resize(size);
  		semi_bdt_factor[i].resize(size);
  		comb_bdt_factor[i].resize(size);
  	}
  }
  acceptance_sys.resize(2);
  mass_scale_sys.resize(2);
  kaon_track_sys.resize(2);
  trigger_sys.resize(2);
  muon_id_sys.resize(2);

  parse_external_numbers(numbers_filename);
  parse_efficiency_numbers();
  if (channels == 4) parse_efficiency_numbers(2);
  if ((simul_bdt_ || simul_all_ || bdt_fit_) && !final_) bdt_effs();

  cout << red_color_bold << "channels = " << channels << "; channels_bdt = " << channels_bdt << "; channels_all = " << channels_all << default_console_color << endl;

  one_over_BRBR_val = 1./ (Bu2JpsiK_BF_val * Jpsi2MuMu_BF_val);
  one_over_BRBR_err = one_over_BRBR_val * sqrt(pow(Bu2JpsiK_BF_err/Bu2JpsiK_BF_val, 2) + pow(Jpsi2MuMu_BF_err/Jpsi2MuMu_BF_val, 2));
  cout << "one_over_BRBR = " << one_over_BRBR_val << " +/- " << one_over_BRBR_err << endl;
  effratio_bs_val.resize(channels);
  effratio_bd_val.resize(channels);
  effratio_bs_err.resize(channels);
  effratio_bd_err.resize(channels);
  for ( int i = 0; i < channels; i++) {
  	int size = 1;
  	if (simul_bdt_ || simul_all_) size = bdt_boundaries[i].size() - 1;
    effratio_bs_val[i].resize(size);
    effratio_bd_val[i].resize(size);
    effratio_bs_err[i].resize(size);
    effratio_bd_err[i].resize(size);
  }
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      effratio_bs_val[i][j] = eff_bs_val[i][j] / eff_bu_val[i][j];
      effratio_bd_val[i][j] = eff_bd_val[i][j] / eff_bu_val[i][j];
      effratio_bs_err[i][j] = eff_rel_err[i][j] * effratio_bs_val[i][j];
      effratio_bd_err[i][j] = eff_rel_err[i][j] * effratio_bd_val[i][j];
      cout << "effratio_bs ("<< i << "," << j << ") = " << effratio_bs_val[i][j] << " +/- " << effratio_bs_err[i][j] << endl;
      cout << "effratio_bd ("<< i << "," << j << ") = " << effratio_bd_val[i][j] << " +/- " << effratio_bd_err[i][j] << endl;
    }
  }
}

void pdf_analysis::parse_efficiency_numbers(int offset) {

	vector < string> bdt_s(6);
	bdt_s[0] = "-cat20";
	bdt_s[1] = "-cat21";
	bdt_s[2] = "-cat20";
	bdt_s[3] = "-cat21";
	bdt_s[4] = "-cat22";
	bdt_s[5] = "-cat23";

//	bdt_s[0] = "-0.10-0.29";
//	bdt_s[1] = "-0.29-1.00";
//	bdt_s[2] = "-cat10";
//	bdt_s[3] = "-cat11";
//	bdt_s[4] = "-cat12";
//	bdt_s[5] = "-cat13";

	int channels_;
	vector <vector <double> > eff_ratio_err;
	vector <vector <double> > eff_ratio_val;
	eff_ratio_err.resize(channels);
	eff_ratio_val.resize(channels);
	for (int i = 0; i < channels; i++) {
		int size = 1;
		if (simul_bdt_ || simul_all_) size = bdt_boundaries[i].size() - 1;
		eff_ratio_err[i].resize(size);
		eff_ratio_val[i].resize(size);
	}

	int kkkk_max = 1;
	if ((simul_all_ || simul_bdt_) && final_) kkkk_max = 4;
	for (int kkkk = 0; kkkk < kkkk_max; kkkk++) {
		if (offset == 0 && (simul_all_ || simul_bdt_) && kkkk > 1) continue; /// 2011 gets 2, 2012 gets 4

		/// total efficiencies
		string filename = "anaBmm.plotResults.";
		if (offset == 0) filename += "2011";
		else filename += "2012";
		if (bdt_fit_ && final_ && offset == 0) filename += "-0.1";
		if (bdt_fit_ && final_ && offset == 2) filename += "-0.2";
		if ((simul_all_ || simul_bdt_) && final_) filename += bdt_s[kkkk + offset];
		filename += ".tex";

		string file_address = "./input/";

		if (offset == 0) {
			file_address += "2011/" + filename;
			channels_ = channels;
		}
		else if (offset == 2) {
			file_address += "2012/" + filename;
			channels_ = channels / 2;
		}
		else {
			cout << "not ready " << endl;
			exit(1);
		}
		cout << "parsing " << file_address << endl;
		FILE *file = fopen(file_address.c_str(), "r");
		if (!file) {cout << "file " << file_address << " does not exist" << endl; exit(1);}

		char buffer[2048];
		char left[1024];
		double number = 0;

		vector < pair<string, string> > end_bd(channels_);
		vector < pair<string, string> > end_bs(channels_);
		vector < pair<string, string> > end_bu(channels_);
		vector < pair<string, string> > end_Nbu(channels_);
		vector < pair<string, string> > end_bs_ratio(channels_);
		int ii = -1;
		for (int i = 0; i < channels_; i++) {
			if (simul_) ii = i;
			else ii = ch_i_;
			end_bd[i] = make_pair(Form("N-EFF-TOT-BDMM%d:val", ii), Form("N-EFF-TOT-BDMM%d:tot", ii));
			end_bs[i] = make_pair(Form("N-EFF-TOT-BSMM%d:val", ii), Form("N-EFF-TOT-BSMM%d:tot", ii));
			end_bu[i] = make_pair(Form("N-EFF-TOT-BPLUS%d:val", ii), Form("N-EFF-TOT-BPLUS%d:tot", ii));
			end_Nbu[i] = make_pair(Form("N-OBS-BPLUS%d:val", ii), Form("N-OBS-BPLUS%d:tot", ii));
			end_bs_ratio[i] = make_pair(Form("N-EFFRATIO-TOT-BSMM%d:val", ii), Form("N-EFFRATIO-TOT-BSMM%d:err", ii));
		}

		while (fgets(buffer, sizeof(buffer), file)) {
			if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
			if (buffer[0] == '\045') continue;
			if (buffer[0] == '\040') continue;
			sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);
			string left_s(left);
			for (int i = 0; i < channels_; i++) {
				size_t found;
				found = left_s.find(end_bd[i].first);
				if (found != string::npos) eff_bd_val[i+offset][kkkk] = number;
				found = left_s.find(end_bd[i].second);
				if (found != string::npos)eff_bd_err[i+offset][kkkk] = number;
				found = left_s.find(end_bs[i].first);
				if (found != string::npos) eff_bs_val[i+offset][kkkk] = number;
				found = left_s.find(end_bs[i].second);
				if (found != string::npos) eff_bs_err[i+offset][kkkk] = number;
				found = left_s.find(end_bu[i].first);
				if (found != string::npos) eff_bu_val[i+offset][kkkk] = number;
				found = left_s.find(end_bu[i].second);
				if (found != string::npos) eff_bu_err[i+offset][kkkk] = number;
				found = left_s.find(end_Nbu[i].first);
				if (found != string::npos) N_bu_val[i+offset][kkkk] = number;
				found = left_s.find(end_Nbu[i].second);
				if (found != string::npos) {
					if (number > 1000000) {
						cout << "something wrong with Bu error: " <<i+offset << " " << kkkk << " " << number << endl;
						N_bu_err[i+offset][kkkk] = N_bu_val[i+offset][kkkk]*5./100.;
					}
					else N_bu_err[i+offset][kkkk] = number;
				}
      	found = left_s.find(end_bs_ratio[i].first);
      	if (found != string::npos) eff_ratio_val[i+offset][kkkk] = number;
      	found = left_s.find(end_bs_ratio[i].second);
      	if (found != string::npos) eff_ratio_err[i+offset][kkkk] = number;
			}
    }
		fclose(file);
  }

	if (!final_) {
		/// we do not know eff and yields for the different bdt categories for the normalization channel
		for (int i = 0; i < channels_; i++) {
			for (int j = 1; j < bdt_index_max(i+offset); j++) {
				eff_bd_val[i+offset][j] = eff_bd_val[i+offset][0];
				eff_bd_err[i+offset][j] = eff_bd_err[i+offset][0];
				eff_bs_val[i+offset][j] = eff_bs_val[i+offset][0];
				eff_bs_err[i+offset][j] = eff_bs_err[i+offset][0];
				eff_bu_val[i+offset][j] = eff_bu_val[i+offset][0];
				eff_bu_err[i+offset][j] = eff_bu_err[i+offset][0];
				N_bu_val[i+offset][j] = N_bu_val[i+offset][0];
				N_bu_err[i+offset][j] = N_bu_err[i+offset][0];
				eff_ratio_val[i+offset][j] = eff_ratio_val[i+offset][0];
				eff_ratio_err[i+offset][j] = eff_ratio_err[i+offset][0];
			}
		}
	}

	for (int i = 0; i < channels_; i++) {
				for (int j = 1; j < bdt_index_max(i+offset); j++) {

				}
	}

	/// eff systs
  vector <double> eff_rel_err_sq(channels_, 0.);
  for (int i = 0; i < channels_; i++) {
  	eff_rel_err_sq[i] =   pow(acceptance_sys[i], 2)
												+ pow(mass_scale_sys[i], 2)
  											+ pow(trigger_sys[i], 2)
    										+ pow(muon_id_sys[i], 2)
    										+ pow(kaon_track_sys[i], 2);
  }

  /// BDT efficiency errors
  string filename = "anaBmm.plotReducedOverlays.";
  if (offset == 0) filename += "2011";
  else filename += "2012";
  if (bdt_fit_ && false) filename += "-0.1"; //////////////////// not yet done!
  filename += ".tex";

  string file_address = "../uml/input/";
  if (offset == 0) {
    file_address += "2011/" + filename;
  }
  else if (offset == 2) {
    file_address += "2012/" + filename;
  }
  else {
    cout << "not ready " << endl;
    exit(1);
  }
  cout << "parsing " << file_address << endl;
  FILE *file_bdt = fopen(file_address.c_str(), "r");
  if (!file_bdt) {cout << "file " << file_address << " does not exist" << endl; exit(1);}

  vector < string > NO_err(channels_);
  vector < string > CS_err(channels_);
  vector < string > MC_err(channels_);
  int ii = -1;
  for (int i = 0; i < channels_; i++) {
    if (simul_) ii = i;
    else ii = ch_i_;
    NO_err[i] = Form("relDeltaEpsNoDataNoMcchan%d:val", ii);
    CS_err[i] = Form("relDeltaEpsCsDataCsMcchan%d:val", ii);
    MC_err[i] = Form("relDeltaEpsSgMcCsMcchan%d:val", ii);
  }

  vector <double> NO_err_d(channels_);
  vector <double> CS_err_d(channels_);
  vector <double> MC_err_d(channels_);

  char buffer[2048];
  char left[1024];
  double number = 0;

  while (fgets(buffer, sizeof(buffer), file_bdt)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    if (buffer[0] == '\040') continue;
    sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);
    string left_s(left);
    for (int i = 0; i < channels_; i++) {
      size_t found;
      found = left_s.find(NO_err[i]);
      if (found != string::npos) NO_err_d[i] = number;

      found = left_s.find(CS_err[i]);
      if (found != string::npos) CS_err_d[i] = number;

      found = left_s.find(MC_err[i]);
      if (found != string::npos) MC_err_d[i] = number;
    }
  }
  for (int i = 0; i < channels_; i++) {
    eff_rel_err_sq[i] += pow(NO_err_d[i], 2)
    										+ pow(CS_err_d[i], 2)
    										/* + pow(MC_err_d[i], 2)*/
    									  ;

  }
  fclose(file_bdt);

  for (int i = 0; i < channels_; i++) {
    eff_rel_err[i+offset][0] = sqrt(eff_rel_err_sq[i]);
  }
  ///
  /// using only NO result table instead
  ///
//  for (int i = 0; i < channels_; i++) {
//  	eff_ratio_err[i+offset][0] = eff_ratio_err[i+offset][0] / (eff_ratio_val[i+offset][0] * eff_ratio_val[i+offset][0]);
//  	eff_rel_err[i+offset][0] = eff_ratio_err[i+offset][0] * eff_ratio_val[i+offset][0];
//  }

  for (int i = 0; i < channels_; i++) {
    for (int j = 1; j < bdt_index_max(i+offset); j++) {
      eff_rel_err[i+offset][j] = eff_rel_err[i+offset][0];
    }
  }

  for (int i = 0; i < channels_ && i < 2; i++) {
  	for (int j = 0; j < bdt_index_max(i + offset); j++) {
  		if (simul_) ii = i;
  		else ii = ch_i_;
  		cout << "etacat " << ii + offset;
  		if ((simul_bdt_ || simul_all_) && final_) cout << " - bdtcat " << j;
  		cout	<< ":" << endl;
  		cout << "bd eff = " << eff_bd_val[i+offset][j] << " \\pm " << eff_bd_err[i+offset][j] << endl;
  		cout << "bs eff = " << eff_bs_val[i+offset][j] << " \\pm " << eff_bs_err[i+offset][j] << endl;
  		cout << "bu eff = " << eff_bu_val[i+offset][j] << " \\pm " << eff_bu_err[i+offset][j] << endl;
  		cout << "N bu   = " << N_bu_val[i+offset][j] << " \\pm " << N_bu_err[i+offset][j] << endl;
  		cout << "eff ratio relative error " << eff_rel_err[i+offset][j] << endl;
  		cout << endl;
  	}
  }
}

void pdf_analysis::parse_external_numbers(string filename) {
  cout << "parsing " << filename << endl;
  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {cout << "file " << filename << " does not exist" << endl; exit(1);}
  char buffer[1024];
  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    sscanf(buffer, "fsfu\t%lf\t%lf", &fs_over_fu_val, &fs_over_fu_err);
    sscanf(buffer, "BuToKpsiK\t%lf\t%lf", &Bu2JpsiK_BF_val, &Bu2JpsiK_BF_err);
    sscanf(buffer, "JpsiToMuMu\t%lf\t%lf", &Jpsi2MuMu_BF_val, &Jpsi2MuMu_BF_err);
    sscanf(buffer, "Bs2MuMu_BF\t%lf\t%lf", &Bs2MuMu_SM_BF_val, &Bs2MuMu_SM_BF_err);
    sscanf(buffer, "Bd2MuMu_BF\t%lf\t%lf", &Bd2MuMu_SM_BF_val, &Bd2MuMu_SM_BF_err);

    sscanf(buffer, "Acceptance_0\t%lf", &acceptance_sys[0]);
    sscanf(buffer, "Acceptance_1\t%lf", &acceptance_sys[1]);

    sscanf(buffer, "MassScale_0\t%lf", &mass_scale_sys[0]);
    sscanf(buffer, "MassScale_1\t%lf", &mass_scale_sys[1]);

    sscanf(buffer, "KaonTrack_0\t%lf", &kaon_track_sys[0]);
    sscanf(buffer, "KaonTrack_1\t%lf", &kaon_track_sys[1]);

    sscanf(buffer, "Trigger_0\t%lf", &trigger_sys[0]);
    sscanf(buffer, "Trigger_1\t%lf", &trigger_sys[1]);

    sscanf(buffer, "MuonID_0\t%lf", &muon_id_sys[0]);
    sscanf(buffer, "MuonID_1\t%lf", &muon_id_sys[1]);

  }
  cout << "Bs2MuMu_SM_BF " <<  Bs2MuMu_SM_BF_val << " \\pm " << Bs2MuMu_SM_BF_err << endl;
  cout << "Bd2MuMu_SM_BF " <<  Bd2MuMu_SM_BF_val << " \\pm " << Bd2MuMu_SM_BF_err << endl;
  cout << "fs/fu        " <<  fs_over_fu_val << " \\pm " << fs_over_fu_err << endl;
  cout << "Bu2JpsiK_BF  " <<  Bu2JpsiK_BF_val << " \\pm " << Bu2JpsiK_BF_err << endl;
  cout << "Jpsi2MuMu_BF " <<  Jpsi2MuMu_BF_val << " \\pm " << Jpsi2MuMu_BF_err << endl;
  for (int i = 0; i < 2; i++) {
  	cout << "Acceptance " << i << " = " <<  acceptance_sys[i] << endl;
  	cout << "Mass scale " << i << " = " <<  mass_scale_sys[i] << endl;
  	cout << "Kaon track " << i << " = " <<  kaon_track_sys[i] << endl;
  	cout << "Trigger " << i << " = " <<  trigger_sys[i] << endl;
  	cout << "Muon ID " << i << " = " <<  muon_id_sys[i] << endl;
  }
}

void pdf_analysis::bdt_effs() {
	string input_file = "input/eff_bdtbins_fact.txt";
	if (bdt_fit_) input_file = "input/eff_2d_fact.txt";
  FILE *file = fopen(input_file.c_str(), "r");
  if (!file) {cout << "file " << input_file << " does not exist" << endl; exit(1);}
  char buffer[1024];
  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
        ostringstream bu_oss, bs_oss, bd_oss, peak_oss, semi_oss, comb_oss;
        bs_oss << "bdt_bs_" << i << "_" << j << "\t%lf";
        bd_oss << "bdt_bd_" << i << "_" << j << "\t%lf";
        bu_oss << "bdt_bu_" << i << "_" << j << "\t%lf";
        peak_oss << "bdt_peak_" << i << "_" << j << "\t%lf";
        semi_oss << "bdt_semi_" << i << "_" << j << "\t%lf";
        comb_oss << "bdt_comb_" << i << "_" << j << "\t%lf";
        sscanf(buffer, bs_oss.str().c_str(), &bs_bdt_factor[i][j]);
        sscanf(buffer, bd_oss.str().c_str(), &bd_bdt_factor[i][j]);
        sscanf(buffer, bu_oss.str().c_str(), &bu_bdt_factor[i][j]);
        sscanf(buffer, peak_oss.str().c_str(), &peak_bdt_factor[i][j]);
        sscanf(buffer, semi_oss.str().c_str(), &semi_bdt_factor[i][j]);
        sscanf(buffer, comb_oss.str().c_str(), &comb_bdt_factor[i][j]);
      }
    }
  }
  fclose(file);

  cout << "efficiencies and bu yields for every category:" << endl;
  int allchannels = 0;
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
    	allchannels++;
    	eff_bs_val[i][j] = bs_bdt_factor[i][j] * eff_bs_val[i][j];
    	eff_bd_val[i][j] = bd_bdt_factor[i][j] * eff_bd_val[i][j];
    	eff_bu_val[i][j] = bu_bdt_factor[i][j] * eff_bu_val[i][j];
    	N_bu_val[i][j] = bu_bdt_factor[i][j] * N_bu_val[i][j];
    	N_bu_err[i][j] = bu_bdt_factor[i][j] * N_bu_err[i][j];
      cout << "[" << i << "," << j << "]  bu eff = " <<  eff_bu_val[i][j] << " bs eff = " << eff_bs_val[i][j] << " bd eff = " << eff_bd_val[i][j] << " bu yield = " << N_bu_val[i][j] << " +/- " << N_bu_err[i][j] << endl;
      cout << "\t peak factor " << peak_bdt_factor[i][j] << " semi factor " << semi_bdt_factor[i][j] << " comb factor " << comb_bdt_factor[i][j] << endl;
    }
  }
  channels_all = allchannels;
}

void pdf_analysis::bdt_fit_effs() {
//  vector <double> eff_Bs_bdt01;
//  vector <double> eff_Bd_bdt01;
//  eff_Bs_bdt01.push_back(37547./24029.);
//  eff_Bs_bdt01.push_back(21375./11095.);
//  eff_Bs_bdt01.push_back(19498./15225.);
//  eff_Bs_bdt01.push_back(9336./6590.);
//  eff_Bd_bdt01.push_back(1960./1241.);
//  eff_Bd_bdt01.push_back(1039./563.);
//  eff_Bd_bdt01.push_back(2007./1597.);
//  eff_Bd_bdt01.push_back(889./629.);
//  for (int i = 0; i < channels; i++) {
//  	eff_bs_val[i][0] *= eff_Bs_bdt01[i];
//  	eff_bs_err[i][0] *= eff_Bs_bdt01[i];
//  	eff_bd_val[i][0] *= eff_Bd_bdt01[i];
//  	eff_bd_err[i][0] *= eff_Bd_bdt01[i];
//  }
	cout << "not used" << endl;
	exit(1);
}

void pdf_analysis::setSBslope(string pdf, RooAbsData *sb_data) {

  RooAbsData* subdata;
  if (!simul_bdt_ && !simul_all_) subdata = sb_data->reduce(Form("etacat==etacat::etacat_%d", channel));
  else subdata = sb_data->reduce(Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", channel, channel_bdt));

  RooExponential expo_temp("expo_temp", "expo_temp", *ws_->var("Mass"), *ws_->var(name("exp_comb", channel, channel_bdt)));
  cout << "comb fitting " << subdata->numEntries() << endl;
  subdata->Print();
  cout << "with " << endl;
  expo_temp.Print();
  ws_->var(name("exp_comb", channel, channel_bdt))->Print();
  cout << ws_->var(name("exp_comb", channel, channel_bdt))->getVal() << endl;
  expo_temp.fitTo(*subdata, Range("sb_lo,sb_hi"));
//  expo_temp.fitTo(*subdata);
  if (print_ || 1) {
    TCanvas c("c", "c", 600, 600);
    RooPlot* rp = ws_->var("Mass")->frame();
    subdata->plotOn(rp, Binning(40));
    expo_temp.plotOn(rp, Range("sb_lo,sb_hi"));
    expo_temp.paramOn(rp);
    rp->Draw();
    c.Print((get_address(pdf, "", true) + ".gif").c_str());
    c.Print((get_address(pdf, "", true) + ".pdf").c_str());
    delete rp;
  }
}

string pdf_analysis::get_title(int i) {
	string output(" ");
	if (i == 0) output = "CMS - L = 5 fb^{-1} #sqrt{s} = 7 TeV - Barrel 2011";
	if (i == 1) output = "CMS - L = 5 fb^{-1} #sqrt{s} = 7 TeV - Endcap 2011";
	if (i == 2) output = "CMS - L = 20 fb^{-1} #sqrt{s} = 8 TeV - Barrel 2012";
	if (i == 3) output = "CMS - L = 20 fb^{-1} #sqrt{s} = 8 TeV - Endcap 2012";
	return output;
}

int pdf_analysis::bdt_index(int eta_ch, double bdt) {
	if (bdt < bdt_boundaries[eta_ch][0]) return -1;
	for (int index = 0; index < (int)bdt_boundaries[eta_ch].size(); index++){
		if (bdt > bdt_boundaries[eta_ch][index+1]) continue;
		else return index;
	}
	return -1;
}

int pdf_analysis::bdt_index_max(int eta_ch) {
  if (!simul_bdt_ && !simul_all_) return 1;
  return bdt_boundaries[eta_ch].size() - 1;
}

int pdf_analysis::super_index(int eta_ch, int bdt_ch) {
	int superindex = -1;
	for (int i = 0; i < (int)bdt_boundaries.size(); i++) {
		for (int j = 0; j < (int)bdt_boundaries[i].size() - 1; j++) {
			superindex++;
			if (eta_ch == i && bdt_ch == j) return superindex;
		}
	}
	return superindex;
}

vector <int> pdf_analysis::get_EtaBdt_bins(int index) {
  vector <int> indexes(2, -1);

  int superindex = -1;

	for (int i = 0; i < (int)bdt_boundaries.size(); i++) {
		indexes[0]++;
		for (int j = 0; j < (int)bdt_boundaries[i].size() - 1; j++) {
			indexes[1]++;
			superindex++;
			if (index == superindex) return indexes;
		}
		indexes[1] = -1;
	}
	return indexes;
}

void pdf_analysis::fill_bdt_boundaries() {
	string input_file("input/bdtbins.txt");
	if (bdt_fit_) input_file = "input/2dbins.txt";
	FILE *file = fopen(input_file.c_str(), "r");
	char buffer[1024];
	int in = 0;
	while (fgets(buffer, sizeof(buffer), file)) {
		cout << "eta category " << in << "; binning = " << buffer;
		string binning = buffer;
		bool not_ended = true;
		while (not_ended) {
			size_t found;
			found = binning.find_first_of(",");
			string sub = binning.substr(0, found);
			binning.erase(0, found+1);
			ostringstream number_oss;
			number_oss << sub;
			double n = atof(number_oss.str().c_str());
			bdt_boundaries[in].push_back(n);
			if (n > 0.999 || in >= channels) not_ended = false;
		}
		in++;
		if (in >= channels) break;
	}
	fclose(file);
}

void pdf_analysis::get_bkg_yields(string filename, string dir, int offset, int bdtcats) {

  string peakdecays[] = {"BgPeakLo", "BgPeakBd", "BgPeakBs", "BgPeakHi"};
//  string semidecays[] = {"BgRslsLo", "BgRslsBd", "BgRslsBs", "BgRslsHi"};
  string semidecays[] = {"BgRslLo", "BgRslBd", "BgRslBs", "BgRslHi"};
  string combdecays[] = {"BgCombLo", "BgCombBd", "BgCombBs", "BgCombHi"};

  Double_t peak_exp_t[4][2];
  Double_t semi_exp_t[4][2];
  Double_t comb_exp_t[4][2];
  Double_t peak_syst_err_t[4][2];
  Double_t semi_syst_err_t[4][2];
  Double_t comb_syst_err_t[4][2];
  Double_t peak_stat_err_t[4][2];
  Double_t semi_stat_err_t[4][2];
  Double_t comb_stat_err_t[4][2];

  string full_address = dir + filename;
  FILE *file = fopen(full_address.c_str(), "r");
  if (!file) {cout << "file " << full_address << " does not exist"; exit(1);}

  char buffer[1024];
  char left[1024];
  double number;
  int peak_n = sizeof(peakdecays)/sizeof(string);
  int semi_n = sizeof(semidecays)/sizeof(string);
  int comb_n = sizeof(combdecays)/sizeof(string);
  vector <double> peak_exp(2, 0);
  vector <double> semi_exp(2, 0);
  vector <double> comb_exp(2, 0);
  vector <double> peak_syst_err(2, 0);
  vector <double> semi_syst_err(2, 0);
  vector <double> comb_syst_err(2, 0);
  vector <double> peak_stat_err(2, 0);
  vector <double> semi_stat_err(2, 0);
  vector <double> comb_stat_err(2, 0);
  string end_0("0:val}");
  string end_1("1:val}");
  string err_stat_0("0:e1}");
  string err_stat_1("1:e1}");
  string err_syst_0("0:e2}");
  string err_syst_1("1:e2}");

  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);
    string left_s(left);
    for (int i = 0; i < peak_n; i++) { /// peak

      size_t found = left_s.find(peakdecays[i]);
      if (found != string::npos) {
        found = left_s.find(end_0);
        if (found != string::npos) {
        	peak_exp_t[i][0] = number;
        }
        found = left_s.find(end_1);
        if (found != string::npos) {
          peak_exp_t[i][1] = number;
        }
        found = left_s.find(err_syst_0);
        if (found != string::npos) {
          peak_syst_err_t[i][0] = number;
        }
        found = left_s.find(err_syst_1);
        if (found != string::npos) {
          peak_syst_err_t[i][1] = number;
        }
        found = left_s.find(err_stat_0);
        if (found != string::npos) {
        	peak_stat_err_t[i][0] = number;
        }
        found = left_s.find(err_stat_1);
        if (found != string::npos) {
        	peak_stat_err_t[i][1] = number;
        }
      }
    }
    for (int i = 0; i < semi_n; i++) { /// semi
      size_t found = left_s.find(semidecays[i]);
      if (found != string::npos) {
        found = left_s.find(end_0);
        if (found != string::npos) {
          semi_exp_t[i][0] = number;
        }
        found = left_s.find(end_1);
        if (found != string::npos) {
          semi_exp_t[i][1] = number;
        }
        found = left_s.find(err_syst_0);
        if (found != string::npos) {
        	semi_syst_err_t[i][0] = number;
        }
        found = left_s.find(err_syst_1);
        if (found != string::npos) {
        	semi_syst_err_t[i][1] = number;
        }
        found = left_s.find(err_stat_0);
        if (found != string::npos) {
        	semi_stat_err_t[i][0] = number;
        }
        found = left_s.find(err_stat_1);
        if (found != string::npos) {
        	semi_stat_err_t[i][1] = number;
        }
      }
    }
    for (int i = 0; i < comb_n; i++) { /// semi
    	size_t found = left_s.find(combdecays[i]);
    	if (found != string::npos) {
    		found = left_s.find(end_0);
    		if (found != string::npos) {
    			comb_exp_t[i][0] = number;
    		}
    		found = left_s.find(end_1);
    		if (found != string::npos) {
    			comb_exp_t[i][1] = number;
    		}
    		found = left_s.find(err_syst_0);
    		if (found != string::npos) {
    			comb_syst_err_t[i][0] = number;
    		}
    		found = left_s.find(err_syst_1);
    		if (found != string::npos) {
    			comb_syst_err_t[i][1] = number;
    		}
    		found = left_s.find(err_stat_0);
    		if (found != string::npos) {
    			comb_stat_err_t[i][0] = number;
    		}
    		found = left_s.find(err_stat_1);
    		if (found != string::npos) {
    			comb_stat_err_t[i][1] = number;
    		}
    	}
    }
  }
  fclose(file);

  for (int i = 0; i < 4; i++) {
  	for (int j = 0; j < 2; j++) {
  		peak_exp[j] += peak_exp_t[i][j];
  		peak_syst_err[j] += peak_syst_err_t[i][j];
  		peak_stat_err[j] += peak_stat_err_t[i][j];
  		semi_exp[j] += semi_exp_t[i][j];
  		semi_syst_err[j] += semi_syst_err_t[i][j];
  		semi_stat_err[j] += semi_stat_err_t[i][j];
  		comb_exp[j] += comb_exp_t[i][j];
  		comb_syst_err[j] += comb_syst_err_t[i][j];
  		comb_stat_err[j] += comb_stat_err_t[i][j];
  	}
  }

  vector <double> peak_err(2, 0);
  vector <double> semi_err(2, 0);
  vector <double> comb_err(2, 0);
  for (int i = 0; i < 2; i++) {
  	peak_err[i] = sqrt(pow(peak_stat_err[i], 2) + pow(peak_syst_err[i], 2));
  	semi_err[i] = sqrt(pow(semi_stat_err[i], 2) + pow(semi_syst_err[i], 2));
  	comb_err[i] = sqrt(pow(comb_stat_err[i], 2) + pow(comb_syst_err[i], 2));
  }

  string outdir("./input/2011/");
  if (offset == 2) outdir = "./input/2012/";
  string tail("/bkg_yields.txt");
  if (final_) {
  	if (bdt_fit_) tail = "/bkg_yields_2D.txt";
  	if (simul_bdt_ || simul_all_) tail = Form("/bkg_yields_bdt_%d.txt", bdtcats);
  }


  string full_output = outdir + tail;
  FILE* file_out = fopen(full_output.c_str(), "w");
  if (final_ && (simul_bdt_ || simul_all_)) {
  	int bdt_index = -1;
  	if (bdtcats < 2) bdt_index = bdtcats;
  	else bdt_index = bdtcats - 2;
  	for (int i = 0; i < 2; i++) {
  		fprintf(file_out, "N_peak_%d_%d\t%lf\t%lf\n", i+offset, bdt_index, peak_exp[i], peak_err[i]);
  		fprintf(file_out, "N_semi_%d_%d\t%lf\t%lf\n", i+offset, bdt_index, semi_exp[i], semi_err[i]);
//  		fprintf(file_out, "N_comb_%d_%d\t%lf\t%lf\n", i+offset, bdt_index, comb_exp[i] == 0 ? 1 : comb_exp[i], comb_err[i]); /// HACK
  		fprintf(file_out, "N_comb_%d_%d\t%lf\t%lf\n", i+offset, bdt_index, comb_exp[i], comb_err[i]);
  	}
  }
  else {
  	for (int i = 0; i < 2; i++) {
  		fprintf(file_out, "N_peak_%d\t%lf\t%lf\n", i+offset, peak_exp[i], peak_err[i]);
  		fprintf(file_out, "N_semi_%d\t%lf\t%lf\n", i+offset, semi_exp[i], semi_err[i]);
//  		fprintf(file_out, "N_comb_%d\t%lf\t%lf\n", i+offset, comb_exp[i] == 0 ? 1 : comb_exp[i], comb_err[i]); /// HACK
  		fprintf(file_out, "N_comb_%d\t%lf\t%lf\n", i+offset, comb_exp[i], comb_err[i]);
  	}
  }
  fprintf(file_out, "######\n");
  fclose(file_out);
  cout << "in file : " << full_output << endl;
  system(Form("cat %s", full_output.c_str()));
}

void pdf_analysis::get_bkg_from_tex() {

	vector <string> input_dir(2);
	input_dir[0] = "./input/2011/";
	input_dir[1] = "./input/2012/";

	vector <string> input_file(2);
	input_file[0] = "anaBmm.plotResults.2011.tex";
	input_file[1] = "anaBmm.plotResults.2012.tex";
	if (bdt_fit_ && final_) {
		input_file[0] = "anaBmm.plotResults.2011-0.1.tex";
		input_file[1] = "anaBmm.plotResults.2012-0.2.tex";
	}
	if ((simul_all_ || simul_bdt_) && final_) {
		input_file[0] = "anaBmm.plotResults.2011-cat20.tex";
		input_file[1] = "anaBmm.plotResults.2011-cat21.tex";
		input_file.push_back("anaBmm.plotResults.2012-cat20.tex");
		input_file.push_back("anaBmm.plotResults.2012-cat21.tex");
		input_file.push_back("anaBmm.plotResults.2012-cat22.tex");
		input_file.push_back("anaBmm.plotResults.2012-cat23.tex");
//		input_file[0] = "anaBmm.plotResults.2011-0.10-0.29.tex";
//		input_file[1] = "anaBmm.plotResults.2011-0.29-1.00.tex";
//		input_file.push_back("anaBmm.plotResults.2012-cat10.tex");
//		input_file.push_back("anaBmm.plotResults.2012-cat11.tex");
//		input_file.push_back("anaBmm.plotResults.2012-cat12.tex");
//		input_file.push_back("anaBmm.plotResults.2012-cat13.tex");
	}

	for (int i = 0; i < input_file.size(); i++) {
		int offset = 2*i;
		int bdtcats = -1;
		int i_dir = i;
		if ((simul_all_ || simul_bdt_) && final_) {
			bdtcats = i;
			if (i < 2) {
				offset = 0;
				i_dir = 0;
			}
			else {
				offset = 2;
				i_dir = 1;
			}

		}
		get_bkg_yields(input_file[i], input_dir[i_dir], offset, bdtcats);
	}

	string command("rm input/bkg_yields.txt; cat input/2011/bkg_yields.txt >> input/bkg_yields.txt; cat input/2012/bkg_yields.txt >> input/bkg_yields.txt;");
	if (final_) {
		if (bdt_fit_) command = "rm input/bkg_yields_2D.txt; cat input/2011/bkg_yields_2D.txt >> input/bkg_yields_2D.txt; cat input/2012/bkg_yields_2D.txt >> input/bkg_yields_2D.txt;";
		if (simul_all_ || simul_bdt_) {
			command = "rm input/bkg_yields_bdt.txt;";
			for (int i = 0; i < input_file.size(); i++) {
				if (i < 2) command += Form("cat input/2011/bkg_yields_bdt_%d.txt >> input/bkg_yields_bdt.txt;", i);
				else command += Form("cat input/2012/bkg_yields_bdt_%d.txt >> input/bkg_yields_bdt.txt;", i);
			}
		}
	}

	system(command.c_str());
	string all_bgk("input/bkg_yields.txt");
	if (final_) {
		if (bdt_fit_) all_bgk = "input/bkg_yields_2D.txt";
		if (simul_all_ || simul_bdt_) all_bgk = "input/bkg_yields_bdt.txt";
	}
  set_bkg_normalization(all_bgk);
}

void pdf_analysis::set_bdt_min(vector<double> cuts) {
	cout << "setting minimum bdt = ";
	bdt_cuts = cuts;
	for (int i = 0; i < channels; i++) {
		ws_->var(name("bdt", i))->setMin(bdt_cuts[i]);
//		ws_->var(name("bdt", i))->setMin(-1);
	}
  cout << ws_->var(name("bdt", 0))->getMin() << endl;
}

void pdf_analysis::hack_comb_slope() {
	cout << "hackering exp_comb" << endl;
	if (!bdt_fit_ && !simul_all_ && !simul_bdt_) {
		ws_->var(name("exp_comb", 0))->setVal(-0.626);
		ws_->var(name("exp_comb", 1))->setVal(0.54);
		ws_->var(name("exp_comb", 2))->setVal(0.04);
		ws_->var(name("exp_comb", 3))->setVal(-0.482);
	}
	if (bdt_fit_ && !simul_all_ && !simul_bdt_) {
		ws_->var(name("exp_comb", 0))->setVal(-6.26);
		ws_->var(name("exp_comb", 1))->setVal(0.2);
		ws_->var(name("exp_comb", 2))->setVal(-1.827);
		ws_->var(name("exp_comb", 3))->setVal(-7.24);
	}
}


