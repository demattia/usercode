#ifndef PDF_ANALYSIS_H
#define PDF_ANALYSIS_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TError.h"
#include "RooGlobalFunc.h"

#include "Math/WrappedMultiTF1.h"
#include "Fit/DataOptions.h"
#include "Fit/DataRange.h"
#include "Fit/BinData.h"
#include "HFitInterface.h"
#include "Fit/Chi2FCN.h"
#include "Fit/Fitter.h"

#include "RooWorkspace.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooUniform.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooSimWSTool.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooExponential.h"
#include "RooSimultaneous.h"
#include "RooProduct.h"
#include "RooGaussModel.h"
#include "RooFFTConvPdf.h"
#include "RooSimWSTool.h"
#include "RooSuperCategory.h"
#include "RooMultiCategory.h"
#include "RooIntegralMorph.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooMsgService.h"
#include "RooBernstein.h"

using namespace std;
using namespace RooFit;

class pdf_analysis {
public:
  pdf_analysis(bool print, string ch_s = "0", string range = "all", int BF = 0, bool SM = false, bool bd_constr = false, int simul = 1, int simulbdt = 1, int simulall = 1, bool pee_ = false, bool bdt_fit = false, bool final = false);
  ~pdf_analysis();
  void set_ws(RooWorkspace *ws) {ws_ = ws;}
  RooWorkspace* get_ws() {return ws_;}
  
  void initialize();
  TH1D* define_massRes_pdf(RooDataSet *rds, string name, bool rkeys);
  TH1D* define_semi_pdf(RooDataSet *rds, string name, bool rkeys);
  TH1D* define_bdt_pdf(RooDataSet *rds, string name, TFile* bdt_syst_f = 0, bool rkeys = false, Double_t bdt_min = -1.);

  void define_N();
  void define_pdfs(bool semi);

  void set_SMratio(double ratio) {ratio_ = ratio;}

  double getformulaerror(RooFormulaVar* formula, int vars);
  
  string define_pdf_sum(string name);
  void define_total_fractional(int i, int j); // final pdf with fractional components, and also extended

  void getBFnumbers(string numbers_filename);

  void fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error = true, bool hesse = true, bool setconstant = true);
  void set_pdf_constant(string pdf);


  void get_bkg_from_tex();
  void get_bkg_yields(string filename, string dir, int offset = 0, int bdtcats = -1);

  string pdf_name;

  int channel;
  int channel_bdt;
  bool SM_;
  bool bd_constr_;
  bool simul_;
  bool simul_bdt_;
  bool simul_all_;

  int BF_;
  RooRealVar* Mass;
  RooRealVar* MassRes;
  RooRealVar* ReducedMassRes;
  RooRealVar* bdt;
  RooRealVar* bdt_0;
  RooRealVar* bdt_1;
  RooRealVar* bdt_2;
  RooRealVar* bdt_3;
  RooRealVar* eta;
  RooRealVar* m1eta;
  RooRealVar* m2eta;
  RooRealVar* weight;

  RooCategory* channels_cat;
  RooCategory* bdt_cat;
  RooCategory* all_cat;
  RooSuperCategory* super_cat;

  int channels;
  int channels_bdt;
  int channels_all;
  string range_;
  RooFitResult* RFR;

  int verbosity;

  void simsplit();

  double getErrorHigh(RooRealVar* var);
  double getErrorLow(RooRealVar* var);

  bool old_tree;
  bool pee;
  bool bdt_fit_;
  bool no_legend;

  void gen_and_fit(string pdfname);
  bool print_;

  const char* name(string name, int i, int j = -1);

  bool newcomb_;
  void setSBslope(string pdf, RooAbsData* sb_data);
  bool rare_constr_;

  int bdt_index(int eta_ch, double bdt);
  int bdt_index_max(int eta_ch) ;
  int super_index(int eta_ch, int bdt_ch);
  int supercatdim;

  void set_bdt_min(vector <double> cuts);

  bool bdtsplit;
  void hack_comb_slope();
  bool minos;

protected:
  string default_console_color;
  string red_color_bold;
  string green_color_bold;
  string meth_;
  string ch_s_;
  int ch_i_;
  RooWorkspace* ws_;

  double ratio_;

  string get_address(string name, string pdf = "", bool channeling = true);
  string get_title(int i);

  RooArgSet* obs;

  vector < string > source;

  void parse_external_numbers(string filename);
  void parse_efficiency_numbers(int offset = 0);

  double Bs2MuMu_SM_BF_val;
  double Bs2MuMu_SM_BF_err;
  double Bd2MuMu_SM_BF_val;
  double Bd2MuMu_SM_BF_err;
  double fs_over_fu_val;
  double Jpsi2MuMu_BF_val;
  double Bu2JpsiK_BF_val;
  double fs_over_fu_err;
  double Jpsi2MuMu_BF_err;
  double Bu2JpsiK_BF_err;

  double one_over_BRBR_val;
  double one_over_BRBR_err;

  vector <vector <double> > effratio_bs_val;
  vector <vector <double> > effratio_bd_val;
  vector <vector <double> > eff_rel_err;
  vector <vector <double> > effratio_bs_err;
  vector <vector <double> > effratio_bd_err;

  vector <vector <double> > eff_bd_val;
  vector <vector <double> > eff_bs_val;
  vector <vector <double> > eff_bu_val;
  vector <vector <double> > N_bu_val;
  vector <vector <double> > eff_bd_err;
  vector <vector <double> > eff_bs_err;
  vector <vector <double> > eff_bu_err;
  vector <vector <double> > N_bu_err;

  vector <vector <double> > bs_bdt_factor;
  vector <vector <double> > bd_bdt_factor;
  vector <vector <double> > bu_bdt_factor;
  vector <vector <double> > peak_bdt_factor;
  vector <vector <double> > semi_bdt_factor;
  vector <vector <double> > comb_bdt_factor;

  vector <double> acceptance_sys;
  vector <double> mass_scale_sys;
  vector <double> kaon_track_sys;
  vector <double> trigger_sys;
  vector <double> muon_id_sys;

  bool syst;
  bool randomsyst;

  vector <string> massrange_names;
  vector <int> get_EtaBdt_bins(int index);
  vector < vector <double> > bdt_boundaries;
  void fill_bdt_boundaries();

  vector <double> bdt_cuts;
  bool final_;
  void set_bkg_normalization(string input); //set peak fraction parameter to Bu2JpsiK

  void define_comb2(int i, int j);


private:
  void define_bs(int i, int j);
  void define_bd(int i, int j);
  void define_peak(int i, int j);
  void define_semi(int i, int j, bool semi);
  void define_comb(int i, int j);
  void define_signals(int i, int j);
  void define_rare(int i, int j);
  void define_rare3(int i, int j);
  void define_bf(int i, int j);
  void print_pdf(RooAbsPdf* pdf, RooRealVar * var);
  void bdt_effs();
  void bdt_fit_effs();



  void print(RooAbsData *data, string output = "");
};

#endif // PDF_ANALYSIS_H
