#ifndef PDF_FITDATA_H
#define PDF_FITDATA_H

#include <sstream>
#include <iostream>
#include <iomanip>

#include "pdf_analysis.h"

#include "TTree.h"
#include "TH2D.h"

#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooRandom.h"
#include "RooGamma.h"
#include "RooMinuit.h"
#include "RooLognormal.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TPolyMarker.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/RooStatsUtils.h"

using namespace RooStats;

class pdf_fitData : public pdf_analysis {
  public:
    pdf_fitData(bool print, string input_estimates = "", string range = "all", int BF = 0, bool SM = false, bool bd_constr = false, int simul = 1, int simulbdt = 1, int simulall = 1, bool pee_ = false , bool bdt_fit = false, bool final = false, string ch_s = "0", int sig = -1, bool asimov = false, bool syste = false, bool randomsyste = false, bool rare_constr = false, int nexp = 10, bool bd = false, string years = "0");
    ~pdf_fitData();

    void set_estimate();
    void parse_estimate(string input = "");
    void print_estimate();

    void print();
    void print_each_channel(string var = "Mass", string output = "", RooWorkspace *ws = 0, RooDataSet *rds_ = 0);

    void define_dataset();
    void make_dataset(bool cut_b, vector<double> cut_, double bdt_cut, TTree *tree, int offset = 0);
    void make_pdf_input(string root_s);
    void set_starting_N();
    void hack_ws(string frozen_ws_address);

    RooDataSet* global_data;
    RooDataHist* global_datahist;
    RooSimultaneous* simul_pdf;

    RooFitResult* fit_pdf(bool do_not_import = false, string pdf_name = "", int strategy = 2, bool extconstr = true);
    void significance();
    void save();

    double lumi;
    bool random;
    void setnewlumi();
    void set_syst();

    int proof;

    void extract_N_inRanges();
    void profile_NLL();
    void doublescan();

    bool SMIsNull;

    void set_final_pdf();
    void reset_minmax();

    void print_gaussian_constraints();
    void tweak_pdf(int free = 0);

    bool berns_;
    bool freeze;
    bool stat_only_;
    bool doubleNull;

    void save_for_cls();

  protected:

    void define_constraints(int i, int j);
    void define_simul();
    void define_perchannel_pdf();
    void randomize_constraints(RooWorkspace* ws);
    void define_total_extended(int i, int j, bool h = false); // final pdf with all extended components

    string input_estimates_;
    vector <double> estimate_bs;
    vector <double> estimate_bd;
    vector <double> estimate_peak;
    vector <double> estimate_semi;
    vector <double> estimate_comb;

    vector <vector <double> > estimate2D_bs;
    vector <vector <double> > estimate2D_bd;
    vector <vector <double> > estimate2D_peak;
    vector <vector <double> > estimate2D_semi;
    vector <vector <double> > estimate2D_comb;
    vector <vector <double> > estimate2D_channel;

    bool parse(char *cutName, float cut);
    string input_cuts_;
    bool asimov_;
    int sign;

    string pdfname;
    int NExp;
    bool Bd;
    string years_;

    void freeze_norm (bool set);
    void stat_error(bool stat_only);

    RooHistPdf * pdf_to_hist(RooAbsPdf* kpdf);
    void merge_mass_to_hist();

  private:

    TFile* ws_file_input;
    RooWorkspace* ws_input;

    void fill_dataset(RooDataSet* dataset, bool cut_b, vector<double> cut_, double bdt_cut, TTree *tree, int offset);
    void changeName(RooWorkspace *ws, int str);

    Double_t sig_hand();
    void sig_plhc();
    void sig_plhts();
    void sig_hybrid_plhts();
    void sig_hybrid_roplhts();
    void make_prior();
    void make_models();

    void plot_hypotest(HypoTestResult * hts);

};

#endif // PDF_FITDATA_H
