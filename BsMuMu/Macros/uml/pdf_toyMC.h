/* 
 * File:   pdf_toyMC.h
 * Author: lucamartini
 *
 * Created on 26 aprile 2012, 11.03
 */

#ifndef PDF_TOYMC_H
#define	PDF_TOYMC_H

#include "pdf_fitData.h"

#include "TPaveStats.h"
#include "TMath.h"

#include "RooArgSet.h"
#include "RooMCStudy.h"
#include "RooDLLSignificanceMCSModule.h"

class pdf_toyMC : public pdf_fitData {
public:
  pdf_toyMC(bool print, string input_estimates = "", string range = "all", int BF = 0, bool SM = false, bool bd_constr = false, int simul = 1, int simulbdt = 1, int simulall = 1, bool pee_ = false, bool bdt_fit = false, bool final = false, string ch_s = "0", int sig = -1, bool asimov = false, bool syste = false, bool randomsyste = false, bool rare_constr = false, int nexp = 10, bool bd = false, string years = "0", string bias = "no");
  ~pdf_toyMC();

  void generate(string pdf_toy, string pdf_test = "total");
  void mcstudy(string pdf_toy, string pdf_test = "total");

private:
  string bias_;
  
  vector <vector <RooDataSet*> > residual_rds_bs;
  vector <vector <RooDataSet*> > residual_rds_bd;
//  vector <vector <RooDataSet*> > residual_rds_peak;
//  vector <vector <RooDataSet*> > residual_rds_semi;
  vector <vector <RooDataSet*> > pull_rds_bs;
  vector <vector <RooDataSet*> > pull_rds_bd;
  vector <vector <RooDataSet*> > pull_rds_peak;
  vector <vector <RooDataSet*> > pull_rds_semi;
  vector <vector <RooDataSet*> > pull_rds_comb;

  vector <vector <RooRealVar*> > residual_bs;
  vector <vector <RooRealVar*> > residual_bd;
//  vector <vector <RooRealVar*> > residual_peak;
//  vector <vector <RooRealVar*> > residual_semi;
  vector <vector <RooRealVar*> > pull_bs;
  vector <vector <RooRealVar*> > pull_bd;
  vector <vector <RooRealVar*> > pull_peak;
  vector <vector <RooRealVar*> > pull_semi;
  vector <vector <RooRealVar*> > pull_comb;
  RooRealVar* pull_BF_bs;
  RooDataSet* pull_rds_BF_bs;
  RooRealVar* pull_BF_bd;
  RooDataSet* pull_rds_BF_bd;
  RooRealVar * pull_BF_bs_minos;
  RooDataSet* pull_rds_BF_bs_minos;

  void fit_pulls(RooRealVar *pull, RooDataSet *rds, int i, int j);
  void print_histos(TH1D* histos, int i, int j);
  RooFitResult* fit_pdf (RooAbsData* data, int printlevel = -1, RooWorkspace *ws = 0);
  Double_t sig_hand(RooAbsData *data, int printlevel, RooWorkspace *ws);
  void do_bias(RooWorkspace* ws);
  Double_t median(TH1D * h);

};

#endif	/* PDF_TOYMC_H */

