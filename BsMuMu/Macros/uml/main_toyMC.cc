/* 
 * File:   main_toyMC.cc
 * Author: lucamartini
 *
 * Created on 24 aprile 2012, 18.25
 */
#include "pdf_analysis.h"
#include "pdf_toyMC.h"

#include <stdio.h>
#include <stdlib.h>

#include "CommonFun.h"

#include "TFile.h"

#include "RooWorkspace.h"

/*
 * 
 */

int main(int argc, char** argv) {

  parse_options(argc, argv);
  if (!simul) {
    if (!input || !pdf) help();
    parse_input(input_ws);
  }
//  if (!(!bias_s.compare("no") || !bias_s.compare("c+") || !bias_s.compare("c-") || !bias_s.compare("tau+") || !bias_s.compare("tau-") || !bias_s.compare("sig+") || !bias_s.compare("sig-"))) { cout << "I don't understand what to bias: please enter -bias c+, c-, tau+, tau-, sig+, sig-" << endl; exit(EXIT_FAILURE);}
  if (!pdf_test_b) pdf_test = pdf_toy;
  pdf_toyMC toy1(print, input_estimates, "all", BF, SM, bd_const, inputs, (!simul_all) ? inputs_bdt : 1, inputs_all, pee, bdt_fit, final, ch_s, sig_meth, asimov, syst, randomsyst, rare_constr, NExp, Bd, years_opt, bias_s);
  toy1.random = true;
  TFile* input_f = new TFile(input_ws.c_str());
  RooWorkspace* ws = (RooWorkspace*)input_f->Get("ws");
  vector <double> cuts_v(inputs, -1);
  if (cuts_f_b) cuts_v = cut_bdt_file();
  if (bdt_fit) toy1.set_bdt_min(cuts_v);
  toy1.set_ws(ws);

  toy1.get_bkg_from_tex();
  toy1.set_estimate();
  toy1.parse_estimate();
  toy1.set_starting_N();
  toy1.setnewlumi();
  if (berns) toy1.berns_ = true;
  toy1.set_final_pdf();

  if (simul_bdt || simul_all) toy1.parse_estimate("input/estimates_5_20_bdt_semi.txt");
  else toy1.parse_estimate("input/estimates_5_20_semi.txt");
  toy1.set_starting_N();

  if (hack_semi2011) toy1.hack_ws("output/frozen/ws_simul4_bdt_BF2_PEE.root");
  toy1.set_syst();
  toy1.tweak_pdf(free_rare);
  toy1.print_estimate();
  if (!simul_all && !simul_bdt && hack) toy1.hack_comb_slope();
  if (roomcs) toy1.mcstudy(pdf_toy, pdf_test);
  if (minos) toy1.minos = true;
  if (!roomcs) toy1.generate(pdf_toy, pdf_test);
  cout << "finishing time ";
  system("date");
  return EXIT_SUCCESS;
}

