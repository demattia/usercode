#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <cmath>

#include "TH1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"

using namespace std;

/// options
static string input_name;
static string input_estimates;
static string input_ws;
static string meth = "bdt";
static string ch_s = "-1";
static string ch_bdt_s = "-1";
static string pdf_toy = "total";
static string pdf_test;
static string tree_name = "SgData_bdt";
static string bias_s = "no";
static string cuts_f = "no";
static string rare_f = "no";
static bool print = false;
static bool simul = false;
static bool simul_bdt = false;
static int BF = 0;
static int NExp = 1;
static int ch_i = -1;
static int ch_bdt_i = -1;
static int inputs = 1;
static int inputs_bdt = 1;
static int inputs_all = 1;
static int sig_meth = -1;
static int proof = 1;
bool syst = false;
bool randomsyst = false;
bool shapesyst = false;
//static string cuts = "bdt>-10.";
static string years_opt = "0";
bool input = false, output = false, channel = false, estimate = false, pdf = false, roomcs = false, SM = false, bd_const = false, pdf_test_b = false, bias = false, SB = false, pee = false, no_legend = false, bdt_fit = false, cuts_b = false, cuts_f_b = false, channel_bdt = false, asimov = false, rare_constr = false, rkeys = false;
static bool newcomb = false;
static bool method = true;
static bool Bd = false;
static bool SMIsNull = false;
static bool LLprofile = false;
static bool hack_semi2011 = false;
bool simul_all = false;
int free_rare = 3;
bool make_bdt_binning_inputs = false;
static string bdtbinnings_s = "input/bdtbinnings.root";
double bdt_cut = -1;
bool semi = false;
bool bdtsplit = false;
bool final = false;
bool hack = false;
bool minos = false;
bool berns = false;
bool freeze = false;
bool stat_only = false;
bool null = false;
bool LLdouble = false;

void help() {

  cout << endl;
  cout << ">>>>>>>>> main_make_pdf.o: makes pdf workspace" << endl;
  cout << "choose one between:" << endl;
  cout << "\t -cha {0, 1} \t barrel OR endcap input" << endl;
  cout << "\t -simul # \t simultaneous fit of # eta channels" << endl;
  cout << "-bdt_fit \t bdt_fit" << endl;
  cout << "-simul_bdt # \t simultaneous fit of # bdt channels" << endl;
  cout << "-simul_all # \t simultaneous fit of # eta and bdt channels" << endl;
  cout << "-BF {1,2,3} \t imposing the same BF in each channel for bs only (1) or for bs and bd (2) or for bs/bd and bd (3)" << endl;
  cout << "-SM \t with SM constraints (incompatible with -bd_const)" << endl;
  cout << "-bd_const \t with Bd constrainted to Bs, over all different channels (incompatible with -SM)" << endl;
  cout << "-print \t save the fits to gif and pdf if -no_legend without parameters on canvas" << endl;
//  cout << "-cuts #cutstring \t string cut for small tree: i.e. \"bdt>0.&&bdt<0.18\"; default is \"" << cuts << "\"" << endl;
  cout << "-cuts_file \t file containing bdt cuts for small tree" << endl;
  cout << "-pee \t per-event-error" << endl;
  cout << "-y {0,1,all} \t year 2011, 2012 or both (this last works only with simul)" << endl;
  cout << "-newcomb \t exponential combinatorial bkg study" << endl;
  cout << "-rkeys \t RooKeysPdf for MassRes and BDT" << endl;
  cout << "-bdtbins \t produces only inputs for evaluation of bdt binnings" << endl;
  cout << "-bdt_cut # \t minimum bdt" << endl;
  cout << "-semi \t semi mass shape from mc data" << endl;
  cout << "-bdtsplit \t splits bdt in 1,2,3,4" << endl;
  cout << "-final \t don't change input numbers, take the ones from AN" << endl;
  cout << endl;
  cout << ">>>>>>>>> make_bdt_uml_inputs.o: make bdt effs" << endl;
  cout << "-bins ### 1 \t ### is a sequence of numbers describing the bdt binning vector. e.g. -bins -0.2 0.2 0.5 1" << endl;
  cout << endl;
  cout << ">>>>>>>>> main_fitData.o: fits events with pdf given by main_pdf_choise or main_simul_maker" << endl;
  cout << "-i #filename \t input for fitting events (MANDATORY)" << endl;
  cout << "-t treename (default bdt)" << endl;
  cout << "-cha {0, 1} \t barrel OR endcap input, incompatible with -simul" << endl;
  cout << "\t -simul # \t simultaneous fit of # eta channels" << endl;
  cout << "\t -simul_bdt # \t simultaneous fit of # bdt channels" << endl;
  cout << "-simul_all # \t simultaneous fit of # eta and bdt channels" << endl;
  cout << "-BF {1,2,3} \t imposing the same BF in each channel for bs only (1) or for bs and bd (2) or for bs/bd and bd (3)" << endl;
  cout << "-SM \t with SM constraints (incompatible with -bd_const)" << endl;
  cout << "-bd_const \t with Bd constrainted to Bs, over all different channels (incompatible with -SM)" << endl;
  cout << "-print \t save the fits to gif and pdf --> -no_legend without parameters on canvas" << endl;
  cout << "-SB \t fit side-bands only" << endl; /// test
  cout << "-sig # \t evaluate significance with method:" << endl << "\t\t 0 by hand; " << endl << "\t\t 1 ProfileLikelihoodCalculator; " << endl << "\t\t 2 frequentist ProfileLikelihoodTestStat; " << endl << "\t\t 3 hybrid ProfileLikelihoodTestStat" << endl << "\t\t 4 hybrid RatioOfProfiledLikelihoodsTestStat" << endl;
  cout << "-e #filename \t estimates file (useful for significance)" << endl;
  cout << "-pee \t per-event-error" << endl;
  cout << "-bdt_fit \t bdt_fit" << endl;
  cout << "-free #" << endl;
//  cout << "-cuts #cutstring \t string cut for small tree: i.e. \"bdt>0.&&bdt<0.18\"; default is \"" << cuts << "\"" << endl;
  cout << "-cuts_file \t file containing bdt cuts for small tree" << endl;
  cout << "-y {0,1,all} \t year 2011, 2012 or both (this last works only with simul)" << endl;
  cout << "-asimov \t asimov dataset for significance estimation" << endl;
  cout << "-syst \t adding syst constraints" << endl;
  cout << "-randomsyst \t syst constraints are randomized" << endl;
  cout << "-rare_constr \t rare yield is constraint" << endl;
  cout << "-proof # \t enable PROOF with # workers" << endl;
  cout << "-nexp # \t number of experiments (default 1)" << endl;
  cout << "-Bd # \t significance for Bd" << endl;
  cout << "-SMIsNull # \t significance against SM (and not against 0)" << endl;
  cout << "-LLprofile" << endl;
  cout << "-hack2011 \t hack 2011 semi to 2012" << endl;
  cout << "-ws \t input workspace" << endl;
  cout << "-free # \t rare decays are not constant: 1 = semi, 2 = peak, 3 = both" << endl;
  cout << "-bdt_cut # \t minimum bdt" << endl;
  cout << "-final \t don't change input numbers, take the ones from AN" << endl;
  cout << "-berns \t linear" << endl;
  cout << "-freeze \t freeze components for significance" << endl;
  cout << "-stat \t only statistical uncertainties" << endl;
  cout << "-null \t Bs and Bd = 0" << endl;
  cout << "-LLdouble \t double profile scan" << endl;
  cout << endl;
  cout << ">>>>>>>>> main_toyMC.o: studies the pdf given by main_pdf_choise or main_simul_maker" << endl;
  cout << "-e #filename \t estimates of events file (MANDATORY)" << endl;
  cout << "-ws \t input workspace" << endl;
  cout << "-roomcs \t toy mc with RooMCStudy, otherwise by hand" << endl;
  cout << "-nexp # \t number of experiments (default 1)" << endl;
  cout << "-BF {1,2,3} \t imposing the same BF in each channel for bs only (1) or for bs and bd (2) or for bs/bd and bd (3)" << endl;
  cout << "-SM \t with SM constraints (incompatible with -bd_const)" << endl;
  cout << "-bd_const \t with Bd constrainted to Bs, over all different channels (incompatible with -SM)" << endl;
  cout << "if simultaneous: " << endl;
  cout << "\t -simul # \t simultaneous fit of # eta channels" << endl;
  cout << "\t -simul_bdt # \t simultaneous fit of # bdt channels" << endl;
  cout << "if NOT simultaneous" << endl;
  cout << "\t -pdf {bs, bd, rare, comb, total} \t combination of pdf names, for generating" << endl;
  cout << "\t -test {bs, bd, rare, comb, total} \t fitting pdf, if different from pdf" << endl;
  cout << "-bias [c+,c-,tau+,tau-,sig+,sig-]\t biasing rare pdf parameters (it works without -roomcs)" << endl;
  cout << "-pee \t per-event-error" << endl;
  cout << "-bdt_fit \t bdt_fit" << endl;
  cout << "-sig # \t evaluate significance with method:" << endl << "\t\t 0 by hand; " << endl;
  cout << "\t -Bd # \t significance for Bd" << endl;
  cout << "-y {0,1,all} \t year 2011, 2012 or both (this last works only with simul)" << endl;
  cout << "-syst \t adding syst constraints" << endl;
  cout << "-randomsyst \t syst constraints are randomized" << endl;
  cout << "-rare_constr \t rare yield is constraint" << endl;
  cout << "-hack2011 \t hack 2011 semi to 2012" << endl;
  cout << "-free # \t rare decays are not constant: 1 = semi, 2 = peak, 3 = both" << endl;
  cout << "-final \t don't change input numbers, take the ones from AN" << endl;
  cout << "-minos \t toyMc with MINOS" << endl;
  cout << "-berns \t linear" << endl;
  cout << endl;

  exit(EXIT_SUCCESS);
}

void parse_options(int argc, char* argv[]){
  
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-meth")) {
      if (!strcmp(argv[i+1],"cnc")) {
        meth = "cnc";
        method = true;
      }
      if (!strcmp(argv[i+1],"bdt")) {
        meth = "bdt";
        method = true;
      }
      cout << "method: " << meth << endl;
    }
    if (!strcmp(argv[i],"-cha")) {
      ch_s = argv[i+1];
      ch_i = atoi(ch_s.c_str());
      channel = true;
      cout << "channel: " << ch_s << endl;
    }
    if (!strcmp(argv[i],"-cha_bdt")) {
      ch_bdt_s = argv[i+1];
      ch_bdt_i = atoi(ch_bdt_s.c_str());
      channel_bdt = true;
      cout << "channel bdt: " << ch_s << endl;
    }
    if (!strcmp(argv[i],"-print")) {
      cout << "print plots" << endl;
      print = true;
    }
    if (!strcmp(argv[i],"-nexp")) {
      NExp = atoi(argv[i+1]);
      cout << "number of experiments: " << NExp << endl;
    }
    if (!strcmp(argv[i],"-i")) {
      input_name = argv[i+1];
      cout << "input = " << input_name << endl;
      input = true;
    }
    if (!strcmp(argv[i],"-ws")) {
      input_ws = argv[i+1];
      cout << "workspace = " << input_ws << endl;
    }
    if (!strcmp(argv[i],"-e")) {
      input_estimates = argv[i+1];
      cout << "estimate file = " << input_estimates << endl;
      estimate = true;
    }
    if (!strcmp(argv[i],"-pdf")) {
      pdf_toy = argv[i+1];
      cout << "pdf = " << pdf_toy << endl;
      pdf = true;
    }
    if (!strcmp(argv[i],"-test")) {
      pdf_test = argv[i+1];
      cout << "testing pdf = " << pdf_test << endl;
      pdf_test_b = true;
    }
    if (!strcmp(argv[i],"-bias")) {
      bias_s = argv[i+1];
      cout << "biasing " << bias_s << endl;
      bias = true;
    }
    if (!strcmp(argv[i],"-roomcs")) {
      cout << "using RooMCStudy" << endl;
      roomcs = true;
    }
    if (!strcmp(argv[i],"-BF")) {
      BF = atoi(argv[i+1]);
      cout << "Branching fraction option " << BF << endl;
    }
    if (!strcmp(argv[i],"-SM")) {
      cout << "SM constraints" << endl;
      SM = true;
    }
    if (!strcmp(argv[i],"-bd_const")) {
      cout << "Bd constrainted" << endl;
      bd_const = true;
    }
    if (!strcmp(argv[i],"-simul")) {
      inputs = atoi(argv[i+1]);
      cout << "simultaneous fits of " << inputs << " eta channels" << endl;
      simul = true;
    }
    if (!strcmp(argv[i],"-simul_all")) {
      inputs_all = atoi(argv[i+1]);
      cout << "simultaneous fits of " << inputs_all << " eta and bdt channels" << endl;
      simul_all = true;
    }
    if (!strcmp(argv[i],"-simul_bdt")) {
      inputs_bdt = atoi(argv[i+1]);
      cout << "simultaneous fits of " << inputs_bdt << " bdt channels" << endl;
      simul_bdt = true;
    }
    if (!strcmp(argv[i],"-t")) {
      tree_name = argv[i+1];
      cout << "tree name = " << tree_name << endl;
    }
    if (!strcmp(argv[i],"-SB")) {
      SB = true;
      cout << "fitting only the side-bands" << endl;
    }
    if (!strcmp(argv[i],"-sig")) {
      sig_meth = atoi(argv[i+1]);
      cout << "significance with method " << sig_meth << endl;
    }
//    if (!strcmp(argv[i],"-cuts")) {
//      cuts_b = true;
//      cuts = argv[i+1];
//      cout << "cut string = " << cuts << endl;
//    }
    if (!strcmp(argv[i],"-cuts_file")) {
      cuts_f_b = true;
      cuts_f = argv[i+1];
      cout << "cuts file = " << cuts_f << endl;
    }
    if (!strcmp(argv[i],"-pee")) {
      pee = true;
      cout << "per-event-error" << endl;
    }
    if (!strcmp(argv[i],"-no_legend")) {
      no_legend = true;
      cout << "no legend on canvas" << endl;
    }
    if (!strcmp(argv[i],"-rare")) {
      rare_f = argv[i+1];
      cout << "rare file = " << rare_f << endl;
    }
    if (!strcmp(argv[i],"-bdt_fit")) {
      bdt_fit = true;
      cout << "2D fit with mass and bdt" << endl;
    }
    if (!strcmp(argv[i],"-y")) {
      years_opt = argv[i+1];
      cout << "years: " << years_opt << endl;
    }
    if (!strcmp(argv[i],"-asimov")) {
      asimov = true;
      cout << "Asimov dataset" << endl;
    }
    if (!strcmp(argv[i],"-syst")) {
      syst = true;
      cout << "systematic constraints" << endl;
    }
    if (!strcmp(argv[i],"-randomsyst")) {
      randomsyst = true;
      cout << "constraints are randomized" << endl;
    }
    if (!strcmp(argv[i],"-newcomb")) {
      newcomb = true;
      cout << "new combinatorial bkg" << endl;
    }
    if (!strcmp(argv[i],"-shapesyst")) {
      shapesyst = true;
      cout << "adding shape systematics" << endl;
    }
    if (!strcmp(argv[i],"-proof")) {
      proof = atoi(argv[i+1]);
      cout << "PROOF with " << proof << " workers" << endl;
    }
    if (!strcmp(argv[i],"-Bd")) {
      Bd = true;
      cout << "significance for Bd hypothesis" << endl;
    }
    if (!strcmp(argv[i],"-SMIsNull")) {
      SMIsNull = true;
      cout << "SM is the null hypothesis" << endl;
    }
    if (!strcmp(argv[i],"-LLprofile")) {
      LLprofile = true;
      cout << "LLprofile" << endl;
    }
    if (!strcmp(argv[i],"-hack2011")) {
      hack_semi2011 = true;
      cout << "hack2011" << endl;
    }
    if (!strcmp(argv[i],"-hack")) {
      hack = true;
      cout << "hack" << endl;
    }
    if (!strcmp(argv[i],"-rare_constr")) {
    	rare_constr = true;
      cout << "rare yield is constrained, so it is with systematics" << endl;
      syst = true;
    }
    if (!strcmp(argv[i],"-rkeys")) {
    	rkeys = true;
      cout << "RooKeysPdf for MassRes and BDT" << endl;
    }
    if (!strcmp(argv[i],"-free")) {
      free_rare = atoi(argv[i+1]);
      cout << "free rare decays, option " << free_rare << endl;
    }
    if (!strcmp(argv[i],"-bdtbins")) {
    	make_bdt_binning_inputs = true;
      cout << "saves bdt histos to " << bdtbinnings_s << endl;
    }
    if (!strcmp(argv[i],"-bdt_cut")) {
    	bdt_cut = atof(argv[i+1]);
    	cout << "bdt_cut = " << bdt_cut << endl;
    }
    if (!strcmp(argv[i],"-semi")) {
    	semi = true;
      cout << "semi mass pdf from mc data" << endl;
    }
    if (!strcmp(argv[i],"-bdtsplit")) {
    	bdtsplit = true;
      cout << "bdt splitted in 1,2,3,4" << endl;
    }
    if (!strcmp(argv[i],"-final")) {
    	final = true;
    	cout << "numbers from AN" << endl;
    }
    if (!strcmp(argv[i],"-minos")) {
    	minos = true;
    	cout << "MINOS" << endl;
    }
    if (!strcmp(argv[i],"-berns")) {
    	berns = true;
    	cout << "bernstein" << endl;
    }
    if (!strcmp(argv[i],"-freeze")) {
    	freeze = true;
    	cout << "freezes components for significance" << endl;
    }
    if (!strcmp(argv[i],"-stat")) {
    	stat_only = true;
    	cout << "only statistica uncertainties" << endl;
    }
    if (!strcmp(argv[i],"-null")) {
    	null = true;
    	cout << "Bs and Bd = 0" << endl;
    }
    if (!strcmp(argv[i],"-LLdouble")) {
    	LLdouble = true;
    	cout << "double profile" << endl;
    }
    if (!strcmp(argv[i],"-h")) help();
  }
}

void parse_input (string input) {
  string a_meth_s[2] = {"cnc", "bdt"};
  string a_ch_s[2] = {"0", "1"};
  size_t found;
  for (int i = 0; i < 2; i++) {
    found = input.find(a_meth_s[i]);
    if (found!=string::npos) {
      meth = a_meth_s[i];
      cout << "meth = " << meth << endl;
    }
  }
  for (int i = 0; i < 2; i++) {
    found = input.find(a_ch_s[i]);
    if (found!=string::npos) {
      ch_s = a_ch_s[i];
      cout << "channel = " << ch_s << endl;
    }
  }
  found = input.find("SM");
  if (found!=string::npos) {
    SM = true;
    cout << "SM" << endl;
  }
  found = input.find("BdConst");
  if (found!=string::npos) {
    bd_const = true;
    cout << "bd constrained" << endl;
  }
  found = input.find("simul");
  if (found!=string::npos) {
    simul = true;
    size_t found2;
    found2 = input.find_first_of("0123456789");
    ostringstream number;
    number << input[found2];
    inputs = atoi(number.str().c_str());
    cout << "simultaneous " << inputs << endl;
  }
}

vector <double> parse_bdt_cuts() {

	vector <double> cuts_v(4, -1);
	for (int y = 0; y < 2; y++) {
		int year = y == 0 ? 2011 : 2012;
		string filename(Form("../uml/input/%d/anaBmm.plotResults.%d.tex", year, year));
		FILE *file = fopen(filename.c_str(), "r");
		if (!file) {cout << "file " << filename << " does not exist"; exit(1);}
		char buffer[1024];
		char left[1024];
		double number;

		string lefty[2];
		lefty[0] = Form("\\vdef{%d:bdt:0}", year);
		lefty[1] = Form("\\vdef{%d:bdt:1}", year);


		while (fgets(buffer, sizeof(buffer), file)) {
			if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
			if (buffer[0] == '\045') continue;
			sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);

			string left_s(left);
			for (int i = 0; i < 2; i++) {
				size_t found = left_s.find(lefty[i]);
				if (found != string::npos) {
					cuts_v[2*y+i] = number;
        }
      }
    }
		fclose(file);
  }
  cout << " bdt cuts:" << endl;
  for (int i = 0; i < 4; i++) cout << "channel " << i << "; bdt > " << cuts_v[i] << endl;
  return cuts_v;
}

vector <double> cut_bdt_file() {
  FILE *file = fopen(cuts_f.c_str(), "r");
  char buffer[1024];
  vector <double> bdt_(4, -1);
  if (!file) {
  	cout << "file " << cuts_f << " does not exist, parsing" << endl;
  	return parse_bdt_cuts();
  }
  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '#') continue;
    sscanf(buffer, "bdt_0\t%lf", &bdt_[0]);
    sscanf(buffer, "bdt_1\t%lf", &bdt_[1]);
    sscanf(buffer, "bdt_2\t%lf", &bdt_[2]);
    sscanf(buffer, "bdt_3\t%lf", &bdt_[3]);
  }
  cout << "bdt_0 cut = " << bdt_[0] << "; bdt_1 cut = " << bdt_[1] << endl;
  cout << "bdt_2 cut = " << bdt_[2] << "; bdt_3 cut = " << bdt_[3] << endl;
  fclose(file);
  return bdt_;
}



vector < double > get_singlerare_normalization(string filename, int endcap, int size) {

  string peakdecays[] = {"bgBs2KK", "bgBs2KPi", "bgBs2PiPi", "bgBd2KK", "bgBd2KPi", "bgBd2PiPi", "bgLb2PiP", "bgLb2KP"};
  string semidecays[] = {"bgBs2KMuNu", "bgBd2PiMuNu", "bgLb2PMuNu", "bgBu2PiMuMu", "bgBu2KMuMu", "bgBd2Pi0MuMu", "bgBd2K0MuMu", "bgBd2MuMuGamma", "bgBs2MuMuGamma"};


  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {cout << "file " << filename << " does not exist"; exit(1);}

  char buffer[1024];
  char left[1024];
  double number;
  int peak_n = sizeof(peakdecays)/sizeof(string);
  int semi_n = sizeof(semidecays)/sizeof(string);

  Double_t peak_t[peak_n][3];
  Double_t semi_t[semi_n][3];

  string end[3] = {Form("bsRare%d}", endcap), Form("bdRare%d}", endcap), Form("loSideband%d:val}", endcap)};

  vector <double> output(size, 0.);

  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);
    string left_s(left);
    for (int i = 0; i < peak_n; i++) {
      size_t found = left_s.find(peakdecays[i]);
      if (found != string::npos) {
        for (int j = 0; j < 3; j++) {
          found = left_s.find(end[j]);
          if (found != string::npos) {
//            output[i+2] += number;
            peak_t[i][j] = number;
          }
        }
      }
    }
    for (int i = 0; i < semi_n; i++) {
      size_t found = left_s.find(semidecays[i]);
      if (found != string::npos) {
        for (int j = 0; j < 3; j++) {
          found = left_s.find(end[j]);
          if (found != string::npos) {
//            output[i+2+peak_n] += number;
            semi_t[i][j] = number;
          }
        }
      }
    }
  }
  fclose(file);

  for (int i = 0; i < peak_n; i++) {
  	for (int j = 0; j < 3; j++) {
  		output[i+2] += peak_t[i][j];
  	}
  }
  for (int i = 0; i < semi_n; i++) {
  	for (int j = 0; j < 3; j++) {
  		output[i+2+peak_n] += semi_t[i][j];
  	}
  }

  output[0] = 1.;
  output[1] = 1.;
  output[size-1] = 1.;
  cout << " rare expected:" << endl;
  for (int i = 0; i < peak_n; i++) cout << peakdecays[i] << " = " << output[i+2] << "; ";
  cout << endl;
  for (int i = 0; i < semi_n; i++) cout << semidecays[i] << " = " << output[i+2+peak_n] << "; ";
  cout << endl << endl;
  return output;
}

