/*
 * make_bdt_uml_inputs.cc
 *
 *  Created on: 12/giu/2013
 *      Author: lucamartini
 */

#include "CommonFun.h"

using namespace std;

int main(int argc, char* argv[]) {

	parse_options(argc, argv);

	if (cuts_f=="no") {
		cout << "insert initial cut file: -cuts_file input/channel_cuts.txt" << endl;
		return 1;
	}
	if (!simul_all && !bdt_fit) {
		cout << "choose between -simul_all 12 or -bdt_fit" << endl;
		return 1;
	}

	string input_s("input/bdtbins.txt");
	if (bdt_fit) input_s = "input/2dbins.txt";

	FILE *file = fopen(input_s.c_str(), "r");
	if (!file) {
		cout << input_s << " does not exist" << endl;
		return 1;
	}
	char buffer[1024];

	vector < vector <double> > binnings(4);

	int in = 0;
	 while (fgets(buffer, sizeof(buffer), file)) {
		 cout << "binning " << buffer << endl;
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
			 binnings[in].push_back(n);
			 if (n > 0.999) not_ended = false;
		 }
		 in++;
	 }
	 fclose(file);

	 for (unsigned int i = 0; i < binnings.size(); i++) {
		 for (unsigned int j = 0; j < binnings[i].size(); j++) {
			 cout << binnings[i][j] << " ";
		 }
		 cout << endl;
	 }

	TFile* bdt_f = new TFile("../uml/input/bdtbinnings.root");
  vector <string> type_s;
  type_s.push_back("bs");
  type_s.push_back("bd");
  type_s.push_back("peak");
  type_s.push_back("semi");
  type_s.push_back("comb");
  type_s.push_back("bu");

  vector <double> cuts_v(binnings.size(), -1);
  if (cuts_f_b) cuts_v = cut_bdt_file();

  string full_output = "input/eff_";
  string ini_output = "input/eff_";
  string fact_output = "input/eff_";
  if (bdt_fit) {
  	full_output += "2d";
  	ini_output += "2d";
  	fact_output += "2d";
  }
  else {
  	full_output += "bdtbins";
  	ini_output += "bdtbins";
  	fact_output += "bdtbins";
  }
	full_output += "_full.txt";
	ini_output += "_ini.txt";
	fact_output += "_fact.txt";

  FILE* file_out = fopen(full_output.c_str(), "w");
  FILE* file_ini = fopen(ini_output.c_str(), "w");
  FILE* file_fact = fopen(fact_output.c_str(), "w");
  for (unsigned int i = 0; i < binnings.size(); i++) {
  	for (unsigned int j = 0; j < type_s.size(); j++) {
  		TH1D* bdt_histo = (TH1D*)bdt_f->Get(Form("bdt_%s_%d_h", type_s[j].c_str(), i));
  		Double_t tot = bdt_histo->Integral();
  		vector <double> eff(binnings[i].size() - 1);
  		for (unsigned int k = 0; k < binnings[i].size() - 1; k++) {
  			Int_t bdt0 = bdt_histo->FindBin(binnings[i][k]);
  			Int_t bdt1 = bdt_histo->FindBin(binnings[i][k+1]);
  			eff[k] = bdt_histo->Integral(bdt0, bdt1) / tot;
  		}
  		cout << bdt_histo->GetName() << "\t";
  		for (unsigned int k = 0; k < binnings[i].size() - 1; k++) {
  			cout << "[" << binnings[i][k] << " - " << binnings[i][k+1] << "] --> " << eff[k] << "\t";
  			fprintf(file_out, "bdt_%s_%d_%d\t%lf\n", type_s[j].c_str(), i, k, eff[k]);
  		}
  		Int_t bdt_ini = bdt_histo->FindBin(cuts_v[i]);
  		Double_t eff_ini = bdt_histo->Integral(bdt_ini, bdt_histo->GetNbinsX()) / tot;
  		fprintf(file_ini, "bdt_%s_%d\t%lf\n", type_s[j].c_str(), i, eff_ini);
  		cout << endl;

  		for (unsigned int k = 0; k < binnings[i].size() - 1; k++) {
  			Double_t fact = eff[k]/eff_ini;
  			fprintf(file_fact, "bdt_%s_%d_%d\t%lf\n", type_s[j].c_str(), i, k, fact);
  		}
  	}
  }
  bdt_f->Close();

  fclose(file_out);
  cout << endl <<  full_output << " contains" << endl;
  system(Form("cat %s", full_output.c_str()));

  fclose(file_ini);
  cout << endl <<  ini_output << " contains" << endl;
  system(Form("cat %s", ini_output.c_str()));

  fclose(file_fact);
  cout << endl <<  fact_output << " contains" << endl;
  system(Form("cat %s", fact_output.c_str()));

  return EXIT_SUCCESS;
}


