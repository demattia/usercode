/*
 * get_quantiles.cc
 *
 *  Created on: 04/lug/2013
 *      Author: lucamartini
 */

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

using namespace std;

void get_quantiles() {

	TFile *_file0 = TFile::Open("input/2011/small-SgMc.root");
	TTree *_tree0 = (TTree*)_file0->Get("SgMc_bdt");
	TH1D  *_h0 = new TH1D("_h0", "_h0", 900, 0.1, 1);
	_tree0->Draw("bdt>>_h0", "m>4.9 && m<5.9 && bdt > 0.1 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4)");

	const Int_t nq = 4;
	Double_t xq[nq];
	Double_t yq[nq];
	for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
	_h0->GetQuantiles(nq,yq,xq);

	for (int i = 0; i< nq; i++) cout << yq[i] << endl;
	cout << endl;

	Double_t zq[nq+1];
	zq[0] = 0;
	for (int i = 0; i< nq; i++) zq[i+1] = yq[i];
	for (int i = 0; i< nq; i++) cout << _h0->Integral(_h0->FindBin(zq[i]), _h0->FindBin(zq[i+1])) / _h0->Integral() << endl;

	TFile *_file1 = TFile::Open("input/2012/small-SgData.root");
	TTree *_tree1 = (TTree*)_file1->Get("SgData_bdt");

	TH1D  *_h1 = new TH1D("_h1", "_h1", 900, 0.1, 1);
	TH1D  *_h2 = new TH1D("_h2", "_h2", 900, 0.1, 1);

	Double_t factor1 = _tree1->Draw("m", "m>4.9 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.1");
	Double_t factor2 = _tree1->Draw("m", "m>4.9 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.1");
//	Double_t large = 5.9 - 5.45 + 5.2 - 4.9;
	Double_t large = 1.;
	Double_t small = 5.9 - 5.45;
	Double_t wind = large / small;
	Double_t factor3 = factor1/factor2;
	Double_t totfact = factor3 * wind;

//	string cat1[4];
//	cat1[0] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.1 && bdt < 0.22";
//	cat1[1] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.22 && bdt < 0.33";
//	cat1[2] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.33 && bdt < 0.45";
//	cat1[3] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.45 && bdt < 1";

//	string cat2[4];
//	cat2[0] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.1 && bdt < 0.22";
//	cat2[1] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.22 && bdt < 0.33";
//	cat2[2] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.33 && bdt < 0.45";
//	cat2[3] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.45 && bdt < 1";

//	Double_t comb_new[4];
//	Double_t comb_true[4];
//	for (int i = 0; i < 4; i++) {
//		comb_new[i] = _tree1->Draw("m", cat2[i].c_str()) * totfact;
//		comb_true[i] = _tree1->Draw("m", cat1[i].c_str()) * wind;
//		cout << comb_new[i] << "  " << comb_true[i] << endl;
//	}


	string cat1[2];
	cat1[0] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.1 && bdt < 0.26";
	cat1[1] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && muid && bdt > 0.38 && bdt < 1";

	string cat2[2];
	cat2[0] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.1 && bdt < 0.26";
	cat2[1] = "m>5.45 && m<5.9 && !(abs(m1eta)<1.4&&abs(m2eta)<1.4) && m1bdt > -1 && m2bdt > -1 && !muid && bdt > 0.38 && bdt < 1";

	Double_t comb_new[2];
	Double_t comb_true[2];
	for (int i = 0; i < 2; i++) {
		comb_new[i] = _tree1->Draw("m", cat2[i].c_str()) * totfact;
		comb_true[i] = _tree1->Draw("m", cat1[i].c_str()) * wind;
		cout << comb_new[i] << "  " << comb_true[i] << endl;
	}

}


