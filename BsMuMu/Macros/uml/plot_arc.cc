/*
 * plot_arc.cc
 *
 *  Created on: 28/giu/2013
 *      Author: lucamartini
 */

#include "CommonFun.h"

#include "TPaveStats.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPDF.h"
#include "TObjString.h"
#include "TH2D.h"

using namespace std;

void print_2D(TH2D h) {
	TCanvas  canvas("canvas");
	h.Draw("COLZ");
	canvas.Print(Form("./plots/%s.pdf", h.GetName()));
	canvas.Print(Form("./plots/%s.gif", h.GetName()));
}

void print_prof(TH1D * prof, double r_min, double r_max) {
	TCanvas  canvas("canvas");
//	prof->Fit("pol1", "", "", r_min, r_max);
	prof->Fit("pol0", "", "", r_min, r_max);
	prof->Draw();
	canvas.Update();
  TPaveStats *st = (TPaveStats*)prof->FindObject("stats");
  st->SetY1NDC(0.87);
	canvas.Print(Form("./plots/%s.pdf", prof->GetName()));
	canvas.Print(Form("./plots/%s.gif", prof->GetName()));
}

int main(int argc, char* argv[]) {

	parse_options(argc, argv);

	TH2D m_bdt("m_bdt", "m_bdt;m;bdt", 100, 4.9, 5.9, 100, -0.7, 0.5);
	TH2D bdt_me("bdt_me", "bdt_me;bdt;me", 100, -0.7, 0.5, 100, 0.01755, 0.1);
	TH2D bdt_rme("bdt_rme", "bdt_rme;bdt;rme", 100, -0.7, 0.5, 100, 0.003, 0.0175);
	TH2D m_me("m_me", "m_me;m;me", 100, 4.9, 5.9, 100, .025, .1);
	TH2D m_rme("m_rme", "m_rme;m;rme", 100, 4.9, 5.9, 100, 0.003, 0.0175);

	TFile smalltree_f("../uml/input/2012/small-SgData.root");
	TTree* reduced_tree = (TTree*)smalltree_f.Get("SgData_bdt");
	Double_t m_t, eta_t, m1eta_t, m2eta_t, bdt_t, me_t, m1bdt, m2bdt;
	Bool_t muid_t;
	reduced_tree->SetBranchAddress("m",     &m_t);
	reduced_tree->SetBranchAddress("bdt",   &bdt_t);
	reduced_tree->SetBranchAddress("eta",   &eta_t);
	reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
	reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
	reduced_tree->SetBranchAddress("me",    &me_t);
	reduced_tree->SetBranchAddress("muid",  &muid_t);
	reduced_tree->SetBranchAddress("m1bdt",  &m1bdt);
	reduced_tree->SetBranchAddress("m2bdt",  &m2bdt);
	double entries = reduced_tree->GetEntries();
	for (int j = 0; j < entries; j++) {
		reduced_tree->GetEntry(j);
		Bool_t barrel = fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4;
//		if (barrel) continue;
		if (m_t > 5.20 && m_t < 5.45) continue; // skip signal windows for comb bkg
		if (m_t < 4.90 || m_t > 5.90 ) continue; // skip outside range
		if (me_t < 0.0 || me_t > 0.2) continue; //skip wrong mass error
		if (bdt_t < bdt_cut) continue;
		if (m1bdt < 0.360 || m2bdt < 0.360) continue;
		Double_t rme_t = me_t/m_t;
		m_bdt.Fill(m_t, bdt_t);
		bdt_me.Fill(bdt_t, me_t);
		m_me.Fill(m_t, me_t);
		bdt_rme.Fill(bdt_t, rme_t);
		m_rme.Fill(m_t, rme_t);
	}

	TH1D *bdt_ID = new TH1D("bdt_ID", "bdt_ID", 50, -1, 1);
	TH1D *bdt_antiID = new TH1D("bdt_antiID", "bdt_antiID", 50, -1, 1);
	TH1D *bdt_ratioID = new TH1D("bdt_ratioID", "bdt_ratioID", 50, -1, 1);
	bdt_ID->Sumw2();
	bdt_antiID->Sumw2();

	TH1D *me_ID = new TH1D("me_ID", "me_ID", 50, 0.01755, 0.1);
	TH1D *me_antiID = new TH1D("me_antiID", "me_antiID", 50, 0.01755, 0.1);
	TH1D *me_ratioID = new TH1D("me_ratioID", "me_ratioID", 50, 0.01755, 0.1);
	me_ID->Sumw2();
	me_antiID->Sumw2();
	TH1D *me_inbdt_ID = new TH1D("me_inbdt_ID", "me_inbdt_ID", 50, 0.01755, 0.1);
	TH1D *me_inbdt_antiID = new TH1D("me_inbdt_antiID", "me_inbdt_antiID", 50, 0.01755, 0.1);
	TH1D *me_inbdt_ratioID = new TH1D("me_inbdt_ratioID", "me_inbdt_ratioID", 50, 0.01755, 0.1);
	me_inbdt_ID->Sumw2();
	me_inbdt_antiID->Sumw2();

	TH1D *rme_ID = new TH1D("rme_ID", "rme_ID", 50, 0.003, 0.0175);
	TH1D *rme_antiID = new TH1D("rme_antiID", "rme_antiID", 50, 0.003, 0.0175);
	TH1D *rme_ratioID = new TH1D("rme_ratioID", "rme_ratioID", 50, 0.003, 0.0175);
	rme_ID->Sumw2();
	rme_antiID->Sumw2();
	TH1D *rme_inbdt_ID = new TH1D("rme_inbdt_ID", "rme_inbdt_ID", 50, 0.003, 0.0175);
	TH1D *rme_inbdt_antiID = new TH1D("rme_inbdt_antiID", "rme_inbdt_antiID", 50, 0.003, 0.0175);
	TH1D *rme_inbdt_ratioID = new TH1D("rme_inbdt_ratioID", "rme_inbdt_ratioID", 50, 0.003, 0.0175);
	rme_inbdt_ID->Sumw2();
	rme_inbdt_antiID->Sumw2();


	TH1D *m_ID = new TH1D("m_ID", "m_ID", 50, 4.9, 5.9);
	TH1D *m_antiID = new TH1D("m_antiID", "m_antiID", 50, 4.9, 5.9);
	TH1D *m_ratioID = new TH1D("m_ratioID", "m_ratioID", 50, 4.9, 5.9);
	m_ID->Sumw2();
	m_antiID->Sumw2();

	for (int j = 0; j < entries; j++) {
		reduced_tree->GetEntry(j);
		Bool_t barrel = fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4;
		if (m_t > 5.20 && m_t < 5.45) continue; // skip signal windows for comb bkg
		if (m_t < 4.90 || m_t > 5.90 ) continue; // skip outside range
		if (me_t < 0.0 || me_t > 0.2) continue; //skip wrong mass error
		Double_t rme_t = me_t/m_t;
		if (bdt_t > bdt_cut) {
			if (m1bdt < 0.20 || m2bdt < 0.20) continue;
			if (m1bdt < 0.360 || m2bdt < 0.360) {
				bdt_antiID->Fill(bdt_t);
				me_antiID->Fill(me_t);
				m_antiID->Fill(m_t);
				rme_antiID->Fill(rme_t);
			}
			else {
				bdt_ID->Fill(bdt_t);
				me_ID->Fill(me_t);
				m_ID->Fill(m_t);
				rme_ID->Fill(rme_t);
			}
		}
		if (m1bdt > 0.360 && m2bdt > 0.360) {
			if (bdt_t < -0.2) continue;
			if (bdt_t > bdt_cut) {
				me_inbdt_ID->Fill(me_t);
				rme_inbdt_ID->Fill(rme_t);
			}
			else {
				me_inbdt_antiID->Fill(me_t);
				rme_inbdt_antiID->Fill(rme_t);
			}
		}
	}

//	gStyle->SetOptFit(0111);
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(0);

	print_2D(m_bdt);
	print_2D(bdt_me);
	print_2D(m_me);
	print_2D(bdt_rme);
	print_2D(m_rme);

	TProfile * m_bdt_x = m_bdt.ProfileX();
	print_prof(m_bdt_x, 4.9, 5.9);
	TProfile * m_bdt_y = m_bdt.ProfileY();
	print_prof(m_bdt_y, -1., 1.);

	TProfile * bdt_me_x = bdt_me.ProfileX();
	print_prof(bdt_me_x, -.4, .5);
//	print_prof(bdt_me_x, -1, 1);
	TProfile * bdt_me_y = bdt_me.ProfileY();
	print_prof(bdt_me_y, 0.01755, .1);

	TProfile * m_me_x = m_me.ProfileX();
	print_prof(m_me_x, 4.9, 5.9);
	TProfile * m_me_y = m_me.ProfileY();
	print_prof(m_me_y, 0.01755, .1);

	TProfile * m_rme_x = m_rme.ProfileX();
	print_prof(m_rme_x, 4.9, 5.9);
	TProfile * m_rme_y = m_rme.ProfileY();
	print_prof(m_rme_y, 0.003, 0.0175);

	bdt_ID->Scale(1./bdt_ID->Integral());
	bdt_antiID->Scale(1./bdt_antiID->Integral());
	bdt_ratioID->Divide(bdt_ID, bdt_antiID);
	print_prof(bdt_ratioID, -1, 1);
	TCanvas c_bdt("c_bdt");
	bdt_ID->Draw();
	bdt_antiID->SetLineColor(kRed);
	bdt_antiID->Draw("same");
	c_bdt.BuildLegend();
	c_bdt.Print("plots/bdt_antiID.gif");
	c_bdt.Print("plots/bdt_antiID.pdf");

	me_ID->Scale(1./me_ID->Integral());
	me_antiID->Scale(1./me_antiID->Integral());
	me_ratioID->Divide(me_ID, me_antiID);
	print_prof(me_ratioID, 0.01755, 0.1);
	TCanvas c_me("c_me");
	me_ID->Draw();
	me_antiID->SetLineColor(kRed);
	me_antiID->Draw("same");
	c_me.BuildLegend();
	c_me.Print("plots/me_antiID.gif");
	c_me.Print("plots/me_antiID.pdf");

	me_inbdt_ID->Scale(1./me_inbdt_ID->Integral());
	me_inbdt_antiID->Scale(1./me_inbdt_antiID->Integral());
	me_inbdt_ratioID->Divide(me_inbdt_ID, me_inbdt_antiID);
	print_prof(me_inbdt_ratioID, 0.01755, 0.1);
	TCanvas c_inbdt_me("c_inbdt_me");
	me_inbdt_ID->Draw();
	me_inbdt_antiID->SetLineColor(kRed);
	me_inbdt_antiID->Draw("same");
	c_inbdt_me.BuildLegend();
	c_inbdt_me.Print("plots/me_inbdt_antiID.gif");
	c_inbdt_me.Print("plots/me_inbdt_antiID.pdf");


	m_ID->Scale(1./m_ID->Integral());
	m_antiID->Scale(1./m_antiID->Integral());
	m_ratioID->Divide(m_ID, m_antiID);
	print_prof(m_ratioID, 4.9, 5.9);
	TCanvas c_m("c_m");
	m_ID->Draw();
	m_antiID->SetLineColor(kRed);
	m_antiID->Draw("same");
	c_m.BuildLegend(0.7, 0.3, 0.9, 0.4);
	c_m.Print("plots/m_antiID.gif");
	c_m.Print("plots/m_antiID.pdf");


	rme_inbdt_ID->Scale(1./rme_inbdt_ID->Integral());
	rme_inbdt_antiID->Scale(1./rme_inbdt_antiID->Integral());
	rme_inbdt_ratioID->Divide(rme_inbdt_ID, rme_inbdt_antiID);
	print_prof(rme_inbdt_ratioID, 0.0045, 0.0175);
	TCanvas c_inbdt_rme("c_inbdt_rme");
	rme_inbdt_ID->Draw();
	rme_inbdt_antiID->SetLineColor(kRed);
	rme_inbdt_antiID->Draw("same");
	c_inbdt_rme.BuildLegend();
	c_inbdt_rme.Print("plots/rme_inbdt_antiID.gif");
	c_inbdt_rme.Print("plots/rme_inbdt_antiID.pdf");

	int colors[4] = {kBlue, kRed, kBlack, kGreen};

	int bins[5] = {42, 47, 50, 100};

	TH1D * m_bdt_px[4];
	TCanvas * m_bdt_px_c = new TCanvas("m_bdt_px_c");
	for (int i = 0; i < 3; i++) {
		m_bdt_px[i]  = m_bdt.ProjectionX(Form("_px%d", i), bins[i], bins[i+1], "e");
		m_bdt_px[i]->SetLineColor(colors[i]);
		m_bdt_px[i]->Rebin(4);
		m_bdt_px[i]->Scale(1./m_bdt_px[i]->Integral());
		m_bdt_px[i]->Draw(i == 0 ? "" : "same");
		if (i > 0) cout << "KS " << m_bdt_px[i]->KolmogorovTest(m_bdt_px[i-1]) << endl;
		cout << "bdt range = " << m_bdt.GetYaxis()->GetBinCenter(bins[i]) <<  "  " << m_bdt.GetYaxis()->GetBinCenter(bins[i+1]) << endl;
	}
	m_bdt_px_c->Print("./plots/m_bdt_px.pdf");

	TH1D * bdt_rme_px[4];
	TCanvas * bdt_rme_px_c = new TCanvas("bdt_rme_px_c");
	for (int i = 0; i < 3; i++) {
		bdt_rme_px[i]  = bdt_rme.ProjectionX(Form("_px%d", i), bins[i], bins[i+1], "e");
		bdt_rme_px[i]->SetLineColor(colors[i]);
		bdt_rme_px[i]->Rebin(4);
		bdt_rme_px[i]->Scale(1./bdt_rme_px[i]->Integral());
		bdt_rme_px[i]->Draw(i == 0 ? "" : "same");
		if (i > 0) cout << "KS " << bdt_rme_px[i]->KolmogorovTest(bdt_rme_px[i-1]) << endl;
	}
	bdt_rme_px_c->Print("./plots/bdt_rme_px.pdf");

	TH1D * bdt_rme_py[4];
	TCanvas * bdt_rme_py_c = new TCanvas("bdt_rme_py_c");
	for (int i = 0; i < 3; i++) {
		bdt_rme_py[i]  = bdt_rme.ProjectionY(Form("_py%d", i), bins[i], bins[i+1], "e");
		bdt_rme_py[i]->SetLineColor(colors[i]);
		bdt_rme_py[i]->Rebin(4);
		bdt_rme_py[i]->Scale(1./bdt_rme_py[i]->Integral());
		bdt_rme_py[i]->Draw(i == 0 ? "" : "same");
		if (i > 0) cout << "KS " << bdt_rme_py[i]->KolmogorovTest(bdt_rme_py[i-1]) << endl;
		cout << "bdt range = " << bdt_rme.GetXaxis()->GetBinCenter(bins[i]) <<  "  " << bdt_rme.GetXaxis()->GetBinCenter(bins[i+1]) << endl;
	}
	bdt_rme_py_c->Print("./plots/bdt_rme_py.pdf");

	TH1D * m_rme_px[4];
	TCanvas * m_rme_px_c = new TCanvas("m_rme_px_c");
	for (int i = 0; i < 3; i++) {
		m_rme_px[i]  = m_rme.ProjectionX(Form("_px%d", i), bins[i], bins[i+1], "e");
		m_rme_px[i]->SetLineColor(colors[i]);
		m_rme_px[i]->Rebin(4);
		m_rme_px[i]->Scale(1./m_rme_px[i]->Integral());
		m_rme_px[i]->Draw(i == 0 ? "" : "same");
		if (i > 0) cout << "KS " << m_rme_px[i]->KolmogorovTest(m_rme_px[i-1]) << endl;
	}
	m_rme_px_c->Print("./plots/m_rme_px.pdf");

//	cout << rme_inbdt_ID->GetEntries() << endl;
//	cout << rme_inbdt_antiID->GetEntries() << endl;
	return EXIT_SUCCESS;
}

