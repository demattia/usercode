#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF2.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector < string> file_name(2);
  file_name[0] = "input/2011/small-SgData.root";
  file_name[1] = "input/2011/small-SgMc.root";

  vector <TFile*> file_f(2);
  int channels = 2;
  vector <TH2D*> mass_vs_bdt_data_h(2);
  vector <TH2D*> mass_vs_bdt_mc_h(2);
  for (int i = 0; i < channels; i ++) {
    mass_vs_bdt_data_h[i] = new TH2D(Form("mass_vs_bdt_data_h_%d", i), Form("mass_vs_bdt_data_h_%d;mass[GeV];bdt", i), 100, 4.9, 5.9, 100, 0., .5);
    mass_vs_bdt_mc_h[i] = new TH2D(Form("mass_vs_bdt_mc_h_%d", i), Form("mass_vs_bdt_mc_h_%d;mass[GeV];bdt", i), 100, 4.9, 5.9, 100, 0., .5);
    file_f[i] = new TFile(file_name[i].c_str());
  }
  vector < string> tree_name(2);
  tree_name[0] = "SgData_bdt";
  tree_name[1] = "SgMc_bdt";

  /// mass vs bdt
  for (int h = 0; h < 2; h++) {
    TTree* tree = (TTree*)file_f[h]->Get(tree_name[h].c_str());
    double m1eta_t, m2eta_t, m_t, bdt_t;
    tree->SetBranchAddress("m1eta", &m1eta_t);
    tree->SetBranchAddress("m2eta", &m2eta_t);
    tree->SetBranchAddress("m", &m_t);
    tree->SetBranchAddress("bdt", &bdt_t);
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if (m_t > 4.9 && m_t < 5.9) {
        if (bdt_t > 0.) {
          if (fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4) {
            if (tree_name[h] == "SgData_bdt") {
              mass_vs_bdt_data_h[0]->Fill(m_t, bdt_t);
            }
            if (tree_name[h] == "SgMc_bdt") {
              mass_vs_bdt_mc_h[0]->Fill(m_t, bdt_t);
            }
          }
          else {
            if (tree_name[h] == "SgData_bdt") {
              mass_vs_bdt_data_h[1]->Fill(m_t, bdt_t);
            }
            if (tree_name[h] == "SgMc_bdt") {
              mass_vs_bdt_mc_h[1]->Fill(m_t, bdt_t);
            }
          }
        }
      }
    }
  }



  for (int i = 0; i < channels; i ++) {
    TCanvas * mass_vs_bdt_c = new TCanvas("mass_vs_bdt_c", "mass_vs_bdt_c", 600, 600);

    mass_vs_bdt_mc_h[i]->SetStats(0);
    mass_vs_bdt_mc_h[i]->SetMarkerStyle(7);
    mass_vs_bdt_mc_h[i]->SetMarkerColor(kRed);
    mass_vs_bdt_mc_h[i]->SetFillColor(kRed);

    mass_vs_bdt_data_h[i]->SetStats(0);
    mass_vs_bdt_data_h[i]->SetMarkerStyle(7);
    mass_vs_bdt_data_h[i]->SetMarkerColor(kBlack);
    mass_vs_bdt_data_h[i]->SetFillColor(kBlack);

    mass_vs_bdt_data_h[i]->SetMaximum(10);

    mass_vs_bdt_mc_h[i]->Draw("");
//    mass_vs_bdt_data_h[i]->Draw("same");

    TLine *line;
    if (i == 0)   {
      line = new TLine(4.9,0.18,5.9,0.18);
    }
    else line = new TLine(4.9,0.19,5.9,0.19);
    line->SetLineWidth(3);
    line->SetLineColor(kViolet);
    line->Draw();
    mass_vs_bdt_c->Print(Form("fig/mass_vs_bdt_%d.gif", i));
    mass_vs_bdt_c->Print(Form("fig/mass_vs_bdt_%d.pdf", i));
//    mass_vs_bdt_c->Print(Form("fig/mass_vs_bdt_%d.root", i));
    delete mass_vs_bdt_c;
  }

}
