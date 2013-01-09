#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>

TH1F * h1 = new TH1F("h1","P_{T} muon",30,0.,30.);
TH1F * h2 = new TH1F("h2","#eta muon",30,-8.,8.);
TH1F * h3 = new TH1F("h3","P_{T} Ups",30,0.,30.);
TH1F * h4 = new TH1F("h4","#eta Ups",30,-8.,8.);

void main(){

  TFile *file1 = TFile::Open("/tmp/aafke/onia2mumu_upsmumu.root");
  TTree *tree1=(TTree*) file1->Get("T1");
  fill_hists(tree1);
  makeplots();
 
}


void fill_hists( TTree* t){

  double n = t->GetEntries();
  cout << "n = " << n << endl;
  TLorentzVector mu1;
  TLorentzVector mu2;
  TClonesArray* Mc_mu_4mom=new TClonesArray("TLorentzVector");
  int Mc_mu_size;
  t->SetBranchAddress("Mc_mu_4mom", &Mc_mu_4mom);
  t->SetBranchAddress("Mc_mu_size", &Mc_mu_size);
  for (int ev=0;ev<n;ev++){
    t->GetEntry(ev);
    if(Mc_mu_size==2){
      for(int i=0; i<Mc_mu_size; i++){
   	mu1= &((TLorentzVector *)Mc_mu_4mom->At(i));
	for(int j=i+1; j<Mc_mu_size; j++){
          mu2=&((TLorentzVector *)Mc_mu_4mom->At(j));
          TLorentzVector aa=mu1;
          TLorentzVector bb=mu2;
	  /if(charge1!=charge2){
	    // Muon Pt
	    h1->Fill(aa.Pt());
	    h1->Fill(bb.Pt());
	    // Muon Eta
	    h2->Fill(aa.Eta());
	    h2->Fill(bb.Eta());
	    
	    TLorentzVector *cc=new TLorentzVector;
	    (*cc)=aa+bb;
	    // Ups Pt (NB: can also get this from Mc_QQ_4mom!)
	    h3->Fill(cc->Pt());
	    //  Ups Eta
	    h4->Fill(cc->Eta());
	  }
        }
      }
    }
  }
}

void makeplots(){
  TCanvas *c1 = new TCanvas("c1","title",800,600);
  c1->cd();
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);
  h1->SetStats(0);
  h1->GetXaxis()->SetTitle("P_{T}^{#mu} (GeV)");
  h1->GetYaxis()->SetTitle("# Entries");
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->SetLineColor(2);
  h1->SetLineStyle(1);
  h1->SetLineWidth(3);
  h1->Draw("");
  c1->Print("plot1.eps");
 
  TCanvas *c2 = new TCanvas("c2","title",800,600);
  c2->cd();
  c2->SetBottomMargin(0.15);
  c2->SetLeftMargin(0.15);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("#eta^{#mu} (GeV)");
  h2->GetYaxis()->SetTitle("# Entries");
  h2->GetYaxis()->SetTitleOffset(1.3);
  h2->SetLineColor(2);
  h2->SetLineStyle(1);
  h2->SetLineWidth(3);
  h2->Draw("same");
  c2->Print("plot2.eps");

  TCanvas *c3 = new TCanvas("c3","title",800,600);
  c3->cd();
  c3->SetBottomMargin(0.15);
  c3->SetLeftMargin(0.15);
  h3->SetStats(0);
  h3->GetXaxis()->SetTitle("P_{T}^{#mu} (GeV)");
  h3->GetYaxis()->SetTitle("# Entries");
  h3->GetYaxis()->SetTitleOffset(1.3);
  h3->SetLineColor(2);
  h3->SetLineStyle(1);
  h3->SetLineWidth(3);
  h3->Draw("");
  c3->Print("plot3.eps");
 
  TCanvas *c4 = new TCanvas("c4","title",800,600);
  c4->cd();
  c4->SetBottomMargin(0.15);
  c4->SetLeftMargin(0.15);
  h4->SetStats(0);
  h4->GetXaxis()->SetTitle("#eta^{#mu} (GeV)");
  h4->GetYaxis()->SetTitle("# Entries");
  h4->GetYaxis()->SetTitleOffset(1.3);
  h4->SetLineColor(2);
  h4->SetLineStyle(1);
  h4->SetLineWidth(3);
  h4->Draw("same");
  c4->Print("plot4.eps");

  

}

