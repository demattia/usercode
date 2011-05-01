#include <vector>
#include "TClonesArray.h"
#include "TFile.h"
#include "TCanvas.h"

void Convert3pbPPtreeToHItree_HIcuts(){
	TFile *f = new TFile("/Users/demattia/HeavyIons/MuOnia.root");
	TDirectory *dir = (TDirectory*)f->Get("detailedDimuonTree");
	TTree *t = (TTree*)dir->Get("probe_tree");
	Long64_t nentries = t->GetEntries();

	std::cout << "nentries = " << nentries << std::endl;

	const int NMAX=100;
	Int_t Centrality;
	Int_t HLTriggers;
	Int_t Reco_QQ_size;
	Int_t Reco_QQ_trig[NMAX];
	Int_t Reco_QQ_type[NMAX];
	Int_t Reco_QQ_sign[NMAX];
	Float_t Reco_QQ_VtxProb[NMAX];
	TClonesArray *Reco_QQ_4mom;
	TClonesArray *Reco_QQ_mupl_4mom;
	TClonesArray *Reco_QQ_mumi_4mom;
	// Int_t Reco_mu_size;
	// Int_t Reco_mu_type[NMAX];
	// Int_t Reco_mu_charge[NMAX];
	// TClonesArray *Reco_mu_4mom = 0;
	//Int_t nReco_QQ;
	Float_t invariantMass;
	Int_t QQtrig;
	Int_t QQsign;
	Float_t upsPt;
	Float_t upsEta;
	Float_t upsRapidity;	
	Float_t vProb;
	Float_t muPlusPt;
	Float_t muMinusPt;
	Float_t muPlusEta;
	Float_t muMinusEta;

	float QQmuPlusDxy[100];
	float QQmuPlusDz[100];
	float QQmuPlusNhits[100];
        float QQmuPlusNPixelhits[100];
	int QQmuPlusNPxlLayers[100];
	float QQmuPlusInnerChi[100];
	float QQmuPlusGlobalChi[100];
	int QQmuPlusNMuonhits[100];
	float QQmuPlusDB[100];

	float QQmuMinusDxy[100];
	float QQmuMinusDz[100];
	float QQmuMinusNhits[100];
	float QQmuMinusNPixelhits[100];
	int QQmuMinusNPxlLayers[100];
	float QQmuMinusInnerChi[100];
	float QQmuMinusGlobalChi[100];
	int QQmuMinusNMuonhits[100];
	float QQmuMinusDB[100];

	int mu1_HLT_L1DoubleMuOpen = false;
	float mu1_pt, mu1_eta, mu1_phi;
	int mu2_HLT_L1DoubleMuOpen = false;
	float mu2_pt, mu2_eta, mu2_phi;
	int mu1_GM, mu2_GM;
	int mu1_charge, mu2_charge;

	float muPlusDxy;
	float muPlusDz;
	int muPlusNhits;
	int muPlusNPixelhits;
	int muPlusNPxlLayers;
	float muPlusInnerChi;
	float muPlusGlobalChi;
	int muPlusNMuonhits;
	float muPlusDB;

	float muMinusDxy;
	float muMinusDz;
	int muMinusNhits;
	int muMinusNPixelhits;
	int muMinusNPxlLayers;
	float muMinusInnerChi;
	float muMinusGlobalChi;
	int muMinusNMuonhits;
	float muMinusDB;

	t->SetBranchAddress("Centrality",    &Centrality    );
	t->SetBranchAddress("HLTriggers",    &HLTriggers    );
	t->SetBranchAddress("Reco_QQ_size",  &Reco_QQ_size  );
	t->SetBranchAddress("Reco_QQ_trig",  &Reco_QQ_trig  );
	t->SetBranchAddress("Reco_QQ_type",  Reco_QQ_type   );
	t->SetBranchAddress("Reco_QQ_sign",  Reco_QQ_sign   );
	t->SetBranchAddress("vProb",  Reco_QQ_VtxProb);
	t->SetBranchAddress("Reco_QQ_4mom",  &Reco_QQ_4mom  );
	t->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom);
	t->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom);

	// t->SetBranchAddress("Centrality",    &Centrality    );
	// t->SetBranchAddress("HLTriggers",    &HLTriggers    );
	// t->SetBranchAddress("Reco_QQ_size",  &Reco_QQ_size  );
	// t->SetBranchAddress("Reco_QQ_trig",  &Reco_QQ_trig  );
	// t->SetBranchAddress("Reco_QQ_type",  Reco_QQ_type   );
	// t->SetBranchAddress("Reco_QQ_sign",  Reco_QQ_sign   );
	// t->SetBranchAddress("Reco_QQ_VtxProb",  Reco_QQ_VtxProb);
	// t->SetBranchAddress("Reco_QQ_4mom",  &Reco_QQ_4mom  );
	// t->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom);
	// t->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom);

	t->SetBranchAddress("mu1_dxy",QQmuPlusDxy);
	t->SetBranchAddress("mu1_dz",QQmuPlusDz);
	t->SetBranchAddress("mu1_nTrHits",QQmuPlusNhits); // Check if it is Trk+Pix
	t->SetBranchAddress("mu1_nPixHits",QQmuPlusNPixelhits);
	t->SetBranchAddress("muPlusNPxlLayers",QQmuPlusNPxlLayers); // Not there
	t->SetBranchAddress("mu1_nChi2",QQmuPlusInnerChi);
	t->SetBranchAddress("muPlusGlobalChi",QQmuPlusGlobalChi); // Not there or check if it is the above
	t->SetBranchAddress("muPlusNMuonhits",QQmuPlusNMuonhits); // Not there. There is nMuSegs though
	t->SetBranchAddress("mu1_dB",QQmuPlusDB);
	t->SetBranchAddress("mu1_HLT_L1DoubleMuOpen", &mu1_HLT_L1DoubleMuOpen);
	t->SetBranchAddress("mu1_pt", &mu1_pt);
	t->SetBranchAddress("mu1_eta", &mu1_eta);
	t->SetBranchAddress("mu1_phi", &mu1_phi);
	t->SetBranchAddress("mu1_GM", &mu1_GM);
	t->SetBranchAddress("mu1_charge", &mu1_charge);

	t->SetBranchAddress("mu2_dxy",QQmuMinusDxy);
	t->SetBranchAddress("mu2_dz",QQmuMinusDz);
	t->SetBranchAddress("mu2_nTrHits",QQmuMinusNhits); // Check if it is Trk+Pix
	t->SetBranchAddress("mu2_nPixHits",QQmuMinusNPixelhits);
	t->SetBranchAddress("muMinusNPxlLayers",QQmuMinusNPxlLayers); // Not there
	t->SetBranchAddress("mu2_nChi2",QQmuMinusInnerChi);
	t->SetBranchAddress("muMinusGlobalChi",QQmuMinusGlobalChi); // Not there or check if it is the above
	t->SetBranchAddress("muMinusNMuonhits",QQmuMinusNMuonhits); // Not there. There is nMuSegs though
	t->SetBranchAddress("mu2_dB",QQmuMinusDB);
	t->SetBranchAddress("mu2_HLT_L1DoubleMuOpen", &mu2_HLT_L1DoubleMuOpen);
	t->SetBranchAddress("mu2_pt", &mu2_pt);
	t->SetBranchAddress("mu2_eta", &mu2_eta);
	t->SetBranchAddress("mu2_phi", &mu2_phi);
	t->SetBranchAddress("mu2_GM", &mu2_GM);
	t->SetBranchAddress("mu2_charge", &mu2_charge);

	// t->SetBranchAddress("muMinusDxy",QQmuMinusDxy);
	// t->SetBranchAddress("muMinusDz",QQmuMinusDz);
	// t->SetBranchAddress("muMinusNhits",QQmuMinusNhits);
	// t->SetBranchAddress("muMinusNPixelhits",QQmuMinusNPixelhits);
	// t->SetBranchAddress("muMinusNPxlLayers",QQmuMinusNPxlLayers);
	// t->SetBranchAddress("muMinusInnerChi",QQmuMinusInnerChi);
	// t->SetBranchAddress("muMinusGlobalChi",QQmuMinusGlobalChi);
	// t->SetBranchAddress("muMinusNMuonhits",QQmuMinusNMuonhits);
	// t->SetBranchAddress("muMinusDB",QQmuMinusDB);

	// t->SetBranchAddress("Reco_mu_size",  &Reco_mu_size  );
	// t->SetBranchAddress("Reco_mu_type",  Reco_mu_type   );
	// t->SetBranchAddress("Reco_mu_charge",Reco_mu_charge );
	// t->SetBranchAddress("Reco_mu_4mom",  &Reco_mu_4mom);

	TH1F *h_QQ_mass = new TH1F("h_QQ_mass","",60,7,13);
	TH1F *h_QQ_mass_1 = new TH1F("h_QQ_mass_1","",60,7,13);
	TH1F *h_QQ_mass_2 = new TH1F("h_QQ_mass_2","",60,7,13);	

	TFile *f1 = new TFile("GlbMuon_hi_MassTree.root","RECREATE");
	TTree *UpsilonTree = new TTree("UpsilonTree","UpsilonTree");
	//UpsilonTree->Branch("nReco_QQ", &nReco_QQ, "nReco_QQ/I");

	// UpsilonTree->Branch("Centrality",    &Centrality,    "Centrality/I");
	UpsilonTree->Branch("HLTriggers",    &HLTriggers,    "HLTriggers/I");
	UpsilonTree->Branch("QQtrig",   &QQtrig,    "QQtrig/I");
	UpsilonTree->Branch("QQsign",   &QQsign,    "QQsign/I");
	UpsilonTree->Branch("vProb",    &vProb,     "vProb/F");
	UpsilonTree->Branch("invariantMass", &invariantMass, "invariantMass/F");
	UpsilonTree->Branch("upsPt", &upsPt, "upsPt/F");
	UpsilonTree->Branch("upsEta", &upsEta, "upsEta/F");
	UpsilonTree->Branch("upsRapidity", &upsRapidity, "upsRapidity/F");
	UpsilonTree->Branch("muPlusPt", &muPlusPt, "muPlusPt/F");
	UpsilonTree->Branch("muMinusPt", &muMinusPt, "muMinusPt/F");
	UpsilonTree->Branch("muPlusEta", &muPlusEta, "muPlusEta/F");
	UpsilonTree->Branch("muMinusEta", &muMinusEta, "muMinusEta/F");

	UpsilonTree->Branch("muPlusDxy", &muPlusDxy,"muPlusDxy/F");
	UpsilonTree->Branch("muPlusDz", &muPlusDz,"muPlusDz/F");
	UpsilonTree->Branch("muPlusNhits", &muPlusNhits,"muPlusNhits/I");
	UpsilonTree->Branch("muPlusNPixelhits", &muPlusNPixelhits,"muPlusNPixelhits/I");
	UpsilonTree->Branch("muPlusNPxlLayers", &muPlusNPxlLayers,"muPlusNPxlLayers/I");
	UpsilonTree->Branch("muPlusInnerChi", &muPlusInnerChi,"muPlusInnerChi/F");
	UpsilonTree->Branch("muPlusGlobalChi", &muPlusGlobalChi,"muPlusGlobalChi/F");
	UpsilonTree->Branch("muPlusNMuonhits", &muPlusNMuonhits,"muPlusNMuonhits/I");
	UpsilonTree->Branch("muPlusDB", &muPlusDB,"muPlusDB/F");

	UpsilonTree->Branch("muMinusDxy", &muMinusDxy,"muMinusDxy/F");
	UpsilonTree->Branch("muMinusDz", &muMinusDz,"muMinusDz/F");
	UpsilonTree->Branch("muMinusNhits", &muMinusNhits,"muMinusNhits/I");
	UpsilonTree->Branch("muMinusNPixelhits", &muMinusNPixelhits,"muMinusNPixelhits/I");
	UpsilonTree->Branch("muMinusNPxlLayers", &muMinusNPxlLayers,"muMinusNPxlLayers/I");
	UpsilonTree->Branch("muMinusInnerChi", &muMinusInnerChi,"muMinusInnerChi/F");
	UpsilonTree->Branch("muMinusGlobalChi", &muMinusGlobalChi,"muMinusGlobalChi/F");
	UpsilonTree->Branch("muMinusNMuonhits", &muMinusNMuonhits,"muMinusNMuonhits/I");
	UpsilonTree->Branch("muMinusDB", &muMinusDB,"muMinusDB/F");


	std::cout << "Number of entries = " << nentries << std::endl;

	for (int i=0; i<nentries; i++) {
		t->GetEntry(i);

		// if( i > 100000 ) break;

	        if( mu1_pt == 0 || mu2_pt == 0 ) continue;

		// std::cout << "mu1_pt = " << mu1_pt << std::endl;

		//nReco_QQ=0;
		if (i%1000==0) cout<<i<<endl;
		// for(int iQQ = 0; iQQ < Reco_QQ_size; iQQ++){
		// Taking only the first candidate to check
		for(int iQQ = 0; iQQ < 1; iQQ++){
			//cout<<iQQ<<endl;
			vProb = Reco_QQ_VtxProb[iQQ];
			QQtrig = Reco_QQ_trig[iQQ];
			// QQsign = Reco_QQ_sign[iQQ];
			// TLorentzVector *Reco_QQ = (TLorentzVector *) Reco_QQ_4mom->At(iQQ);
			// TLorentzVector *Reco_QQ_mupl = (TLorentzVector *) Reco_QQ_mupl_4mom->At(iQQ);
			// TLorentzVector *Reco_QQ_mumi = (TLorentzVector *) Reco_QQ_mumi_4mom->At(iQQ);

			QQsign = 1;
			if( mu1_charge != mu2_charge ) QQsign = 0;
			TLorentzVector *Reco_QQ_mupl = new TLorentzVector;
			Reco_QQ_mupl->SetPtEtaPhiM(mu1_pt, mu1_eta, mu1_phi, 0.1056583);
			TLorentzVector *Reco_QQ_mumi = new TLorentzVector;
			Reco_QQ_mumi->SetPtEtaPhiM(mu2_pt, mu2_eta, mu2_phi, 0.1056583);
			TLorentzVector *Reco_QQ = new TLorentzVector(*Reco_QQ_mupl + *Reco_QQ_mumi);

			invariantMass=Reco_QQ->M();
			upsPt=Reco_QQ->Pt();
			upsEta=Reco_QQ->Eta();
			upsRapidity=Reco_QQ->Rapidity();
			muMinusPt=Reco_QQ_mumi->Pt();
			muMinusEta=Reco_QQ_mumi->Eta();
			muPlusPt=Reco_QQ_mupl->Pt();
			muPlusEta=Reco_QQ_mupl->Eta();

			muPlusDxy=QQmuPlusDxy[iQQ];
			muPlusDz=QQmuPlusDz[iQQ];
			muPlusNhits=QQmuPlusNhits[iQQ];
			muPlusNPixelhits=QQmuPlusNPixelhits[iQQ];
			muPlusNPxlLayers=QQmuPlusNPxlLayers[iQQ];
			muPlusInnerChi=QQmuPlusInnerChi[iQQ];
			muPlusGlobalChi=QQmuPlusGlobalChi[iQQ];
			muPlusNMuonhits=QQmuPlusNMuonhits[iQQ];
			muPlusDB=QQmuPlusDB[iQQ];

			muMinusDxy=QQmuMinusDxy[iQQ];
			muMinusDz=QQmuMinusDz[iQQ];
			muMinusNhits=QQmuMinusNhits[iQQ];
			muMinusNPixelhits=QQmuMinusNPixelhits[iQQ];
			muMinusNPxlLayers=QQmuMinusNPxlLayers[iQQ];
			muMinusInnerChi=QQmuMinusInnerChi[iQQ];
			muMinusGlobalChi=QQmuMinusGlobalChi[iQQ];
			muMinusNMuonhits=QQmuMinusNMuonhits[iQQ];
			muMinusDB=QQmuMinusDB[iQQ];

			// std::cout << "muPlusDxy = " << muPlusDxy << std::endl;
			// std::cout << "muPlusDz = " << muPlusDz << std::endl;
			// std::cout << "muPlusNhits = " << muPlusNhits << std::endl;
			// std::cout << "muPlusNPxlLayers = " << muPlusNPxlLayers << std::endl;
			// std::cout << "muPlusInnerChi = " << muPlusInnerChi << std::endl;
			// std::cout << "mu1_pt = " << mu1_pt << std::endl;

			// std::cout << "muMinusDxy = " << muMinusDxy << std::endl;
			// std::cout << "muMinusDz = " << muMinusDz << std::endl;
			// std::cout << "muMinusNhits = " << muMinusNhits << std::endl;
			// std::cout << "muMinusNPxlLayers = " << muMinusNPxlLayers << std::endl;
			// std::cout << "muMinusInnerChi = " << muMinusInnerChi << std::endl;
			// std::cout << "mu2_pt = " << mu2_pt << std::endl;

			if( !(fabs(muPlusDxy)<3 && fabs(muPlusDz)<15
			      && muPlusNhits>10
			      // && muPlusNPxlLayers>0
			      && muPlusNPixelhits>0 // replaced the above, should be checked
			      && muPlusInnerChi<4
			      // && muPlusGlobalChi<20
			      // && mu1_pt > 4
			      && mu1_GM
			      && fabs(muMinusDxy)<3 && fabs(muMinusDz)<15 
			      && muMinusNhits>10
			      // && muMinusNPxlLayers>0
			      && muMinusNPixelhits>0 // replaced the above, should be checked
			      && muMinusInnerChi<4
			      // && muMinusGlobalChi<20 
			      // && mu2_pt > 4
			      && mu2_GM
			      && vProb>0.01) ) continue;

			// if( !(fabs(muPlusDxy)<3 && fabs(muPlusDz)<15 
			//       && muPlusNhits>10 && muPlusNPxlLayers>0 && muPlusNMuonhits>12 
			//       && muPlusInnerChi<2.1 && muPlusGlobalChi<2.7
			//       && fabs(muMinusDxy)<3 && fabs(muMinusDz)<15 
			//       && muMinusNhits>10 && muMinusNPxlLayers>0 && muMinusNMuonhits>12 
			//       && muMinusInnerChi<2.1 && muMinusGlobalChi<2.7 
			//       && vProb>0.19) ) continue;

			// if (HLTriggers&1==1 && Reco_QQ_trig[iQQ]&1==1
			if( mu1_HLT_L1DoubleMuOpen && mu2_HLT_L1DoubleMuOpen )
			{
				if (QQsign==0) {
					UpsilonTree->Fill();

					if (Reco_QQ->M()>7 && Reco_QQ->M()<14) {
						h_QQ_mass->Fill(Reco_QQ->M());
					}
					if (Reco_QQ_mupl->Pt()>=4 && Reco_QQ_mumi->Pt()>=4 && Reco_QQ->M()>7 && Reco_QQ->M()<14) {
						h_QQ_mass_1->Fill(Reco_QQ->M());
					}
				}
				else if (Reco_QQ_mupl->Pt()>=4 && Reco_QQ_mumi->Pt()>=4 && Reco_QQ->M()>7 && Reco_QQ->M()<14) {
					h_QQ_mass_2->Fill(Reco_QQ->M());
				}

			}
			//nReco_QQ++;
		}
		//	UpsilonTree->Fill();
	}

	TCanvas *c1 = new TCanvas("c1","c1");
	h_QQ_mass->GetYaxis()->SetRangeUser(0,300);
	h_QQ_mass->SetMarkerColor(kRed);
	h_QQ_mass->SetMarkerStyle(22);
	h_QQ_mass->SetXTitle("#mu#mu invariant mass [GeV/c^{2}]");
	h_QQ_mass->SetYTitle("Events/(0.1 GeV/c^{2})");
	h_QQ_mass->Draw("PE");

	h_QQ_mass_1->SetMarkerColor(kBlue);
	h_QQ_mass_1->SetMarkerStyle(20);
	h_QQ_mass_1->Draw("samePE");

	h_QQ_mass_2->SetMarkerColor(kRed);
	h_QQ_mass_2->SetMarkerStyle(24);
	h_QQ_mass_2->Draw("samePE");

	TLegend *legend = new TLegend(.45,.75,.99,.99);
	legend->AddEntry(h_QQ_mass,"#mu^{+}#mu^{-} pairs, no #mu p_{T} cut","P");
	legend->AddEntry(h_QQ_mass_1,"#mu^{+}#mu^{-} pairs, #mu p_{T} > 4 GeV/c","P");
	legend->AddEntry(h_QQ_mass_2,"same sign #mu pairs, #mu p_{T} > 4 GeV/c","P");
	legend->Draw();

	h_QQ_mass->Write();
	h_QQ_mass_1->Write();
	h_QQ_mass_2->Write();
	c1->SaveAs("paperPlots/ups_invmass_muptcut_hi_nocuts.pdf");
	c1->SaveAs("paperPlots/ups_invmass_muptcut_hi_nocuts.gif");
	c1->Write();
	f1->Write();

}
