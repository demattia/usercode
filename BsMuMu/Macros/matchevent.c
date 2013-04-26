void makeplot(TH1F *h1, TH1F *h2, TString name, TCanvas *c, int n) {
	c->cd(n);
	h1->SetLineColor(kRed);
	h1->GetXaxis()->SetTitle(name);
	h2->SetLineColor(kBlue);
	h2->SetLineStyle(kDashed);
	h1->Draw();
	h2->Draw("same");
}

void matchevent(){
	TFile *fxchk = new TFile("Barrel_preselection.root");
	//TFile *fmain = new TFile("data_main_endcaps.root");
	TFile *fmain = new TFile("MainAnalysis/data_afterCuts_0.root");

	TTree *Txchk = (TTree*)fxchk->Get("probe_tree");
	TTree *Tmain = (TTree*)fmain->Get("events");

	Float_t mu1_pt, mu2_pt, l3d, l3dsig, mass, dca, delta3d, cosAlpha3D, minDca, isolation;
	Float_t ntrk20;
	UInt_t event;
	Double_t m1pt, m2pt, fl3d, fls3d, m, maxdoca, pvip, cosa, docatrk, iso;
	Int_t closetrk;
	Long64_t evt;

	Txchk->SetBranchAddress("mu1_pt",&mu1_pt);
	Txchk->SetBranchAddress("mu2_pt",&mu2_pt);
	Txchk->SetBranchAddress("l3d",&l3d);
	Txchk->SetBranchAddress("l3dsig",&l3dsig);
	Txchk->SetBranchAddress("mass",&mass);
	Txchk->SetBranchAddress("ntrk20",&ntrk20);
	Txchk->SetBranchAddress("dca",&dca);
	Txchk->SetBranchAddress("delta3d",&delta3d);
	Txchk->SetBranchAddress("cosAlpha3D",&cosAlpha3D);
	Txchk->SetBranchAddress("minDca",&minDca);
	Txchk->SetBranchAddress("isolation",&isolation);
	Txchk->SetBranchAddress("event",&event);

	Tmain->SetBranchAddress("m1pt",&m1pt);
	Tmain->SetBranchAddress("m2pt",&m2pt);
	Tmain->SetBranchAddress("fl3d",&fl3d);
	Tmain->SetBranchAddress("fls3d",&fls3d);
	Tmain->SetBranchAddress("m",&m);
	Tmain->SetBranchAddress("closetrk",&closetrk);
	Tmain->SetBranchAddress("maxdoca",&maxdoca);
	Tmain->SetBranchAddress("pvip",&pvip);
	Tmain->SetBranchAddress("cosa",&cosa);
	Tmain->SetBranchAddress("docatrk",&docatrk);
	Tmain->SetBranchAddress("iso",&iso);
	Tmain->SetBranchAddress("evt",&evt);

	TH1F *hmu1_pt = new TH1F("mu1_pt","mu1_pt",40,0,50);
	TH1F *hmu2_pt = new TH1F("mu2_pt","mu2_pt",40,0,50);
	TH1F *hl3d = new TH1F("l3d","l3d",40,0,5);
	TH1F *hl3dsig = new TH1F("l3dsig","l3dsig",40,0,10);
	TH1F *hmass = new TH1F("mass","mass",40,4.9,5.9);
	TH1F *hntrk20 = new TH1F("ntrk20","ntrk20",40,0,20);
	TH1F *hdca = new TH1F("dca","dca",40,0,0.1);
	TH1F *hdelta3d = new TH1F("delta3d","delta3d",40,0,0.5);
	TH1F *hcosAlpha3D = new TH1F("cosAlpha3D","cosAlpha3D",40,0,1);
	TH1F *hminDca = new TH1F("minDca","minDca",40,0,0.2);
	TH1F *hisolation = new TH1F("isolation","isolation",40,0,1);

	TH1F *hm1pt = new TH1F("m1pt","m1pt",40,0,50);
	TH1F *hm2pt = new TH1F("m2pt","m2pt",40,0,50);
	TH1F *hfl3d = new TH1F("fl3d","fl3d",40,0,5);
	TH1F *hfls3d = new TH1F("fls3d","fls3d",40,0,10);
	TH1F *hm = new TH1F("m","m",40,4.9,5.9);
	TH1F *hclosetrk = new TH1F("closetrk","closetrk",40,0,20);
	TH1F *hmaxdoca = new TH1F("maxdoca","maxdoca",40,0,0.1);
	TH1F *hpvip = new TH1F("pvip","pvip",40,0,0.5);
	TH1F *hcosa = new TH1F("cosa","cosa",40,0,1);
	TH1F *hdocatrk = new TH1F("docatrk","docatrk",40,0,0.2);
	TH1F *hiso = new TH1F("iso","iso",40,0,1);

	Int_t nentries_xchk = (Int_t)Txchk->GetEntries();
	Int_t nentries_main = (Int_t)Tmain->GetEntries();

	for (int i=0; i<100; i++) {
		Txchk->GetEntry(i);
		cout<<event<<" "<<mu1_pt<<" "<<mu2_pt<<" "<<l3d<<" "<<l3dsig<<" "<<mass<<" "<<dca<<" "<<delta3d<<" "<<cosAlpha3D<<" "<<minDca<<" "<<isolation<<" "<<ntrk20<<endl;
		hmu1_pt->Fill(mu1_pt);
		hmu2_pt->Fill(mu2_pt);
		hl3d->Fill(l3d);
		hl3dsig->Fill(l3dsig);
		hmass->Fill(mass);
		hntrk20->Fill(ntrk20);
		hdca->Fill(dca);
		hdelta3d->Fill(delta3d);
		hcosAlpha3D->Fill(cosAlpha3D);
		hminDca->Fill(minDca);
		hisolation->Fill(isolation);


		for (int j=0; j<nentries_main; j++) {
			Tmain->GetEntry(j);
			bool match = false; 
			if (evt==event) {
				//cout<<"found the matched event: "<<endl;
				cout<<evt<<" "<<m1pt<<" "<<m2pt<<" "<<fl3d<<" "<<fls3d<<" "<<m<<" "<<maxdoca<<" "<<pvip<<" "<<cosa<<" "<<docatrk<<" "<<iso<<" "<<closetrk<<endl<<endl;
				hm1pt->Fill(m1pt);
				hm2pt->Fill(m2pt);
				hfl3d->Fill(fl3d);
				hfls3d->Fill(fls3d);
				hm->Fill(m);
				hclosetrk->Fill(closetrk);
				hmaxdoca->Fill(maxdoca);
				hpvip->Fill(pvip);
				hcosa->Fill(cosa);
				hdocatrk->Fill(docatrk);
				hiso->Fill(iso);

				match = true;
				break;
			}
		}

		if (!match)
			cout<<"!!!!!***** No Matched Event in the Main Analysis Tree *****!!!!!"<<endl<<endl;
	}

	TCanvas *c = new TCanvas("c","c",1600,1200);
	c->Divide(4,3);

	makeplot(hmu1_pt,hm1pt,"mu1_pt",c,1);
	makeplot(hmu2_pt,hm2pt,"mu2_pt",c,2);
	makeplot(hl3d,hfl3d,"l3d",c,3);
	makeplot(hl3dsig,hfls3d,"l3dsig",c,4);
	makeplot(hmass,hm,"mass",c,5);
	makeplot(hntrk20,hclosetrk,"ntrk20",c,6);
	makeplot(hdca,hmaxdoca,"dca",c,7);
	makeplot(hdelta3d,hpvip,"delta3d",c,8);
	makeplot(hcosAlpha3D,hcosa,"cosAlpha3D",c,9);
	makeplot(hminDca,hdocatrk,"minDca",c,10);
	makeplot(hisolation,hiso,"isolation",c,11);

	c->SaveAs("plots/data_barrel.gif");
}
