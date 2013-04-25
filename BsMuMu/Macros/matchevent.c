{
	TFile *fxchk = new TFile("BsMC12_barrel_preselection.root");
	TFile *fmain = new TFile("MC_main_barrel.root");

	TTree *Txchk = (TTree*)fxchk->Get("probe_tree");
	TTree *Tmain = (TTree*)fmain->Get("events");

	Float_t mu1_pt, mu2_pt, l3dsig;
	UInt_t event;
	Double_t m1pt, m2pt, fls3d;
	Long64_t evt;

	Txchk->SetBranchAddress("mu1_pt",&mu1_pt);
	Txchk->SetBranchAddress("mu2_pt",&mu2_pt);
	Txchk->SetBranchAddress("l3dsig",&l3dsig);
	Txchk->SetBranchAddress("event",&event);
	Tmain->SetBranchAddress("m1pt",&m1pt);
	Tmain->SetBranchAddress("m2pt",&m2pt);
	Tmain->SetBranchAddress("fls3d",&fls3d);
	Tmain->SetBranchAddress("evt",&evt);

	TH1F *hmu1_pt = new TH1F("hmu1_pt","mu1_pt",50,0,50);
	TH1F *hl3dsig = new TH1F("hl3dsig","l3dsig",50,0,10);
	TH1F *hm1pt = new TH1F("hm1pt","m1pt",50,0,50);
	TH1F *hfls3d = new TH1F("hfls3d","fls3d",50,0,10);

	Int_t nentries_xchk = (Int_t)Txchk->GetEntries();
	Int_t nentries_main = (Int_t)Tmain->GetEntries();

	for (int i=0; i<100; i++) {
		Txchk->GetEntry(i);
		cout<<event<<" "<<mu1_pt<<" "<<mu2_pt<<" "<<l3dsig<<" "<<endl;
		hmu1_pt->Fill(mu1_pt);
		hl3dsig->Fill(l3dsig);

		for (int j=0; j<nentries_main; j++) {
			Tmain->GetEntry(j);
			bool match = false; 
			if (evt==event) {
				cout<<"found the matched event: "<<endl;
				cout<<evt<<" "<<m1pt<<" "<<m2pt<<" "<<fls3d<<endl<<endl;
				hm1pt->Fill(m1pt);
				hfls3d->Fill(fls3d);
				match = true;
				break;
			}
		}

		if (!match)
			cout<<"!!!!!***** No Matched Event in the Main Analysis Tree *****!!!!!"<<endl<<endl;
	}

	TCanvas *cmu1_pt = new TCanvas("cmu1_pt ","cmu1_pt");
	hmu1_pt->SetLineColor(kRed);
	hm1pt->SetLineColor(kBlue);
	hm1pt->SetLineStyle(kDashed);
	hmu1_pt->Draw();
	hm1pt->Draw("same");

    TCanvas *cl3dsig = new TCanvas("cl3dsig ","cl3dsig");
    hl3dsig->SetLineColor(kRed);
    hfls3d->SetLineColor(kBlue);
	hfls3d->SetLineStyle(kDashed);
    hl3dsig->Draw();
    hfls3d->Draw("same");
}
