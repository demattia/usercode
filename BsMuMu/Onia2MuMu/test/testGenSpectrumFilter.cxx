
void makePlots(TString name, bool normalize, const char *dir1, const char *dir2=0, const char *dir3 = 0, const char *dir4 = 0) {
    THStack *stk = new THStack("st","st");

    TH1 *h1 = (TH1*) gFile->Get(dir1+("/"+name));
    if (h1 == 0) { std::cerr << "Failed to find histogram 1:" << dir1+("/"+name) << std::endl; return; }
    if (normalize) { h1->Sumw2(); h1->Scale(1.0/h1->Integral()); }
    h1->SetLineWidth(2); h1->SetLineColor(kBlue);
    stk->Add(h1);

    if (dir2) {
        TH1 *h2 = (TH1*) gFile->Get(dir2+("/"+name));
        if (h2 == 0) { std::cerr << "Failed to find histogram 2:" << dir2+("/"+name) << std::endl; return; }
        if (normalize) { h2->Sumw2(); h2->Scale(1.0/h2->Integral()); }
        h2->SetLineWidth(2); h2->SetLineColor(kRed);
        stk->Add(h2);
    }

    if (dir3) {
        TH1 *h3 = (TH1*) gFile->Get(dir3+("/"+name));
        if (h3 == 0) { std::cerr << "Failed to find histogram 3:" << dir3+("/"+name) << std::endl; return; }
        if (normalize) { h3->Sumw2(); h3->Scale(1.0/h3->Integral()); }
        h3->SetLineWidth(2); h3->SetLineColor(209);
        stk->Add(h3);
    }

    if (dir4) {
        TH1 *h4 = (TH1*) gFile->Get(dir4+("/"+name));
        if (h4 == 0) { std::cerr << "Failed to find histogram 4:" << dir4+("/"+name) << std::endl; return; }
        if (normalize) { h4->Sumw2(); h4->Scale(1.0/h4->Integral()); }
        h4->SetLineWidth(2); h4->SetLineColor(222);
        stk->Add(h4);
    }

    stk->Draw("HIST NOSTACK");
    c1->Print("plots/"+TString(dir1)+"_"+name+".png");
}
void testGenSpectrumFilter() {
    makePlots("pt", 1,"plotJPsi", "lowPtPlots","highPtPlots");
    makePlots("eta",1,"plotJPsi", "lowEtaPlots","highEtaPlots");
    makePlots("y",  1,"plotJPsi", "lowEtaPlots","highEtaPlots");
    makePlots("pt", 1,"plotRecoJPsi", "lowPtRecoPlots","highPtRecoPlots");
    makePlots("eta",1,"plotRecoJPsi", "lowEtaRecoPlots","highEtaRecoPlots");
    makePlots("y",  1,"plotRecoJPsi", "lowEtaRecoPlots","highEtaRecoPlots");
}

