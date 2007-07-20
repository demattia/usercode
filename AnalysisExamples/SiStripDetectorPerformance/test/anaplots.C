void anaplots(){

  //TIB: subid==3
  //TOB: subid==5
  //TEC: subid==  



 //cluster positions

 THStack Position("Position","Hits on tracks Position");
 TH1F* pos1 = new TH1F("pos1","pos1",400, 0, 800);
 TH1F* pos2 = new TH1F("pos2","pos2",400, 0, 800);
 TH1F* pos1tk = new TH1F("pos1tk","pos1tk",400, 0, 800);
 TH1F* pos2tk = new TH1F("pos2tk","pos2tk",400, 0, 800);

 nt->Draw("Clu_1strip >> pos1");
 pos1->SetFillColor(kBlue);
 Position.Add(pos1);
 nt->Draw("Clu_bar >> pos1");
 pos2->SetFillColor(kBlue);
 Position.Add(pos2);


 Position.Draw();


 THStack allchg("allchg", "Total Charge");
 TH1F* chrg1 = new TH1F("chrg1", "chrg1",  60, 0., 300.);
 TH1F* chrg2 = new TH1F("chrg2", "chrg2",  60, 0., 300.);
 TH1F* chrg3 = new TH1F("chrg3", "chrg3",  60, 0., 300.);
 TH1F* chrg4 = new TH1F("chrg4", "chrg4",  60, 0., 300.);
 nt->Draw("Clu_ch >> chrg1","Clu_size==1");
 chrg1->SetFillColor(kBlue);
 allchg.Add(chrg1);
 nt->Draw("Clu_ch >> chrg2","Clu_size==2");
 chrg2->SetFillColor(kRed);
 allchg.Add(chrg2);
 nt->Draw("Clu_ch >> chrg3","Clu_size==3");
 chrg3->SetFillColor(kGreen);
 allchg.Add(chrg3);
 nt->Draw("Clu_ch >> chrg4","Clu_size>=4");
 chrg4->SetFillColor(kYellow);
 allchg.Add(chrg4);


 //now plots

//track parameters
TCanvas* c1 = new TCanvas("c1");
                                                                                                          
 c1->Divide(2,2);
 c1->cd(1);
 nt->Draw("p_tk","p_tk<100");
 c1->cd(2);
 nt->Draw("pt_tk","p_tk<100");
 c1->cd(3);
 nt->Draw("phi_tk","p_tk<100");
 c1->cd(4);
 nt->Draw("eta_tk","p_tk<100");


 TCanvas* c2 = new TCanvas("c2");
 c2->Divide(2,2);
 c2->cd(1);
 Position.Draw();
 c2->cd(2);
 nt->Draw("Clu_size");
 c2->cd(3);
 nt->Draw("Clu_ang");
 c2->cd(4);
 allchg.Draw(); 
}
