#include "setdirs.h"

TString fileName(const TString & method, const TString & region, const TString & index);
void significanceCuts(const double & nsig, const double & nbkg, TString method="BDT",TString region="barrel", TString index="", const int subdir = 1);

void significance(const double & nsig, const double & nbkg, TString method="BDT",TString region="barrel", TString index="0", const int subdir = 1) {

  cout << "significance: processing " << method << " for " << region << "(sub sample:"<< index << ")" << endl;
  
  if(method.Contains("Cuts")) {
    significanceCuts(nsig,nbkg,method,region,index,1);
    return; 
  }

  TString fnameA = fileName(method, region, index);
  TFile* inputA = TFile::Open(fnameA);

  if(subdir) gDirectory->Cd(fnameA+":/Method_"+method+"/"+method);
  gDirectory->pwd();
  TH1F* tmva_s = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_S_high"));
  TH1F* tmva_b = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_B_high"));

  const int rb =1;//250;
  tmva_s->Rebin(rb);
  tmva_b->Rebin(rb);

  int nbins   = tmva_s -> GetNbinsX();
  double xmin = tmva_s -> GetBinLowEdge(1);
  double xmax = tmva_s -> GetBinLowEdge(nbins+1);

  double normS = tmva_s -> Integral(1,nbins);
  double normB = tmva_b -> Integral(1,nbins);
  tmva_s->Scale(nsig/normS);
  tmva_b->Scale(nbkg/normB);
  normS = tmva_s -> Integral(1,nbins);
  normB = tmva_b -> Integral(1,nbins);

  double sumS(0), sumB(0), signi(0), signi1(0), signi2(0); 

  //printf("check bins: %d  %f %f %f %f \n", nbins, xmin, xmax, sumS, normS);

  TH1F* effS = new TH1F("effS","effS",nbins,xmin,xmax);
  TH1F* effB = new TH1F("effB","effB",nbins,xmin,xmax);
  TH1F* rejB = new TH1F("rejB","rejB",nbins,xmin,xmax);
  TH1F* sign = new TH1F("sign","sign",nbins,xmin,xmax);
  TH1F* sign1 = new TH1F("sign1","sign1",nbins,xmin,xmax);
  TH1F* sign2 = new TH1F("sign2","sign2",nbins,xmin,xmax);

  TH2F* roc  = new TH2F("roc","roc",nbins,0,1.02,nbins,0,1.02);

  for(int i=0; i<nbins+1; i++) {

    sumS  = tmva_s -> Integral(i,nbins);
    sumB  = tmva_b -> Integral(i,nbins);

    rejB->SetBinContent(i, 1.-sumB);
    effS->SetBinContent(i, normS?sumS/normS:0);
    effB->SetBinContent(i, normB?sumB/normB:0);
    rejB->SetBinContent(i, normB?1.-sumB/normB:1);

    roc->Fill(normS?sumS/normS:0,normB?1.-sumB/normB:1);

    signi = sumS/sqrt(sumS + sumB);
    sign->SetBinContent(i, signi);

    signi1 = sumB?sumS/sqrt(sumB):0;
    sign1->SetBinContent(i, signi1);

    signi2 = sumS/(sqrt(sumB)+0.5);
    sign2->SetBinContent(i, signi2);

    //printf("--> signif : %f  %f %f\n", signi, signi1, signi2 );
	   

  }


  int sig_max_bin = sign->GetMaximumBin(); 
  double sig_max_mva = effS->GetBinCenter(sig_max_bin);
  
  TF1 *pol = new TF1("fa1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)",sig_max_mva-0.1,sig_max_mva+0.1);
  pol->SetParameters(1,1,1,1);
  //pol->FixParameter(3,0); pol->SetRange(0.1,0.5);
  sign->Fit(pol, "rq");

  TF1 *pol1 = (TF1*) pol->Clone();
  sign1 -> Scale(sign->Integral(1,nbins)/sign1->Integral(1,nbins));
  sign1->Fit(pol1, "rq");
  sign1->Scale(pol->GetMaximum()/pol1->GetMaximum());
  sign1->Fit(pol1, "rq");

  TF1 *pol2 = (TF1*) pol->Clone();
  sign2 -> Scale(sign->Integral(1,nbins)/sign2->Integral(1,nbins));
  sign2->Fit(pol2, "rq");
  sign2->Scale(pol->GetMaximum()/pol2->GetMaximum());
  sign2->Fit(pol2, "rq");
  
  sign->Fit(pol, "rq");

  double sig_max      = sign->GetBinContent(sig_max_bin);
  double sig_max_effS = effS->GetBinContent(sig_max_bin);
  double sig_max_effB = effB->GetBinContent(sig_max_bin);
  //printf("maximum significance: %f bin:%d mva:%f %f %f \n", sig_max, sig_max_bin, sig_max_mva, sig_max_effS, sig_max_effB);

  double max_sig_pol_val  = pol ->GetMaximum();
  double max_sig_pol_mva  = pol ->GetMaximumX();
  double max_sig_pol1_mva = pol1->GetMaximumX();
  double max_sig_pol2_mva = pol2->GetMaximumX();
  //double max_sig_pol_bin  = sign->GetBinCenter(sign->FindBin(max_sig_pol_mva));
  double max_sig_pol_effS = effS->GetBinContent(sign->FindBin(max_sig_pol_mva));
  double max_sig_pol_effB = effB->GetBinContent(sign->FindBin(max_sig_pol_mva));

  printf("pol fit to sign: max %f, mva %f effb %f effs %f \n",
	 max_sig_pol_val, max_sig_pol_mva, max_sig_pol_effB, max_sig_pol_effS);
  printf("histogram sign:  max %f  mva %f effb %f effs %f \n", 
	 sig_max, sig_max_mva, sig_max_effB, sig_max_effS);

   char ss[100];
   sprintf(ss, "%s,  max. significance: %3.1f, %s>%5.4f, #epsilon_{S}=%4.3f, #epsilon_{B}=%5.4f",
	   region.Data(), max_sig_pol_val, method.Data(), pol->GetMaximumX(), max_sig_pol_effS, max_sig_pol_effB);
   printf("%s \n",ss);

   ofstream outputTxt(logsDir+"maxsignificance_"+method+"_"+region+".txt");
   outputTxt << "s/sqrt(s+b): " << max_sig_pol_mva << "\t s/sqrt(b): " << max_sig_pol1_mva << "\t s/sqrt(b)+0.5: "  << max_sig_pol2_mva << " (mva cut maximizing significance estimators)";
   outputTxt.close();

  //setTDRStyle(false);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(0); 

  TCanvas c;
  c.SetGrid();
  effS->SetTitle("");
  effS->GetYaxis()->SetTitle("#epsilon_{S}   |   1-#epsilon_{B}  ");
  effS->GetXaxis()->SetTitle(method+ " >     ");
  effS->SetLineColor(kBlue);
  effS->Draw("c");
  rejB->SetLineColor(kRed);
  rejB->Draw("same c");
  c.Update();  

  Float_t rightmax = 1.02*sign1->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;

  sign ->Scale(scale);
  sign1->Scale(scale);
  sign2->Scale(scale);

  sign ->SetLineColor(8);
  sign1->SetLineColor(38);
  sign2->SetLineColor(42);

  sign ->SetLineWidth(2);
  sign1->SetLineWidth(2);
  sign2->SetLineWidth(2);

  sign ->SetLineStyle(4);
  sign1->SetLineStyle(3);
  sign2->SetLineStyle(2);

  sign ->Draw("same c");
  sign1->Draw("same c");
  sign2->Draw("same c");

  pol ->SetLineColor(8);
  pol1->SetLineColor(38);
  pol2->SetLineColor(42);

  pol ->SetLineWidth(2);
  pol1->SetLineWidth(1);
  pol2->SetLineWidth(1);

  //pol1->SetLineStyle(3);
  //pol2->SetLineStyle(2);

  pol ->SetParameter(0, pol ->GetParameter(0)*scale);
  pol1->SetParameter(0, pol1->GetParameter(0)*scale);
  pol2->SetParameter(0, pol2->GetParameter(0)*scale);

  pol ->Draw("same");
  pol1->Draw("same");
  pol2->Draw("same");


  TLegend* leg = new TLegend(0.12,0.65,0.45,0.75);
  //leg->SetHeader("The Legend Title");
  char ss0[50], ss1[50],ss2[50];
  sprintf(ss0, "S/#sqrt{S+B},     max.@ %s>%5.4f",     method.Data(), max_sig_pol_mva);
  sprintf(ss1, "S/#sqrt{B},          max.@ %s>%5.4f",  method.Data(), max_sig_pol1_mva);
  sprintf(ss2, "S/(#sqrt{B}+0.5), max.@ %s>%5.4f",     method.Data(), max_sig_pol2_mva);
  leg->AddEntry(sign ,ss0,"l");
  leg->AddEntry(sign1,ss1,"l");
  leg->AddEntry(sign2,ss2,"l");
  leg->Draw();


  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
  	  gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor (8);
   axis->SetLabelColor(8);
   axis->SetTextColor (8);
   axis->SetTitle("Significance");
   axis->Draw();

   TPaveText *tp = new TPaveText(0.2,0.9,0.9,0.95,"brNDC");
   tp->SetBorderSize(0);
   tp->AddText(ss);
   tp->Draw("same");
   
   TLine* tl = new TLine(sig_max_mva,0, max_sig_pol_mva, gPad->GetUymax());
   tl->SetLineStyle(4);
   tl->Draw("same");

   TString ext(".gif");


  TString name = method + "_" + region;
  if(index!="") name += "_"+index;

   c.SaveAs(TString(plotsDir+name+"_eff"+ext));

   TString lpath(figuresDir);
   if(method == "BDT")  lpath += "bdt/";
   if(method == "MLP")  lpath += "mlp/";
   if(method.Contains("Cuts")) lpath += "cnt/";
   c.SaveAs(TString(lpath+name+"_eff.pdf"));

   TCanvas c2;
   roc->SetTitle("");
   roc->GetXaxis()->SetTitle("#epsilon_{S}");
   roc->GetYaxis()->SetTitle("1-#epsilon_{B}");
   //roc->SetMarkerStyle(7);
   roc->Draw();
   //tp->Draw();
   TMarker* tm = new TMarker(sig_max_effS,1.-sig_max_effB, 28);
   tm->SetMarkerSize(1.5);
   tm->SetMarkerColor(2);
   tm->Draw("same");
   c2.SaveAs(TString(plotsDir+name+"_roc"+ext));
   c2.SaveAs(TString(lpath+name+"_roc.pdf"));

   TCanvas c3;
   roc->GetXaxis()->SetRangeUser(0.,1);
   roc->GetYaxis()->SetRangeUser(0.99,1.);
   roc->Draw();
   //tp->Draw("same");
   tm->Draw("same");
   c3.SaveAs(TString(plotsDir+name+"_roc_zoom"+ext));
   c3.SaveAs(TString(lpath+name+"_roc_zoom.pdf"));

   TGraph* mark = new TGraph(1);
   mark->SetPoint(1,sig_max_effS,1.-sig_max_effB);
   mark->SetName("marker");
   TString outfileName(plotsDir+"signif_"+name+".root");
   TFile* outputFile = TFile::Open(outfileName,"RECREATE");
   roc->Write();   
   mark->Write();   
   outputFile->Close();

   
   inputA->Close();
}

void significanceCuts(const double & nsig, const double & nbkg, TString method,TString region, TString index, const int subdir) {

  assert(method=="CutsSA");
  //TString method("CutsSA");

  TString fname = fileName(method, region, index);
  TFile* input = TFile::Open(fname);
  if(subdir) gDirectory->Cd(fname+":/Method_Cuts/"+method);
  //gDirectory->pwd();
  //gDirectory->ls();
  TH1F* effB = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_effB"));
  TH1F* rejB = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_rejBvsS"));

  int nbins   = effB -> GetNbinsX();
  double xmin = effB -> GetBinLowEdge(1);
  double xmax = effB -> GetBinLowEdge(nbins+1);

  const int sfac = 1;
  //nbins /= sfac;

  TH1F* sign = new TH1F("sign","sign",nbins/sfac,xmin,xmax);
  TH1F* sign1 = new TH1F("sign1","sign1",nbins/sfac,xmin,xmax);
  TH1F* sign2 = new TH1F("sign2","sign2",nbins/sfac,xmin,xmax);

  double effb(0), effs(0), signi(0),  signi1(0), signi2(0);

  for(int i=0; i<nbins/sfac; i++) {

    int ibin= i; 
    //int ibint = ibin*sfac + int(sfac/2) +1;
    effs = effB->GetBinCenter(ibin*sfac + int(sfac/2) +1);

    effb=0;
    for(int j=0; j<sfac; j++) 
      effb+=effB->GetBinContent(ibin*sfac+j);    
    effb=effb/sfac;

    signi = effs*nsig/sqrt(effs*nsig + effb*nbkg);
    if(effb==0) signi=0;
    sign->SetBinContent(ibin, signi);

    signi1 = effb?effs*nsig/sqrt(effb*nbkg):0;
    sign1->SetBinContent(i, signi1);

    signi2 = effb?effs*nsig/(sqrt(effb*nbkg) + 0.5):0;
    sign2->SetBinContent(i, signi2);
    //printf("== bin %d effs %f effb %f  signi %f signi1 %f signi2 %f \n", ibin, effs, effb, signi, signi1, signi2);
    //printf("== ibin %d ibinav %d effs %f  effb %f   sign %f\n",ibin, ibint, effs, effb, signi);
  }
  
  //double sig_max = sign->GetMaximum();
  double yoffset = 0;
  for(int i=0; i<nbins/sfac; i++) {
    sign->SetBinContent(i, sign->GetBinContent(i)+yoffset);
  }


  TF1 *pol = new TF1("fa1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)",0.0,0.8);
  pol->SetParameters(1,1,1,1);
  //pol->FixParameter(3,0); pol->SetRange(0.1,0.5);
  sign->Fit(pol, "rq");

  TF1 *pol1 = (TF1*) pol->Clone();
  sign1 -> Scale(sign->Integral(1,nbins)/sign1->Integral(1,nbins));
  sign1->Fit(pol1, "rq");
  sign1->Scale(pol->GetMaximum()/pol1->GetMaximum());
  sign1->Fit(pol1, "rq");

  TF1 *pol2 = (TF1*) pol->Clone();
  sign2 -> Scale(sign->Integral(1,nbins)/sign2->Integral(1,nbins));
  sign2->Fit(pol2, "rq");
  sign2->Scale(pol->GetMaximum()/pol2->GetMaximum());
  sign2->Fit(pol2, "rq");

  sign->Fit(pol, "rq");

  //int sig_max_bin = sign->GetMaximumBin(); 
  //double sig_max = sign->GetBinContent(sig_max_bin);
  //double sig_max_effS = effB->GetBinCenter(sig_max_bin);
  //double sig_max_effB = effB->GetBinContent(sig_max_bin);

  double max_sig_pol_val = pol->GetMaximum();
  double max_sig_pol_x   = pol->GetMaximumX();
  double max_sig_pol1_x  = pol1->GetMaximumX();
  double max_sig_pol2_x  = pol2->GetMaximumX();
  //double max_sig_pol_bin = sign->GetBinCenter(sign->FindBin(max_sig_pol_x));
  //double max_sig_pol_effS = effB->GetBinCenter(sign->FindBin(max_sig_pol_x));
  double max_sig_pol_effB = effB->GetBinContent(sign->FindBin(max_sig_pol_x));
  
  //printf("pol fit to sign: max %f  max-x %f binval: %f\n ", max_sig_pol_val, max_sig_pol_x, max_sig_pol_bin);
  
  char ss[100];
  sprintf(ss, "%s,  max. significance: %3.1f, %s, #epsilon_{S}=%4.3f, #epsilon_{B}=%5.4f",region.Data(),max_sig_pol_val,method.Data(), pol->GetMaximumX(), max_sig_pol_effB);
  printf("%s \n",ss);
  
  ofstream outputTxt(logsDir+"maxsignificance_"+method+"_"+region+".txt");
  outputTxt << "s/sqrt(s+b): " << max_sig_pol_x << "\t s/sqrt(b): " <<  max_sig_pol1_x << "\t s/sqrt(b)+0.5: "  <<  max_sig_pol2_x << " (mva cut maximizing significance estimators)"; 
  outputTxt.close();
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(0); 
  
  
  TCanvas c;
  c.SetGrid();

  sign->GetXaxis()->SetRangeUser(0.,0.76);

  //double ymarg = 0.9;
  //Float_t rightmax = sign->GetMaximum();
  //Float_t scale = gPad->GetUymax()/rightmax;
  //scale *= (1-sig_max*(1.-ymarg)/rightmax);
  //sign->Scale(scale);
  //sign->Scale(1-sig_max*(1.-ymarg)/rightmax);
  //sign->Sumw2();

  sign->GetYaxis()->SetTitle("significance");
  sign->GetXaxis()->SetTitle("#epsilon_{S}");
  //sign->GetXaxis()->SetTextColor(8);

  sign->SetTitle("");
  sign->SetLineWidth(1);
  sign->SetLineColor(8);
  sign->Draw("hist c");

  pol  -> SetRange(pol ->GetMaximumX()-0.2,pol ->GetMaximumX()+0.2);
  pol1 -> SetRange(pol1->GetMaximumX()-0.2,pol1->GetMaximumX()+0.2);
  pol2 -> SetRange(pol2->GetMaximumX()-0.2,pol2->GetMaximumX()+0.2);

  pol ->SetLineColor(kBlue);
  pol1->SetLineColor(38);
  pol2->SetLineColor(42);

  pol ->SetLineWidth(2);
  pol1->SetLineWidth(2);
  pol2->SetLineWidth(2);

  //pol1->SetLineStyle(3);
  //pol2->SetLineStyle(2);

  pol ->Draw("same");
  pol1->Draw("same");
  pol2->Draw("same");

  TLegend* leg = new TLegend(0.3,0.2,0.7,0.3);
  char ss0[50], ss1[50],ss2[50];
  sprintf(ss0, "S/#sqrt{S+B},     max.@ #epsilon_{S}>%5.3f"   ,pol  ->GetMaximumX());
  sprintf(ss1, "S/#sqrt{B},          max.@ #epsilon_{S}>%5.3f",pol1 ->GetMaximumX());
  sprintf(ss2, "S/(#sqrt{B}+0.5), max.@ #epsilon_{S}>%5.3f"   ,pol2 ->GetMaximumX());
  leg->AddEntry(pol ,ss0,"l");
  leg->AddEntry(pol1,ss1,"l");
  leg->AddEntry(pol2,ss2,"l");
  leg->Draw();
  
  c.Update();  

  rejB->Scale(gPad->GetUymax()/1.1);
  rejB->SetTitle("");
  rejB->SetLineColor(kRed);
  rejB->Draw("c same");


  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),0, 1.1,510,"+L");
  axis->SetLineColor (kRed);
  //axis->SetLabelColor(8);
  //axis->SetTextColor (kRed);
  axis->SetTitle("1 - #epsilon_{B}");
  axis->Draw("same");

   TPaveText *tp = new TPaveText(0.2,0.9,0.9,0.95,"brNDC");
   tp->SetBorderSize(0);
   tp->AddText(ss);
   tp->Draw("same");

  TLine* tl = new TLine(pol->GetMaximumX(),0,pol->GetMaximumX(), gPad->GetUymax());
  tl->SetLineStyle(4);
  tl->Draw("same");
  
  TString ext(".gif");
  c.SaveAs(TString("plots/Cuts_"+region+"_eff"+ext));
  TString name = method + "_" + region;
  if(index!="") name += "_"+index;
  TString lpath(figuresDir+"cnt/");
  c.SaveAs(TString(lpath+name+"_eff.pdf"));

  input->Close();
}


TString fileName(const TString & method, const TString & region, const TString & index) {
  TString fnameA = "TMVA_" + region;
  if(index=="merged") 
    fnameA += "_"+method;
  if(index!="") 
    fnameA += "_"+index;
  fnameA += ".root";
  return rootDir+fnameA;
}
