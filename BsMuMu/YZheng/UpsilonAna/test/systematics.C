#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TPaveLabel.h"

using namespace std;

void systematics() {

  const int nsys = 1;

  ///open files
  TString fname[nsys+1]={"xsection.root", "xsection_syst1.root"};
  TFile* file[nsys+1];
  for(int i=0; i<nsys+1; i++) {
    file[i] = TFile::Open(fname[i],"READ");
  }
  
  /// nominal results

  gDirectory->Cd(fname[0]+":/plots");
  TGraphAsymmErrors *tot0 = (TGraphAsymmErrors*)gROOT->FindObject("xsection_total");
  TGraphAsymmErrors *dif0 = (TGraphAsymmErrors*)gROOT->FindObject("xsection_pt");
  //tot0->Draw("ap");
  //dif0->Draw("ap");

  // get pt bins
  int nn = dif0->GetN();
  const int nptb = nn;

  double *pt, *pt_el, *pt_eh;
  pt = dif0->GetX();
  pt_el = dif0->GetEXlow();
  pt_eh = dif0->GetEXhigh();

  double systematics[nptb+1][nsys+2];
  for(int i=0;i<nptb+1;i++) for(int j=0; j<nsys+1; j++) systematics[i][j]=0.;

  systematics[0][0] = (tot0->GetY())    [0]; // total|central
  systematics[0][1] = (tot0->GetEYlow())[0]; // total|statistical

  for(int i=0; i<nptb; i++) {
    systematics[i+1][0] = (dif0->GetY())    [i]; // differential|central
    systematics[i+1][1] = (dif0->GetEYlow())[i]; // differential|statistical
  }

  printf("\nNOMINAL MEASURED CROSS SECTION\n");
  printf("Total value: %5.2f +- %5.2f (stat) nb-1\n", 
	 systematics[0][0],systematics[0][1]);
  printf("Differential:\n");
  for(int i=0; i<nptb; i++) {
    printf("\t(pt=%4.1f):  %5.2f +- %5.2f (stat) nb-1\n", 
	   pt[i],systematics[i+1][0], systematics[i+1][1]);
  }


  printf("\nVARIED RESULTS\n\n");

  TGraphAsymmErrors *tot1, *dif1;

  for (int isys=0; isys<nsys; isys++) {

    gDirectory->Cd(fname[isys+1]+":/plots");
    tot1 = (TGraphAsymmErrors*)gROOT->FindObject("xsection_total");
    dif1 = (TGraphAsymmErrors*)gROOT->FindObject("xsection_pt");
    
    tot1->Draw("ap");
    dif1->Draw("ap");
    
    double fake = 1.5;
    systematics[0][isys+2] 
      = fabs( fake*(tot0->GetY())[0] - systematics[0][0]); // total|systematics-isys  
    
    for(int i=0; i<nptb; i++) {
      systematics[i+1][isys+2] 
	= fabs( fake*(dif0->GetY())[i] - systematics[i+1][0]); // differential|central
    }
    
    printf("  SYSTEMATICS %d  (%s)\n", isys, fname[isys+1].Data());
    printf("\tTotal value: %5.2f\n", systematics[0][isys+2]);
    printf("\tDifferential:\n"); 
    for(int i=0; i<nptb; i++) {
      printf("\t\t(pt=%4.1f):  val:%5.2f\n",
	     pt[i],systematics[i+1][isys+2]);
    }
  }
  
  /// systematics table
  for (int ibin=0; ibin<nptb+1; ibin++) {//row

    if(ibin==0) {
      printf("\n\npt-bin\t\tcentral\tstat.");
    
      for (int isys=0; isys<nsys; isys++) {
	printf("\tsyst.%d",isys+1);
      }
      printf("\tsyst.all\tstat+syst\t");
      cout << endl << flush;
    }
    
    if(ibin==0)
      printf("pt-integrated");
    else
      printf("[%4.1f-%4.1f]",pt[ibin-1]-pt_el[ibin-1],pt[ibin-1]+pt_eh[ibin-1]);

    double binw(0),tot_sys(0.), tot_err(0);;
    for (int isys=0; isys<nsys+2; isys++) {
      binw = (ibin==0)?1.:(pt_eh[ibin-1]+pt_el[ibin-1]);
      double val = systematics[ibin][isys]*binw;
      printf("\t%5.2f",val);
      if(isys>1) tot_sys += sqrt(pow(tot_sys,2)+pow(val,2));
      if(isys>0) tot_err += sqrt(pow(tot_err,2)+pow(val,2));
    }
    binw*=0.01;
    printf("\t%5.2f (%2.0d%%)\t%5.2f (%2.0d%%)",
	   tot_sys, tot_sys/systematics[ibin][0]/binw ,tot_err,tot_err/systematics[ibin][0]/binw);
    cout << endl << flush;


    //if(ibin==nptb) {
    //  printf("[%4.1f-%4.1f]",pt[0]-pt_el[0],pt[nptb-1]+pt_eh[nptb-1]);
    //  cout << endl << flush;
    //}

  }
}
