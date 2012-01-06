#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include <fstream>
using namespace RooFit;

const char *filenames[4] = {"limits/masses_data_muon.txt", "limits/masses_backgroundMC_muon.txt",
 			    "limits/masses_data_electron.txt", "limits/masses_backgroundMC_electron.txt"};
const char *outnames[4] = {"data.root", "bkgdata.root", "data.root", "bkgdata.root"};

void makeDataFile(int which, float minmass = 20.0) {
  RooRealVar mass("mass", "mass", 0);
  RooRealVar weight("weight", "weight", 0);

  RooDataSet d("modelData", "modelData", RooArgSet(mass,weight), WeightVar(weight));

  std::cout << "Reading in" << std::endl;

  std::ifstream infile(filenames[which]);

  double m,w;
  while ( true ) {

    infile >> m >> w;
    
    if ( infile.eof() ) break;

    if (m >= minmass) {
      mass = m;
      weight = w;
      d.add(RooArgSet(mass,weight));
    }
  }


  d.Print("v");
  cout << endl;

  // S a v i n g   a n d   l o a d i n g   f r o m   f i l e 
  // -------------------------------------------------------

  // Datasets can be persisted with ROOT I/O
  cout << endl << ">> Persisting d via ROOT I/O" << endl;
  TFile f(outnames[which], "RECREATE");
  d.Write();
  f.ls();

}
