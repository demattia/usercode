#ifndef EVALCUTS_CC
#define EVALCUTS_CC

#include <iostream>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TDirectory.h"

using namespace std;

class EvalCuts {

  ofstream * multijetPixelEffFile;
  ofstream * allgoodcuts;

 public:

  EvalCuts() {
    // File to store the output
    multijetPixelEffFile = new ofstream( "EvalCuts.txt" );
    allgoodcuts = new ofstream( "EvalAllGoodCuts.txt" );
  }

  ~EvalCuts() {
    delete multijetPixelEffFile;
    delete allgoodcuts;
  }

  void eval( ) {
    // Do not display anything on screen
    gROOT->SetBatch(kTRUE);

    // Rate from the standard trigger (Hz)
    float standardRateMulti = 1576.30;
    float standardRateMEtJet = 4000.;

    TFile* file1 = new TFile("QCD.root");
    TFile* file2 = new TFile("Ana_TTH_120_tk3.root");

    TObject *obj1;
    TObject *obj2;
    TKey *key1;
    TKey *key2;
    TIter next1( file1->GetListOfKeys());
    TIter next2( file2->GetListOfKeys());
    while ((key1 = (TKey *) next1())) {
      key2 = (TKey *) next2();
      obj1 = file1->Get(key1->GetName()); // copy object one to memory
      obj2 = file2->Get(key2->GetName()); // copy object two to memory

      // Select histograms
      string className = key1->GetClassName();
      string TH1FName = "TH1F";
      bool comp = (className == TH1FName);

      if ( comp ) {

        TH1F * myhist1 = dynamic_cast<TH1F*>(obj1);
        TH1F * myhist2 = dynamic_cast<TH1F*>(obj2);

        string histoName( myhist1->GetName() );

        // Evaluate only on efficiency histograms
        string eff("Eff");

        //        cout << "histoName.find(\"Eff\") = " << histoName.find("Eff") << endl;

        if ( histoName.find( eff ) != string::npos ) {

          // New histogram
          TString division("_div");
          // If they do not have the same number of bins quit
          if ( myhist1->GetSize() != myhist2->GetSize() ) exit(1);

          // Take a map<float, int>, key is the efficiency and value is the bin index
          // the last element in the map will automatically be the highest efficiency.
          map<float, pair<int, float> > effMap;

          *allgoodcuts << histoName << endl;
          *allgoodcuts << "----------------" << endl;
          for (Int_t i=0; i<myhist1->GetSize(); ++i ) {

            // Use the rate condition on qcd to select corresponding bins for ttH
            float standardRate = 0.;
            if ( histoName.find( "Multi" ) != string::npos ) standardRate = standardRateMulti;
            else standardRate = standardRateMEtJet;

            if ( (*myhist1)[i] < standardRate ) {
              effMap.insert( make_pair( (*myhist2)[i], make_pair(i, (*myhist1)[i]) ) );
              *allgoodcuts << "eff = " << (*myhist2)[i] << " for rate = " << (*myhist1)[i] << " in bin = " << i << endl;
            }
          }
          *allgoodcuts << endl;

          map<float, pair<int, float> >::const_reverse_iterator rit = effMap.rbegin();

          *multijetPixelEffFile << histoName << endl;
          *multijetPixelEffFile << "-------------------" << endl;
          *multijetPixelEffFile << "highest efficiency = " << rit->first  << endl;
          *multijetPixelEffFile << "corresponding to rate = " << rit->second.second  << endl;
          *multijetPixelEffFile << "for i = " << rit->second.first << endl;

          if ( histoName.find( "Multi" ) != string::npos ) {
            // Determine the bin
            int multijetCount = 0;
            for ( int Et1=205; Et1<255; Et1+=5 ) {
              for ( int Et2=155; Et2<205; Et2+=5 ) {
                for ( int Et3=55; Et3<105; Et3+=5 ) {
                  for ( int Et4=35; Et4<85; Et4+=5 ) {
                    // Increment before, since the bins start from 1
                    ++multijetCount;
                    if ( multijetCount == rit->second.first ) {
                      *multijetPixelEffFile << "Et1 = " << Et1 << endl;
                      *multijetPixelEffFile << "Et2 = " << Et2 << endl;
                      *multijetPixelEffFile << "Et3 = " << Et3 << endl;
                      *multijetPixelEffFile << "Et4 = " << Et4 << endl;
                    }
                  }
                }
              }
            }
          }
          else {
            *multijetPixelEffFile << "MEt = " << endl;
            *multijetPixelEffFile << "Jet Et = " << endl;
          }
          *multijetPixelEffFile << endl;

        } // end if Eff
      } // end if TH1F
    } // end loop on keys
    multijetPixelEffFile->close();
    allgoodcuts->close();
  }
};

#endif
