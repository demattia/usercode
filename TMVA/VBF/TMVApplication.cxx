/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVApplication                                                     *
 *                                                                                *
 * This exectutable provides a simple example on how to use the trained           *
 * classifiers within a C++ analysis module                                       *
 *                                                                                *
 * ------------------------------------------------------------------------------ *
 * see also the alternative (slightly faster) way to retrieve the MVA values in   *
 * examples/TMVApplicationAlternative.cxx                                         *
 * ------------------------------------------------------------------------------ *
 **********************************************************************************/

#include <iostream> // Stream declarations
#include <algorithm>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "EventsReader.h"

using namespace std;

// ---------------------------------------------------------------
// choose MVA methods to be applied
Bool_t Use_Cuts            = 0;
Bool_t Use_CutsD           = 0;
Bool_t Use_CutsGA          = 0;
Bool_t Use_Likelihood      = 0;
Bool_t Use_LikelihoodD     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
Bool_t Use_LikelihoodPCA   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
Bool_t Use_PDERS           = 0;
Bool_t Use_PDERSD          = 0;
Bool_t Use_PDERSPCA        = 0;
Bool_t Use_KNN             = 0;
Bool_t Use_HMatrix         = 0;
Bool_t Use_Fisher          = 1;
Bool_t Use_FDA_GA          = 0;
Bool_t Use_FDA_MT          = 0;
Bool_t Use_MLP             = 0; // this is the recommended ANN
Bool_t Use_CFMlpANN        = 0; 
Bool_t Use_TMlpANN         = 0; 
Bool_t Use_SVM_Gauss       = 0;
Bool_t Use_SVM_Poly        = 0;
Bool_t Use_SVM_Lin         = 0;
Bool_t Use_BDT             = 1;
Bool_t Use_BDTD            = 0;
Bool_t Use_RuleFit         = 0;
// ---------------------------------------------------------------

int main( void ) 
{
  cout << endl;
  cout << "==> Start TMVApplication" << endl;
 
  //
  // create the Reader object
  //
  TMVA::Reader *reader = new TMVA::Reader("!Color");    

  int variablesNumber = 48;

  vector<TString> eventVariablesNamesVector;
  eventVariablesNamesVector.push_back("etq1"    );
  eventVariablesNamesVector.push_back("etq2"    );
  eventVariablesNamesVector.push_back("etaq1"   );
  eventVariablesNamesVector.push_back("etaq2"   );
  eventVariablesNamesVector.push_back("phiq1"   );
  eventVariablesNamesVector.push_back("phiq2"   );
  eventVariablesNamesVector.push_back("dphiqq"  );
  eventVariablesNamesVector.push_back("detaqq"  );
  eventVariablesNamesVector.push_back("ptminqq" );
  eventVariablesNamesVector.push_back("ptmaxqq" );
  eventVariablesNamesVector.push_back("etaminqq");
  eventVariablesNamesVector.push_back("etamaxqq");
  // Z jets variables
  eventVariablesNamesVector.push_back("eth1"    );
  eventVariablesNamesVector.push_back("eth2"    );
  eventVariablesNamesVector.push_back("etah1"   );
  eventVariablesNamesVector.push_back("etah2"   );
  eventVariablesNamesVector.push_back("phih1"   );
  eventVariablesNamesVector.push_back("phih2"   );
  eventVariablesNamesVector.push_back("dphihh"  );
  eventVariablesNamesVector.push_back("detahh"  );
  eventVariablesNamesVector.push_back("ptminhh" );
  eventVariablesNamesVector.push_back("ptmaxhh" );
  eventVariablesNamesVector.push_back("etaminhh");
  eventVariablesNamesVector.push_back("etamaxhh");
  // tag jets system variables
  eventVariablesNamesVector.push_back("ptqq" );
  eventVariablesNamesVector.push_back("mqq"  );
  eventVariablesNamesVector.push_back("etaqq");
  // Z jets system variables
  eventVariablesNamesVector.push_back("pthh" );
  eventVariablesNamesVector.push_back("mhh"  );
  eventVariablesNamesVector.push_back("etahh");
  // Z leptons system variables
  eventVariablesNamesVector.push_back("ptll" );
  eventVariablesNamesVector.push_back("mll"  );
  eventVariablesNamesVector.push_back("etall");
  // 2-system variables
  eventVariablesNamesVector.push_back("dphiTjetZjet");
  eventVariablesNamesVector.push_back("dphiTjetZlep");
  eventVariablesNamesVector.push_back("dphiminTZ"   );
  eventVariablesNamesVector.push_back("detaTjetZjet");
  eventVariablesNamesVector.push_back("detaTjetZlep");
  eventVariablesNamesVector.push_back("detaminTZ"   );
  eventVariablesNamesVector.push_back("dphiZjetZlep");
  eventVariablesNamesVector.push_back("detaZjetZlep");
  eventVariablesNamesVector.push_back("massTjetZjet");
  eventVariablesNamesVector.push_back("massTjetZlep");
  eventVariablesNamesVector.push_back("massZjetZlep");
  // 3-system variables
  eventVariablesNamesVector.push_back("massTZZ"     );
  eventVariablesNamesVector.push_back("etaTZZ"      );
  eventVariablesNamesVector.push_back("ptTZZ"       );
  eventVariablesNamesVector.push_back("eventNumber" );

  // Initialize the readers:
  // the last variable must be the event number.
  EventsReader signalReader( eventVariablesNamesVector, "../../tmvaSIGNAL200.root", "tree_SIGNAL200", variablesNumber );
  EventsReader backgroundReader( eventVariablesNamesVector, "../../tmvaCOMBINATORIAL200.root", "tree_COMBINATORIAL200", variablesNumber );





//   TMVA::Reader * reader_ = new TMVA::Reader("!Color");    
//   Float_t * variables_ = new Float_t[variablesNumber];
//   for( int iVar = 0; iVar < variablesNumber - 1; ++iVar ) {
//     reader_->AddVariable( eventVariablesNamesVector[iVar], &(variables_[iVar]) );
//   }
//   string dir    = "weights/";
//   string prefix = "TMVAnalysis";
//   reader_->BookMVA( "BDT method", dir + prefix + "_BDT.weights.txt" );

//   TTree * theTree_ = signalReader.tree();
//   vector<TString>::const_iterator varName = eventVariablesNamesVector.begin();
//   int iVar = 0;
//   for( ; varName != eventVariablesNamesVector.end(); ++varName, ++iVar ) {
//     theTree_->SetBranchAddress( *varName, &(variables_[iVar]) );
//   }
//   for (Long64_t ievt=0; ievt<theTree_->GetEntries();ievt++) {
//     theTree_->GetEntry(ievt);
// //     cout << "variable[46] inside = " << (signalReader.variables())[46] << endl;
// //     cout << "signalReader[46] inside = " << signalReader[46] << endl;
// //     cout << "variable[46] address = " << &(signalReader.variables())[46] << endl;
// //     cout << "variable[46] outside = " << variables_[46] << endl;
//     cout << "variable[46] before = " << variables_[46] << endl;
//     cout << "discriminant = " << reader_->EvaluateMVA( "BDT method" ) << endl;
//     cout << "variable[46] after = " << variables_[46] << endl;
//     cout << "discriminant = " << reader_->EvaluateMVA( "BDT method" ) << endl;
//     cout << "variable[46] after second call = " << variables_[46] << endl;
// //     cout << "discriminant = " << signalReader.discriminant() << endl;
//   }



  // Efficiency counter
  float efficiency = 0;


  // Loop on the tree
  // ----------------

  TTree * backgroundTree = backgroundReader.tree();
  TTree * signalTree = signalReader.tree();
  cout << "--- Processing: " << backgroundTree->GetEntries() << " events" << endl;
  TStopwatch sw;
  sw.Start();

  // Store the event number
  int previousEventNum = 0;
  int eventNum = 0;
  vector<double> backgroundDiscriminants;
  double signalDiscriminant;
  Long64_t signalTreeNum = 0;
  signalTree->GetEntries();
  for (Long64_t ievt=0; ievt<backgroundTree->GetEntries();ievt++) {

    if (ievt%1000 == 0)
      cout << "--- ... Processing event: " << ievt << endl;

    backgroundTree->GetEntry(ievt);

    // cout << "background variable 46 = " << backgroundReader[46] << endl;

    eventNum = int(backgroundReader[variablesNumber-1]);
    if ( previousEventNum == eventNum ) {
      // Still in the same event, fill the discriminant
      backgroundDiscriminants.push_back(backgroundReader.discriminant());
    }
    else {
      previousEventNum = eventNum;
      // combinatorics for this events finished. Now process the signal combination
      cout << "--- ... Processing combinatorials for event: " << backgroundReader[variablesNumber-1] << " of entry: " << signalTreeNum << endl;
      signalTree->GetEntry(signalTreeNum);
      // We can now increase the event number
      ++signalTreeNum;
//       cout << "variable 46 = " << signalReader[46] << endl;
      signalDiscriminant = signalReader.discriminant();
//       cout << "signalDiscriminant = " << signalDiscriminant << endl;

      // Compare signal to combinatorics and compute the efficiency to take the signal combination when
      // taking the one with the biggest discriminant.
      sort( backgroundDiscriminants.rbegin(), backgroundDiscriminants.rend() );
      int i = 0;
//       for( vector<double>::const_iterator it = backgroundDiscriminants.begin(); it != backgroundDiscriminants.end(); ++it, ++i ) {
//         cout << "double["<<i<<"] = " << *it << endl;
//       }
      // Compare the discrimant for the signal with the one for the combinatorial to comput the efficiency
      if( !(backgroundDiscriminants.empty()) ) {
        if( signalDiscriminant > backgroundDiscriminants.front() ) {
          ++efficiency;
        }
      }

      // Empty the vector before starting with the combinatiorial for a new event
      backgroundDiscriminants.clear();
    }
  }
  // Now process the last event, which is skipped by the previous loop
  // -----------------------------------------------------------------
  cout << "--- ... Processing combinatorials for event: " << backgroundReader[variablesNumber-1] << " of entry: " << signalTreeNum << endl;
  signalTree->GetEntry(signalTreeNum);
  // We can now increase the event number
  ++signalTreeNum;
  signalDiscriminant = signalReader.discriminant();
  // Compare signal to combinatorics and compute the efficiency to take the signal combination when
  // taking the one with the biggest discriminant.
  sort( backgroundDiscriminants.rbegin(), backgroundDiscriminants.rend() );
  int i = 0;
  // Compare the discrimant for the signal with the one for the combinatorial to comput the efficiency
  if( !(backgroundDiscriminants.empty()) ) {
    if( signalDiscriminant > backgroundDiscriminants.front() ) {
      ++efficiency;
    }
  }

  sw.Stop();
  cout << "--- End of event loop: "; sw.Print();

  if( signalTreeNum > 0 ) cout << "efficiency = " << efficiency/signalTreeNum << endl;
  else cout << "Error: events = " << signalTreeNum << endl;

  //
  // write histograms
  //
  TFile *target  = new TFile( "TMVApp.root","RECREATE" );
  signalReader.write();


  target->Close();

  cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << endl;

  cout << "==> TMVApplication is done!" << endl << endl;
} 
