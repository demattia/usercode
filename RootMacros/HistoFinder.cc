/**
 * This macro takes a file and the name of a histogram and look for
 * a histogram of the same name in all the subdirectories.
 * If it finds a match, it returns a pointer to it (but a TObject pointer
 * needs a cast to use it outside), if not returns
 * a null pointer and also does a cout saying it did not find it.
 * ATTENTION: if the same name matches more than one histogram in different
 * subdirectories, only the last match will be returned.
 */

// If you need to compile this macro, uncomment the following includes:
// --------------------------------------------------------------------
// Needed to use gROOT in a compiled macro
// #include "TROOT.h"
// #include <string.h>
// #include "TH1.h"
// #include "TChain.h"
// #include "TFile.h"
// #include "TTree.h"
// #include "TKey.h"
// #include "Riostream.h"
// #include "TCanvas.h"

class HistoFinder {
 public:
  HistoFinder() { histoFound_ = 0; }
  TObject * operator()( const TString & name, TDirectory * input ) {
    histoFound_ = 0;
    findHisto_( name, input );
    if( histoFound_ == 0 ) cout << "WARNING: histogram " << name << " not found" << endl;
    return histoFound_;
  }
 protected:
  void findHisto_( const TString & name, TDirectory * input ) {

    TIter nextkey( input->GetListOfKeys() );
    TKey *key, *oldkey=0;
    while ( (key = (TKey*)nextkey())) {

      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

      // read object from first source file
      TObject *obj = key->ReadObj();

      if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
        // If it matches the name, return it
        if( name == obj->GetName() ) histoFound_ = obj;
      }
      else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
        // it's a subdirectory
        // cout << "Found subdirectory " << obj->GetName() << endl;
        findHisto_( name, (TDirectory*)obj );
      }
      else {
        // object is of no type that we know or can handle
        cout << "Unknown object type, name: "
             << obj->GetName() << " title: " << obj->GetTitle() << endl;
      }
    }
  }
  TObject * histoFound_;
};
