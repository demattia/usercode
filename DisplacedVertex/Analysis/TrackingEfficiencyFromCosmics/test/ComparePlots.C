/**
 * This macro takes two histogram files and compare them.
 * It produces a "compared.root" file with the corresponding
 * histograms superimposed (2nd in red).
 * It also gives an output with the result of the difference between
 * the corresponding histograms. If everything is zero, the two files
 * have the same histograms.
 */

// Needed to use gROOT in a compiled macro
#include "TROOT.h"

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "THStack.h"
#include "TCanvas.h"

#include "tdrstyle.C"

TList *FileList;
TFile *Target;
unsigned int countDifferences = 0;

void CompareRootfile( TDirectory *target, TList *sourcelist );

void ComparePlots() {

  setTDRStyle();

  // in an interactive ROOT session, edit the file names
  // Target and FileList, then
  // root > .L hadd.C
  // root > hadd()

  gROOT->SetBatch(true);

  Target = TFile::Open( "compared.root", "RECREATE" );

  FileList = new TList();

  // ************************************************************
  // List of Files
  FileList->Add( TFile::Open("TrackingEfficiencyFromCosmics_Data.root") );    // 1
  FileList->Add( TFile::Open("TrackingEfficiencyFromCosmics_Sim.root") );    // 2

  CompareRootfile( Target, FileList );

  cout << endl << "Number of different histograms: " << countDifferences << endl;
}

void CompareRootfile( TDirectory *target, TList *sourcelist ) {

  // Create a THStack to draw the histograms superimposed
  THStack * comp = 0;

  //  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile *first_source = (TFile*)sourcelist->First();

  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> merge it

      // cout << "Comparing histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;

      Double_t norm = h1->GetEntries();
      if (norm) h1->Scale(1/norm);

      comp = new THStack(TString(h1->GetName()) + "Stack", TString(h1->GetTitle()) + " stack");

      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      comp->Add(h1);

      bool found = false;
      while ( nextsource ) {
        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
        if (key2) {
          found = true;
          TH1 *h2 = (TH1*)key2->ReadObj();
          h2->SetLineColor(kRed);
	  norm = h2->GetEntries();
	  if (norm) h2->Scale(1/norm);
          comp->Add(h2);

          if ( h1->GetNbinsX() != h2->GetNbinsX() ) {
            cout << "ERROR: histograms " << h1->GetName() << " have different bins on the x axis" << endl;
            exit(1);
          }
          // 0 is the underflow and nBins+1 is the overflow.
          double diff = 0.;
          for( int iBin = 1; iBin <= h1->GetNbinsX(); ++iBin ) {

            // Using long double for precision. This could still fail for approximations
            long double binH1 = h1->GetBinContent(iBin);
            long double binH2 = h2->GetBinContent(iBin);

            // if ( string(h1->GetName()).find("hLikeVSMu_LikelihoodVSPt_prof") != string::npos ) {
            //   cout << "h1["<<iBin<<"] = " << binH1 << ", h2["<<iBin<<"] = " << binH2 << endl;
            // }

            // Sum only if the bin contents are different (to avoid approximation errors)
            //if( binH1 != binH2 ) {
              diff = binH1 - binH2;
              // cout << "binH1 = " << binH1 << ", binH2 = " << binH2 << endl;
            //}

          }
          cout << h1->GetName() << " difference = " << diff << endl;
          if( diff != 0 ) ++countDifferences;
        }
        nextsource = (TFile*)sourcelist->After( nextsource );
      }
      if (!found) cout << "ERROR: no matching histogram found in the second file" << endl;
    }
    else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      CompareRootfile( newdir, sourcelist );
    }
    else {
      // object is of no type that we know or can handle
      cout << "Unknown object type, name: "
           << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }

    // now write the compared histograms (which are "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();

      if( obj->IsA()->InheritsFrom( "TH1" ) ) {
        // Write the superimposed histograms to file
	//obj->Write( key->GetName() );
        TCanvas canvas( TString(obj->GetName())+"Canvas", TString(obj->GetName())+" canvas", 1000, 800 );
        comp->Draw("nostack");
	comp->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
	comp->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
        canvas.Write();
	canvas.SaveAs("comparedPlots/"+TString(obj->GetName())+".pdf");
	canvas.SaveAs("comparedPlots/"+TString(obj->GetName())+".png");
        comp->Write();
      }
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);
}
