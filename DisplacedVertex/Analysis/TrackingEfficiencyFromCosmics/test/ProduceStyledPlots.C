/**
 * This macro loops on all the histograms in a given file and
 * draws them again with some modified settings.
 */

// Needed to use gROOT in a compiled macro
#include "TROOT.h"
// #include "TStyle.h"

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
// #include "gstyle.C"

TList *FileList;
TFile *Target;
void draw( TDirectory *target, TList *sourcelist );

void ProduceStyledPlots(const TString & inputFile = "TrackingEfficiencyFromCosmics.root", const TString & type = "")
{
  TString underscoreType("");
  if( type != "" ) underscoreType = "_"+type;
  gROOT->SetBatch(true);
  //   gROOT->SetStyle("Plain");
  //   gStyle->SetCanvasColor(kWhite);
  //   gStyle->SetCanvasBorderMode(0);
  //   gStyle->SetPadBorderMode(0);
  //   gStyle->SetTitleFillColor(kWhite);
  //   gStyle->SetTitleColor(kWhite);

  setTDRStyle();

  Target = TFile::Open( "redrawed.root", "RECREATE" );

  FileList = new TList();

  // ************************************************************
  // List of Files
  FileList->Add( TFile::Open(inputFile) );    // 1

  draw( Target, FileList, underscoreType );
}

void draw( TDirectory *target, TList *sourcelist, const TString & type )
{
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
      // descendant of TH1 -> redraw it

      // TH1 *h1 = (TH1*)obj;

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
      draw( newdir, sourcelist, type );
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
        // obj->Write( key->GetName() );
	TString objName(obj->GetName());
        TCanvas canvas( objName+"Canvas", TString(obj->GetName())+" canvas", 1000, 800 );
        obj->Draw();
        // If the update is not called, the gStyle attributes are not used.
        canvas.Update();
        canvas.Write();
	canvas.SaveAs("plots"+type+"/"+TString(obj->GetName())+".pdf");
	canvas.SaveAs("plots"+type+"/"+TString(obj->GetName())+".png");
      }
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);
}
