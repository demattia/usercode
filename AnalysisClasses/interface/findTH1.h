/**

  This macro will add histograms multiplied by their cross-section
  from a list of root files and write them
  to a target root file. The target file is newly created and must not be
  identical to one of the source files.


  Author: Sven A. Schmidt, sven.schmidt@cern.ch
  Date:   13.2.2001

  Editing Author: Michael B. Anderson, mbanderson@hep.wisc.edu
  Date:  July 12, 2007

  This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
  which had a problem with directories more than one level deep.
  (see macro hadd_old.C for this previous implementation).

  The macro from Sven has been enhanced by
     Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
   to automatically add Trees (via a chain of trees).

  To use this macro, modify the file names in function hadd.

  NB: This macro is provided as a tutorial.
      Use $ROOTSYS/bin/hadd to merge many histogram files

 */


#include <string.h>
#include <sstream>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

void findTH1( std::vector<TH1*>    & outputVector,
	      TList                * sourcelist, 
	      TString                variableName, 
	      std::vector<TString> & sourceSuffixVector,
	      std::vector<double>  * xSecVector = 0  
	      ) {
  
  TFile *first_source = (TFile*)sourcelist->First();

  TString path( (char*)strstr( first_source->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  // loop over all keys in this directory
  TChain *globChain = 0;
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

      bool foundName = false;
      if ( variableName == obj->GetName() ) 
	foundName=true;
      
      if ( !foundName ) continue;
      cout << "finding histogram " << variableName << endl;
      TH1 *h1 = (TH1*)obj->Clone();

      std::ostringstream sourceIndex;
      TString h1Name = h1->GetName();

      // Scale by the cross-section factor
      if (xSecVector) h1->Scale((*xSecVector)[0]);
      sourceIndex << 0;
      //      h1->SetName(h1Name+"_"+sourceIndex.str());
      h1->SetName(h1Name+"_"+sourceSuffixVector[0]);
      outputVector.push_back(h1);
      sourceIndex.str("");

      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );

      int q = 1; // This keeps track of which
                 // cross section factor to use

      while ( nextsource ) {
        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(variableName);
        if (key2) {
           TH1 *h2 = (TH1*)(key2->ReadObj())->Clone();
	   // Scale by the cross section factor
           // before adding.
           if (xSecVector) h2->Scale((*xSecVector)[q]);
	   sourceIndex << q;
	   h2->SetName(variableName+"_"+sourceSuffixVector[q]);
	   outputVector.push_back(h2);
	   sourceIndex.str("");
           q++;
        }
        nextsource = (TFile*)sourcelist->After( nextsource );
      }
    }
    else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
      
      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();
      
      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
      while ( nextsource ) {
	
	globChain->Add(nextsource->GetName());
	nextsource = (TFile*)sourcelist->After( nextsource );
      }
    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory
      cout << "Found subdirectory " << obj->GetName() << endl;
      
      // newdir is now the starting point of another round of finding
      findTH1( outputVector, 
	       sourcelist, 
	       variableName, 
	       sourceSuffixVector,
	       xSecVector
	       );
      
    } else {
      
      // object is of no type that we know or can handle
      cout << "Unknown object type, name: "
	   << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }
  }
}
