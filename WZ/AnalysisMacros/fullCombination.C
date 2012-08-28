#include <map>
#include <TString.h>
#include <iostream>
#include <TH1F.h>
#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TClass.h>
#include <TChain.h>
#include "FilesAndWeights.h"
#include "MergeFiles.h"
#include "GetCategory.h"

void fullCombination(const bool electrons = false)
{
  TString type("_muons");
  if( electrons ) type = "_electrons";

  std::map<TString, TList*> cathegories;
  // Loop over all files
  std::map<TString, double> fw(filesAndWeightsMap( electrons ));
  std::map<TString, double>::const_iterator it = fw.begin();
  for( ; it != fw.end(); ++it ) {
    TString cathegory(getCategory(it->first));
    if( cathegories[cathegory] == 0 ) cathegories[cathegory] = new TList();
    std::cout << "fileName = " << std::string(it->first).substr(std::string(it->first).find_last_of("/")+1)+"weighted"+type+".root" << std::endl;
    cathegories[cathegory]->Add( TFile::Open(std::string(it->first).substr(std::string(it->first).find_last_of("/")+1)+"weighted"+type+".root") );
  }

  TList * combinationsList = new TList();
  std::map<TString, TList*>::const_iterator it2 = cathegories.begin();
  for( ; it2 != cathegories.end(); ++it2 ) {
    TFile *Target = TFile::Open( it2->first+"_combined"+type+".root", "RECREATE" );
    std::cout << "fileName = " << it2->first << std::endl;
    combinationsList->Add( Target );
    mergeFiles( Target, it2->second );
  }
}
