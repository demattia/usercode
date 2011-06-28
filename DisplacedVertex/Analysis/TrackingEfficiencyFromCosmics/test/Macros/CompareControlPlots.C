#include <TFile.h>
#include <TH1F.h>
#include <TObject.h>
#include <TList.h>
#include <TKey.h>
#include <TClass.h>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>
#include <TCanvas.h>

/**
  * Loop over all keys in a file and save pointers to objects inheriting from TH1* to the map.
  * Go recursively in subdirectories.
  */
struct SaveHistogramList
{
  void saveAllHistograms(TDirectory * inputDir, const TString & identifier = "")
  {
    TIter next(inputDir->GetListOfKeys());
    TKey * key;
    while( (key = (TKey*)next()) ) {
      TObject * obj = key->ReadObj();
      if( obj->IsA()->InheritsFrom("TH1") ) {
        TString nameForMap(obj->GetName());
        if( identifier == "" || nameForMap.BeginsWith(identifier) ) {
          if( identifier != "" ) {
            nameForMap.Remove(nameForMap.Index(identifier), identifier.Length()+1);
          }
          histoMap[nameForMap].push_back((TH1*)obj);
        }
      }
      else if( (obj != inputDir) && (obj->IsA()->InheritsFrom(TDirectory::Class())) ) {
        saveAllHistograms((TDirectory*)obj, identifier);
      }
    }
  }
  std::map<TString, std::vector<TH1*> > histoMap;
};

void CompareControlPlots(const TString & inputFileName = "/home/demattia/TrackingEfficiencyFromCosmics/Simulation/TrackingEfficiencyFromCosmics.root")
{
  TFile * inputFile = new TFile(inputFileName);

  SaveHistogramList saveHistogramList;
  saveHistogramList.saveAllHistograms(inputFile, "generalTracks");
  saveHistogramList.saveAllHistograms(inputFile, "standAloneMuons");

  std::map<TString, std::vector<TH1*> >::const_iterator vit = saveHistogramList.histoMap.begin();
  for( ; vit != saveHistogramList.histoMap.end(); ++vit ) {
    TCanvas * canvas = new TCanvas;
    canvas->Divide(1,2);
    int index = 1;
    for( std::vector<TH1*>::const_iterator it = vit->second.begin(); it != vit->second.end(); ++it, ++index ) {
      canvas->cd(index);
      (*it)->Draw();
    }
    canvas->Print(vit->first+".pdf");
  }
}

