#include "treeAnalyzer.C"
//#include "muonAnalyzer.C"
#include <map>
#include <TString.h>
#include "FilesAndWeights.h"

void fullAnalyzer(const bool electrons = false, const bool trackAnalysis = true, const bool useSAMuons = false)
{
  // passing file name, weight and if it is electrons or muons.
  // The weights are the cross sections in pb
  std::map<TString, double> fw(filesAndWeightsMap( electrons ));
  std::map<TString, double>::const_iterator it = fw.begin();
  for( ; it != fw.end(); ++it ) {
    std::cout << "processing: " << it->first << std::endl;

    if ( trackAnalysis ) {
      treeAnalyzer analyzer(it->first+"/histograms.root", it->second, electrons);
      analyzer.Loop();
    }
//    else {
////      muonAnalyzer analyzer(it->first+"/histograms.root", it->second, useSAMuons, electrons);
////      analyzer.Loop();
//      if ( useSAMuons ) {}
//    }
  }
}
