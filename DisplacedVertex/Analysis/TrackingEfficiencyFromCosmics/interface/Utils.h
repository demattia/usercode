#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <map>

namespace utils
{
  TH1F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name,
		       const TString & nameForTitle, const TString axisTitle, const TString & unit,
		       const int bins, const double & min, const double & max)
  {
    TString underscore("");
    if( nameForTitle != "" ) underscore = "_";

    TString spaceAndUnit("");
    if( unit != "" ) spaceAndUnit = " " + unit;

    TH1F * histo = fileService->make<TH1F>(name+underscore+nameForTitle, name+" "+nameForTitle, bins, min, max);
    histo->GetXaxis()->SetTitle(axisTitle+spaceAndUnit);
    std::stringstream ss;
    ss << histo->GetBinWidth(1);
    histo->GetYaxis()->SetTitle("entries / ( "+ss.str() + spaceAndUnit + " )");

    return histo;
  }
}

#endif
