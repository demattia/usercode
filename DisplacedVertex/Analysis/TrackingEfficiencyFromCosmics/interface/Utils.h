#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <map>

namespace utils
{
  TH1F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name,
		       const TString & nameForTitle, const TString axisTitle, const TString & unit,
		       const int bins, const double min, const double max)
  {
    TString underscore("");
    if( nameForTitle != "" ) underscore = "_";

    TH1F * histo = fileService->make<TH1F>(name+underscore+nameForTitle, name+" "+nameForTitle,500,0,500);
    histo->GetXaxis()->SetTitle(axisTitle+" "+unit);
    std::stringstream ss;
    ss << histo->GetBinWidth(1);
    histo->GetYaxis()->SetTitle("entries / ( "+ss.str() + unit + " )");

    return histo;
  }
}

#endif
