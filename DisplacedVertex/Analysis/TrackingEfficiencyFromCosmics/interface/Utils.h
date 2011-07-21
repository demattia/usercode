#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <map>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>

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

  TH2F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name, const TString & nameForTitle,
		       const TString xAxisTitle, const TString & xUnit,
		       const int xBins, const double & xMin, const double & xMax,
		       const TString yAxisTitle, const TString & yUnit,
		       const int yBins, const double & yMin, const double & yMax)
  {
    TString underscore("");
    if( nameForTitle != "" ) underscore = "_";

    TString spaceAndUnit("");
    if( xUnit != "" ) spaceAndUnit = " " + xUnit;

    TH2F * histo = fileService->make<TH2F>(name+underscore+nameForTitle, name+" "+nameForTitle, xBins, xMin, xMax, yBins, yMin, yMax);
    histo->GetXaxis()->SetTitle(xAxisTitle+spaceAndUnit);
    // std::stringstream ss;
    // ss << histo->GetXaxis()->GetXBinWidth(1);
    // histo->GetYaxis()->SetTitle("entries / ( "+ss.str() + spaceAndUnit + " )");
    // ss.str("");
    // ss << histo->GetYaxis()->GetYBinWidth(1);
    if( yUnit != "" ) spaceAndUnit = " " + yUnit;
    histo->GetYaxis()->SetTitle(yAxisTitle+spaceAndUnit);

    return histo;
  }

  /// Polynomial function fitted to the dxy_error vs Pt distribution
  double polynomial( const double & pt )
  {
    // return ( 459.572 - 1222.24*pt + 1423.85*pt*pt -756.519*pt*pt*pt + 150.481*pt*pt*pt*pt );
    return ( 0.1+1.91364 - 0.0211496*pt + 0.0000906055*pt*pt - 0.000000130650*pt*pt*pt );
  }

  /// Returns the maximum dxy error value allowed for the cuts
  double dxyErrMax( const double & pt )
  {
    double dxyErrMax = 1.;
    if(pt < 200 ) dxyErrMax = polynomial(pt);
    else dxyErrMax = polynomial(200);
    return std::min(dxyErrMax, 1.);
  }
}

#endif
