#ifndef MuScleFitLikelihoodPdf_h
#define MuScleFitLikelihoodPdf_h

#include "CondFormats/PhysicsToolsObjects/interface/Histogram2D.h"
#include <vector>
#include <string>

struct MuScleFitLikelihoodPdf
{
  std::vector<PhysicsTools::Calibration::HistogramD2D> histograms;
  std::vector<std::string> names;
  std::vector<int> xBins;
  std::vector<int> yBins;
};

#endif // MuScleFitLikelihoodPdf
