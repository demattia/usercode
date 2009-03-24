#ifndef ClusterizerUnitTesterTest_h
#define ClusterizerUnitTesterTest_h
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithm.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithmFactory.h"
#include "RecoLocalTracker/SiStripClusterizer/test/UnitTester.h"
#include "RecoLocalTracker/SiStripClusterizer/test/TestObjectWrappers.h"

#include <map>
#include <string>

using namespace std;

class ClusterizerUnitTesterTest : public UnitTester {
  
  typedef edm::ParameterSet PSet;
  typedef std::vector<PSet> VPSet;
  typedef VPSet::const_iterator iter_t;
  typedef edmNew::DetSetVector<SiStripCluster> output_t;
  typedef std::auto_ptr<StripClusterizerAlgorithm> algo_t;

 public:
  
  ClusterizerUnitTesterTest(const PSet& conf) : UnitTester( conf ) {}
  ~ClusterizerUnitTesterTest() {}

 protected:

  virtual baseObjectWrapper * objectFactory(const PSet& test, const edm::EventSetup& es);
};

#endif
