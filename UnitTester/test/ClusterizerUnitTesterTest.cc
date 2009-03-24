#include "RecoLocalTracker/SiStripClusterizer/test/ClusterizerUnitTesterTest.h"

#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include <functional>
#include <numeric>
#include <vector>
#include <iostream>
#include <sstream>

baseObjectWrapper * ClusterizerUnitTesterTest::objectFactory(const PSet& test, const edm::EventSetup& es)
{
  // PSet must contain: string label, VPSet parameters and alternatively a base type, another VPSet or nothing.
  // In the case of a basic type it is not necessary to provide an objectWrapper.

  string type = test.getParameter<string>("Label");

  cout << "type = " << type << endl;

  if( type == "Cluster" ) return new clusterWrapper(test);
  if( type == "Algo" ) return new algoWrapper(test, es);
  throw(string("unidentifiedType"));

  return 0;
}
