#include "RecoLocalTracker/SiStripClusterizer/test/TestObjectWrappers.h"

bool ClusterWrapper::clusterDSVsIdenticalTest(const output_t& L, const output_t& R)
{
//   cout << "L.size() = " << L.size() << ", R.size() = " << R.size() << endl;
  return
    L.size() == R.size() &&
    inner_product(L.begin(), L.end(), R.begin(),
		  bool(true), std::logical_and<bool>(), clusterDetSetsIdenticalTest );  
}

bool ClusterWrapper::clusterDetSetsIdenticalTest(const edmNew::DetSet<SiStripCluster>& L, const edmNew::DetSet<SiStripCluster>& R)
{
  return 
    L.size() == R.size() &&
    inner_product(L.begin(), L.end(), R.begin(),
		  bool(true), std::logical_and<bool>(), clustersIdenticalTest );
}

bool ClusterWrapper::clustersIdenticalTest(const SiStripCluster& L, const SiStripCluster& R)
{
//   cout << "L.geographicalId() = " << L.geographicalId() << ", R.geographicalId() = " << R.geographicalId() << endl;
//   cout << "L.firstStrip() = " << L.firstStrip() << ", R.firstStrip() = " << R.firstStrip() << endl;
//   cout << "L.amplitudes().size() = " << L.amplitudes().size() << ", R.amplitudes().size() = " << R.amplitudes().size() << endl;
//   cout << "inner_product(L.amplitudes().begin(), L.amplitudes().end(), R.amplitudes().begin(), bool(true), std::logical_and<bool>(), std::equal_to<uint16_t>() ) = " << inner_product(L.amplitudes().begin(), L.amplitudes().end(), R.amplitudes().begin(), bool(true), std::logical_and<bool>(), std::equal_to<uint16_t>() ) << endl;
  return
    L.geographicalId() == R.geographicalId()
    && L.firstStrip() == R.firstStrip() 
    && L.amplitudes().size() == R.amplitudes().size()
    && inner_product(L.amplitudes().begin(), L.amplitudes().end(), R.amplitudes().begin(), 
		     bool(true), std::logical_and<bool>(), std::equal_to<uint16_t>() );
}
