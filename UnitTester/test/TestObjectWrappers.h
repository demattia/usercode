#ifndef TestObjectWrappers_h
#define TestObjectWrappers_h

#include <iostream>
#include <boost/any.hpp>

/**
 * boost::any is used to pass the objects contained in classes derived from baseObjectWrapper. To use the method
 * we need to define it in the base classe, but the return type depends on the template parameter. Thus we return
 * a boost::any object containing the actual object. The operator== willing to access the object (and knowing its
 * type) can do the conversion back with boost::any_cast<type>(anyObject).
 */

using namespace std;

namespace {
  typedef edm::ParameterSet PSet;
  typedef std::vector<PSet> VPSet;
  typedef VPSet::const_iterator iter_t;
}

/// Empty base class only used to have a polymorphic pointer
class baseObjectWrapper {
public:
  virtual bool operator==( const baseObjectWrapper & other ) const = 0;
  virtual bool operator!=( const baseObjectWrapper & other ) const { return !(*this == other); }
  virtual baseObjectWrapper * eval(const PSet & parameters) { return this; }
  virtual boost::any object() const = 0;
};

/**
 * This class will hold the tested object
 * The object must be assignable
 */
template <class T>
class objectWrapper : public baseObjectWrapper
{
public:
  /// Not initializing empty constructor used by derived classes filling the object_ on their own.
  objectWrapper() {}
  objectWrapper(const T & test) : object_(test) {}
  virtual ~objectWrapper() {}
  virtual bool operator==( const baseObjectWrapper & other ) const
  {
    cout << "Warning: using the default operator== which always returns false. It must be overridden in derived classes" << endl;
    return false;
  }
//   {
//     return object_ == boost::any_cast<T>(other.object());
//   }
  boost::any object() const { return boost::any(object_); }
protected:
  T object_;
};



// Specializations
#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithm.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithmFactory.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
namespace {
  typedef edmNew::DetSetVector<SiStripCluster> output_t;
  typedef boost::shared_ptr<StripClusterizerAlgorithm> algo_t;
}

// For the expected cluster
class clusterWrapper : public objectWrapper<output_t>
{
public:
  clusterWrapper(const PSet & test)
  {
    constructClusters(test.getParameter<PSet>("Parameters"), object_);
  }
  virtual ~clusterWrapper() {}

  virtual bool operator==(const baseObjectWrapper & other) const
  {
    cout << "clusterWrapper operator==" << endl;
    output_t otherClusterSet = boost::any_cast<output_t>(other.object());

    // Check if they have the same number of clusters
    if( otherClusterSet.size() != object_.size() ) return false;

    // Check if the clusters are the same
    output_t::const_iterator otherCluster = otherClusterSet.begin();
    output_t::const_iterator thisCluster = object_.begin();
    for( ; otherCluster != otherClusterSet.end(); ++otherCluster, ++thisCluster ) {
      if( otherCluster->detId() != thisCluster->detId() ) return false;
    }
    return true;
  }
protected:
  void constructClusters(const PSet& parameters, output_t& clusters) {
    int detId = parameters.getParameter<uint32_t>("DetId");
    cout << "Building a clusterSet with detId = " << detId << endl;
    VPSet clusterset = parameters.getParameter<VPSet>("ClusterSet");
    output_t::FastFiller clustersFF(clusters, detId);
    for(iter_t c = clusterset.begin(); c<clusterset.end(); ++c) {
      uint16_t firststrip =  c->getParameter<unsigned>("FirstStrip");
      std::vector<unsigned> amplitudes =  c->getParameter<std::vector<unsigned> >("Amplitudes");
      std::vector<uint16_t> a16(amplitudes.begin(),amplitudes.end());
      clustersFF.push_back(SiStripCluster(detId, firststrip, a16.begin(),a16.end()));
    }
    if(clustersFF.empty()) clustersFF.abort();
  }
};

// For the algorithm
class algoWrapper : public objectWrapper<algo_t>
{
public:
  algoWrapper(const PSet & test, const edm::EventSetup & es)
  {
    // Cannot store an auto_ptr inside the boost::any container. We extract the pointer and save it into
    // a shared_ptr (which is copiable).
    object_.reset(StripClusterizerAlgorithmFactory::create(test.getParameter<PSet>("Parameters")).release());
    object_->initialize(es);
  }
  virtual ~algoWrapper() {}

  virtual bool operator==(const baseObjectWrapper & other) const
  {
    cout << "Algo wrapper operator ==" << endl;
    return true;
  }

  // This is a fixture which also needs to receive parameters for the tests. We redefine the
  // eval methods accordingly.
  virtual baseObjectWrapper * eval(const PSet & parameters)
  {
    return 0;
  }
protected:
  void constructDigis(const VPSet& stripset, edmNew::DetSetVector<SiStripDigi>& digis, const int detId)
  {
    edmNew::DetSetVector<SiStripDigi>::FastFiller digisFF(digis, detId);
    for(iter_t strip = stripset.begin(); strip < stripset.end(); strip++) {
      digisFF.push_back( SiStripDigi(strip->getParameter<unsigned>("Strip"),
                                     strip->getParameter<unsigned>("ADC") ));
    }
    if(digisFF.empty()) digisFF.abort();
  }
};

#endif
