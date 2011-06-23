#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <vector>
#include <memory>
#include "boost/shared_array.hpp"
#include "boost/shared_ptr.hpp"

/**
  * Used to compute efficiency.
  * An arbitrary number of variables can be defined and it provides methods to project the results over a subset.
  * The efficiencies are computed in an N dimensional matrix of size S = Prod(s_i) for i=1...N. Where s_i is the
  * number of bins of the i-th variable. Note that the last bin is alwyas treated as the overflow bin.
  * Internally it uses a linear representation of the matrix and it converts the indexes to the linearIndex.
  * Note that the methods operate with arrays passed by a shared_array smart pointer and require that the size
  * matches the one specified in the constructor. No checks are done for this consistency.
  */
class Efficiency
{
public:

  /**
    * Struct to pass parameters for the initialization of the Efficiency objects.
    */
  struct Parameters
  {
    Parameters(const unsigned int inputBins, const double & inputMin, const double & inputMax) :
      min(inputMin),
      max(inputMax),
      bins(inputBins),
    {}
    double min;
    double max;
    unsigned int bins;
  };

  /**
    * The constructor receives a vector with the sizes of all the variables.
    * The order of these variables must be preserved when they are filled later on.
    */
  Efficiency(const std::vector<Parameters> & vPars) :
    N_(vPars.size()),
    Nminus1_(vPars.size()-1),
    S_(1),
    vIndexesTempPtr_(0)
  {
    vSizes_.reset(new unsigned int[N_]);
    vBinSizes_.reset(new double[N_]);
    vMin_.reset(new double[N_]);
    vMax_.reset(new double[N_]);
    std::vector<Parameters>::const_iterator it = vPars.begin();
    unsigned int i = 0;
    for( ; it != vPars.end(); ++it, ++i ) {
      vSizes_[i] = it->bins;
      S_ *= it->bins;
      vMin_[i] = it->min;
      vMax_[i] = it->max;
      vBinSizes_[i] = (vMax_[i] - vMin_[i])/vSizes_[i];
    }
    values_.reset(new std::pair<unsigned int, unsigned int>[S_]);
    // This initialization should not be necessary as when the pair is built the integers are default initialized to 0.
    // for( unsigned int i=0; i<S_; ++i ) {
    //   values_[i] = std::pair<unsigned int, unsigned int>(0,0);
    // }
  }

  unsigned int getLinearIndex(const boost::shared_array<unsigned int> & vIndexes)
  {
    vIndexesTempPtr_ = vIndexes.get();
    return linearIndex(0, 0);
  }

  /**
    * Projects the values over a subset of the dimensions (reduces matrix size).
    * Accepts an array where the int is -1 for a variable to drop.
    * This method builds a new Efficiency object with the new sizes and values and returns a shared_ptr to it.
    */
  boost::shared_ptr<Efficiency> project(const boost::shared_array<int> & vKeep)
  {
    // Find the sizes to be kept
    std::vector<Parameters> newPars;
    for( unsigned int i=0; i<N_; ++i ) {
      if( vKeep[i] != -1 ) {
        newPars.push_back(Parameters(vSizes_[i], vMin_[i], vMax_[i]));
      }
    }
    boost::shared_ptr<Efficiency> newEff(new Efficiency(newPars));
    // Loop on all the values and fill them in the new Efficiency object.
    for( int i=0; i<S_; ++i ) {
      newEff.fill();
    }
    return newEff;
  }

  /**
    * Fills the counters for the efficiency.
    * The first counter is increased for all fills. The second counter is increased if accept == true.
    * At the end the efficiency is given by second/first.
    */
  void fill(const boost::shared_array<double> & variables, const bool accept)
  {
    // Compute the index for the given variable
    for( unsigned int i=0; i<N_; ++i ) {
      if( variables[i] >= vMax_[i] ) {
        values_[i].first += vSizes_[i]-1;
        if( accept ) values_[i].second += vSizes_[i]-1;
      }
      else {
        unsigned int index = (unsigned int)((variables[i] - vMin_[i])/vBinSizes_[i]);
        values_[i].first += index;
        if( accept ) values_[i].second += index;
      }
    }
  }

  /**
    * Alternate method taking an array of indexes and the values to add.
    * Used when reducing the matrix size by projecting.
    */
  void fill(const boost::shared_array<unsigned int> & variables, const std::pair<unsigned int, unsigned int> & values)
  {
    // Compute the index for the given variable
    for( unsigned int i=0; i<N_; ++i ) {
      if( variables[i] >= vMax_[i] ) {
        vIndexes_[i] = vSizes_[i];
        // values_[i].first += vSizes_[i]-1;
        // if( accept ) values_[i].second += vSizes_[i]-1;
      }
      else {
        unsigned int index = (unsigned int)((variables[i] - vMin_[i])/vBinSizes_[i]);
        vIndexes_[i] = index;
        // values_[i].first += index;
        // if( accept ) values_[i].second += index;
      }
    }
  }

  inline unsigned int linearSize() const {return S_;}
  inline unsigned int bins(unsigned int i) const {return vSizes_[i];}
  inline double binsSize(unsigned int i) const {return vBinSizes_[i];}
  inline double min(unsigned int i) const {return vMin_[i];}
  inline double max(unsigned int i) const {return vMax_[i];}

protected:
  /// Compute the linearIndex for a given vector of indexes
  unsigned int linearIndex(unsigned int i, unsigned int previous) const
  {
    if( i < Nminus1_ ) return( linearIndex(i+1, vIndexesTempPtr_[i] + vSizes_[i]*previous) );
    return (vIndexesTempPtr_[i] + vSizes_[i]*previous);
  }

  unsigned int N_;
  unsigned int Nminus1_;
  // Size of the linear representation
  unsigned int S_;

  boost::shared_array<unsigned int> vSizes_;
  boost::shared_array<double> vBinSizes_;
  boost::shared_array<double> vMin_;
  boost::shared_array<double> vMax_;
  unsigned int * vIndexesTempPtr_;
  boost::shared_array<unsigned int> vIndexes_;
  // Linear representation of the N-dimensional matrix
  boost::shared_array<std::pair<unsigned int, unsigned int> > values_;
};

#endif // EFFICIENCY_H
