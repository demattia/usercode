#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <vector>
#include <memory>
#include <cmath>
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
      bins(inputBins)
    {}
    double min;
    double max;
    unsigned int bins;
  };

  /// This is used only when building from a SerializedEfficiency
  Efficiency() {}

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
    vIndexes_.reset(new unsigned int[N_]);
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
    return computeLinearIndex(0, 0);
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
    for( unsigned int i=0; i<S_; ++i ) {
      // if( values_[i].first != 0 ) {
      //   std::cout << "Filling i = " << i << " with (" << values_[i].first << ", " << values_[i].second << ")" << std::endl;
      //   for( unsigned int i=0; i<newPars.size(); ++i) {
      //     std::cout << "newIndex("<<i<<") = " << getIndexes(i, vKeep, newPars.size()) << std::endl;
      //   }
      // }
      newEff->fill(getIndexes(i, vKeep, newPars.size()), values_[i]);
    }
    return newEff;
  }

  /**
    * Method to compute the indexes starting from the linearIndex. It is used when reducing the matrix dimensions by projecting.
    * The computation stops when the last variable needed is reached.
    * It takes:
    * - the value of the linearIndex
    * - the array with the selection of variables to keep (see project method)
    * - the number of variables to keep
    */
  boost::shared_array<unsigned int> getIndexes(const unsigned int linearIndex, const boost::shared_array<int> & vKeep, const unsigned int newN)
  {
    boost::shared_array<unsigned int> vNewIndexes(new unsigned int[newN]);
    unsigned int tempS = S_;
    unsigned int newCounter = 0;
    unsigned int counter = 0;
    computeIndexes(linearIndex, tempS, vKeep, newN, vNewIndexes, counter, newCounter);
    return vNewIndexes;
  }

  /**
    * Compute the indexes only for the variables to be kept and save them in the vNewIndexes.
    * It is used when reducing the matrix dimensions by projecting.
    */
  void computeIndexes(const unsigned int L, unsigned int & S, const boost::shared_array<int> & vKeep,
                      const unsigned int newN, boost::shared_array<unsigned int> & vNewIndexes,
                      unsigned int & counter, unsigned int & newCounter)
  {
    if( newCounter == newN ) return;

    S /= vSizes_[counter];

    if( vKeep[counter] != -1 ) {
      vNewIndexes[newCounter] = (unsigned int)(L/S);
      ++newCounter;
    }
    ++counter;
    computeIndexes(L%S, S, vKeep, newN, vNewIndexes, counter, newCounter);
  }

  /**
    * Fills the counters for the efficiency.
    * The first counter is increased for all fills. The second counter is increased if accept == true.
    * At the end the efficiency is given by second/first.
    */
  void fill(const boost::shared_array<double> & variables, const bool accept)
  {
    // Compute the indexes for the given values
    for( unsigned int i=0; i<N_; ++i ) {
      if( variables[i] >= vMax_[i] ) {
        vIndexes_[i] = vSizes_[i] - 1;
        // std::cout << "vIndexes_["<<i<<"] = " << vIndexes_[i] << std::endl;
      }
      else if( variables[i] < vMin_[i] ) {
        vIndexes_[i] = 0;
      }
      else {
        vIndexes_[i] = (unsigned int)((variables[i] - vMin_[i])/vBinSizes_[i]);
        // std::cout << "vIndexes_["<<i<<"] = " << vIndexes_[i] << std::endl;
      }
    }
    values_[getLinearIndex(vIndexes_)].first += 1;
    if( accept ) values_[getLinearIndex(vIndexes_)].second += 1;
  }

  /**
    * Alternate method taking an array of indexes and the values to add.
    * Used when reducing the matrix size by projecting.
    */
  inline void fill(const boost::shared_array<unsigned int> & indexes, const std::pair<unsigned int, unsigned int> & values)
  {
    values_[getLinearIndex(indexes)].first += values.first;
    values_[getLinearIndex(indexes)].second += values.second;
  }

  double getEff(const boost::shared_array<unsigned int> & vIndexes)
  {
    unsigned int linearIndex = getLinearIndex(vIndexes);
    if( (values_[linearIndex].first) == 0 ) return -1;
    return( double(values_[linearIndex].second)/double(values_[linearIndex].first) );
  }

  double getEffError(const boost::shared_array<unsigned int> & vIndexes)
  {
    unsigned int linearIndex = getLinearIndex(vIndexes);
    if( (values_[linearIndex].first) == 0 ) return 0;
    double p = double(values_[linearIndex].second)/double(values_[linearIndex].first);
    return( sqrt(p*(1-p)/(double(values_[linearIndex].first))) );
  }

  /// Overloaded for the case of a single variable. No check is performed and the index is assumed == linearIndex.
  double getEff(unsigned int index)
  {
    if( (values_[index].first) == 0 ) return -1;
    return( double(values_[index].second)/double(values_[index].first) );
  }

  /// Overloaded for the case of a single variable. No check is performed and the index is assumed == linearIndex.
  double getEffError(unsigned int index)
  {
    if( (values_[index].first) == 0 ) return 0;
    double p = (values_[index].second)/double(values_[index].first);
    return( sqrt(p*(1-p)/(double(values_[index].first))) );
  }

  inline unsigned int getLinearSize() const {return S_;}
  inline unsigned int getN() const {return N_;}
  inline unsigned int getVIndexes(const unsigned int i) const {return vIndexes_[i];}
  inline std::pair<unsigned int, unsigned int> getValues(const unsigned int i) const {return values_[i];}
  inline unsigned int bins(const unsigned int i) const {return vSizes_[i];}
  inline double binsSize(const unsigned int i) const {return vBinSizes_[i];}
  inline double min(const unsigned int i) const {return vMin_[i];}
  inline double max(const unsigned int i) const {return vMax_[i];}

  // Setter methods for building from serialized objects
  inline void setLinearSize(const unsigned int S) {S_ = S;}
  inline void setN(const unsigned int N) {N_ = N;}
  inline void setVIndexes(boost::shared_array<unsigned int> & vIndexes) {vIndexes_ = vIndexes;}
  inline void setVSizes(boost::shared_array<unsigned int> & vSizes) {vSizes_ = vSizes;}
  inline void setVMin(boost::shared_array<double> & vMin) {vMin_ = vMin;}
  inline void setVMax(boost::shared_array<double> & vMax) {vMax_ = vMax;}
  inline void setVBinSizes(boost::shared_array<double> & vBinSizes) {vBinSizes_ = vBinSizes;}
  inline void setValues(boost::shared_array< std::pair<unsigned int, unsigned int> > & values) {values_ = values;}

protected:
  /// Compute the linearIndex for a given vector of indexes
  unsigned int computeLinearIndex(unsigned int i, unsigned int previous) const
  {
    if( i < Nminus1_ ) return( computeLinearIndex(i+1, vIndexesTempPtr_[i] + vSizes_[i]*previous) );
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
