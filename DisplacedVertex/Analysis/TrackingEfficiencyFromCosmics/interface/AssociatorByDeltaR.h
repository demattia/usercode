#ifndef ASSOCIATORBYDELTAR_H
#define ASSOCIATORBYDELTAR_H

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

/**
  * The matching is done in the following way:
  * - loop on the firstCollection
  * - find the closest element in the secondCollection and take it as matching if deltaR < maxDeltaR
  * - remove this element from the secondCollection and move to the next in the first collection
  */
class AssociatorByDeltaR
{
public:
  AssociatorByDeltaR(const double & maxDeltaR, const bool removeMatches = true, const bool takeOppositeTracks = false) :
    maxDeltaR_(maxDeltaR),
    removeMatches_(removeMatches),
  takeOppositeTracks_(takeOppositeTracks)
  {}

  template <class T1, class T2>
  void fillAssociationMap( const std::vector<T1> & firstCollection,
                           const std::vector<T2> & secondCollection,
                           std::map<const T1 *, const T2 *> & matchesMap,
                           TH1 * hMinDeltaR = 0 )
  {
    // Clone the second collection
    typename std::vector<const T2 *> clonedSecond;
    typename std::vector<T2>::const_iterator itSecond = secondCollection.begin();
    for( ; itSecond != secondCollection.end(); ++itSecond ) {
      clonedSecond.push_back(&*itSecond);
    }

    // Loop over the first collection
    typename std::vector<T1>::const_iterator itFirst = firstCollection.begin();
    for( ; itFirst != firstCollection.end(); ++itFirst ) {
      double minDeltaR = maxDeltaR_;
      typename std::vector<const T2 *>::iterator matched = clonedSecond.end();
      double deltaR = maxDeltaR_;
      typename std::vector<const T2 *>::iterator itSecond = clonedSecond.begin();
      for( ; itSecond != clonedSecond.end(); ++itSecond ) {
        deltaR = reco::deltaR(*itFirst, **itSecond);
        if( deltaR < minDeltaR ) {
          minDeltaR = deltaR;
          matched = itSecond;
        }
        else if( takeOppositeTracks_ ) {
          if( (deltaR > (3.14 - minDeltaR)) && (deltaR < 3.14 + minDeltaR) ) {
            minDeltaR = deltaR;
            matched = itSecond;
          }
        }
      }
      if( matched != clonedSecond.end() ) {
        matchesMap.insert(std::make_pair(&*itFirst, *matched));

        if( removeMatches_ ) clonedSecond.erase(matched);

        if( hMinDeltaR != 0 ) hMinDeltaR->Fill(deltaR);
      }
      else {
        T2 * emptyPointer = 0;
        matchesMap.insert(std::make_pair(&*itFirst, emptyPointer));
      }
    }
  }

  /// Specialization to cope with SimTracks...
  template <class T2>
  void fillAssociationMap( const edm::SimTrackContainer & firstCollection,
                           const std::vector<T2> & secondCollection,
                           std::map<const math::XYZTLorentzVectorD *, const T2 *> & matchesMap,
                           TH1 * hMinDeltaR = 0 ) {
    std::vector<math::XYZTLorentzVectorD> simTrackMomenta;
    edm::SimTrackContainer::const_iterator it = firstCollection.begin();
    for( ; it != firstCollection.end(); ++it ) {
      simTrackMomenta.push_back(it->momentum());
    }
    fillAssociationMap(simTrackMomenta, secondCollection, matchesMap, hMinDeltaR);
  }

  /// Specialization to cope with SimTracks...
  template <class T1>
  void fillAssociationMap( const std::vector<T1> & firstCollection,
                           const edm::SimTrackContainer & secondCollection,
                           std::map<const T1 *, const math::XYZTLorentzVectorD *> & matchesMap,
                           TH1 * hMinDeltaR = 0 ) {
    std::vector<const math::XYZTLorentzVectorD *> simTrackMomenta;
    edm::SimTrackContainer::const_iterator it = secondCollection.begin();
    for( ; it != secondCollection.end(); ++it ) {
      simTrackMomenta.push_back(&(it->momentum()));
    }
    fillAssociationMap(firstCollection, simTrackMomenta, matchesMap, hMinDeltaR);
  }

protected:
  double maxDeltaR_;
  bool removeMatches_;
  bool takeOppositeTracks_;
};

#endif // ASSOCIATORBYDELTAR_H
