#ifndef ASSOCIATORET_HH
#define ASSOCIATORET_HH

/**
 *
 * A DeltaR associator
 *
 * Author M. De Mattia - 9/8/2007
 *
 * It receives two vectors of objects which must have the following methods:
 * eta(), phi().
 * It is an exclusive associator, but does not associate always the closest.
 * It loops on the first vector, and associates the closest element to it.
 * (This does not ensure that it is closer to one of the following elements
 * of the first vector).
 * It then removes this element from the second vector and iterates.
 * The first vector should be provided ordered (for example in Pt) to
 * justify this approach. This way the higher pt elements are compared
 * first and have higher priority on being associated.
 *
 * The association is done within a cone of radius R (which can be set).
 *
 * Returns a map with for key the pointer to the const elements of the first
 * vector and value those to the const elements of the second vector.
 *
 * (Possibly return also a map with (DeltaR, key of the other map) to be
 * able to recover the DeltaR).
 *
 * Modifications on 6/1/2008:
 * Added an Associator method receiving directly vectors of pointers.
 *
 * This method modifies the vectors of pointers passed to it.
 */

#include <memory>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>

template <class T1, class T2>
bool AssociatorEt_PtSort( const T1* first, const T2* second ) {

  return first->et() < second->et();
}

template <class T1, class T2>
class AssociatorEt {
public:
  AssociatorEt( double CONER_CUT ) {
    PI_ = 3.141593;
    ConeR_cut_ = CONER_CUT;
  }
  std::auto_ptr<std::map<const T1*, const T2*> > Associate( const std::vector<T1> & v_T1, const std::vector<T2> & v_T2 );
  std::auto_ptr<std::map<const T1*, const T2*> > Associate( std::vector<const T1*> & v_T1, const std::vector<T2> & v_T2 );
  std::auto_ptr<std::map<const T1*, const T2*> > Associate( std::vector<const T1*> & v_T1, std::vector<const T2*> & v_T2 );

 private:
  double PI_;
  double ConeR_cut_;
};

template <class T1, class T2>
std::auto_ptr<std::map<const T1*, const T2*> > AssociatorEt<T1, T2>::Associate( const std::vector<T1> & v_T1, 
                                                                                const std::vector<T2> & v_T2 ) {
  using namespace std;

  // Store vector of pointers if are not passed
  // ------------------------------------------
  vector<const T1*> v_T1_ptr;
  typename vector<T1>::const_iterator temp_it1 = v_T1.begin();
  for ( ; temp_it1 != v_T1.end(); ++temp_it1 ) {
    v_T1_ptr.push_back( &(*temp_it1) );
  }
  vector<const T2*> v_T2_ptr;
  typename vector<T2>::const_iterator temp_it2 = v_T2.begin();
  for ( ; temp_it2 != v_T2.end(); ++temp_it2 ) {
    v_T2_ptr.push_back( &(*temp_it2) );
  }
  // Call the method with the vectors of pointers
  return( Associate( v_T1_ptr, v_T2_ptr ) );
}

template <class T1, class T2>
std::auto_ptr<std::map<const T1*, const T2*> > AssociatorEt<T1, T2>::Associate( std::vector<const T1*> & v_T1_ptr, 
                                                                                const std::vector<T2> & v_T2 ) {
  using namespace std;

  // Store vector of pointers if are not passed
  // ------------------------------------------
  // Call the method with the vectors of pointers
  vector<const T2*> v_T2_ptr;
  typename vector<T2>::const_iterator temp_it2 = v_T2.begin();
  for ( ; temp_it2 != v_T2.end(); ++temp_it2 ) {
    v_T2_ptr.push_back( &(*temp_it2) );
  }
  return( Associate( v_T1_ptr, v_T2_ptr ) );
}

template <class T1, class T2>
std::auto_ptr<std::map<const T1*, const T2*> > AssociatorEt<T1, T2>::Associate( std::vector<const T1*> & v_T1_ptr, 
                                                                                std::vector<const T2*> & v_T2_ptr ) {

  using namespace std;

  sort( v_T1_ptr.begin(), v_T1_ptr.end(), AssociatorEt_PtSort<T1, T1> );
//  sort( v_T2_ptr.begin(), v_T2_ptr.end(), AssociatorEt_PtSort<T2, T2> );

  // Map to be returned
  auto_ptr<map<const T1*, const T2*> > AssocMap_ptr( new map<const T1*, const T2*> );

  typename vector<const T1*>::const_iterator it1 = v_T1_ptr.begin();
  typename vector<const T2*>::const_iterator it2;

  // Loop on the first vector
  for( ; it1 != v_T1_ptr.end(); ++it1 ) {
    double seedEta = (*it1)->eta();
    double seedPhi = (*it1)->phi();

    // Map to evaluate the closest to it1
    map<double, const T2*> temp_map;

    // Loop on the second vector and associate the closest
    it2 = v_T2_ptr.begin();
    for( ; it2 != v_T2_ptr.end(); ++it2 ) {
      double candEta = (*it2)->eta();
      double candPhi = (*it2)->phi();

      // Evaluate deltaR
      double deltaPhi = PI_ - fabs( fabs( seedPhi - candPhi ) - PI_ );
      double deltaEta = seedEta - candEta;
      double deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

      // If it is smaller then the association radius
      if ( deltaR < ConeR_cut_ ) {

	temp_map.insert( make_pair( deltaR, *it2 ) );

      }
    }

    // If it was associated
    if ( temp_map.size() != 0 ) {
      // Select the closest
      const T2* assoc_ptr = temp_map.begin()->second;
      // Write it in the map to be returned
      AssocMap_ptr->insert( make_pair( *it1, assoc_ptr ) );
      // and remove it from the v_T2.
      typename vector<const T2*>::iterator erase_it = find( v_T2_ptr.begin(), v_T2_ptr.end(), assoc_ptr );
      v_T2_ptr.erase( erase_it );
    }
  } // end loop on v_T1_ptr

  return AssocMap_ptr;

}

#endif // ASSOCIATOR_HH
