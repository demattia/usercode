//#define DEBUG

#ifndef DIVISIVE_ALGORITHM_H
#define DIVISIVE_ALGORITHM_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

#include "Vertex.h"

using namespace std;

template <class T1>
struct Sort_Greater_Z {
  bool operator() ( const T1* a, const T1* b ) {
    return a->z() < b->z();
  }
  bool operator() ( const T1& a, const T1& b ) {
    return a.z() < b.z();
  }
};
// bool Sort_Greater_Z ( const T1* a, const T1* b ) {
//   return a->z() < b->z();
// }

template <class T1>
class DivisiveAlgorithm {
 public:
  DivisiveAlgorithm( double DZ = 0.4 ) {
    dz_ = DZ;
    // Create the auto_ptr to the vector of verteces
    vec_vertex_aptr_ = auto_ptr<vector<Vertex<T1> > > ( new vector<Vertex<T1> > );
    // Put the first, empty, vertex in the vector
    vec_vertex_aptr_->push_back( Vertex<T1>() );
  }
  void Fill( vector<const T1*> & vec_obj );
  auto_ptr<vector<Vertex<T1> > > VecVertex() {
    return vec_vertex_aptr_;
  }
 private:
  auto_ptr<vector<Vertex<T1> > > vec_vertex_aptr_;
//  vector<Vertex<T1> > vec_vertex_;
  double dz_;
};

template <class T1>
void DivisiveAlgorithm<T1>::Fill( vector<const T1*> & vec_obj ) {
  if ( vec_obj.size() > 0 ) {
    sort( vec_obj.begin(), vec_obj.end(), Sort_Greater_Z<T1>() );
    //  sort( vec_obj.begin(), vec_obj.end(), Sort_Greater_Z<T1> );

#ifdef DEBUG
    // Check that they are ordered in ascending z
    typename vector<const T1*>::const_iterator obj_it_test (vec_obj.begin()+1);
    int test_id = 0;
    for ( ; obj_it_test != vec_obj.end(); ++obj_it_test, ++test_id ) {
      //    cout << "pointer["<<test_id<<"] = " << *obj_it_test << endl;
      cout << "z["<<test_id<<"] = " << (*obj_it_test)->z() << endl;
    }
#endif

    // Create a vertex with the first obj
    (vec_vertex_aptr_->back()).Fill( *(vec_obj.begin()) );

    // Start looping on the remaining obj (loop from the second)
    if ( vec_obj.size() > 1 ) {
      typename vector<const T1*>::const_iterator obj_it (vec_obj.begin()+1);
      for ( ; obj_it != vec_obj.end(); ++obj_it ) {
        // If it is close to the last one add it to the vertex
        if ( fabs( (*obj_it)->z() - (*(obj_it-1))->z() ) < dz_  ) {
          (vec_vertex_aptr_->back()).Fill( *(obj_it) );
        }
        // Else, create a new vertex and add the obj to it
        else {
          // First close the last vertex. This divides the values by sumpt, so that
          // the weighted mean is completed.
          (vec_vertex_aptr_->back()).Close();

          vec_vertex_aptr_->push_back( Vertex<T1>() );
          (vec_vertex_aptr_->back()).Fill( *(obj_it) );
        }
      }
    }
    // Close the last vertex
    (vec_vertex_aptr_->back()).Close();
  }
}

#endif // DIVISIVE_ALGORITHM_H
