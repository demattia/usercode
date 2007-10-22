#ifndef L1PIXELTRIG_H
#define L1PIXELTRIG_H

#include "DivisiveAlgorithm.h"

template <class T1>
class L1PixelTrig {
 public:
  L1PixelTrig( int PJNUMBER ) {
    pjnumber_ = PJNUMBER;
    response_ = false;
  }
//  void Fill( vector<const PixelJet> );
  void Fill( const vector<T1> & vec_pj );
  bool Response() const {
    return response_;
  }
 private:
  int pjnumber_;
  bool response_;
  vector<const T1 *> vec_pj_ptr;
};

template <class T1>
//void L1PixelTrig::Fill ( vector<const PixelJet> & vec_pj ) {
void L1PixelTrig<T1>::Fill ( const vector<T1> & vec_pj ) {

  // Create a new vector<pointers to obj>, so that the passed vector
  // is not modified in any way
  typename vector<T1>::const_iterator vec_pj_it = vec_pj.begin();
  for ( ; vec_pj_it != vec_pj.end(); ++vec_pj_it ) {
    vec_pj_ptr.push_back( &*(vec_pj_it) );
  }

//  DivisiveAlgorithm<const PixelJet> div;
  DivisiveAlgorithm<T1> div;
  div.Fill( vec_pj_ptr );

  // Get the vector of verteces
//  auto_ptr<vector<Vertex<PixelJet> > > vertex_vec_aptr( div.VecVertex() );
  auto_ptr<vector<Vertex<T1> > > vertex_vec_aptr( div.VecVertex() );

  sort( vertex_vec_aptr->begin(), vertex_vec_aptr->end() );

  int c = 0;
  for( typename vector<Vertex<T1> >::const_iterator it_t = vertex_vec_aptr->begin(); it_t != vertex_vec_aptr->end(); ++it_t, ++c ) {
    cout << "trig pt["<<c<<"] = " << it_t->pt() << endl;
    cout << "trig eta["<<c<<"] = " << it_t->eta() << endl;
    cout << "trig phi["<<c<<"] = " << it_t->phi() << endl;
    cout << "trig z["<<c<<"] = " << it_t->z() << endl;
    cout << "trig number["<<c<<"] = " << it_t->number() << endl;
  }

  if ( (vertex_vec_aptr->back()).number() >= pjnumber_ ) {
    response_ = true;
  }
}

#endif // L1PIXELTRIG_H
