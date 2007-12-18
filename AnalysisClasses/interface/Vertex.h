#ifndef VERTEX_H
#define VERTEX_H

/**
 * Simple Vertex template class
 * The Fill method requires a pointer to the template object. It
 * evaluates the weighted mean by pt of the eta, phi and z values.
 *
 * IMPORTANT:
 * the Close method needs to be called after the fill
 * is complete, or the eta, phi and z will still be multiplied by the
 * total pt. Close divides them, givin the correct weighted mean values.
 *
 * Overloads the operator< to allow sorting the verteces by Pt.
 *
 * Author M. De Mattia - 19/10/2007
 */

#include <iostream>
#include <vector>

template <class T1>
class Vertex {

 public:
  Vertex() {
    pt_ = 0;
    eta_ = 0;
    phi_ = 0;
    z_ = 0;
    objnumber_ = 0;
  }

  /// Fill the Vertex with the object, which has pt(), eta(), phi() and z() methods
  void Fill( const T1* obj_ptr );

  /// Close the Vertex. The correct eta, phi and z are evaluated
  /** eta, phi and z are evaluated as the weighted mean in pt of the corresponding
   *  corresponding values from the constituents.
   */
  void Close();

  double pt() const {
    return pt_;
  }
  double eta() const {
    return eta_;
  }
  double phi() const {
    return phi_;
  }
  double z() const {
    return z_;
  }
  int number() const {
    return objnumber_;
  }

  /// Return the vector of pointers to the objects used to build this vertex
  std::vector<const T1*> constituents() const {
    return vec_obj_ptr;
  }

  /// To allow ordering in Pt
  bool operator< ( const Vertex<T1> & comp_vertex ) const {
    return pt() < comp_vertex.pt();
  }

 private:
  double pt_;
  double eta_;
  double phi_;
  double z_;
  int objnumber_;
  std::vector<const T1*> vec_obj_ptr;
};

template <class T1>
void Vertex<T1>::Fill( const T1* obj_ptr ) {
  double pt_temp = obj_ptr->pt();

  ++objnumber_;

  pt_ += pt_temp;
  eta_ += obj_ptr->eta()*pt_temp;
  phi_ += obj_ptr->phi()*pt_temp;
  z_ += obj_ptr->z()*pt_temp;
  vec_obj_ptr.push_back( obj_ptr );
}

template <class T1>
void Vertex<T1>::Close() {
  if ( pt_ != 0 ) {
    eta_ = eta_/pt_;
    phi_ = phi_/pt_;
    z_ = z_/pt_;
  }
  else {
    std::cout << "Total pt of Vertex = 0, fill it before calling close or check the object passed to it" << std::endl;
  }
}

#endif // VERTEX_H
