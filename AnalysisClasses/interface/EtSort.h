#ifndef ETSORT_HH
#define ETSORT_HH

/**
 * Can be used to sort any class which has the method et()
 * The two overloaded methods correspond to value and pointer template type.
 *
 * Author M. De Mattia - 5/1/2008
 */

template <class T1>
struct EtSort {
  bool operator() ( const T1 first, const T1 second ) {
    return first.et() < second.et();
  }
  bool operator() ( const T1 * first, const T1 * second ) {
    return first->et() < second->et();
  }
};

#endif // ETSORT_HH
