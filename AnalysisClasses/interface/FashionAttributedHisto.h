#ifndef FASHIONATTRIBUTEDHISTO_H
#define FASHIONATTRIBUTEDHISTO_H

#include "TH1.h"

namespace anaobj {
  
  template <class T>
    class FashionAttributedHisto : 
    public T { 
    
    public:
    FashionAttributedHisto( const T*     HISTO, 
			    const int    COLOR       = 1,
			    const int    FILLSTYLE   = 1,
			    const int    MARKERSTYLE = 1,
			    const double MARKERSIZE  = 1.
			    ) : 
      T( *HISTO ) { 
      
      T::fLineColor   = COLOR;
      T::fFillColor   = COLOR;
      T::fFillStyle   = FILLSTYLE;
      T::fMarkerColor = COLOR;
      T::fMarkerStyle = MARKERSTYLE;
      T::fMarkerSize  = MARKERSIZE;

    }

  };

}
#endif // FASHIONATTRIBUTEDHISTO_H
