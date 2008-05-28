#ifndef LISTFASHIONATTRIBUTEDHISTO_H
#define LISTFASHIONATTRIBUTEDHISTO_H

#include <iostream>
#include "AnalysisExamples/AnalysisClasses/interface/FashionAttributedHisto.h"

namespace anaobj {

  template <class T>  
    class H160 : public FashionAttributedHisto<T> {

    public:
    H160( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 86, 0, 20, 0.4 ) { 
      std::cout << "entering in H160" << std::endl;      
    }
  };
  
  template <class T>
    class H200 : public FashionAttributedHisto<T> {
    public:
    H200( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 84, 3001, 20, 0.6 ) { }
  };
  
  template <class T>
    class H400 : public FashionAttributedHisto<T> {
    public:
    H400( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 211, 3009, 20, 0.8 ) { }
  };
  
  template <class T>
    class H800 : public FashionAttributedHisto<T> {
    public:
    H800( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 210, 3015, 20, 1.0 ) { }
  };
  
  template <class T>
    class ZZNjets : public FashionAttributedHisto<T> {
    public:
    ZZNjets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 65, 0, 23, 0.5 ) { }
  };
  
  template <class T>
    class ZZ0jets : public FashionAttributedHisto<T> {
    public:
    ZZ0jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 67, 3001, 23, 0.4 ) { }
  };
  
  template <class T>
    class ZZ1jets : public FashionAttributedHisto<T> {
    public:
    ZZ1jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 65, 3009, 23, 0.6 ) { }
  };
  
  template <class T>
    class ZZ2jets : public FashionAttributedHisto<T> {
    public:
    ZZ2jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 63, 3015, 23, 0.8 ) { }
  };
  
  template <class T>
    class ZbbNjets : public FashionAttributedHisto<T> {
    public:
    ZbbNjets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 215, 0, 22, 0.5 ) { }
  };

  template <class T>
    class Zbb0jets : public FashionAttributedHisto<T> {
    public:
    Zbb0jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 216, 3001, 22, 0.4 ) { }
  };

  template <class T>
    class Zbb1jets : public FashionAttributedHisto<T> {
    public:
    Zbb1jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 215, 3009, 22, 0.6 ) { }
  };

  template <class T>
    class Zbb2jets : public FashionAttributedHisto<T> {
    public:
    Zbb2jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 214, 3015, 22, 0.8 ) { }
  };

  template <class T>
    class WZNjets : public FashionAttributedHisto<T> {
    public:
    WZNjets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 98, 0, 29, 0.5 ) { }
  };

  template <class T>
    class WZ0jets : public FashionAttributedHisto<T> {
    public:
    WZ0jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 96, 3001, 29, 0.4 ) { }
  };

  template <class T>
    class WZ1jets : public FashionAttributedHisto<T> {
    public:
    WZ1jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 98, 3009, 29, 0.6 ) {
    }
    
  };

  template <class T>
    class WZ2jets : public FashionAttributedHisto<T> {
    public:
    WZ2jets( T* HISTO 
	     ) : FashionAttributedHisto<T>( HISTO, 100, 3015, 29, 0.8 ) { }
  };

  template <class T>
    class ttNjets : public FashionAttributedHisto<T> {
    public:
    ttNjets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 92, 0, 21, 0.5 ) { }
  };

  template <class T>
    class tt0jets : public FashionAttributedHisto<T> {
    public:
    tt0jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 91, 3001, 21, 0.4 ) { }
  };
  
  template <class T>
    class tt1jets : public FashionAttributedHisto<T> {
    public:
    tt1jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 92, 3009, 21, 0.6 ) { }
  };

  template <class T>
    class tt2jets : public FashionAttributedHisto<T> {
    public:
    tt2jets( T* HISTO ) : 
      FashionAttributedHisto<T>( HISTO, 93, 3015, 21, 0.8 ) { }
  };

}
#endif // LISTFASHIONATTRIBUTEDHISTO_H
