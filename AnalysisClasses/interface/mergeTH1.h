/**

  This macro will add histograms multiplied by their cross-section
  from a list of root files and write them
  to a target root file. The target file is newly created and must not be
  identical to one of the source files.


  Author: Sven A. Schmidt, sven.schmidt@cern.ch
  Date:   13.2.2001

  Editing Author: Michael B. Anderson, mbanderson@hep.wisc.edu
  Date:  July 12, 2007

  This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
  which had a problem with directories more than one level deep.
  (see macro hadd_old.C for this previous implementation).

  The macro from Sven has been enhanced by
     Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
   to automatically add Trees (via a chain of trees).

  To use this macro, modify the file names in function hadd.

  NB: This macro is provided as a tutorial.
      Use $ROOTSYS/bin/hadd to merge many histogram files

 */


#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

template <class T1>
T1* mergeObj( TString               outputTH1name,
	       std::vector<T1*>   & inputVector,
	       std::vector<double> & xSecVector,  
	       bool                  scale = true
	       ) {



  std::vector<T1*>::const_iterator inputVector_itr  = inputVector.begin();
  std::vector<double>::const_iterator xSecVector_itr  = xSecVector.begin();
  T1 * outputTH1 = (T1*)(*inputVector_itr)->Clone();
  outputTH1->SetName(outputTH1name);
  if (scale) outputTH1->Scale(*xSecVector_itr);
  ++inputVector_itr;
  ++xSecVector_itr;
  for ( ; inputVector_itr != inputVector.end(); ++inputVector_itr,
	                                          ++xSecVector_itr) { 
    T1 * histo = (T1*)(*inputVector_itr)->Clone();
    if (scale) histo->Scale(*xSecVector_itr);
    outputTH1->Add(histo);
    delete histo;
  }
  return outputTH1;
}
