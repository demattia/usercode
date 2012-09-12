#ifndef CANDIDATES_H
#define CANDIDATES_H

#include <TObject.h>
#include <TLorentzVector.h>
#include <string>
#include <vector>
#include "TreeCandidate.h"
#include "TreeLepton.h"

struct Candidates : public TObject
{
  Candidates()
  {}

  void clear()
  {
    run = -999;
    event = -999;
    numPV = -999;
    nvtx_m1 = -999;
    nvtx_0 = -999;
    nvtx_p1 = -999;
    nvtx_true = -999;
    MET = -999.;
    METPhi = -999.;
    genCTau_.clear();
    // triggers.clear();
    leptons_.clear();
    candidates_.clear();
  }

  // Event related variables
  Int_t run;
  Int_t event;
  Int_t numPV;
  Int_t nvtx_m1;
  Int_t nvtx_0;
  Int_t nvtx_p1;
  Float_t nvtx_true;
  Float_t MET;
  Float_t METPhi;
  // HLT paths triggered on in a particular event
  // std::vector<std::string> triggers;
  std::vector<double> genCTau_;

  std::vector<TreeLepton> leptons_;
  std::vector<TreeCandidate> candidates_;

  ClassDef(Candidates, 1)
};
ClassImp(Candidates)

#endif // CANDIDATES_H
