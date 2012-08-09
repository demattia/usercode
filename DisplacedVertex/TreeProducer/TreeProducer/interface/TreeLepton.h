#ifndef TREELEPTON_H
#define TREELEPTON_H

struct TreeLepton
{
  TreeLepton() :
    pt(-999.), eta(-999.), phi(-999.), charge(-999), d0(-999.), absD0Significance(-999.), iso(-999.),
    isStandAlone(-999), isGlobalMuon(-999), hasCaloMatch(-999), triggerMatch(-999), index(-999),
    vx(-999.), vy(-999.), vz(-999.)
  {}
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Int_t charge;
  Float_t d0;
  Float_t absD0Significance;
  Float_t iso;
  Int_t isStandAlone;
  Int_t isGlobalMuon;
  Int_t isTrackerMuon;
  Int_t hasCaloMatch;
  Int_t triggerMatch;
  Int_t index;
  Float_t vx;
  Float_t vy;
  Float_t vz;
  ClassDef(TreeLepton, 1)
};
ClassImp(TreeLepton)

#endif
