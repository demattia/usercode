#ifndef TREELEPTON_H
#define TREELEPTON_H

struct TreeLepton
{
  TreeLepton() :
    pt(-999.), eta(-999.), phi(-999.), charge(-999), d0(-999.), absD0Significance(-999.), iso(-999.),
    hasCaloMatch(-999), scEnergy(-999), scEt(-999), scEta(-999), scPhi(-999),
    isStandAlone(-999), isGlobalMuon(-999), isTrackerMuon(-999),
    triggerMatch(-999), index(-999),
    vx(-999.), vy(-999.), vz(-999.)
  {}
  // Pt, eta, phi etc. of original lepton tracks
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Int_t charge;
  Float_t d0;
  Float_t absD0Significance;
  Float_t iso;

  // Info on matched SC if any
  Int_t hasCaloMatch;
  Float_t scEnergy;
  Float_t scEt;
  Float_t scEta;
  Float_t scPhi;

  // What type of lepton
  Int_t isStandAlone;
  Int_t isGlobalMuon;
  Int_t isTrackerMuon;

  // Info on trigger match if any
  Int_t triggerMatch;

  Int_t index;

  // Reference point on track
  Float_t vx;
  Float_t vy;
  Float_t vz;

  ClassDef(TreeLepton, 1)
};
ClassImp(TreeLepton)

#endif
