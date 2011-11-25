#ifndef TREETRACK_H
#define TREETRACK_H

#include <TObject.h>
#include <TLorentzVector.h>
#ifndef __CINT__
#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"
#endif

/**
 * Simple class used to save a track in a root tree. <br>
 * Includes the information on the genTrack.
 */

///// Simple structure used for genTracks (also the base for TreeTrackStruct that is used for reco tracks)
//struct TreeBaseTrackStruct
//{
//  TreeBaseTrackStruct() :
//    dxy(0.), dxyErr(0.), dz(0.), dzErr(0.),
//    charge(0),
//    vertexX(0.), vertexY(0.), vertexZ(0.)
//  {}

//  // virtual ~TreeBaseTrackStruct() {}

//#ifndef __CINT__
//  template <class T, class S>
//  TreeBaseTrackStruct(const T & inputTrack, const S & vertex) :
//    dxy(inputTrack.dxy()), dxyErr(inputTrack.dxyError()), dz(inputTrack.dz()), dzErr(inputTrack.dzError()),
//    charge(inputTrack.charge()),
//    vertexX(vertex.x()), vertexY(vertex.y()), vertexZ(vertex.z())
//  {
//    // Muon mass used
//    track.SetPtEtaPhiM(inputTrack.pt(), inputTrack.eta(), inputTrack.phi(), 0.105658);
//  }

//  template <class T>
//  TreeBaseTrackStruct(const SmartPropagatorWithIP::IP & ip, const int charge, const T & vertex) :
//    dxy(ip.dxyValue), dxyErr(ip.dxyError), dz(ip.dzValue), dzErr(ip.dzError),
//    charge(charge),
//    vertexX(vertex.x()), vertexY(vertex.y()), vertexZ(vertex.z())
//  {
//    // Muon mass used
//    track.SetPtEtaPhiM(ip.pt, ip.eta, ip.phi, 0.105658);
//  }
//#endif

//  double dxy, dxyErr, dz, dzErr;
//  int charge;
//  double vertexX, vertexY, vertexZ;
//  // Default initialized to (0,0,0,0)
//  TLorentzVector track;
//};

///// Simple structure used for reco tracks
//struct TreeTrackStruct : public TreeBaseTrackStruct
//{
//  TreeTrackStruct() :
//    referencePointRadius(0.), referencePointZ(0.),
//    nHits(0), nValidHits(0), nValidPlusInvalidHits(0),
//    innermostHitRadius(0.), innermostHitZ(0.),
//    chi2(0.), normalizedChi2(0.)
//  {}

//  // virtual ~TreeTrackStruct() {}

//#ifndef __CINT__
//  template <class T>
//  TreeTrackStruct(const T & inputTrack) :
//    TreeBaseTrackStruct(inputTrack, inputTrack.vertex()),
//    referencePointRadius(sqrt(pow(inputTrack.referencePoint().x(),2)+pow(inputTrack.referencePoint().y(),2))),
//    referencePointZ(inputTrack.referencePoint().z()),
//    nHits(inputTrack.recHitsSize()), nValidHits(inputTrack.found()), nValidPlusInvalidHits(inputTrack.found() + inputTrack.lost()),
//    innermostHitRadius(inputTrack.innerPosition().r()), innermostHitZ(inputTrack.innerPosition().z()),
//    chi2(inputTrack.chi2()), normalizedChi2(inputTrack.normalizedChi2())
//  {}

//  template <class T>
//  TreeTrackStruct(const T & inputTrack, const SmartPropagatorWithIP::IP & ip) :
//    TreeBaseTrackStruct(ip, inputTrack.charge(), inputTrack.vertex()),
//    referencePointRadius(sqrt(pow(inputTrack.referencePoint().x(),2)+pow(inputTrack.referencePoint().y(),2))),
//    referencePointZ(inputTrack.referencePoint().z()),
//    nHits(inputTrack.recHitsSize()), nValidHits(inputTrack.found()), nValidPlusInvalidHits(inputTrack.found() + inputTrack.lost()),
//    innermostHitRadius(inputTrack.innerPosition().r()), innermostHitZ(inputTrack.innerPosition().z()),
//    chi2(inputTrack.chi2()), normalizedChi2(inputTrack.normalizedChi2())
//  {}
//#endif

//  double referencePointRadius, referencePointZ;
//  int nHits, nValidHits, nValidPlusInvalidHits;
//  double innermostHitRadius, innermostHitZ;
//  double chi2, normalizedChi2;
//};

/// Class to be saved in the tree
class TreeTrack : public TObject
{
public:

//  TreeTrack() {}
//  // virtual ~TreeTrack() {}

////  template<class T, class S, class U>
////  TreeTrack(const T & track, const S genTrack, const U genVertex) :
////    recoTrack_(track), genTrack_(genTrack, genVertex)
////  {}

//  // Do not let CINT parse this or it will complain...
//#ifndef __CINT__
//  template <class T>
//  TreeTrack(const T & track) :
//    recoTrack(track)
//  {}

//  template <class T>
//  TreeTrack(const T & track, const SmartPropagatorWithIP::IP & ip) :
//    recoTrack(track, ip)
//  {}
//#endif

//  template <class T, class S, class U>
//  void setGen(const T & inputGenTrack, const S & inputGenVertex, const U & inputGenIP)
//  {
//    genTrack.dxy = inputGenIP.dxyValue;
//    genTrack.dxyErr = inputGenIP.dxyError;
//    genTrack.dz = inputGenIP.dzValue;
//    genTrack.dzErr = inputGenIP.dzError;
//    genTrack.charge = inputGenTrack.charge();
//    genTrack.track.SetPtEtaPhiM(inputGenTrack.pt(), inputGenTrack.eta(), inputGenTrack.phi(), 0.105658);
//    genTrack.vertexX = inputGenVertex.x();
//    genTrack.vertexY = inputGenVertex.y();
//    genTrack.vertexZ = inputGenVertex.z();
//  }

//protected:
//  TreeTrackStruct recoTrack;
//  TreeBaseTrackStruct genTrack;

  TreeTrack() {}

  double pt, ptError, eta, etaError, phi, phiError;
  int charge;
  double dxy, dxyError, dz, dzError;
  double vx, vy, vz;
  double chi2, normalizedChi2;
  double referencePointRadius, referencePointZ;
  int nHits, nValidHits, nValidPlusInvalidHits;
  double innermostHitRadius, innermostHitZ;

  double genPt, genEta, genPhi;
  int genCharge;
  double genDxy, genDxyError, genDz, genDzError;
  double genVx, genVy, genVz;

  ClassDef(TreeTrack, 1)
};
ClassImp(TreeTrack)

#endif // TREETRACK_H
