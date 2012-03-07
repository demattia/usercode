#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <map>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Analysis/TrackingEfficiencyFromCosmics/interface/Utils.h"
#include "Analysis/SmartPropagatorWithIP/interface/SmartPropagatorWithIP.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/TreeTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/**
  * Namespace containing some utility functions. They are declared static to limit their scope to this file and avoid linking errors
  * when including it in multiple plugins.
  */
namespace utils
{
  TH1F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name,
                       const TString & nameForTitle, const TString axisTitle, const TString & unit,
                       const int bins, const double & min, const double & max)
  {
    TString underscore("");
    if( nameForTitle != "" ) underscore = "_";

    TString spaceAndUnit("");
    if( unit != "" ) spaceAndUnit = " " + unit;

    TH1F * histo = fileService->make<TH1F>(name+underscore+nameForTitle, name+" "+nameForTitle, bins, min, max);
    histo->GetXaxis()->SetTitle(axisTitle+spaceAndUnit);
    std::stringstream ss;
    ss << histo->GetBinWidth(1);
    histo->GetYaxis()->SetTitle("entries / ( "+ss.str() + spaceAndUnit + " )");

    return histo;
  }

  TH2F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name, const TString & nameForTitle,
                       const TString xAxisTitle, const TString & xUnit,
                       const int xBins, const double & xMin, const double & xMax,
                       const TString yAxisTitle, const TString & yUnit,
                       const int yBins, const double & yMin, const double & yMax)
  {
    TString underscore("");
    if( nameForTitle != "" ) underscore = "_";

    TString spaceAndUnit("");
    if( xUnit != "" ) spaceAndUnit = " " + xUnit;

    TH2F * histo = fileService->make<TH2F>(name+underscore+nameForTitle, name+" "+nameForTitle, xBins, xMin, xMax, yBins, yMin, yMax);
    histo->GetXaxis()->SetTitle(xAxisTitle+spaceAndUnit);
    if( yUnit != "" ) spaceAndUnit = " " + yUnit;
    histo->GetYaxis()->SetTitle(yAxisTitle+spaceAndUnit);

    return histo;
  }

  /// Polynomial function fitted to the dxy_error vs Pt distribution
  double polynomial( const double & pt )
  {
    return ( 0.1+1.91364 - 0.0211496*pt + 0.0000906055*pt*pt - 0.000000130650*pt*pt*pt );
  }

  /// Returns the maximum dxy error value allowed for the cuts
  double dxyErrMax( const double & pt )
  {
    double dxyErrMax = 1.;
    if(pt < 200 ) dxyErrMax = polynomial(pt);
    else dxyErrMax = polynomial(200);
    return std::min(dxyErrMax, 1.);
  }

  void fillTrackToTreeTrack( TreeTrack & treeTrack, const reco::Track & track )
  {
    treeTrack.pt = track.pt();
    treeTrack.ptError = track.ptError();
    treeTrack.eta = track.eta();
    treeTrack.etaError = track.etaError();
    treeTrack.phi = track.phi();
    treeTrack.phiError = track.phiError();
    treeTrack.charge = track.charge();
    treeTrack.vx = track.vx();
    treeTrack.vy = track.vy();
    treeTrack.vz = track.vz();
    treeTrack.chi2 = track.chi2();
    treeTrack.normalizedChi2 = track.normalizedChi2();
    treeTrack.referencePointRadius = std::sqrt(std::pow(track.referencePoint().x(),2)+std::pow(track.referencePoint().y(),2));
    treeTrack.referencePointZ = track.referencePoint().z();
    treeTrack.nHits = track.recHitsSize();
    treeTrack.nValidHits = track.found();
    treeTrack.nValidPlusInvalidHits = track.found()+track.lost();
    treeTrack.innermostHitRadius = track.innerPosition().r();
    treeTrack.innermostHitZ = track.innerPosition().z();
  }

  void fillTrackToTreeTrack(TreeTrack & treeTrack, const reco::Track & track, const SmartPropagatorWithIP::IP & ip)
  {
    fillTrackToTreeTrack(treeTrack, track);
    treeTrack.pt = ip.pt;
    treeTrack.ptError = ip.ptError;
    treeTrack.eta = ip.eta;
    treeTrack.etaError = ip.etaError;
    treeTrack.phi = ip.phi;
    treeTrack.phiError = ip.phiError;
    treeTrack.dxy = ip.dxyValue;
    treeTrack.dxyError = ip.dxyError;
    treeTrack.dz = ip.dzValue;
    treeTrack.dzError = ip.dzError;
  }

  void fillGenToTreeTrack(TreeTrack & treeTrack, const reco::GenParticle & genTrack, const SmartPropagatorWithIP::IP & genIp)
  {
    treeTrack.genPt = genTrack.pt();
    treeTrack.genEta = genTrack.eta();
    treeTrack.genPhi = genTrack.phi();
    treeTrack.genCharge = genTrack.charge();
    treeTrack.genDxy = genIp.dxyValue;
    treeTrack.genDxyError = genIp.dxyError;
    treeTrack.genDz = genIp.dzValue;
    treeTrack.genDzError = genIp.dzError;
    treeTrack.genVx = genTrack.vx();
    treeTrack.genVy = genTrack.vy();
    treeTrack.genVz = genTrack.vz();
  }
}

#endif
