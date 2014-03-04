// -*- C++ -*-
//
// Package:    EfficiencyMerger
// Class:      EfficiencyMerger
//
/**\class EfficiencyMerger EfficiencyMerger.cc Analysis/EfficiencyMerger/src/EfficiencyMerger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,42 R-23,
//         Created:  Mon Aug 3 16:20:0 CEST 2011
// $Id: EfficiencyMerger.cc,v 1.1 2011/08/11 10:37:05 demattia Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Analysis/TrackingEfficiencyFromCosmics/interface/EfficiencyTree.h"

#include <vector>
#include <boost/foreach.hpp>

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

//
// class declaration
//

class EfficiencyMerger : public edm::EDAnalyzer {
public:
  explicit EfficiencyMerger(const edm::ParameterSet&);
  ~EfficiencyMerger();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void fillHistogram(const TString & name, const TString & title, const boost::shared_ptr<Efficiency> & eff );

  // ----------member data ---------------------------
  std::auto_ptr<Efficiency> efficiency_;
  std::vector<std::string> inputFileNames_;
};

EfficiencyMerger::EfficiencyMerger(const edm::ParameterSet& iConfig) :
  inputFileNames_(iConfig.getParameter<std::vector<std::string> >("InputFileNames"))
{
  if( inputFileNames_.size() < 2 ) {
    std::cout << "Not enought objects to merge. Object provied are " << inputFileNames_.size() << std::endl;
  }

  efficiency_.reset(new Efficiency);
  EfficiencyTree tree;
  std::vector<std::string>::const_iterator it = inputFileNames_.begin();
  tree.readTree(*it, &*efficiency_);
  ++it;
  for( ; it != inputFileNames_.end(); ++it ) {
    std::auto_ptr<Efficiency> eff;
    eff.reset(new Efficiency);
    tree.readTree(*it, &*eff);
    efficiency_->add(*eff);
  }

  tree.writeTree("MergedEfficiency.root", &*efficiency_);
}

EfficiencyMerger::~EfficiencyMerger() {}

void EfficiencyMerger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just before starting event loop  ------------
void EfficiencyMerger::beginJob()
{
  // edm::Service<TFileService> fileService;
}

// ------------ method called once each job just after ending the event loop  ------------
void EfficiencyMerger::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void EfficiencyMerger::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void EfficiencyMerger::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void EfficiencyMerger::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void EfficiencyMerger::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EfficiencyMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyMerger);
