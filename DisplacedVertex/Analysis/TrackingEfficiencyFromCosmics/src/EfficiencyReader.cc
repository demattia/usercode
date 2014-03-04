// -*- C++ -*-
//
// Package:    EfficiencyReader
// Class:      EfficiencyReader
//
/**\class EfficiencyReader EfficiencyReader.cc Analysis/EfficiencyReader/src/EfficiencyReader.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,42 R-23,
//         Created:  Mon Jul 4 18:38:0 CEST 2011
// $Id: EfficiencyReader.cc,v 1.1 2011/07/04 17:01:11 demattia Exp $
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

#include <boost/foreach.hpp>

//
// class declaration
//

class EfficiencyReader : public edm::EDAnalyzer {
public:
  explicit EfficiencyReader(const edm::ParameterSet&);
  ~EfficiencyReader();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  // double maxDeltaR_;
  std::auto_ptr<Efficiency> efficiency_;
  std::string inputFileName_;
};

EfficiencyReader::EfficiencyReader(const edm::ParameterSet& iConfig) :
  inputFileName_(iConfig.getParameter<std::string>("InputFileName"))
{
  efficiency_.reset(new Efficiency);

  EfficiencyTree tree;
  tree.readTree(inputFileName_, &*efficiency_);

  unsigned int S = efficiency_->getLinearSize();
  for( unsigned int i=0; i<S; ++i ) {
    std::cout << "reco eff["<<i<<"] = " << efficiency_->getEff(i) << " +/- " << efficiency_->getEffError(i) << std::endl;
  }
}

EfficiencyReader::~EfficiencyReader() {}

void EfficiencyReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just before starting event loop  ------------
void EfficiencyReader::beginJob()
{
  edm::Service<TFileService> fileService;
}

// ------------ method called once each job just after ending the event loop  ------------
void EfficiencyReader::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void EfficiencyReader::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void EfficiencyReader::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void EfficiencyReader::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void EfficiencyReader::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EfficiencyReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyReader);
