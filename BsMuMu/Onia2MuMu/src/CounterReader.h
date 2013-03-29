#ifndef COUNTERREADER_H
#define COUNTERREADER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//
// class decleration
//

class CounterReader : public edm::EDAnalyzer {
public:
  explicit CounterReader(const edm::ParameterSet&);
  ~CounterReader();

private:
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
  virtual void endJob();

  // ----------member data ---------------------------
  // unsigned int numProcessedEvents_;
  // unsigned int numEventsPassingFilter_;
  TH1D *numProcessedEventsHisto_;
  TH1D *numEventsPassingFilterHisto_;
};

#endif // COUNTERREADER_H
