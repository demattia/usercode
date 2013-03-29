#include "CounterReader.h"

void CounterReader::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & setup)
{
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  lumi.getByLabel("nEventsTotal", nEventsTotalCounter);
  // numProcessedEvents_ += nEventsTotalCounter->value;
  std::cout << "nEventsTotalCounter = " << nEventsTotalCounter->value << std::endl;
  numProcessedEventsHisto_->Fill(1, nEventsTotalCounter->value);

  // Number of events after filters
  edm::Handle<edm::MergeableCounter> nEventsAfterFilter;
  lumi.getByLabel("nEventsAfterFilter", nEventsAfterFilter);
  // numEventsPassingFilter_ += nEventsAfterFilter->value;
  std::cout << "nEventsAfterFilter = " << nEventsAfterFilter->value << std::endl;
  numEventsPassingFilterHisto_->Fill(1, nEventsAfterFilter->value);
}

CounterReader::CounterReader(const edm::ParameterSet& iConfig)
{
  // I want to make a histogram of number of tracks in an event
  edm::Service<TFileService> fs;
  numProcessedEventsHisto_ = fs->make<TH1D>("numProcessedEvents" , "Number of processed events" , 1, 0, 2 );
  numEventsPassingFilterHisto_ = fs->make<TH1D>("numEventsPassingFilter" , "Number of events passing the filters" , 1, 0, 2 );
}

CounterReader::~CounterReader()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called to for each event  ------------
void CounterReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just before starting event loop  ------------
void CounterReader::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void CounterReader::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(CounterReader);
