#include <Analysis/RootTreeProducers/interface/RootTreeHandler.h>

RootTreeHandler::RootTreeHandler(const TString & fileName, const TString & treeName)
{
  file_ = new TFile(fileName, "RECREATE");
  tree_ = new TTree("T", treeName);
  tree_->Branch("event", &eventNumber_);
  tree_->Branch("run", &runNumber_);
}

// RootTreeHandler::~RootTreeHandler()
void RootTreeHandler::writeTree()
{
  file_->Write();
  file_->Close();
}

void RootTreeHandler::saveToTree( const unsigned int eventNumber,
                                  const unsigned int runNumber )
{
  eventNumber_ = eventNumber;
  runNumber_ = runNumber;
  tree_->Fill();
}
