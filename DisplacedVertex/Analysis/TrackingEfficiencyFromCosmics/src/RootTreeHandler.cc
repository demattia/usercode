#include <Analysis/TrackingEfficiencyFromCosmics/interface/RootTreeHandler.h>

RootTreeHandler::RootTreeHandler(const TString & fileName)
{
  file_ = new TFile(fileName, "RECREATE");
  tree_ = new TTree("T", "Muons");
  tree_->Branch("event", &eventNumber_);
  tree_->Branch("run", &runNumber_);
  tree_->Branch("tracks", "std::vector<TreeTrack>", &tracks_);
}

// RootTreeHandler::~RootTreeHandler()
void RootTreeHandler::writeTree()
{
  file_->Write();
  file_->Close();
}

void RootTreeHandler::saveToTree( const std::vector<TreeTrack> & tracks, const unsigned int eventNumber, const unsigned int runNumber )
{
  std::cout << "Filling tree with " << tracks.size() << " tracks" << std::endl;
  eventNumber_ = eventNumber;
  runNumber_ = runNumber;
  tracks_ = tracks;
  tree_->Fill();
}
