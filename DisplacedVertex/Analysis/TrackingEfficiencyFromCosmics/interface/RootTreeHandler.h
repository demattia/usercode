#ifndef ROOTTREEHANDLER_H
#define ROOTTREEHANDLER_H

#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include <Analysis/TrackingEfficiencyFromCosmics/interface/TreeTrack.h>
#include <TH1F.h>
#include <stdlib.h>

/**
 * This class can be used to save the muons (and gen muons if any) to a root tree. <br>
 * The writeTree method gets the name of the file to store the tree and the muons to save. <br>
 */

class RootTreeHandler
{
public:
  RootTreeHandler(const TString & fileName);
  // ~RootTreeHandler();

  void saveToTree( const std::vector<TreeTrack> & tracks, const unsigned int eventNumber, const unsigned int runNumber );
  void writeTree();

protected:
  TFile * file_;
  TTree * tree_;
  unsigned int eventNumber_;
  unsigned int runNumber_;
  std::vector<TreeTrack> tracks_;
};

#endif // ROOTTREEHANDLER_H
