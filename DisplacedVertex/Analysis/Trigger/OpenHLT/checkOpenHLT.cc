#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

void checkOpenHLT()
{
   gROOT->Reset();

//   Connect file generated in $ROOTSYS/test
   TFile f("Event.root");

//   Create an histogram
   TH1F *hnseg = new TH1F("hnseg","Number of segments for selected tracks",5000,0,5000);

//   Start main loop on all events
   TClonesArray *tracks = 0;
   Event *event = new Event();   //object must be created before setting the branch address

   TBranch *bntrack = T.GetBranch("fNtrack");
   TBranch *branch  = T.GetBranch("event");
   branch->SetAddress(&event);
   Int_t nevent = T.GetEntries();
   Int_t nselected = 0;
   for (Int_t i=0;i<nevent;i++) {
      bntrack->GetEvent(i);
      if (event->GetNtrack() > 587)continue;
      nselected++;
      T.GetEvent(i);                  //read complete accepted event in memory
      hnseg->Fill(event->GetNseg());  //Fill histogram with number of segments
      tracks = event->GetTracks();    //get pointer to the TClonesArray object
      tracks->Clear();                //clear it
   }
//  Delete ClonesArray and histogram objects
   event->Finish();

//  Draw the histogram
   hnseg->Draw();
}
