// -*- C++ -*-
//
// Package:    TrackingEfficiencyFromCosmics
// Class:      TrackingEfficiencyFromCosmics
// 
/**\class TrackingEfficiencyFromCosmics TrackingEfficiencyFromCosmics.cc Analysis/TrackingEfficiencyFromCosmics/src/TrackingEfficiencyFromCosmics.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,40 3-B32,+41227671551,
//         Created:  Wed May 25 16:44:02 CEST 2011
// $Id: TrackingEfficiencyFromCosmics.cc,v 1.1 2011/05/26 10:50:27 demattia Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class TrackingEfficiencyFromCosmics : public edm::EDAnalyzer {
   public:
      explicit TrackingEfficiencyFromCosmics(const edm::ParameterSet&);
      ~TrackingEfficiencyFromCosmics();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackingEfficiencyFromCosmics::TrackingEfficiencyFromCosmics(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


TrackingEfficiencyFromCosmics::~TrackingEfficiencyFromCosmics()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackingEfficiencyFromCosmics::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);

  if( tracks->size() > 0 ) {
    LogInfo("Demo") << "number of tracks "<<tracks->size();

    unsigned int trackNumber = 0;
    std::vector<reco::Track>::const_iterator ittrk = tracks->begin();
    for( ; ittrk != tracks->end(); ++ittrk, ++trackNumber ) {
      LogInfo("Demo") << "track["<<trackNumber<<"] eta = " << ittrk->eta();
      LogInfo("Demo") << "track["<<trackNumber<<"] number of hits = " << ittrk->recHitsSize() << std::endl;
    }
  }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackingEfficiencyFromCosmics::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackingEfficiencyFromCosmics::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TrackingEfficiencyFromCosmics::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackingEfficiencyFromCosmics::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackingEfficiencyFromCosmics::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackingEfficiencyFromCosmics::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackingEfficiencyFromCosmics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingEfficiencyFromCosmics);
