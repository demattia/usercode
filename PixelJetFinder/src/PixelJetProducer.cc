// #define DEBUG
// #define SIM

#include <memory>
//#include <string>
//#include <iostream>
//#include <fstream>

// Definitions 
// -----------
#include "AnalysisExamples/PixelJetFinder/interface/PixelJetProducer.h"

using namespace std;

// Constructor
// -----------
PixelJetProducer::PixelJetProducer(edm::ParameterSet const& conf) : 
  conf_(conf), 
  filename_(conf.getParameter<std::string>("fileName")),
  pixeljet_(conf.getParameter<std::string>("pixeljet")),
  eta_cut_(conf_.getParameter<double>( "Eta_cut" )),
  ConeR_cut_(conf_.getParameter<double>( "ConeR_cut" )),
  NumTk_cut_(conf_.getParameter<int>( "Minimum_PixelTracks" )) {

  produces<PixelJetCollection>( pixeljet_ );
  Eta_cut_Eff_ = 0.;
  eventcounter_ = 0;
}


// BeginJob
// --------
void PixelJetProducer::beginJob(const edm::EventSetup& es) {
}

// Virtual destructor needed
// -------------------------

PixelJetProducer::~PixelJetProducer() {
}

// Binary predicate to sort the vector<pair<const reco::Track*, unsigned int> >
// this will sort the pixeltracks in decreasing pt
bool PixelJetProducer_PtSort( const std::pair<const reco::Track*, unsigned int> & pair_1, const std::pair<const reco::Track*, unsigned int> & pair_2 ) {
  return (pair_1).first->pt() > (pair_2).first->pt();
}

// Producer: Function that gets called by framework every event
// ------------------------------------------------------------
void PixelJetProducer::produce(edm::Event& e, const edm::EventSetup& es) {

  using namespace edm;

  // Vector to store the AnalyzedCluster structs:
  std::auto_ptr<PixelJetCollection> v_pj_ptr(new PixelJetCollection);

  // TrackCollection
  // -------------------
  edm::Handle<reco::TrackCollection> trackCollection;
  e.getByLabel( conf_.getParameter<std::string>( "PixelTracksLabel" ), trackCollection);

#ifdef SIM
  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  reco::RecoToSimCollection myRecoToSimCollection;
  if ( SIM_ ) {
    // RECO/MC ASSOCIATION:
    e.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
    // Associator setup, change per event
    es.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
    associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();
    //Associate reco to sim tracks
    myRecoToSimCollection = associatorByHits->associateRecoToSim(trackCollection, TPCollectionH, &e);
  }
#endif

  // TrackCollection
  // ---------------
  const reco::TrackCollection *tracks=trackCollection.product();
  // Take the number of tracks in the event
  // --------------------------------------
  unsigned int numberoftracks = tracks->size();


  // Initializations per event
  // -------------------------
  
  eventcounter_++;



#ifdef DEBUG
  std::cout << "Event number " << eventcounter_ << std::endl;
#endif

  // Perform track study
  // -------------------
  if(numberoftracks>0){



    double PI = 3.141593;

    // First fill a vector of pointers to PixelTracks,
    // this vector will be looped on and will be ordered in decreasing Pt.
    // Inside the loop to build the PixelJet remove the used tracks. This way
    // the loop will become shorter with each iteration.

    vector<pair<const reco::Track*,unsigned int> > v_ptrpixeltracks;

    // Loop on track collection
    // ------------------------
    unsigned int TrackNum = 0;
    unsigned int TrackNum_Eta_cut = 0;
    reco::TrackCollection::const_iterator trackIter;
    for( trackIter=trackCollection->begin(); trackIter!=trackCollection->end(); ++trackIter ){

      // Store the pointers to the pixel tracks and their indeces in a vector.
      // The indeces are needed to create the Refs which need the index in the original collection.
      // Select only barrel tracks
      if ( fabs( trackIter->eta() ) < eta_cut_ ) {
	v_ptrpixeltracks.push_back( std::make_pair( &(*trackIter), TrackNum ) );
	++TrackNum_Eta_cut;
      }
      ++TrackNum;
    }

    Eta_cut_Eff_ += float(TrackNum_Eta_cut)/float(TrackNum);

    // Sort the vector in decreasing pt
    // This way the higher pt tracks will be used as seeds for the pixeljets
    // --------------------------------
    sort( v_ptrpixeltracks.begin(), v_ptrpixeltracks.end(), PixelJetProducer_PtSort );

#ifdef DEBUG
    int precount = 0;
    vector<pair<const reco::Track*, unsigned int> >::iterator preseedIter = v_ptrpixeltracks.begin();
    for( ; preseedIter != v_ptrpixeltracks.end(); ++preseedIter ) {
      ++precount;
      std::cout << "preloop: pt["<<precount<<"] = " << (*preseedIter).first->pt() << std::endl;
    }

    int postcount = 0;
#endif

    // Now loop on this vector and if the pixeltrack is assigned to a PixelJet, remove it from the vector

    // Pixel Jet construction
    // ----------------------

    bool erase = false;
    vector<pair<const reco::Track*, unsigned int> >::iterator seedIter = v_ptrpixeltracks.begin();
    for( ; seedIter != v_ptrpixeltracks.end(); ++seedIter ) {

//      double seedPt = (*seedIter).first->pt();
      double seedEta = (*seedIter).first->eta();
      double seedPhi = (*seedIter).first->phi();

      // construct an empty pixel jet
      PixelJet * PJ = 0;
      std::auto_ptr<PixelJet> PixelJet_auto_ptr; 

#ifdef DEBUG
      ++postcount;
      std::cout << std::endl;
      std::cout << "post pt["<<postcount<<"] = " << (*seedIter).first->pt() << std::endl;

      int candcount = 0;
#endif

      // Iterate on the remaining pixeltracks
      vector<pair<const reco::Track*, unsigned int> >::iterator pixeltkIter = seedIter+1;
      for( ; pixeltkIter != v_ptrpixeltracks.end(); ) {

	// Initialize the erase value
	erase = false;

#ifdef DEBUG
	++candcount;
	std::cout << "cand pt["<<candcount<<"] = " << (*pixeltkIter).first->pt() << std::endl;
#endif

//	double candPt = (*pixeltkIter).first->pt();
	double candEta = (*pixeltkIter).first->eta();
	double candPhi = (*pixeltkIter).first->phi();

	// Evaluate deltaR
	double deltaPhi = PI - fabs( fabs( seedPhi - candPhi ) - PI );
	double deltaEta = seedEta - candEta;
	double deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

	// If a compatible track is found
	if ( deltaR < ConeR_cut_ ) {

#ifdef DEBUG
	  std::cout << "erase this track["<<candcount<<"] with pt = " << (*pixeltkIter).first->pt() << std::endl;
#endif

	  // If there is no Pixel Jet construct a new one
	  if ( PJ == 0 ) {
	    PJ = new PixelJet;
	    // Put it in the auto_ptr;
	    PixelJet_auto_ptr.reset( PJ );
	    reco::TrackRef seedpixeltrackref=reco::TrackRef( trackCollection, (*seedIter).second );
	    PJ->PixelTrackRef( seedpixeltrackref );

#ifdef DEBUG
	    // seedpixeltrack
	    std::cout << "seedpixeltrackref->pt() = " << seedpixeltrackref->pt() << std::endl;
	    std::cout << "seedpixeltrackref->eta() = " << seedpixeltrackref->eta() << std::endl;
	    std::cout << "seedpixeltrackref->phi() = " << seedpixeltrackref->phi() << std::endl;
	    // PJ
	    std::cout << "PJ->pt() = " << PJ->pt() << std::endl;
	    std::cout << "PJ->eta() = " << PJ->eta() << std::endl;
	    std::cout << "PJ->phi() = " << PJ->phi() << std::endl;
#endif
	  }

	  // Ref to PixelTrack
	  reco::TrackRef pixeltrackref=reco::TrackRef( trackCollection, (*pixeltkIter).second );
	  PJ->PixelTrackRef( pixeltrackref );

#ifdef DEBUG
	    std::cout << "pixeltrackref->pt() = " << pixeltrackref->pt() << std::endl;
	    std::cout << "pixeltrackref->eta() = " << pixeltrackref->eta() << std::endl;
	    std::cout << "pixeltrackref->phi() = " << pixeltrackref->phi() << std::endl;
	    // PJ
	    std::cout << "PJ->pt() = " << PJ->pt() << std::endl;
	    std::cout << "PJ->eta() = " << PJ->eta() << std::endl;
	    std::cout << "PJ->phi() = " << PJ->phi() << std::endl;
#endif

	  // Remove this track from the vector
	  erase = true;
	}

	// This is to avoid to get an invalid iterator when the vector is modifyed.
	// The erase method returns the first valid iterator after the deleted one or end().
	// The iterator is only incremented if the vector is not modifyed.
	if ( erase ) {
	  pixeltkIter = v_ptrpixeltracks.erase(pixeltkIter);
	}
	else {
	  ++pixeltkIter;
	}

      } // end loop on remaining pixeltracks


//      ++trackcounter;  // only count processed tracks (not tracks->size())

//       momentum     = trackIter->p();
//       pt           = trackIter->pt();
//       charge       = trackIter->charge();
//       eta          = trackIter->eta();
//       phi          = trackIter->phi();
//       hitspertrack = trackIter->recHitsSize();
//       normchi2     = trackIter->normalizedChi2();
//       chi2         = trackIter->chi2();
//       ndof         = trackIter->ndof();
//       d0           = trackIter->d0();
//       vx           = trackIter->vx();
//       vy           = trackIter->vy();
//       vz           = trackIter->vz();
//       outerPt      = trackIter->outerPt();
//       found        = trackIter->found();
//       lost         = trackIter->lost();



#ifdef SIM
      if ( SIM_ ) {
	// PID (associate a TrackingParticle to the Track):
	std::vector<std::pair<TrackingParticleRef, double> > myTrackingParticleRefMap1 = myRecoToSimCollection[trackref];
	if ( myTrackingParticleRefMap1.size() != 0 ) {
	  std::vector<std::pair<TrackingParticleRef, double> >::const_iterator myTrackingParticleIt1 = myTrackingParticleRefMap1.begin();
	  TrackingParticleRef myTrackingParticle1 = myTrackingParticleIt1->first;
#ifdef DEBUG
	  cout << "~~~ Reco Track " << trackref.index() << " pT: " << trackref->pt() 
	       <<  " matched to " << myTrackingParticleRefMap1.size() << " MC Tracks" << std::endl;
#endif

	  // ---------------------------------------------------------------------------------------------------------
	  // Taking the iterator and overwriting the value. The track will have the values of the last associated
	  // mctrack. Not a problem in TIF simulation since are single muon events and usually there is only one match
	  // but should be changed
	  // ---------------------------------------------------------------------------------------------------------
	  for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = myTrackingParticleRefMap1.begin(); 
	       it != myTrackingParticleRefMap1.end(); ++it) {
	    TrackingParticleRef tpr = it->first;

	    trackingparticleP = tpr->p();
	    trackingparticlePt = tpr->pt();
	    trackingparticleEta = tpr->eta();
	    trackingparticlePhi = tpr->phi();
	    trackingparticleAssocChi2 = it->second;
	    trackingparticlepdgId = tpr->pdgId();

	    // Loop on simtracks
	    TrackingParticle::g4t_iterator g4T = tpr->g4Track_begin();
	    if ( g4T != tpr->g4Track_end() ) {
	      trackingparticletrackId = g4T->trackId();
	    }

	    // 	  for (TrackingParticle::g4t_iterator g4T = tpr->g4Track_begin();
	    // 	       g4T !=  tpr->g4Track_end(); ++g4T) {
	    // 	    std::cout << "simtrackid = " << g4T->trackId() << std::endl;
	    // 	  }

	  }
	}
      } // end if SIM_
#endif


      // Put the PixelJet in the collection
      // ----------------------------------
      // ----------------------------------

      // If a pixeljet was created
      if ( PJ != 0 ) {

	// First close the PJ, this evaluates the correct pt, eta and phi (divides by the number of pixeltracks)
	PJ->Close();

	// If the PJ has enough tracks, save it
	if ( PJ->tkNum() >= NumTk_cut_ ) {
	  // Put the PixelJet in the collection
	  v_pj_ptr->push_back( *PJ );
	}
      }

    } // Loop on track collection [trackIter]
  } // numberoftracks > 0
  
  // Fill track tree
  // ---------------

//  int numPJ = v_pj_ptr->size();

#ifdef DEBUG
  std::cout << "putting a vector of size = " << v_pj_ptr->size() << std::endl;
#endif

  e.put(v_pj_ptr, pixeljet_);

} // end produce


// EndJob
// ------
void PixelJetProducer::endJob(){

  std::cout << "Mean efficiency of Eta < "<<eta_cut_<< " cut on pixeltracks = " << Eta_cut_Eff_/float(eventcounter_) << std::endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelJetProducer);
