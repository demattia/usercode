// -*- C++ -*-
//
// Package:    IsoTrackPrefilter
// Class:      IsoTrackPrefilter
// 
/**\class IsoTrackPrefilter IsoTrackPrefilter.cc HarderAnalysis/IsoTrackPrefilter/src/IsoTrackPrefilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kristian Harder
//         Created:  Tue Feb 22 21:16:02 GMT 2011
// $Id: IsoTrackPrefilter.cc,v 1.9 2011/07/15 21:17:20 harder Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DisplacedLeptons/Samples/interface/GenEventProperties.h"
#include "HarderAnalysisTools/Histograms/interface/HistMap.h"

//
// class declaration
//

class IsoTrackPrefilter : public edm::EDFilter {
public:
  explicit IsoTrackPrefilter(const edm::ParameterSet&);
  ~IsoTrackPrefilter();
  
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  
  // ----------member data ---------------------------

  edm::InputTag trackTag_;
  edm::InputTag muonTag_;
  edm::InputTag generatorTag_;
  edm::InputTag pileupTag_;
  edm::InputTag genEventInfoTag_;

  double minPt_;
  double coneSize_;
  double isolationThreshold_;

  int signalPDGId_;

  unsigned numEvents_;
  unsigned numEvents_pass_;

  HistMap* filterHistos_;

};

//
// constructors and destructor
//
IsoTrackPrefilter::IsoTrackPrefilter(const edm::ParameterSet& iConfig) :
  trackTag_(iConfig.getParameter<edm::InputTag>("trackSrc")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonSrc")),
  generatorTag_(iConfig.getParameter<edm::InputTag>("generatorSrc")),
  pileupTag_(edm::InputTag("addPileupInfo")),
  genEventInfoTag_(edm::InputTag("generator")),
  minPt_(iConfig.getParameter<double>("minPt")),
  coneSize_(iConfig.getParameter<double>("coneSize")),
  isolationThreshold_(iConfig.getParameter<double>("isolationThreshold")),
  signalPDGId_(iConfig.getParameter<int>("signalPDGId"))
{
  numEvents_=0;
  numEvents_pass_=0;

  filterHistos_ = new HistMap();

  std::cout << "IsoTrackPrefilter: p_t cut is "
	    << minPt_ << std::endl;
  std::cout << "IsoTrackPrefilter: isolation cone size: "
	    << coneSize_ << std::endl;
  std::cout << "IsoTrackPrefilter: isolation threshold: "
	    << isolationThreshold_ << std::endl;
  std::cout << "IsoTrackPrefilter: signal PDG is "
	    << signalPDGId_ << std::endl;
}


IsoTrackPrefilter::~IsoTrackPrefilter()
{
  std::cout << "IsoTrackPrefilter: num events = "
	    << numEvents_ << std::endl;
  std::cout << "IsoTrackPrefilter: num events passed = " 
	    << numEvents_pass_ << std::endl;
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
IsoTrackPrefilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // access track collection
   Handle<RefVector<reco::TrackCollection> > tracks;
   iEvent.getByLabel(trackTag_,tracks);

   // sort tracks into high p_t and low p_t categories
   std::vector<const reco::Track*> highPtTracks,lowPtTracks;
   highPtTracks.clear();
   lowPtTracks.clear();
   for (unsigned i=0; i<tracks->size(); i++) {
     Ref<reco::TrackCollection> trk = (*tracks)[i];
     if (trk->pt()>minPt_) highPtTracks.push_back(&(*trk)); else lowPtTracks.push_back(&(*trk));
   }

   // add high p_t standalone muons to list of high_pt tracks
   Handle<reco::TrackCollection> muons;
   iEvent.getByLabel(muonTag_,muons);
   for (reco::TrackCollection::const_iterator trk=muons->begin(); trk!=muons->end(); trk++) {
     if (trk->pt()>minPt_) highPtTracks.push_back(&(*trk));
   }
   unsigned nHighPt=highPtTracks.size();
   bool pass_nHighPt=(nHighPt>=2);

   // calculate tracker isolation for high p_t tracks, counting only low p_t tracks
   std::vector<double> isolation;
   unsigned nIsolHighPt=0;
   isolation.clear();
   for (unsigned iHighPt=0; iHighPt<highPtTracks.size(); iHighPt++) {
     double sumPt=0;
     for (unsigned iLowPt=0; iLowPt<lowPtTracks.size(); iLowPt++) {
       double dR=deltaR<const reco::Track>(*(highPtTracks[iHighPt]),*(lowPtTracks[iLowPt]));
       if (dR<coneSize_) sumPt+=lowPtTracks[iLowPt]->pt();
     }
     isolation.push_back(sumPt);
     if (sumPt<isolationThreshold_) ++nIsolHighPt;
   }
   bool pass_nIsolHighPt=(nIsolHighPt>=2);

   ++numEvents_;
   if (pass_nIsolHighPt) ++numEvents_pass_;

   // histograms to evaluate usefulness of selection criteria
   filterHistos_->fill("numHighPtTracks",nHighPt,1,
		       "number of high p_t tracks in event",20,0,20);
   filterHistos_->fill("numIsolHighPtTracks",nIsolHighPt,1,
		       "number of isolated high p_t tracks in event",20,0,20);
   for (unsigned i=0; i<isolation.size(); i++) {
     filterHistos_->fill("isolationDistr",isolation[i],1,
			 "sum p_t in cone around high p_t track,"
			 " excluding other high p_t tracks",
			 100,0,10);
     filterHistos_->fill("isolationDistrWide",isolation[i],1,
			 "sum p_t in cone around high p_t track,"
			 " excluding other high p_t tracks",
			 100,0,100);
   }

   GenEventProperties genProp(iEvent,signalPDGId_,generatorTag_,
			      pileupTag_,genEventInfoTag_);
   
   // pile-up information
   filterHistos_->fill("pileup_1bx",genProp.numPileupInTime(),1,
		       "in-time pile-up distribution",30,0,30);
   filterHistos_->fill("pileup_3bx",genProp.aveNumPileup3BX(),1,
		       "3BX average pile-up distribution",30,0,30);
   if (pass_nIsolHighPt) {
     filterHistos_->fill("pileup_1bx_pass",genProp.numPileupInTime(),1,
			 "in-time pile-up distribution",30,0,30);
     filterHistos_->fill("pileup_3bx_pass",genProp.aveNumPileup3BX(),1,
			 "3BX average pile-up distribution",30,0,30);
   }
   
   // PDF information
   filterHistos_->fill("pdf_id1",genProp.pdfId(0),1,"PDF ID 1",10,0,10);
   filterHistos_->fill("pdf_id2",genProp.pdfId(1),1,"PDF ID 2",10,0,10);
   filterHistos_->fill("pdf_x1",genProp.pdfX(0),1,"PDF x 1",100,0,1);
   filterHistos_->fill("pdf_x2",genProp.pdfX(1),1,"PDF x 2",100,0,1);
   filterHistos_->fill("pdf_pdfval1",genProp.pdfXPDF(0),1,"PDF value 1",100,0,2);
   filterHistos_->fill("pdf_pdfval2",genProp.pdfXPDF(1),1,"PDF value 2",100,0,2);
   filterHistos_->fill("pdf_scale",genProp.pdfScale(),1,"PDF scale",100,0,10000);
   
   // generator information (if any)
   // to evaluate efficiency for specific signal channels
   int n_elec=0, i_elec=0;
   int n_muon=0, i_muon=0;
   int n_tau=0, i_tau=0;
   for (unsigned i=0; i<genProp.numDecays(); i++) {
     if (genProp.decayMode(i)==11) { ++n_elec; i_elec=i; }
     if (genProp.decayMode(i)==13) { ++n_muon; i_muon=i; }
     if (genProp.decayMode(i)==15) { ++n_tau;  i_tau=i;  }
   }

   // histograms: number of events and event type before filter
   filterHistos_->fill("numSignal",genProp.numDecays(),1,
		       "number of signal particles in event",
		       10,0,10);
   filterHistos_->fill("numSignalE",n_elec,1,
		       "number of signal particles in event decaying to ee",
		       10,0,10);
   filterHistos_->fill("numSignalMu",n_muon,1,
		       "number of signal particles in event decaying to mumu",
		       10,0,10);
   filterHistos_->fill("numSignalTau",n_tau,1,
		       "number of signal particles in event decaying to tautau",
		       10,0,10);
   if (genProp.numDecays()==2) {
     filterHistos_->fill("numSignalEMu",n_elec,n_muon,1,
			 "number of signal particles in event decaying to"
			 " electrons (x) and muons (y)",10,0,10,10,0,10);
   }
   if (n_elec==1) {
     filterHistos_->fill("decayLength2D_1elec",genProp.decayLength2D(i_elec),1,
			 "2d decay length, electron channel, 1 electron decay",
			 100,0,100);
   } else if (n_elec==2) {
     for (unsigned i=0; i<2; i++) 
       filterHistos_->fill("decayLength2D_2elec",genProp.decayLength2D(i),1,
			   "2d decay length, electron channel, 2 electron decays",
			   100,0,100);
   }
   if (n_muon==1) {
     filterHistos_->fill("decayLength2D_1muon",genProp.decayLength2D(i_muon),1,
			 "2d decay length, muon channel, 1 muon decay",
			 100,0,100);
   } else if (n_muon==2) {
     for (unsigned i=0; i<2; i++)
       filterHistos_->fill("decayLength2D_2muon",genProp.decayLength2D(i),1,
			   "2d decay length, muon channel, 2 muon decays",
			   100,0,100);
   }
   if (n_tau==1) {
     filterHistos_->fill("decayLength2D_1tau",genProp.decayLength2D(i_tau),1,
			 "2d decay length, tau channel, 1 tau decay",
			 100,0,100);
   } else if (n_tau==2) {
     for (unsigned i=0; i<2; i++)
       filterHistos_->fill("decayLength2D_2tau",genProp.decayLength2D(i),1,
			   "2d decay length, tau channel, 2 tau decays",
			   100,0,100);
   }
   for (unsigned i=0; i<genProp.numDecays(); i++) {
     filterHistos_->fill("ctau",genProp.ctau(i),1,"ctau distribution",100,0,100);
   }
   
   // histograms: number of events and event type after track requirement
   if (pass_nHighPt) {
     filterHistos_->fill("numSignalPassNTrk",genProp.numDecays(),1,
			 "number of signal particles in event,"
			 " pass nHighPt requirement in filter",10,0,10);
     
     filterHistos_->fill("numSignalEPassNTrk",n_elec,1,
			 "number of signal particles in event decaying"
			 " to ee, pass nHighPt requirement in filter",
			 10,0,10);
     filterHistos_->fill("numSignalMuPassNTrk",n_muon,1,
			 "number of signal particles in event decaying"
			 " to mumu, pass nHighPt requirement in filter",
			 10,0,10);
     filterHistos_->fill("numSignalTauPassNTrk",n_tau,1,
			 "number of signal particles in event decaying"
			 " to tautau, pass nHighPt requirement in filter",
			 10,0,10);
     if (genProp.numDecays()==2) {
       filterHistos_->fill("numSignalEMuPassNTrk",n_elec,n_muon,1,
			   "number of signal particles in event"
			   " decaying to electrons (x) and muons (y),"
			   " pass nHighPt requirement in filter",
			   10,0,10,10,0,10);
     }
   }

   // histograms: number of events and event type after isolation requirement
   if (pass_nIsolHighPt) {
     filterHistos_->fill("numSignalPass",genProp.numDecays(),1,
			 "number of signal particles in event, pass filter",
			 10,0,10);
     
     filterHistos_->fill("numSignalEPass",n_elec,1,
			 "number of signal particles in event"
			 " decaying to ee, pass filter",
			 10,0,10);
     filterHistos_->fill("numSignalMuPass",n_muon,1,
			 "number of signal particles in event"
			 " decaying to mumu, pass filter",
			 10,0,10);
     filterHistos_->fill("numSignalTauPass",n_tau,1,
			 "number of signal particles in event"
			 " decaying to tautau, pass filter",
			 10,0,10);
     if (genProp.numDecays()==2) {
       filterHistos_->fill("numSignalEMuPass",n_elec,n_muon,1,
			   "number of signal particles in event decaying to"
			   " electrons (x) and muons (y)",10,0,10,10,0,10);
     }
     if (n_elec==1) {
       filterHistos_->fill("decayLength2D_1elec_pass",genProp.decayLength2D(i_elec),1,
			   "2d decay length, electron channel, 1 electron decay",
			   100,0,100);
     } else if (n_elec==2) {
       for (unsigned i=0; i<2; i++) 
	 filterHistos_->fill("decayLength2D_2elec_pass",genProp.decayLength2D(i),1,
			     "2d decay length, electron channel, 2 electron decays",
			     100,0,100);
     }
     if (n_muon==1) {
       filterHistos_->fill("decayLength2D_1muon_pass",genProp.decayLength2D(i_muon),1,
			   "2d decay length, muon channel, 1 muon decay",
			   100,0,100);
     } else if (n_muon==2) {
       for (unsigned i=0; i<2; i++)
	 filterHistos_->fill("decayLength2D_2muon_pass",genProp.decayLength2D(i),1,
			     "2d decay length, muon channel, 2 muon decays",
			     100,0,100);
     }
     if (n_tau==1) {
       filterHistos_->fill("decayLength2D_1tau_pass",genProp.decayLength2D(i_tau),1,
			   "2d decay length, tau channel, 1 tau decay",
			   100,0,100);
     } else if (n_tau==2) {
       for (unsigned i=0; i<2; i++)
	 filterHistos_->fill("decayLength2D_2tau_pass",genProp.decayLength2D(i),1,
			     "2d decay length, tau channel, 2 tau decays",
			     100,0,100);
     }
     for (unsigned i=0; i<genProp.numDecays(); i++) {
       filterHistos_->fill("ctau_pass",genProp.ctau(i),1,"ctau distribution",100,0,100);
     }
   }
   
   return pass_nIsolHighPt;
}


//define this as a plug-in
DEFINE_FWK_MODULE(IsoTrackPrefilter);
