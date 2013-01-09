// -*- C++ -*-
//
// Package:    GenSpectrumFilter
// Class:      GenSpectrumFilter
// 
/**\class GenSpectrumFilter GenSpectrumFilter.cc PhysicsTools/PatAlgos/src/GenSpectrumFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nov 16 16:12 (lxplus231.cern.ch)
//         Created:  Sun Nov 16 16:14:09 CET 2008
// $Id: GenSpectrumFilter.cc,v 1.1 2010/03/07 21:05:33 gpetrucc Exp $
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

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// class decleration
class GenSpectrumFilter : public edm::EDFilter {
    public:
        explicit GenSpectrumFilter(const edm::ParameterSet&);
        ~GenSpectrumFilter();

    private:
        virtual bool filter(edm::Event&, const edm::EventSetup&);

        /// The GEN objects
        edm::InputTag src_;

        /// The pdgId of the particle
        std::vector<int32_t> pdgIds_;

        typedef StringObjectFunction<reco::GenParticle> Function;
        /// The Probability
        Function probability_;

        /// The Random number generator
        std::auto_ptr<CLHEP::RandFlat> flatDistribution_;
};

GenSpectrumFilter::GenSpectrumFilter(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    pdgIds_(iConfig.getParameter<std::vector<int32_t> >("pdgIds")),
    probability_(iConfig.getParameter<std::string>("probability")),
    flatDistribution_(0)
{
    edm::Service<edm::RandomNumberGenerator> rng;
    if ( ! rng.isAvailable()) {
        throw cms::Exception("Configuration")
            << "XXXXXXX requires the RandomNumberGeneratorService\n"
            "which is not present in the configuration file.  You must add the service\n"
            "in the configuration file or remove the modules that require it.";
    }

    CLHEP::HepRandomEngine& engine = rng->getEngine();

    // engine MUST be a reference here, if a pointer is used the
    // distribution will destroy the engine in its destructor, a major
    // problem because the service owns the engine and will destroy it 
    flatDistribution_.reset(new CLHEP::RandFlat(engine, 0, 1));
}

GenSpectrumFilter::~GenSpectrumFilter() 
{
}

bool
GenSpectrumFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> src; 
    iEvent.getByLabel(src_, src);

    foreach (int32_t pdgId, pdgIds_) {
        foreach (const reco::GenParticle &gen, *src) {
            if (gen.pdgId() == pdgId) {
                return flatDistribution_->fire() <= probability_(gen);
            }
        } 
    }

    return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSpectrumFilter);
