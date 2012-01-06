#include "DisplacedLeptons/Samples/interface/GenEventProperties.h"


#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

GenEventProperties::GenEventProperties() {

  // initialize data members with unphysical values
  decayModes_.clear();
  decayLength2D_.clear();
  decayLength3D_.clear();
  ctau_.clear();
  nvtx_        = -999;
  ave_nvtx_    = -999;
  pdf_id_[0]   = -999;
  pdf_id_[1]   = -999;
  pdf_x_[0]    = -999;
  pdf_x_[1]    = -999;
  pdf_xPDF_[0] = -999;
  pdf_xPDF_[1] = -999;
  pdf_scale_   = -999;
}


GenEventProperties::GenEventProperties(const edm::Event& iEvent,
				       int signalPDGId,
				       edm::InputTag& generatorTag,
				       edm::InputTag& pileupTag,
				       edm::InputTag& genEventInfoTag) {

  // initialize data members with unphysical values
  decayModes_.clear();
  decayLength2D_.clear();
  decayLength3D_.clear();
  ctau_.clear();
  nvtx_        = -999;
  ave_nvtx_    = -999;
  pdf_id_[0]   = -999;
  pdf_id_[1]   = -999;
  pdf_x_[0]    = -999;
  pdf_x_[1]    = -999;
  pdf_xPDF_[0] = -999;
  pdf_xPDF_[1] = -999;
  pdf_scale_   = -999;

  // look for long-lived particles and evaluate decay length
  edm::Handle<edm::View<reco::GenParticle> > mcParticles;
  iEvent.getByLabel(generatorTag,mcParticles);
  if (mcParticles.failedToGet()) return;

  for (unsigned imc=0; imc<mcParticles->size(); imc++) {
    if ((abs((*mcParticles)[imc].pdgId())==signalPDGId) &&
	((*mcParticles)[imc].status()==3)) {
      // long-lived particle and its production vertex
      const reco::GenParticle part=(*mcParticles)[imc];
      reco::Candidate::Point productionVertex = part.vertex();
      // decay products and decay vertex
      reco::Candidate::Point decayVertex=productionVertex;
      int decay_pid=0;
      for (unsigned i=0; i<part.numberOfDaughters(); i++) {
	const reco::Candidate* daughter = part.daughter(i);
	int pid = abs(daughter->pdgId());
	if (pid==11 || pid==13 || pid==15) {
	  decay_pid=pid;
	  // note that original MC signal particle and decay products
	  // are in documentation lines (Pythia status 3) and the actual
	  // decay length is only available in the actual decay tree
	  // (Pythia status 2). The corresponding status 2 particle for each
	  // status 3 particle is referenced by a daughter link.
	  if (daughter->numberOfDaughters()>0) {
	    decayVertex=daughter->daughter(0)->vertex();
	  }
	}
      }
      decayModes_.push_back(decay_pid);
      decayLength2D_.push_back((decayVertex-productionVertex).rho());
      double distance=(decayVertex-productionVertex).r();
      decayLength3D_.push_back(distance);
      ctau_.push_back(distance*part.p4().M()/part.p4().P());
    }
  }
  
  // pile-up information (if any)
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(pileupTag, PupInfo);
  if (!PupInfo.failedToGet()) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float sum_nvtx_3bx = 0;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      unsigned npv = PVI->getPU_NumInteractions();
      sum_nvtx_3bx += float(npv);
      int BX = PVI->getBunchCrossing();
      if (BX == 0) nvtx_=npv; 
    }
    ave_nvtx_ = sum_nvtx_3bx/3.;
  }
  
  // PDF information (if any)
  edm::Handle<GenEventInfoProduct>  genEventInfo;
  iEvent.getByLabel(genEventInfoTag, genEventInfo);
  if (!genEventInfo.failedToGet()) {
    if (genEventInfo->hasPDF()) {
      const gen::PdfInfo* pdf = genEventInfo->pdf();
      pdf_id_[0]=pdf->id.first;
      pdf_id_[1]=pdf->id.second;
      pdf_x_[0]=pdf->x.first;
      pdf_x_[1]=pdf->x.second;
      pdf_xPDF_[0]=pdf->xPDF.first;
      pdf_xPDF_[1]=pdf->xPDF.second;
      pdf_scale_=pdf->scalePDF;
    }
  }
}
