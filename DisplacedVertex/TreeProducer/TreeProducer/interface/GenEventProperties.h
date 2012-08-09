#ifndef GenEventProperties_h
#define GenEventProperties_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class GenEventProperties
{
 private:

  // long-lived decay properties
  unsigned numDecays_;
  std::vector<int> decayModes_;
  std::vector<double> decayLength2D_;
  std::vector<double> decayLength3D_;
  std::vector<double> ctau_;

  // pile-up vertices
  float nvtx_m1_;
  float nvtx_0_;
  float nvtx_p1_;
  // float ave_nvtx_;

  // PDF description
  int pdf_id_[2];
  double pdf_x_[2];
  double pdf_xPDF_[2];
  double pdf_scale_;

 public:

  struct DecayLengthAndType
  {
    DecayLengthAndType(const int inputDecayPid,
                       const double & inputDecayLength2D,
                       const double & inputDecayLength3D,
                       const double & inputctau
                       ) :
      decayPid(inputDecayPid),
      decayLength2D(inputDecayLength2D),
      decayLength3D(inputDecayLength3D),
      ctau(inputctau)
    {}

    int decayPid;
    double decayLength2D;
    double decayLength3D;
    double ctau;
  };

  GenEventProperties();
  GenEventProperties(const edm::Event&,
                     int signalPDGId,
                     edm::InputTag& generatorTag,
                     edm::InputTag& pileupTag,
                     edm::InputTag& genEventInfoTag);

  static DecayLengthAndType getDecayLengthAndType(const reco::GenParticle & part);

  // accessors for decay properties
  unsigned numDecays() { return decayModes_.size(); }
  int decayMode(unsigned i) { if (i<2) return decayModes_[i]; else return -999; }
  double decayLength2D(unsigned i) { if (i<2) return decayLength2D_[i]; else return -999; }
  double decayLength3D(unsigned i) { if (i<2) return decayLength3D_[i]; else return -999; }
  double ctau(unsigned i) { if (i<2) return ctau_[i]; else return -999; }

  // accessors for pile-up
  float numPileupInTime() { return nvtx_0_; }
  float numPileup( int BX ) { if (BX!=0) return (BX==-1) ? nvtx_m1_ : nvtx_p1_; else return nvtx_0_; }
  // float aveNumPileup3BX() { return ave_nvtx_; }

  // accessors for PDF
  int pdfId(unsigned i) { if (i<2) return pdf_id_[i]; else return -999; }
  double pdfX(unsigned i) { if (i<2) return pdf_x_[i]; else return -999; }
  double pdfXPDF(unsigned i) { if (i<2) return pdf_xPDF_[i]; else return -999; }
  double pdfScale() { return pdf_scale_; }
};

#endif
