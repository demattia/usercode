#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HarderAnalysisTools/Histograms/interface/HistMap.h"
#include <string>

class SecondaryParticles : public edm::EDAnalyzer {

public:
  SecondaryParticles(const edm::ParameterSet&);
  ~SecondaryParticles();
  
private:

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void writeHistogram(std::string);

  HistMap* histos_;
};
