// -*- C++ -*-
//
// Package:    EfficiencyRatioProducer
// Class:      EfficiencyRatioProducer
//
/**\class EfficiencyRatioProducer EfficiencyRatioProducer.cc Analysis/EfficiencyRatioProducer/src/EfficiencyRatioProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,42 R-23,
//         Created:  Mon Jul 4 18:38:0 CEST 2011
// $Id: EfficiencyRatioProducer.cc,v 1.2 2011/07/15 15:53:04 demattia Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Analysis/TrackingEfficiencyFromCosmics/interface/EfficiencyTree.h"

#include <boost/foreach.hpp>

#include <TCanvas.h>
#include <TGraphErrors.h>

//
// class declaration
//

class EfficiencyRatioProducer : public edm::EDAnalyzer {
public:
  explicit EfficiencyRatioProducer(const edm::ParameterSet&);
  ~EfficiencyRatioProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void fillHistogram(const TString & name, const TString & title,
                     const boost::shared_ptr<Efficiency> & effRatioNumerator,
                     const boost::shared_ptr<Efficiency> & effRatioDenominator);

  // ----------member data ---------------------------
  std::auto_ptr<Efficiency> efficiencyRatioNumerator_;
  std::auto_ptr<Efficiency> efficiencyRatioDenominator_;
  std::string inputFileNameNumerator_, inputFileNameDenominator_;
  unsigned int rebin_;
};

void EfficiencyRatioProducer::fillHistogram(const TString & name, const TString & title,
                                            const boost::shared_ptr<Efficiency> & effRatioNumerator,
                                            const boost::shared_ptr<Efficiency> & effRatioDenominator )
{
  TCanvas *c1 = new TCanvas("c1"+name,"A Simple Graph Example",200,10,700,500);
  Int_t n = effRatioNumerator->bins(0);
  std::cout << "number of bins = " << n << std::endl;
  Double_t x[100], y[100];
  Double_t ex[100];
  Double_t ey[100];
  for (Int_t i=0;i<n;++i) {
    double num = effRatioNumerator->getEff(i);
    double numE = effRatioNumerator->getEffError(i);
    if( numE == 0 ) {
      numE = 0.0001;
      std::cout << "error on the numerator = 0! Setting it to " << numE << std::endl;
    }
    double den = effRatioDenominator->getEff(i);
    double denE = effRatioDenominator->getEffError(i);
    if( denE == 0 ) {
      denE = 0.0001;
      std::cout << "error on the denominator = 0! Setting it to " << denE << std::endl;
    }
    std::cout << "num = " << num << std::endl;
    std::cout << "den = " << den << std::endl;
    x[i] = effRatioNumerator->min(0) + (i+0.5)*effRatioNumerator->binsSize(0);
    if( den == 0 ) {
      std::cout << "denominator = 0! Setting the ratio to 0" << std::endl;
      y[i] = 0;
    }
    else {
      y[i] = num/den;
    }
    ex[i] = effRatioNumerator->binsSize(0)/2.;
    ey[i] = y[i]*sqrt(pow(numE/num,2)+pow(denE/den,2));
    std::cout << "x["<<i<<"] = " << x[i] << std::endl;
    std::cout << "ex["<<i<<"] = " << ex[i] << std::endl;
    std::cout << "y["<<i<<"] = " << y[i] << std::endl;
    std::cout << "ey["<<i<<"] = " << ey[i] << std::endl;
  }
  std::cout << "Building the TGraphErrors" << std::endl;
  TGraphErrors * hEffRatio = new TGraphErrors(n, x, y, ex, ey);
  hEffRatio->SetName(name);
  hEffRatio->SetTitle(name);
  hEffRatio->Draw("AP");
  c1->Draw();
  c1->SaveAs(name+".root");
}

EfficiencyRatioProducer::EfficiencyRatioProducer(const edm::ParameterSet& iConfig) :
  inputFileNameNumerator_(iConfig.getParameter<std::string>("InputFileNameNumerator")),
  inputFileNameDenominator_(iConfig.getParameter<std::string>("InputFileNameDenominator")),
  rebin_(iConfig.getParameter<unsigned int>("Rebin"))
{
  efficiencyRatioNumerator_.reset(new Efficiency);
  EfficiencyTree treeRatioNumerator;
  treeRatioNumerator.readTree(inputFileNameNumerator_, &*efficiencyRatioNumerator_);
  boost::shared_array<unsigned int> vKeep(new unsigned int[3]);
  vKeep[0] = rebin_;
  vKeep[1] = 0;
  vKeep[2] = 0;
  boost::shared_ptr<Efficiency> effVsDxyRatioNumerator(efficiencyRatioNumerator_->projectAndRebin(vKeep));

  efficiencyRatioDenominator_.reset(new Efficiency);
  EfficiencyTree treeRatioDenominator;
  treeRatioDenominator.readTree(inputFileNameDenominator_, &*efficiencyRatioDenominator_);
  vKeep[0] = rebin_;
  vKeep[1] = 0;
  vKeep[2] = 0;
  boost::shared_ptr<Efficiency> effVsDxyRatioDenominator(efficiencyRatioDenominator_->projectAndRebin(vKeep));

  TFile * outputFile = new TFile("EfficiencyRatioProducer_1.root", "RECREATE");
  outputFile->cd();

  fillHistogram("EffVsDxyRatio", "Ratio of efficiencies vs absolute transverse impact parameter", effVsDxyRatioNumerator, effVsDxyRatioDenominator);


  //  vKeep[0] = 0;
  //  vKeep[1] = 1;
  //  vKeep[2] = 0;
  //  boost::shared_ptr<Efficiency> effVsDz(efficiency_->projectAndRebin(vKeep));
  //  fillHistogram("EffVsDz", "Efficiency vs absolute longitudinal impact parameter", effVsDz);

  //  vKeep[0] = 0;
  //  vKeep[1] = 0;
  //  vKeep[2] = 1;
  //  boost::shared_ptr<Efficiency> effVsPt(efficiency_->projectAndRebin(vKeep));
  //  fillHistogram("EffVsPt", "Efficiency vs Pt", effVsPt);

}

EfficiencyRatioProducer::~EfficiencyRatioProducer() {}

void EfficiencyRatioProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just before starting event loop  ------------
void EfficiencyRatioProducer::beginJob()
{
  // edm::Service<TFileService> fileService;
}

// ------------ method called once each job just after ending the event loop  ------------
void EfficiencyRatioProducer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void EfficiencyRatioProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void EfficiencyRatioProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void EfficiencyRatioProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void EfficiencyRatioProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EfficiencyRatioProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyRatioProducer);


