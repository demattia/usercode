// -*- C++ -*-
//
// Package:    EfficiencyAnalyzer
// Class:      EfficiencyAnalyzer
//
/**\class EfficiencyAnalyzer EfficiencyAnalyzer.cc Analysis/EfficiencyAnalyzer/src/EfficiencyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,42 R-23,
//         Created:  Mon Jul 4 18:38:0 CEST 2011
// $Id: EfficiencyAnalyzer.cc,v 1.5 2011/10/18 16:18:24 demattia Exp $
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
#include <TGraphAsymmErrors.h>

//
// class declaration
//

class EfficiencyAnalyzer : public edm::EDAnalyzer {
public:
  explicit EfficiencyAnalyzer(const edm::ParameterSet&);
  ~EfficiencyAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void fillHistogram(const TString & name, const TString & title, const boost::shared_ptr<Efficiency> & eff );

  // ----------member data ---------------------------
  std::auto_ptr<Efficiency> efficiency_;
  std::string inputFileName_;
  unsigned int rebin_;
  // edm::Service<TFileService> fileService_;
};

void EfficiencyAnalyzer::fillHistogram(const TString & name, const TString & title, const boost::shared_ptr<Efficiency> & eff )
{
  // TFile *outputFile = new TFile(name+".root", "RECREATE");
  // outputFile->cd();
  // TCanvas *c1 = new TCanvas("c1"+name,"A Simple Graph Example",200,10,700,500);
  unsigned int n = eff->bins(0);
  std::cout << "Filling " << n << " bins" << std::endl;
  Double_t x[100], y[100];
  Double_t exl[100], exh[100];
  Double_t eyl[100], eyh[100];
  for (unsigned i=0;i<n;++i) {
    // std::cout << eff->min(0) << " "<<i<<"*"<<eff->binsSize(0)<<" = "<< i*(eff->binsSize(0)) << std::endl;
    x[i] = eff->min(0) + (i+0.5)*eff->binsSize(0);
    y[i] = eff->getEff(i);
    exl[i] = eff->binsSize(0)/2.;
    exh[i] = eff->binsSize(0)/2.;
    eyl[i] = std::min(y[i], eff->getEffError(i));
    eyh[i] = std::min(1-y[i], eff->getEffError(i));
  }
  // TGraph * hEff = fileService_->make<TGraph>(name, title, eff->bins(0), eff->min(0), eff->min(0));
  TGraphAsymmErrors * hEff = new TGraphAsymmErrors(n, x, y, exl, exh, eyl, eyh);
  hEff->SetName(name);
  hEff->Draw("AP");
  // c1->Draw();
  // c1->Write();
  // outputFile->Write();
  // outputFile->Close();
  // c1->SaveAs(name+".root");
  hEff->SaveAs(name+".root");
  // c1->Save();

//  for( unsigned int i=0; i<eff->getLinearSize(); ++i ) {
//    std::cout << "eff["<<i<<"] = " << eff->getEff(i) << std::endl;
//    hEff->SetBinContent(i+1, eff->getEff(i));
//  }
}

EfficiencyAnalyzer::EfficiencyAnalyzer(const edm::ParameterSet& iConfig) :
  inputFileName_(iConfig.getParameter<std::string>("InputFileName")),
  rebin_(iConfig.getParameter<unsigned int>("Rebin"))
{
  efficiency_.reset(new Efficiency);
  EfficiencyTree tree;
  tree.readTree(inputFileName_, &*efficiency_);

  TFile * outputFile = new TFile("EfficiencyAnalyzer_1.root", "RECREATE");
  outputFile->cd();

  boost::shared_array<unsigned int> vKeep(new unsigned int[4]);
  vKeep[0] = rebin_;
  vKeep[1] = 0;
  vKeep[2] = 0;
  vKeep[3] = 0;

  boost::shared_ptr<Efficiency> effClonedForDxy(efficiency_->clone());
  // Apply a cut on dz
  effClonedForDxy->cut(1, 0, 10);
  // boost::shared_ptr<Efficiency> effVsDxy(efficiency_->projectAndRebin(vKeep));
  boost::shared_ptr<Efficiency> effVsDxy(effClonedForDxy->projectAndRebin(vKeep));
  fillHistogram("EffVsDxy", "Efficiency vs absolute transverse impact parameter", effVsDxy );

  unsigned int S = effVsDxy->getLinearSize();
  for( unsigned int i=0; i<S; ++i ) {
    std::cout << "reco eff vs Dxy ["<<i<<"] = " << effVsDxy->getValues(i).second<<"/"<<effVsDxy->getValues(i).first << " = "
	      << effVsDxy->getEff(i) << " +/- " << effVsDxy->getEffError(i) << std::endl;
  }

  vKeep[0] = 0;
  vKeep[1] = rebin_;
  vKeep[2] = 0;
  vKeep[3] = 0;

  boost::shared_ptr<Efficiency> effClonedForDz(efficiency_->clone());
  // Apply a cut on dxy  
  effClonedForDz->cut(0, 0, 4);
  // boost::shared_ptr<Efficiency> effVsDz(efficiency_->projectAndRebin(vKeep));
  boost::shared_ptr<Efficiency> effVsDz(effClonedForDz->projectAndRebin(vKeep));
  fillHistogram("EffVsDz", "Efficiency vs absolute longitudinal impact parameter", effVsDz);

  vKeep[0] = 0;
  vKeep[1] = 0;
  vKeep[2] = rebin_;
  vKeep[3] = 0;
  boost::shared_ptr<Efficiency> effVsPt(efficiency_->projectAndRebin(vKeep));
  fillHistogram("EffVsPt", "Efficiency vs Pt", effVsPt);


  vKeep[0] = 0;
  vKeep[1] = 0;
  vKeep[2] = 0;
  vKeep[3] = rebin_;
  boost::shared_ptr<Efficiency> effClonedForEta(efficiency_->clone());
  effClonedForEta->cut(0, 0, 4);
  effClonedForEta->cut(1, 0, 10);
  boost::shared_ptr<Efficiency> effVsEta(effClonedForEta->projectAndRebin(vKeep));
  fillHistogram("EffVsEta", "Efficiency vs #eta", effVsEta);

  unsigned int Seta = effVsEta->getLinearSize();
  for( unsigned int i=0; i<Seta; ++i ) {
    std::cout << "reco eff vs Eta ["<<i<<"] = " << effVsEta->getValues(i).second<<"/"<<effVsEta->getValues(i).first << " = "
	      << effVsEta->getEff(i) << " +/- " << effVsEta->getEffError(i) << std::endl;
  }


  // unsigned int S = efficiency_->getLinearSize();
  // for( unsigned int i=0; i<S; ++i ) {
  //   std::cout << "reco eff["<<i<<"] = " << efficiency_->getEff(i) << " +/- " << efficiency_->getEffError(i) << std::endl;
  // }
}

EfficiencyAnalyzer::~EfficiencyAnalyzer() {}

void EfficiencyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just before starting event loop  ------------
void EfficiencyAnalyzer::beginJob()
{
  // edm::Service<TFileService> fileService;
}

// ------------ method called once each job just after ending the event loop  ------------
void EfficiencyAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void EfficiencyAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void EfficiencyAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void EfficiencyAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void EfficiencyAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EfficiencyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyAnalyzer);

