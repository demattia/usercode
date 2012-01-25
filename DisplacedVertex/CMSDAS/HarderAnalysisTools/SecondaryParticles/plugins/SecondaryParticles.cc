#include "HarderAnalysisTools/SecondaryParticles/plugins/SecondaryParticles.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include <cmath>
#include <fstream>


SecondaryParticles::SecondaryParticles(const edm::ParameterSet& iConfig)
{
}


SecondaryParticles::~SecondaryParticles()
{
}


void SecondaryParticles::beginJob(const edm::EventSetup&)
{
  histos_ = new HistMap();
}


void SecondaryParticles::writeHistogram(std::string histname)
{
  TH1* hist = histos_->getHist1D(histname);
  if (!hist) {
    std::cout << "ERROR: histogram " << histname << " not found" << std::endl;
    return;
  }

  ofstream fout((histname+".inc").c_str());
  if (fout.is_open()) {
    fout << "unsigned " << histname << "[" << hist->GetXaxis()->GetNbins() << "] = {" << std::endl;
    for (int i=1; i<hist->GetXaxis()->GetNbins(); i++) {
      fout << hist->GetBinContent(i) << "," << std::endl;
    }
    fout << hist->GetBinContent(hist->GetXaxis()->GetNbins()) << "};" <<std::endl;
    fout << "double " << histname << "_min=" << hist->GetXaxis()->GetXmin() << ";" << std::endl;
    fout << "double " << histname << "_max=" << hist->GetXaxis()->GetXmax() << ";" << std::endl;
    fout.close();
  } else {
    std::cout << "ERROR: could not open output file " << histname << ".inc" << std::endl;
  }
}


void SecondaryParticles::endJob()
{
  // dump histograms into files
  writeHistogram("dphi");
  writeHistogram("deta");
  writeHistogram("mom_e");
  writeHistogram("mom_pi");
  writeHistogram("mom_p");

  // distribution of particle types in secondary particles
  std::cout << std::endl << "CHARGED SECONDARY PARTICLE CONTRIBUTIONS:" << std::endl;
  double sum_pepi=histos_->getHist1D("mom_e")->GetEntries()
    +histos_->getHist1D("mom_pi")->GetEntries()
    +histos_->getHist1D("mom_p")->GetEntries();
  double sum_all=sum_pepi
    +histos_->getHist1D("mom_mu")->GetEntries()
    +histos_->getHist1D("mom_other")->GetEntries();
  if (sum_pepi==0) {
    std::cout << "no particles found!" << std::endl;
    return;
  }
  std::cout << "fraction of e+pi+p among all secondaries: " << sum_pepi/sum_all << std::endl;
  std::cout << "within e+pi+p only: e  fraction "
	    << histos_->getHist1D("mom_e")->GetEntries()/sum_pepi << std::endl;
  std::cout << "                    pi fraction "
	    << histos_->getHist1D("mom_pi")->GetEntries()/sum_pepi << std::endl;
  std::cout << "                    p  fraction "
	    << histos_->getHist1D("mom_p")->GetEntries()/sum_pepi << std::endl;
}


void SecondaryParticles::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::vector<edm::Handle<TrackingParticleCollection> >  TPCollections ;
  iEvent.getManyByType(TPCollections);

  if (TPCollections.size()>0) {
    const TrackingParticleCollection tpc   = *(TPCollections[0].product());
    for (size_t i=0; i<tpc.size(); i++) {
      const TrackingParticle& part = tpc[i];

      // skip particles coming from the MC generator. that leaves us mostly with secondaries
      // from material interactions
      if (part.genParticle().size()>0) continue;

      // also, we are only interested in charged particles
      if (part.charge()==0) continue;

      double z = part.vertex().z();
      double r = std::sqrt(part.vertex().x()*part.vertex().x()
			   +part.vertex().y()*part.vertex().y());
      histos_->fill("vertex_rz",z,r,1,"secondary particle vertex distribution r vs z",
		    100,-300,300,100,0,100);

      histos_->fill("phi",part.vertex().phi(),part.momentum().phi(),1,
		    "vertex phi coordinate vs momentum phi direction",
		    100,-3.142,3.142,100,-3.142,3.142);
      histos_->fill("dphi",part.vertex().phi()-part.momentum().phi(),1,
		    "vertex phi coordinate minus momentum phi direction",
		    100,-1,1);

      histos_->fill("eta",part.vertex().eta(),part.momentum().eta(),1,
		    "vertex eta coordinate vs momentum eta direction",
		    100,-5,5,100,-5,5);
      histos_->fill("deta",fabs(part.vertex().eta())-fabs(part.momentum().eta()),1,
		    "vertex eta coordinate minus momentum eta direction",
		    100,-1,5);

      double mom=std::sqrt(part.energy()*part.energy()-part.mass()*part.mass());
      histos_->fill("mass",part.mass(),1,"secondary particle mass",100,0,1.1);
      if (abs(part.pdgId())==211) {
	histos_->fill("mom_pi",mom,1,"secondary particle 3-momentum, pion",100,0,5);
      } else if (abs(part.pdgId())==11) {
	histos_->fill("mom_e",mom,1,"secondary particle 3-momentum, electron",100,0,5);
      } else if (abs(part.pdgId())==13) {
	histos_->fill("mom_mu",mom,1,"secondary particle 3-momentum, muon",100,0,5);
      } else if (abs(part.pdgId())==2212) {
	histos_->fill("mom_p",mom,1,"secondary particle 3-momentum, proton",100,0,5);
	histos_->fill("energy_p",part.energy(),1,"secondary particle energy, proton",100,0,5);
      } else {
	histos_->fill("mom_other",mom,1,"secondary particle 3-momentum, other",100,0,5);
      }
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SecondaryParticles);
