/**
 * This macro shows how to access FW objects in an uncompiled,
 * nameless FWLite macro. It uses directly the fwlite::ChainEvent method.
 * As of now, the MC information is not accessed.
 */
{
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");

  gROOT->SetBatch(true);

  TH2F * ptVSetaGen = new TH2F("ptVSetaGen","genMuon pt vs eta", 100, -3.2, 3.2, 100, 0, 120);
  TProfile * ptVSetaGenProfile = new TProfile("ptVSetaGenProfile", "genMuon pt vs eta", 100, -3.2, 3.2, 0, 120);

  TH2F * ptVSeta = new TH2F("ptVSeta","muon pt vs eta", 100, -3.2, 3.2, 100, 0, 120);
  TProfile * ptVSetaProfile = new TProfile("ptVSetaProfile", "muon pt vs eta", 100, -3.2, 3.2, 0, 120);

#include "DataFormats/FWLite/interface/Handle.h"
  vector<string> fileNames;
  fileNames.push_back("/data2/demattia/Data/MuonGun/MuonGun_0.root");

  fwlite::ChainEvent ev(fileNames);

  for( ev.toBegin(); !ev.atEnd(); ++ev ) {
    fwlite::Handle<std::vector<reco::Muon> > muons;
    muons.getByLabel(ev,"muons");
    // now can access data
    // Either of the two is ok.
    // std::cout <<" size "<<muons.ptr()->size()<<std::endl;
    // std::cout <<" size "<<muons->size()<<std::endl;
    for(unsigned int muonId=0; muonId<muons->size(); ++muonId) {
      reco::Muon * muon = &((*muons)[muonId]);
      // cout << "(*muons)["<<muonId<<"] = " << muon << endl;
      ptVSeta->Fill(muon->eta(), muon->pt());
      ptVSetaProfile->Fill(muon->eta(), muon->pt());
    }
  }

  TCanvas recoMuonCanvas("ptVSeta","reco muon pt VS eta", 1000, 800);
  recoMuonCanvas.cd();
  ptVSeta->Draw();
  recoMuonCanvas->Print("ptVSetaReco.pdf");
  TCanvas recoMuonCanvasProfile("ptVSeta","reco muon pt VS eta", 1000, 800);
  recoMuonCanvasProfile.cd();
  ptVSetaProfile->Draw();
  recoMuonCanvasProfile.Print("ptVSetaRecoProfile.pdf");
}
