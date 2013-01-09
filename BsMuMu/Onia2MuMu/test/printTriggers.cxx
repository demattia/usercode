/** Debug macro to print Onia2MuMuPAT events and all the triggers fired by the event and all its muons */
#include "DataFormats/FWLite/interface/Handle.h"
#include <TFile.h>
#include <iostream>

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h" 
#endif

void printTriggers(const char *process = "HLT") {
    fwlite::Event ev(gFile);
    fwlite::Handle<edm::TriggerResults> hTriggerResults;
    fwlite::Handle<std::vector<pat::CompositeCandidate> > hOnia;
    fwlite::Handle<std::vector<l1extra::L1MuonParticle> > hL1;
#if !defined(__CINT__) && !defined(__MAKECINT__)
    fwlite::Handle<pat::TriggerObjectStandAloneCollection> hHLTObjs;
#endif

    int iEvent = 0;
    for (ev.toBegin(); ! ev.atEnd(); ++ev) {
        ++iEvent;
        hTriggerResults.getByLabel(ev,"TriggerResults","",process);
        edm::TriggerNames const&  triggerNames = ev.triggerNames(*hTriggerResults);
        std::cout << "\n\n\n === EVENT " << iEvent << ": HLT RESULTS ====" << std::endl;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            if (strstr(triggerNames.triggerName(i).c_str(), "Mu")) {
                std::cout << "\t" << triggerNames.triggerName(i) << ": " << (hTriggerResults->accept(i) ? "PASS" : "FAIL" ) << std::endl;
            }
        }

        hL1.getByLabel(ev,"l1extraParticles");
        if (!hL1.failedToGet()) {
            std::cout << " === EVENT " << iEvent << ": L1 Mu Particles ====" << std::endl;
            for (unsigned i = 0; i < hL1->size(); ++i) {
                const l1extra::L1MuonParticle &l1 = (*hL1)[i];
                std::cout << "\tL1 "<< i <<": pt " << l1.pt() << ", eta " << l1.eta() << ", phi " << l1.phi() << ", quality " << l1.gmtMuonCand().quality();
                if (l1.gmtMuonCand().quality() == 7) {
                    std::cout << ", Matched (DT or CSC) + RPC" << std::endl;
                } else if (l1.gmtMuonCand().isRPC()) {
                    std::cout << ", RPC" << std::endl;
                } else {
                    std::cout << (l1.gmtMuonCand().isFwd() ? ", CSC" : ", DT")  << std::endl;
                }
            }
        }

#if !defined(__CINT__) && !defined(__MAKECINT__)
        hHLTObjs.getByLabel(ev,"patTrigger");
        if (!hHLTObjs.failedToGet()) {
            std::cout << " === EVENT " << iEvent << ": ALL HLT Objects ====" << std::endl;
            const pat::TriggerObjectStandAloneCollection & objs = *hHLTObjs;
            for (unsigned i = 0; i < objs.size(); ++i) {
                const pat::TriggerObjectStandAlone &o = objs[i];
                std::cout << "\tTrigger object " << i << ", pt " << o.pt() << ", eta " << o.eta() << ", phi " << o.phi() << std::endl;
                std::cout << "\t   Collection: " << o.collection() << std::endl;
                std::cout << "\t   Type IDs:   ";
                for (unsigned h = 0; h < o.filterIds().size(); ++h) std::cout << " " << o.filterIds()[h] ;
                std::cout << std::endl;
                std::cout << "\t   Filters:    ";
                for (unsigned h = 0; h < o.filterLabels().size(); ++h) std::cout << " " << o.filterLabels()[h];
                std::cout << std::endl;
                std::cout << "\t   Paths  :    ";
                for (unsigned h = 0; h < o.pathNames().size(); ++h) std::cout << " " << o.pathNames()[h];
                std::cout << std::endl;
            }
        }
#endif

        std::cout << " === EVENT " << iEvent << ": ONIA ====" << std::endl;
        hOnia.getByLabel(ev,"onia2MuMuPatTrkTrk");
        const std::vector<pat::CompositeCandidate> & onias = *hOnia;
        for (unsigned i = 0; i < onias.size(); ++i) {
            std::cout << "  onia " << i << ", mass " << onias[i].mass() << ", charge " << onias[i].charge() << std::endl;
            for (unsigned j = 0; j < 2; ++j) {
                const pat::Muon *mu = dynamic_cast<const pat::Muon *>(onias[i].daughter(j));
                std::cout << "     muon " << j << " pt " << mu->pt() << ", eta " << mu->eta() << ", global? " << mu->isGlobalMuon() << std::endl;
                if (!mu->triggerObjectMatchesByFilter("propagatedToM2").empty()) std::cout << "\t\t- Reco can be propagated at muon station 2\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltL1MuOpenL1Filtered0").empty()) std::cout << "\t\t- L1MuOpen\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltL2Mu0L2Filtered0").empty()) std::cout << "\t\t- L2Mu0\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltSingleMu3L2Filtered3").empty()) std::cout << "\t\t- L2Mu3\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltSingleMu3L3Filtered3").empty()) std::cout << "\t\t- Mu3\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltSingleMu5L3Filtered5").empty()) std::cout << "\t\t- Mu5\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()) std::cout << "\t\t- L1DoubleMuOpen\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDiMuonL2PreFiltered0").empty()) std::cout << "\t\t- L2DoubleMu0\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty()) std::cout << "\t\t- DoubleMu0\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) std::cout << "\t\t- DoubleMu3\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu0L1MuOpenL3Filtered0").empty()) std::cout << "\t\t- Mu0+L1MuOpen (L3)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu3L1MuOpenL3Filtered3").empty()) std::cout << "\t\t- Mu3+L1MuOpen (L3)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu5L1MuOpenL3Filtered5").empty()) std::cout << "\t\t- Mu5+L1MuOpen (L3)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu0L1MuOpenL1Filtered0").empty()) std::cout << "\t\t- Mu0+L1MuOpen (L1)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu3L1MuOpenL1Filtered0").empty()) std::cout << "\t\t- Mu3+L1MuOpen (L1)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu5L1MuOpenL1Filtered0").empty()) std::cout << "\t\t- Mu5+L1MuOpen (L1)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu0L2Mu0L3Filtered0").empty()) std::cout << "\t\t- Mu0+L2Mu0 (L3)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu3L2Mu0L3Filtered3").empty()) std::cout << "\t\t- Mu3+L2Mu0 (L3)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltMu5L2Mu0L3Filtered5").empty()) std::cout << "\t\t- Mu5+L2Mu0 (L3)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDiMuonL2PreFiltered0").empty()) std::cout << "\t\t- Mu0+L2Mu0 (L2; same as L2DoubleMu0)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDiMuonL2PreFiltered0").empty()) std::cout << "\t\t- Mu3+L2Mu0 (L2; same as L2DoubleMu0)\n"; 
                if (!mu->triggerObjectMatchesByFilter("hltDiMuonL2PreFiltered0").empty()) std::cout << "\t\t- Mu5+L2Mu0 (L2; same as L2DoubleMu0)\n"; 
                /// MuX + Track works only compiled, not in CINT
#if !defined(__CINT__) && !defined(__MAKECINT__)
                pat::TriggerObjectStandAloneCollection mu0tkMatch = mu->triggerObjectMatchesByFilter("hltMu0TrackJpsiTrackMassFiltered");
                bool mu0tkMu = false, mu0tkTk = false;
                for (unsigned k = 0; k < mu0tkMatch.size(); ++k) {
                    if (mu0tkMatch[k].collection() == "hltL3MuonCandidates::HLT") mu0tkMu = true;
                    if (mu0tkMatch[k].collection() == "hltMuTrackJpsiCtfTrackCands::HLT") mu0tkTk = true;
                }
                if (mu0tkMu) std::cout << "\t\t- Mu0+Track (Mu)\n";
                if (mu0tkTk) std::cout << "\t\t- Mu0+Track (Tk)\n";
                pat::TriggerObjectStandAloneCollection mu3tkMatch = mu->triggerObjectMatchesByFilter("hltMu3TrackJpsiTrackMassFiltered");
                bool mu3tkMu = false, mu3tkTk = false;
                for (unsigned k = 0; k < mu3tkMatch.size(); ++k) {
                    if (mu3tkMatch[k].collection() == "hltL3MuonCandidates::HLT") mu3tkMu = true;
                    if (mu3tkMatch[k].collection() == "hltMuTrackJpsiCtfTrackCands::HLT") mu3tkTk = true;
                }
                if (mu3tkMu) std::cout << "\t\t- Mu3+Track (Mu)\n";
                if (mu3tkTk) std::cout << "\t\t- Mu3+Track (Tk)\n";
                pat::TriggerObjectStandAloneCollection mu5tkMatch = mu->triggerObjectMatchesByFilter("hltMu5TrackJpsiTrackMassFiltered");
                bool mu5tkMu = false, mu5tkTk = false;
                for (unsigned k = 0; k < mu5tkMatch.size(); ++k) {
                    if (mu5tkMatch[k].collection() == "hltL3MuonCandidates::HLT") mu5tkMu = true;
                    if (mu5tkMatch[k].collection() == "hltMuTrackJpsiCtfTrackCands::HLT") mu5tkTk = true;
                }
                if (mu5tkMu) std::cout << "\t\t- Mu5+Track (Mu)\n";
                if (mu5tkTk) std::cout << "\t\t- Mu5+Track (Tk)\n";
#endif
            } // muon  
        } // onia
        //if (iEvent >= 10) break;
    } // event
} // macro
