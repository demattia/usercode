#ifndef TREECANDIDATE_H
#define TREECANDIDATE_H

#include <TObject.h>
#include <TLorentzVector.h>

struct TreeCandidate
{
  TreeCandidate() :
    leptonPtL(-999.), leptonPtH(-999.), leptonEtaL(-999.), leptonEtaH(-999.), leptonChargeL(-999), leptonChargeH(-999), leptonIndexL(-999), leptonIndexH(-999),
    leptonD0L(-999.), leptonD0H(-999.), leptonAbsD0SignificanceL(-999.), leptonAbsD0SignificanceH(-999.), leptonIsoL(-999.), leptonIsoH(-999.),
    triggerMatchL(-999), triggerMatchH(-999), cosine(-999.), deltaR(-999.), cosThetaStar(-999.), mass(-999.), corrDileptonMass(-999.),
    triggerCorrDileptonMass(-999.), caloCorrMass(-999.), scaleCorrMass(-999.), eta(-999.), phi(-999.), etaCorr(-999.), phiCorr(-999.),
    etaTriggerCorr(-999.), phiTriggerCorr(-999.), etaCaloCorr(-999.), phiCaloCorr(-999.), decayLength(-999.), decayLengthSignificance(-999.),
    dPhi(-999.), dPhiCorr(-999.), dPhiTriggerCorr(-999.), dPhiCaloCorr(-999.), hitsBeforeVertexL(-999.), hitsBeforeVertexH(-999.),
    missedLayersAfterVertexL(-999.), missedLayersAfterVertexH(-999.), isStandAloneL(-999), isStandAloneH(-999), isGlobalMuonL(-999), isGlobalMuonH(-999), hasCaloMatchL(-999), hasCaloMatchH(-999),
    originPdgIdL(-999), originPdgIdH(-999), pdgIdL(-999), pdgIdH(-999), genctau(-999.), genDecayLength2D(-999.), genDecayLength3D(-999.),
    validVertex(-999), vertexChi2(-999.), vx(-999.), vy(-999.), vz(-999.)
  {}

  // Single lepton variables
  Float_t leptonPtL;
  Float_t leptonPtH;
  Float_t leptonEtaL;
  Float_t leptonEtaH;
  Int_t leptonChargeL;
  Int_t leptonChargeH;
  Int_t leptonIndexL;
  Int_t leptonIndexH;
  Float_t leptonD0L;
  Float_t leptonD0H;
  Float_t leptonAbsD0SignificanceL;
  Float_t leptonAbsD0SignificanceH;
  Float_t leptonIsoL;
  Float_t leptonIsoH;
  Int_t triggerMatchL;
  Int_t triggerMatchH;
  // Di-lepton variables
  Float_t cosine; // to reject cosmics. Leave unbiased by taking the values before refit to vertex
  Float_t deltaR; // for trigger inefficiencies. Leave unbiased as above.
  Float_t cosThetaStar; // angle between positive lepton momentum in dilepton rest frame and dilepton momentum
  // In the following:
  // - the default values are from the leptons
  // - the "corr" values are from the leptons after refit to common vertex
  // - the "triggerCorr" values are from the trigger matches
  // - the "caloCorr" values are from the calo matches
  Float_t mass;
  Float_t corrDileptonMass;
  Float_t triggerCorrDileptonMass;
  Float_t caloCorrMass;
  Float_t scaleCorrMass; // corrected so that deltaR between vertex flight direction and di-lepton momentum is 0.
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t ptCorr;
  Float_t etaCorr;
  Float_t phiCorr;
  Float_t ptTriggerCorr;
  Float_t etaTriggerCorr;
  Float_t phiTriggerCorr;
  Float_t ptCaloCorr;
  Float_t etaCaloCorr;
  Float_t phiCaloCorr;

  Float_t decayLength;
  Float_t decayLengthSignificance;
  Float_t dPhi; // Angle in transverse plane between vertex (secondary-primay) flight direction and di-lepton momentum
  Float_t dPhiCorr; // Same as above but using leptons refitted to vertex to compute di-lepton momentum
  Float_t dPhiTriggerCorr; // Same as above using trigger matches to compute di-lepton momentum
  Float_t dPhiCaloCorr;

  Float_t hitsBeforeVertexL;
  Float_t hitsBeforeVertexH;
  Float_t missedLayersAfterVertexL;
  Float_t missedLayersAfterVertexH;
  Int_t isStandAloneL;
  Int_t isStandAloneH;
  Int_t isGlobalMuonL;
  Int_t isGlobalMuonH;
  Int_t isTrackerMuonL;
  Int_t isTrackerMuonH;
  Int_t hasCaloMatchL;
  Int_t hasCaloMatchH;
  Int_t originPdgIdL;
  Int_t originPdgIdH;
  Int_t pdgIdL;
  Int_t pdgIdH;
  Float_t genctau;
  Float_t genDecayLength2D;
  Float_t genDecayLength3D;

  Int_t validVertex;
  Float_t vertexChi2;
  Float_t vx;
  Float_t vy;
  Float_t vz;

  ClassDef(TreeCandidate, 1)
};
ClassImp(TreeCandidate)

#endif
