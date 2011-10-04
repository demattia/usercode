//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep 25 21:04:12 2011 by ROOT version 5.28/00b
// from TTree HltTree/
// found on file: openhlt_merge.root
//////////////////////////////////////////////////////////

#ifndef checkOpenHLT_h
#define checkOpenHLT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <sstream>

class checkOpenHLT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           NoRecoPFTausSignal;
   Int_t           signalTrToPFTauMatch[1];   //[NoRecoPFTausSignal]
   Float_t         recoPFTauSignalTrDz[1];   //[NoRecoPFTausSignal]
   Float_t         recoPFTauSignalTrPt[1];   //[NoRecoPFTausSignal]
   Int_t           NoRecoPFTausIso;
   Int_t           isoTrToPFTauMatch[1];   //[NoRecoPFTausIso]
   Float_t         recoPFTauIsoTrDz[1];   //[NoRecoPFTausIso]
   Float_t         recoPFTauIsoTrPt[1];   //[NoRecoPFTausIso]
   Int_t           NoHLTPFTausSignal;
   Int_t           hltpftauSignalTrToPFTauMatch[54];   //[NoHLTPFTausSignal]
   Float_t         HLTPFTauSignalTrDz[54];   //[NoHLTPFTausSignal]
   Float_t         HLTPFTauSignalTrPt[54];   //[NoHLTPFTausSignal]
   Int_t           NoHLTPFTausIso;
   Int_t           hltpftauIsoTrToPFTauMatch[35];   //[NoHLTPFTausIso]
   Float_t         HLTPFTauIsoTrDz[35];   //[NoHLTPFTausIso]
   Float_t         HLTPFTauIsoTrPt[35];   //[NoHLTPFTausIso]
   Int_t           NrecoJetGen;
   Int_t           NrecoTowCal;
   Int_t           NrecoJetCal;
   Float_t         recoJetCalPt[1];   //[NrecoJetCal]
   Float_t         recoJetCalPhi[1];   //[NrecoJetCal]
   Float_t         recoJetCalEta[1];   //[NrecoJetCal]
   Float_t         recoJetCalE[1];   //[NrecoJetCal]
   Float_t         recoJetCalEMF[1];   //[NrecoJetCal]
   Float_t         recoJetCalN90[1];   //[NrecoJetCal]
   Float_t         recoJetCalN90hits[1];   //[NrecoJetCal]
   Int_t           NrecoJetCorCal;
   Float_t         recoJetCorCalPt[1];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalPhi[1];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalEta[1];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalE[1];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalEMF[1];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalN90[1];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalN90hits[1];   //[NrecoJetCorCal]
   Int_t           NohJetCal;
   Float_t         ohJetCalPt[82];   //[NohJetCal]
   Float_t         ohJetCalPhi[82];   //[NohJetCal]
   Float_t         ohJetCalEta[82];   //[NohJetCal]
   Float_t         ohJetCalE[82];   //[NohJetCal]
   Float_t         ohJetCalEMF[82];   //[NohJetCal]
   Float_t         ohJetCalN90[82];   //[NohJetCal]
   Float_t         ohJetCalN90hits[82];   //[NohJetCal]
   Int_t           NohJetCorCal;
   Float_t         ohJetCorCalPt[80];   //[NohJetCorCal]
   Float_t         ohJetCorCalPhi[80];   //[NohJetCorCal]
   Float_t         ohJetCorCalEta[80];   //[NohJetCorCal]
   Float_t         ohJetCorCalE[80];   //[NohJetCorCal]
   Float_t         ohJetCorCalEMF[80];   //[NohJetCorCal]
   Float_t         ohJetCorCalN90[80];   //[NohJetCorCal]
   Float_t         ohJetCorCalN90hits[80];   //[NohJetCorCal]
   Float_t         recoJetGenPt[1];   //[NrecoJetGen]
   Float_t         recoJetGenPhi[1];   //[NrecoJetGen]
   Float_t         recoJetGenEta[1];   //[NrecoJetGen]
   Float_t         recoJetGenE[1];   //[NrecoJetGen]
   Float_t         recoTowEt[665];   //[NrecoTowCal]
   Float_t         recoTowEta[665];   //[NrecoTowCal]
   Float_t         recoTowPhi[665];   //[NrecoTowCal]
   Float_t         recoTowE[665];   //[NrecoTowCal]
   Float_t         recoTowEm[665];   //[NrecoTowCal]
   Float_t         recoTowHad[665];   //[NrecoTowCal]
   Float_t         recoTowOE[665];   //[NrecoTowCal]
   Float_t         recoMetCal;
   Float_t         recoMetCalPhi;
   Float_t         recoMetCalSum;
   Float_t         recoMetGen;
   Float_t         recoMetGenPhi;
   Float_t         recoMetGenSum;
   Float_t         recoHTCal;
   Float_t         recoHTCalPhi;
   Float_t         recoHTCalSum;
   Float_t         recoMetPF;
   Float_t         recoMetPFSum;
   Float_t         recoMetPFPhi;
   Int_t           NohTau;
   Float_t         ohTauEta[18];   //[NohTau]
   Float_t         ohTauPhi[18];   //[NohTau]
   Float_t         ohTauPt[18];   //[NohTau]
   Float_t         ohTauEiso[18];   //[NohTau]
   Float_t         ohTauL25Tpt[18];   //[NohTau]
   Int_t           ohTauL3Tiso[18];   //[NohTau]
   Int_t           NohpfTau;
   Float_t         ohpfTauPt[137];   //[NohpfTau]
   Int_t           ohpfTauProngs[137];   //[NohpfTau]
   Float_t         ohpfTauEta[137];   //[NohpfTau]
   Float_t         ohpfTauPhi[137];   //[NohpfTau]
   Float_t         ohpfTauLeadTrackPt[137];   //[NohpfTau]
   Float_t         ohpfTauLeadPionPt[137];   //[NohpfTau]
   Float_t         ohpfTauTrkIso[137];   //[NohpfTau]
   Float_t         ohpfTauGammaIso[137];   //[NohpfTau]
   Float_t         ohpfTauJetPt[137];   //[NohpfTau]
   Int_t           NohpfTauTightCone;
   Float_t         ohpfTauTightConePt[137];   //[NohpfTauTightCone]
   Int_t           ohpfTauTightConeProngs[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConeEta[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConePhi[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConeLeadTrackPt[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConeLeadPionPt[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConeTrkIso[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConeGammaIso[137];   //[NohpfTauTightCone]
   Float_t         ohpfTauTightConeJetPt[137];   //[NohpfTauTightCone]
   Int_t           NRecoPFTau;
   Float_t         recopfTauPt[1];   //[NRecoPFTau]
   Float_t         recopfTauEta[1];   //[NRecoPFTau]
   Float_t         recopfTauPhi[1];   //[NRecoPFTau]
   Float_t         recopfTauLeadTrackPt[1];   //[NRecoPFTau]
   Float_t         recopfTauLeadPionPt[1];   //[NRecoPFTau]
   Int_t           recopfTauTrkIso[1];   //[NRecoPFTau]
   Int_t           recopfTauGammaIso[1];   //[NRecoPFTau]
   Float_t         recopfTauJetPt[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrByTancOnePercent[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrByTancHalfPercent[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrByTancQuarterPercent[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrByTancTenthPercent[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrByIso[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrAgainstMuon[1];   //[NRecoPFTau]
   Float_t         recopfTauDiscrAgainstElec[1];   //[NRecoPFTau]
   Float_t         pfHT;
   Float_t         pfMHT;
   Int_t           NohPFJet;
   Float_t         pfJetPt[137];   //[NohPFJet]
   Float_t         pfJetE[137];   //[NohPFJet]
   Float_t         pfJetEta[137];   //[NohPFJet]
   Float_t         pfJetPhi[137];   //[NohPFJet]
   Int_t           nrpj;
   Float_t         recopfJetpt[1];   //[nrpj]
   Float_t         recopfJete[1];   //[nrpj]
   Float_t         recopfJetphi[1];   //[nrpj]
   Float_t         recopfJeteta[1];   //[nrpj]
   Float_t         recopfJetneutralHadronFraction[1];   //[nrpj]
   Float_t         recopfJetneutralEMFraction[1];   //[nrpj]
   Float_t         recopfJetchargedHadronFraction[1];   //[nrpj]
   Float_t         recopfJetchargedEMFraction[1];   //[nrpj]
   Int_t           recopfJetneutralMultiplicity[1];   //[nrpj]
   Int_t           recopfJetchargedMultiplicity[1];   //[nrpj]
   Int_t           NohBJetL2;
   Float_t         ohBJetL2Energy[10];   //[NohBJetL2]
   Float_t         ohBJetL2Et[10];   //[NohBJetL2]
   Float_t         ohBJetL2Pt[10];   //[NohBJetL2]
   Float_t         ohBJetL2Eta[10];   //[NohBJetL2]
   Float_t         ohBJetL2Phi[10];   //[NohBJetL2]
   Int_t           NohBJetL2Corrected;
   Float_t         ohBJetL2CorrectedEnergy[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetL2CorrectedEt[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetL2CorrectedPt[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetL2CorrectedEta[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetL2CorrectedPhi[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetIPL25Tag[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetIPL3Tag[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetIPL25TagSingleTrack[10];   //[NohBJetL2Corrected]
   Float_t         ohBJetIPL3TagSingleTrack[10];   //[NohBJetL2Corrected]
   Int_t           ohBJetPerfL25Tag[10];   //[NohBJetL2Corrected]
   Int_t           ohBJetPerfL3Tag[10];   //[NohBJetL2Corrected]
   Int_t           NrecoElec;
   Float_t         recoElecPt[1];   //[NrecoElec]
   Float_t         recoElecPhi[1];   //[NrecoElec]
   Float_t         recoElecEta[1];   //[NrecoElec]
   Float_t         recoElecEt[1];   //[NrecoElec]
   Float_t         recoElecE[1];   //[NrecoElec]
   Int_t           recoElecEleID[1];   //[NrecoElec]
   Float_t         recoElecIP[1];   //[NrecoElec]
   Int_t           recoElecNLostHits[1];   //[NrecoElec]
   Float_t         recoElecChi2NDF[1];   //[NrecoElec]
   Float_t         recoElecTrkIsoR03[1];   //[NrecoElec]
   Float_t         recoElecECaloIsoR03[1];   //[NrecoElec]
   Float_t         recoElecHCaloIsoR03[1];   //[NrecoElec]
   Bool_t          recoElecIsEcalDriven[1];   //[NrecoElec]
   Float_t         recoElecFbrem[1];   //[NrecoElec]
   Int_t           recoElecmishits[1];   //[NrecoElec]
   Float_t         recoElecdist[1];   //[NrecoElec]
   Float_t         recoElecdcot[1];   //[NrecoElec]
   Float_t         recoElectrkiso[1];   //[NrecoElec]
   Float_t         recoElececaliso[1];   //[NrecoElec]
   Float_t         recoElechcaliso[1];   //[NrecoElec]
   Float_t         recoElecsigmaietaieta[1];   //[NrecoElec]
   Float_t         recoElecdeltaPhiIn[1];   //[NrecoElec]
   Float_t         recoElecdeltaEtaIn[1];   //[NrecoElec]
   Float_t         recoElechOverE[1];   //[NrecoElec]
   Float_t         recoElecscEt[1];   //[NrecoElec]
   Float_t         recoElecd0corr[1];   //[NrecoElec]
   Bool_t          recoElecqGsfCtfScPixConsistent[1];   //[NrecoElec]
   Int_t           NrecoPhot;
   Float_t         recoPhotPt[1];   //[NrecoPhot]
   Float_t         recoPhotPhi[1];   //[NrecoPhot]
   Float_t         recoPhotEta[1];   //[NrecoPhot]
   Float_t         recoPhotEt[1];   //[NrecoPhot]
   Float_t         recoPhotE[1];   //[NrecoPhot]
   Float_t         recoPhotTiso[1];   //[NrecoPhot]
   Float_t         recoPhotEiso[1];   //[NrecoPhot]
   Float_t         recoPhotHiso[1];   //[NrecoPhot]
   Float_t         recoPhotHoverE[1];   //[NrecoPhot]
   Float_t         recoPhotClusShap[1];   //[NrecoPhot]
   Float_t         recoPhotR9ID[1];   //[NrecoPhot]
   Int_t           NohPhot;
   Float_t         ohPhotEt[5];   //[NohPhot]
   Float_t         ohPhotEta[5];   //[NohPhot]
   Float_t         ohPhotPhi[5];   //[NohPhot]
   Float_t         ohPhotEiso[5];   //[NohPhot]
   Float_t         ohPhotHiso[5];   //[NohPhot]
   Float_t         ohPhotTiso[5];   //[NohPhot]
   Int_t           ohPhotL1iso[5];   //[NohPhot]
   Float_t         ohPhotClusShap[5];   //[NohPhot]
   Float_t         ohPhotR9[5];   //[NohPhot]
   Float_t         ohPhotHforHoverE[5];   //[NohPhot]
   Float_t         ohPhotR9ID[5];   //[NohPhot]
   Int_t           NohEcalActiv;
   Float_t         ohEcalActivEt[7];   //[NohEcalActiv]
   Float_t         ohEcalActivEta[7];   //[NohEcalActiv]
   Float_t         ohEcalActivPhi[7];   //[NohEcalActiv]
   Float_t         ohEcalActivEiso[7];   //[NohEcalActiv]
   Float_t         ohEcalActivHiso[7];   //[NohEcalActiv]
   Float_t         ohEcalActivTiso[7];   //[NohEcalActiv]
   Int_t           ohEcalActivL1iso[7];   //[NohEcalActiv]
   Float_t         ohEcalActivClusShap[7];   //[NohEcalActiv]
   Float_t         ohEcalActivR9[7];   //[NohEcalActiv]
   Float_t         ohEcalActivHforHoverE[7];   //[NohEcalActiv]
   Float_t         ohEcalActivR9ID[7];   //[NohEcalActiv]
   Int_t           NohEle;
   Float_t         ohEleEt[6];   //[NohEle]
   Float_t         ohEleEta[6];   //[NohEle]
   Float_t         ohElePhi[6];   //[NohEle]
   Float_t         ohEleVtxZ[6];   //[NohEle]
   Float_t         ohEleE[6];   //[NohEle]
   Float_t         ohEleP[6];   //[NohEle]
   Float_t         ohEleHiso[6];   //[NohEle]
   Float_t         ohEleTiso[6];   //[NohEle]
   Float_t         ohEleEiso[6];   //[NohEle]
   Int_t           ohEleL1iso[6];   //[NohEle]
   Int_t           ohElePixelSeeds[6];   //[NohEle]
   Int_t           ohEleNewSC[6];   //[NohEle]
   Float_t         ohEleClusShap[6];   //[NohEle]
   Float_t         ohEleDeta[6];   //[NohEle]
   Float_t         ohEleDphi[6];   //[NohEle]
   Float_t         ohEleR9[6];   //[NohEle]
   Float_t         ohEleHforHoverE[6];   //[NohEle]
   Float_t         ohEleR9ID[6];   //[NohEle]
   Int_t           NohHFEle;
   Float_t         ohHFElePt[1];   //[NohHFEle]
   Float_t         ohHFEleEta[1];   //[NohHFEle]
   Int_t           NohHFECALClus;
   Float_t         ohHFEleClustere9e25[1];   //[NohHFECALClus]
   Float_t         ohHFEleClustere1e9[1];   //[NohHFECALClus]
   Float_t         ohHFEleClustereCOREe9[1];   //[NohHFECALClus]
   Float_t         ohHFEleClustereSeL[1];   //[NohHFECALClus]
   Float_t         ohHFEleCluster2Dcut[1];   //[NohHFECALClus]
   Float_t         ohHFEleClusterEta[1];   //[NohHFECALClus]
   Float_t         ohHFEleClusterPhi[1];   //[NohHFECALClus]
   Int_t           NrecoMuon;
   Float_t         recoMuonPt[1];   //[NrecoMuon]
   Float_t         recoMuonPhi[1];   //[NrecoMuon]
   Float_t         recoMuonEta[1];   //[NrecoMuon]
   Float_t         recoMuonEt[1];   //[NrecoMuon]
   Float_t         recoMuonE[1];   //[NrecoMuon]
   Float_t         recoMuonChi2NDF[1];   //[NrecoMuon]
   Float_t         recoMuonCharge[1];   //[NrecoMuon]
   Float_t         recoMuonTrkIsoR03[1];   //[NrecoMuon]
   Float_t         recoMuonECalIsoR03[1];   //[NrecoMuon]
   Float_t         recoMuonHCalIsoR03[1];   //[NrecoMuon]
   Float_t         recoMuonD0[1];   //[NrecoMuon]
   Int_t           recoMuonType[1];   //[NrecoMuon]
   Int_t           recoMuonNValidTrkHits[1];   //[NrecoMuon]
   Int_t           recoMuonNValidMuonHits[1];   //[NrecoMuon]
   Int_t           NohMuL2;
   Float_t         ohMuL2Pt[3];   //[NohMuL2]
   Float_t         ohMuL2Phi[3];   //[NohMuL2]
   Float_t         ohMuL2Eta[3];   //[NohMuL2]
   Int_t           ohMuL2Chg[3];   //[NohMuL2]
   Float_t         ohMuL2PtErr[3];   //[NohMuL2]
   Int_t           ohMuL2Iso[3];   //[NohMuL2]
   Float_t         ohMuL2Dr[3];   //[NohMuL2]
   Float_t         ohMuL2Dz[3];   //[NohMuL2]
   Float_t         ohMuL2VtxZ[3];   //[NohMuL2]
   Int_t           ohMuL2Nhits[3];   //[NohMuL2]
   Int_t           ohMuL2Nchambers[3];   //[NohMuL2]
   Int_t           ohMuL2Nstat[3];   //[NohMuL2]
   Int_t           ohMuL2L1idx[3];   //[NohMuL2]
   Int_t           NohMuL3;
   Float_t         ohMuL3Pt[3];   //[NohMuL3]
   Float_t         ohMuL3Phi[3];   //[NohMuL3]
   Float_t         ohMuL3Eta[3];   //[NohMuL3]
   Int_t           ohMuL3Chg[3];   //[NohMuL3]
   Float_t         ohMuL3PtErr[3];   //[NohMuL3]
   Int_t           ohMuL3Iso[3];   //[NohMuL3]
   Int_t           ohMuL3Trk10Iso[3];   //[NohMuL3]
   Float_t         ohMuL3Dr[3];   //[NohMuL3]
   Float_t         ohMuL3Dz[3];   //[NohMuL3]
   Float_t         ohMuL3VtxZ[3];   //[NohMuL3]
   Int_t           ohMuL3Nhits[3];   //[NohMuL3]
   Float_t         ohMuL3NormChi2[3];   //[NohMuL3]
   Int_t           ohMuL3Ntrackerhits[3];   //[NohMuL3]
   Int_t           ohMuL3Nmuonhits[3];   //[NohMuL3]
   Int_t           ohMuL3L2idx[3];   //[NohMuL3]
   Int_t           NohOniaPixel;
   Float_t         ohOniaPixelPt[49];   //[NohOniaPixel]
   Float_t         ohOniaPixelPhi[49];   //[NohOniaPixel]
   Float_t         ohOniaPixelEta[49];   //[NohOniaPixel]
   Int_t           ohOniaPixelChg[49];   //[NohOniaPixel]
   Float_t         ohOniaPixelDr[49];   //[NohOniaPixel]
   Float_t         ohOniaPixelDz[49];   //[NohOniaPixel]
   Int_t           ohOniaPixelHits[49];   //[NohOniaPixel]
   Float_t         ohOniaPixelNormChi2[49];   //[NohOniaPixel]
   Int_t           NohOniaTrack;
   Float_t         ohOniaTrackPt[21];   //[NohOniaTrack]
   Float_t         ohOniaTrackPhi[21];   //[NohOniaTrack]
   Float_t         ohOniaTrackEta[21];   //[NohOniaTrack]
   Int_t           ohOniaTrackChg[21];   //[NohOniaTrack]
   Float_t         ohOniaTrackDr[21];   //[NohOniaTrack]
   Float_t         ohOniaTrackDz[21];   //[NohOniaTrack]
   Int_t           ohOniaTrackHits[21];   //[NohOniaTrack]
   Float_t         ohOniaTrackNormChi2[21];   //[NohOniaTrack]
   Int_t           NohMuL2NoVtx;
   Float_t         ohMuL2NoVtxPt[3];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxPhi[3];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxEta[3];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxChg[3];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxPtErr[3];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxDr[3];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxDz[3];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxNhits[3];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxNchambers[3];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxL1idx[3];   //[NohMuL2NoVtx]
   Int_t           NohDiMu;
   Float_t         ohDiMuDCA[3];   //[NohDiMu]
   Int_t           ohDiMu1st[3];   //[NohDiMu]
   Int_t           ohDiMu2nd[3];   //[NohDiMu]
   Int_t           NohDiMuVtx;
   Int_t           ohDiMuVtx1st[1];   //[NohDiMuVtx]
   Int_t           ohDiMuVtx2nd[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxChi2[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxR[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxRSig[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxROverSig[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxCosAlpha[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxMu2DIpMax[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxMu2DIpMin[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxMu2DIpSigMax[1];   //[NohDiMuVtx]
   Float_t         ohDiMuVtxMu2DIpSigMin[1];   //[NohDiMuVtx]
   Int_t           NohIsoPixelTrackL3;
   Float_t         ohIsoPixelTrackL3Pt[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL3Eta[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL3Phi[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL3MaxPtPxl[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL3Energy[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL2pt[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL2eta[1];   //[NohIsoPixelTrackL3]
   Float_t         ohIsoPixelTrackL2dXY[1];   //[NohIsoPixelTrackL3]
   Int_t           NohPixelTracksL3;
   Float_t         ohPixelTracksL3Pt[1];   //[NohPixelTracksL3]
   Float_t         ohPixelTracksL3Eta[1];   //[NohPixelTracksL3]
   Float_t         ohPixelTracksL3Phi[1];   //[NohPixelTracksL3]
   Float_t         ohPixelTracksL3Vz[1];   //[NohPixelTracksL3]
   Int_t           ohPixelFEDSize;
   Int_t           NohPixelClusters;
   Int_t           NMCpart;
   Int_t           MCpid[1];   //[NMCpart]
   Int_t           MCstatus[1];   //[NMCpart]
   Float_t         MCvtxX[1];   //[NMCpart]
   Float_t         MCvtxY[1];   //[NMCpart]
   Float_t         MCvtxZ[1];   //[NMCpart]
   Float_t         MCpt[1];   //[NMCpart]
   Float_t         MCeta[1];   //[NMCpart]
   Float_t         MCphi[1];   //[NMCpart]
   Float_t         MCPtHat;
   Int_t           MCmu3;
   Int_t           MCel3;
   Int_t           MCbb;
   Int_t           MCab;
   Int_t           MCWenu;
   Int_t           MCWmunu;
   Int_t           MCZee;
   Int_t           MCZmumu;
   Float_t         MCptEleMax;
   Float_t         MCptMuMax;
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[4];   //[NL1IsolEm]
   Float_t         L1IsolEmE[4];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[4];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[4];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[4];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[4];   //[NL1Mu]
   Float_t         L1MuE[4];   //[NL1Mu]
   Float_t         L1MuEta[4];   //[NL1Mu]
   Float_t         L1MuPhi[4];   //[NL1Mu]
   Int_t           L1MuIsol[4];   //[NL1Mu]
   Int_t           L1MuMip[4];   //[NL1Mu]
   Int_t           L1MuFor[4];   //[NL1Mu]
   Int_t           L1MuRPC[4];   //[NL1Mu]
   Int_t           L1MuQal[4];   //[NL1Mu]
   Int_t           L1MuChg[4];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[4];   //[NL1CenJet]
   Float_t         L1CenJetE[4];   //[NL1CenJet]
   Float_t         L1CenJetEta[4];   //[NL1CenJet]
   Float_t         L1CenJetPhi[4];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[4];   //[NL1ForJet]
   Float_t         L1ForJetE[4];   //[NL1ForJet]
   Float_t         L1ForJetEta[4];   //[NL1ForJet]
   Float_t         L1ForJetPhi[4];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[4];   //[NL1Tau]
   Float_t         L1TauE[4];   //[NL1Tau]
   Float_t         L1TauEta[4];   //[NL1Tau]
   Float_t         L1TauPhi[4];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1EtTot;
   Float_t         L1Mht;
   Float_t         L1MhtPhi;
   Float_t         L1EtHad;
   Int_t           L1HfRing1EtSumPositiveEta;
   Int_t           L1HfRing2EtSumPositiveEta;
   Int_t           L1HfRing1EtSumNegativeEta;
   Int_t           L1HfRing2EtSumNegativeEta;
   Int_t           L1HfTowerCountPositiveEtaRing1;
   Int_t           L1HfTowerCountNegativeEtaRing1;
   Int_t           L1HfTowerCountPositiveEtaRing2;
   Int_t           L1HfTowerCountNegativeEtaRing2;
   Int_t           recoNVrtOffline0;
   Float_t         recoVrtXOffline0[1];   //[NVrtx]
   Float_t         recoVrtYOffline0[1];   //[NVrtx]
   Float_t         recoVrtZOffline0[1];   //[NVrtx]
   Int_t           recoVrtNtrkOffline0[1];   //[NVrtx]
   Float_t         recoVrtChi2Offline0[1];   //[NVrtx]
   Float_t         recoVrtNdofOffline0[1];   //[NVrtx]
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiBlock;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           HLT_Activity_Ecal_SC7_v1;
   Int_t           HLT_Activity_Ecal_SC7_v1_Prescl;
   Int_t           HLT_L1SingleJet36_v1;
   Int_t           HLT_L1SingleJet36_v1_Prescl;
   Int_t           HLT_Jet30_v1;
   Int_t           HLT_Jet30_v1_Prescl;
   Int_t           HLT_Jet60_v1;
   Int_t           HLT_Jet60_v1_Prescl;
   Int_t           HLT_Jet80_v1;
   Int_t           HLT_Jet80_v1_Prescl;
   Int_t           HLT_Jet110_v1;
   Int_t           HLT_Jet110_v1_Prescl;
   Int_t           HLT_Jet150_v1;
   Int_t           HLT_Jet150_v1_Prescl;
   Int_t           HLT_Jet190_v1;
   Int_t           HLT_Jet190_v1_Prescl;
   Int_t           HLT_Jet240_v1;
   Int_t           HLT_Jet240_v1_Prescl;
   Int_t           HLT_Jet370_v1;
   Int_t           HLT_Jet370_v1_Prescl;
   Int_t           HLT_Jet370_NoJetID_v1;
   Int_t           HLT_Jet370_NoJetID_v1_Prescl;
   Int_t           HLT_DiJetAve15U_v4;
   Int_t           HLT_DiJetAve15U_v4_Prescl;
   Int_t           HLT_DiJetAve30U_v4;
   Int_t           HLT_DiJetAve30U_v4_Prescl;
   Int_t           HLT_DiJetAve50U_v4;
   Int_t           HLT_DiJetAve50U_v4_Prescl;
   Int_t           HLT_DiJetAve70U_v4;
   Int_t           HLT_DiJetAve70U_v4_Prescl;
   Int_t           HLT_DiJetAve100U_v4;
   Int_t           HLT_DiJetAve100U_v4_Prescl;
   Int_t           HLT_DiJetAve140U_v4;
   Int_t           HLT_DiJetAve140U_v4_Prescl;
   Int_t           HLT_DiJetAve180U_v4;
   Int_t           HLT_DiJetAve180U_v4_Prescl;
   Int_t           HLT_DiJetAve300U_v4;
   Int_t           HLT_DiJetAve300U_v4_Prescl;
   Int_t           HLT_DoubleJet30_ForwardBackward_v2;
   Int_t           HLT_DoubleJet30_ForwardBackward_v2_Prescl;
   Int_t           HLT_DoubleJet60_ForwardBackward_v2;
   Int_t           HLT_DoubleJet60_ForwardBackward_v2_Prescl;
   Int_t           HLT_DoubleJet70_ForwardBackward_v2;
   Int_t           HLT_DoubleJet70_ForwardBackward_v2_Prescl;
   Int_t           HLT_DoubleJet80_ForwardBackward_v2;
   Int_t           HLT_DoubleJet80_ForwardBackward_v2_Prescl;
   Int_t           HLT_CentralJet80_MET65_v1;
   Int_t           HLT_CentralJet80_MET65_v1_Prescl;
   Int_t           HLT_CentralJet80_MET80_v1;
   Int_t           HLT_CentralJet80_MET80_v1_Prescl;
   Int_t           HLT_CentralJet80_MET100_v1;
   Int_t           HLT_CentralJet80_MET100_v1_Prescl;
   Int_t           HLT_CentralJet80_MET160_v1;
   Int_t           HLT_CentralJet80_MET160_v1_Prescl;
   Int_t           HLT_DiJet60_MET45_v1;
   Int_t           HLT_DiJet60_MET45_v1_Prescl;
   Int_t           HLT_DiJet70_PT70_v1;
   Int_t           HLT_DiJet70_PT70_v1_Prescl;
   Int_t           HLT_DiJet100_PT100_v1;
   Int_t           HLT_DiJet100_PT100_v1_Prescl;
   Int_t           HLT_DiJet130_PT130_v1;
   Int_t           HLT_DiJet130_PT130_v1_Prescl;
   Int_t           HLT_QuadJet40_v2;
   Int_t           HLT_QuadJet40_v2_Prescl;
   Int_t           HLT_QuadJet40_IsoPFTau40_v1;
   Int_t           HLT_QuadJet40_IsoPFTau40_v1_Prescl;
   Int_t           HLT_QuadJet50_BTagIP_v1;
   Int_t           HLT_QuadJet50_BTagIP_v1_Prescl;
   Int_t           HLT_QuadJet50_Jet40_v1;
   Int_t           HLT_QuadJet50_Jet40_v1_Prescl;
   Int_t           HLT_QuadJet60_v1;
   Int_t           HLT_QuadJet60_v1_Prescl;
   Int_t           HLT_QuadJet70_v1;
   Int_t           HLT_QuadJet70_v1_Prescl;
   Int_t           HLT_ExclDiJet60_HFOR_v1;
   Int_t           HLT_ExclDiJet60_HFOR_v1_Prescl;
   Int_t           HLT_ExclDiJet60_HFAND_v1;
   Int_t           HLT_ExclDiJet60_HFAND_v1_Prescl;
   Int_t           HLT_JetE30_NoBPTX_v2;
   Int_t           HLT_JetE30_NoBPTX_v2_Prescl;
   Int_t           HLT_JetE30_NoBPTX_NoHalo_v4;
   Int_t           HLT_JetE30_NoBPTX_NoHalo_v4_Prescl;
   Int_t           HLT_JetE30_NoBPTX3BX_NoHalo_v4;
   Int_t           HLT_JetE30_NoBPTX3BX_NoHalo_v4_Prescl;
   Int_t           HLT_HT150_v2;
   Int_t           HLT_HT150_v2_Prescl;
   Int_t           HLT_HT150_AlphaT0p60_v1;
   Int_t           HLT_HT150_AlphaT0p60_v1_Prescl;
   Int_t           HLT_HT150_AlphaT0p70_v1;
   Int_t           HLT_HT150_AlphaT0p70_v1_Prescl;
   Int_t           HLT_HT200_v2;
   Int_t           HLT_HT200_v2_Prescl;
   Int_t           HLT_HT200_AlphaT0p60_v1;
   Int_t           HLT_HT200_AlphaT0p60_v1_Prescl;
   Int_t           HLT_HT200_AlphaT0p65_v1;
   Int_t           HLT_HT200_AlphaT0p65_v1_Prescl;
   Int_t           HLT_HT250_v2;
   Int_t           HLT_HT250_v2_Prescl;
   Int_t           HLT_HT250_AlphaT0p55_v1;
   Int_t           HLT_HT250_AlphaT0p55_v1_Prescl;
   Int_t           HLT_HT250_AlphaT0p62_v1;
   Int_t           HLT_HT250_AlphaT0p62_v1_Prescl;
   Int_t           HLT_HT250_DoubleDisplacedJet60_v1;
   Int_t           HLT_HT250_DoubleDisplacedJet60_v1_Prescl;
   Int_t           HLT_HT250_MHT60_v2;
   Int_t           HLT_HT250_MHT60_v2_Prescl;
   Int_t           HLT_HT300_v3;
   Int_t           HLT_HT300_v3_Prescl;
   Int_t           HLT_HT300_MHT75_v3;
   Int_t           HLT_HT300_MHT75_v3_Prescl;
   Int_t           HLT_HT300_AlphaT0p52_v1;
   Int_t           HLT_HT300_AlphaT0p52_v1_Prescl;
   Int_t           HLT_HT300_AlphaT0p54_v1;
   Int_t           HLT_HT300_AlphaT0p54_v1_Prescl;
   Int_t           HLT_HT350_v2;
   Int_t           HLT_HT350_v2_Prescl;
   Int_t           HLT_HT350_AlphaT0p51_v1;
   Int_t           HLT_HT350_AlphaT0p51_v1_Prescl;
   Int_t           HLT_HT350_AlphaT0p53_v1;
   Int_t           HLT_HT350_AlphaT0p53_v1_Prescl;
   Int_t           HLT_HT400_v2;
   Int_t           HLT_HT400_v2_Prescl;
   Int_t           HLT_HT400_AlphaT0p51_v1;
   Int_t           HLT_HT400_AlphaT0p51_v1_Prescl;
   Int_t           HLT_HT450_v2;
   Int_t           HLT_HT450_v2_Prescl;
   Int_t           HLT_HT500_v2;
   Int_t           HLT_HT500_v2_Prescl;
   Int_t           HLT_HT550_v2;
   Int_t           HLT_HT550_v2_Prescl;
   Int_t           HLT_PFMHT150_v2;
   Int_t           HLT_PFMHT150_v2_Prescl;
   Int_t           HLT_MET100_v1;
   Int_t           HLT_MET100_v1_Prescl;
   Int_t           HLT_MET120_v1;
   Int_t           HLT_MET120_v1_Prescl;
   Int_t           HLT_MET200_v1;
   Int_t           HLT_MET200_v1_Prescl;
   Int_t           HLT_Meff440_v2;
   Int_t           HLT_Meff440_v2_Prescl;
   Int_t           HLT_Meff520_v2;
   Int_t           HLT_Meff520_v2_Prescl;
   Int_t           HLT_Meff640_v2;
   Int_t           HLT_Meff640_v2_Prescl;
   Int_t           HLT_MR100_v1;
   Int_t           HLT_MR100_v1_Prescl;
   Int_t           HLT_R032_v1;
   Int_t           HLT_R032_v1_Prescl;
   Int_t           HLT_R032_MR100_v1;
   Int_t           HLT_R032_MR100_v1_Prescl;
   Int_t           HLT_R035_MR100_v1;
   Int_t           HLT_R035_MR100_v1_Prescl;
   Int_t           HLT_L1SingleMuOpen_v1;
   Int_t           HLT_L1SingleMuOpen_v1_Prescl;
   Int_t           HLT_L1SingleMuOpen_DT_v1;
   Int_t           HLT_L1SingleMuOpen_DT_v1_Prescl;
   Int_t           HLT_L1SingleMu10_v1;
   Int_t           HLT_L1SingleMu10_v1_Prescl;
   Int_t           HLT_L1SingleMu20_v1;
   Int_t           HLT_L1SingleMu20_v1_Prescl;
   Int_t           HLT_L1DoubleMu0_v1;
   Int_t           HLT_L1DoubleMu0_v1_Prescl;
   Int_t           HLT_L2Mu10_v1;
   Int_t           HLT_L2Mu10_v1_Prescl;
   Int_t           HLT_L2Mu20_v1;
   Int_t           HLT_L2Mu20_v1_Prescl;
   Int_t           HLT_L2DoubleMu0_v2;
   Int_t           HLT_L2DoubleMu0_v2_Prescl;
   Int_t           HLT_Mu3_v3;
   Int_t           HLT_Mu3_v3_Prescl;
   Int_t           HLT_Mu5_v3;
   Int_t           HLT_Mu5_v3_Prescl;
   Int_t           HLT_Mu8_v1;
   Int_t           HLT_Mu8_v1_Prescl;
   Int_t           HLT_Mu12_v1;
   Int_t           HLT_Mu12_v1_Prescl;
   Int_t           HLT_Mu15_v2;
   Int_t           HLT_Mu15_v2_Prescl;
   Int_t           HLT_Mu20_v1;
   Int_t           HLT_Mu20_v1_Prescl;
   Int_t           HLT_Mu24_v1;
   Int_t           HLT_Mu24_v1_Prescl;
   Int_t           HLT_Mu30_v1;
   Int_t           HLT_Mu30_v1_Prescl;
   Int_t           HLT_IsoMu12_v1;
   Int_t           HLT_IsoMu12_v1_Prescl;
   Int_t           HLT_IsoMu15_v5;
   Int_t           HLT_IsoMu15_v5_Prescl;
   Int_t           HLT_IsoMu17_v5;
   Int_t           HLT_IsoMu17_v5_Prescl;
   Int_t           HLT_IsoMu24_v1;
   Int_t           HLT_IsoMu24_v1_Prescl;
   Int_t           HLT_IsoMu30_v1;
   Int_t           HLT_IsoMu30_v1_Prescl;
   Int_t           HLT_L2DoubleMu23_NoVertex_v1;
   Int_t           HLT_L2DoubleMu23_NoVertex_v1_Prescl;
   Int_t           HLT_DoubleMu3_v3;
   Int_t           HLT_DoubleMu3_v3_Prescl;
   Int_t           HLT_DoubleMu6_v1;
   Int_t           HLT_DoubleMu6_v1_Prescl;
   Int_t           HLT_DoubleMu7_v1;
   Int_t           HLT_DoubleMu7_v1_Prescl;
   Int_t           HLT_DoubleMu2_Bs_v1;
   Int_t           HLT_DoubleMu2_Bs_v1_Prescl;
   Int_t           HLT_DoubleMu3_Jpsi_v2;
   Int_t           HLT_DoubleMu3_Jpsi_v2_Prescl;
   Int_t           HLT_DoubleMu3_Quarkonium_v2;
   Int_t           HLT_DoubleMu3_Quarkonium_v2_Prescl;
   Int_t           HLT_DoubleMu3_Upsilon_v1;
   Int_t           HLT_DoubleMu3_Upsilon_v1_Prescl;
   Int_t           HLT_DoubleMu3_LowMass_v1;
   Int_t           HLT_DoubleMu3_LowMass_v1_Prescl;
   Int_t           HLT_DoubleMu4_Acoplanarity03_v1;
   Int_t           HLT_DoubleMu4_Acoplanarity03_v1_Prescl;
   Int_t           HLT_TripleMu5_v2;
   Int_t           HLT_TripleMu5_v2_Prescl;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1_Prescl;
   Int_t           HLT_Mu5_L2Mu2_v2;
   Int_t           HLT_Mu5_L2Mu2_v2_Prescl;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v2;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v2_Prescl;
   Int_t           HLT_Mu3_Track3_Jpsi_v5;
   Int_t           HLT_Mu3_Track3_Jpsi_v5_Prescl;
   Int_t           HLT_Mu5_Track2_Jpsi_v1;
   Int_t           HLT_Mu5_Track2_Jpsi_v1_Prescl;
   Int_t           HLT_Mu7_Track5_Jpsi_v2;
   Int_t           HLT_Mu7_Track5_Jpsi_v2_Prescl;
   Int_t           HLT_Mu7_Track7_Jpsi_v2;
   Int_t           HLT_Mu7_Track7_Jpsi_v2_Prescl;
   Int_t           HLT_Photon20_CaloIdVL_IsoL_v1;
   Int_t           HLT_Photon20_CaloIdVL_IsoL_v1_Prescl;
   Int_t           HLT_Photon20_R9Id_Photon18_R9Id_v2;
   Int_t           HLT_Photon20_R9Id_Photon18_R9Id_v2_Prescl;
   Int_t           HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2;
   Int_t           HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl;
   Int_t           HLT_Photon26_Photon18_v2;
   Int_t           HLT_Photon26_Photon18_v2_Prescl;
   Int_t           HLT_Photon26_IsoVL_Photon18_v2;
   Int_t           HLT_Photon26_IsoVL_Photon18_v2_Prescl;
   Int_t           HLT_Photon26_IsoVL_Photon18_IsoVL_v2;
   Int_t           HLT_Photon26_IsoVL_Photon18_IsoVL_v2_Prescl;
   Int_t           HLT_Photon26_CaloIdL_IsoVL_Photon18_v2;
   Int_t           HLT_Photon26_CaloIdL_IsoVL_Photon18_v2_Prescl;
   Int_t           HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1;
   Int_t           HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1_Prescl;
   Int_t           HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2;
   Int_t           HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2_Prescl;
   Int_t           HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1;
   Int_t           HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1_Prescl;
   Int_t           HLT_Photon30_CaloIdVL_v2;
   Int_t           HLT_Photon30_CaloIdVL_v2_Prescl;
   Int_t           HLT_Photon30_CaloIdVL_IsoL_v2;
   Int_t           HLT_Photon30_CaloIdVL_IsoL_v2_Prescl;
   Int_t           HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2;
   Int_t           HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2_Prescl;
   Int_t           HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1;
   Int_t           HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1_Prescl;
   Int_t           HLT_Photon50_CaloIdVL_IsoL_v1;
   Int_t           HLT_Photon50_CaloIdVL_IsoL_v1_Prescl;
   Int_t           HLT_Photon60_CaloIdL_HT200_v2;
   Int_t           HLT_Photon60_CaloIdL_HT200_v2_Prescl;
   Int_t           HLT_Photon70_CaloIdL_HT200_v2;
   Int_t           HLT_Photon70_CaloIdL_HT200_v2_Prescl;
   Int_t           HLT_Photon70_CaloIdL_HT300_v2;
   Int_t           HLT_Photon70_CaloIdL_HT300_v2_Prescl;
   Int_t           HLT_Photon70_CaloIdL_MHT30_v2;
   Int_t           HLT_Photon70_CaloIdL_MHT30_v2_Prescl;
   Int_t           HLT_Photon70_CaloIdL_MHT50_v2;
   Int_t           HLT_Photon70_CaloIdL_MHT50_v2_Prescl;
   Int_t           HLT_Photon75_CaloIdVL_v2;
   Int_t           HLT_Photon75_CaloIdVL_v2_Prescl;
   Int_t           HLT_Photon75_CaloIdVL_IsoL_v2;
   Int_t           HLT_Photon75_CaloIdVL_IsoL_v2_Prescl;
   Int_t           HLT_Photon125_NoSpikeFilter_v2;
   Int_t           HLT_Photon125_NoSpikeFilter_v2_Prescl;
   Int_t           HLT_DoublePhoton33_v2;
   Int_t           HLT_DoublePhoton33_v2_Prescl;
   Int_t           HLT_DoublePhoton5_IsoVL_CEP_v1;
   Int_t           HLT_DoublePhoton5_IsoVL_CEP_v1_Prescl;
   Int_t           HLT_L1SingleEG5_v1;
   Int_t           HLT_L1SingleEG5_v1_Prescl;
   Int_t           HLT_L1SingleEG12_v1;
   Int_t           HLT_L1SingleEG12_v1_Prescl;
   Int_t           HLT_Ele8_v2;
   Int_t           HLT_Ele8_v2_Prescl;
   Int_t           HLT_Ele8_CaloIdL_CaloIsoVL_v2;
   Int_t           HLT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl;
   Int_t           HLT_Ele8_CaloIdL_TrkIdVL_v2;
   Int_t           HLT_Ele8_CaloIdL_TrkIdVL_v2_Prescl;
   Int_t           HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2;
   Int_t           HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl;
   Int_t           HLT_Ele17_CaloIdL_CaloIsoVL_v2;
   Int_t           HLT_Ele17_CaloIdL_CaloIsoVL_v2_Prescl;
   Int_t           HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2;
   Int_t           HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2_Prescl;
   Int_t           HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2;
   Int_t           HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2_Prescl;
   Int_t           HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2;
   Int_t           HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2_Prescl;
   Int_t           HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2;
   Int_t           HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2_Prescl;
   Int_t           HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2;
   Int_t           HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl;
   Int_t           HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1;
   Int_t           HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1_Prescl;
   Int_t           HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2;
   Int_t           HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2_Prescl;
   Int_t           HLT_Ele45_CaloIdVT_TrkIdT_v2;
   Int_t           HLT_Ele45_CaloIdVT_TrkIdT_v2_Prescl;
   Int_t           HLT_Ele90_NoSpikeFilter_v2;
   Int_t           HLT_Ele90_NoSpikeFilter_v2_Prescl;
   Int_t           HLT_IsoPFTau35_Trk20_MET45_v2;
   Int_t           HLT_IsoPFTau35_Trk20_MET45_v2_Prescl;
   Int_t           HLT_DoubleIsoPFTau20_Trk5_v2;
   Int_t           HLT_DoubleIsoPFTau20_Trk5_v2_Prescl;
   Int_t           HLT_BTagMu_DiJet20_Mu5_v2;
   Int_t           HLT_BTagMu_DiJet20_Mu5_v2_Prescl;
   Int_t           HLT_BTagMu_DiJet60_Mu7_v2;
   Int_t           HLT_BTagMu_DiJet60_Mu7_v2_Prescl;
   Int_t           HLT_BTagMu_DiJet80_Mu9_v2;
   Int_t           HLT_BTagMu_DiJet80_Mu9_v2_Prescl;
   Int_t           HLT_BTagMu_DiJet100_Mu9_v2;
   Int_t           HLT_BTagMu_DiJet100_Mu9_v2_Prescl;
   Int_t           HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3;
   Int_t           HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3_Prescl;
   Int_t           HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3;
   Int_t           HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3_Prescl;
   Int_t           HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3;
   Int_t           HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3_Prescl;
   Int_t           HLT_Mu5_DoubleEle8_v3;
   Int_t           HLT_Mu5_DoubleEle8_v3_Prescl;
   Int_t           HLT_Mu5_HT200_v4;
   Int_t           HLT_Mu5_HT200_v4_Prescl;
   Int_t           HLT_Mu8_HT200_v3;
   Int_t           HLT_Mu8_HT200_v3_Prescl;
   Int_t           HLT_Mu8_Ele17_CaloIdL_v2;
   Int_t           HLT_Mu8_Ele17_CaloIdL_v2_Prescl;
   Int_t           HLT_Mu8_Photon20_CaloIdVT_IsoT_v2;
   Int_t           HLT_Mu8_Photon20_CaloIdVT_IsoT_v2_Prescl;
   Int_t           HLT_Mu8_Jet40_v3;
   Int_t           HLT_Mu8_Jet40_v3_Prescl;
   Int_t           HLT_Mu10_Ele10_CaloIdL_v3;
   Int_t           HLT_Mu10_Ele10_CaloIdL_v3_Prescl;
   Int_t           HLT_Mu15_Photon20_CaloIdL_v3;
   Int_t           HLT_Mu15_Photon20_CaloIdL_v3_Prescl;
   Int_t           HLT_Mu15_DoublePhoton15_CaloIdL_v3;
   Int_t           HLT_Mu15_DoublePhoton15_CaloIdL_v3_Prescl;
   Int_t           HLT_Mu15_LooseIsoPFTau20_v2;
   Int_t           HLT_Mu15_LooseIsoPFTau20_v2_Prescl;
   Int_t           HLT_Mu17_CentralJet30_v2;
   Int_t           HLT_Mu17_CentralJet30_v2_Prescl;
   Int_t           HLT_Mu17_DiCentralJet30_v2;
   Int_t           HLT_Mu17_DiCentralJet30_v2_Prescl;
   Int_t           HLT_Mu17_TriCentralJet30_v2;
   Int_t           HLT_Mu17_TriCentralJet30_v2_Prescl;
   Int_t           HLT_Mu17_Ele8_CaloIdL_v2;
   Int_t           HLT_Mu17_Ele8_CaloIdL_v2_Prescl;
   Int_t           HLT_Mu17_CentralJet40_BTagIP_v2;
   Int_t           HLT_Mu17_CentralJet40_BTagIP_v2_Prescl;
   Int_t           HLT_IsoMu12_LooseIsoPFTau10_v2;
   Int_t           HLT_IsoMu12_LooseIsoPFTau10_v2_Prescl;
   Int_t           HLT_IsoMu17_CentralJet40_BTagIP_v2;
   Int_t           HLT_IsoMu17_CentralJet40_BTagIP_v2_Prescl;
   Int_t           HLT_DoubleMu3_HT160_v3;
   Int_t           HLT_DoubleMu3_HT160_v3_Prescl;
   Int_t           HLT_DoubleMu3_HT200_v3;
   Int_t           HLT_DoubleMu3_HT200_v3_Prescl;
   Int_t           HLT_DoubleMu5_Ele8_v3;
   Int_t           HLT_DoubleMu5_Ele8_v3_Prescl;
   Int_t           HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3;
   Int_t           HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3_Prescl;
   Int_t           HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2;
   Int_t           HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2_Prescl;
   Int_t           HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3;
   Int_t           HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3_Prescl;
   Int_t           HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3;
   Int_t           HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3_Prescl;
   Int_t           HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2;
   Int_t           HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2_Prescl;
   Int_t           HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2;
   Int_t           HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2_Prescl;
   Int_t           HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2;
   Int_t           HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2_Prescl;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2_Prescl;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2_Prescl;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2_Prescl;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2;
   Int_t           HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2_Prescl;
   Int_t           HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3;
   Int_t           HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3_Prescl;
   Int_t           HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3;
   Int_t           HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3_Prescl;
   Int_t           HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2;
   Int_t           HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2_Prescl;
   Int_t           HLT_TripleEle10_CaloIdL_TrkIdVL_v2;
   Int_t           HLT_TripleEle10_CaloIdL_TrkIdVL_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity80_v2;
   Int_t           HLT_PixelTracks_Multiplicity80_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity100_v2;
   Int_t           HLT_PixelTracks_Multiplicity100_v2_Prescl;
   Int_t           HLT_BeamGas_HF_v2;
   Int_t           HLT_BeamGas_HF_v2_Prescl;
   Int_t           HLT_BeamGas_BSC_v2;
   Int_t           HLT_BeamGas_BSC_v2_Prescl;
   Int_t           HLT_BeamHalo_v2;
   Int_t           HLT_BeamHalo_v2_Prescl;
   Int_t           HLT_L1Tech_BSC_minBias_threshold1_v1;
   Int_t           HLT_L1Tech_BSC_minBias_threshold1_v1_Prescl;
   Int_t           HLT_L1Tech_BSC_halo_v1;
   Int_t           HLT_L1Tech_BSC_halo_v1_Prescl;
   Int_t           HLT_L1Tech_CASTOR_HaloMuon_v1;
   Int_t           HLT_L1Tech_CASTOR_HaloMuon_v1_Prescl;
   Int_t           HLT_L1_PreCollisions_v1;
   Int_t           HLT_L1_PreCollisions_v1_Prescl;
   Int_t           HLT_L1_Interbunch_BSC_v1;
   Int_t           HLT_L1_Interbunch_BSC_v1_Prescl;
   Int_t           HLT_IsoTrackHE_v3;
   Int_t           HLT_IsoTrackHE_v3_Prescl;
   Int_t           HLT_IsoTrackHB_v2;
   Int_t           HLT_IsoTrackHB_v2_Prescl;
   Int_t           HLT_HcalPhiSym_v3;
   Int_t           HLT_HcalPhiSym_v3_Prescl;
   Int_t           HLT_HcalNZS_v3;
   Int_t           HLT_HcalNZS_v3_Prescl;
   Int_t           HLT_GlobalRunHPDNoise_v2;
   Int_t           HLT_GlobalRunHPDNoise_v2_Prescl;
   Int_t           HLT_L1Tech_HBHEHO_totalOR_v1;
   Int_t           HLT_L1Tech_HBHEHO_totalOR_v1_Prescl;
   Int_t           HLT_ZeroBias_v1;
   Int_t           HLT_ZeroBias_v1_Prescl;
   Int_t           HLT_Physics_v1;
   Int_t           HLT_Physics_v1_Prescl;
   Int_t           HLT_Physics_NanoDST_v1;
   Int_t           HLT_Physics_NanoDST_v1_Prescl;
   Int_t           HLT_Calibration_v1;
   Int_t           HLT_Calibration_v1_Prescl;
   Int_t           HLT_EcalCalibration_v1;
   Int_t           HLT_EcalCalibration_v1_Prescl;
   Int_t           HLT_HcalCalibration_v1;
   Int_t           HLT_HcalCalibration_v1_Prescl;
   Int_t           HLT_TrackerCalibration_v1;
   Int_t           HLT_TrackerCalibration_v1_Prescl;
   Int_t           HLT_Random_v1;
   Int_t           HLT_Random_v1_Prescl;
   Int_t           HLT_L1SingleMuOpen_AntiBPTX_v1;
   Int_t           HLT_L1SingleMuOpen_AntiBPTX_v1_Prescl;
   Int_t           HLT_L1TrackerCosmics_v2;
   Int_t           HLT_L1TrackerCosmics_v2_Prescl;
   Int_t           HLT_RegionalCosmicTracking_v1;
   Int_t           HLT_RegionalCosmicTracking_v1_Prescl;
   Int_t           HLT_L3MuonsCosmicTracking_v1;
   Int_t           HLT_L3MuonsCosmicTracking_v1_Prescl;
   Int_t           HLT_LogMonitor_v1;
   Int_t           HLT_LogMonitor_v1_Prescl;
   Int_t           HLT_DTErrors_v1;
   Int_t           HLT_DTErrors_v1_Prescl;
   Int_t           AlCa_EcalPi0_v4;
   Int_t           AlCa_EcalPi0_v4_Prescl;
   Int_t           AlCa_EcalEta_v3;
   Int_t           AlCa_EcalEta_v3_Prescl;
   Int_t           AlCa_EcalPhiSym_v2;
   Int_t           AlCa_EcalPhiSym_v2_Prescl;
   Int_t           AlCa_RPCMuonNoTriggers_v2;
   Int_t           AlCa_RPCMuonNoTriggers_v2_Prescl;
   Int_t           AlCa_RPCMuonNoHits_v2;
   Int_t           AlCa_RPCMuonNoHits_v2_Prescl;
   Int_t           AlCa_RPCMuonNormalisation_v2;
   Int_t           AlCa_RPCMuonNormalisation_v2_Prescl;
   Int_t           DQM_FEDIntegrity_v3;
   Int_t           DQM_FEDIntegrity_v3_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           L1_BeamGas_Bsc;
   Int_t           L1_BeamGas_Bsc_Prescl;
   Int_t           L1_BeamGas_Bsc_5bx;
   Int_t           L1_BeamGas_Hf;
   Int_t           L1_BeamGas_Hf_Prescl;
   Int_t           L1_BeamGas_Hf_5bx;
   Int_t           L1_BeamHalo;
   Int_t           L1_BeamHalo_Prescl;
   Int_t           L1_BeamHalo_5bx;
   Int_t           L1_BptxMinus_NotBptxPlus;
   Int_t           L1_BptxMinus_NotBptxPlus_Prescl;
   Int_t           L1_BptxMinus_NotBptxPlus_5bx;
   Int_t           L1_BptxPlus_NotBptxMinus;
   Int_t           L1_BptxPlus_NotBptxMinus_Prescl;
   Int_t           L1_BptxPlus_NotBptxMinus_5bx;
   Int_t           L1_Bsc2Minus_BptxMinus;
   Int_t           L1_Bsc2Minus_BptxMinus_Prescl;
   Int_t           L1_Bsc2Minus_BptxMinus_5bx;
   Int_t           L1_Bsc2Plus_BptxPlus;
   Int_t           L1_Bsc2Plus_BptxPlus_Prescl;
   Int_t           L1_Bsc2Plus_BptxPlus_5bx;
   Int_t           L1_BscMinBiasOR_BptxPlusANDMinus;
   Int_t           L1_BscMinBiasOR_BptxPlusANDMinus_Prescl;
   Int_t           L1_BscMinBiasOR_BptxPlusANDMinus_5bx;
   Int_t           L1_DoubleEG10;
   Int_t           L1_DoubleEG10_Prescl;
   Int_t           L1_DoubleEG10_5bx;
   Int_t           L1_DoubleEG2_FwdVeto;
   Int_t           L1_DoubleEG2_FwdVeto_Prescl;
   Int_t           L1_DoubleEG2_FwdVeto_5bx;
   Int_t           L1_DoubleEG3;
   Int_t           L1_DoubleEG3_Prescl;
   Int_t           L1_DoubleEG3_5bx;
   Int_t           L1_DoubleEG5;
   Int_t           L1_DoubleEG5_Prescl;
   Int_t           L1_DoubleEG5_5bx;
   Int_t           L1_DoubleEG5_HTT50;
   Int_t           L1_DoubleEG5_HTT50_Prescl;
   Int_t           L1_DoubleEG5_HTT50_5bx;
   Int_t           L1_DoubleEG5_HTT75;
   Int_t           L1_DoubleEG5_HTT75_Prescl;
   Int_t           L1_DoubleEG5_HTT75_5bx;
   Int_t           L1_DoubleEG8;
   Int_t           L1_DoubleEG8_Prescl;
   Int_t           L1_DoubleEG8_5bx;
   Int_t           L1_DoubleEG_12_5;
   Int_t           L1_DoubleEG_12_5_Prescl;
   Int_t           L1_DoubleEG_12_5_5bx;
   Int_t           L1_DoubleEG_12_5_Eta1p39;
   Int_t           L1_DoubleEG_12_5_Eta1p39_Prescl;
   Int_t           L1_DoubleEG_12_5_Eta1p39_5bx;
   Int_t           L1_DoubleForJet32_EtaOpp;
   Int_t           L1_DoubleForJet32_EtaOpp_Prescl;
   Int_t           L1_DoubleForJet32_EtaOpp_5bx;
   Int_t           L1_DoubleForJet44_EtaOpp;
   Int_t           L1_DoubleForJet44_EtaOpp_Prescl;
   Int_t           L1_DoubleForJet44_EtaOpp_5bx;
   Int_t           L1_DoubleIsoEG10;
   Int_t           L1_DoubleIsoEG10_Prescl;
   Int_t           L1_DoubleIsoEG10_5bx;
   Int_t           L1_DoubleJet36_Central;
   Int_t           L1_DoubleJet36_Central_Prescl;
   Int_t           L1_DoubleJet36_Central_5bx;
   Int_t           L1_DoubleJet44_Central;
   Int_t           L1_DoubleJet44_Central_Prescl;
   Int_t           L1_DoubleJet44_Central_5bx;
   Int_t           L1_DoubleJet52;
   Int_t           L1_DoubleJet52_Prescl;
   Int_t           L1_DoubleJet52_5bx;
   Int_t           L1_DoubleMu0;
   Int_t           L1_DoubleMu0_Prescl;
   Int_t           L1_DoubleMu0_5bx;
   Int_t           L1_DoubleMu0_HighQ;
   Int_t           L1_DoubleMu0_HighQ_Prescl;
   Int_t           L1_DoubleMu0_HighQ_5bx;
   Int_t           L1_DoubleMu0_HighQ_EtaCuts;
   Int_t           L1_DoubleMu0_HighQ_EtaCuts_Prescl;
   Int_t           L1_DoubleMu0_HighQ_EtaCuts_5bx;
   Int_t           L1_DoubleMu3;
   Int_t           L1_DoubleMu3_Prescl;
   Int_t           L1_DoubleMu3_5bx;
   Int_t           L1_DoubleMu3_EG5;
   Int_t           L1_DoubleMu3_EG5_Prescl;
   Int_t           L1_DoubleMu3_EG5_5bx;
   Int_t           L1_DoubleMu3p5;
   Int_t           L1_DoubleMu3p5_Prescl;
   Int_t           L1_DoubleMu3p5_5bx;
   Int_t           L1_DoubleMu5_v1;
   Int_t           L1_DoubleMu5_v1_Prescl;
   Int_t           L1_DoubleMu5_v1_5bx;
   Int_t           L1_DoubleTauJet28;
   Int_t           L1_DoubleTauJet28_Prescl;
   Int_t           L1_DoubleTauJet28_5bx;
   Int_t           L1_DoubleTauJet32;
   Int_t           L1_DoubleTauJet32_Prescl;
   Int_t           L1_DoubleTauJet32_5bx;
   Int_t           L1_DoubleTauJet36;
   Int_t           L1_DoubleTauJet36_Prescl;
   Int_t           L1_DoubleTauJet36_5bx;
   Int_t           L1_DoubleTauJet40;
   Int_t           L1_DoubleTauJet40_Prescl;
   Int_t           L1_DoubleTauJet40_5bx;
   Int_t           L1_EG10_Jet24_Central_deltaPhi1;
   Int_t           L1_EG10_Jet24_Central_deltaPhi1_Prescl;
   Int_t           L1_EG10_Jet24_Central_deltaPhi1_5bx;
   Int_t           L1_EG12_Jet24_Central_deltaPhi1;
   Int_t           L1_EG12_Jet24_Central_deltaPhi1_Prescl;
   Int_t           L1_EG12_Jet24_Central_deltaPhi1_5bx;
   Int_t           L1_EG12_TauJet20_deltaPhi1;
   Int_t           L1_EG12_TauJet20_deltaPhi1_Prescl;
   Int_t           L1_EG12_TauJet20_deltaPhi1_5bx;
   Int_t           L1_EG5_HTT100;
   Int_t           L1_EG5_HTT100_Prescl;
   Int_t           L1_EG5_HTT100_5bx;
   Int_t           L1_EG5_HTT125;
   Int_t           L1_EG5_HTT125_Prescl;
   Int_t           L1_EG5_HTT125_5bx;
   Int_t           L1_EG5_HTT75;
   Int_t           L1_EG5_HTT75_Prescl;
   Int_t           L1_EG5_HTT75_5bx;
   Int_t           L1_EG5_Jet36_deltaPhi1;
   Int_t           L1_EG5_Jet36_deltaPhi1_Prescl;
   Int_t           L1_EG5_Jet36_deltaPhi1_5bx;
   Int_t           L1_EG8_Jet20_Central_deltaPhi1;
   Int_t           L1_EG8_Jet20_Central_deltaPhi1_Prescl;
   Int_t           L1_EG8_Jet20_Central_deltaPhi1_5bx;
   Int_t           L1_ETM100;
   Int_t           L1_ETM100_Prescl;
   Int_t           L1_ETM100_5bx;
   Int_t           L1_ETM20;
   Int_t           L1_ETM20_Prescl;
   Int_t           L1_ETM20_5bx;
   Int_t           L1_ETM30;
   Int_t           L1_ETM30_Prescl;
   Int_t           L1_ETM30_5bx;
   Int_t           L1_ETM50;
   Int_t           L1_ETM50_Prescl;
   Int_t           L1_ETM50_5bx;
   Int_t           L1_ETM70;
   Int_t           L1_ETM70_Prescl;
   Int_t           L1_ETM70_5bx;
   Int_t           L1_ETT220;
   Int_t           L1_ETT220_Prescl;
   Int_t           L1_ETT220_5bx;
   Int_t           L1_ETT260_EG5;
   Int_t           L1_ETT260_EG5_Prescl;
   Int_t           L1_ETT260_EG5_5bx;
   Int_t           L1_ETT300_EG5;
   Int_t           L1_ETT300_EG5_Prescl;
   Int_t           L1_ETT300_EG5_5bx;
   Int_t           L1_HTM50;
   Int_t           L1_HTM50_Prescl;
   Int_t           L1_HTM50_5bx;
   Int_t           L1_HTT100;
   Int_t           L1_HTT100_Prescl;
   Int_t           L1_HTT100_5bx;
   Int_t           L1_HTT150;
   Int_t           L1_HTT150_Prescl;
   Int_t           L1_HTT150_5bx;
   Int_t           L1_HTT50;
   Int_t           L1_HTT50_Prescl;
   Int_t           L1_HTT50_5bx;
   Int_t           L1_HTT50_HTM30;
   Int_t           L1_HTT50_HTM30_Prescl;
   Int_t           L1_HTT50_HTM30_5bx;
   Int_t           L1_HTT50_HTM50;
   Int_t           L1_HTT50_HTM50_Prescl;
   Int_t           L1_HTT50_HTM50_5bx;
   Int_t           L1_HTT75;
   Int_t           L1_HTT75_Prescl;
   Int_t           L1_HTT75_5bx;
   Int_t           L1_InterBunch_Bsc;
   Int_t           L1_InterBunch_Bsc_Prescl;
   Int_t           L1_InterBunch_Bsc_5bx;
   Int_t           L1_InterBunch_Hf;
   Int_t           L1_InterBunch_Hf_Prescl;
   Int_t           L1_InterBunch_Hf_5bx;
   Int_t           L1_Mu0_HTT50;
   Int_t           L1_Mu0_HTT50_Prescl;
   Int_t           L1_Mu0_HTT50_5bx;
   Int_t           L1_Mu0_HTT75;
   Int_t           L1_Mu0_HTT75_Prescl;
   Int_t           L1_Mu0_HTT75_5bx;
   Int_t           L1_Mu10_Jet36_Central;
   Int_t           L1_Mu10_Jet36_Central_Prescl;
   Int_t           L1_Mu10_Jet36_Central_5bx;
   Int_t           L1_Mu12_EG5;
   Int_t           L1_Mu12_EG5_Prescl;
   Int_t           L1_Mu12_EG5_5bx;
   Int_t           L1_Mu3_DoubleEG5;
   Int_t           L1_Mu3_DoubleEG5_Prescl;
   Int_t           L1_Mu3_DoubleEG5_5bx;
   Int_t           L1_Mu3_EG5;
   Int_t           L1_Mu3_EG5_Prescl;
   Int_t           L1_Mu3_EG5_5bx;
   Int_t           L1_Mu3_Jet16_Central;
   Int_t           L1_Mu3_Jet16_Central_Prescl;
   Int_t           L1_Mu3_Jet16_Central_5bx;
   Int_t           L1_Mu3_Jet20_Central;
   Int_t           L1_Mu3_Jet20_Central_Prescl;
   Int_t           L1_Mu3_Jet20_Central_5bx;
   Int_t           L1_Mu3_Jet28_Central;
   Int_t           L1_Mu3_Jet28_Central_Prescl;
   Int_t           L1_Mu3_Jet28_Central_5bx;
   Int_t           L1_Mu5_EG12;
   Int_t           L1_Mu5_EG12_Prescl;
   Int_t           L1_Mu5_EG12_5bx;
   Int_t           L1_Mu7_EG5;
   Int_t           L1_Mu7_EG5_Prescl;
   Int_t           L1_Mu7_EG5_5bx;
   Int_t           L1_Mu7_Jet20_Central;
   Int_t           L1_Mu7_Jet20_Central_Prescl;
   Int_t           L1_Mu7_Jet20_Central_5bx;
   Int_t           L1_Mu7_TauJet16;
   Int_t           L1_Mu7_TauJet16_Prescl;
   Int_t           L1_Mu7_TauJet16_5bx;
   Int_t           L1_MuOpen_EG12;
   Int_t           L1_MuOpen_EG12_Prescl;
   Int_t           L1_MuOpen_EG12_5bx;
   Int_t           L1_MuOpen_EG5;
   Int_t           L1_MuOpen_EG5_Prescl;
   Int_t           L1_MuOpen_EG5_5bx;
   Int_t           L1_PreCollisions;
   Int_t           L1_PreCollisions_Prescl;
   Int_t           L1_PreCollisions_5bx;
   Int_t           L1_QuadJet20_Central;
   Int_t           L1_QuadJet20_Central_Prescl;
   Int_t           L1_QuadJet20_Central_5bx;
   Int_t           L1_QuadJet28_Central;
   Int_t           L1_QuadJet28_Central_Prescl;
   Int_t           L1_QuadJet28_Central_5bx;
   Int_t           L1_SingleEG12;
   Int_t           L1_SingleEG12_Prescl;
   Int_t           L1_SingleEG12_5bx;
   Int_t           L1_SingleEG12_Eta1p39;
   Int_t           L1_SingleEG12_Eta1p39_Prescl;
   Int_t           L1_SingleEG12_Eta1p39_5bx;
   Int_t           L1_SingleEG12_Eta2p17;
   Int_t           L1_SingleEG12_Eta2p17_Prescl;
   Int_t           L1_SingleEG12_Eta2p17_5bx;
   Int_t           L1_SingleEG15;
   Int_t           L1_SingleEG15_Prescl;
   Int_t           L1_SingleEG15_5bx;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG20_5bx;
   Int_t           L1_SingleEG30;
   Int_t           L1_SingleEG30_Prescl;
   Int_t           L1_SingleEG30_5bx;
   Int_t           L1_SingleEG5;
   Int_t           L1_SingleEG5_Prescl;
   Int_t           L1_SingleEG5_5bx;
   Int_t           L1_SingleIsoEG12;
   Int_t           L1_SingleIsoEG12_Prescl;
   Int_t           L1_SingleIsoEG12_5bx;
   Int_t           L1_SingleIsoEG12_Eta1p39;
   Int_t           L1_SingleIsoEG12_Eta1p39_Prescl;
   Int_t           L1_SingleIsoEG12_Eta1p39_5bx;
   Int_t           L1_SingleIsoEG12_Eta2p17;
   Int_t           L1_SingleIsoEG12_Eta2p17_Prescl;
   Int_t           L1_SingleIsoEG12_Eta2p17_5bx;
   Int_t           L1_SingleJet128;
   Int_t           L1_SingleJet128_Prescl;
   Int_t           L1_SingleJet128_5bx;
   Int_t           L1_SingleJet16;
   Int_t           L1_SingleJet16_Prescl;
   Int_t           L1_SingleJet16_5bx;
   Int_t           L1_SingleJet20_NotBptxOR;
   Int_t           L1_SingleJet20_NotBptxOR_Prescl;
   Int_t           L1_SingleJet20_NotBptxOR_5bx;
   Int_t           L1_SingleJet20_NotBptxOR_NotMuBeamHalo;
   Int_t           L1_SingleJet20_NotBptxOR_NotMuBeamHalo_Prescl;
   Int_t           L1_SingleJet20_NotBptxOR_NotMuBeamHalo_5bx;
   Int_t           L1_SingleJet32_NotBptxOR_NotMuBeamHalo;
   Int_t           L1_SingleJet32_NotBptxOR_NotMuBeamHalo_Prescl;
   Int_t           L1_SingleJet32_NotBptxOR_NotMuBeamHalo_5bx;
   Int_t           L1_SingleJet36;
   Int_t           L1_SingleJet36_Prescl;
   Int_t           L1_SingleJet36_5bx;
   Int_t           L1_SingleJet36_FwdVeto;
   Int_t           L1_SingleJet36_FwdVeto_Prescl;
   Int_t           L1_SingleJet36_FwdVeto_5bx;
   Int_t           L1_SingleJet52;
   Int_t           L1_SingleJet52_Prescl;
   Int_t           L1_SingleJet52_5bx;
   Int_t           L1_SingleJet68;
   Int_t           L1_SingleJet68_Prescl;
   Int_t           L1_SingleJet68_5bx;
   Int_t           L1_SingleJet80_Central;
   Int_t           L1_SingleJet80_Central_Prescl;
   Int_t           L1_SingleJet80_Central_5bx;
   Int_t           L1_SingleJet92;
   Int_t           L1_SingleJet92_Prescl;
   Int_t           L1_SingleJet92_5bx;
   Int_t           L1_SingleJet92_Central;
   Int_t           L1_SingleJet92_Central_Prescl;
   Int_t           L1_SingleJet92_Central_5bx;
   Int_t           L1_SingleMu10;
   Int_t           L1_SingleMu10_Prescl;
   Int_t           L1_SingleMu10_5bx;
   Int_t           L1_SingleMu12;
   Int_t           L1_SingleMu12_Prescl;
   Int_t           L1_SingleMu12_5bx;
   Int_t           L1_SingleMu12_Debug;
   Int_t           L1_SingleMu12_Debug_Prescl;
   Int_t           L1_SingleMu12_Debug_5bx;
   Int_t           L1_SingleMu16;
   Int_t           L1_SingleMu16_Prescl;
   Int_t           L1_SingleMu16_5bx;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu20_Prescl;
   Int_t           L1_SingleMu20_5bx;
   Int_t           L1_SingleMu25;
   Int_t           L1_SingleMu25_Prescl;
   Int_t           L1_SingleMu25_5bx;
   Int_t           L1_SingleMu3;
   Int_t           L1_SingleMu3_Prescl;
   Int_t           L1_SingleMu3_5bx;
   Int_t           L1_SingleMu5_Eta1p5_Q80_v1;
   Int_t           L1_SingleMu5_Eta1p5_Q80_v1_Prescl;
   Int_t           L1_SingleMu5_Eta1p5_Q80_v1_5bx;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu7_Prescl;
   Int_t           L1_SingleMu7_5bx;
   Int_t           L1_SingleMu7_Barrel;
   Int_t           L1_SingleMu7_Barrel_Prescl;
   Int_t           L1_SingleMu7_Barrel_5bx;
   Int_t           L1_SingleMu7_Eta2p1;
   Int_t           L1_SingleMu7_Eta2p1_Prescl;
   Int_t           L1_SingleMu7_Eta2p1_5bx;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuBeamHalo_5bx;
   Int_t           L1_SingleMuOpen;
   Int_t           L1_SingleMuOpen_Prescl;
   Int_t           L1_SingleMuOpen_5bx;
   Int_t           L1_SingleTauJet52;
   Int_t           L1_SingleTauJet52_Prescl;
   Int_t           L1_SingleTauJet52_5bx;
   Int_t           L1_SingleTauJet68;
   Int_t           L1_SingleTauJet68_Prescl;
   Int_t           L1_SingleTauJet68_5bx;
   Int_t           L1_SingleTauJet80;
   Int_t           L1_SingleTauJet80_Prescl;
   Int_t           L1_SingleTauJet80_5bx;
   Int_t           L1_TripleEG5;
   Int_t           L1_TripleEG5_Prescl;
   Int_t           L1_TripleEG5_5bx;
   Int_t           L1_TripleEG7;
   Int_t           L1_TripleEG7_Prescl;
   Int_t           L1_TripleEG7_5bx;
   Int_t           L1_TripleEG_8_5_5;
   Int_t           L1_TripleEG_8_5_5_Prescl;
   Int_t           L1_TripleEG_8_5_5_5bx;
   Int_t           L1_TripleEG_8_8_5;
   Int_t           L1_TripleEG_8_8_5_Prescl;
   Int_t           L1_TripleEG_8_8_5_5bx;
   Int_t           L1_TripleJet28_Central;
   Int_t           L1_TripleJet28_Central_Prescl;
   Int_t           L1_TripleJet28_Central_5bx;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;
   Int_t           L1_ZeroBias_5bx;
   Int_t           L1Tech_BPTX_minus_v0;
   Int_t           L1Tech_BPTX_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_v0_5bx;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_5bx;
   Int_t           L1Tech_BPTX_plus_v0;
   Int_t           L1Tech_BPTX_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_v0_5bx;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_5bx;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_5bx;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_5bx;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_5bx;
   Int_t           L1Tech_BPTX_quiet_v0;
   Int_t           L1Tech_BPTX_quiet_v0_Prescl;
   Int_t           L1Tech_BPTX_quiet_v0_5bx;
   Int_t           L1Tech_BSC_HighMultiplicity_v0;
   Int_t           L1Tech_BSC_HighMultiplicity_v0_Prescl;
   Int_t           L1Tech_BSC_HighMultiplicity_v0_5bx;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0_5bx;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0_5bx;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0_5bx;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0_5bx;
   Int_t           L1Tech_BSC_minBias_OR_v0;
   Int_t           L1Tech_BSC_minBias_OR_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_OR_v0_5bx;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0_5bx;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0_5bx;
   Int_t           L1Tech_BSC_minBias_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold1_v0_5bx;
   Int_t           L1Tech_BSC_minBias_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold2_v0_5bx;
   Int_t           L1Tech_BSC_splash_beam1_v0;
   Int_t           L1Tech_BSC_splash_beam1_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam1_v0_5bx;
   Int_t           L1Tech_BSC_splash_beam2_v0;
   Int_t           L1Tech_BSC_splash_beam2_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam2_v0_5bx;
   Int_t           L1Tech_CASTOR_HaloMuon_v0;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_Prescl;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_5bx;
   Int_t           L1Tech_DT_GlobalOR_v0;
   Int_t           L1Tech_DT_GlobalOR_v0_Prescl;
   Int_t           L1Tech_DT_GlobalOR_v0_5bx;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_5bx;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v0;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v0_5bx;
   Int_t           L1Tech_HCAL_HF_MM_or_PP_or_PM_v0;
   Int_t           L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_5bx;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v1;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v1_Prescl;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v1_5bx;
   Int_t           L1Tech_HCAL_HO_totalOR_v0;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_5bx;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0_5bx;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0_5bx;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0_5bx;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_5bx;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_5bx;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0_5bx;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_5bx;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_5bx;
   Int_t           L1Tech_ZDC_Scint_loose_vertex_v0;
   Int_t           L1Tech_ZDC_Scint_loose_vertex_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_loose_vertex_v0_5bx;
   Int_t           L1Tech_ZDC_Scint_minus_v0;
   Int_t           L1Tech_ZDC_Scint_minus_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_minus_v0_5bx;
   Int_t           L1Tech_ZDC_Scint_plus_v0;
   Int_t           L1Tech_ZDC_Scint_plus_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_plus_v0_5bx;
   Int_t           L1Tech_ZDC_Scint_tight_vertex_v0;
   Int_t           L1Tech_ZDC_Scint_tight_vertex_v0_Prescl;
   Int_t           L1Tech_ZDC_Scint_tight_vertex_v0_5bx;

   // List of branches
   TBranch        *b_NoRecoPFTausSignal;   //!
   TBranch        *b_signalTrToPFTauMatch;   //!
   TBranch        *b_recoPFTauSignalTrDz;   //!
   TBranch        *b_recoPFTauSignalTrPt;   //!
   TBranch        *b_NoRecoPFTausIso;   //!
   TBranch        *b_isoTrToPFTauMatch;   //!
   TBranch        *b_recoPFTauIsoTrDz;   //!
   TBranch        *b_recoPFTauIsoTrPt;   //!
   TBranch        *b_NoHLTPFTausSignal;   //!
   TBranch        *b_hltpftauSignalTrToPFTauMatch;   //!
   TBranch        *b_HLTPFTauSignalTrDz;   //!
   TBranch        *b_HLTPFTauSignalTrPt;   //!
   TBranch        *b_NoHLTPFTausIso;   //!
   TBranch        *b_hltpftauIsoTrToPFTauMatch;   //!
   TBranch        *b_HLTPFTauIsoTrDz;   //!
   TBranch        *b_HLTPFTauIsoTrPt;   //!
   TBranch        *b_NrecoJetGen;   //!
   TBranch        *b_NrecoTowCal;   //!
   TBranch        *b_NrecoJetCal;   //!
   TBranch        *b_recoJetCalPt;   //!
   TBranch        *b_recoJetCalPhi;   //!
   TBranch        *b_recoJetCalEta;   //!
   TBranch        *b_recoJetCalE;   //!
   TBranch        *b_recoJetCalEMF;   //!
   TBranch        *b_recoJetCalN90;   //!
   TBranch        *b_recoJetCalN90hits;   //!
   TBranch        *b_NrecoJetCorCal;   //!
   TBranch        *b_recoJetCorCalPt;   //!
   TBranch        *b_recoJetCorCalPhi;   //!
   TBranch        *b_recoJetCorCalEta;   //!
   TBranch        *b_recoJetCorCalE;   //!
   TBranch        *b_recoJetCorCalEMF;   //!
   TBranch        *b_recoJetCorCalN90;   //!
   TBranch        *b_recoJetCorCalN90hits;   //!
   TBranch        *b_NohJetCal;   //!
   TBranch        *b_ohJetCalPt;   //!
   TBranch        *b_ohJetCalPhi;   //!
   TBranch        *b_ohJetCalEta;   //!
   TBranch        *b_ohJetCalE;   //!
   TBranch        *b_ohJetCalEMF;   //!
   TBranch        *b_ohJetCalN90;   //!
   TBranch        *b_ohJetCalN90hits;   //!
   TBranch        *b_NohJetCorCal;   //!
   TBranch        *b_ohJetCorCalPt;   //!
   TBranch        *b_ohJetCorCalPhi;   //!
   TBranch        *b_ohJetCorCalEta;   //!
   TBranch        *b_ohJetCorCalE;   //!
   TBranch        *b_ohJetCorCalEMF;   //!
   TBranch        *b_ohJetCorCalN90;   //!
   TBranch        *b_ohJetCorCalN90hits;   //!
   TBranch        *b_recoJetGenPt;   //!
   TBranch        *b_recoJetGenPhi;   //!
   TBranch        *b_recoJetGenEta;   //!
   TBranch        *b_recoJetGenE;   //!
   TBranch        *b_recoTowEt;   //!
   TBranch        *b_recoTowEta;   //!
   TBranch        *b_recoTowPhi;   //!
   TBranch        *b_recoTowE;   //!
   TBranch        *b_recoTowEm;   //!
   TBranch        *b_recoTowHad;   //!
   TBranch        *b_recoTowOE;   //!
   TBranch        *b_recoMetCal;   //!
   TBranch        *b_recoMetCalPhi;   //!
   TBranch        *b_recoMetCalSum;   //!
   TBranch        *b_recoMetGen;   //!
   TBranch        *b_recoMetGenPhi;   //!
   TBranch        *b_recoMetGenSum;   //!
   TBranch        *b_recoHTCal;   //!
   TBranch        *b_recoHTCalPhi;   //!
   TBranch        *b_recoHTCalSum;   //!
   TBranch        *b_recoMetPF;   //!
   TBranch        *b_recoMetPFSum;   //!
   TBranch        *b_recoMetPFPhi;   //!
   TBranch        *b_NohTau;   //!
   TBranch        *b_ohTauEta;   //!
   TBranch        *b_ohTauPhi;   //!
   TBranch        *b_ohTauPt;   //!
   TBranch        *b_ohTauEiso;   //!
   TBranch        *b_ohTauL25Tpt;   //!
   TBranch        *b_ohTauL3Tiso;   //!
   TBranch        *b_NohpfTau;   //!
   TBranch        *b_ohpfTauPt;   //!
   TBranch        *b_ohpfTauProngs;   //!
   TBranch        *b_ohpfTauEta;   //!
   TBranch        *b_ohpfTauPhi;   //!
   TBranch        *b_ohpfTauLeadTrackPt;   //!
   TBranch        *b_ohpfTauLeadPionPt;   //!
   TBranch        *b_ohpfTauTrkIso;   //!
   TBranch        *b_ohpfTauGammaIso;   //!
   TBranch        *b_ohpfTauJetPt;   //!
   TBranch        *b_NohpfTauTightCone;   //!
   TBranch        *b_ohpfTauTightConePt;   //!
   TBranch        *b_ohpfTauTightConeProngs;   //!
   TBranch        *b_ohpfTauTightConeEta;   //!
   TBranch        *b_ohpfTauTightConePhi;   //!
   TBranch        *b_ohpfTauTightConeLeadTrackPt;   //!
   TBranch        *b_ohpfTauTightConeLeadPionPt;   //!
   TBranch        *b_ohpfTauTightConeTrkIso;   //!
   TBranch        *b_ohpfTauTightConeGammaIso;   //!
   TBranch        *b_ohpfTauTightConeJetPt;   //!
   TBranch        *b_NRecoPFTau;   //!
   TBranch        *b_recopfTauPt;   //!
   TBranch        *b_recopfTauEta;   //!
   TBranch        *b_recopfTauPhi;   //!
   TBranch        *b_recopfTauLeadTrackPt;   //!
   TBranch        *b_recopfTauLeadPionPt;   //!
   TBranch        *b_recopfTauTrkIso;   //!
   TBranch        *b_recopfTauGammaIso;   //!
   TBranch        *b_recopfTauJetPt;   //!
   TBranch        *b_recopfTauDiscrByTancOnePercent;   //!
   TBranch        *b_recopfTauDiscrByTancHalfPercent;   //!
   TBranch        *b_recopfTauDiscrByTancQuarterPercent;   //!
   TBranch        *b_recopfTauDiscrByTancTenthPercent;   //!
   TBranch        *b_recopfTauDiscrByIso;   //!
   TBranch        *b_recopfTauDiscrAgainstMuon;   //!
   TBranch        *b_recopfTauDiscrAgainstElec;   //!
   TBranch        *b_pfHT;   //!
   TBranch        *b_pfMHT;   //!
   TBranch        *b_NohPFJet;   //!
   TBranch        *b_pfJetPt;   //!
   TBranch        *b_pfJetE;   //!
   TBranch        *b_pfJetEta;   //!
   TBranch        *b_pfJetPhi;   //!
   TBranch        *b_nrpj;   //!
   TBranch        *b_recopfJetpt;   //!
   TBranch        *b_recopfJete;   //!
   TBranch        *b_recopfJetphi;   //!
   TBranch        *b_recopfJeteta;   //!
   TBranch        *b_recopfJetneutralHadronFraction;   //!
   TBranch        *b_recopfJetneutralEMFraction;   //!
   TBranch        *b_recopfJetchargedHadronFraction;   //!
   TBranch        *b_recopfJetchargedEMFraction;   //!
   TBranch        *b_recopfJetneutralMultiplicity;   //!
   TBranch        *b_recopfJetchargedMultiplicity;   //!
   TBranch        *b_NohBJetL2;   //!
   TBranch        *b_ohBJetL2Energy;   //!
   TBranch        *b_ohBJetL2Et;   //!
   TBranch        *b_ohBJetL2Pt;   //!
   TBranch        *b_ohBJetL2Eta;   //!
   TBranch        *b_ohBJetL2Phi;   //!
   TBranch        *b_NohBJetL2Corrected;   //!
   TBranch        *b_ohBJetL2CorrectedEnergy;   //!
   TBranch        *b_ohBJetL2CorrectedEt;   //!
   TBranch        *b_ohBJetL2CorrectedPt;   //!
   TBranch        *b_ohBJetL2CorrectedEta;   //!
   TBranch        *b_ohBJetL2CorrectedPhi;   //!
   TBranch        *b_ohBJetIPL25Tag;   //!
   TBranch        *b_ohBJetIPL3Tag;   //!
   TBranch        *b_ohBJetIPL25TagSingleTrack;   //!
   TBranch        *b_ohBJetIPL3TagSingleTrack;   //!
   TBranch        *b_ohBJetPerfL25Tag;   //!
   TBranch        *b_ohBJetPerfL3Tag;   //!
   TBranch        *b_NrecoElec;   //!
   TBranch        *b_recoElecPt;   //!
   TBranch        *b_recoElecPhi;   //!
   TBranch        *b_recoElecEta;   //!
   TBranch        *b_recoElecEt;   //!
   TBranch        *b_recoElecE;   //!
   TBranch        *b_recoElecEleID;   //!
   TBranch        *b_recoElecIP;   //!
   TBranch        *b_recoElecNLostHits;   //!
   TBranch        *b_recoElecChi2NDF;   //!
   TBranch        *b_recoElecTrkIsoR03;   //!
   TBranch        *b_recoElecECaloIsoR03;   //!
   TBranch        *b_recoElecHCaloIsoR03;   //!
   TBranch        *b_recoElecIsEcalDriven;   //!
   TBranch        *b_recoElecFbrem;   //!
   TBranch        *b_recoElecmishits;   //!
   TBranch        *b_recoElecdist;   //!
   TBranch        *b_recoElecdcot;   //!
   TBranch        *b_recoElectrkiso;   //!
   TBranch        *b_recoElececaliso;   //!
   TBranch        *b_recoElechcaliso;   //!
   TBranch        *b_recoElecsigmaietaieta;   //!
   TBranch        *b_recoElecdeltaPhiIn;   //!
   TBranch        *b_recoElecdeltaEtaIn;   //!
   TBranch        *b_recoElechOverE;   //!
   TBranch        *b_recoElecscEt;   //!
   TBranch        *b_recoElecd0corr;   //!
   TBranch        *b_recoElecqGsfCtfScPixConsistent;   //!
   TBranch        *b_NrecoPhot;   //!
   TBranch        *b_recoPhotPt;   //!
   TBranch        *b_recoPhotPhi;   //!
   TBranch        *b_recoPhotEta;   //!
   TBranch        *b_recoPhotEt;   //!
   TBranch        *b_recoPhotE;   //!
   TBranch        *b_recoPhotTiso;   //!
   TBranch        *b_recoPhotEiso;   //!
   TBranch        *b_recoPhotHiso;   //!
   TBranch        *b_recoPhotHoverE;   //!
   TBranch        *b_recoPhotClusShap;   //!
   TBranch        *b_recoPhotR9ID;   //!
   TBranch        *b_NohPhot;   //!
   TBranch        *b_ohPhotEt;   //!
   TBranch        *b_ohPhotEta;   //!
   TBranch        *b_ohPhotPhi;   //!
   TBranch        *b_ohPhotEiso;   //!
   TBranch        *b_ohPhotHiso;   //!
   TBranch        *b_ohPhotTiso;   //!
   TBranch        *b_ohPhotL1iso;   //!
   TBranch        *b_ohPhotClusShap;   //!
   TBranch        *b_ohPhotR9;   //!
   TBranch        *b_ohPhotHforHoverE;   //!
   TBranch        *b_ohPhotR9ID;   //!
   TBranch        *b_NohEcalActiv;   //!
   TBranch        *b_ohEcalActivEt;   //!
   TBranch        *b_ohEcalActivEta;   //!
   TBranch        *b_ohEcalActivPhi;   //!
   TBranch        *b_ohEcalActivEiso;   //!
   TBranch        *b_ohEcalActivHiso;   //!
   TBranch        *b_ohEcalActivTiso;   //!
   TBranch        *b_ohEcalActivL1iso;   //!
   TBranch        *b_ohEcalActivClusShap;   //!
   TBranch        *b_ohEcalActivR9;   //!
   TBranch        *b_ohEcalActivHforHoverE;   //!
   TBranch        *b_ohEcalActivR9ID;   //!
   TBranch        *b_NohEle;   //!
   TBranch        *b_ohEleEt;   //!
   TBranch        *b_ohEleEta;   //!
   TBranch        *b_ohElePhi;   //!
   TBranch        *b_ohEleVtxZ;   //!
   TBranch        *b_ohEleE;   //!
   TBranch        *b_ohEleP;   //!
   TBranch        *b_ohEleHiso;   //!
   TBranch        *b_ohEleTiso;   //!
   TBranch        *b_ohEleEiso;   //!
   TBranch        *b_ohEleL1iso;   //!
   TBranch        *b_ohElePixelSeeds;   //!
   TBranch        *b_ohEleNewSC;   //!
   TBranch        *b_ohEleClusShap;   //!
   TBranch        *b_ohEleDeta;   //!
   TBranch        *b_ohEleDphi;   //!
   TBranch        *b_ohEleR9;   //!
   TBranch        *b_ohEleHforHoverE;   //!
   TBranch        *b_ohEleR9ID;   //!
   TBranch        *b_NohHFEle;   //!
   TBranch        *b_ohHFElePt;   //!
   TBranch        *b_ohHFEleEta;   //!
   TBranch        *b_NohHFECALClus;   //!
   TBranch        *b_ohHFEleClustere9e25;   //!
   TBranch        *b_ohHFEleClustere1e9;   //!
   TBranch        *b_ohHFEleClustereCOREe9;   //!
   TBranch        *b_ohHFEleClustereSeL;   //!
   TBranch        *b_ohHFEleCluster2Dcut;   //!
   TBranch        *b_ohHFEleClusterEta;   //!
   TBranch        *b_ohHFEleClusterPhi;   //!
   TBranch        *b_NrecoMuon;   //!
   TBranch        *b_recoMuonPt;   //!
   TBranch        *b_recoMuonPhi;   //!
   TBranch        *b_recoMuonEta;   //!
   TBranch        *b_recoMuonEt;   //!
   TBranch        *b_recoMuonE;   //!
   TBranch        *b_recoMuonChi2NDF;   //!
   TBranch        *b_recoMuonCharge;   //!
   TBranch        *b_recoMuonTrkIsoR03;   //!
   TBranch        *b_recoMuonECalIsoR03;   //!
   TBranch        *b_recoMuonHCalIsoR03;   //!
   TBranch        *b_recoMuonD0;   //!
   TBranch        *b_recoMuonType;   //!
   TBranch        *b_recoMuonNValidTrkHits;   //!
   TBranch        *b_recoMuonNValidMuonHits;   //!
   TBranch        *b_NohMuL2;   //!
   TBranch        *b_ohMuL2Pt;   //!
   TBranch        *b_ohMuL2Phi;   //!
   TBranch        *b_ohMuL2Eta;   //!
   TBranch        *b_ohMuL2Chg;   //!
   TBranch        *b_ohMuL2PtErr;   //!
   TBranch        *b_ohMuL2Iso;   //!
   TBranch        *b_ohMuL2Dr;   //!
   TBranch        *b_ohMuL2Dz;   //!
   TBranch        *b_ohMuL2VtxZ;   //!
   TBranch        *b_ohMuL2Nhits;   //!
   TBranch        *b_ohMuL2Nchambers;   //!
   TBranch        *b_ohMuL2Nstat;   //!
   TBranch        *b_ohMuL2L1idx;   //!
   TBranch        *b_NohMuL3;   //!
   TBranch        *b_ohMuL3Pt;   //!
   TBranch        *b_ohMuL3Phi;   //!
   TBranch        *b_ohMuL3Eta;   //!
   TBranch        *b_ohMuL3Chg;   //!
   TBranch        *b_ohMuL3PtErr;   //!
   TBranch        *b_ohMuL3Iso;   //!
   TBranch        *b_ohMuL3Trk10Iso;   //!
   TBranch        *b_ohMuL3Dr;   //!
   TBranch        *b_ohMuL3Dz;   //!
   TBranch        *b_ohMuL3VtxZ;   //!
   TBranch        *b_ohMuL3Nhits;   //!
   TBranch        *b_ohMuL3NormChi2;   //!
   TBranch        *b_ohMuL3Ntrackerhits;   //!
   TBranch        *b_ohMuL3Nmuonhits;   //!
   TBranch        *b_ohMuL3L2idx;   //!
   TBranch        *b_NohOniaPixel;   //!
   TBranch        *b_ohOniaPixelPt;   //!
   TBranch        *b_ohOniaPixelPhi;   //!
   TBranch        *b_ohOniaPixelEta;   //!
   TBranch        *b_ohOniaPixelChg;   //!
   TBranch        *b_ohOniaPixelDr;   //!
   TBranch        *b_ohOniaPixelDz;   //!
   TBranch        *b_ohOniaPixelHits;   //!
   TBranch        *b_ohOniaPixelNormChi2;   //!
   TBranch        *b_NohOniaTrack;   //!
   TBranch        *b_ohOniaTrackPt;   //!
   TBranch        *b_ohOniaTrackPhi;   //!
   TBranch        *b_ohOniaTrackEta;   //!
   TBranch        *b_ohOniaTrackChg;   //!
   TBranch        *b_ohOniaTrackDr;   //!
   TBranch        *b_ohOniaTrackDz;   //!
   TBranch        *b_ohOniaTrackHits;   //!
   TBranch        *b_ohOniaTrackNormChi2;   //!
   TBranch        *b_NohMuL2NoVtx;   //!
   TBranch        *b_ohMuL2NoVtxPt;   //!
   TBranch        *b_ohMuL2NoVtxPhi;   //!
   TBranch        *b_ohMuL2NoVtxEta;   //!
   TBranch        *b_ohMuL2NoVtxChg;   //!
   TBranch        *b_ohMuL2NoVtxPtErr;   //!
   TBranch        *b_ohMuL2NoVtxDr;   //!
   TBranch        *b_ohMuL2NoVtxDz;   //!
   TBranch        *b_ohMuL2NoVtxNhits;   //!
   TBranch        *b_ohMuL2NoVtxNchambers;   //!
   TBranch        *b_ohMuL2NoVtxL1idx;   //!
   TBranch        *b_NohDiMu;   //!
   TBranch        *b_ohDiMuDCA;   //!
   TBranch        *b_ohDiMu1st;   //!
   TBranch        *b_ohDiMu2nd;   //!
   TBranch        *b_NohDiMuVtx;   //!
   TBranch        *b_ohDiMuVtx1st;   //!
   TBranch        *b_ohDiMuVtx2nd;   //!
   TBranch        *b_ohDiMuVtxChi2;   //!
   TBranch        *b_ohDiMuVtxR;   //!
   TBranch        *b_ohDiMuVtxRSig;   //!
   TBranch        *b_ohDiMuVtxROverSig;   //!
   TBranch        *b_ohDiMuVtxCosAlpha;   //!
   TBranch        *b_ohDiMuVtxMu2DIpMax;   //!
   TBranch        *b_ohDiMuVtxMu2DIpMin;   //!
   TBranch        *b_ohDiMuVtxMu2DIpSigMax;   //!
   TBranch        *b_ohDiMuVtxMu2DIpSigMin;   //!
   TBranch        *b_NohIsoPixelTrackL3;   //!
   TBranch        *b_ohIsoPixelTrackL3Pt;   //!
   TBranch        *b_ohIsoPixelTrackL3Eta;   //!
   TBranch        *b_ohIsoPixelTrackL3Phi;   //!
   TBranch        *b_ohIsoPixelTrackL3MaxPtPxl;   //!
   TBranch        *b_ohIsoPixelTrackL3Energy;   //!
   TBranch        *b_ohIsoPixelTrackL2pt;   //!
   TBranch        *b_ohIsoPixelTrackL2eta;   //!
   TBranch        *b_ohIsoPixelTrackL2dXY;   //!
   TBranch        *b_NohPixelTracksL3;   //!
   TBranch        *b_ohPixelTracksL3Pt;   //!
   TBranch        *b_ohPixelTracksL3Eta;   //!
   TBranch        *b_ohPixelTracksL3Phi;   //!
   TBranch        *b_ohPixelTracksL3Vz;   //!
   TBranch        *b_ohPixelFEDSize;   //!
   TBranch        *b_NohPixelClusters;   //!
   TBranch        *b_NMCpart;   //!
   TBranch        *b_MCpid;   //!
   TBranch        *b_MCstatus;   //!
   TBranch        *b_MCvtxX;   //!
   TBranch        *b_MCvtxY;   //!
   TBranch        *b_MCvtxZ;   //!
   TBranch        *b_MCpt;   //!
   TBranch        *b_MCeta;   //!
   TBranch        *b_MCphi;   //!
   TBranch        *b_MCPtHat;   //!
   TBranch        *b_MCmu3;   //!
   TBranch        *b_MCel3;   //!
   TBranch        *b_MCbb;   //!
   TBranch        *b_MCab;   //!
   TBranch        *b_MCWenu;   //!
   TBranch        *b_MCmunu;   //!
   TBranch        *b_MCZee;   //!
   TBranch        *b_MCZmumu;   //!
   TBranch        *b_MCptEleMax;   //!
   TBranch        *b_MCptMuMax;   //!
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_L1MuChg;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1EtTot;   //!
   TBranch        *b_L1Mht;   //!
   TBranch        *b_L1MhtPhi;   //!
   TBranch        *b_L1EtHad;   //!
   TBranch        *b_L1HfRing1EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing2EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing1EtSumNegativeEta;   //!
   TBranch        *b_L1HfRing2EtSumNegativeEta;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing1;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing1;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing2;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing2;   //!
   TBranch        *b_NVrtx;   //!
   TBranch        *b_recoVrtXOffline0;   //!
   TBranch        *b_recoVrtYOffline0;   //!
   TBranch        *b_recoVrtZOffline0;   //!
   TBranch        *b_recoVrtNtrkOffline0;   //!
   TBranch        *b_recoVrtChi2Offline0;   //!
   TBranch        *b_recoVrtNdofOffline0;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_HLT_Activity_Ecal_SC7_v1;   //!
   TBranch        *b_HLT_Activity_Ecal_SC7_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleJet36_v1;   //!
   TBranch        *b_HLT_L1SingleJet36_v1_Prescl;   //!
   TBranch        *b_HLT_Jet30_v1;   //!
   TBranch        *b_HLT_Jet30_v1_Prescl;   //!
   TBranch        *b_HLT_Jet60_v1;   //!
   TBranch        *b_HLT_Jet60_v1_Prescl;   //!
   TBranch        *b_HLT_Jet80_v1;   //!
   TBranch        *b_HLT_Jet80_v1_Prescl;   //!
   TBranch        *b_HLT_Jet110_v1;   //!
   TBranch        *b_HLT_Jet110_v1_Prescl;   //!
   TBranch        *b_HLT_Jet150_v1;   //!
   TBranch        *b_HLT_Jet150_v1_Prescl;   //!
   TBranch        *b_HLT_Jet190_v1;   //!
   TBranch        *b_HLT_Jet190_v1_Prescl;   //!
   TBranch        *b_HLT_Jet240_v1;   //!
   TBranch        *b_HLT_Jet240_v1_Prescl;   //!
   TBranch        *b_HLT_Jet370_v1;   //!
   TBranch        *b_HLT_Jet370_v1_Prescl;   //!
   TBranch        *b_HLT_Jet370_NoJetID_v1;   //!
   TBranch        *b_HLT_Jet370_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_DiJetAve15U_v4;   //!
   TBranch        *b_HLT_DiJetAve15U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve30U_v4;   //!
   TBranch        *b_HLT_DiJetAve30U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve50U_v4;   //!
   TBranch        *b_HLT_DiJetAve50U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve70U_v4;   //!
   TBranch        *b_HLT_DiJetAve70U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve100U_v4;   //!
   TBranch        *b_HLT_DiJetAve100U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve140U_v4;   //!
   TBranch        *b_HLT_DiJetAve140U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve180U_v4;   //!
   TBranch        *b_HLT_DiJetAve180U_v4_Prescl;   //!
   TBranch        *b_HLT_DiJetAve300U_v4;   //!
   TBranch        *b_HLT_DiJetAve300U_v4_Prescl;   //!
   TBranch        *b_HLT_DoubleJet30_ForwardBackward_v2;   //!
   TBranch        *b_HLT_DoubleJet30_ForwardBackward_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleJet60_ForwardBackward_v2;   //!
   TBranch        *b_HLT_DoubleJet60_ForwardBackward_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleJet70_ForwardBackward_v2;   //!
   TBranch        *b_HLT_DoubleJet70_ForwardBackward_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleJet80_ForwardBackward_v2;   //!
   TBranch        *b_HLT_DoubleJet80_ForwardBackward_v2_Prescl;   //!
   TBranch        *b_HLT_CentralJet80_MET65_v1;   //!
   TBranch        *b_HLT_CentralJet80_MET65_v1_Prescl;   //!
   TBranch        *b_HLT_CentralJet80_MET80_v1;   //!
   TBranch        *b_HLT_CentralJet80_MET80_v1_Prescl;   //!
   TBranch        *b_HLT_CentralJet80_MET100_v1;   //!
   TBranch        *b_HLT_CentralJet80_MET100_v1_Prescl;   //!
   TBranch        *b_HLT_CentralJet80_MET160_v1;   //!
   TBranch        *b_HLT_CentralJet80_MET160_v1_Prescl;   //!
   TBranch        *b_HLT_DiJet60_MET45_v1;   //!
   TBranch        *b_HLT_DiJet60_MET45_v1_Prescl;   //!
   TBranch        *b_HLT_DiJet70_PT70_v1;   //!
   TBranch        *b_HLT_DiJet70_PT70_v1_Prescl;   //!
   TBranch        *b_HLT_DiJet100_PT100_v1;   //!
   TBranch        *b_HLT_DiJet100_PT100_v1_Prescl;   //!
   TBranch        *b_HLT_DiJet130_PT130_v1;   //!
   TBranch        *b_HLT_DiJet130_PT130_v1_Prescl;   //!
   TBranch        *b_HLT_QuadJet40_v2;   //!
   TBranch        *b_HLT_QuadJet40_v2_Prescl;   //!
   TBranch        *b_HLT_QuadJet40_IsoPFTau40_v1;   //!
   TBranch        *b_HLT_QuadJet40_IsoPFTau40_v1_Prescl;   //!
   TBranch        *b_HLT_QuadJet50_BTagIP_v1;   //!
   TBranch        *b_HLT_QuadJet50_BTagIP_v1_Prescl;   //!
   TBranch        *b_HLT_QuadJet50_Jet40_v1;   //!
   TBranch        *b_HLT_QuadJet50_Jet40_v1_Prescl;   //!
   TBranch        *b_HLT_QuadJet60_v1;   //!
   TBranch        *b_HLT_QuadJet60_v1_Prescl;   //!
   TBranch        *b_HLT_QuadJet70_v1;   //!
   TBranch        *b_HLT_QuadJet70_v1_Prescl;   //!
   TBranch        *b_HLT_ExclDiJet60_HFOR_v1;   //!
   TBranch        *b_HLT_ExclDiJet60_HFOR_v1_Prescl;   //!
   TBranch        *b_HLT_ExclDiJet60_HFAND_v1;   //!
   TBranch        *b_HLT_ExclDiJet60_HFAND_v1_Prescl;   //!
   TBranch        *b_HLT_JetE30_NoBPTX_v2;   //!
   TBranch        *b_HLT_JetE30_NoBPTX_v2_Prescl;   //!
   TBranch        *b_HLT_JetE30_NoBPTX_NoHalo_v4;   //!
   TBranch        *b_HLT_JetE30_NoBPTX_NoHalo_v4_Prescl;   //!
   TBranch        *b_HLT_JetE30_NoBPTX3BX_NoHalo_v4;   //!
   TBranch        *b_HLT_JetE30_NoBPTX3BX_NoHalo_v4_Prescl;   //!
   TBranch        *b_HLT_HT150_v2;   //!
   TBranch        *b_HLT_HT150_v2_Prescl;   //!
   TBranch        *b_HLT_HT150_AlphaT0p60_v1;   //!
   TBranch        *b_HLT_HT150_AlphaT0p60_v1_Prescl;   //!
   TBranch        *b_HLT_HT150_AlphaT0p70_v1;   //!
   TBranch        *b_HLT_HT150_AlphaT0p70_v1_Prescl;   //!
   TBranch        *b_HLT_HT200_v2;   //!
   TBranch        *b_HLT_HT200_v2_Prescl;   //!
   TBranch        *b_HLT_HT200_AlphaT0p60_v1;   //!
   TBranch        *b_HLT_HT200_AlphaT0p60_v1_Prescl;   //!
   TBranch        *b_HLT_HT200_AlphaT0p65_v1;   //!
   TBranch        *b_HLT_HT200_AlphaT0p65_v1_Prescl;   //!
   TBranch        *b_HLT_HT250_v2;   //!
   TBranch        *b_HLT_HT250_v2_Prescl;   //!
   TBranch        *b_HLT_HT250_AlphaT0p55_v1;   //!
   TBranch        *b_HLT_HT250_AlphaT0p55_v1_Prescl;   //!
   TBranch        *b_HLT_HT250_AlphaT0p62_v1;   //!
   TBranch        *b_HLT_HT250_AlphaT0p62_v1_Prescl;   //!
   TBranch        *b_HLT_HT250_DoubleDisplacedJet60_v1;   //!
   TBranch        *b_HLT_HT250_DoubleDisplacedJet60_v1_Prescl;   //!
   TBranch        *b_HLT_HT250_MHT60_v2;   //!
   TBranch        *b_HLT_HT250_MHT60_v2_Prescl;   //!
   TBranch        *b_HLT_HT300_v3;   //!
   TBranch        *b_HLT_HT300_v3_Prescl;   //!
   TBranch        *b_HLT_HT300_MHT75_v3;   //!
   TBranch        *b_HLT_HT300_MHT75_v3_Prescl;   //!
   TBranch        *b_HLT_HT300_AlphaT0p52_v1;   //!
   TBranch        *b_HLT_HT300_AlphaT0p52_v1_Prescl;   //!
   TBranch        *b_HLT_HT300_AlphaT0p54_v1;   //!
   TBranch        *b_HLT_HT300_AlphaT0p54_v1_Prescl;   //!
   TBranch        *b_HLT_HT350_v2;   //!
   TBranch        *b_HLT_HT350_v2_Prescl;   //!
   TBranch        *b_HLT_HT350_AlphaT0p51_v1;   //!
   TBranch        *b_HLT_HT350_AlphaT0p51_v1_Prescl;   //!
   TBranch        *b_HLT_HT350_AlphaT0p53_v1;   //!
   TBranch        *b_HLT_HT350_AlphaT0p53_v1_Prescl;   //!
   TBranch        *b_HLT_HT400_v2;   //!
   TBranch        *b_HLT_HT400_v2_Prescl;   //!
   TBranch        *b_HLT_HT400_AlphaT0p51_v1;   //!
   TBranch        *b_HLT_HT400_AlphaT0p51_v1_Prescl;   //!
   TBranch        *b_HLT_HT450_v2;   //!
   TBranch        *b_HLT_HT450_v2_Prescl;   //!
   TBranch        *b_HLT_HT500_v2;   //!
   TBranch        *b_HLT_HT500_v2_Prescl;   //!
   TBranch        *b_HLT_HT550_v2;   //!
   TBranch        *b_HLT_HT550_v2_Prescl;   //!
   TBranch        *b_HLT_PFMHT150_v2;   //!
   TBranch        *b_HLT_PFMHT150_v2_Prescl;   //!
   TBranch        *b_HLT_MET100_v1;   //!
   TBranch        *b_HLT_MET100_v1_Prescl;   //!
   TBranch        *b_HLT_MET120_v1;   //!
   TBranch        *b_HLT_MET120_v1_Prescl;   //!
   TBranch        *b_HLT_MET200_v1;   //!
   TBranch        *b_HLT_MET200_v1_Prescl;   //!
   TBranch        *b_HLT_Meff440_v2;   //!
   TBranch        *b_HLT_Meff440_v2_Prescl;   //!
   TBranch        *b_HLT_Meff520_v2;   //!
   TBranch        *b_HLT_Meff520_v2_Prescl;   //!
   TBranch        *b_HLT_Meff640_v2;   //!
   TBranch        *b_HLT_Meff640_v2_Prescl;   //!
   TBranch        *b_HLT_MR100_v1;   //!
   TBranch        *b_HLT_MR100_v1_Prescl;   //!
   TBranch        *b_HLT_R032_v1;   //!
   TBranch        *b_HLT_R032_v1_Prescl;   //!
   TBranch        *b_HLT_R032_MR100_v1;   //!
   TBranch        *b_HLT_R032_MR100_v1_Prescl;   //!
   TBranch        *b_HLT_R035_MR100_v1;   //!
   TBranch        *b_HLT_R035_MR100_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleMuOpen_v1;   //!
   TBranch        *b_HLT_L1SingleMuOpen_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleMuOpen_DT_v1;   //!
   TBranch        *b_HLT_L1SingleMuOpen_DT_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleMu10_v1;   //!
   TBranch        *b_HLT_L1SingleMu10_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleMu20_v1;   //!
   TBranch        *b_HLT_L1SingleMu20_v1_Prescl;   //!
   TBranch        *b_HLT_L1DoubleMu0_v1;   //!
   TBranch        *b_HLT_L1DoubleMu0_v1_Prescl;   //!
   TBranch        *b_HLT_L2Mu10_v1;   //!
   TBranch        *b_HLT_L2Mu10_v1_Prescl;   //!
   TBranch        *b_HLT_L2Mu20_v1;   //!
   TBranch        *b_HLT_L2Mu20_v1_Prescl;   //!
   TBranch        *b_HLT_L2DoubleMu0_v2;   //!
   TBranch        *b_HLT_L2DoubleMu0_v2_Prescl;   //!
   TBranch        *b_HLT_Mu3_v3;   //!
   TBranch        *b_HLT_Mu3_v3_Prescl;   //!
   TBranch        *b_HLT_Mu5_v3;   //!
   TBranch        *b_HLT_Mu5_v3_Prescl;   //!
   TBranch        *b_HLT_Mu8_v1;   //!
   TBranch        *b_HLT_Mu8_v1_Prescl;   //!
   TBranch        *b_HLT_Mu12_v1;   //!
   TBranch        *b_HLT_Mu12_v1_Prescl;   //!
   TBranch        *b_HLT_Mu15_v2;   //!
   TBranch        *b_HLT_Mu15_v2_Prescl;   //!
   TBranch        *b_HLT_Mu20_v1;   //!
   TBranch        *b_HLT_Mu20_v1_Prescl;   //!
   TBranch        *b_HLT_Mu24_v1;   //!
   TBranch        *b_HLT_Mu24_v1_Prescl;   //!
   TBranch        *b_HLT_Mu30_v1;   //!
   TBranch        *b_HLT_Mu30_v1_Prescl;   //!
   TBranch        *b_HLT_IsoMu12_v1;   //!
   TBranch        *b_HLT_IsoMu12_v1_Prescl;   //!
   TBranch        *b_HLT_IsoMu15_v5;   //!
   TBranch        *b_HLT_IsoMu15_v5_Prescl;   //!
   TBranch        *b_HLT_IsoMu17_v5;   //!
   TBranch        *b_HLT_IsoMu17_v5_Prescl;   //!
   TBranch        *b_HLT_IsoMu24_v1;   //!
   TBranch        *b_HLT_IsoMu24_v1_Prescl;   //!
   TBranch        *b_HLT_IsoMu30_v1;   //!
   TBranch        *b_HLT_IsoMu30_v1_Prescl;   //!
   TBranch        *b_HLT_L2DoubleMu23_NoVertex_v1;   //!
   TBranch        *b_HLT_L2DoubleMu23_NoVertex_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_v3;   //!
   TBranch        *b_HLT_DoubleMu3_v3_Prescl;   //!
   TBranch        *b_HLT_DoubleMu6_v1;   //!
   TBranch        *b_HLT_DoubleMu6_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu7_v1;   //!
   TBranch        *b_HLT_DoubleMu7_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v1;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_Jpsi_v2;   //!
   TBranch        *b_HLT_DoubleMu3_Jpsi_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_Quarkonium_v2;   //!
   TBranch        *b_HLT_DoubleMu3_Quarkonium_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_Upsilon_v1;   //!
   TBranch        *b_HLT_DoubleMu3_Upsilon_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_LowMass_v1;   //!
   TBranch        *b_HLT_DoubleMu3_LowMass_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu4_Acoplanarity03_v1;   //!
   TBranch        *b_HLT_DoubleMu4_Acoplanarity03_v1_Prescl;   //!
   TBranch        *b_HLT_TripleMu5_v2;   //!
   TBranch        *b_HLT_TripleMu5_v2_Prescl;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_v2;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_v2_Prescl;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v2_Prescl;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v5;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v5_Prescl;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v1_Prescl;   //!
   TBranch        *b_HLT_Mu7_Track5_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu7_Track5_Jpsi_v2_Prescl;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v2_Prescl;   //!
   TBranch        *b_HLT_Photon20_CaloIdVL_IsoL_v1;   //!
   TBranch        *b_HLT_Photon20_CaloIdVL_IsoL_v1_Prescl;   //!
   TBranch        *b_HLT_Photon20_R9Id_Photon18_R9Id_v2;   //!
   TBranch        *b_HLT_Photon20_R9Id_Photon18_R9Id_v2_Prescl;   //!
   TBranch        *b_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2;   //!
   TBranch        *b_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon26_Photon18_v2;   //!
   TBranch        *b_HLT_Photon26_Photon18_v2_Prescl;   //!
   TBranch        *b_HLT_Photon26_IsoVL_Photon18_v2;   //!
   TBranch        *b_HLT_Photon26_IsoVL_Photon18_v2_Prescl;   //!
   TBranch        *b_HLT_Photon26_IsoVL_Photon18_IsoVL_v2;   //!
   TBranch        *b_HLT_Photon26_IsoVL_Photon18_IsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon26_CaloIdL_IsoVL_Photon18_v2;   //!
   TBranch        *b_HLT_Photon26_CaloIdL_IsoVL_Photon18_v2_Prescl;   //!
   TBranch        *b_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1;   //!
   TBranch        *b_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1_Prescl;   //!
   TBranch        *b_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2;   //!
   TBranch        *b_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1;   //!
   TBranch        *b_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1_Prescl;   //!
   TBranch        *b_HLT_Photon30_CaloIdVL_v2;   //!
   TBranch        *b_HLT_Photon30_CaloIdVL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon30_CaloIdVL_IsoL_v2;   //!
   TBranch        *b_HLT_Photon30_CaloIdVL_IsoL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2;   //!
   TBranch        *b_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1;   //!
   TBranch        *b_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1_Prescl;   //!
   TBranch        *b_HLT_Photon50_CaloIdVL_IsoL_v1;   //!
   TBranch        *b_HLT_Photon50_CaloIdVL_IsoL_v1_Prescl;   //!
   TBranch        *b_HLT_Photon60_CaloIdL_HT200_v2;   //!
   TBranch        *b_HLT_Photon60_CaloIdL_HT200_v2_Prescl;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_HT200_v2;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_HT200_v2_Prescl;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_HT300_v2;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_HT300_v2_Prescl;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_MHT30_v2;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_MHT30_v2_Prescl;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_MHT50_v2;   //!
   TBranch        *b_HLT_Photon70_CaloIdL_MHT50_v2_Prescl;   //!
   TBranch        *b_HLT_Photon75_CaloIdVL_v2;   //!
   TBranch        *b_HLT_Photon75_CaloIdVL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon75_CaloIdVL_IsoL_v2;   //!
   TBranch        *b_HLT_Photon75_CaloIdVL_IsoL_v2_Prescl;   //!
   TBranch        *b_HLT_Photon125_NoSpikeFilter_v2;   //!
   TBranch        *b_HLT_Photon125_NoSpikeFilter_v2_Prescl;   //!
   TBranch        *b_HLT_DoublePhoton33_v2;   //!
   TBranch        *b_HLT_DoublePhoton33_v2_Prescl;   //!
   TBranch        *b_HLT_DoublePhoton5_IsoVL_CEP_v1;   //!
   TBranch        *b_HLT_DoublePhoton5_IsoVL_CEP_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleEG5_v1;   //!
   TBranch        *b_HLT_L1SingleEG5_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleEG12_v1;   //!
   TBranch        *b_HLT_L1SingleEG12_v1_Prescl;   //!
   TBranch        *b_HLT_Ele8_v2;   //!
   TBranch        *b_HLT_Ele8_v2_Prescl;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_CaloIsoVL_v2;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_TrkIdVL_v2;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_TrkIdVL_v2_Prescl;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_CaloIsoVL_v2;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_CaloIsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2;   //!
   TBranch        *b_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2_Prescl;   //!
   TBranch        *b_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2;   //!
   TBranch        *b_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2_Prescl;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2_Prescl;   //!
   TBranch        *b_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2;   //!
   TBranch        *b_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl;   //!
   TBranch        *b_HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1;   //!
   TBranch        *b_HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1_Prescl;   //!
   TBranch        *b_HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2;   //!
   TBranch        *b_HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2_Prescl;   //!
   TBranch        *b_HLT_Ele45_CaloIdVT_TrkIdT_v2;   //!
   TBranch        *b_HLT_Ele45_CaloIdVT_TrkIdT_v2_Prescl;   //!
   TBranch        *b_HLT_Ele90_NoSpikeFilter_v2;   //!
   TBranch        *b_HLT_Ele90_NoSpikeFilter_v2_Prescl;   //!
   TBranch        *b_HLT_IsoPFTau35_Trk20_MET45_v2;   //!
   TBranch        *b_HLT_IsoPFTau35_Trk20_MET45_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleIsoPFTau20_Trk5_v2;   //!
   TBranch        *b_HLT_DoubleIsoPFTau20_Trk5_v2_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet20_Mu5_v2;   //!
   TBranch        *b_HLT_BTagMu_DiJet20_Mu5_v2_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet60_Mu7_v2;   //!
   TBranch        *b_HLT_BTagMu_DiJet60_Mu7_v2_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet80_Mu9_v2;   //!
   TBranch        *b_HLT_BTagMu_DiJet80_Mu9_v2_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet100_Mu9_v2;   //!
   TBranch        *b_HLT_BTagMu_DiJet100_Mu9_v2_Prescl;   //!
   TBranch        *b_HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3;   //!
   TBranch        *b_HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3_Prescl;   //!
   TBranch        *b_HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3;   //!
   TBranch        *b_HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3_Prescl;   //!
   TBranch        *b_HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3;   //!
   TBranch        *b_HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3_Prescl;   //!
   TBranch        *b_HLT_Mu5_DoubleEle8_v3;   //!
   TBranch        *b_HLT_Mu5_DoubleEle8_v3_Prescl;   //!
   TBranch        *b_HLT_Mu5_HT200_v4;   //!
   TBranch        *b_HLT_Mu5_HT200_v4_Prescl;   //!
   TBranch        *b_HLT_Mu8_HT200_v3;   //!
   TBranch        *b_HLT_Mu8_HT200_v3_Prescl;   //!
   TBranch        *b_HLT_Mu8_Ele17_CaloIdL_v2;   //!
   TBranch        *b_HLT_Mu8_Ele17_CaloIdL_v2_Prescl;   //!
   TBranch        *b_HLT_Mu8_Photon20_CaloIdVT_IsoT_v2;   //!
   TBranch        *b_HLT_Mu8_Photon20_CaloIdVT_IsoT_v2_Prescl;   //!
   TBranch        *b_HLT_Mu8_Jet40_v3;   //!
   TBranch        *b_HLT_Mu8_Jet40_v3_Prescl;   //!
   TBranch        *b_HLT_Mu10_Ele10_CaloIdL_v3;   //!
   TBranch        *b_HLT_Mu10_Ele10_CaloIdL_v3_Prescl;   //!
   TBranch        *b_HLT_Mu15_Photon20_CaloIdL_v3;   //!
   TBranch        *b_HLT_Mu15_Photon20_CaloIdL_v3_Prescl;   //!
   TBranch        *b_HLT_Mu15_DoublePhoton15_CaloIdL_v3;   //!
   TBranch        *b_HLT_Mu15_DoublePhoton15_CaloIdL_v3_Prescl;   //!
   TBranch        *b_HLT_Mu15_LooseIsoPFTau20_v2;   //!
   TBranch        *b_HLT_Mu15_LooseIsoPFTau20_v2_Prescl;   //!
   TBranch        *b_HLT_Mu17_CentralJet30_v2;   //!
   TBranch        *b_HLT_Mu17_CentralJet30_v2_Prescl;   //!
   TBranch        *b_HLT_Mu17_DiCentralJet30_v2;   //!
   TBranch        *b_HLT_Mu17_DiCentralJet30_v2_Prescl;   //!
   TBranch        *b_HLT_Mu17_TriCentralJet30_v2;   //!
   TBranch        *b_HLT_Mu17_TriCentralJet30_v2_Prescl;   //!
   TBranch        *b_HLT_Mu17_Ele8_CaloIdL_v2;   //!
   TBranch        *b_HLT_Mu17_Ele8_CaloIdL_v2_Prescl;   //!
   TBranch        *b_HLT_Mu17_CentralJet40_BTagIP_v2;   //!
   TBranch        *b_HLT_Mu17_CentralJet40_BTagIP_v2_Prescl;   //!
   TBranch        *b_HLT_IsoMu12_LooseIsoPFTau10_v2;   //!
   TBranch        *b_HLT_IsoMu12_LooseIsoPFTau10_v2_Prescl;   //!
   TBranch        *b_HLT_IsoMu17_CentralJet40_BTagIP_v2;   //!
   TBranch        *b_HLT_IsoMu17_CentralJet40_BTagIP_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_HT160_v3;   //!
   TBranch        *b_HLT_DoubleMu3_HT160_v3_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_HT200_v3;   //!
   TBranch        *b_HLT_DoubleMu3_HT200_v3_Prescl;   //!
   TBranch        *b_HLT_DoubleMu5_Ele8_v3;   //!
   TBranch        *b_HLT_DoubleMu5_Ele8_v3_Prescl;   //!
   TBranch        *b_HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3;   //!
   TBranch        *b_HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3_Prescl;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2_Prescl;   //!
   TBranch        *b_HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3;   //!
   TBranch        *b_HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3_Prescl;   //!
   TBranch        *b_HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3;   //!
   TBranch        *b_HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3_Prescl;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2_Prescl;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2_Prescl;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2;   //!
   TBranch        *b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2_Prescl;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2_Prescl;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2_Prescl;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2_Prescl;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2;   //!
   TBranch        *b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3_Prescl;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3_Prescl;   //!
   TBranch        *b_HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2;   //!
   TBranch        *b_HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2_Prescl;   //!
   TBranch        *b_HLT_TripleEle10_CaloIdL_TrkIdVL_v2;   //!
   TBranch        *b_HLT_TripleEle10_CaloIdL_TrkIdVL_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity80_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity80_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity100_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity100_v2_Prescl;   //!
   TBranch        *b_HLT_BeamGas_HF_v2;   //!
   TBranch        *b_HLT_BeamGas_HF_v2_Prescl;   //!
   TBranch        *b_HLT_BeamGas_BSC_v2;   //!
   TBranch        *b_HLT_BeamGas_BSC_v2_Prescl;   //!
   TBranch        *b_HLT_BeamHalo_v2;   //!
   TBranch        *b_HLT_BeamHalo_v2_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_minBias_threshold1_v1;   //!
   TBranch        *b_HLT_L1Tech_BSC_minBias_threshold1_v1_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_halo_v1;   //!
   TBranch        *b_HLT_L1Tech_BSC_halo_v1_Prescl;   //!
   TBranch        *b_HLT_L1Tech_CASTOR_HaloMuon_v1;   //!
   TBranch        *b_HLT_L1Tech_CASTOR_HaloMuon_v1_Prescl;   //!
   TBranch        *b_HLT_L1_PreCollisions_v1;   //!
   TBranch        *b_HLT_L1_PreCollisions_v1_Prescl;   //!
   TBranch        *b_HLT_L1_Interbunch_BSC_v1;   //!
   TBranch        *b_HLT_L1_Interbunch_BSC_v1_Prescl;   //!
   TBranch        *b_HLT_IsoTrackHE_v3;   //!
   TBranch        *b_HLT_IsoTrackHE_v3_Prescl;   //!
   TBranch        *b_HLT_IsoTrackHB_v2;   //!
   TBranch        *b_HLT_IsoTrackHB_v2_Prescl;   //!
   TBranch        *b_HLT_HcalPhiSym_v3;   //!
   TBranch        *b_HLT_HcalPhiSym_v3_Prescl;   //!
   TBranch        *b_HLT_HcalNZS_v3;   //!
   TBranch        *b_HLT_HcalNZS_v3_Prescl;   //!
   TBranch        *b_HLT_GlobalRunHPDNoise_v2;   //!
   TBranch        *b_HLT_GlobalRunHPDNoise_v2_Prescl;   //!
   TBranch        *b_HLT_L1Tech_HBHEHO_totalOR_v1;   //!
   TBranch        *b_HLT_L1Tech_HBHEHO_totalOR_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_v1;   //!
   TBranch        *b_HLT_ZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_Physics_v1;   //!
   TBranch        *b_HLT_Physics_v1_Prescl;   //!
   TBranch        *b_HLT_Physics_NanoDST_v1;   //!
   TBranch        *b_HLT_Physics_NanoDST_v1_Prescl;   //!
   TBranch        *b_HLT_Calibration_v1;   //!
   TBranch        *b_HLT_Calibration_v1_Prescl;   //!
   TBranch        *b_HLT_EcalCalibration_v1;   //!
   TBranch        *b_HLT_EcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HcalCalibration_v1;   //!
   TBranch        *b_HLT_HcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_TrackerCalibration_v1;   //!
   TBranch        *b_HLT_TrackerCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_Random_v1;   //!
   TBranch        *b_HLT_Random_v1_Prescl;   //!
   TBranch        *b_HLT_L1SingleMuOpen_AntiBPTX_v1;   //!
   TBranch        *b_HLT_L1SingleMuOpen_AntiBPTX_v1_Prescl;   //!
   TBranch        *b_HLT_L1TrackerCosmics_v2;   //!
   TBranch        *b_HLT_L1TrackerCosmics_v2_Prescl;   //!
   TBranch        *b_HLT_RegionalCosmicTracking_v1;   //!
   TBranch        *b_HLT_RegionalCosmicTracking_v1_Prescl;   //!
   TBranch        *b_HLT_L3MuonsCosmicTracking_v1;   //!
   TBranch        *b_HLT_L3MuonsCosmicTracking_v1_Prescl;   //!
   TBranch        *b_HLT_LogMonitor_v1;   //!
   TBranch        *b_HLT_LogMonitor_v1_Prescl;   //!
   TBranch        *b_HLT_DTErrors_v1;   //!
   TBranch        *b_HLT_DTErrors_v1_Prescl;   //!
   TBranch        *b_AlCa_EcalPi0_v4;   //!
   TBranch        *b_AlCa_EcalPi0_v4_Prescl;   //!
   TBranch        *b_AlCa_EcalEta_v3;   //!
   TBranch        *b_AlCa_EcalEta_v3_Prescl;   //!
   TBranch        *b_AlCa_EcalPhiSym_v2;   //!
   TBranch        *b_AlCa_EcalPhiSym_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v2;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v2;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v2;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v2_Prescl;   //!
   TBranch        *b_DQM_FEDIntegrity_v3;   //!
   TBranch        *b_DQM_FEDIntegrity_v3_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_L1_BeamGas_Bsc;   //!
   TBranch        *b_L1_BeamGas_Bsc_Prescl;   //!
   TBranch        *b_L1_BeamGas_Bsc_5bx;   //!
   TBranch        *b_L1_BeamGas_Hf;   //!
   TBranch        *b_L1_BeamGas_Hf_Prescl;   //!
   TBranch        *b_L1_BeamGas_Hf_5bx;   //!
   TBranch        *b_L1_BeamHalo;   //!
   TBranch        *b_L1_BeamHalo_Prescl;   //!
   TBranch        *b_L1_BeamHalo_5bx;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_5bx;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_5bx;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus_Prescl;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus_5bx;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus_Prescl;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus_5bx;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusANDMinus;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusANDMinus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusANDMinus_5bx;   //!
   TBranch        *b_L1_DoubleEG10;   //!
   TBranch        *b_L1_DoubleEG10_Prescl;   //!
   TBranch        *b_L1_DoubleEG10_5bx;   //!
   TBranch        *b_L1_DoubleEG2_FwdVeto;   //!
   TBranch        *b_L1_DoubleEG2_FwdVeto_Prescl;   //!
   TBranch        *b_L1_DoubleEG2_FwdVeto_5bx;   //!
   TBranch        *b_L1_DoubleEG3;   //!
   TBranch        *b_L1_DoubleEG3_Prescl;   //!
   TBranch        *b_L1_DoubleEG3_5bx;   //!
   TBranch        *b_L1_DoubleEG5;   //!
   TBranch        *b_L1_DoubleEG5_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_5bx;   //!
   TBranch        *b_L1_DoubleEG5_HTT50;   //!
   TBranch        *b_L1_DoubleEG5_HTT50_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_HTT50_5bx;   //!
   TBranch        *b_L1_DoubleEG5_HTT75;   //!
   TBranch        *b_L1_DoubleEG5_HTT75_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_HTT75_5bx;   //!
   TBranch        *b_L1_DoubleEG8;   //!
   TBranch        *b_L1_DoubleEG8_Prescl;   //!
   TBranch        *b_L1_DoubleEG8_5bx;   //!
   TBranch        *b_L1_DoubleEG_12_5;   //!
   TBranch        *b_L1_DoubleEG_12_5_Prescl;   //!
   TBranch        *b_L1_DoubleEG_12_5_5bx;   //!
   TBranch        *b_L1_DoubleEG_12_5_Eta1p39;   //!
   TBranch        *b_L1_DoubleEG_12_5_Eta1p39_Prescl;   //!
   TBranch        *b_L1_DoubleEG_12_5_Eta1p39_5bx;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp_5bx;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp_5bx;   //!
   TBranch        *b_L1_DoubleIsoEG10;   //!
   TBranch        *b_L1_DoubleIsoEG10_Prescl;   //!
   TBranch        *b_L1_DoubleIsoEG10_5bx;   //!
   TBranch        *b_L1_DoubleJet36_Central;   //!
   TBranch        *b_L1_DoubleJet36_Central_Prescl;   //!
   TBranch        *b_L1_DoubleJet36_Central_5bx;   //!
   TBranch        *b_L1_DoubleJet44_Central;   //!
   TBranch        *b_L1_DoubleJet44_Central_Prescl;   //!
   TBranch        *b_L1_DoubleJet44_Central_5bx;   //!
   TBranch        *b_L1_DoubleJet52;   //!
   TBranch        *b_L1_DoubleJet52_Prescl;   //!
   TBranch        *b_L1_DoubleJet52_5bx;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_5bx;   //!
   TBranch        *b_L1_DoubleMu0_HighQ;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_5bx;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_EtaCuts;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_EtaCuts_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_HighQ_EtaCuts_5bx;   //!
   TBranch        *b_L1_DoubleMu3;   //!
   TBranch        *b_L1_DoubleMu3_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_5bx;   //!
   TBranch        *b_L1_DoubleMu3_EG5;   //!
   TBranch        *b_L1_DoubleMu3_EG5_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_EG5_5bx;   //!
   TBranch        *b_L1_DoubleMu3p5;   //!
   TBranch        *b_L1_DoubleMu3p5_Prescl;   //!
   TBranch        *b_L1_DoubleMu3p5_5bx;   //!
   TBranch        *b_L1_DoubleMu5_v1;   //!
   TBranch        *b_L1_DoubleMu5_v1_Prescl;   //!
   TBranch        *b_L1_DoubleMu5_v1_5bx;   //!
   TBranch        *b_L1_DoubleTauJet28;   //!
   TBranch        *b_L1_DoubleTauJet28_Prescl;   //!
   TBranch        *b_L1_DoubleTauJet28_5bx;   //!
   TBranch        *b_L1_DoubleTauJet32;   //!
   TBranch        *b_L1_DoubleTauJet32_Prescl;   //!
   TBranch        *b_L1_DoubleTauJet32_5bx;   //!
   TBranch        *b_L1_DoubleTauJet36;   //!
   TBranch        *b_L1_DoubleTauJet36_Prescl;   //!
   TBranch        *b_L1_DoubleTauJet36_5bx;   //!
   TBranch        *b_L1_DoubleTauJet40;   //!
   TBranch        *b_L1_DoubleTauJet40_Prescl;   //!
   TBranch        *b_L1_DoubleTauJet40_5bx;   //!
   TBranch        *b_L1_EG10_Jet24_Central_deltaPhi1;   //!
   TBranch        *b_L1_EG10_Jet24_Central_deltaPhi1_Prescl;   //!
   TBranch        *b_L1_EG10_Jet24_Central_deltaPhi1_5bx;   //!
   TBranch        *b_L1_EG12_Jet24_Central_deltaPhi1;   //!
   TBranch        *b_L1_EG12_Jet24_Central_deltaPhi1_Prescl;   //!
   TBranch        *b_L1_EG12_Jet24_Central_deltaPhi1_5bx;   //!
   TBranch        *b_L1_EG12_TauJet20_deltaPhi1;   //!
   TBranch        *b_L1_EG12_TauJet20_deltaPhi1_Prescl;   //!
   TBranch        *b_L1_EG12_TauJet20_deltaPhi1_5bx;   //!
   TBranch        *b_L1_EG5_HTT100;   //!
   TBranch        *b_L1_EG5_HTT100_Prescl;   //!
   TBranch        *b_L1_EG5_HTT100_5bx;   //!
   TBranch        *b_L1_EG5_HTT125;   //!
   TBranch        *b_L1_EG5_HTT125_Prescl;   //!
   TBranch        *b_L1_EG5_HTT125_5bx;   //!
   TBranch        *b_L1_EG5_HTT75;   //!
   TBranch        *b_L1_EG5_HTT75_Prescl;   //!
   TBranch        *b_L1_EG5_HTT75_5bx;   //!
   TBranch        *b_L1_EG5_Jet36_deltaPhi1;   //!
   TBranch        *b_L1_EG5_Jet36_deltaPhi1_Prescl;   //!
   TBranch        *b_L1_EG5_Jet36_deltaPhi1_5bx;   //!
   TBranch        *b_L1_EG8_Jet20_Central_deltaPhi1;   //!
   TBranch        *b_L1_EG8_Jet20_Central_deltaPhi1_Prescl;   //!
   TBranch        *b_L1_EG8_Jet20_Central_deltaPhi1_5bx;   //!
   TBranch        *b_L1_ETM100;   //!
   TBranch        *b_L1_ETM100_Prescl;   //!
   TBranch        *b_L1_ETM100_5bx;   //!
   TBranch        *b_L1_ETM20;   //!
   TBranch        *b_L1_ETM20_Prescl;   //!
   TBranch        *b_L1_ETM20_5bx;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM30_Prescl;   //!
   TBranch        *b_L1_ETM30_5bx;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM50_Prescl;   //!
   TBranch        *b_L1_ETM50_5bx;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM70_Prescl;   //!
   TBranch        *b_L1_ETM70_5bx;   //!
   TBranch        *b_L1_ETT220;   //!
   TBranch        *b_L1_ETT220_Prescl;   //!
   TBranch        *b_L1_ETT220_5bx;   //!
   TBranch        *b_L1_ETT260_EG5;   //!
   TBranch        *b_L1_ETT260_EG5_Prescl;   //!
   TBranch        *b_L1_ETT260_EG5_5bx;   //!
   TBranch        *b_L1_ETT300_EG5;   //!
   TBranch        *b_L1_ETT300_EG5_Prescl;   //!
   TBranch        *b_L1_ETT300_EG5_5bx;   //!
   TBranch        *b_L1_HTM50;   //!
   TBranch        *b_L1_HTM50_Prescl;   //!
   TBranch        *b_L1_HTM50_5bx;   //!
   TBranch        *b_L1_HTT100;   //!
   TBranch        *b_L1_HTT100_Prescl;   //!
   TBranch        *b_L1_HTT100_5bx;   //!
   TBranch        *b_L1_HTT150;   //!
   TBranch        *b_L1_HTT150_Prescl;   //!
   TBranch        *b_L1_HTT150_5bx;   //!
   TBranch        *b_L1_HTT50;   //!
   TBranch        *b_L1_HTT50_Prescl;   //!
   TBranch        *b_L1_HTT50_5bx;   //!
   TBranch        *b_L1_HTT50_HTM30;   //!
   TBranch        *b_L1_HTT50_HTM30_Prescl;   //!
   TBranch        *b_L1_HTT50_HTM30_5bx;   //!
   TBranch        *b_L1_HTT50_HTM50;   //!
   TBranch        *b_L1_HTT50_HTM50_Prescl;   //!
   TBranch        *b_L1_HTT50_HTM50_5bx;   //!
   TBranch        *b_L1_HTT75;   //!
   TBranch        *b_L1_HTT75_Prescl;   //!
   TBranch        *b_L1_HTT75_5bx;   //!
   TBranch        *b_L1_InterBunch_Bsc;   //!
   TBranch        *b_L1_InterBunch_Bsc_Prescl;   //!
   TBranch        *b_L1_InterBunch_Bsc_5bx;   //!
   TBranch        *b_L1_InterBunch_Hf;   //!
   TBranch        *b_L1_InterBunch_Hf_Prescl;   //!
   TBranch        *b_L1_InterBunch_Hf_5bx;   //!
   TBranch        *b_L1_Mu0_HTT50;   //!
   TBranch        *b_L1_Mu0_HTT50_Prescl;   //!
   TBranch        *b_L1_Mu0_HTT50_5bx;   //!
   TBranch        *b_L1_Mu0_HTT75;   //!
   TBranch        *b_L1_Mu0_HTT75_Prescl;   //!
   TBranch        *b_L1_Mu0_HTT75_5bx;   //!
   TBranch        *b_L1_Mu10_Jet36_Central;   //!
   TBranch        *b_L1_Mu10_Jet36_Central_Prescl;   //!
   TBranch        *b_L1_Mu10_Jet36_Central_5bx;   //!
   TBranch        *b_L1_Mu12_EG5;   //!
   TBranch        *b_L1_Mu12_EG5_Prescl;   //!
   TBranch        *b_L1_Mu12_EG5_5bx;   //!
   TBranch        *b_L1_Mu3_DoubleEG5;   //!
   TBranch        *b_L1_Mu3_DoubleEG5_Prescl;   //!
   TBranch        *b_L1_Mu3_DoubleEG5_5bx;   //!
   TBranch        *b_L1_Mu3_EG5;   //!
   TBranch        *b_L1_Mu3_EG5_Prescl;   //!
   TBranch        *b_L1_Mu3_EG5_5bx;   //!
   TBranch        *b_L1_Mu3_Jet16_Central;   //!
   TBranch        *b_L1_Mu3_Jet16_Central_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet16_Central_5bx;   //!
   TBranch        *b_L1_Mu3_Jet20_Central;   //!
   TBranch        *b_L1_Mu3_Jet20_Central_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet20_Central_5bx;   //!
   TBranch        *b_L1_Mu3_Jet28_Central;   //!
   TBranch        *b_L1_Mu3_Jet28_Central_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet28_Central_5bx;   //!
   TBranch        *b_L1_Mu5_EG12;   //!
   TBranch        *b_L1_Mu5_EG12_Prescl;   //!
   TBranch        *b_L1_Mu5_EG12_5bx;   //!
   TBranch        *b_L1_Mu7_EG5;   //!
   TBranch        *b_L1_Mu7_EG5_Prescl;   //!
   TBranch        *b_L1_Mu7_EG5_5bx;   //!
   TBranch        *b_L1_Mu7_Jet20_Central;   //!
   TBranch        *b_L1_Mu7_Jet20_Central_Prescl;   //!
   TBranch        *b_L1_Mu7_Jet20_Central_5bx;   //!
   TBranch        *b_L1_Mu7_TauJet16;   //!
   TBranch        *b_L1_Mu7_TauJet16_Prescl;   //!
   TBranch        *b_L1_Mu7_TauJet16_5bx;   //!
   TBranch        *b_L1_MuOpen_EG12;   //!
   TBranch        *b_L1_MuOpen_EG12_Prescl;   //!
   TBranch        *b_L1_MuOpen_EG12_5bx;   //!
   TBranch        *b_L1_MuOpen_EG5;   //!
   TBranch        *b_L1_MuOpen_EG5_Prescl;   //!
   TBranch        *b_L1_MuOpen_EG5_5bx;   //!
   TBranch        *b_L1_PreCollisions;   //!
   TBranch        *b_L1_PreCollisions_Prescl;   //!
   TBranch        *b_L1_PreCollisions_5bx;   //!
   TBranch        *b_L1_QuadJet20_Central;   //!
   TBranch        *b_L1_QuadJet20_Central_Prescl;   //!
   TBranch        *b_L1_QuadJet20_Central_5bx;   //!
   TBranch        *b_L1_QuadJet28_Central;   //!
   TBranch        *b_L1_QuadJet28_Central_Prescl;   //!
   TBranch        *b_L1_QuadJet28_Central_5bx;   //!
   TBranch        *b_L1_SingleEG12;   //!
   TBranch        *b_L1_SingleEG12_Prescl;   //!
   TBranch        *b_L1_SingleEG12_5bx;   //!
   TBranch        *b_L1_SingleEG12_Eta1p39;   //!
   TBranch        *b_L1_SingleEG12_Eta1p39_Prescl;   //!
   TBranch        *b_L1_SingleEG12_Eta1p39_5bx;   //!
   TBranch        *b_L1_SingleEG12_Eta2p17;   //!
   TBranch        *b_L1_SingleEG12_Eta2p17_Prescl;   //!
   TBranch        *b_L1_SingleEG12_Eta2p17_5bx;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG15_Prescl;   //!
   TBranch        *b_L1_SingleEG15_5bx;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG20_5bx;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG30_Prescl;   //!
   TBranch        *b_L1_SingleEG30_5bx;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG5_Prescl;   //!
   TBranch        *b_L1_SingleEG5_5bx;   //!
   TBranch        *b_L1_SingleIsoEG12;   //!
   TBranch        *b_L1_SingleIsoEG12_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG12_5bx;   //!
   TBranch        *b_L1_SingleIsoEG12_Eta1p39;   //!
   TBranch        *b_L1_SingleIsoEG12_Eta1p39_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG12_Eta1p39_5bx;   //!
   TBranch        *b_L1_SingleIsoEG12_Eta2p17;   //!
   TBranch        *b_L1_SingleIsoEG12_Eta2p17_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG12_Eta2p17_5bx;   //!
   TBranch        *b_L1_SingleJet128;   //!
   TBranch        *b_L1_SingleJet128_Prescl;   //!
   TBranch        *b_L1_SingleJet128_5bx;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet16_Prescl;   //!
   TBranch        *b_L1_SingleJet16_5bx;   //!
   TBranch        *b_L1_SingleJet20_NotBptxOR;   //!
   TBranch        *b_L1_SingleJet20_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJet20_NotBptxOR_5bx;   //!
   TBranch        *b_L1_SingleJet20_NotBptxOR_NotMuBeamHalo;   //!
   TBranch        *b_L1_SingleJet20_NotBptxOR_NotMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleJet20_NotBptxOR_NotMuBeamHalo_5bx;   //!
   TBranch        *b_L1_SingleJet32_NotBptxOR_NotMuBeamHalo;   //!
   TBranch        *b_L1_SingleJet32_NotBptxOR_NotMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleJet32_NotBptxOR_NotMuBeamHalo_5bx;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet36_5bx;   //!
   TBranch        *b_L1_SingleJet36_FwdVeto;   //!
   TBranch        *b_L1_SingleJet36_FwdVeto_Prescl;   //!
   TBranch        *b_L1_SingleJet36_FwdVeto_5bx;   //!
   TBranch        *b_L1_SingleJet52;   //!
   TBranch        *b_L1_SingleJet52_Prescl;   //!
   TBranch        *b_L1_SingleJet52_5bx;   //!
   TBranch        *b_L1_SingleJet68;   //!
   TBranch        *b_L1_SingleJet68_Prescl;   //!
   TBranch        *b_L1_SingleJet68_5bx;   //!
   TBranch        *b_L1_SingleJet80_Central;   //!
   TBranch        *b_L1_SingleJet80_Central_Prescl;   //!
   TBranch        *b_L1_SingleJet80_Central_5bx;   //!
   TBranch        *b_L1_SingleJet92;   //!
   TBranch        *b_L1_SingleJet92_Prescl;   //!
   TBranch        *b_L1_SingleJet92_5bx;   //!
   TBranch        *b_L1_SingleJet92_Central;   //!
   TBranch        *b_L1_SingleJet92_Central_Prescl;   //!
   TBranch        *b_L1_SingleJet92_Central_5bx;   //!
   TBranch        *b_L1_SingleMu10;   //!
   TBranch        *b_L1_SingleMu10_Prescl;   //!
   TBranch        *b_L1_SingleMu10_5bx;   //!
   TBranch        *b_L1_SingleMu12;   //!
   TBranch        *b_L1_SingleMu12_Prescl;   //!
   TBranch        *b_L1_SingleMu12_5bx;   //!
   TBranch        *b_L1_SingleMu12_Debug;   //!
   TBranch        *b_L1_SingleMu12_Debug_Prescl;   //!
   TBranch        *b_L1_SingleMu12_Debug_5bx;   //!
   TBranch        *b_L1_SingleMu16;   //!
   TBranch        *b_L1_SingleMu16_Prescl;   //!
   TBranch        *b_L1_SingleMu16_5bx;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20_Prescl;   //!
   TBranch        *b_L1_SingleMu20_5bx;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleMu25_Prescl;   //!
   TBranch        *b_L1_SingleMu25_5bx;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu3_Prescl;   //!
   TBranch        *b_L1_SingleMu3_5bx;   //!
   TBranch        *b_L1_SingleMu5_Eta1p5_Q80_v1;   //!
   TBranch        *b_L1_SingleMu5_Eta1p5_Q80_v1_Prescl;   //!
   TBranch        *b_L1_SingleMu5_Eta1p5_Q80_v1_5bx;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_Prescl;   //!
   TBranch        *b_L1_SingleMu7_5bx;   //!
   TBranch        *b_L1_SingleMu7_Barrel;   //!
   TBranch        *b_L1_SingleMu7_Barrel_Prescl;   //!
   TBranch        *b_L1_SingleMu7_Barrel_5bx;   //!
   TBranch        *b_L1_SingleMu7_Eta2p1;   //!
   TBranch        *b_L1_SingleMu7_Eta2p1_Prescl;   //!
   TBranch        *b_L1_SingleMu7_Eta2p1_5bx;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo_5bx;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_5bx;   //!
   TBranch        *b_L1_SingleTauJet52;   //!
   TBranch        *b_L1_SingleTauJet52_Prescl;   //!
   TBranch        *b_L1_SingleTauJet52_5bx;   //!
   TBranch        *b_L1_SingleTauJet68;   //!
   TBranch        *b_L1_SingleTauJet68_Prescl;   //!
   TBranch        *b_L1_SingleTauJet68_5bx;   //!
   TBranch        *b_L1_SingleTauJet80;   //!
   TBranch        *b_L1_SingleTauJet80_Prescl;   //!
   TBranch        *b_L1_SingleTauJet80_5bx;   //!
   TBranch        *b_L1_TripleEG5;   //!
   TBranch        *b_L1_TripleEG5_Prescl;   //!
   TBranch        *b_L1_TripleEG5_5bx;   //!
   TBranch        *b_L1_TripleEG7;   //!
   TBranch        *b_L1_TripleEG7_Prescl;   //!
   TBranch        *b_L1_TripleEG7_5bx;   //!
   TBranch        *b_L1_TripleEG_8_5_5;   //!
   TBranch        *b_L1_TripleEG_8_5_5_Prescl;   //!
   TBranch        *b_L1_TripleEG_8_5_5_5bx;   //!
   TBranch        *b_L1_TripleEG_8_8_5;   //!
   TBranch        *b_L1_TripleEG_8_8_5_Prescl;   //!
   TBranch        *b_L1_TripleEG_8_8_5_5bx;   //!
   TBranch        *b_L1_TripleJet28_Central;   //!
   TBranch        *b_L1_TripleJet28_Central_Prescl;   //!
   TBranch        *b_L1_TripleJet28_Central_5bx;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!
   TBranch        *b_L1_ZeroBias_5bx;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_5bx;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0_5bx;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0_5bx;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_5bx;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_5bx;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_5bx;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v0_5bx;   //!
   TBranch        *b_L1Tech_HCAL_HF_MM_or_PP_or_PM_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_5bx;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v1;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v1_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v1_5bx;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_Scint_loose_vertex_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_loose_vertex_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_loose_vertex_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_Scint_minus_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_minus_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_Scint_plus_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_plus_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_Scint_tight_vertex_v0;   //!
   TBranch        *b_L1Tech_ZDC_Scint_tight_vertex_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_Scint_tight_vertex_v0_5bx;   //!

   checkOpenHLT(TTree *tree=0);
   virtual ~checkOpenHLT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void prepareHistograms(const TString & name, const int bins, const double & min, const double & max, const TString & axisTitle);
   void prepareAllHistograms(const TString & name, TFile * outputFile);
   void fillNameArray(std::string * nameArray);
   void applyCuts(const int arraySize, const bool selectOnChambers, const double & parallelDiff, const bool selectOnParallelism, bool * selectionArray);

   template <class T>
   void fillHistograms(const TString & name, const T * valueArray, const int valueArraySize, const bool * selectionArray)
   {
     std::locale loc;
     std::string nameArray[4];
     fillNameArray(nameArray);

     for( int i=0; i<valueArraySize; ++i ) {
       std::string nameCopy(nameArray[i]);
       TString namePart(name + std::toupper(nameArray[i][0], loc) + nameCopy.erase(0,1));
       if( selectionArray[i] ) {
         histoMap_[namePart]->Fill(valueArray[i]);
       }
       for( int j=i+1; j<valueArraySize; ++j ) {
         std::stringstream ss;
         ss << i << "_" << j;
         TString correlationName(name+"Correlation_"+ss.str());
         if( selectionArray[i] && selectionArray[j] ) {
           histoMap_[correlationName]->Fill(valueArray[i], valueArray[j]);
         }
       }
     }
   }

   void fillAllHistograms(const TString & name, const int arraySize, const bool * selectionArray);
   void saveHistograms(const TString & name);
   void saveAllHistograms(const TString & name);
   void saveHistogram(TH1F * histo);
   TLorentzVector fromPtEtaPhiToPxPyPz( const double & pt, const double & eta, const double & phi );

   std::map<TString, TH1*> histoMap_;
   TString dir_;
   double parallelDiff_;
   bool defaultTriggerCuts_;
};

#endif

#ifdef checkOpenHLT_cxx
checkOpenHLT::checkOpenHLT(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("openhlt_merge.root");
      if (!f) {
         f = new TFile("openhlt_merge.root");
      }
      tree = (TTree*)gDirectory->Get("HltTree");

   }
   Init(tree);
   Loop();
}

checkOpenHLT::~checkOpenHLT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t checkOpenHLT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t checkOpenHLT::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void checkOpenHLT::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NoRecoPFTausSignal", &NoRecoPFTausSignal, &b_NoRecoPFTausSignal);
   fChain->SetBranchAddress("signalTrToPFTauMatch", &signalTrToPFTauMatch, &b_signalTrToPFTauMatch);
   fChain->SetBranchAddress("recoPFTauSignalTrDz", &recoPFTauSignalTrDz, &b_recoPFTauSignalTrDz);
   fChain->SetBranchAddress("recoPFTauSignalTrPt", &recoPFTauSignalTrPt, &b_recoPFTauSignalTrPt);
   fChain->SetBranchAddress("NoRecoPFTausIso", &NoRecoPFTausIso, &b_NoRecoPFTausIso);
   fChain->SetBranchAddress("isoTrToPFTauMatch", &isoTrToPFTauMatch, &b_isoTrToPFTauMatch);
   fChain->SetBranchAddress("recoPFTauIsoTrDz", &recoPFTauIsoTrDz, &b_recoPFTauIsoTrDz);
   fChain->SetBranchAddress("recoPFTauIsoTrPt", &recoPFTauIsoTrPt, &b_recoPFTauIsoTrPt);
   fChain->SetBranchAddress("NoHLTPFTausSignal", &NoHLTPFTausSignal, &b_NoHLTPFTausSignal);
   fChain->SetBranchAddress("hltpftauSignalTrToPFTauMatch", hltpftauSignalTrToPFTauMatch, &b_hltpftauSignalTrToPFTauMatch);
   fChain->SetBranchAddress("HLTPFTauSignalTrDz", HLTPFTauSignalTrDz, &b_HLTPFTauSignalTrDz);
   fChain->SetBranchAddress("HLTPFTauSignalTrPt", HLTPFTauSignalTrPt, &b_HLTPFTauSignalTrPt);
   fChain->SetBranchAddress("NoHLTPFTausIso", &NoHLTPFTausIso, &b_NoHLTPFTausIso);
   fChain->SetBranchAddress("hltpftauIsoTrToPFTauMatch", hltpftauIsoTrToPFTauMatch, &b_hltpftauIsoTrToPFTauMatch);
   fChain->SetBranchAddress("HLTPFTauIsoTrDz", HLTPFTauIsoTrDz, &b_HLTPFTauIsoTrDz);
   fChain->SetBranchAddress("HLTPFTauIsoTrPt", HLTPFTauIsoTrPt, &b_HLTPFTauIsoTrPt);
   fChain->SetBranchAddress("NrecoJetGen", &NrecoJetGen, &b_NrecoJetGen);
   fChain->SetBranchAddress("NrecoTowCal", &NrecoTowCal, &b_NrecoTowCal);
   fChain->SetBranchAddress("NrecoJetCal", &NrecoJetCal, &b_NrecoJetCal);
   fChain->SetBranchAddress("recoJetCalPt", &recoJetCalPt, &b_recoJetCalPt);
   fChain->SetBranchAddress("recoJetCalPhi", &recoJetCalPhi, &b_recoJetCalPhi);
   fChain->SetBranchAddress("recoJetCalEta", &recoJetCalEta, &b_recoJetCalEta);
   fChain->SetBranchAddress("recoJetCalE", &recoJetCalE, &b_recoJetCalE);
   fChain->SetBranchAddress("recoJetCalEMF", &recoJetCalEMF, &b_recoJetCalEMF);
   fChain->SetBranchAddress("recoJetCalN90", &recoJetCalN90, &b_recoJetCalN90);
   fChain->SetBranchAddress("recoJetCalN90hits", &recoJetCalN90hits, &b_recoJetCalN90hits);
   fChain->SetBranchAddress("NrecoJetCorCal", &NrecoJetCorCal, &b_NrecoJetCorCal);
   fChain->SetBranchAddress("recoJetCorCalPt", &recoJetCorCalPt, &b_recoJetCorCalPt);
   fChain->SetBranchAddress("recoJetCorCalPhi", &recoJetCorCalPhi, &b_recoJetCorCalPhi);
   fChain->SetBranchAddress("recoJetCorCalEta", &recoJetCorCalEta, &b_recoJetCorCalEta);
   fChain->SetBranchAddress("recoJetCorCalE", &recoJetCorCalE, &b_recoJetCorCalE);
   fChain->SetBranchAddress("recoJetCorCalEMF", &recoJetCorCalEMF, &b_recoJetCorCalEMF);
   fChain->SetBranchAddress("recoJetCorCalN90", &recoJetCorCalN90, &b_recoJetCorCalN90);
   fChain->SetBranchAddress("recoJetCorCalN90hits", &recoJetCorCalN90hits, &b_recoJetCorCalN90hits);
   fChain->SetBranchAddress("NohJetCal", &NohJetCal, &b_NohJetCal);
   fChain->SetBranchAddress("ohJetCalPt", ohJetCalPt, &b_ohJetCalPt);
   fChain->SetBranchAddress("ohJetCalPhi", ohJetCalPhi, &b_ohJetCalPhi);
   fChain->SetBranchAddress("ohJetCalEta", ohJetCalEta, &b_ohJetCalEta);
   fChain->SetBranchAddress("ohJetCalE", ohJetCalE, &b_ohJetCalE);
   fChain->SetBranchAddress("ohJetCalEMF", ohJetCalEMF, &b_ohJetCalEMF);
   fChain->SetBranchAddress("ohJetCalN90", ohJetCalN90, &b_ohJetCalN90);
   fChain->SetBranchAddress("ohJetCalN90hits", ohJetCalN90hits, &b_ohJetCalN90hits);
   fChain->SetBranchAddress("NohJetCorCal", &NohJetCorCal, &b_NohJetCorCal);
   fChain->SetBranchAddress("ohJetCorCalPt", ohJetCorCalPt, &b_ohJetCorCalPt);
   fChain->SetBranchAddress("ohJetCorCalPhi", ohJetCorCalPhi, &b_ohJetCorCalPhi);
   fChain->SetBranchAddress("ohJetCorCalEta", ohJetCorCalEta, &b_ohJetCorCalEta);
   fChain->SetBranchAddress("ohJetCorCalE", ohJetCorCalE, &b_ohJetCorCalE);
   fChain->SetBranchAddress("ohJetCorCalEMF", ohJetCorCalEMF, &b_ohJetCorCalEMF);
   fChain->SetBranchAddress("ohJetCorCalN90", ohJetCorCalN90, &b_ohJetCorCalN90);
   fChain->SetBranchAddress("ohJetCorCalN90hits", ohJetCorCalN90hits, &b_ohJetCorCalN90hits);
   fChain->SetBranchAddress("recoJetGenPt", &recoJetGenPt, &b_recoJetGenPt);
   fChain->SetBranchAddress("recoJetGenPhi", &recoJetGenPhi, &b_recoJetGenPhi);
   fChain->SetBranchAddress("recoJetGenEta", &recoJetGenEta, &b_recoJetGenEta);
   fChain->SetBranchAddress("recoJetGenE", &recoJetGenE, &b_recoJetGenE);
   fChain->SetBranchAddress("recoTowEt", recoTowEt, &b_recoTowEt);
   fChain->SetBranchAddress("recoTowEta", recoTowEta, &b_recoTowEta);
   fChain->SetBranchAddress("recoTowPhi", recoTowPhi, &b_recoTowPhi);
   fChain->SetBranchAddress("recoTowE", recoTowE, &b_recoTowE);
   fChain->SetBranchAddress("recoTowEm", recoTowEm, &b_recoTowEm);
   fChain->SetBranchAddress("recoTowHad", recoTowHad, &b_recoTowHad);
   fChain->SetBranchAddress("recoTowOE", recoTowOE, &b_recoTowOE);
   fChain->SetBranchAddress("recoMetCal", &recoMetCal, &b_recoMetCal);
   fChain->SetBranchAddress("recoMetCalPhi", &recoMetCalPhi, &b_recoMetCalPhi);
   fChain->SetBranchAddress("recoMetCalSum", &recoMetCalSum, &b_recoMetCalSum);
   fChain->SetBranchAddress("recoMetGen", &recoMetGen, &b_recoMetGen);
   fChain->SetBranchAddress("recoMetGenPhi", &recoMetGenPhi, &b_recoMetGenPhi);
   fChain->SetBranchAddress("recoMetGenSum", &recoMetGenSum, &b_recoMetGenSum);
   fChain->SetBranchAddress("recoHTCal", &recoHTCal, &b_recoHTCal);
   fChain->SetBranchAddress("recoHTCalPhi", &recoHTCalPhi, &b_recoHTCalPhi);
   fChain->SetBranchAddress("recoHTCalSum", &recoHTCalSum, &b_recoHTCalSum);
   fChain->SetBranchAddress("recoMetPF", &recoMetPF, &b_recoMetPF);
   fChain->SetBranchAddress("recoMetPFSum", &recoMetPFSum, &b_recoMetPFSum);
   fChain->SetBranchAddress("recoMetPFPhi", &recoMetPFPhi, &b_recoMetPFPhi);
   fChain->SetBranchAddress("NohTau", &NohTau, &b_NohTau);
   fChain->SetBranchAddress("ohTauEta", ohTauEta, &b_ohTauEta);
   fChain->SetBranchAddress("ohTauPhi", ohTauPhi, &b_ohTauPhi);
   fChain->SetBranchAddress("ohTauPt", ohTauPt, &b_ohTauPt);
   fChain->SetBranchAddress("ohTauEiso", ohTauEiso, &b_ohTauEiso);
   fChain->SetBranchAddress("ohTauL25Tpt", ohTauL25Tpt, &b_ohTauL25Tpt);
   fChain->SetBranchAddress("ohTauL3Tiso", ohTauL3Tiso, &b_ohTauL3Tiso);
   fChain->SetBranchAddress("NohpfTau", &NohpfTau, &b_NohpfTau);
   fChain->SetBranchAddress("ohpfTauPt", ohpfTauPt, &b_ohpfTauPt);
   fChain->SetBranchAddress("ohpfTauProngs", ohpfTauProngs, &b_ohpfTauProngs);
   fChain->SetBranchAddress("ohpfTauEta", ohpfTauEta, &b_ohpfTauEta);
   fChain->SetBranchAddress("ohpfTauPhi", ohpfTauPhi, &b_ohpfTauPhi);
   fChain->SetBranchAddress("ohpfTauLeadTrackPt", ohpfTauLeadTrackPt, &b_ohpfTauLeadTrackPt);
   fChain->SetBranchAddress("ohpfTauLeadPionPt", ohpfTauLeadPionPt, &b_ohpfTauLeadPionPt);
   fChain->SetBranchAddress("ohpfTauTrkIso", ohpfTauTrkIso, &b_ohpfTauTrkIso);
   fChain->SetBranchAddress("ohpfTauGammaIso", ohpfTauGammaIso, &b_ohpfTauGammaIso);
   fChain->SetBranchAddress("ohpfTauJetPt", ohpfTauJetPt, &b_ohpfTauJetPt);
   fChain->SetBranchAddress("NohpfTauTightCone", &NohpfTauTightCone, &b_NohpfTauTightCone);
   fChain->SetBranchAddress("ohpfTauTightConePt", ohpfTauTightConePt, &b_ohpfTauTightConePt);
   fChain->SetBranchAddress("ohpfTauTightConeProngs", ohpfTauTightConeProngs, &b_ohpfTauTightConeProngs);
   fChain->SetBranchAddress("ohpfTauTightConeEta", ohpfTauTightConeEta, &b_ohpfTauTightConeEta);
   fChain->SetBranchAddress("ohpfTauTightConePhi", ohpfTauTightConePhi, &b_ohpfTauTightConePhi);
   fChain->SetBranchAddress("ohpfTauTightConeLeadTrackPt", ohpfTauTightConeLeadTrackPt, &b_ohpfTauTightConeLeadTrackPt);
   fChain->SetBranchAddress("ohpfTauTightConeLeadPionPt", ohpfTauTightConeLeadPionPt, &b_ohpfTauTightConeLeadPionPt);
   fChain->SetBranchAddress("ohpfTauTightConeTrkIso", ohpfTauTightConeTrkIso, &b_ohpfTauTightConeTrkIso);
   fChain->SetBranchAddress("ohpfTauTightConeGammaIso", ohpfTauTightConeGammaIso, &b_ohpfTauTightConeGammaIso);
   fChain->SetBranchAddress("ohpfTauTightConeJetPt", ohpfTauTightConeJetPt, &b_ohpfTauTightConeJetPt);
   fChain->SetBranchAddress("NRecoPFTau", &NRecoPFTau, &b_NRecoPFTau);
   fChain->SetBranchAddress("recopfTauPt", &recopfTauPt, &b_recopfTauPt);
   fChain->SetBranchAddress("recopfTauEta", &recopfTauEta, &b_recopfTauEta);
   fChain->SetBranchAddress("recopfTauPhi", &recopfTauPhi, &b_recopfTauPhi);
   fChain->SetBranchAddress("recopfTauLeadTrackPt", &recopfTauLeadTrackPt, &b_recopfTauLeadTrackPt);
   fChain->SetBranchAddress("recopfTauLeadPionPt", &recopfTauLeadPionPt, &b_recopfTauLeadPionPt);
   fChain->SetBranchAddress("recopfTauTrkIso", &recopfTauTrkIso, &b_recopfTauTrkIso);
   fChain->SetBranchAddress("recopfTauGammaIso", &recopfTauGammaIso, &b_recopfTauGammaIso);
   fChain->SetBranchAddress("recopfTauJetPt", &recopfTauJetPt, &b_recopfTauJetPt);
   fChain->SetBranchAddress("recopfTauDiscrByTancOnePercent", &recopfTauDiscrByTancOnePercent, &b_recopfTauDiscrByTancOnePercent);
   fChain->SetBranchAddress("recopfTauDiscrByTancHalfPercent", &recopfTauDiscrByTancHalfPercent, &b_recopfTauDiscrByTancHalfPercent);
   fChain->SetBranchAddress("recopfTauDiscrByTancQuarterPercent", &recopfTauDiscrByTancQuarterPercent, &b_recopfTauDiscrByTancQuarterPercent);
   fChain->SetBranchAddress("recopfTauDiscrByTancTenthPercent", &recopfTauDiscrByTancTenthPercent, &b_recopfTauDiscrByTancTenthPercent);
   fChain->SetBranchAddress("recopfTauDiscrByIso", &recopfTauDiscrByIso, &b_recopfTauDiscrByIso);
   fChain->SetBranchAddress("recopfTauDiscrAgainstMuon", &recopfTauDiscrAgainstMuon, &b_recopfTauDiscrAgainstMuon);
   fChain->SetBranchAddress("recopfTauDiscrAgainstElec", &recopfTauDiscrAgainstElec, &b_recopfTauDiscrAgainstElec);
   fChain->SetBranchAddress("pfHT", &pfHT, &b_pfHT);
   fChain->SetBranchAddress("pfMHT", &pfMHT, &b_pfMHT);
   fChain->SetBranchAddress("NohPFJet", &NohPFJet, &b_NohPFJet);
   fChain->SetBranchAddress("pfJetPt", pfJetPt, &b_pfJetPt);
   fChain->SetBranchAddress("pfJetE", pfJetE, &b_pfJetE);
   fChain->SetBranchAddress("pfJetEta", pfJetEta, &b_pfJetEta);
   fChain->SetBranchAddress("pfJetPhi", pfJetPhi, &b_pfJetPhi);
   fChain->SetBranchAddress("nrpj", &nrpj, &b_nrpj);
   fChain->SetBranchAddress("recopfJetpt", &recopfJetpt, &b_recopfJetpt);
   fChain->SetBranchAddress("recopfJete", &recopfJete, &b_recopfJete);
   fChain->SetBranchAddress("recopfJetphi", &recopfJetphi, &b_recopfJetphi);
   fChain->SetBranchAddress("recopfJeteta", &recopfJeteta, &b_recopfJeteta);
   fChain->SetBranchAddress("recopfJetneutralHadronFraction", &recopfJetneutralHadronFraction, &b_recopfJetneutralHadronFraction);
   fChain->SetBranchAddress("recopfJetneutralEMFraction", &recopfJetneutralEMFraction, &b_recopfJetneutralEMFraction);
   fChain->SetBranchAddress("recopfJetchargedHadronFraction", &recopfJetchargedHadronFraction, &b_recopfJetchargedHadronFraction);
   fChain->SetBranchAddress("recopfJetchargedEMFraction", &recopfJetchargedEMFraction, &b_recopfJetchargedEMFraction);
   fChain->SetBranchAddress("recopfJetneutralMultiplicity", &recopfJetneutralMultiplicity, &b_recopfJetneutralMultiplicity);
   fChain->SetBranchAddress("recopfJetchargedMultiplicity", &recopfJetchargedMultiplicity, &b_recopfJetchargedMultiplicity);
   fChain->SetBranchAddress("NohBJetL2", &NohBJetL2, &b_NohBJetL2);
   fChain->SetBranchAddress("ohBJetL2Energy", ohBJetL2Energy, &b_ohBJetL2Energy);
   fChain->SetBranchAddress("ohBJetL2Et", ohBJetL2Et, &b_ohBJetL2Et);
   fChain->SetBranchAddress("ohBJetL2Pt", ohBJetL2Pt, &b_ohBJetL2Pt);
   fChain->SetBranchAddress("ohBJetL2Eta", ohBJetL2Eta, &b_ohBJetL2Eta);
   fChain->SetBranchAddress("ohBJetL2Phi", ohBJetL2Phi, &b_ohBJetL2Phi);
   fChain->SetBranchAddress("NohBJetL2Corrected", &NohBJetL2Corrected, &b_NohBJetL2Corrected);
   fChain->SetBranchAddress("ohBJetL2CorrectedEnergy", ohBJetL2CorrectedEnergy, &b_ohBJetL2CorrectedEnergy);
   fChain->SetBranchAddress("ohBJetL2CorrectedEt", ohBJetL2CorrectedEt, &b_ohBJetL2CorrectedEt);
   fChain->SetBranchAddress("ohBJetL2CorrectedPt", ohBJetL2CorrectedPt, &b_ohBJetL2CorrectedPt);
   fChain->SetBranchAddress("ohBJetL2CorrectedEta", ohBJetL2CorrectedEta, &b_ohBJetL2CorrectedEta);
   fChain->SetBranchAddress("ohBJetL2CorrectedPhi", ohBJetL2CorrectedPhi, &b_ohBJetL2CorrectedPhi);
   fChain->SetBranchAddress("ohBJetIPL25Tag", ohBJetIPL25Tag, &b_ohBJetIPL25Tag);
   fChain->SetBranchAddress("ohBJetIPL3Tag", ohBJetIPL3Tag, &b_ohBJetIPL3Tag);
   fChain->SetBranchAddress("ohBJetIPL25TagSingleTrack", ohBJetIPL25TagSingleTrack, &b_ohBJetIPL25TagSingleTrack);
   fChain->SetBranchAddress("ohBJetIPL3TagSingleTrack", ohBJetIPL3TagSingleTrack, &b_ohBJetIPL3TagSingleTrack);
   fChain->SetBranchAddress("ohBJetPerfL25Tag", ohBJetPerfL25Tag, &b_ohBJetPerfL25Tag);
   fChain->SetBranchAddress("ohBJetPerfL3Tag", ohBJetPerfL3Tag, &b_ohBJetPerfL3Tag);
   fChain->SetBranchAddress("NrecoElec", &NrecoElec, &b_NrecoElec);
   fChain->SetBranchAddress("recoElecPt", &recoElecPt, &b_recoElecPt);
   fChain->SetBranchAddress("recoElecPhi", &recoElecPhi, &b_recoElecPhi);
   fChain->SetBranchAddress("recoElecEta", &recoElecEta, &b_recoElecEta);
   fChain->SetBranchAddress("recoElecEt", &recoElecEt, &b_recoElecEt);
   fChain->SetBranchAddress("recoElecE", &recoElecE, &b_recoElecE);
   fChain->SetBranchAddress("recoElecEleID", &recoElecEleID, &b_recoElecEleID);
   fChain->SetBranchAddress("recoElecIP", &recoElecIP, &b_recoElecIP);
   fChain->SetBranchAddress("recoElecNLostHits", &recoElecNLostHits, &b_recoElecNLostHits);
   fChain->SetBranchAddress("recoElecChi2NDF", &recoElecChi2NDF, &b_recoElecChi2NDF);
   fChain->SetBranchAddress("recoElecTrkIsoR03", &recoElecTrkIsoR03, &b_recoElecTrkIsoR03);
   fChain->SetBranchAddress("recoElecECaloIsoR03", &recoElecECaloIsoR03, &b_recoElecECaloIsoR03);
   fChain->SetBranchAddress("recoElecHCaloIsoR03", &recoElecHCaloIsoR03, &b_recoElecHCaloIsoR03);
   fChain->SetBranchAddress("recoElecIsEcalDriven", &recoElecIsEcalDriven, &b_recoElecIsEcalDriven);
   fChain->SetBranchAddress("recoElecFbrem", &recoElecFbrem, &b_recoElecFbrem);
   fChain->SetBranchAddress("recoElecmishits", &recoElecmishits, &b_recoElecmishits);
   fChain->SetBranchAddress("recoElecdist", &recoElecdist, &b_recoElecdist);
   fChain->SetBranchAddress("recoElecdcot", &recoElecdcot, &b_recoElecdcot);
   fChain->SetBranchAddress("recoElectrkiso", &recoElectrkiso, &b_recoElectrkiso);
   fChain->SetBranchAddress("recoElececaliso", &recoElececaliso, &b_recoElececaliso);
   fChain->SetBranchAddress("recoElechcaliso", &recoElechcaliso, &b_recoElechcaliso);
   fChain->SetBranchAddress("recoElecsigmaietaieta", &recoElecsigmaietaieta, &b_recoElecsigmaietaieta);
   fChain->SetBranchAddress("recoElecdeltaPhiIn", &recoElecdeltaPhiIn, &b_recoElecdeltaPhiIn);
   fChain->SetBranchAddress("recoElecdeltaEtaIn", &recoElecdeltaEtaIn, &b_recoElecdeltaEtaIn);
   fChain->SetBranchAddress("recoElechOverE", &recoElechOverE, &b_recoElechOverE);
   fChain->SetBranchAddress("recoElecscEt", &recoElecscEt, &b_recoElecscEt);
   fChain->SetBranchAddress("recoElecd0corr", &recoElecd0corr, &b_recoElecd0corr);
   fChain->SetBranchAddress("recoElecqGsfCtfScPixConsistent", &recoElecqGsfCtfScPixConsistent, &b_recoElecqGsfCtfScPixConsistent);
   fChain->SetBranchAddress("NrecoPhot", &NrecoPhot, &b_NrecoPhot);
   fChain->SetBranchAddress("recoPhotPt", &recoPhotPt, &b_recoPhotPt);
   fChain->SetBranchAddress("recoPhotPhi", &recoPhotPhi, &b_recoPhotPhi);
   fChain->SetBranchAddress("recoPhotEta", &recoPhotEta, &b_recoPhotEta);
   fChain->SetBranchAddress("recoPhotEt", &recoPhotEt, &b_recoPhotEt);
   fChain->SetBranchAddress("recoPhotE", &recoPhotE, &b_recoPhotE);
   fChain->SetBranchAddress("recoPhotTiso", &recoPhotTiso, &b_recoPhotTiso);
   fChain->SetBranchAddress("recoPhotEiso", &recoPhotEiso, &b_recoPhotEiso);
   fChain->SetBranchAddress("recoPhotHiso", &recoPhotHiso, &b_recoPhotHiso);
   fChain->SetBranchAddress("recoPhotHoverE", &recoPhotHoverE, &b_recoPhotHoverE);
   fChain->SetBranchAddress("recoPhotClusShap", &recoPhotClusShap, &b_recoPhotClusShap);
   fChain->SetBranchAddress("recoPhotR9ID", &recoPhotR9ID, &b_recoPhotR9ID);
   fChain->SetBranchAddress("NohPhot", &NohPhot, &b_NohPhot);
   fChain->SetBranchAddress("ohPhotEt", ohPhotEt, &b_ohPhotEt);
   fChain->SetBranchAddress("ohPhotEta", ohPhotEta, &b_ohPhotEta);
   fChain->SetBranchAddress("ohPhotPhi", ohPhotPhi, &b_ohPhotPhi);
   fChain->SetBranchAddress("ohPhotEiso", ohPhotEiso, &b_ohPhotEiso);
   fChain->SetBranchAddress("ohPhotHiso", ohPhotHiso, &b_ohPhotHiso);
   fChain->SetBranchAddress("ohPhotTiso", ohPhotTiso, &b_ohPhotTiso);
   fChain->SetBranchAddress("ohPhotL1iso", ohPhotL1iso, &b_ohPhotL1iso);
   fChain->SetBranchAddress("ohPhotClusShap", ohPhotClusShap, &b_ohPhotClusShap);
   fChain->SetBranchAddress("ohPhotR9", ohPhotR9, &b_ohPhotR9);
   fChain->SetBranchAddress("ohPhotHforHoverE", ohPhotHforHoverE, &b_ohPhotHforHoverE);
   fChain->SetBranchAddress("ohPhotR9ID", ohPhotR9ID, &b_ohPhotR9ID);
   fChain->SetBranchAddress("NohEcalActiv", &NohEcalActiv, &b_NohEcalActiv);
   fChain->SetBranchAddress("ohEcalActivEt", ohEcalActivEt, &b_ohEcalActivEt);
   fChain->SetBranchAddress("ohEcalActivEta", ohEcalActivEta, &b_ohEcalActivEta);
   fChain->SetBranchAddress("ohEcalActivPhi", ohEcalActivPhi, &b_ohEcalActivPhi);
   fChain->SetBranchAddress("ohEcalActivEiso", ohEcalActivEiso, &b_ohEcalActivEiso);
   fChain->SetBranchAddress("ohEcalActivHiso", ohEcalActivHiso, &b_ohEcalActivHiso);
   fChain->SetBranchAddress("ohEcalActivTiso", ohEcalActivTiso, &b_ohEcalActivTiso);
   fChain->SetBranchAddress("ohEcalActivL1iso", ohEcalActivL1iso, &b_ohEcalActivL1iso);
   fChain->SetBranchAddress("ohEcalActivClusShap", ohEcalActivClusShap, &b_ohEcalActivClusShap);
   fChain->SetBranchAddress("ohEcalActivR9", ohEcalActivR9, &b_ohEcalActivR9);
   fChain->SetBranchAddress("ohEcalActivHforHoverE", ohEcalActivHforHoverE, &b_ohEcalActivHforHoverE);
   fChain->SetBranchAddress("ohEcalActivR9ID", ohEcalActivR9ID, &b_ohEcalActivR9ID);
   fChain->SetBranchAddress("NohEle", &NohEle, &b_NohEle);
   fChain->SetBranchAddress("ohEleEt", ohEleEt, &b_ohEleEt);
   fChain->SetBranchAddress("ohEleEta", ohEleEta, &b_ohEleEta);
   fChain->SetBranchAddress("ohElePhi", ohElePhi, &b_ohElePhi);
   fChain->SetBranchAddress("ohEleVtxZ", ohEleVtxZ, &b_ohEleVtxZ);
   fChain->SetBranchAddress("ohEleE", ohEleE, &b_ohEleE);
   fChain->SetBranchAddress("ohEleP", ohEleP, &b_ohEleP);
   fChain->SetBranchAddress("ohEleHiso", ohEleHiso, &b_ohEleHiso);
   fChain->SetBranchAddress("ohEleTiso", ohEleTiso, &b_ohEleTiso);
   fChain->SetBranchAddress("ohEleEiso", ohEleEiso, &b_ohEleEiso);
   fChain->SetBranchAddress("ohEleL1iso", ohEleL1iso, &b_ohEleL1iso);
   fChain->SetBranchAddress("ohElePixelSeeds", ohElePixelSeeds, &b_ohElePixelSeeds);
   fChain->SetBranchAddress("ohEleNewSC", ohEleNewSC, &b_ohEleNewSC);
   fChain->SetBranchAddress("ohEleClusShap", ohEleClusShap, &b_ohEleClusShap);
   fChain->SetBranchAddress("ohEleDeta", ohEleDeta, &b_ohEleDeta);
   fChain->SetBranchAddress("ohEleDphi", ohEleDphi, &b_ohEleDphi);
   fChain->SetBranchAddress("ohEleR9", ohEleR9, &b_ohEleR9);
   fChain->SetBranchAddress("ohEleHforHoverE", ohEleHforHoverE, &b_ohEleHforHoverE);
   fChain->SetBranchAddress("ohEleR9ID", ohEleR9ID, &b_ohEleR9ID);
   fChain->SetBranchAddress("NohHFEle", &NohHFEle, &b_NohHFEle);
   fChain->SetBranchAddress("ohHFElePt", ohHFElePt, &b_ohHFElePt);
   fChain->SetBranchAddress("ohHFEleEta", ohHFEleEta, &b_ohHFEleEta);
   fChain->SetBranchAddress("NohHFECALClus", &NohHFECALClus, &b_NohHFECALClus);
   fChain->SetBranchAddress("ohHFEleClustere9e25", ohHFEleClustere9e25, &b_ohHFEleClustere9e25);
   fChain->SetBranchAddress("ohHFEleClustere1e9", ohHFEleClustere1e9, &b_ohHFEleClustere1e9);
   fChain->SetBranchAddress("ohHFEleClustereCOREe9", ohHFEleClustereCOREe9, &b_ohHFEleClustereCOREe9);
   fChain->SetBranchAddress("ohHFEleClustereSeL", ohHFEleClustereSeL, &b_ohHFEleClustereSeL);
   fChain->SetBranchAddress("ohHFEleCluster2Dcut", ohHFEleCluster2Dcut, &b_ohHFEleCluster2Dcut);
   fChain->SetBranchAddress("ohHFEleClusterEta", ohHFEleClusterEta, &b_ohHFEleClusterEta);
   fChain->SetBranchAddress("ohHFEleClusterPhi", ohHFEleClusterPhi, &b_ohHFEleClusterPhi);
   fChain->SetBranchAddress("NrecoMuon", &NrecoMuon, &b_NrecoMuon);
   fChain->SetBranchAddress("recoMuonPt", &recoMuonPt, &b_recoMuonPt);
   fChain->SetBranchAddress("recoMuonPhi", &recoMuonPhi, &b_recoMuonPhi);
   fChain->SetBranchAddress("recoMuonEta", &recoMuonEta, &b_recoMuonEta);
   fChain->SetBranchAddress("recoMuonEt", &recoMuonEt, &b_recoMuonEt);
   fChain->SetBranchAddress("recoMuonE", &recoMuonE, &b_recoMuonE);
   fChain->SetBranchAddress("recoMuonChi2NDF", &recoMuonChi2NDF, &b_recoMuonChi2NDF);
   fChain->SetBranchAddress("recoMuonCharge", &recoMuonCharge, &b_recoMuonCharge);
   fChain->SetBranchAddress("recoMuonTrkIsoR03", &recoMuonTrkIsoR03, &b_recoMuonTrkIsoR03);
   fChain->SetBranchAddress("recoMuonECalIsoR03", &recoMuonECalIsoR03, &b_recoMuonECalIsoR03);
   fChain->SetBranchAddress("recoMuonHCalIsoR03", &recoMuonHCalIsoR03, &b_recoMuonHCalIsoR03);
   fChain->SetBranchAddress("recoMuonD0", &recoMuonD0, &b_recoMuonD0);
   fChain->SetBranchAddress("recoMuonType", &recoMuonType, &b_recoMuonType);
   fChain->SetBranchAddress("recoMuonNValidTrkHits", &recoMuonNValidTrkHits, &b_recoMuonNValidTrkHits);
   fChain->SetBranchAddress("recoMuonNValidMuonHits", &recoMuonNValidMuonHits, &b_recoMuonNValidMuonHits);
   fChain->SetBranchAddress("NohMuL2", &NohMuL2, &b_NohMuL2);
   fChain->SetBranchAddress("ohMuL2Pt", ohMuL2Pt, &b_ohMuL2Pt);
   fChain->SetBranchAddress("ohMuL2Phi", ohMuL2Phi, &b_ohMuL2Phi);
   fChain->SetBranchAddress("ohMuL2Eta", ohMuL2Eta, &b_ohMuL2Eta);
   fChain->SetBranchAddress("ohMuL2Chg", ohMuL2Chg, &b_ohMuL2Chg);
   fChain->SetBranchAddress("ohMuL2PtErr", ohMuL2PtErr, &b_ohMuL2PtErr);
   fChain->SetBranchAddress("ohMuL2Iso", ohMuL2Iso, &b_ohMuL2Iso);
   fChain->SetBranchAddress("ohMuL2Dr", ohMuL2Dr, &b_ohMuL2Dr);
   fChain->SetBranchAddress("ohMuL2Dz", ohMuL2Dz, &b_ohMuL2Dz);
   fChain->SetBranchAddress("ohMuL2VtxZ", ohMuL2VtxZ, &b_ohMuL2VtxZ);
   fChain->SetBranchAddress("ohMuL2Nhits", ohMuL2Nhits, &b_ohMuL2Nhits);
   fChain->SetBranchAddress("ohMuL2Nchambers", ohMuL2Nchambers, &b_ohMuL2Nchambers);
   fChain->SetBranchAddress("ohMuL2Nstat", ohMuL2Nstat, &b_ohMuL2Nstat);
   fChain->SetBranchAddress("ohMuL2L1idx", ohMuL2L1idx, &b_ohMuL2L1idx);
   fChain->SetBranchAddress("NohMuL3", &NohMuL3, &b_NohMuL3);
   fChain->SetBranchAddress("ohMuL3Pt", ohMuL3Pt, &b_ohMuL3Pt);
   fChain->SetBranchAddress("ohMuL3Phi", ohMuL3Phi, &b_ohMuL3Phi);
   fChain->SetBranchAddress("ohMuL3Eta", ohMuL3Eta, &b_ohMuL3Eta);
   fChain->SetBranchAddress("ohMuL3Chg", ohMuL3Chg, &b_ohMuL3Chg);
   fChain->SetBranchAddress("ohMuL3PtErr", ohMuL3PtErr, &b_ohMuL3PtErr);
   fChain->SetBranchAddress("ohMuL3Iso", ohMuL3Iso, &b_ohMuL3Iso);
   fChain->SetBranchAddress("ohMuL3Trk10Iso", ohMuL3Trk10Iso, &b_ohMuL3Trk10Iso);
   fChain->SetBranchAddress("ohMuL3Dr", ohMuL3Dr, &b_ohMuL3Dr);
   fChain->SetBranchAddress("ohMuL3Dz", ohMuL3Dz, &b_ohMuL3Dz);
   fChain->SetBranchAddress("ohMuL3VtxZ", ohMuL3VtxZ, &b_ohMuL3VtxZ);
   fChain->SetBranchAddress("ohMuL3Nhits", ohMuL3Nhits, &b_ohMuL3Nhits);
   fChain->SetBranchAddress("ohMuL3NormChi2", ohMuL3NormChi2, &b_ohMuL3NormChi2);
   fChain->SetBranchAddress("ohMuL3Ntrackerhits", ohMuL3Ntrackerhits, &b_ohMuL3Ntrackerhits);
   fChain->SetBranchAddress("ohMuL3Nmuonhits", ohMuL3Nmuonhits, &b_ohMuL3Nmuonhits);
   fChain->SetBranchAddress("ohMuL3L2idx", ohMuL3L2idx, &b_ohMuL3L2idx);
   fChain->SetBranchAddress("NohOniaPixel", &NohOniaPixel, &b_NohOniaPixel);
   fChain->SetBranchAddress("ohOniaPixelPt", ohOniaPixelPt, &b_ohOniaPixelPt);
   fChain->SetBranchAddress("ohOniaPixelPhi", ohOniaPixelPhi, &b_ohOniaPixelPhi);
   fChain->SetBranchAddress("ohOniaPixelEta", ohOniaPixelEta, &b_ohOniaPixelEta);
   fChain->SetBranchAddress("ohOniaPixelChg", ohOniaPixelChg, &b_ohOniaPixelChg);
   fChain->SetBranchAddress("ohOniaPixelDr", ohOniaPixelDr, &b_ohOniaPixelDr);
   fChain->SetBranchAddress("ohOniaPixelDz", ohOniaPixelDz, &b_ohOniaPixelDz);
   fChain->SetBranchAddress("ohOniaPixelHits", ohOniaPixelHits, &b_ohOniaPixelHits);
   fChain->SetBranchAddress("ohOniaPixelNormChi2", ohOniaPixelNormChi2, &b_ohOniaPixelNormChi2);
   fChain->SetBranchAddress("NohOniaTrack", &NohOniaTrack, &b_NohOniaTrack);
   fChain->SetBranchAddress("ohOniaTrackPt", ohOniaTrackPt, &b_ohOniaTrackPt);
   fChain->SetBranchAddress("ohOniaTrackPhi", ohOniaTrackPhi, &b_ohOniaTrackPhi);
   fChain->SetBranchAddress("ohOniaTrackEta", ohOniaTrackEta, &b_ohOniaTrackEta);
   fChain->SetBranchAddress("ohOniaTrackChg", ohOniaTrackChg, &b_ohOniaTrackChg);
   fChain->SetBranchAddress("ohOniaTrackDr", ohOniaTrackDr, &b_ohOniaTrackDr);
   fChain->SetBranchAddress("ohOniaTrackDz", ohOniaTrackDz, &b_ohOniaTrackDz);
   fChain->SetBranchAddress("ohOniaTrackHits", ohOniaTrackHits, &b_ohOniaTrackHits);
   fChain->SetBranchAddress("ohOniaTrackNormChi2", ohOniaTrackNormChi2, &b_ohOniaTrackNormChi2);
   fChain->SetBranchAddress("NohMuL2NoVtx", &NohMuL2NoVtx, &b_NohMuL2NoVtx);
   fChain->SetBranchAddress("ohMuL2NoVtxPt", ohMuL2NoVtxPt, &b_ohMuL2NoVtxPt);
   fChain->SetBranchAddress("ohMuL2NoVtxPhi", ohMuL2NoVtxPhi, &b_ohMuL2NoVtxPhi);
   fChain->SetBranchAddress("ohMuL2NoVtxEta", ohMuL2NoVtxEta, &b_ohMuL2NoVtxEta);
   fChain->SetBranchAddress("ohMuL2NoVtxChg", ohMuL2NoVtxChg, &b_ohMuL2NoVtxChg);
   fChain->SetBranchAddress("ohMuL2NoVtxPtErr", ohMuL2NoVtxPtErr, &b_ohMuL2NoVtxPtErr);
   fChain->SetBranchAddress("ohMuL2NoVtxDr", ohMuL2NoVtxDr, &b_ohMuL2NoVtxDr);
   fChain->SetBranchAddress("ohMuL2NoVtxDz", ohMuL2NoVtxDz, &b_ohMuL2NoVtxDz);
   fChain->SetBranchAddress("ohMuL2NoVtxNhits", ohMuL2NoVtxNhits, &b_ohMuL2NoVtxNhits);
   fChain->SetBranchAddress("ohMuL2NoVtxNchambers", ohMuL2NoVtxNchambers, &b_ohMuL2NoVtxNchambers);
   fChain->SetBranchAddress("ohMuL2NoVtxL1idx", ohMuL2NoVtxL1idx, &b_ohMuL2NoVtxL1idx);
   fChain->SetBranchAddress("NohDiMu", &NohDiMu, &b_NohDiMu);
   fChain->SetBranchAddress("ohDiMuDCA", ohDiMuDCA, &b_ohDiMuDCA);
   fChain->SetBranchAddress("ohDiMu1st", ohDiMu1st, &b_ohDiMu1st);
   fChain->SetBranchAddress("ohDiMu2nd", ohDiMu2nd, &b_ohDiMu2nd);
   fChain->SetBranchAddress("NohDiMuVtx", &NohDiMuVtx, &b_NohDiMuVtx);
   fChain->SetBranchAddress("ohDiMuVtx1st", &ohDiMuVtx1st, &b_ohDiMuVtx1st);
   fChain->SetBranchAddress("ohDiMuVtx2nd", &ohDiMuVtx2nd, &b_ohDiMuVtx2nd);
   fChain->SetBranchAddress("ohDiMuVtxChi2", &ohDiMuVtxChi2, &b_ohDiMuVtxChi2);
   fChain->SetBranchAddress("ohDiMuVtxR", &ohDiMuVtxR, &b_ohDiMuVtxR);
   fChain->SetBranchAddress("ohDiMuVtxRSig", &ohDiMuVtxRSig, &b_ohDiMuVtxRSig);
   fChain->SetBranchAddress("ohDiMuVtxROverSig", &ohDiMuVtxROverSig, &b_ohDiMuVtxROverSig);
   fChain->SetBranchAddress("ohDiMuVtxCosAlpha", &ohDiMuVtxCosAlpha, &b_ohDiMuVtxCosAlpha);
   fChain->SetBranchAddress("ohDiMuVtxMu2DIpMax", &ohDiMuVtxMu2DIpMax, &b_ohDiMuVtxMu2DIpMax);
   fChain->SetBranchAddress("ohDiMuVtxMu2DIpMin", &ohDiMuVtxMu2DIpMin, &b_ohDiMuVtxMu2DIpMin);
   fChain->SetBranchAddress("ohDiMuVtxMu2DIpSigMax", &ohDiMuVtxMu2DIpSigMax, &b_ohDiMuVtxMu2DIpSigMax);
   fChain->SetBranchAddress("ohDiMuVtxMu2DIpSigMin", &ohDiMuVtxMu2DIpSigMin, &b_ohDiMuVtxMu2DIpSigMin);
   fChain->SetBranchAddress("NohIsoPixelTrackL3", &NohIsoPixelTrackL3, &b_NohIsoPixelTrackL3);
   fChain->SetBranchAddress("ohIsoPixelTrackL3Pt", &ohIsoPixelTrackL3Pt, &b_ohIsoPixelTrackL3Pt);
   fChain->SetBranchAddress("ohIsoPixelTrackL3Eta", &ohIsoPixelTrackL3Eta, &b_ohIsoPixelTrackL3Eta);
   fChain->SetBranchAddress("ohIsoPixelTrackL3Phi", &ohIsoPixelTrackL3Phi, &b_ohIsoPixelTrackL3Phi);
   fChain->SetBranchAddress("ohIsoPixelTrackL3MaxPtPxl", &ohIsoPixelTrackL3MaxPtPxl, &b_ohIsoPixelTrackL3MaxPtPxl);
   fChain->SetBranchAddress("ohIsoPixelTrackL3Energy", &ohIsoPixelTrackL3Energy, &b_ohIsoPixelTrackL3Energy);
   fChain->SetBranchAddress("ohIsoPixelTrackL2pt", &ohIsoPixelTrackL2pt, &b_ohIsoPixelTrackL2pt);
   fChain->SetBranchAddress("ohIsoPixelTrackL2eta", &ohIsoPixelTrackL2eta, &b_ohIsoPixelTrackL2eta);
   fChain->SetBranchAddress("ohIsoPixelTrackL2dXY", &ohIsoPixelTrackL2dXY, &b_ohIsoPixelTrackL2dXY);
   fChain->SetBranchAddress("NohPixelTracksL3", &NohPixelTracksL3, &b_NohPixelTracksL3);
   fChain->SetBranchAddress("ohPixelTracksL3Pt", &ohPixelTracksL3Pt, &b_ohPixelTracksL3Pt);
   fChain->SetBranchAddress("ohPixelTracksL3Eta", &ohPixelTracksL3Eta, &b_ohPixelTracksL3Eta);
   fChain->SetBranchAddress("ohPixelTracksL3Phi", &ohPixelTracksL3Phi, &b_ohPixelTracksL3Phi);
   fChain->SetBranchAddress("ohPixelTracksL3Vz", &ohPixelTracksL3Vz, &b_ohPixelTracksL3Vz);
   fChain->SetBranchAddress("ohPixelFEDSize", &ohPixelFEDSize, &b_ohPixelFEDSize);
   fChain->SetBranchAddress("NohPixelClusters", &NohPixelClusters, &b_NohPixelClusters);
   fChain->SetBranchAddress("NMCpart", &NMCpart, &b_NMCpart);
   fChain->SetBranchAddress("MCpid", &MCpid, &b_MCpid);
   fChain->SetBranchAddress("MCstatus", &MCstatus, &b_MCstatus);
   fChain->SetBranchAddress("MCvtxX", &MCvtxX, &b_MCvtxX);
   fChain->SetBranchAddress("MCvtxY", &MCvtxY, &b_MCvtxY);
   fChain->SetBranchAddress("MCvtxZ", &MCvtxZ, &b_MCvtxZ);
   fChain->SetBranchAddress("MCpt", &MCpt, &b_MCpt);
   fChain->SetBranchAddress("MCeta", &MCeta, &b_MCeta);
   fChain->SetBranchAddress("MCphi", &MCphi, &b_MCphi);
   fChain->SetBranchAddress("MCPtHat", &MCPtHat, &b_MCPtHat);
   fChain->SetBranchAddress("MCmu3", &MCmu3, &b_MCmu3);
   fChain->SetBranchAddress("MCel3", &MCel3, &b_MCel3);
   fChain->SetBranchAddress("MCbb", &MCbb, &b_MCbb);
   fChain->SetBranchAddress("MCab", &MCab, &b_MCab);
   fChain->SetBranchAddress("MCWenu", &MCWenu, &b_MCWenu);
   fChain->SetBranchAddress("MCWmunu", &MCWmunu, &b_MCmunu);
   fChain->SetBranchAddress("MCZee", &MCZee, &b_MCZee);
   fChain->SetBranchAddress("MCZmumu", &MCZmumu, &b_MCZmumu);
   fChain->SetBranchAddress("MCptEleMax", &MCptEleMax, &b_MCptEleMax);
   fChain->SetBranchAddress("MCptMuMax", &MCptMuMax, &b_MCptMuMax);
   fChain->SetBranchAddress("NL1IsolEm", &NL1IsolEm, &b_NL1IsolEm);
   fChain->SetBranchAddress("L1IsolEmEt", L1IsolEmEt, &b_L1IsolEmEt);
   fChain->SetBranchAddress("L1IsolEmE", L1IsolEmE, &b_L1IsolEmE);
   fChain->SetBranchAddress("L1IsolEmEta", L1IsolEmEta, &b_L1IsolEmEta);
   fChain->SetBranchAddress("L1IsolEmPhi", L1IsolEmPhi, &b_L1IsolEmPhi);
   fChain->SetBranchAddress("NL1NIsolEm", &NL1NIsolEm, &b_NL1NIsolEm);
   fChain->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt, &b_L1NIsolEmEt);
   fChain->SetBranchAddress("L1NIsolEmE", L1NIsolEmE, &b_L1NIsolEmE);
   fChain->SetBranchAddress("L1NIsolEmEta", L1NIsolEmEta, &b_L1NIsolEmEta);
   fChain->SetBranchAddress("L1NIsolEmPhi", L1NIsolEmPhi, &b_L1NIsolEmPhi);
   fChain->SetBranchAddress("NL1Mu", &NL1Mu, &b_NL1Mu);
   fChain->SetBranchAddress("L1MuPt", L1MuPt, &b_L1MuPt);
   fChain->SetBranchAddress("L1MuE", L1MuE, &b_L1MuE);
   fChain->SetBranchAddress("L1MuEta", L1MuEta, &b_L1MuEta);
   fChain->SetBranchAddress("L1MuPhi", L1MuPhi, &b_L1MuPhi);
   fChain->SetBranchAddress("L1MuIsol", L1MuIsol, &b_L1MuIsol);
   fChain->SetBranchAddress("L1MuMip", L1MuMip, &b_L1MuMip);
   fChain->SetBranchAddress("L1MuFor", L1MuFor, &b_L1MuFor);
   fChain->SetBranchAddress("L1MuRPC", L1MuRPC, &b_L1MuRPC);
   fChain->SetBranchAddress("L1MuQal", L1MuQal, &b_L1MuQal);
   fChain->SetBranchAddress("L1MuChg", L1MuChg, &b_L1MuChg);
   fChain->SetBranchAddress("NL1CenJet", &NL1CenJet, &b_NL1CenJet);
   fChain->SetBranchAddress("L1CenJetEt", L1CenJetEt, &b_L1CenJetEt);
   fChain->SetBranchAddress("L1CenJetE", L1CenJetE, &b_L1CenJetE);
   fChain->SetBranchAddress("L1CenJetEta", L1CenJetEta, &b_L1CenJetEta);
   fChain->SetBranchAddress("L1CenJetPhi", L1CenJetPhi, &b_L1CenJetPhi);
   fChain->SetBranchAddress("NL1ForJet", &NL1ForJet, &b_NL1ForJet);
   fChain->SetBranchAddress("L1ForJetEt", L1ForJetEt, &b_L1ForJetEt);
   fChain->SetBranchAddress("L1ForJetE", L1ForJetE, &b_L1ForJetE);
   fChain->SetBranchAddress("L1ForJetEta", L1ForJetEta, &b_L1ForJetEta);
   fChain->SetBranchAddress("L1ForJetPhi", L1ForJetPhi, &b_L1ForJetPhi);
   fChain->SetBranchAddress("NL1Tau", &NL1Tau, &b_NL1Tau);
   fChain->SetBranchAddress("L1TauEt", L1TauEt, &b_L1TauEt);
   fChain->SetBranchAddress("L1TauE", L1TauE, &b_L1TauE);
   fChain->SetBranchAddress("L1TauEta", L1TauEta, &b_L1TauEta);
   fChain->SetBranchAddress("L1TauPhi", L1TauPhi, &b_L1TauPhi);
   fChain->SetBranchAddress("L1Met", &L1Met, &b_L1Met);
   fChain->SetBranchAddress("L1MetPhi", &L1MetPhi, &b_L1MetPhi);
   fChain->SetBranchAddress("L1EtTot", &L1EtTot, &b_L1EtTot);
   fChain->SetBranchAddress("L1Mht", &L1Mht, &b_L1Mht);
   fChain->SetBranchAddress("L1MhtPhi", &L1MhtPhi, &b_L1MhtPhi);
   fChain->SetBranchAddress("L1EtHad", &L1EtHad, &b_L1EtHad);
   fChain->SetBranchAddress("L1HfRing1EtSumPositiveEta", &L1HfRing1EtSumPositiveEta, &b_L1HfRing1EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing2EtSumPositiveEta", &L1HfRing2EtSumPositiveEta, &b_L1HfRing2EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing1EtSumNegativeEta", &L1HfRing1EtSumNegativeEta, &b_L1HfRing1EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfRing2EtSumNegativeEta", &L1HfRing2EtSumNegativeEta, &b_L1HfRing2EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing1", &L1HfTowerCountPositiveEtaRing1, &b_L1HfTowerCountPositiveEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing1", &L1HfTowerCountNegativeEtaRing1, &b_L1HfTowerCountNegativeEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing2", &L1HfTowerCountPositiveEtaRing2, &b_L1HfTowerCountPositiveEtaRing2);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing2", &L1HfTowerCountNegativeEtaRing2, &b_L1HfTowerCountNegativeEtaRing2);
   fChain->SetBranchAddress("recoNVrtOffline0", &recoNVrtOffline0, &b_NVrtx);
   fChain->SetBranchAddress("recoVrtXOffline0", &recoVrtXOffline0, &b_recoVrtXOffline0);
   fChain->SetBranchAddress("recoVrtYOffline0", &recoVrtYOffline0, &b_recoVrtYOffline0);
   fChain->SetBranchAddress("recoVrtZOffline0", &recoVrtZOffline0, &b_recoVrtZOffline0);
   fChain->SetBranchAddress("recoVrtNtrkOffline0", &recoVrtNtrkOffline0, &b_recoVrtNtrkOffline0);
   fChain->SetBranchAddress("recoVrtChi2Offline0", &recoVrtChi2Offline0, &b_recoVrtChi2Offline0);
   fChain->SetBranchAddress("recoVrtNdofOffline0", &recoVrtNdofOffline0, &b_recoVrtNdofOffline0);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Bx", &Bx, &b_Bx);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("AvgInstDelLumi", &AvgInstDelLumi, &b_AvgInstDelLumi);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC7_v1", &HLT_Activity_Ecal_SC7_v1, &b_HLT_Activity_Ecal_SC7_v1);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC7_v1_Prescl", &HLT_Activity_Ecal_SC7_v1_Prescl, &b_HLT_Activity_Ecal_SC7_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleJet36_v1", &HLT_L1SingleJet36_v1, &b_HLT_L1SingleJet36_v1);
   fChain->SetBranchAddress("HLT_L1SingleJet36_v1_Prescl", &HLT_L1SingleJet36_v1_Prescl, &b_HLT_L1SingleJet36_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet30_v1", &HLT_Jet30_v1, &b_HLT_Jet30_v1);
   fChain->SetBranchAddress("HLT_Jet30_v1_Prescl", &HLT_Jet30_v1_Prescl, &b_HLT_Jet30_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet60_v1", &HLT_Jet60_v1, &b_HLT_Jet60_v1);
   fChain->SetBranchAddress("HLT_Jet60_v1_Prescl", &HLT_Jet60_v1_Prescl, &b_HLT_Jet60_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet80_v1", &HLT_Jet80_v1, &b_HLT_Jet80_v1);
   fChain->SetBranchAddress("HLT_Jet80_v1_Prescl", &HLT_Jet80_v1_Prescl, &b_HLT_Jet80_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet110_v1", &HLT_Jet110_v1, &b_HLT_Jet110_v1);
   fChain->SetBranchAddress("HLT_Jet110_v1_Prescl", &HLT_Jet110_v1_Prescl, &b_HLT_Jet110_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet150_v1", &HLT_Jet150_v1, &b_HLT_Jet150_v1);
   fChain->SetBranchAddress("HLT_Jet150_v1_Prescl", &HLT_Jet150_v1_Prescl, &b_HLT_Jet150_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet190_v1", &HLT_Jet190_v1, &b_HLT_Jet190_v1);
   fChain->SetBranchAddress("HLT_Jet190_v1_Prescl", &HLT_Jet190_v1_Prescl, &b_HLT_Jet190_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet240_v1", &HLT_Jet240_v1, &b_HLT_Jet240_v1);
   fChain->SetBranchAddress("HLT_Jet240_v1_Prescl", &HLT_Jet240_v1_Prescl, &b_HLT_Jet240_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet370_v1", &HLT_Jet370_v1, &b_HLT_Jet370_v1);
   fChain->SetBranchAddress("HLT_Jet370_v1_Prescl", &HLT_Jet370_v1_Prescl, &b_HLT_Jet370_v1_Prescl);
   fChain->SetBranchAddress("HLT_Jet370_NoJetID_v1", &HLT_Jet370_NoJetID_v1, &b_HLT_Jet370_NoJetID_v1);
   fChain->SetBranchAddress("HLT_Jet370_NoJetID_v1_Prescl", &HLT_Jet370_NoJetID_v1_Prescl, &b_HLT_Jet370_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve15U_v4", &HLT_DiJetAve15U_v4, &b_HLT_DiJetAve15U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve15U_v4_Prescl", &HLT_DiJetAve15U_v4_Prescl, &b_HLT_DiJetAve15U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve30U_v4", &HLT_DiJetAve30U_v4, &b_HLT_DiJetAve30U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve30U_v4_Prescl", &HLT_DiJetAve30U_v4_Prescl, &b_HLT_DiJetAve30U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve50U_v4", &HLT_DiJetAve50U_v4, &b_HLT_DiJetAve50U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve50U_v4_Prescl", &HLT_DiJetAve50U_v4_Prescl, &b_HLT_DiJetAve50U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve70U_v4", &HLT_DiJetAve70U_v4, &b_HLT_DiJetAve70U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve70U_v4_Prescl", &HLT_DiJetAve70U_v4_Prescl, &b_HLT_DiJetAve70U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve100U_v4", &HLT_DiJetAve100U_v4, &b_HLT_DiJetAve100U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve100U_v4_Prescl", &HLT_DiJetAve100U_v4_Prescl, &b_HLT_DiJetAve100U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve140U_v4", &HLT_DiJetAve140U_v4, &b_HLT_DiJetAve140U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve140U_v4_Prescl", &HLT_DiJetAve140U_v4_Prescl, &b_HLT_DiJetAve140U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve180U_v4", &HLT_DiJetAve180U_v4, &b_HLT_DiJetAve180U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve180U_v4_Prescl", &HLT_DiJetAve180U_v4_Prescl, &b_HLT_DiJetAve180U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve300U_v4", &HLT_DiJetAve300U_v4, &b_HLT_DiJetAve300U_v4);
   fChain->SetBranchAddress("HLT_DiJetAve300U_v4_Prescl", &HLT_DiJetAve300U_v4_Prescl, &b_HLT_DiJetAve300U_v4_Prescl);
   fChain->SetBranchAddress("HLT_DoubleJet30_ForwardBackward_v2", &HLT_DoubleJet30_ForwardBackward_v2, &b_HLT_DoubleJet30_ForwardBackward_v2);
   fChain->SetBranchAddress("HLT_DoubleJet30_ForwardBackward_v2_Prescl", &HLT_DoubleJet30_ForwardBackward_v2_Prescl, &b_HLT_DoubleJet30_ForwardBackward_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleJet60_ForwardBackward_v2", &HLT_DoubleJet60_ForwardBackward_v2, &b_HLT_DoubleJet60_ForwardBackward_v2);
   fChain->SetBranchAddress("HLT_DoubleJet60_ForwardBackward_v2_Prescl", &HLT_DoubleJet60_ForwardBackward_v2_Prescl, &b_HLT_DoubleJet60_ForwardBackward_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleJet70_ForwardBackward_v2", &HLT_DoubleJet70_ForwardBackward_v2, &b_HLT_DoubleJet70_ForwardBackward_v2);
   fChain->SetBranchAddress("HLT_DoubleJet70_ForwardBackward_v2_Prescl", &HLT_DoubleJet70_ForwardBackward_v2_Prescl, &b_HLT_DoubleJet70_ForwardBackward_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleJet80_ForwardBackward_v2", &HLT_DoubleJet80_ForwardBackward_v2, &b_HLT_DoubleJet80_ForwardBackward_v2);
   fChain->SetBranchAddress("HLT_DoubleJet80_ForwardBackward_v2_Prescl", &HLT_DoubleJet80_ForwardBackward_v2_Prescl, &b_HLT_DoubleJet80_ForwardBackward_v2_Prescl);
   fChain->SetBranchAddress("HLT_CentralJet80_MET65_v1", &HLT_CentralJet80_MET65_v1, &b_HLT_CentralJet80_MET65_v1);
   fChain->SetBranchAddress("HLT_CentralJet80_MET65_v1_Prescl", &HLT_CentralJet80_MET65_v1_Prescl, &b_HLT_CentralJet80_MET65_v1_Prescl);
   fChain->SetBranchAddress("HLT_CentralJet80_MET80_v1", &HLT_CentralJet80_MET80_v1, &b_HLT_CentralJet80_MET80_v1);
   fChain->SetBranchAddress("HLT_CentralJet80_MET80_v1_Prescl", &HLT_CentralJet80_MET80_v1_Prescl, &b_HLT_CentralJet80_MET80_v1_Prescl);
   fChain->SetBranchAddress("HLT_CentralJet80_MET100_v1", &HLT_CentralJet80_MET100_v1, &b_HLT_CentralJet80_MET100_v1);
   fChain->SetBranchAddress("HLT_CentralJet80_MET100_v1_Prescl", &HLT_CentralJet80_MET100_v1_Prescl, &b_HLT_CentralJet80_MET100_v1_Prescl);
   fChain->SetBranchAddress("HLT_CentralJet80_MET160_v1", &HLT_CentralJet80_MET160_v1, &b_HLT_CentralJet80_MET160_v1);
   fChain->SetBranchAddress("HLT_CentralJet80_MET160_v1_Prescl", &HLT_CentralJet80_MET160_v1_Prescl, &b_HLT_CentralJet80_MET160_v1_Prescl);
   fChain->SetBranchAddress("HLT_DiJet60_MET45_v1", &HLT_DiJet60_MET45_v1, &b_HLT_DiJet60_MET45_v1);
   fChain->SetBranchAddress("HLT_DiJet60_MET45_v1_Prescl", &HLT_DiJet60_MET45_v1_Prescl, &b_HLT_DiJet60_MET45_v1_Prescl);
   fChain->SetBranchAddress("HLT_DiJet70_PT70_v1", &HLT_DiJet70_PT70_v1, &b_HLT_DiJet70_PT70_v1);
   fChain->SetBranchAddress("HLT_DiJet70_PT70_v1_Prescl", &HLT_DiJet70_PT70_v1_Prescl, &b_HLT_DiJet70_PT70_v1_Prescl);
   fChain->SetBranchAddress("HLT_DiJet100_PT100_v1", &HLT_DiJet100_PT100_v1, &b_HLT_DiJet100_PT100_v1);
   fChain->SetBranchAddress("HLT_DiJet100_PT100_v1_Prescl", &HLT_DiJet100_PT100_v1_Prescl, &b_HLT_DiJet100_PT100_v1_Prescl);
   fChain->SetBranchAddress("HLT_DiJet130_PT130_v1", &HLT_DiJet130_PT130_v1, &b_HLT_DiJet130_PT130_v1);
   fChain->SetBranchAddress("HLT_DiJet130_PT130_v1_Prescl", &HLT_DiJet130_PT130_v1_Prescl, &b_HLT_DiJet130_PT130_v1_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet40_v2", &HLT_QuadJet40_v2, &b_HLT_QuadJet40_v2);
   fChain->SetBranchAddress("HLT_QuadJet40_v2_Prescl", &HLT_QuadJet40_v2_Prescl, &b_HLT_QuadJet40_v2_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet40_IsoPFTau40_v1", &HLT_QuadJet40_IsoPFTau40_v1, &b_HLT_QuadJet40_IsoPFTau40_v1);
   fChain->SetBranchAddress("HLT_QuadJet40_IsoPFTau40_v1_Prescl", &HLT_QuadJet40_IsoPFTau40_v1_Prescl, &b_HLT_QuadJet40_IsoPFTau40_v1_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet50_BTagIP_v1", &HLT_QuadJet50_BTagIP_v1, &b_HLT_QuadJet50_BTagIP_v1);
   fChain->SetBranchAddress("HLT_QuadJet50_BTagIP_v1_Prescl", &HLT_QuadJet50_BTagIP_v1_Prescl, &b_HLT_QuadJet50_BTagIP_v1_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet50_Jet40_v1", &HLT_QuadJet50_Jet40_v1, &b_HLT_QuadJet50_Jet40_v1);
   fChain->SetBranchAddress("HLT_QuadJet50_Jet40_v1_Prescl", &HLT_QuadJet50_Jet40_v1_Prescl, &b_HLT_QuadJet50_Jet40_v1_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet60_v1", &HLT_QuadJet60_v1, &b_HLT_QuadJet60_v1);
   fChain->SetBranchAddress("HLT_QuadJet60_v1_Prescl", &HLT_QuadJet60_v1_Prescl, &b_HLT_QuadJet60_v1_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet70_v1", &HLT_QuadJet70_v1, &b_HLT_QuadJet70_v1);
   fChain->SetBranchAddress("HLT_QuadJet70_v1_Prescl", &HLT_QuadJet70_v1_Prescl, &b_HLT_QuadJet70_v1_Prescl);
   fChain->SetBranchAddress("HLT_ExclDiJet60_HFOR_v1", &HLT_ExclDiJet60_HFOR_v1, &b_HLT_ExclDiJet60_HFOR_v1);
   fChain->SetBranchAddress("HLT_ExclDiJet60_HFOR_v1_Prescl", &HLT_ExclDiJet60_HFOR_v1_Prescl, &b_HLT_ExclDiJet60_HFOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_ExclDiJet60_HFAND_v1", &HLT_ExclDiJet60_HFAND_v1, &b_HLT_ExclDiJet60_HFAND_v1);
   fChain->SetBranchAddress("HLT_ExclDiJet60_HFAND_v1_Prescl", &HLT_ExclDiJet60_HFAND_v1_Prescl, &b_HLT_ExclDiJet60_HFAND_v1_Prescl);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX_v2", &HLT_JetE30_NoBPTX_v2, &b_HLT_JetE30_NoBPTX_v2);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX_v2_Prescl", &HLT_JetE30_NoBPTX_v2_Prescl, &b_HLT_JetE30_NoBPTX_v2_Prescl);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX_NoHalo_v4", &HLT_JetE30_NoBPTX_NoHalo_v4, &b_HLT_JetE30_NoBPTX_NoHalo_v4);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX_NoHalo_v4_Prescl", &HLT_JetE30_NoBPTX_NoHalo_v4_Prescl, &b_HLT_JetE30_NoBPTX_NoHalo_v4_Prescl);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX3BX_NoHalo_v4", &HLT_JetE30_NoBPTX3BX_NoHalo_v4, &b_HLT_JetE30_NoBPTX3BX_NoHalo_v4);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX3BX_NoHalo_v4_Prescl", &HLT_JetE30_NoBPTX3BX_NoHalo_v4_Prescl, &b_HLT_JetE30_NoBPTX3BX_NoHalo_v4_Prescl);
   fChain->SetBranchAddress("HLT_HT150_v2", &HLT_HT150_v2, &b_HLT_HT150_v2);
   fChain->SetBranchAddress("HLT_HT150_v2_Prescl", &HLT_HT150_v2_Prescl, &b_HLT_HT150_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT150_AlphaT0p60_v1", &HLT_HT150_AlphaT0p60_v1, &b_HLT_HT150_AlphaT0p60_v1);
   fChain->SetBranchAddress("HLT_HT150_AlphaT0p60_v1_Prescl", &HLT_HT150_AlphaT0p60_v1_Prescl, &b_HLT_HT150_AlphaT0p60_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT150_AlphaT0p70_v1", &HLT_HT150_AlphaT0p70_v1, &b_HLT_HT150_AlphaT0p70_v1);
   fChain->SetBranchAddress("HLT_HT150_AlphaT0p70_v1_Prescl", &HLT_HT150_AlphaT0p70_v1_Prescl, &b_HLT_HT150_AlphaT0p70_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT200_v2", &HLT_HT200_v2, &b_HLT_HT200_v2);
   fChain->SetBranchAddress("HLT_HT200_v2_Prescl", &HLT_HT200_v2_Prescl, &b_HLT_HT200_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT200_AlphaT0p60_v1", &HLT_HT200_AlphaT0p60_v1, &b_HLT_HT200_AlphaT0p60_v1);
   fChain->SetBranchAddress("HLT_HT200_AlphaT0p60_v1_Prescl", &HLT_HT200_AlphaT0p60_v1_Prescl, &b_HLT_HT200_AlphaT0p60_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT200_AlphaT0p65_v1", &HLT_HT200_AlphaT0p65_v1, &b_HLT_HT200_AlphaT0p65_v1);
   fChain->SetBranchAddress("HLT_HT200_AlphaT0p65_v1_Prescl", &HLT_HT200_AlphaT0p65_v1_Prescl, &b_HLT_HT200_AlphaT0p65_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT250_v2", &HLT_HT250_v2, &b_HLT_HT250_v2);
   fChain->SetBranchAddress("HLT_HT250_v2_Prescl", &HLT_HT250_v2_Prescl, &b_HLT_HT250_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT250_AlphaT0p55_v1", &HLT_HT250_AlphaT0p55_v1, &b_HLT_HT250_AlphaT0p55_v1);
   fChain->SetBranchAddress("HLT_HT250_AlphaT0p55_v1_Prescl", &HLT_HT250_AlphaT0p55_v1_Prescl, &b_HLT_HT250_AlphaT0p55_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT250_AlphaT0p62_v1", &HLT_HT250_AlphaT0p62_v1, &b_HLT_HT250_AlphaT0p62_v1);
   fChain->SetBranchAddress("HLT_HT250_AlphaT0p62_v1_Prescl", &HLT_HT250_AlphaT0p62_v1_Prescl, &b_HLT_HT250_AlphaT0p62_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT250_DoubleDisplacedJet60_v1", &HLT_HT250_DoubleDisplacedJet60_v1, &b_HLT_HT250_DoubleDisplacedJet60_v1);
   fChain->SetBranchAddress("HLT_HT250_DoubleDisplacedJet60_v1_Prescl", &HLT_HT250_DoubleDisplacedJet60_v1_Prescl, &b_HLT_HT250_DoubleDisplacedJet60_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT250_MHT60_v2", &HLT_HT250_MHT60_v2, &b_HLT_HT250_MHT60_v2);
   fChain->SetBranchAddress("HLT_HT250_MHT60_v2_Prescl", &HLT_HT250_MHT60_v2_Prescl, &b_HLT_HT250_MHT60_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT300_v3", &HLT_HT300_v3, &b_HLT_HT300_v3);
   fChain->SetBranchAddress("HLT_HT300_v3_Prescl", &HLT_HT300_v3_Prescl, &b_HLT_HT300_v3_Prescl);
   fChain->SetBranchAddress("HLT_HT300_MHT75_v3", &HLT_HT300_MHT75_v3, &b_HLT_HT300_MHT75_v3);
   fChain->SetBranchAddress("HLT_HT300_MHT75_v3_Prescl", &HLT_HT300_MHT75_v3_Prescl, &b_HLT_HT300_MHT75_v3_Prescl);
   fChain->SetBranchAddress("HLT_HT300_AlphaT0p52_v1", &HLT_HT300_AlphaT0p52_v1, &b_HLT_HT300_AlphaT0p52_v1);
   fChain->SetBranchAddress("HLT_HT300_AlphaT0p52_v1_Prescl", &HLT_HT300_AlphaT0p52_v1_Prescl, &b_HLT_HT300_AlphaT0p52_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT300_AlphaT0p54_v1", &HLT_HT300_AlphaT0p54_v1, &b_HLT_HT300_AlphaT0p54_v1);
   fChain->SetBranchAddress("HLT_HT300_AlphaT0p54_v1_Prescl", &HLT_HT300_AlphaT0p54_v1_Prescl, &b_HLT_HT300_AlphaT0p54_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT350_v2", &HLT_HT350_v2, &b_HLT_HT350_v2);
   fChain->SetBranchAddress("HLT_HT350_v2_Prescl", &HLT_HT350_v2_Prescl, &b_HLT_HT350_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT350_AlphaT0p51_v1", &HLT_HT350_AlphaT0p51_v1, &b_HLT_HT350_AlphaT0p51_v1);
   fChain->SetBranchAddress("HLT_HT350_AlphaT0p51_v1_Prescl", &HLT_HT350_AlphaT0p51_v1_Prescl, &b_HLT_HT350_AlphaT0p51_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT350_AlphaT0p53_v1", &HLT_HT350_AlphaT0p53_v1, &b_HLT_HT350_AlphaT0p53_v1);
   fChain->SetBranchAddress("HLT_HT350_AlphaT0p53_v1_Prescl", &HLT_HT350_AlphaT0p53_v1_Prescl, &b_HLT_HT350_AlphaT0p53_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT400_v2", &HLT_HT400_v2, &b_HLT_HT400_v2);
   fChain->SetBranchAddress("HLT_HT400_v2_Prescl", &HLT_HT400_v2_Prescl, &b_HLT_HT400_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT400_AlphaT0p51_v1", &HLT_HT400_AlphaT0p51_v1, &b_HLT_HT400_AlphaT0p51_v1);
   fChain->SetBranchAddress("HLT_HT400_AlphaT0p51_v1_Prescl", &HLT_HT400_AlphaT0p51_v1_Prescl, &b_HLT_HT400_AlphaT0p51_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT450_v2", &HLT_HT450_v2, &b_HLT_HT450_v2);
   fChain->SetBranchAddress("HLT_HT450_v2_Prescl", &HLT_HT450_v2_Prescl, &b_HLT_HT450_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT500_v2", &HLT_HT500_v2, &b_HLT_HT500_v2);
   fChain->SetBranchAddress("HLT_HT500_v2_Prescl", &HLT_HT500_v2_Prescl, &b_HLT_HT500_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT550_v2", &HLT_HT550_v2, &b_HLT_HT550_v2);
   fChain->SetBranchAddress("HLT_HT550_v2_Prescl", &HLT_HT550_v2_Prescl, &b_HLT_HT550_v2_Prescl);
   fChain->SetBranchAddress("HLT_PFMHT150_v2", &HLT_PFMHT150_v2, &b_HLT_PFMHT150_v2);
   fChain->SetBranchAddress("HLT_PFMHT150_v2_Prescl", &HLT_PFMHT150_v2_Prescl, &b_HLT_PFMHT150_v2_Prescl);
   fChain->SetBranchAddress("HLT_MET100_v1", &HLT_MET100_v1, &b_HLT_MET100_v1);
   fChain->SetBranchAddress("HLT_MET100_v1_Prescl", &HLT_MET100_v1_Prescl, &b_HLT_MET100_v1_Prescl);
   fChain->SetBranchAddress("HLT_MET120_v1", &HLT_MET120_v1, &b_HLT_MET120_v1);
   fChain->SetBranchAddress("HLT_MET120_v1_Prescl", &HLT_MET120_v1_Prescl, &b_HLT_MET120_v1_Prescl);
   fChain->SetBranchAddress("HLT_MET200_v1", &HLT_MET200_v1, &b_HLT_MET200_v1);
   fChain->SetBranchAddress("HLT_MET200_v1_Prescl", &HLT_MET200_v1_Prescl, &b_HLT_MET200_v1_Prescl);
   fChain->SetBranchAddress("HLT_Meff440_v2", &HLT_Meff440_v2, &b_HLT_Meff440_v2);
   fChain->SetBranchAddress("HLT_Meff440_v2_Prescl", &HLT_Meff440_v2_Prescl, &b_HLT_Meff440_v2_Prescl);
   fChain->SetBranchAddress("HLT_Meff520_v2", &HLT_Meff520_v2, &b_HLT_Meff520_v2);
   fChain->SetBranchAddress("HLT_Meff520_v2_Prescl", &HLT_Meff520_v2_Prescl, &b_HLT_Meff520_v2_Prescl);
   fChain->SetBranchAddress("HLT_Meff640_v2", &HLT_Meff640_v2, &b_HLT_Meff640_v2);
   fChain->SetBranchAddress("HLT_Meff640_v2_Prescl", &HLT_Meff640_v2_Prescl, &b_HLT_Meff640_v2_Prescl);
   fChain->SetBranchAddress("HLT_MR100_v1", &HLT_MR100_v1, &b_HLT_MR100_v1);
   fChain->SetBranchAddress("HLT_MR100_v1_Prescl", &HLT_MR100_v1_Prescl, &b_HLT_MR100_v1_Prescl);
   fChain->SetBranchAddress("HLT_R032_v1", &HLT_R032_v1, &b_HLT_R032_v1);
   fChain->SetBranchAddress("HLT_R032_v1_Prescl", &HLT_R032_v1_Prescl, &b_HLT_R032_v1_Prescl);
   fChain->SetBranchAddress("HLT_R032_MR100_v1", &HLT_R032_MR100_v1, &b_HLT_R032_MR100_v1);
   fChain->SetBranchAddress("HLT_R032_MR100_v1_Prescl", &HLT_R032_MR100_v1_Prescl, &b_HLT_R032_MR100_v1_Prescl);
   fChain->SetBranchAddress("HLT_R035_MR100_v1", &HLT_R035_MR100_v1, &b_HLT_R035_MR100_v1);
   fChain->SetBranchAddress("HLT_R035_MR100_v1_Prescl", &HLT_R035_MR100_v1_Prescl, &b_HLT_R035_MR100_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_v1", &HLT_L1SingleMuOpen_v1, &b_HLT_L1SingleMuOpen_v1);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_v1_Prescl", &HLT_L1SingleMuOpen_v1_Prescl, &b_HLT_L1SingleMuOpen_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_DT_v1", &HLT_L1SingleMuOpen_DT_v1, &b_HLT_L1SingleMuOpen_DT_v1);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_DT_v1_Prescl", &HLT_L1SingleMuOpen_DT_v1_Prescl, &b_HLT_L1SingleMuOpen_DT_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleMu10_v1", &HLT_L1SingleMu10_v1, &b_HLT_L1SingleMu10_v1);
   fChain->SetBranchAddress("HLT_L1SingleMu10_v1_Prescl", &HLT_L1SingleMu10_v1_Prescl, &b_HLT_L1SingleMu10_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleMu20_v1", &HLT_L1SingleMu20_v1, &b_HLT_L1SingleMu20_v1);
   fChain->SetBranchAddress("HLT_L1SingleMu20_v1_Prescl", &HLT_L1SingleMu20_v1_Prescl, &b_HLT_L1SingleMu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1DoubleMu0_v1", &HLT_L1DoubleMu0_v1, &b_HLT_L1DoubleMu0_v1);
   fChain->SetBranchAddress("HLT_L1DoubleMu0_v1_Prescl", &HLT_L1DoubleMu0_v1_Prescl, &b_HLT_L1DoubleMu0_v1_Prescl);
   fChain->SetBranchAddress("HLT_L2Mu10_v1", &HLT_L2Mu10_v1, &b_HLT_L2Mu10_v1);
   fChain->SetBranchAddress("HLT_L2Mu10_v1_Prescl", &HLT_L2Mu10_v1_Prescl, &b_HLT_L2Mu10_v1_Prescl);
   fChain->SetBranchAddress("HLT_L2Mu20_v1", &HLT_L2Mu20_v1, &b_HLT_L2Mu20_v1);
   fChain->SetBranchAddress("HLT_L2Mu20_v1_Prescl", &HLT_L2Mu20_v1_Prescl, &b_HLT_L2Mu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_L2DoubleMu0_v2", &HLT_L2DoubleMu0_v2, &b_HLT_L2DoubleMu0_v2);
   fChain->SetBranchAddress("HLT_L2DoubleMu0_v2_Prescl", &HLT_L2DoubleMu0_v2_Prescl, &b_HLT_L2DoubleMu0_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_v3", &HLT_Mu3_v3, &b_HLT_Mu3_v3);
   fChain->SetBranchAddress("HLT_Mu3_v3_Prescl", &HLT_Mu3_v3_Prescl, &b_HLT_Mu3_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_v3", &HLT_Mu5_v3, &b_HLT_Mu5_v3);
   fChain->SetBranchAddress("HLT_Mu5_v3_Prescl", &HLT_Mu5_v3_Prescl, &b_HLT_Mu5_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu8_v1", &HLT_Mu8_v1, &b_HLT_Mu8_v1);
   fChain->SetBranchAddress("HLT_Mu8_v1_Prescl", &HLT_Mu8_v1_Prescl, &b_HLT_Mu8_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu12_v1", &HLT_Mu12_v1, &b_HLT_Mu12_v1);
   fChain->SetBranchAddress("HLT_Mu12_v1_Prescl", &HLT_Mu12_v1_Prescl, &b_HLT_Mu12_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu15_v2", &HLT_Mu15_v2, &b_HLT_Mu15_v2);
   fChain->SetBranchAddress("HLT_Mu15_v2_Prescl", &HLT_Mu15_v2_Prescl, &b_HLT_Mu15_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu20_v1", &HLT_Mu20_v1, &b_HLT_Mu20_v1);
   fChain->SetBranchAddress("HLT_Mu20_v1_Prescl", &HLT_Mu20_v1_Prescl, &b_HLT_Mu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu24_v1", &HLT_Mu24_v1, &b_HLT_Mu24_v1);
   fChain->SetBranchAddress("HLT_Mu24_v1_Prescl", &HLT_Mu24_v1_Prescl, &b_HLT_Mu24_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu30_v1", &HLT_Mu30_v1, &b_HLT_Mu30_v1);
   fChain->SetBranchAddress("HLT_Mu30_v1_Prescl", &HLT_Mu30_v1_Prescl, &b_HLT_Mu30_v1_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu12_v1", &HLT_IsoMu12_v1, &b_HLT_IsoMu12_v1);
   fChain->SetBranchAddress("HLT_IsoMu12_v1_Prescl", &HLT_IsoMu12_v1_Prescl, &b_HLT_IsoMu12_v1_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu15_v5", &HLT_IsoMu15_v5, &b_HLT_IsoMu15_v5);
   fChain->SetBranchAddress("HLT_IsoMu15_v5_Prescl", &HLT_IsoMu15_v5_Prescl, &b_HLT_IsoMu15_v5_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu17_v5", &HLT_IsoMu17_v5, &b_HLT_IsoMu17_v5);
   fChain->SetBranchAddress("HLT_IsoMu17_v5_Prescl", &HLT_IsoMu17_v5_Prescl, &b_HLT_IsoMu17_v5_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu24_v1", &HLT_IsoMu24_v1, &b_HLT_IsoMu24_v1);
   fChain->SetBranchAddress("HLT_IsoMu24_v1_Prescl", &HLT_IsoMu24_v1_Prescl, &b_HLT_IsoMu24_v1_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu30_v1", &HLT_IsoMu30_v1, &b_HLT_IsoMu30_v1);
   fChain->SetBranchAddress("HLT_IsoMu30_v1_Prescl", &HLT_IsoMu30_v1_Prescl, &b_HLT_IsoMu30_v1_Prescl);
   fChain->SetBranchAddress("HLT_L2DoubleMu23_NoVertex_v1", &HLT_L2DoubleMu23_NoVertex_v1, &b_HLT_L2DoubleMu23_NoVertex_v1);
   fChain->SetBranchAddress("HLT_L2DoubleMu23_NoVertex_v1_Prescl", &HLT_L2DoubleMu23_NoVertex_v1_Prescl, &b_HLT_L2DoubleMu23_NoVertex_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_v3", &HLT_DoubleMu3_v3, &b_HLT_DoubleMu3_v3);
   fChain->SetBranchAddress("HLT_DoubleMu3_v3_Prescl", &HLT_DoubleMu3_v3_Prescl, &b_HLT_DoubleMu3_v3_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu6_v1", &HLT_DoubleMu6_v1, &b_HLT_DoubleMu6_v1);
   fChain->SetBranchAddress("HLT_DoubleMu6_v1_Prescl", &HLT_DoubleMu6_v1_Prescl, &b_HLT_DoubleMu6_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu7_v1", &HLT_DoubleMu7_v1, &b_HLT_DoubleMu7_v1);
   fChain->SetBranchAddress("HLT_DoubleMu7_v1_Prescl", &HLT_DoubleMu7_v1_Prescl, &b_HLT_DoubleMu7_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v1", &HLT_DoubleMu2_Bs_v1, &b_HLT_DoubleMu2_Bs_v1);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v1_Prescl", &HLT_DoubleMu2_Bs_v1_Prescl, &b_HLT_DoubleMu2_Bs_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_Jpsi_v2", &HLT_DoubleMu3_Jpsi_v2, &b_HLT_DoubleMu3_Jpsi_v2);
   fChain->SetBranchAddress("HLT_DoubleMu3_Jpsi_v2_Prescl", &HLT_DoubleMu3_Jpsi_v2_Prescl, &b_HLT_DoubleMu3_Jpsi_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_Quarkonium_v2", &HLT_DoubleMu3_Quarkonium_v2, &b_HLT_DoubleMu3_Quarkonium_v2);
   fChain->SetBranchAddress("HLT_DoubleMu3_Quarkonium_v2_Prescl", &HLT_DoubleMu3_Quarkonium_v2_Prescl, &b_HLT_DoubleMu3_Quarkonium_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_Upsilon_v1", &HLT_DoubleMu3_Upsilon_v1, &b_HLT_DoubleMu3_Upsilon_v1);
   fChain->SetBranchAddress("HLT_DoubleMu3_Upsilon_v1_Prescl", &HLT_DoubleMu3_Upsilon_v1_Prescl, &b_HLT_DoubleMu3_Upsilon_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_LowMass_v1", &HLT_DoubleMu3_LowMass_v1, &b_HLT_DoubleMu3_LowMass_v1);
   fChain->SetBranchAddress("HLT_DoubleMu3_LowMass_v1_Prescl", &HLT_DoubleMu3_LowMass_v1_Prescl, &b_HLT_DoubleMu3_LowMass_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu4_Acoplanarity03_v1", &HLT_DoubleMu4_Acoplanarity03_v1, &b_HLT_DoubleMu4_Acoplanarity03_v1);
   fChain->SetBranchAddress("HLT_DoubleMu4_Acoplanarity03_v1_Prescl", &HLT_DoubleMu4_Acoplanarity03_v1_Prescl, &b_HLT_DoubleMu4_Acoplanarity03_v1_Prescl);
   fChain->SetBranchAddress("HLT_TripleMu5_v2", &HLT_TripleMu5_v2, &b_HLT_TripleMu5_v2);
   fChain->SetBranchAddress("HLT_TripleMu5_v2_Prescl", &HLT_TripleMu5_v2_Prescl, &b_HLT_TripleMu5_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1_Prescl", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1_Prescl, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v2", &HLT_Mu5_L2Mu2_v2, &b_HLT_Mu5_L2Mu2_v2);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v2_Prescl", &HLT_Mu5_L2Mu2_v2_Prescl, &b_HLT_Mu5_L2Mu2_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v2", &HLT_Mu5_L2Mu2_Jpsi_v2, &b_HLT_Mu5_L2Mu2_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v2_Prescl", &HLT_Mu5_L2Mu2_Jpsi_v2_Prescl, &b_HLT_Mu5_L2Mu2_Jpsi_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v5", &HLT_Mu3_Track3_Jpsi_v5, &b_HLT_Mu3_Track3_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v5_Prescl", &HLT_Mu3_Track3_Jpsi_v5_Prescl, &b_HLT_Mu3_Track3_Jpsi_v5_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v1", &HLT_Mu5_Track2_Jpsi_v1, &b_HLT_Mu5_Track2_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v1_Prescl", &HLT_Mu5_Track2_Jpsi_v1_Prescl, &b_HLT_Mu5_Track2_Jpsi_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v2", &HLT_Mu7_Track5_Jpsi_v2, &b_HLT_Mu7_Track5_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v2_Prescl", &HLT_Mu7_Track5_Jpsi_v2_Prescl, &b_HLT_Mu7_Track5_Jpsi_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v2", &HLT_Mu7_Track7_Jpsi_v2, &b_HLT_Mu7_Track7_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v2_Prescl", &HLT_Mu7_Track7_Jpsi_v2_Prescl, &b_HLT_Mu7_Track7_Jpsi_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon20_CaloIdVL_IsoL_v1", &HLT_Photon20_CaloIdVL_IsoL_v1, &b_HLT_Photon20_CaloIdVL_IsoL_v1);
   fChain->SetBranchAddress("HLT_Photon20_CaloIdVL_IsoL_v1_Prescl", &HLT_Photon20_CaloIdVL_IsoL_v1_Prescl, &b_HLT_Photon20_CaloIdVL_IsoL_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon20_R9Id_Photon18_R9Id_v2", &HLT_Photon20_R9Id_Photon18_R9Id_v2, &b_HLT_Photon20_R9Id_Photon18_R9Id_v2);
   fChain->SetBranchAddress("HLT_Photon20_R9Id_Photon18_R9Id_v2_Prescl", &HLT_Photon20_R9Id_Photon18_R9Id_v2_Prescl, &b_HLT_Photon20_R9Id_Photon18_R9Id_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2", &HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2, &b_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2);
   fChain->SetBranchAddress("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl", &HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl, &b_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_Photon18_v2", &HLT_Photon26_Photon18_v2, &b_HLT_Photon26_Photon18_v2);
   fChain->SetBranchAddress("HLT_Photon26_Photon18_v2_Prescl", &HLT_Photon26_Photon18_v2_Prescl, &b_HLT_Photon26_Photon18_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_IsoVL_Photon18_v2", &HLT_Photon26_IsoVL_Photon18_v2, &b_HLT_Photon26_IsoVL_Photon18_v2);
   fChain->SetBranchAddress("HLT_Photon26_IsoVL_Photon18_v2_Prescl", &HLT_Photon26_IsoVL_Photon18_v2_Prescl, &b_HLT_Photon26_IsoVL_Photon18_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_IsoVL_Photon18_IsoVL_v2", &HLT_Photon26_IsoVL_Photon18_IsoVL_v2, &b_HLT_Photon26_IsoVL_Photon18_IsoVL_v2);
   fChain->SetBranchAddress("HLT_Photon26_IsoVL_Photon18_IsoVL_v2_Prescl", &HLT_Photon26_IsoVL_Photon18_IsoVL_v2_Prescl, &b_HLT_Photon26_IsoVL_Photon18_IsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_CaloIdL_IsoVL_Photon18_v2", &HLT_Photon26_CaloIdL_IsoVL_Photon18_v2, &b_HLT_Photon26_CaloIdL_IsoVL_Photon18_v2);
   fChain->SetBranchAddress("HLT_Photon26_CaloIdL_IsoVL_Photon18_v2_Prescl", &HLT_Photon26_CaloIdL_IsoVL_Photon18_v2_Prescl, &b_HLT_Photon26_CaloIdL_IsoVL_Photon18_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1", &HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1, &b_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1);
   fChain->SetBranchAddress("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1_Prescl", &HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1_Prescl, &b_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2", &HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2, &b_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2);
   fChain->SetBranchAddress("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2_Prescl", &HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2_Prescl, &b_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1", &HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1, &b_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1);
   fChain->SetBranchAddress("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1_Prescl", &HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1_Prescl, &b_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon30_CaloIdVL_v2", &HLT_Photon30_CaloIdVL_v2, &b_HLT_Photon30_CaloIdVL_v2);
   fChain->SetBranchAddress("HLT_Photon30_CaloIdVL_v2_Prescl", &HLT_Photon30_CaloIdVL_v2_Prescl, &b_HLT_Photon30_CaloIdVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon30_CaloIdVL_IsoL_v2", &HLT_Photon30_CaloIdVL_IsoL_v2, &b_HLT_Photon30_CaloIdVL_IsoL_v2);
   fChain->SetBranchAddress("HLT_Photon30_CaloIdVL_IsoL_v2_Prescl", &HLT_Photon30_CaloIdVL_IsoL_v2_Prescl, &b_HLT_Photon30_CaloIdVL_IsoL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2", &HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2, &b_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2);
   fChain->SetBranchAddress("HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2_Prescl", &HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2_Prescl, &b_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1", &HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1, &b_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1);
   fChain->SetBranchAddress("HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1_Prescl", &HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1_Prescl, &b_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon50_CaloIdVL_IsoL_v1", &HLT_Photon50_CaloIdVL_IsoL_v1, &b_HLT_Photon50_CaloIdVL_IsoL_v1);
   fChain->SetBranchAddress("HLT_Photon50_CaloIdVL_IsoL_v1_Prescl", &HLT_Photon50_CaloIdVL_IsoL_v1_Prescl, &b_HLT_Photon50_CaloIdVL_IsoL_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon60_CaloIdL_HT200_v2", &HLT_Photon60_CaloIdL_HT200_v2, &b_HLT_Photon60_CaloIdL_HT200_v2);
   fChain->SetBranchAddress("HLT_Photon60_CaloIdL_HT200_v2_Prescl", &HLT_Photon60_CaloIdL_HT200_v2_Prescl, &b_HLT_Photon60_CaloIdL_HT200_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_HT200_v2", &HLT_Photon70_CaloIdL_HT200_v2, &b_HLT_Photon70_CaloIdL_HT200_v2);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_HT200_v2_Prescl", &HLT_Photon70_CaloIdL_HT200_v2_Prescl, &b_HLT_Photon70_CaloIdL_HT200_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_HT300_v2", &HLT_Photon70_CaloIdL_HT300_v2, &b_HLT_Photon70_CaloIdL_HT300_v2);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_HT300_v2_Prescl", &HLT_Photon70_CaloIdL_HT300_v2_Prescl, &b_HLT_Photon70_CaloIdL_HT300_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_MHT30_v2", &HLT_Photon70_CaloIdL_MHT30_v2, &b_HLT_Photon70_CaloIdL_MHT30_v2);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_MHT30_v2_Prescl", &HLT_Photon70_CaloIdL_MHT30_v2_Prescl, &b_HLT_Photon70_CaloIdL_MHT30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_MHT50_v2", &HLT_Photon70_CaloIdL_MHT50_v2, &b_HLT_Photon70_CaloIdL_MHT50_v2);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdL_MHT50_v2_Prescl", &HLT_Photon70_CaloIdL_MHT50_v2_Prescl, &b_HLT_Photon70_CaloIdL_MHT50_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon75_CaloIdVL_v2", &HLT_Photon75_CaloIdVL_v2, &b_HLT_Photon75_CaloIdVL_v2);
   fChain->SetBranchAddress("HLT_Photon75_CaloIdVL_v2_Prescl", &HLT_Photon75_CaloIdVL_v2_Prescl, &b_HLT_Photon75_CaloIdVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon75_CaloIdVL_IsoL_v2", &HLT_Photon75_CaloIdVL_IsoL_v2, &b_HLT_Photon75_CaloIdVL_IsoL_v2);
   fChain->SetBranchAddress("HLT_Photon75_CaloIdVL_IsoL_v2_Prescl", &HLT_Photon75_CaloIdVL_IsoL_v2_Prescl, &b_HLT_Photon75_CaloIdVL_IsoL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Photon125_NoSpikeFilter_v2", &HLT_Photon125_NoSpikeFilter_v2, &b_HLT_Photon125_NoSpikeFilter_v2);
   fChain->SetBranchAddress("HLT_Photon125_NoSpikeFilter_v2_Prescl", &HLT_Photon125_NoSpikeFilter_v2_Prescl, &b_HLT_Photon125_NoSpikeFilter_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoublePhoton33_v2", &HLT_DoublePhoton33_v2, &b_HLT_DoublePhoton33_v2);
   fChain->SetBranchAddress("HLT_DoublePhoton33_v2_Prescl", &HLT_DoublePhoton33_v2_Prescl, &b_HLT_DoublePhoton33_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoublePhoton5_IsoVL_CEP_v1", &HLT_DoublePhoton5_IsoVL_CEP_v1, &b_HLT_DoublePhoton5_IsoVL_CEP_v1);
   fChain->SetBranchAddress("HLT_DoublePhoton5_IsoVL_CEP_v1_Prescl", &HLT_DoublePhoton5_IsoVL_CEP_v1_Prescl, &b_HLT_DoublePhoton5_IsoVL_CEP_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleEG5_v1", &HLT_L1SingleEG5_v1, &b_HLT_L1SingleEG5_v1);
   fChain->SetBranchAddress("HLT_L1SingleEG5_v1_Prescl", &HLT_L1SingleEG5_v1_Prescl, &b_HLT_L1SingleEG5_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleEG12_v1", &HLT_L1SingleEG12_v1, &b_HLT_L1SingleEG12_v1);
   fChain->SetBranchAddress("HLT_L1SingleEG12_v1_Prescl", &HLT_L1SingleEG12_v1_Prescl, &b_HLT_L1SingleEG12_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele8_v2", &HLT_Ele8_v2, &b_HLT_Ele8_v2);
   fChain->SetBranchAddress("HLT_Ele8_v2_Prescl", &HLT_Ele8_v2_Prescl, &b_HLT_Ele8_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_CaloIsoVL_v2", &HLT_Ele8_CaloIdL_CaloIsoVL_v2, &b_HLT_Ele8_CaloIdL_CaloIsoVL_v2);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl", &HLT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl, &b_HLT_Ele8_CaloIdL_CaloIsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrkIdVL_v2", &HLT_Ele8_CaloIdL_TrkIdVL_v2, &b_HLT_Ele8_CaloIdL_TrkIdVL_v2);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrkIdVL_v2_Prescl", &HLT_Ele8_CaloIdL_TrkIdVL_v2_Prescl, &b_HLT_Ele8_CaloIdL_TrkIdVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2, &b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl, &b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_CaloIsoVL_v2", &HLT_Ele17_CaloIdL_CaloIsoVL_v2, &b_HLT_Ele17_CaloIdL_CaloIsoVL_v2);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_CaloIsoVL_v2_Prescl", &HLT_Ele17_CaloIdL_CaloIsoVL_v2_Prescl, &b_HLT_Ele17_CaloIdL_CaloIsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2", &HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2, &b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2_Prescl", &HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2_Prescl, &b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2", &HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2, &b_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2_Prescl", &HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2_Prescl, &b_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2", &HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2, &b_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2_Prescl", &HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2_Prescl, &b_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2", &HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2, &b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2_Prescl", &HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2_Prescl, &b_HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2", &HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2, &b_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2);
   fChain->SetBranchAddress("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl", &HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl, &b_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1", &HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1, &b_HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1);
   fChain->SetBranchAddress("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1_Prescl", &HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1_Prescl, &b_HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2", &HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2, &b_HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2);
   fChain->SetBranchAddress("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2_Prescl", &HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2_Prescl, &b_HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele45_CaloIdVT_TrkIdT_v2", &HLT_Ele45_CaloIdVT_TrkIdT_v2, &b_HLT_Ele45_CaloIdVT_TrkIdT_v2);
   fChain->SetBranchAddress("HLT_Ele45_CaloIdVT_TrkIdT_v2_Prescl", &HLT_Ele45_CaloIdVT_TrkIdT_v2_Prescl, &b_HLT_Ele45_CaloIdVT_TrkIdT_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele90_NoSpikeFilter_v2", &HLT_Ele90_NoSpikeFilter_v2, &b_HLT_Ele90_NoSpikeFilter_v2);
   fChain->SetBranchAddress("HLT_Ele90_NoSpikeFilter_v2_Prescl", &HLT_Ele90_NoSpikeFilter_v2_Prescl, &b_HLT_Ele90_NoSpikeFilter_v2_Prescl);
   fChain->SetBranchAddress("HLT_IsoPFTau35_Trk20_MET45_v2", &HLT_IsoPFTau35_Trk20_MET45_v2, &b_HLT_IsoPFTau35_Trk20_MET45_v2);
   fChain->SetBranchAddress("HLT_IsoPFTau35_Trk20_MET45_v2_Prescl", &HLT_IsoPFTau35_Trk20_MET45_v2_Prescl, &b_HLT_IsoPFTau35_Trk20_MET45_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleIsoPFTau20_Trk5_v2", &HLT_DoubleIsoPFTau20_Trk5_v2, &b_HLT_DoubleIsoPFTau20_Trk5_v2);
   fChain->SetBranchAddress("HLT_DoubleIsoPFTau20_Trk5_v2_Prescl", &HLT_DoubleIsoPFTau20_Trk5_v2_Prescl, &b_HLT_DoubleIsoPFTau20_Trk5_v2_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20_Mu5_v2", &HLT_BTagMu_DiJet20_Mu5_v2, &b_HLT_BTagMu_DiJet20_Mu5_v2);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20_Mu5_v2_Prescl", &HLT_BTagMu_DiJet20_Mu5_v2_Prescl, &b_HLT_BTagMu_DiJet20_Mu5_v2_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet60_Mu7_v2", &HLT_BTagMu_DiJet60_Mu7_v2, &b_HLT_BTagMu_DiJet60_Mu7_v2);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet60_Mu7_v2_Prescl", &HLT_BTagMu_DiJet60_Mu7_v2_Prescl, &b_HLT_BTagMu_DiJet60_Mu7_v2_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet80_Mu9_v2", &HLT_BTagMu_DiJet80_Mu9_v2, &b_HLT_BTagMu_DiJet80_Mu9_v2);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet80_Mu9_v2_Prescl", &HLT_BTagMu_DiJet80_Mu9_v2_Prescl, &b_HLT_BTagMu_DiJet80_Mu9_v2_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet100_Mu9_v2", &HLT_BTagMu_DiJet100_Mu9_v2, &b_HLT_BTagMu_DiJet100_Mu9_v2);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet100_Mu9_v2_Prescl", &HLT_BTagMu_DiJet100_Mu9_v2_Prescl, &b_HLT_BTagMu_DiJet100_Mu9_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3", &HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3, &b_HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3);
   fChain->SetBranchAddress("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3_Prescl", &HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3_Prescl, &b_HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3", &HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3, &b_HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3);
   fChain->SetBranchAddress("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3_Prescl", &HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3_Prescl, &b_HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3", &HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3, &b_HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3);
   fChain->SetBranchAddress("HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3_Prescl", &HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3_Prescl, &b_HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_DoubleEle8_v3", &HLT_Mu5_DoubleEle8_v3, &b_HLT_Mu5_DoubleEle8_v3);
   fChain->SetBranchAddress("HLT_Mu5_DoubleEle8_v3_Prescl", &HLT_Mu5_DoubleEle8_v3_Prescl, &b_HLT_Mu5_DoubleEle8_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_HT200_v4", &HLT_Mu5_HT200_v4, &b_HLT_Mu5_HT200_v4);
   fChain->SetBranchAddress("HLT_Mu5_HT200_v4_Prescl", &HLT_Mu5_HT200_v4_Prescl, &b_HLT_Mu5_HT200_v4_Prescl);
   fChain->SetBranchAddress("HLT_Mu8_HT200_v3", &HLT_Mu8_HT200_v3, &b_HLT_Mu8_HT200_v3);
   fChain->SetBranchAddress("HLT_Mu8_HT200_v3_Prescl", &HLT_Mu8_HT200_v3_Prescl, &b_HLT_Mu8_HT200_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu8_Ele17_CaloIdL_v2", &HLT_Mu8_Ele17_CaloIdL_v2, &b_HLT_Mu8_Ele17_CaloIdL_v2);
   fChain->SetBranchAddress("HLT_Mu8_Ele17_CaloIdL_v2_Prescl", &HLT_Mu8_Ele17_CaloIdL_v2_Prescl, &b_HLT_Mu8_Ele17_CaloIdL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2", &HLT_Mu8_Photon20_CaloIdVT_IsoT_v2, &b_HLT_Mu8_Photon20_CaloIdVT_IsoT_v2);
   fChain->SetBranchAddress("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2_Prescl", &HLT_Mu8_Photon20_CaloIdVT_IsoT_v2_Prescl, &b_HLT_Mu8_Photon20_CaloIdVT_IsoT_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu8_Jet40_v3", &HLT_Mu8_Jet40_v3, &b_HLT_Mu8_Jet40_v3);
   fChain->SetBranchAddress("HLT_Mu8_Jet40_v3_Prescl", &HLT_Mu8_Jet40_v3_Prescl, &b_HLT_Mu8_Jet40_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu10_Ele10_CaloIdL_v3", &HLT_Mu10_Ele10_CaloIdL_v3, &b_HLT_Mu10_Ele10_CaloIdL_v3);
   fChain->SetBranchAddress("HLT_Mu10_Ele10_CaloIdL_v3_Prescl", &HLT_Mu10_Ele10_CaloIdL_v3_Prescl, &b_HLT_Mu10_Ele10_CaloIdL_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu15_Photon20_CaloIdL_v3", &HLT_Mu15_Photon20_CaloIdL_v3, &b_HLT_Mu15_Photon20_CaloIdL_v3);
   fChain->SetBranchAddress("HLT_Mu15_Photon20_CaloIdL_v3_Prescl", &HLT_Mu15_Photon20_CaloIdL_v3_Prescl, &b_HLT_Mu15_Photon20_CaloIdL_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu15_DoublePhoton15_CaloIdL_v3", &HLT_Mu15_DoublePhoton15_CaloIdL_v3, &b_HLT_Mu15_DoublePhoton15_CaloIdL_v3);
   fChain->SetBranchAddress("HLT_Mu15_DoublePhoton15_CaloIdL_v3_Prescl", &HLT_Mu15_DoublePhoton15_CaloIdL_v3_Prescl, &b_HLT_Mu15_DoublePhoton15_CaloIdL_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu15_LooseIsoPFTau20_v2", &HLT_Mu15_LooseIsoPFTau20_v2, &b_HLT_Mu15_LooseIsoPFTau20_v2);
   fChain->SetBranchAddress("HLT_Mu15_LooseIsoPFTau20_v2_Prescl", &HLT_Mu15_LooseIsoPFTau20_v2_Prescl, &b_HLT_Mu15_LooseIsoPFTau20_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu17_CentralJet30_v2", &HLT_Mu17_CentralJet30_v2, &b_HLT_Mu17_CentralJet30_v2);
   fChain->SetBranchAddress("HLT_Mu17_CentralJet30_v2_Prescl", &HLT_Mu17_CentralJet30_v2_Prescl, &b_HLT_Mu17_CentralJet30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu17_DiCentralJet30_v2", &HLT_Mu17_DiCentralJet30_v2, &b_HLT_Mu17_DiCentralJet30_v2);
   fChain->SetBranchAddress("HLT_Mu17_DiCentralJet30_v2_Prescl", &HLT_Mu17_DiCentralJet30_v2_Prescl, &b_HLT_Mu17_DiCentralJet30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu17_TriCentralJet30_v2", &HLT_Mu17_TriCentralJet30_v2, &b_HLT_Mu17_TriCentralJet30_v2);
   fChain->SetBranchAddress("HLT_Mu17_TriCentralJet30_v2_Prescl", &HLT_Mu17_TriCentralJet30_v2_Prescl, &b_HLT_Mu17_TriCentralJet30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu17_Ele8_CaloIdL_v2", &HLT_Mu17_Ele8_CaloIdL_v2, &b_HLT_Mu17_Ele8_CaloIdL_v2);
   fChain->SetBranchAddress("HLT_Mu17_Ele8_CaloIdL_v2_Prescl", &HLT_Mu17_Ele8_CaloIdL_v2_Prescl, &b_HLT_Mu17_Ele8_CaloIdL_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu17_CentralJet40_BTagIP_v2", &HLT_Mu17_CentralJet40_BTagIP_v2, &b_HLT_Mu17_CentralJet40_BTagIP_v2);
   fChain->SetBranchAddress("HLT_Mu17_CentralJet40_BTagIP_v2_Prescl", &HLT_Mu17_CentralJet40_BTagIP_v2_Prescl, &b_HLT_Mu17_CentralJet40_BTagIP_v2_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu12_LooseIsoPFTau10_v2", &HLT_IsoMu12_LooseIsoPFTau10_v2, &b_HLT_IsoMu12_LooseIsoPFTau10_v2);
   fChain->SetBranchAddress("HLT_IsoMu12_LooseIsoPFTau10_v2_Prescl", &HLT_IsoMu12_LooseIsoPFTau10_v2_Prescl, &b_HLT_IsoMu12_LooseIsoPFTau10_v2_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu17_CentralJet40_BTagIP_v2", &HLT_IsoMu17_CentralJet40_BTagIP_v2, &b_HLT_IsoMu17_CentralJet40_BTagIP_v2);
   fChain->SetBranchAddress("HLT_IsoMu17_CentralJet40_BTagIP_v2_Prescl", &HLT_IsoMu17_CentralJet40_BTagIP_v2_Prescl, &b_HLT_IsoMu17_CentralJet40_BTagIP_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_HT160_v3", &HLT_DoubleMu3_HT160_v3, &b_HLT_DoubleMu3_HT160_v3);
   fChain->SetBranchAddress("HLT_DoubleMu3_HT160_v3_Prescl", &HLT_DoubleMu3_HT160_v3_Prescl, &b_HLT_DoubleMu3_HT160_v3_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_HT200_v3", &HLT_DoubleMu3_HT200_v3, &b_HLT_DoubleMu3_HT200_v3);
   fChain->SetBranchAddress("HLT_DoubleMu3_HT200_v3_Prescl", &HLT_DoubleMu3_HT200_v3_Prescl, &b_HLT_DoubleMu3_HT200_v3_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu5_Ele8_v3", &HLT_DoubleMu5_Ele8_v3, &b_HLT_DoubleMu5_Ele8_v3);
   fChain->SetBranchAddress("HLT_DoubleMu5_Ele8_v3_Prescl", &HLT_DoubleMu5_Ele8_v3_Prescl, &b_HLT_DoubleMu5_Ele8_v3_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3", &HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3, &b_HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3);
   fChain->SetBranchAddress("HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3_Prescl", &HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3_Prescl, &b_HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v3_Prescl);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2", &HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2, &b_HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2_Prescl", &HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2_Prescl, &b_HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3", &HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3, &b_HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3_Prescl", &HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3_Prescl, &b_HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT200_v3_Prescl);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3", &HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3, &b_HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3_Prescl", &HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3_Prescl, &b_HLT_Ele10_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v3_Prescl);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2", &HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2, &b_HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2_Prescl", &HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2_Prescl, &b_HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2, &b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2_Prescl", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2_Prescl, &b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2, &b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2);
   fChain->SetBranchAddress("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2_Prescl", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2_Prescl, &b_HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2", &HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2_Prescl", &HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2_Prescl, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2", &HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2_Prescl", &HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2_Prescl, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralDiJet30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2", &HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2_Prescl", &HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2_Prescl, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2_Prescl);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2", &HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2);
   fChain->SetBranchAddress("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2_Prescl", &HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2_Prescl, &b_HLT_Ele25_CaloIdVT_TrkIdT_CentralJet40_BTagIP_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3", &HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3, &b_HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3_Prescl", &HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3_Prescl, &b_HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3_Prescl);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3", &HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3, &b_HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3_Prescl", &HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3_Prescl, &b_HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3_Prescl);
   fChain->SetBranchAddress("HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2", &HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2, &b_HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2);
   fChain->SetBranchAddress("HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2_Prescl", &HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2_Prescl, &b_HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2_Prescl);
   fChain->SetBranchAddress("HLT_TripleEle10_CaloIdL_TrkIdVL_v2", &HLT_TripleEle10_CaloIdL_TrkIdVL_v2, &b_HLT_TripleEle10_CaloIdL_TrkIdVL_v2);
   fChain->SetBranchAddress("HLT_TripleEle10_CaloIdL_TrkIdVL_v2_Prescl", &HLT_TripleEle10_CaloIdL_TrkIdVL_v2_Prescl, &b_HLT_TripleEle10_CaloIdL_TrkIdVL_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity80_v2", &HLT_PixelTracks_Multiplicity80_v2, &b_HLT_PixelTracks_Multiplicity80_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity80_v2_Prescl", &HLT_PixelTracks_Multiplicity80_v2_Prescl, &b_HLT_PixelTracks_Multiplicity80_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity100_v2", &HLT_PixelTracks_Multiplicity100_v2, &b_HLT_PixelTracks_Multiplicity100_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity100_v2_Prescl", &HLT_PixelTracks_Multiplicity100_v2_Prescl, &b_HLT_PixelTracks_Multiplicity100_v2_Prescl);
   fChain->SetBranchAddress("HLT_BeamGas_HF_v2", &HLT_BeamGas_HF_v2, &b_HLT_BeamGas_HF_v2);
   fChain->SetBranchAddress("HLT_BeamGas_HF_v2_Prescl", &HLT_BeamGas_HF_v2_Prescl, &b_HLT_BeamGas_HF_v2_Prescl);
   fChain->SetBranchAddress("HLT_BeamGas_BSC_v2", &HLT_BeamGas_BSC_v2, &b_HLT_BeamGas_BSC_v2);
   fChain->SetBranchAddress("HLT_BeamGas_BSC_v2_Prescl", &HLT_BeamGas_BSC_v2_Prescl, &b_HLT_BeamGas_BSC_v2_Prescl);
   fChain->SetBranchAddress("HLT_BeamHalo_v2", &HLT_BeamHalo_v2, &b_HLT_BeamHalo_v2);
   fChain->SetBranchAddress("HLT_BeamHalo_v2_Prescl", &HLT_BeamHalo_v2_Prescl, &b_HLT_BeamHalo_v2_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_minBias_threshold1_v1", &HLT_L1Tech_BSC_minBias_threshold1_v1, &b_HLT_L1Tech_BSC_minBias_threshold1_v1);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_minBias_threshold1_v1_Prescl", &HLT_L1Tech_BSC_minBias_threshold1_v1_Prescl, &b_HLT_L1Tech_BSC_minBias_threshold1_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_halo_v1", &HLT_L1Tech_BSC_halo_v1, &b_HLT_L1Tech_BSC_halo_v1);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_halo_v1_Prescl", &HLT_L1Tech_BSC_halo_v1_Prescl, &b_HLT_L1Tech_BSC_halo_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_CASTOR_HaloMuon_v1", &HLT_L1Tech_CASTOR_HaloMuon_v1, &b_HLT_L1Tech_CASTOR_HaloMuon_v1);
   fChain->SetBranchAddress("HLT_L1Tech_CASTOR_HaloMuon_v1_Prescl", &HLT_L1Tech_CASTOR_HaloMuon_v1_Prescl, &b_HLT_L1Tech_CASTOR_HaloMuon_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1_PreCollisions_v1", &HLT_L1_PreCollisions_v1, &b_HLT_L1_PreCollisions_v1);
   fChain->SetBranchAddress("HLT_L1_PreCollisions_v1_Prescl", &HLT_L1_PreCollisions_v1_Prescl, &b_HLT_L1_PreCollisions_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1_Interbunch_BSC_v1", &HLT_L1_Interbunch_BSC_v1, &b_HLT_L1_Interbunch_BSC_v1);
   fChain->SetBranchAddress("HLT_L1_Interbunch_BSC_v1_Prescl", &HLT_L1_Interbunch_BSC_v1_Prescl, &b_HLT_L1_Interbunch_BSC_v1_Prescl);
   fChain->SetBranchAddress("HLT_IsoTrackHE_v3", &HLT_IsoTrackHE_v3, &b_HLT_IsoTrackHE_v3);
   fChain->SetBranchAddress("HLT_IsoTrackHE_v3_Prescl", &HLT_IsoTrackHE_v3_Prescl, &b_HLT_IsoTrackHE_v3_Prescl);
   fChain->SetBranchAddress("HLT_IsoTrackHB_v2", &HLT_IsoTrackHB_v2, &b_HLT_IsoTrackHB_v2);
   fChain->SetBranchAddress("HLT_IsoTrackHB_v2_Prescl", &HLT_IsoTrackHB_v2_Prescl, &b_HLT_IsoTrackHB_v2_Prescl);
   fChain->SetBranchAddress("HLT_HcalPhiSym_v3", &HLT_HcalPhiSym_v3, &b_HLT_HcalPhiSym_v3);
   fChain->SetBranchAddress("HLT_HcalPhiSym_v3_Prescl", &HLT_HcalPhiSym_v3_Prescl, &b_HLT_HcalPhiSym_v3_Prescl);
   fChain->SetBranchAddress("HLT_HcalNZS_v3", &HLT_HcalNZS_v3, &b_HLT_HcalNZS_v3);
   fChain->SetBranchAddress("HLT_HcalNZS_v3_Prescl", &HLT_HcalNZS_v3_Prescl, &b_HLT_HcalNZS_v3_Prescl);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise_v2", &HLT_GlobalRunHPDNoise_v2, &b_HLT_GlobalRunHPDNoise_v2);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise_v2_Prescl", &HLT_GlobalRunHPDNoise_v2_Prescl, &b_HLT_GlobalRunHPDNoise_v2_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_HBHEHO_totalOR_v1", &HLT_L1Tech_HBHEHO_totalOR_v1, &b_HLT_L1Tech_HBHEHO_totalOR_v1);
   fChain->SetBranchAddress("HLT_L1Tech_HBHEHO_totalOR_v1_Prescl", &HLT_L1Tech_HBHEHO_totalOR_v1_Prescl, &b_HLT_L1Tech_HBHEHO_totalOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_v1", &HLT_ZeroBias_v1, &b_HLT_ZeroBias_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_v1_Prescl", &HLT_ZeroBias_v1_Prescl, &b_HLT_ZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_Physics_v1", &HLT_Physics_v1, &b_HLT_Physics_v1);
   fChain->SetBranchAddress("HLT_Physics_v1_Prescl", &HLT_Physics_v1_Prescl, &b_HLT_Physics_v1_Prescl);
   fChain->SetBranchAddress("HLT_Physics_NanoDST_v1", &HLT_Physics_NanoDST_v1, &b_HLT_Physics_NanoDST_v1);
   fChain->SetBranchAddress("HLT_Physics_NanoDST_v1_Prescl", &HLT_Physics_NanoDST_v1_Prescl, &b_HLT_Physics_NanoDST_v1_Prescl);
   fChain->SetBranchAddress("HLT_Calibration_v1", &HLT_Calibration_v1, &b_HLT_Calibration_v1);
   fChain->SetBranchAddress("HLT_Calibration_v1_Prescl", &HLT_Calibration_v1_Prescl, &b_HLT_Calibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_EcalCalibration_v1", &HLT_EcalCalibration_v1, &b_HLT_EcalCalibration_v1);
   fChain->SetBranchAddress("HLT_EcalCalibration_v1_Prescl", &HLT_EcalCalibration_v1_Prescl, &b_HLT_EcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_HcalCalibration_v1", &HLT_HcalCalibration_v1, &b_HLT_HcalCalibration_v1);
   fChain->SetBranchAddress("HLT_HcalCalibration_v1_Prescl", &HLT_HcalCalibration_v1_Prescl, &b_HLT_HcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_TrackerCalibration_v1", &HLT_TrackerCalibration_v1, &b_HLT_TrackerCalibration_v1);
   fChain->SetBranchAddress("HLT_TrackerCalibration_v1_Prescl", &HLT_TrackerCalibration_v1_Prescl, &b_HLT_TrackerCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_Random_v1", &HLT_Random_v1, &b_HLT_Random_v1);
   fChain->SetBranchAddress("HLT_Random_v1_Prescl", &HLT_Random_v1_Prescl, &b_HLT_Random_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_AntiBPTX_v1", &HLT_L1SingleMuOpen_AntiBPTX_v1, &b_HLT_L1SingleMuOpen_AntiBPTX_v1);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_AntiBPTX_v1_Prescl", &HLT_L1SingleMuOpen_AntiBPTX_v1_Prescl, &b_HLT_L1SingleMuOpen_AntiBPTX_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1TrackerCosmics_v2", &HLT_L1TrackerCosmics_v2, &b_HLT_L1TrackerCosmics_v2);
   fChain->SetBranchAddress("HLT_L1TrackerCosmics_v2_Prescl", &HLT_L1TrackerCosmics_v2_Prescl, &b_HLT_L1TrackerCosmics_v2_Prescl);
   fChain->SetBranchAddress("HLT_RegionalCosmicTracking_v1", &HLT_RegionalCosmicTracking_v1, &b_HLT_RegionalCosmicTracking_v1);
   fChain->SetBranchAddress("HLT_RegionalCosmicTracking_v1_Prescl", &HLT_RegionalCosmicTracking_v1_Prescl, &b_HLT_RegionalCosmicTracking_v1_Prescl);
   fChain->SetBranchAddress("HLT_L3MuonsCosmicTracking_v1", &HLT_L3MuonsCosmicTracking_v1, &b_HLT_L3MuonsCosmicTracking_v1);
   fChain->SetBranchAddress("HLT_L3MuonsCosmicTracking_v1_Prescl", &HLT_L3MuonsCosmicTracking_v1_Prescl, &b_HLT_L3MuonsCosmicTracking_v1_Prescl);
   fChain->SetBranchAddress("HLT_LogMonitor_v1", &HLT_LogMonitor_v1, &b_HLT_LogMonitor_v1);
   fChain->SetBranchAddress("HLT_LogMonitor_v1_Prescl", &HLT_LogMonitor_v1_Prescl, &b_HLT_LogMonitor_v1_Prescl);
   fChain->SetBranchAddress("HLT_DTErrors_v1", &HLT_DTErrors_v1, &b_HLT_DTErrors_v1);
   fChain->SetBranchAddress("HLT_DTErrors_v1_Prescl", &HLT_DTErrors_v1_Prescl, &b_HLT_DTErrors_v1_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPi0_v4", &AlCa_EcalPi0_v4, &b_AlCa_EcalPi0_v4);
   fChain->SetBranchAddress("AlCa_EcalPi0_v4_Prescl", &AlCa_EcalPi0_v4_Prescl, &b_AlCa_EcalPi0_v4_Prescl);
   fChain->SetBranchAddress("AlCa_EcalEta_v3", &AlCa_EcalEta_v3, &b_AlCa_EcalEta_v3);
   fChain->SetBranchAddress("AlCa_EcalEta_v3_Prescl", &AlCa_EcalEta_v3_Prescl, &b_AlCa_EcalEta_v3_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v2", &AlCa_EcalPhiSym_v2, &b_AlCa_EcalPhiSym_v2);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v2_Prescl", &AlCa_EcalPhiSym_v2_Prescl, &b_AlCa_EcalPhiSym_v2_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_v2", &AlCa_RPCMuonNoTriggers_v2, &b_AlCa_RPCMuonNoTriggers_v2);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_v2_Prescl", &AlCa_RPCMuonNoTriggers_v2_Prescl, &b_AlCa_RPCMuonNoTriggers_v2_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_v2", &AlCa_RPCMuonNoHits_v2, &b_AlCa_RPCMuonNoHits_v2);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_v2_Prescl", &AlCa_RPCMuonNoHits_v2_Prescl, &b_AlCa_RPCMuonNoHits_v2_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_v2", &AlCa_RPCMuonNormalisation_v2, &b_AlCa_RPCMuonNormalisation_v2);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_v2_Prescl", &AlCa_RPCMuonNormalisation_v2_Prescl, &b_AlCa_RPCMuonNormalisation_v2_Prescl);
   fChain->SetBranchAddress("DQM_FEDIntegrity_v3", &DQM_FEDIntegrity_v3, &b_DQM_FEDIntegrity_v3);
   fChain->SetBranchAddress("DQM_FEDIntegrity_v3_Prescl", &DQM_FEDIntegrity_v3_Prescl, &b_DQM_FEDIntegrity_v3_Prescl);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFinalPath_Prescl", &HLTriggerFinalPath_Prescl, &b_HLTriggerFinalPath_Prescl);
   fChain->SetBranchAddress("L1_BeamGas_Bsc", &L1_BeamGas_Bsc, &b_L1_BeamGas_Bsc);
   fChain->SetBranchAddress("L1_BeamGas_Bsc_Prescl", &L1_BeamGas_Bsc_Prescl, &b_L1_BeamGas_Bsc_Prescl);
   fChain->SetBranchAddress("L1_BeamGas_Bsc_5bx", &L1_BeamGas_Bsc_5bx, &b_L1_BeamGas_Bsc_5bx);
   fChain->SetBranchAddress("L1_BeamGas_Hf", &L1_BeamGas_Hf, &b_L1_BeamGas_Hf);
   fChain->SetBranchAddress("L1_BeamGas_Hf_Prescl", &L1_BeamGas_Hf_Prescl, &b_L1_BeamGas_Hf_Prescl);
   fChain->SetBranchAddress("L1_BeamGas_Hf_5bx", &L1_BeamGas_Hf_5bx, &b_L1_BeamGas_Hf_5bx);
   fChain->SetBranchAddress("L1_BeamHalo", &L1_BeamHalo, &b_L1_BeamHalo);
   fChain->SetBranchAddress("L1_BeamHalo_Prescl", &L1_BeamHalo_Prescl, &b_L1_BeamHalo_Prescl);
   fChain->SetBranchAddress("L1_BeamHalo_5bx", &L1_BeamHalo_5bx, &b_L1_BeamHalo_5bx);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus", &L1_BptxMinus_NotBptxPlus, &b_L1_BptxMinus_NotBptxPlus);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus_Prescl", &L1_BptxMinus_NotBptxPlus_Prescl, &b_L1_BptxMinus_NotBptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus_5bx", &L1_BptxMinus_NotBptxPlus_5bx, &b_L1_BptxMinus_NotBptxPlus_5bx);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus", &L1_BptxPlus_NotBptxMinus, &b_L1_BptxPlus_NotBptxMinus);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus_Prescl", &L1_BptxPlus_NotBptxMinus_Prescl, &b_L1_BptxPlus_NotBptxMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus_5bx", &L1_BptxPlus_NotBptxMinus_5bx, &b_L1_BptxPlus_NotBptxMinus_5bx);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus", &L1_Bsc2Minus_BptxMinus, &b_L1_Bsc2Minus_BptxMinus);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus_Prescl", &L1_Bsc2Minus_BptxMinus_Prescl, &b_L1_Bsc2Minus_BptxMinus_Prescl);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus_5bx", &L1_Bsc2Minus_BptxMinus_5bx, &b_L1_Bsc2Minus_BptxMinus_5bx);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus", &L1_Bsc2Plus_BptxPlus, &b_L1_Bsc2Plus_BptxPlus);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus_Prescl", &L1_Bsc2Plus_BptxPlus_Prescl, &b_L1_Bsc2Plus_BptxPlus_Prescl);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus_5bx", &L1_Bsc2Plus_BptxPlus_5bx, &b_L1_Bsc2Plus_BptxPlus_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusANDMinus", &L1_BscMinBiasOR_BptxPlusANDMinus, &b_L1_BscMinBiasOR_BptxPlusANDMinus);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusANDMinus_Prescl", &L1_BscMinBiasOR_BptxPlusANDMinus_Prescl, &b_L1_BscMinBiasOR_BptxPlusANDMinus_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusANDMinus_5bx", &L1_BscMinBiasOR_BptxPlusANDMinus_5bx, &b_L1_BscMinBiasOR_BptxPlusANDMinus_5bx);
   fChain->SetBranchAddress("L1_DoubleEG10", &L1_DoubleEG10, &b_L1_DoubleEG10);
   fChain->SetBranchAddress("L1_DoubleEG10_Prescl", &L1_DoubleEG10_Prescl, &b_L1_DoubleEG10_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG10_5bx", &L1_DoubleEG10_5bx, &b_L1_DoubleEG10_5bx);
   fChain->SetBranchAddress("L1_DoubleEG2_FwdVeto", &L1_DoubleEG2_FwdVeto, &b_L1_DoubleEG2_FwdVeto);
   fChain->SetBranchAddress("L1_DoubleEG2_FwdVeto_Prescl", &L1_DoubleEG2_FwdVeto_Prescl, &b_L1_DoubleEG2_FwdVeto_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG2_FwdVeto_5bx", &L1_DoubleEG2_FwdVeto_5bx, &b_L1_DoubleEG2_FwdVeto_5bx);
   fChain->SetBranchAddress("L1_DoubleEG3", &L1_DoubleEG3, &b_L1_DoubleEG3);
   fChain->SetBranchAddress("L1_DoubleEG3_Prescl", &L1_DoubleEG3_Prescl, &b_L1_DoubleEG3_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG3_5bx", &L1_DoubleEG3_5bx, &b_L1_DoubleEG3_5bx);
   fChain->SetBranchAddress("L1_DoubleEG5", &L1_DoubleEG5, &b_L1_DoubleEG5);
   fChain->SetBranchAddress("L1_DoubleEG5_Prescl", &L1_DoubleEG5_Prescl, &b_L1_DoubleEG5_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5_5bx", &L1_DoubleEG5_5bx, &b_L1_DoubleEG5_5bx);
   fChain->SetBranchAddress("L1_DoubleEG5_HTT50", &L1_DoubleEG5_HTT50, &b_L1_DoubleEG5_HTT50);
   fChain->SetBranchAddress("L1_DoubleEG5_HTT50_Prescl", &L1_DoubleEG5_HTT50_Prescl, &b_L1_DoubleEG5_HTT50_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5_HTT50_5bx", &L1_DoubleEG5_HTT50_5bx, &b_L1_DoubleEG5_HTT50_5bx);
   fChain->SetBranchAddress("L1_DoubleEG5_HTT75", &L1_DoubleEG5_HTT75, &b_L1_DoubleEG5_HTT75);
   fChain->SetBranchAddress("L1_DoubleEG5_HTT75_Prescl", &L1_DoubleEG5_HTT75_Prescl, &b_L1_DoubleEG5_HTT75_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5_HTT75_5bx", &L1_DoubleEG5_HTT75_5bx, &b_L1_DoubleEG5_HTT75_5bx);
   fChain->SetBranchAddress("L1_DoubleEG8", &L1_DoubleEG8, &b_L1_DoubleEG8);
   fChain->SetBranchAddress("L1_DoubleEG8_Prescl", &L1_DoubleEG8_Prescl, &b_L1_DoubleEG8_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG8_5bx", &L1_DoubleEG8_5bx, &b_L1_DoubleEG8_5bx);
   fChain->SetBranchAddress("L1_DoubleEG_12_5", &L1_DoubleEG_12_5, &b_L1_DoubleEG_12_5);
   fChain->SetBranchAddress("L1_DoubleEG_12_5_Prescl", &L1_DoubleEG_12_5_Prescl, &b_L1_DoubleEG_12_5_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG_12_5_5bx", &L1_DoubleEG_12_5_5bx, &b_L1_DoubleEG_12_5_5bx);
   fChain->SetBranchAddress("L1_DoubleEG_12_5_Eta1p39", &L1_DoubleEG_12_5_Eta1p39, &b_L1_DoubleEG_12_5_Eta1p39);
   fChain->SetBranchAddress("L1_DoubleEG_12_5_Eta1p39_Prescl", &L1_DoubleEG_12_5_Eta1p39_Prescl, &b_L1_DoubleEG_12_5_Eta1p39_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG_12_5_Eta1p39_5bx", &L1_DoubleEG_12_5_Eta1p39_5bx, &b_L1_DoubleEG_12_5_Eta1p39_5bx);
   fChain->SetBranchAddress("L1_DoubleForJet32_EtaOpp", &L1_DoubleForJet32_EtaOpp, &b_L1_DoubleForJet32_EtaOpp);
   fChain->SetBranchAddress("L1_DoubleForJet32_EtaOpp_Prescl", &L1_DoubleForJet32_EtaOpp_Prescl, &b_L1_DoubleForJet32_EtaOpp_Prescl);
   fChain->SetBranchAddress("L1_DoubleForJet32_EtaOpp_5bx", &L1_DoubleForJet32_EtaOpp_5bx, &b_L1_DoubleForJet32_EtaOpp_5bx);
   fChain->SetBranchAddress("L1_DoubleForJet44_EtaOpp", &L1_DoubleForJet44_EtaOpp, &b_L1_DoubleForJet44_EtaOpp);
   fChain->SetBranchAddress("L1_DoubleForJet44_EtaOpp_Prescl", &L1_DoubleForJet44_EtaOpp_Prescl, &b_L1_DoubleForJet44_EtaOpp_Prescl);
   fChain->SetBranchAddress("L1_DoubleForJet44_EtaOpp_5bx", &L1_DoubleForJet44_EtaOpp_5bx, &b_L1_DoubleForJet44_EtaOpp_5bx);
   fChain->SetBranchAddress("L1_DoubleIsoEG10", &L1_DoubleIsoEG10, &b_L1_DoubleIsoEG10);
   fChain->SetBranchAddress("L1_DoubleIsoEG10_Prescl", &L1_DoubleIsoEG10_Prescl, &b_L1_DoubleIsoEG10_Prescl);
   fChain->SetBranchAddress("L1_DoubleIsoEG10_5bx", &L1_DoubleIsoEG10_5bx, &b_L1_DoubleIsoEG10_5bx);
   fChain->SetBranchAddress("L1_DoubleJet36_Central", &L1_DoubleJet36_Central, &b_L1_DoubleJet36_Central);
   fChain->SetBranchAddress("L1_DoubleJet36_Central_Prescl", &L1_DoubleJet36_Central_Prescl, &b_L1_DoubleJet36_Central_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet36_Central_5bx", &L1_DoubleJet36_Central_5bx, &b_L1_DoubleJet36_Central_5bx);
   fChain->SetBranchAddress("L1_DoubleJet44_Central", &L1_DoubleJet44_Central, &b_L1_DoubleJet44_Central);
   fChain->SetBranchAddress("L1_DoubleJet44_Central_Prescl", &L1_DoubleJet44_Central_Prescl, &b_L1_DoubleJet44_Central_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet44_Central_5bx", &L1_DoubleJet44_Central_5bx, &b_L1_DoubleJet44_Central_5bx);
   fChain->SetBranchAddress("L1_DoubleJet52", &L1_DoubleJet52, &b_L1_DoubleJet52);
   fChain->SetBranchAddress("L1_DoubleJet52_Prescl", &L1_DoubleJet52_Prescl, &b_L1_DoubleJet52_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet52_5bx", &L1_DoubleJet52_5bx, &b_L1_DoubleJet52_5bx);
   fChain->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   fChain->SetBranchAddress("L1_DoubleMu0_Prescl", &L1_DoubleMu0_Prescl, &b_L1_DoubleMu0_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_5bx", &L1_DoubleMu0_5bx, &b_L1_DoubleMu0_5bx);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ", &L1_DoubleMu0_HighQ, &b_L1_DoubleMu0_HighQ);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_Prescl", &L1_DoubleMu0_HighQ_Prescl, &b_L1_DoubleMu0_HighQ_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_5bx", &L1_DoubleMu0_HighQ_5bx, &b_L1_DoubleMu0_HighQ_5bx);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_EtaCuts", &L1_DoubleMu0_HighQ_EtaCuts, &b_L1_DoubleMu0_HighQ_EtaCuts);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_EtaCuts_Prescl", &L1_DoubleMu0_HighQ_EtaCuts_Prescl, &b_L1_DoubleMu0_HighQ_EtaCuts_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_HighQ_EtaCuts_5bx", &L1_DoubleMu0_HighQ_EtaCuts_5bx, &b_L1_DoubleMu0_HighQ_EtaCuts_5bx);
   fChain->SetBranchAddress("L1_DoubleMu3", &L1_DoubleMu3, &b_L1_DoubleMu3);
   fChain->SetBranchAddress("L1_DoubleMu3_Prescl", &L1_DoubleMu3_Prescl, &b_L1_DoubleMu3_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3_5bx", &L1_DoubleMu3_5bx, &b_L1_DoubleMu3_5bx);
   fChain->SetBranchAddress("L1_DoubleMu3_EG5", &L1_DoubleMu3_EG5, &b_L1_DoubleMu3_EG5);
   fChain->SetBranchAddress("L1_DoubleMu3_EG5_Prescl", &L1_DoubleMu3_EG5_Prescl, &b_L1_DoubleMu3_EG5_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3_EG5_5bx", &L1_DoubleMu3_EG5_5bx, &b_L1_DoubleMu3_EG5_5bx);
   fChain->SetBranchAddress("L1_DoubleMu3p5", &L1_DoubleMu3p5, &b_L1_DoubleMu3p5);
   fChain->SetBranchAddress("L1_DoubleMu3p5_Prescl", &L1_DoubleMu3p5_Prescl, &b_L1_DoubleMu3p5_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3p5_5bx", &L1_DoubleMu3p5_5bx, &b_L1_DoubleMu3p5_5bx);
   fChain->SetBranchAddress("L1_DoubleMu5_v1", &L1_DoubleMu5_v1, &b_L1_DoubleMu5_v1);
   fChain->SetBranchAddress("L1_DoubleMu5_v1_Prescl", &L1_DoubleMu5_v1_Prescl, &b_L1_DoubleMu5_v1_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu5_v1_5bx", &L1_DoubleMu5_v1_5bx, &b_L1_DoubleMu5_v1_5bx);
   fChain->SetBranchAddress("L1_DoubleTauJet28", &L1_DoubleTauJet28, &b_L1_DoubleTauJet28);
   fChain->SetBranchAddress("L1_DoubleTauJet28_Prescl", &L1_DoubleTauJet28_Prescl, &b_L1_DoubleTauJet28_Prescl);
   fChain->SetBranchAddress("L1_DoubleTauJet28_5bx", &L1_DoubleTauJet28_5bx, &b_L1_DoubleTauJet28_5bx);
   fChain->SetBranchAddress("L1_DoubleTauJet32", &L1_DoubleTauJet32, &b_L1_DoubleTauJet32);
   fChain->SetBranchAddress("L1_DoubleTauJet32_Prescl", &L1_DoubleTauJet32_Prescl, &b_L1_DoubleTauJet32_Prescl);
   fChain->SetBranchAddress("L1_DoubleTauJet32_5bx", &L1_DoubleTauJet32_5bx, &b_L1_DoubleTauJet32_5bx);
   fChain->SetBranchAddress("L1_DoubleTauJet36", &L1_DoubleTauJet36, &b_L1_DoubleTauJet36);
   fChain->SetBranchAddress("L1_DoubleTauJet36_Prescl", &L1_DoubleTauJet36_Prescl, &b_L1_DoubleTauJet36_Prescl);
   fChain->SetBranchAddress("L1_DoubleTauJet36_5bx", &L1_DoubleTauJet36_5bx, &b_L1_DoubleTauJet36_5bx);
   fChain->SetBranchAddress("L1_DoubleTauJet40", &L1_DoubleTauJet40, &b_L1_DoubleTauJet40);
   fChain->SetBranchAddress("L1_DoubleTauJet40_Prescl", &L1_DoubleTauJet40_Prescl, &b_L1_DoubleTauJet40_Prescl);
   fChain->SetBranchAddress("L1_DoubleTauJet40_5bx", &L1_DoubleTauJet40_5bx, &b_L1_DoubleTauJet40_5bx);
   fChain->SetBranchAddress("L1_EG10_Jet24_Central_deltaPhi1", &L1_EG10_Jet24_Central_deltaPhi1, &b_L1_EG10_Jet24_Central_deltaPhi1);
   fChain->SetBranchAddress("L1_EG10_Jet24_Central_deltaPhi1_Prescl", &L1_EG10_Jet24_Central_deltaPhi1_Prescl, &b_L1_EG10_Jet24_Central_deltaPhi1_Prescl);
   fChain->SetBranchAddress("L1_EG10_Jet24_Central_deltaPhi1_5bx", &L1_EG10_Jet24_Central_deltaPhi1_5bx, &b_L1_EG10_Jet24_Central_deltaPhi1_5bx);
   fChain->SetBranchAddress("L1_EG12_Jet24_Central_deltaPhi1", &L1_EG12_Jet24_Central_deltaPhi1, &b_L1_EG12_Jet24_Central_deltaPhi1);
   fChain->SetBranchAddress("L1_EG12_Jet24_Central_deltaPhi1_Prescl", &L1_EG12_Jet24_Central_deltaPhi1_Prescl, &b_L1_EG12_Jet24_Central_deltaPhi1_Prescl);
   fChain->SetBranchAddress("L1_EG12_Jet24_Central_deltaPhi1_5bx", &L1_EG12_Jet24_Central_deltaPhi1_5bx, &b_L1_EG12_Jet24_Central_deltaPhi1_5bx);
   fChain->SetBranchAddress("L1_EG12_TauJet20_deltaPhi1", &L1_EG12_TauJet20_deltaPhi1, &b_L1_EG12_TauJet20_deltaPhi1);
   fChain->SetBranchAddress("L1_EG12_TauJet20_deltaPhi1_Prescl", &L1_EG12_TauJet20_deltaPhi1_Prescl, &b_L1_EG12_TauJet20_deltaPhi1_Prescl);
   fChain->SetBranchAddress("L1_EG12_TauJet20_deltaPhi1_5bx", &L1_EG12_TauJet20_deltaPhi1_5bx, &b_L1_EG12_TauJet20_deltaPhi1_5bx);
   fChain->SetBranchAddress("L1_EG5_HTT100", &L1_EG5_HTT100, &b_L1_EG5_HTT100);
   fChain->SetBranchAddress("L1_EG5_HTT100_Prescl", &L1_EG5_HTT100_Prescl, &b_L1_EG5_HTT100_Prescl);
   fChain->SetBranchAddress("L1_EG5_HTT100_5bx", &L1_EG5_HTT100_5bx, &b_L1_EG5_HTT100_5bx);
   fChain->SetBranchAddress("L1_EG5_HTT125", &L1_EG5_HTT125, &b_L1_EG5_HTT125);
   fChain->SetBranchAddress("L1_EG5_HTT125_Prescl", &L1_EG5_HTT125_Prescl, &b_L1_EG5_HTT125_Prescl);
   fChain->SetBranchAddress("L1_EG5_HTT125_5bx", &L1_EG5_HTT125_5bx, &b_L1_EG5_HTT125_5bx);
   fChain->SetBranchAddress("L1_EG5_HTT75", &L1_EG5_HTT75, &b_L1_EG5_HTT75);
   fChain->SetBranchAddress("L1_EG5_HTT75_Prescl", &L1_EG5_HTT75_Prescl, &b_L1_EG5_HTT75_Prescl);
   fChain->SetBranchAddress("L1_EG5_HTT75_5bx", &L1_EG5_HTT75_5bx, &b_L1_EG5_HTT75_5bx);
   fChain->SetBranchAddress("L1_EG5_Jet36_deltaPhi1", &L1_EG5_Jet36_deltaPhi1, &b_L1_EG5_Jet36_deltaPhi1);
   fChain->SetBranchAddress("L1_EG5_Jet36_deltaPhi1_Prescl", &L1_EG5_Jet36_deltaPhi1_Prescl, &b_L1_EG5_Jet36_deltaPhi1_Prescl);
   fChain->SetBranchAddress("L1_EG5_Jet36_deltaPhi1_5bx", &L1_EG5_Jet36_deltaPhi1_5bx, &b_L1_EG5_Jet36_deltaPhi1_5bx);
   fChain->SetBranchAddress("L1_EG8_Jet20_Central_deltaPhi1", &L1_EG8_Jet20_Central_deltaPhi1, &b_L1_EG8_Jet20_Central_deltaPhi1);
   fChain->SetBranchAddress("L1_EG8_Jet20_Central_deltaPhi1_Prescl", &L1_EG8_Jet20_Central_deltaPhi1_Prescl, &b_L1_EG8_Jet20_Central_deltaPhi1_Prescl);
   fChain->SetBranchAddress("L1_EG8_Jet20_Central_deltaPhi1_5bx", &L1_EG8_Jet20_Central_deltaPhi1_5bx, &b_L1_EG8_Jet20_Central_deltaPhi1_5bx);
   fChain->SetBranchAddress("L1_ETM100", &L1_ETM100, &b_L1_ETM100);
   fChain->SetBranchAddress("L1_ETM100_Prescl", &L1_ETM100_Prescl, &b_L1_ETM100_Prescl);
   fChain->SetBranchAddress("L1_ETM100_5bx", &L1_ETM100_5bx, &b_L1_ETM100_5bx);
   fChain->SetBranchAddress("L1_ETM20", &L1_ETM20, &b_L1_ETM20);
   fChain->SetBranchAddress("L1_ETM20_Prescl", &L1_ETM20_Prescl, &b_L1_ETM20_Prescl);
   fChain->SetBranchAddress("L1_ETM20_5bx", &L1_ETM20_5bx, &b_L1_ETM20_5bx);
   fChain->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   fChain->SetBranchAddress("L1_ETM30_Prescl", &L1_ETM30_Prescl, &b_L1_ETM30_Prescl);
   fChain->SetBranchAddress("L1_ETM30_5bx", &L1_ETM30_5bx, &b_L1_ETM30_5bx);
   fChain->SetBranchAddress("L1_ETM50", &L1_ETM50, &b_L1_ETM50);
   fChain->SetBranchAddress("L1_ETM50_Prescl", &L1_ETM50_Prescl, &b_L1_ETM50_Prescl);
   fChain->SetBranchAddress("L1_ETM50_5bx", &L1_ETM50_5bx, &b_L1_ETM50_5bx);
   fChain->SetBranchAddress("L1_ETM70", &L1_ETM70, &b_L1_ETM70);
   fChain->SetBranchAddress("L1_ETM70_Prescl", &L1_ETM70_Prescl, &b_L1_ETM70_Prescl);
   fChain->SetBranchAddress("L1_ETM70_5bx", &L1_ETM70_5bx, &b_L1_ETM70_5bx);
   fChain->SetBranchAddress("L1_ETT220", &L1_ETT220, &b_L1_ETT220);
   fChain->SetBranchAddress("L1_ETT220_Prescl", &L1_ETT220_Prescl, &b_L1_ETT220_Prescl);
   fChain->SetBranchAddress("L1_ETT220_5bx", &L1_ETT220_5bx, &b_L1_ETT220_5bx);
   fChain->SetBranchAddress("L1_ETT260_EG5", &L1_ETT260_EG5, &b_L1_ETT260_EG5);
   fChain->SetBranchAddress("L1_ETT260_EG5_Prescl", &L1_ETT260_EG5_Prescl, &b_L1_ETT260_EG5_Prescl);
   fChain->SetBranchAddress("L1_ETT260_EG5_5bx", &L1_ETT260_EG5_5bx, &b_L1_ETT260_EG5_5bx);
   fChain->SetBranchAddress("L1_ETT300_EG5", &L1_ETT300_EG5, &b_L1_ETT300_EG5);
   fChain->SetBranchAddress("L1_ETT300_EG5_Prescl", &L1_ETT300_EG5_Prescl, &b_L1_ETT300_EG5_Prescl);
   fChain->SetBranchAddress("L1_ETT300_EG5_5bx", &L1_ETT300_EG5_5bx, &b_L1_ETT300_EG5_5bx);
   fChain->SetBranchAddress("L1_HTM50", &L1_HTM50, &b_L1_HTM50);
   fChain->SetBranchAddress("L1_HTM50_Prescl", &L1_HTM50_Prescl, &b_L1_HTM50_Prescl);
   fChain->SetBranchAddress("L1_HTM50_5bx", &L1_HTM50_5bx, &b_L1_HTM50_5bx);
   fChain->SetBranchAddress("L1_HTT100", &L1_HTT100, &b_L1_HTT100);
   fChain->SetBranchAddress("L1_HTT100_Prescl", &L1_HTT100_Prescl, &b_L1_HTT100_Prescl);
   fChain->SetBranchAddress("L1_HTT100_5bx", &L1_HTT100_5bx, &b_L1_HTT100_5bx);
   fChain->SetBranchAddress("L1_HTT150", &L1_HTT150, &b_L1_HTT150);
   fChain->SetBranchAddress("L1_HTT150_Prescl", &L1_HTT150_Prescl, &b_L1_HTT150_Prescl);
   fChain->SetBranchAddress("L1_HTT150_5bx", &L1_HTT150_5bx, &b_L1_HTT150_5bx);
   fChain->SetBranchAddress("L1_HTT50", &L1_HTT50, &b_L1_HTT50);
   fChain->SetBranchAddress("L1_HTT50_Prescl", &L1_HTT50_Prescl, &b_L1_HTT50_Prescl);
   fChain->SetBranchAddress("L1_HTT50_5bx", &L1_HTT50_5bx, &b_L1_HTT50_5bx);
   fChain->SetBranchAddress("L1_HTT50_HTM30", &L1_HTT50_HTM30, &b_L1_HTT50_HTM30);
   fChain->SetBranchAddress("L1_HTT50_HTM30_Prescl", &L1_HTT50_HTM30_Prescl, &b_L1_HTT50_HTM30_Prescl);
   fChain->SetBranchAddress("L1_HTT50_HTM30_5bx", &L1_HTT50_HTM30_5bx, &b_L1_HTT50_HTM30_5bx);
   fChain->SetBranchAddress("L1_HTT50_HTM50", &L1_HTT50_HTM50, &b_L1_HTT50_HTM50);
   fChain->SetBranchAddress("L1_HTT50_HTM50_Prescl", &L1_HTT50_HTM50_Prescl, &b_L1_HTT50_HTM50_Prescl);
   fChain->SetBranchAddress("L1_HTT50_HTM50_5bx", &L1_HTT50_HTM50_5bx, &b_L1_HTT50_HTM50_5bx);
   fChain->SetBranchAddress("L1_HTT75", &L1_HTT75, &b_L1_HTT75);
   fChain->SetBranchAddress("L1_HTT75_Prescl", &L1_HTT75_Prescl, &b_L1_HTT75_Prescl);
   fChain->SetBranchAddress("L1_HTT75_5bx", &L1_HTT75_5bx, &b_L1_HTT75_5bx);
   fChain->SetBranchAddress("L1_InterBunch_Bsc", &L1_InterBunch_Bsc, &b_L1_InterBunch_Bsc);
   fChain->SetBranchAddress("L1_InterBunch_Bsc_Prescl", &L1_InterBunch_Bsc_Prescl, &b_L1_InterBunch_Bsc_Prescl);
   fChain->SetBranchAddress("L1_InterBunch_Bsc_5bx", &L1_InterBunch_Bsc_5bx, &b_L1_InterBunch_Bsc_5bx);
   fChain->SetBranchAddress("L1_InterBunch_Hf", &L1_InterBunch_Hf, &b_L1_InterBunch_Hf);
   fChain->SetBranchAddress("L1_InterBunch_Hf_Prescl", &L1_InterBunch_Hf_Prescl, &b_L1_InterBunch_Hf_Prescl);
   fChain->SetBranchAddress("L1_InterBunch_Hf_5bx", &L1_InterBunch_Hf_5bx, &b_L1_InterBunch_Hf_5bx);
   fChain->SetBranchAddress("L1_Mu0_HTT50", &L1_Mu0_HTT50, &b_L1_Mu0_HTT50);
   fChain->SetBranchAddress("L1_Mu0_HTT50_Prescl", &L1_Mu0_HTT50_Prescl, &b_L1_Mu0_HTT50_Prescl);
   fChain->SetBranchAddress("L1_Mu0_HTT50_5bx", &L1_Mu0_HTT50_5bx, &b_L1_Mu0_HTT50_5bx);
   fChain->SetBranchAddress("L1_Mu0_HTT75", &L1_Mu0_HTT75, &b_L1_Mu0_HTT75);
   fChain->SetBranchAddress("L1_Mu0_HTT75_Prescl", &L1_Mu0_HTT75_Prescl, &b_L1_Mu0_HTT75_Prescl);
   fChain->SetBranchAddress("L1_Mu0_HTT75_5bx", &L1_Mu0_HTT75_5bx, &b_L1_Mu0_HTT75_5bx);
   fChain->SetBranchAddress("L1_Mu10_Jet36_Central", &L1_Mu10_Jet36_Central, &b_L1_Mu10_Jet36_Central);
   fChain->SetBranchAddress("L1_Mu10_Jet36_Central_Prescl", &L1_Mu10_Jet36_Central_Prescl, &b_L1_Mu10_Jet36_Central_Prescl);
   fChain->SetBranchAddress("L1_Mu10_Jet36_Central_5bx", &L1_Mu10_Jet36_Central_5bx, &b_L1_Mu10_Jet36_Central_5bx);
   fChain->SetBranchAddress("L1_Mu12_EG5", &L1_Mu12_EG5, &b_L1_Mu12_EG5);
   fChain->SetBranchAddress("L1_Mu12_EG5_Prescl", &L1_Mu12_EG5_Prescl, &b_L1_Mu12_EG5_Prescl);
   fChain->SetBranchAddress("L1_Mu12_EG5_5bx", &L1_Mu12_EG5_5bx, &b_L1_Mu12_EG5_5bx);
   fChain->SetBranchAddress("L1_Mu3_DoubleEG5", &L1_Mu3_DoubleEG5, &b_L1_Mu3_DoubleEG5);
   fChain->SetBranchAddress("L1_Mu3_DoubleEG5_Prescl", &L1_Mu3_DoubleEG5_Prescl, &b_L1_Mu3_DoubleEG5_Prescl);
   fChain->SetBranchAddress("L1_Mu3_DoubleEG5_5bx", &L1_Mu3_DoubleEG5_5bx, &b_L1_Mu3_DoubleEG5_5bx);
   fChain->SetBranchAddress("L1_Mu3_EG5", &L1_Mu3_EG5, &b_L1_Mu3_EG5);
   fChain->SetBranchAddress("L1_Mu3_EG5_Prescl", &L1_Mu3_EG5_Prescl, &b_L1_Mu3_EG5_Prescl);
   fChain->SetBranchAddress("L1_Mu3_EG5_5bx", &L1_Mu3_EG5_5bx, &b_L1_Mu3_EG5_5bx);
   fChain->SetBranchAddress("L1_Mu3_Jet16_Central", &L1_Mu3_Jet16_Central, &b_L1_Mu3_Jet16_Central);
   fChain->SetBranchAddress("L1_Mu3_Jet16_Central_Prescl", &L1_Mu3_Jet16_Central_Prescl, &b_L1_Mu3_Jet16_Central_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet16_Central_5bx", &L1_Mu3_Jet16_Central_5bx, &b_L1_Mu3_Jet16_Central_5bx);
   fChain->SetBranchAddress("L1_Mu3_Jet20_Central", &L1_Mu3_Jet20_Central, &b_L1_Mu3_Jet20_Central);
   fChain->SetBranchAddress("L1_Mu3_Jet20_Central_Prescl", &L1_Mu3_Jet20_Central_Prescl, &b_L1_Mu3_Jet20_Central_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet20_Central_5bx", &L1_Mu3_Jet20_Central_5bx, &b_L1_Mu3_Jet20_Central_5bx);
   fChain->SetBranchAddress("L1_Mu3_Jet28_Central", &L1_Mu3_Jet28_Central, &b_L1_Mu3_Jet28_Central);
   fChain->SetBranchAddress("L1_Mu3_Jet28_Central_Prescl", &L1_Mu3_Jet28_Central_Prescl, &b_L1_Mu3_Jet28_Central_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet28_Central_5bx", &L1_Mu3_Jet28_Central_5bx, &b_L1_Mu3_Jet28_Central_5bx);
   fChain->SetBranchAddress("L1_Mu5_EG12", &L1_Mu5_EG12, &b_L1_Mu5_EG12);
   fChain->SetBranchAddress("L1_Mu5_EG12_Prescl", &L1_Mu5_EG12_Prescl, &b_L1_Mu5_EG12_Prescl);
   fChain->SetBranchAddress("L1_Mu5_EG12_5bx", &L1_Mu5_EG12_5bx, &b_L1_Mu5_EG12_5bx);
   fChain->SetBranchAddress("L1_Mu7_EG5", &L1_Mu7_EG5, &b_L1_Mu7_EG5);
   fChain->SetBranchAddress("L1_Mu7_EG5_Prescl", &L1_Mu7_EG5_Prescl, &b_L1_Mu7_EG5_Prescl);
   fChain->SetBranchAddress("L1_Mu7_EG5_5bx", &L1_Mu7_EG5_5bx, &b_L1_Mu7_EG5_5bx);
   fChain->SetBranchAddress("L1_Mu7_Jet20_Central", &L1_Mu7_Jet20_Central, &b_L1_Mu7_Jet20_Central);
   fChain->SetBranchAddress("L1_Mu7_Jet20_Central_Prescl", &L1_Mu7_Jet20_Central_Prescl, &b_L1_Mu7_Jet20_Central_Prescl);
   fChain->SetBranchAddress("L1_Mu7_Jet20_Central_5bx", &L1_Mu7_Jet20_Central_5bx, &b_L1_Mu7_Jet20_Central_5bx);
   fChain->SetBranchAddress("L1_Mu7_TauJet16", &L1_Mu7_TauJet16, &b_L1_Mu7_TauJet16);
   fChain->SetBranchAddress("L1_Mu7_TauJet16_Prescl", &L1_Mu7_TauJet16_Prescl, &b_L1_Mu7_TauJet16_Prescl);
   fChain->SetBranchAddress("L1_Mu7_TauJet16_5bx", &L1_Mu7_TauJet16_5bx, &b_L1_Mu7_TauJet16_5bx);
   fChain->SetBranchAddress("L1_MuOpen_EG12", &L1_MuOpen_EG12, &b_L1_MuOpen_EG12);
   fChain->SetBranchAddress("L1_MuOpen_EG12_Prescl", &L1_MuOpen_EG12_Prescl, &b_L1_MuOpen_EG12_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_EG12_5bx", &L1_MuOpen_EG12_5bx, &b_L1_MuOpen_EG12_5bx);
   fChain->SetBranchAddress("L1_MuOpen_EG5", &L1_MuOpen_EG5, &b_L1_MuOpen_EG5);
   fChain->SetBranchAddress("L1_MuOpen_EG5_Prescl", &L1_MuOpen_EG5_Prescl, &b_L1_MuOpen_EG5_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_EG5_5bx", &L1_MuOpen_EG5_5bx, &b_L1_MuOpen_EG5_5bx);
   fChain->SetBranchAddress("L1_PreCollisions", &L1_PreCollisions, &b_L1_PreCollisions);
   fChain->SetBranchAddress("L1_PreCollisions_Prescl", &L1_PreCollisions_Prescl, &b_L1_PreCollisions_Prescl);
   fChain->SetBranchAddress("L1_PreCollisions_5bx", &L1_PreCollisions_5bx, &b_L1_PreCollisions_5bx);
   fChain->SetBranchAddress("L1_QuadJet20_Central", &L1_QuadJet20_Central, &b_L1_QuadJet20_Central);
   fChain->SetBranchAddress("L1_QuadJet20_Central_Prescl", &L1_QuadJet20_Central_Prescl, &b_L1_QuadJet20_Central_Prescl);
   fChain->SetBranchAddress("L1_QuadJet20_Central_5bx", &L1_QuadJet20_Central_5bx, &b_L1_QuadJet20_Central_5bx);
   fChain->SetBranchAddress("L1_QuadJet28_Central", &L1_QuadJet28_Central, &b_L1_QuadJet28_Central);
   fChain->SetBranchAddress("L1_QuadJet28_Central_Prescl", &L1_QuadJet28_Central_Prescl, &b_L1_QuadJet28_Central_Prescl);
   fChain->SetBranchAddress("L1_QuadJet28_Central_5bx", &L1_QuadJet28_Central_5bx, &b_L1_QuadJet28_Central_5bx);
   fChain->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12, &b_L1_SingleEG12);
   fChain->SetBranchAddress("L1_SingleEG12_Prescl", &L1_SingleEG12_Prescl, &b_L1_SingleEG12_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12_5bx", &L1_SingleEG12_5bx, &b_L1_SingleEG12_5bx);
   fChain->SetBranchAddress("L1_SingleEG12_Eta1p39", &L1_SingleEG12_Eta1p39, &b_L1_SingleEG12_Eta1p39);
   fChain->SetBranchAddress("L1_SingleEG12_Eta1p39_Prescl", &L1_SingleEG12_Eta1p39_Prescl, &b_L1_SingleEG12_Eta1p39_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12_Eta1p39_5bx", &L1_SingleEG12_Eta1p39_5bx, &b_L1_SingleEG12_Eta1p39_5bx);
   fChain->SetBranchAddress("L1_SingleEG12_Eta2p17", &L1_SingleEG12_Eta2p17, &b_L1_SingleEG12_Eta2p17);
   fChain->SetBranchAddress("L1_SingleEG12_Eta2p17_Prescl", &L1_SingleEG12_Eta2p17_Prescl, &b_L1_SingleEG12_Eta2p17_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12_Eta2p17_5bx", &L1_SingleEG12_Eta2p17_5bx, &b_L1_SingleEG12_Eta2p17_5bx);
   fChain->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   fChain->SetBranchAddress("L1_SingleEG15_Prescl", &L1_SingleEG15_Prescl, &b_L1_SingleEG15_Prescl);
   fChain->SetBranchAddress("L1_SingleEG15_5bx", &L1_SingleEG15_5bx, &b_L1_SingleEG15_5bx);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG20_Prescl", &L1_SingleEG20_Prescl, &b_L1_SingleEG20_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20_5bx", &L1_SingleEG20_5bx, &b_L1_SingleEG20_5bx);
   fChain->SetBranchAddress("L1_SingleEG30", &L1_SingleEG30, &b_L1_SingleEG30);
   fChain->SetBranchAddress("L1_SingleEG30_Prescl", &L1_SingleEG30_Prescl, &b_L1_SingleEG30_Prescl);
   fChain->SetBranchAddress("L1_SingleEG30_5bx", &L1_SingleEG30_5bx, &b_L1_SingleEG30_5bx);
   fChain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   fChain->SetBranchAddress("L1_SingleEG5_Prescl", &L1_SingleEG5_Prescl, &b_L1_SingleEG5_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5_5bx", &L1_SingleEG5_5bx, &b_L1_SingleEG5_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG12", &L1_SingleIsoEG12, &b_L1_SingleIsoEG12);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Prescl", &L1_SingleIsoEG12_Prescl, &b_L1_SingleIsoEG12_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG12_5bx", &L1_SingleIsoEG12_5bx, &b_L1_SingleIsoEG12_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Eta1p39", &L1_SingleIsoEG12_Eta1p39, &b_L1_SingleIsoEG12_Eta1p39);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Eta1p39_Prescl", &L1_SingleIsoEG12_Eta1p39_Prescl, &b_L1_SingleIsoEG12_Eta1p39_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Eta1p39_5bx", &L1_SingleIsoEG12_Eta1p39_5bx, &b_L1_SingleIsoEG12_Eta1p39_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Eta2p17", &L1_SingleIsoEG12_Eta2p17, &b_L1_SingleIsoEG12_Eta2p17);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Eta2p17_Prescl", &L1_SingleIsoEG12_Eta2p17_Prescl, &b_L1_SingleIsoEG12_Eta2p17_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Eta2p17_5bx", &L1_SingleIsoEG12_Eta2p17_5bx, &b_L1_SingleIsoEG12_Eta2p17_5bx);
   fChain->SetBranchAddress("L1_SingleJet128", &L1_SingleJet128, &b_L1_SingleJet128);
   fChain->SetBranchAddress("L1_SingleJet128_Prescl", &L1_SingleJet128_Prescl, &b_L1_SingleJet128_Prescl);
   fChain->SetBranchAddress("L1_SingleJet128_5bx", &L1_SingleJet128_5bx, &b_L1_SingleJet128_5bx);
   fChain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   fChain->SetBranchAddress("L1_SingleJet16_Prescl", &L1_SingleJet16_Prescl, &b_L1_SingleJet16_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_5bx", &L1_SingleJet16_5bx, &b_L1_SingleJet16_5bx);
   fChain->SetBranchAddress("L1_SingleJet20_NotBptxOR", &L1_SingleJet20_NotBptxOR, &b_L1_SingleJet20_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJet20_NotBptxOR_Prescl", &L1_SingleJet20_NotBptxOR_Prescl, &b_L1_SingleJet20_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_NotBptxOR_5bx", &L1_SingleJet20_NotBptxOR_5bx, &b_L1_SingleJet20_NotBptxOR_5bx);
   fChain->SetBranchAddress("L1_SingleJet20_NotBptxOR_NotMuBeamHalo", &L1_SingleJet20_NotBptxOR_NotMuBeamHalo, &b_L1_SingleJet20_NotBptxOR_NotMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleJet20_NotBptxOR_NotMuBeamHalo_Prescl", &L1_SingleJet20_NotBptxOR_NotMuBeamHalo_Prescl, &b_L1_SingleJet20_NotBptxOR_NotMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_NotBptxOR_NotMuBeamHalo_5bx", &L1_SingleJet20_NotBptxOR_NotMuBeamHalo_5bx, &b_L1_SingleJet20_NotBptxOR_NotMuBeamHalo_5bx);
   fChain->SetBranchAddress("L1_SingleJet32_NotBptxOR_NotMuBeamHalo", &L1_SingleJet32_NotBptxOR_NotMuBeamHalo, &b_L1_SingleJet32_NotBptxOR_NotMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleJet32_NotBptxOR_NotMuBeamHalo_Prescl", &L1_SingleJet32_NotBptxOR_NotMuBeamHalo_Prescl, &b_L1_SingleJet32_NotBptxOR_NotMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleJet32_NotBptxOR_NotMuBeamHalo_5bx", &L1_SingleJet32_NotBptxOR_NotMuBeamHalo_5bx, &b_L1_SingleJet32_NotBptxOR_NotMuBeamHalo_5bx);
   fChain->SetBranchAddress("L1_SingleJet36", &L1_SingleJet36, &b_L1_SingleJet36);
   fChain->SetBranchAddress("L1_SingleJet36_Prescl", &L1_SingleJet36_Prescl, &b_L1_SingleJet36_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36_5bx", &L1_SingleJet36_5bx, &b_L1_SingleJet36_5bx);
   fChain->SetBranchAddress("L1_SingleJet36_FwdVeto", &L1_SingleJet36_FwdVeto, &b_L1_SingleJet36_FwdVeto);
   fChain->SetBranchAddress("L1_SingleJet36_FwdVeto_Prescl", &L1_SingleJet36_FwdVeto_Prescl, &b_L1_SingleJet36_FwdVeto_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36_FwdVeto_5bx", &L1_SingleJet36_FwdVeto_5bx, &b_L1_SingleJet36_FwdVeto_5bx);
   fChain->SetBranchAddress("L1_SingleJet52", &L1_SingleJet52, &b_L1_SingleJet52);
   fChain->SetBranchAddress("L1_SingleJet52_Prescl", &L1_SingleJet52_Prescl, &b_L1_SingleJet52_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52_5bx", &L1_SingleJet52_5bx, &b_L1_SingleJet52_5bx);
   fChain->SetBranchAddress("L1_SingleJet68", &L1_SingleJet68, &b_L1_SingleJet68);
   fChain->SetBranchAddress("L1_SingleJet68_Prescl", &L1_SingleJet68_Prescl, &b_L1_SingleJet68_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68_5bx", &L1_SingleJet68_5bx, &b_L1_SingleJet68_5bx);
   fChain->SetBranchAddress("L1_SingleJet80_Central", &L1_SingleJet80_Central, &b_L1_SingleJet80_Central);
   fChain->SetBranchAddress("L1_SingleJet80_Central_Prescl", &L1_SingleJet80_Central_Prescl, &b_L1_SingleJet80_Central_Prescl);
   fChain->SetBranchAddress("L1_SingleJet80_Central_5bx", &L1_SingleJet80_Central_5bx, &b_L1_SingleJet80_Central_5bx);
   fChain->SetBranchAddress("L1_SingleJet92", &L1_SingleJet92, &b_L1_SingleJet92);
   fChain->SetBranchAddress("L1_SingleJet92_Prescl", &L1_SingleJet92_Prescl, &b_L1_SingleJet92_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92_5bx", &L1_SingleJet92_5bx, &b_L1_SingleJet92_5bx);
   fChain->SetBranchAddress("L1_SingleJet92_Central", &L1_SingleJet92_Central, &b_L1_SingleJet92_Central);
   fChain->SetBranchAddress("L1_SingleJet92_Central_Prescl", &L1_SingleJet92_Central_Prescl, &b_L1_SingleJet92_Central_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92_Central_5bx", &L1_SingleJet92_Central_5bx, &b_L1_SingleJet92_Central_5bx);
   fChain->SetBranchAddress("L1_SingleMu10", &L1_SingleMu10, &b_L1_SingleMu10);
   fChain->SetBranchAddress("L1_SingleMu10_Prescl", &L1_SingleMu10_Prescl, &b_L1_SingleMu10_Prescl);
   fChain->SetBranchAddress("L1_SingleMu10_5bx", &L1_SingleMu10_5bx, &b_L1_SingleMu10_5bx);
   fChain->SetBranchAddress("L1_SingleMu12", &L1_SingleMu12, &b_L1_SingleMu12);
   fChain->SetBranchAddress("L1_SingleMu12_Prescl", &L1_SingleMu12_Prescl, &b_L1_SingleMu12_Prescl);
   fChain->SetBranchAddress("L1_SingleMu12_5bx", &L1_SingleMu12_5bx, &b_L1_SingleMu12_5bx);
   fChain->SetBranchAddress("L1_SingleMu12_Debug", &L1_SingleMu12_Debug, &b_L1_SingleMu12_Debug);
   fChain->SetBranchAddress("L1_SingleMu12_Debug_Prescl", &L1_SingleMu12_Debug_Prescl, &b_L1_SingleMu12_Debug_Prescl);
   fChain->SetBranchAddress("L1_SingleMu12_Debug_5bx", &L1_SingleMu12_Debug_5bx, &b_L1_SingleMu12_Debug_5bx);
   fChain->SetBranchAddress("L1_SingleMu16", &L1_SingleMu16, &b_L1_SingleMu16);
   fChain->SetBranchAddress("L1_SingleMu16_Prescl", &L1_SingleMu16_Prescl, &b_L1_SingleMu16_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16_5bx", &L1_SingleMu16_5bx, &b_L1_SingleMu16_5bx);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu20_Prescl", &L1_SingleMu20_Prescl, &b_L1_SingleMu20_Prescl);
   fChain->SetBranchAddress("L1_SingleMu20_5bx", &L1_SingleMu20_5bx, &b_L1_SingleMu20_5bx);
   fChain->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   fChain->SetBranchAddress("L1_SingleMu25_Prescl", &L1_SingleMu25_Prescl, &b_L1_SingleMu25_Prescl);
   fChain->SetBranchAddress("L1_SingleMu25_5bx", &L1_SingleMu25_5bx, &b_L1_SingleMu25_5bx);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu3_Prescl", &L1_SingleMu3_Prescl, &b_L1_SingleMu3_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_5bx", &L1_SingleMu3_5bx, &b_L1_SingleMu3_5bx);
   fChain->SetBranchAddress("L1_SingleMu5_Eta1p5_Q80_v1", &L1_SingleMu5_Eta1p5_Q80_v1, &b_L1_SingleMu5_Eta1p5_Q80_v1);
   fChain->SetBranchAddress("L1_SingleMu5_Eta1p5_Q80_v1_Prescl", &L1_SingleMu5_Eta1p5_Q80_v1_Prescl, &b_L1_SingleMu5_Eta1p5_Q80_v1_Prescl);
   fChain->SetBranchAddress("L1_SingleMu5_Eta1p5_Q80_v1_5bx", &L1_SingleMu5_Eta1p5_Q80_v1_5bx, &b_L1_SingleMu5_Eta1p5_Q80_v1_5bx);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu7_Prescl", &L1_SingleMu7_Prescl, &b_L1_SingleMu7_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7_5bx", &L1_SingleMu7_5bx, &b_L1_SingleMu7_5bx);
   fChain->SetBranchAddress("L1_SingleMu7_Barrel", &L1_SingleMu7_Barrel, &b_L1_SingleMu7_Barrel);
   fChain->SetBranchAddress("L1_SingleMu7_Barrel_Prescl", &L1_SingleMu7_Barrel_Prescl, &b_L1_SingleMu7_Barrel_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7_Barrel_5bx", &L1_SingleMu7_Barrel_5bx, &b_L1_SingleMu7_Barrel_5bx);
   fChain->SetBranchAddress("L1_SingleMu7_Eta2p1", &L1_SingleMu7_Eta2p1, &b_L1_SingleMu7_Eta2p1);
   fChain->SetBranchAddress("L1_SingleMu7_Eta2p1_Prescl", &L1_SingleMu7_Eta2p1_Prescl, &b_L1_SingleMu7_Eta2p1_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7_Eta2p1_5bx", &L1_SingleMu7_Eta2p1_5bx, &b_L1_SingleMu7_Eta2p1_5bx);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo", &L1_SingleMuBeamHalo, &b_L1_SingleMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &L1_SingleMuBeamHalo_Prescl, &b_L1_SingleMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_5bx", &L1_SingleMuBeamHalo_5bx, &b_L1_SingleMuBeamHalo_5bx);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_Prescl", &L1_SingleMuOpen_Prescl, &b_L1_SingleMuOpen_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_5bx", &L1_SingleMuOpen_5bx, &b_L1_SingleMuOpen_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet52", &L1_SingleTauJet52, &b_L1_SingleTauJet52);
   fChain->SetBranchAddress("L1_SingleTauJet52_Prescl", &L1_SingleTauJet52_Prescl, &b_L1_SingleTauJet52_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet52_5bx", &L1_SingleTauJet52_5bx, &b_L1_SingleTauJet52_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet68", &L1_SingleTauJet68, &b_L1_SingleTauJet68);
   fChain->SetBranchAddress("L1_SingleTauJet68_Prescl", &L1_SingleTauJet68_Prescl, &b_L1_SingleTauJet68_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet68_5bx", &L1_SingleTauJet68_5bx, &b_L1_SingleTauJet68_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet80", &L1_SingleTauJet80, &b_L1_SingleTauJet80);
   fChain->SetBranchAddress("L1_SingleTauJet80_Prescl", &L1_SingleTauJet80_Prescl, &b_L1_SingleTauJet80_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet80_5bx", &L1_SingleTauJet80_5bx, &b_L1_SingleTauJet80_5bx);
   fChain->SetBranchAddress("L1_TripleEG5", &L1_TripleEG5, &b_L1_TripleEG5);
   fChain->SetBranchAddress("L1_TripleEG5_Prescl", &L1_TripleEG5_Prescl, &b_L1_TripleEG5_Prescl);
   fChain->SetBranchAddress("L1_TripleEG5_5bx", &L1_TripleEG5_5bx, &b_L1_TripleEG5_5bx);
   fChain->SetBranchAddress("L1_TripleEG7", &L1_TripleEG7, &b_L1_TripleEG7);
   fChain->SetBranchAddress("L1_TripleEG7_Prescl", &L1_TripleEG7_Prescl, &b_L1_TripleEG7_Prescl);
   fChain->SetBranchAddress("L1_TripleEG7_5bx", &L1_TripleEG7_5bx, &b_L1_TripleEG7_5bx);
   fChain->SetBranchAddress("L1_TripleEG_8_5_5", &L1_TripleEG_8_5_5, &b_L1_TripleEG_8_5_5);
   fChain->SetBranchAddress("L1_TripleEG_8_5_5_Prescl", &L1_TripleEG_8_5_5_Prescl, &b_L1_TripleEG_8_5_5_Prescl);
   fChain->SetBranchAddress("L1_TripleEG_8_5_5_5bx", &L1_TripleEG_8_5_5_5bx, &b_L1_TripleEG_8_5_5_5bx);
   fChain->SetBranchAddress("L1_TripleEG_8_8_5", &L1_TripleEG_8_8_5, &b_L1_TripleEG_8_8_5);
   fChain->SetBranchAddress("L1_TripleEG_8_8_5_Prescl", &L1_TripleEG_8_8_5_Prescl, &b_L1_TripleEG_8_8_5_Prescl);
   fChain->SetBranchAddress("L1_TripleEG_8_8_5_5bx", &L1_TripleEG_8_8_5_5bx, &b_L1_TripleEG_8_8_5_5bx);
   fChain->SetBranchAddress("L1_TripleJet28_Central", &L1_TripleJet28_Central, &b_L1_TripleJet28_Central);
   fChain->SetBranchAddress("L1_TripleJet28_Central_Prescl", &L1_TripleJet28_Central_Prescl, &b_L1_TripleJet28_Central_Prescl);
   fChain->SetBranchAddress("L1_TripleJet28_Central_5bx", &L1_TripleJet28_Central_5bx, &b_L1_TripleJet28_Central_5bx);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_Prescl", &L1_ZeroBias_Prescl, &b_L1_ZeroBias_Prescl);
   fChain->SetBranchAddress("L1_ZeroBias_5bx", &L1_ZeroBias_5bx, &b_L1_ZeroBias_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0", &L1Tech_BPTX_minus_v0, &b_L1Tech_BPTX_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0_Prescl", &L1Tech_BPTX_minus_v0_Prescl, &b_L1Tech_BPTX_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0_5bx", &L1Tech_BPTX_minus_v0_5bx, &b_L1Tech_BPTX_minus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0", &L1Tech_BPTX_minus_AND_not_plus_v0, &b_L1Tech_BPTX_minus_AND_not_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl", &L1Tech_BPTX_minus_AND_not_plus_v0_Prescl, &b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_5bx", &L1Tech_BPTX_minus_AND_not_plus_v0_5bx, &b_L1Tech_BPTX_minus_AND_not_plus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0", &L1Tech_BPTX_plus_v0, &b_L1Tech_BPTX_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0_Prescl", &L1Tech_BPTX_plus_v0_Prescl, &b_L1Tech_BPTX_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0_5bx", &L1Tech_BPTX_plus_v0_5bx, &b_L1Tech_BPTX_plus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0", &L1Tech_BPTX_plus_AND_NOT_minus_v0, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_5bx", &L1Tech_BPTX_plus_AND_NOT_minus_v0_5bx, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &L1Tech_BPTX_plus_AND_minus_v0, &b_L1Tech_BPTX_plus_AND_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_5bx", &L1Tech_BPTX_plus_AND_minus_v0_5bx, &b_L1Tech_BPTX_plus_AND_minus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0", &L1Tech_BPTX_plus_AND_minus_instance1_v0, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_5bx", &L1Tech_BPTX_plus_AND_minus_instance1_v0_5bx, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &L1Tech_BPTX_plus_OR_minus_v0, &b_L1Tech_BPTX_plus_OR_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_Prescl", &L1Tech_BPTX_plus_OR_minus_v0_Prescl, &b_L1Tech_BPTX_plus_OR_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_5bx", &L1Tech_BPTX_plus_OR_minus_v0_5bx, &b_L1Tech_BPTX_plus_OR_minus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0", &L1Tech_BPTX_quiet_v0, &b_L1Tech_BPTX_quiet_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0_Prescl", &L1Tech_BPTX_quiet_v0_Prescl, &b_L1Tech_BPTX_quiet_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0_5bx", &L1Tech_BPTX_quiet_v0_5bx, &b_L1Tech_BPTX_quiet_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0", &L1Tech_BSC_HighMultiplicity_v0, &b_L1Tech_BSC_HighMultiplicity_v0);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0_Prescl", &L1Tech_BSC_HighMultiplicity_v0_Prescl, &b_L1Tech_BSC_HighMultiplicity_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0_5bx", &L1Tech_BSC_HighMultiplicity_v0_5bx, &b_L1Tech_BSC_HighMultiplicity_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0", &L1Tech_BSC_halo_beam1_inner_v0, &b_L1Tech_BSC_halo_beam1_inner_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0_Prescl", &L1Tech_BSC_halo_beam1_inner_v0_Prescl, &b_L1Tech_BSC_halo_beam1_inner_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0_5bx", &L1Tech_BSC_halo_beam1_inner_v0_5bx, &b_L1Tech_BSC_halo_beam1_inner_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0", &L1Tech_BSC_halo_beam1_outer_v0, &b_L1Tech_BSC_halo_beam1_outer_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0_Prescl", &L1Tech_BSC_halo_beam1_outer_v0_Prescl, &b_L1Tech_BSC_halo_beam1_outer_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0_5bx", &L1Tech_BSC_halo_beam1_outer_v0_5bx, &b_L1Tech_BSC_halo_beam1_outer_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0", &L1Tech_BSC_halo_beam2_inner_v0, &b_L1Tech_BSC_halo_beam2_inner_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0_Prescl", &L1Tech_BSC_halo_beam2_inner_v0_Prescl, &b_L1Tech_BSC_halo_beam2_inner_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0_5bx", &L1Tech_BSC_halo_beam2_inner_v0_5bx, &b_L1Tech_BSC_halo_beam2_inner_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0", &L1Tech_BSC_halo_beam2_outer_v0, &b_L1Tech_BSC_halo_beam2_outer_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0_Prescl", &L1Tech_BSC_halo_beam2_outer_v0_Prescl, &b_L1Tech_BSC_halo_beam2_outer_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0_5bx", &L1Tech_BSC_halo_beam2_outer_v0_5bx, &b_L1Tech_BSC_halo_beam2_outer_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0", &L1Tech_BSC_minBias_OR_v0, &b_L1Tech_BSC_minBias_OR_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0_Prescl", &L1Tech_BSC_minBias_OR_v0_Prescl, &b_L1Tech_BSC_minBias_OR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0_5bx", &L1Tech_BSC_minBias_OR_v0_5bx, &b_L1Tech_BSC_minBias_OR_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0", &L1Tech_BSC_minBias_inner_threshold1_v0, &b_L1Tech_BSC_minBias_inner_threshold1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0_Prescl", &L1Tech_BSC_minBias_inner_threshold1_v0_Prescl, &b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0_5bx", &L1Tech_BSC_minBias_inner_threshold1_v0_5bx, &b_L1Tech_BSC_minBias_inner_threshold1_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0", &L1Tech_BSC_minBias_inner_threshold2_v0, &b_L1Tech_BSC_minBias_inner_threshold2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0_Prescl", &L1Tech_BSC_minBias_inner_threshold2_v0_Prescl, &b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0_5bx", &L1Tech_BSC_minBias_inner_threshold2_v0_5bx, &b_L1Tech_BSC_minBias_inner_threshold2_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0", &L1Tech_BSC_minBias_threshold1_v0, &b_L1Tech_BSC_minBias_threshold1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0_Prescl", &L1Tech_BSC_minBias_threshold1_v0_Prescl, &b_L1Tech_BSC_minBias_threshold1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0_5bx", &L1Tech_BSC_minBias_threshold1_v0_5bx, &b_L1Tech_BSC_minBias_threshold1_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0", &L1Tech_BSC_minBias_threshold2_v0, &b_L1Tech_BSC_minBias_threshold2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0_Prescl", &L1Tech_BSC_minBias_threshold2_v0_Prescl, &b_L1Tech_BSC_minBias_threshold2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0_5bx", &L1Tech_BSC_minBias_threshold2_v0_5bx, &b_L1Tech_BSC_minBias_threshold2_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0", &L1Tech_BSC_splash_beam1_v0, &b_L1Tech_BSC_splash_beam1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0_Prescl", &L1Tech_BSC_splash_beam1_v0_Prescl, &b_L1Tech_BSC_splash_beam1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0_5bx", &L1Tech_BSC_splash_beam1_v0_5bx, &b_L1Tech_BSC_splash_beam1_v0_5bx);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0", &L1Tech_BSC_splash_beam2_v0, &b_L1Tech_BSC_splash_beam2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0_Prescl", &L1Tech_BSC_splash_beam2_v0_Prescl, &b_L1Tech_BSC_splash_beam2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0_5bx", &L1Tech_BSC_splash_beam2_v0_5bx, &b_L1Tech_BSC_splash_beam2_v0_5bx);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0", &L1Tech_CASTOR_HaloMuon_v0, &b_L1Tech_CASTOR_HaloMuon_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_Prescl", &L1Tech_CASTOR_HaloMuon_v0_Prescl, &b_L1Tech_CASTOR_HaloMuon_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_5bx", &L1Tech_CASTOR_HaloMuon_v0_5bx, &b_L1Tech_CASTOR_HaloMuon_v0_5bx);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0", &L1Tech_DT_GlobalOR_v0, &b_L1Tech_DT_GlobalOR_v0);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0_Prescl", &L1Tech_DT_GlobalOR_v0_Prescl, &b_L1Tech_DT_GlobalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0_5bx", &L1Tech_DT_GlobalOR_v0_5bx, &b_L1Tech_DT_GlobalOR_v0_5bx);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0", &L1Tech_HCAL_HBHE_totalOR_v0, &b_L1Tech_HCAL_HBHE_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_Prescl", &L1Tech_HCAL_HBHE_totalOR_v0_Prescl, &b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_5bx", &L1Tech_HCAL_HBHE_totalOR_v0_5bx, &b_L1Tech_HCAL_HBHE_totalOR_v0_5bx);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v0", &L1Tech_HCAL_HF_MMP_or_MPP_v0, &b_L1Tech_HCAL_HF_MMP_or_MPP_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v0_Prescl", &L1Tech_HCAL_HF_MMP_or_MPP_v0_Prescl, &b_L1Tech_HCAL_HF_MMP_or_MPP_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v0_5bx", &L1Tech_HCAL_HF_MMP_or_MPP_v0_5bx, &b_L1Tech_HCAL_HF_MMP_or_MPP_v0_5bx);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MM_or_PP_or_PM.v0", &L1Tech_HCAL_HF_MM_or_PP_or_PM_v0, &b_L1Tech_HCAL_HF_MM_or_PP_or_PM_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MM_or_PP_or_PM.v0_Prescl", &L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_Prescl, &b_L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MM_or_PP_or_PM.v0_5bx", &L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_5bx, &b_L1Tech_HCAL_HF_MM_or_PP_or_PM_v0_5bx);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v1", &L1Tech_HCAL_HF_coincidence_PM_v1, &b_L1Tech_HCAL_HF_coincidence_PM_v1);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v1_Prescl", &L1Tech_HCAL_HF_coincidence_PM_v1_Prescl, &b_L1Tech_HCAL_HF_coincidence_PM_v1_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v1_5bx", &L1Tech_HCAL_HF_coincidence_PM_v1_5bx, &b_L1Tech_HCAL_HF_coincidence_PM_v1_5bx);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0", &L1Tech_HCAL_HO_totalOR_v0, &b_L1Tech_HCAL_HO_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_Prescl", &L1Tech_HCAL_HO_totalOR_v0_Prescl, &b_L1Tech_HCAL_HO_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_5bx", &L1Tech_HCAL_HO_totalOR_v0_5bx, &b_L1Tech_HCAL_HO_totalOR_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0", &L1Tech_RPC_TTU_RB0_Cosmics_v0, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0_5bx", &L1Tech_RPC_TTU_RB0_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0_5bx", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0_5bx", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_5bx", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_5bx", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0", &L1Tech_RPC_TTU_RBst1_collisions_v0, &b_L1Tech_RPC_TTU_RBst1_collisions_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0_Prescl", &L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl, &b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0_5bx", &L1Tech_RPC_TTU_RBst1_collisions_v0_5bx, &b_L1Tech_RPC_TTU_RBst1_collisions_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0", &L1Tech_RPC_TTU_barrel_Cosmics_v0, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_5bx", &L1Tech_RPC_TTU_barrel_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0", &L1Tech_RPC_TTU_pointing_Cosmics_v0, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_5bx", &L1Tech_RPC_TTU_pointing_Cosmics_v0_5bx, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_loose_vertex.v0", &L1Tech_ZDC_Scint_loose_vertex_v0, &b_L1Tech_ZDC_Scint_loose_vertex_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_loose_vertex.v0_Prescl", &L1Tech_ZDC_Scint_loose_vertex_v0_Prescl, &b_L1Tech_ZDC_Scint_loose_vertex_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_loose_vertex.v0_5bx", &L1Tech_ZDC_Scint_loose_vertex_v0_5bx, &b_L1Tech_ZDC_Scint_loose_vertex_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_minus.v0", &L1Tech_ZDC_Scint_minus_v0, &b_L1Tech_ZDC_Scint_minus_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_minus.v0_Prescl", &L1Tech_ZDC_Scint_minus_v0_Prescl, &b_L1Tech_ZDC_Scint_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_minus.v0_5bx", &L1Tech_ZDC_Scint_minus_v0_5bx, &b_L1Tech_ZDC_Scint_minus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_plus.v0", &L1Tech_ZDC_Scint_plus_v0, &b_L1Tech_ZDC_Scint_plus_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_plus.v0_Prescl", &L1Tech_ZDC_Scint_plus_v0_Prescl, &b_L1Tech_ZDC_Scint_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_plus.v0_5bx", &L1Tech_ZDC_Scint_plus_v0_5bx, &b_L1Tech_ZDC_Scint_plus_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_tight_vertex.v0", &L1Tech_ZDC_Scint_tight_vertex_v0, &b_L1Tech_ZDC_Scint_tight_vertex_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_tight_vertex.v0_Prescl", &L1Tech_ZDC_Scint_tight_vertex_v0_Prescl, &b_L1Tech_ZDC_Scint_tight_vertex_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_Scint_tight_vertex.v0_5bx", &L1Tech_ZDC_Scint_tight_vertex_v0_5bx, &b_L1Tech_ZDC_Scint_tight_vertex_v0_5bx);
   Notify();
}

Bool_t checkOpenHLT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void checkOpenHLT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t checkOpenHLT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef checkOpenHLT_cxx
