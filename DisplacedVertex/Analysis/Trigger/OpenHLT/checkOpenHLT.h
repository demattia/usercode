//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  6 12:27:22 2011 by ROOT version 5.28/00b
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
   Int_t           NrecoJetCal;
   Int_t           NrecoJetGen;
   Int_t           NrecoTowCal;
   Float_t         recoJetCalPt[33];   //[NrecoJetCal]
   Float_t         recoJetCalPhi[33];   //[NrecoJetCal]
   Float_t         recoJetCalEta[33];   //[NrecoJetCal]
   Float_t         recoJetCalE[33];   //[NrecoJetCal]
   Float_t         recoJetCalEMF[33];   //[NrecoJetCal]
   Float_t         recoJetCalN90[33];   //[NrecoJetCal]
   Float_t         recoJetGenPt[1];   //[NrecoJetGen]
   Float_t         recoJetGenPhi[1];   //[NrecoJetGen]
   Float_t         recoJetGenEta[1];   //[NrecoJetGen]
   Float_t         recoJetGenE[1];   //[NrecoJetGen]
   Float_t         recoTowEt[563];   //[NrecoTowCal]
   Float_t         recoTowEta[563];   //[NrecoTowCal]
   Float_t         recoTowPhi[563];   //[NrecoTowCal]
   Float_t         recoTowE[563];   //[NrecoTowCal]
   Float_t         recoTowEm[563];   //[NrecoTowCal]
   Float_t         recoTowHad[563];   //[NrecoTowCal]
   Float_t         recoTowOE[563];   //[NrecoTowCal]
   Float_t         recoMetCal;
   Float_t         recoMetCalPhi;
   Float_t         recoMetCalSum;
   Float_t         recoMetGen;
   Float_t         recoMetGenPhi;
   Float_t         recoMetGenSum;
   Float_t         recoHTCal;
   Float_t         recoHTCalPhi;
   Float_t         recoHTCalSum;
   Int_t           NrecoJetCorCal;
   Float_t         recoJetCorCalPt[30];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalPhi[30];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalEta[30];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalE[30];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalEMF[30];   //[NrecoJetCorCal]
   Float_t         recoJetCorCalN90[30];   //[NrecoJetCorCal]
   Int_t           NohTau;
   Float_t         ohTauEta[22];   //[NohTau]
   Float_t         ohTauPhi[22];   //[NohTau]
   Float_t         ohTauPt[22];   //[NohTau]
   Float_t         ohTauEiso[22];   //[NohTau]
   Float_t         ohTauL25Tpt[22];   //[NohTau]
   Int_t           ohTauL3Tiso[22];   //[NohTau]
   Int_t           NohpfTau;
   Float_t         ohpfTauPt[6];   //[NohpfTau]
   Float_t         ohpfTauEta[6];   //[NohpfTau]
   Float_t         ohpfTauPhi[6];   //[NohpfTau]
   Float_t         ohpfTauLeadTrackPt[6];   //[NohpfTau]
   Float_t         ohpfTauLeadPionPt[6];   //[NohpfTau]
   Float_t         ohpfTauTrkIso[6];   //[NohpfTau]
   Float_t         ohpfTauGammaIso[6];   //[NohpfTau]
   Float_t         ohpfTauJetPt[6];   //[NohpfTau]
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
   Float_t         pfMHT;
   Int_t           NohPFJet;
   Float_t         pfJetPt[6];   //[NohPFJet]
   Float_t         pfJetEta[6];   //[NohPFJet]
   Float_t         pfJetPhi[6];   //[NohPFJet]
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
   Float_t         ohBJetIPL25Tag[10];   //[NohBJetL2]
   Float_t         ohBJetIPL3Tag[10];   //[NohBJetL2]
   Float_t         ohBJetIPLooseL25Tag[10];   //[NohBJetL2]
   Float_t         ohBJetIPLooseL3Tag[10];   //[NohBJetL2]
   Int_t           ohBJetMuL25Tag[10];   //[NohBJetL2]
   Float_t         ohBJetMuL3Tag[10];   //[NohBJetL2]
   Int_t           ohBJetPerfL25Tag[10];   //[NohBJetL2]
   Int_t           ohBJetPerfL3Tag[10];   //[NohBJetL2]
   Int_t           NrecoElec;
   Float_t         recoElecPt[1];   //[NrecoElec]
   Float_t         recoElecPhi[1];   //[NrecoElec]
   Float_t         recoElecEta[1];   //[NrecoElec]
   Float_t         recoElecEt[1];   //[NrecoElec]
   Float_t         recoElecE[1];   //[NrecoElec]
   Int_t           recoElecEleID[1];   //[NrecoElec]
   Int_t           NrecoPhot;
   Float_t         recoPhotPt[1];   //[NrecoPhot]
   Float_t         recoPhotPhi[1];   //[NrecoPhot]
   Float_t         recoPhotEta[1];   //[NrecoPhot]
   Float_t         recoPhotEt[1];   //[NrecoPhot]
   Float_t         recoPhotE[1];   //[NrecoPhot]
   Int_t           NohPhot;
   Float_t         ohPhotEt[6];   //[NohPhot]
   Float_t         ohPhotEta[6];   //[NohPhot]
   Float_t         ohPhotPhi[6];   //[NohPhot]
   Float_t         ohPhotEiso[6];   //[NohPhot]
   Float_t         ohPhotHiso[6];   //[NohPhot]
   Float_t         ohPhotTiso[6];   //[NohPhot]
   Int_t           ohPhotL1iso[6];   //[NohPhot]
   Float_t         ohPhotClusShap[6];   //[NohPhot]
   Float_t         ohPhotR9[6];   //[NohPhot]
   Float_t         ohPhotHforHoverE[6];   //[NohPhot]
   Int_t           NohEle;
   Float_t         ohEleEt[9];   //[NohEle]
   Float_t         ohEleEta[9];   //[NohEle]
   Float_t         ohElePhi[9];   //[NohEle]
   Float_t         ohEleE[9];   //[NohEle]
   Float_t         ohEleP[9];   //[NohEle]
   Float_t         ohEleHiso[9];   //[NohEle]
   Float_t         ohEleTiso[9];   //[NohEle]
   Float_t         ohEleEiso[9];   //[NohEle]
   Int_t           ohEleL1iso[9];   //[NohEle]
   Int_t           ohElePixelSeeds[9];   //[NohEle]
   Int_t           ohEleNewSC[9];   //[NohEle]
   Float_t         ohEleClusShap[9];   //[NohEle]
   Float_t         ohEleDeta[9];   //[NohEle]
   Float_t         ohEleDphi[9];   //[NohEle]
   Float_t         ohEleR9[9];   //[NohEle]
   Float_t         ohEleHforHoverE[9];   //[NohEle]
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
   Float_t         ohMuL2Pt[4];   //[NohMuL2]
   Float_t         ohMuL2Phi[4];   //[NohMuL2]
   Float_t         ohMuL2Eta[4];   //[NohMuL2]
   Int_t           ohMuL2Chg[4];   //[NohMuL2]
   Float_t         ohMuL2PtErr[4];   //[NohMuL2]
   Int_t           ohMuL2Iso[4];   //[NohMuL2]
   Float_t         ohMuL2Dr[4];   //[NohMuL2]
   Float_t         ohMuL2Dz[4];   //[NohMuL2]
   Int_t           ohMuL2L1idx[4];   //[NohMuL2]
   Int_t           NohMuL3;
   Float_t         ohMuL3Pt[3];   //[NohMuL3]
   Float_t         ohMuL3Phi[3];   //[NohMuL3]
   Float_t         ohMuL3Eta[3];   //[NohMuL3]
   Int_t           ohMuL3Chg[3];   //[NohMuL3]
   Float_t         ohMuL3PtErr[3];   //[NohMuL3]
   Int_t           ohMuL3Iso[3];   //[NohMuL3]
   Float_t         ohMuL3Dr[3];   //[NohMuL3]
   Float_t         ohMuL3Dz[3];   //[NohMuL3]
   Int_t           ohMuL3L2idx[3];   //[NohMuL3]
   Int_t           NohOniaPixel;
   Float_t         ohOniaPixelPt[31];   //[NohOniaPixel]
   Float_t         ohOniaPixelPhi[31];   //[NohOniaPixel]
   Float_t         ohOniaPixelEta[31];   //[NohOniaPixel]
   Int_t           ohOniaPixelChg[31];   //[NohOniaPixel]
   Float_t         ohOniaPixelDr[31];   //[NohOniaPixel]
   Float_t         ohOniaPixelDz[31];   //[NohOniaPixel]
   Int_t           ohOniaPixelHits[31];   //[NohOniaPixel]
   Float_t         ohOniaPixelNormChi2[31];   //[NohOniaPixel]
   Int_t           NohOniaTrack;
   Float_t         ohOniaTrackPt[10];   //[NohOniaTrack]
   Float_t         ohOniaTrackPhi[10];   //[NohOniaTrack]
   Float_t         ohOniaTrackEta[10];   //[NohOniaTrack]
   Int_t           ohOniaTrackChg[10];   //[NohOniaTrack]
   Float_t         ohOniaTrackDr[10];   //[NohOniaTrack]
   Float_t         ohOniaTrackDz[10];   //[NohOniaTrack]
   Int_t           ohOniaTrackHits[10];   //[NohOniaTrack]
   Float_t         ohOniaTrackNormChi2[10];   //[NohOniaTrack]
   Int_t           NohMuL2NoVtx;
   Float_t         ohMuL2NoVtxPt[4];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxPhi[4];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxEta[4];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxChg[4];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxPtErr[4];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxDr[4];   //[NohMuL2NoVtx]
   Float_t         ohMuL2NoVtxDz[4];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxL1idx[4];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxNhits[4];   //[NohMuL2NoVtx]
   Int_t           ohMuL2NoVtxNchambers[4];   //[NohMuL2NoVtx]
   Float_t         ohHighestEnergyEERecHit;
   Float_t         ohHighestEnergyEBRecHit;
   Float_t         ohHighestEnergyHBHERecHit;
   Float_t         ohHighestEnergyHORecHit;
   Float_t         ohHighestEnergyHFRecHit;
   Int_t           Nalcapi0clusters;
   Float_t         ohAlcapi0ptClusAll[43];   //[Nalcapi0clusters]
   Float_t         ohAlcapi0etaClusAll[43];   //[Nalcapi0clusters]
   Float_t         ohAlcapi0phiClusAll[43];   //[Nalcapi0clusters]
   Float_t         ohAlcapi0s4s9ClusAll[43];   //[Nalcapi0clusters]
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
   Float_t         ohPixelTracksL3Pt[932];   //[NohPixelTracksL3]
   Float_t         ohPixelTracksL3Eta[932];   //[NohPixelTracksL3]
   Float_t         ohPixelTracksL3Phi[932];   //[NohPixelTracksL3]
   Float_t         ohPixelTracksL3Vz[932];   //[NohPixelTracksL3]
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
   Int_t           recoNVrt;
   Float_t         recoVrtX[11];   //[NVrtx]
   Float_t         recoVrtY[11];   //[NVrtx]
   Float_t         recoVrtZ[11];   //[NVrtx]
   Int_t           recoVrtNtrk[11];   //[NVrtx]
   Float_t         recoVrtChi2[11];   //[NVrtx]
   Float_t         recoVrtNdof[11];   //[NVrtx]
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiBlock;
   Int_t           Bx;
   Int_t           Orbit;
   Int_t           HLT_Activity_CSC;
   Int_t           HLT_Activity_CSC_Prescl;
   Int_t           HLT_Activity_DT;
   Int_t           HLT_Activity_DT_Prescl;
   Int_t           HLT_Activity_DT_Tuned;
   Int_t           HLT_Activity_DT_Tuned_Prescl;
   Int_t           HLT_Activity_Ecal_SC7;
   Int_t           HLT_Activity_Ecal_SC7_Prescl;
   Int_t           HLT_Activity_Ecal_SC17;
   Int_t           HLT_Activity_Ecal_SC17_Prescl;
   Int_t           HLT_L1Jet6U;
   Int_t           HLT_L1Jet6U_Prescl;
   Int_t           HLT_L1Jet10U;
   Int_t           HLT_L1Jet10U_Prescl;
   Int_t           HLT_Jet15U;
   Int_t           HLT_Jet15U_Prescl;
   Int_t           HLT_Jet15U_HcalNoiseFiltered;
   Int_t           HLT_Jet15U_HcalNoiseFiltered_Prescl;
   Int_t           HLT_Jet30U;
   Int_t           HLT_Jet30U_Prescl;
   Int_t           HLT_Jet50U;
   Int_t           HLT_Jet50U_Prescl;
   Int_t           HLT_Jet70U_v2;
   Int_t           HLT_Jet70U_v2_Prescl;
   Int_t           HLT_Jet100U_v2;
   Int_t           HLT_Jet100U_v2_Prescl;
   Int_t           HLT_Jet140U_v1;
   Int_t           HLT_Jet140U_v1_Prescl;
   Int_t           HLT_DiJetAve15U;
   Int_t           HLT_DiJetAve15U_Prescl;
   Int_t           HLT_DiJetAve30U;
   Int_t           HLT_DiJetAve30U_Prescl;
   Int_t           HLT_DiJetAve50U;
   Int_t           HLT_DiJetAve50U_Prescl;
   Int_t           HLT_DiJetAve70U_v2;
   Int_t           HLT_DiJetAve70U_v2_Prescl;
   Int_t           HLT_DiJetAve100U_v1;
   Int_t           HLT_DiJetAve100U_v1_Prescl;
   Int_t           HLT_DoubleJet15U_ForwardBackward;
   Int_t           HLT_DoubleJet15U_ForwardBackward_Prescl;
   Int_t           HLT_DoubleJet25U_ForwardBackward;
   Int_t           HLT_DoubleJet25U_ForwardBackward_Prescl;
   Int_t           HLT_ExclDiJet30U_HFAND_v1;
   Int_t           HLT_ExclDiJet30U_HFAND_v1_Prescl;
   Int_t           HLT_ExclDiJet30U_HFOR_v1;
   Int_t           HLT_ExclDiJet30U_HFOR_v1_Prescl;
   Int_t           HLT_QuadJet15U_v2;
   Int_t           HLT_QuadJet15U_v2_Prescl;
   Int_t           HLT_QuadJet20U_v2;
   Int_t           HLT_QuadJet20U_v2_Prescl;
   Int_t           HLT_QuadJet25U_v2;
   Int_t           HLT_QuadJet25U_v2_Prescl;
   Int_t           HLT_L1ETT100;
   Int_t           HLT_L1ETT100_Prescl;
   Int_t           HLT_L1ETT140_v1;
   Int_t           HLT_L1ETT140_v1_Prescl;
   Int_t           HLT_EcalOnly_SumEt160_v2;
   Int_t           HLT_EcalOnly_SumEt160_v2_Prescl;
   Int_t           HLT_L1MET20;
   Int_t           HLT_L1MET20_Prescl;
   Int_t           HLT_MET45;
   Int_t           HLT_MET45_Prescl;
   Int_t           HLT_MET45_HT100U_v1;
   Int_t           HLT_MET45_HT100U_v1_Prescl;
   Int_t           HLT_MET45_HT120U_v1;
   Int_t           HLT_MET45_HT120U_v1_Prescl;
   Int_t           HLT_MET65;
   Int_t           HLT_MET65_Prescl;
   Int_t           HLT_MET80_v1;
   Int_t           HLT_MET80_v1_Prescl;
   Int_t           HLT_MET100_v2;
   Int_t           HLT_MET100_v2_Prescl;
   Int_t           HLT_HT50U_v1;
   Int_t           HLT_HT50U_v1_Prescl;
   Int_t           HLT_HT100U;
   Int_t           HLT_HT100U_Prescl;
   Int_t           HLT_HT120U;
   Int_t           HLT_HT120U_Prescl;
   Int_t           HLT_HT140U;
   Int_t           HLT_HT140U_Prescl;
   Int_t           HLT_HT140U_Eta3_v1;
   Int_t           HLT_HT140U_Eta3_v1_Prescl;
   Int_t           HLT_HT160U_v1;
   Int_t           HLT_HT160U_v1_Prescl;
   Int_t           HLT_HT200U_v1;
   Int_t           HLT_HT200U_v1_Prescl;
   Int_t           HLT_L1MuOpen;
   Int_t           HLT_L1MuOpen_Prescl;
   Int_t           HLT_L1MuOpen_DT;
   Int_t           HLT_L1MuOpen_DT_Prescl;
   Int_t           HLT_L1MuOpen_AntiBPTX;
   Int_t           HLT_L1MuOpen_AntiBPTX_Prescl;
   Int_t           HLT_L1Mu7_v1;
   Int_t           HLT_L1Mu7_v1_Prescl;
   Int_t           HLT_L1Mu20;
   Int_t           HLT_L1Mu20_Prescl;
   Int_t           HLT_L2Mu0_NoVertex;
   Int_t           HLT_L2Mu0_NoVertex_Prescl;
   Int_t           HLT_L2Mu7_v1;
   Int_t           HLT_L2Mu7_v1_Prescl;
   Int_t           HLT_L2Mu30_v1;
   Int_t           HLT_L2Mu30_v1_Prescl;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu3_Prescl;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu5_Prescl;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu7_Prescl;
   Int_t           HLT_Mu9;
   Int_t           HLT_Mu9_Prescl;
   Int_t           HLT_Mu11;
   Int_t           HLT_Mu11_Prescl;
   Int_t           HLT_Mu13_v1;
   Int_t           HLT_Mu13_v1_Prescl;
   Int_t           HLT_Mu15_v1;
   Int_t           HLT_Mu15_v1_Prescl;
   Int_t           HLT_IsoMu9;
   Int_t           HLT_IsoMu9_Prescl;
   Int_t           HLT_IsoMu11_v1;
   Int_t           HLT_IsoMu11_v1_Prescl;
   Int_t           HLT_Mu20_NoVertex;
   Int_t           HLT_Mu20_NoVertex_Prescl;
   Int_t           HLT_L1DoubleMuOpen;
   Int_t           HLT_L1DoubleMuOpen_Prescl;
   Int_t           HLT_L2DoubleMu0;
   Int_t           HLT_L2DoubleMu0_Prescl;
   Int_t           HLT_L2DoubleMu20_NoVertex_v1;
   Int_t           HLT_L2DoubleMu20_NoVertex_v1_Prescl;
   Int_t           HLT_DoubleMu0_Quarkonium_v1;
   Int_t           HLT_DoubleMu0_Quarkonium_v1_Prescl;
   Int_t           HLT_DoubleMu0_Quarkonium_LS_v1;
   Int_t           HLT_DoubleMu0_Quarkonium_LS_v1_Prescl;
   Int_t           HLT_DoubleMu0;
   Int_t           HLT_DoubleMu0_Prescl;
   Int_t           HLT_DoubleMu3_v2;
   Int_t           HLT_DoubleMu3_v2_Prescl;
   Int_t           HLT_DoubleMu5_v1;
   Int_t           HLT_DoubleMu5_v1_Prescl;
   Int_t           HLT_Mu5_L2Mu0;
   Int_t           HLT_Mu5_L2Mu0_Prescl;
   Int_t           HLT_Mu3_Track3_Jpsi;
   Int_t           HLT_Mu3_Track3_Jpsi_Prescl;
   Int_t           HLT_Mu3_Track5_Jpsi_v1;
   Int_t           HLT_Mu3_Track5_Jpsi_v1_Prescl;
   Int_t           HLT_Mu5_Track0_Jpsi;
   Int_t           HLT_Mu5_Track0_Jpsi_Prescl;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Prescl;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1_Prescl;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Prescl;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Prescl;
   Int_t           HLT_L1SingleEG2;
   Int_t           HLT_L1SingleEG2_Prescl;
   Int_t           HLT_L1SingleEG8;
   Int_t           HLT_L1SingleEG8_Prescl;
   Int_t           HLT_Ele10_SW_L1R;
   Int_t           HLT_Ele10_SW_L1R_Prescl;
   Int_t           HLT_Ele12_SW_TightEleId_L1R;
   Int_t           HLT_Ele12_SW_TightEleId_L1R_Prescl;
   Int_t           HLT_Ele12_SW_TighterEleId_L1R_v1;
   Int_t           HLT_Ele12_SW_TighterEleId_L1R_v1_Prescl;
   Int_t           HLT_Ele12_SW_TighterEleIdIsol_L1R_v1;
   Int_t           HLT_Ele12_SW_TighterEleIdIsol_L1R_v1_Prescl;
   Int_t           HLT_Ele17_SW_L1R;
   Int_t           HLT_Ele17_SW_L1R_Prescl;
   Int_t           HLT_Ele17_SW_TightEleId_L1R;
   Int_t           HLT_Ele17_SW_TightEleId_L1R_Prescl;
   Int_t           HLT_Ele17_SW_TighterEleId_L1R_v1;
   Int_t           HLT_Ele17_SW_TighterEleId_L1R_v1_Prescl;
   Int_t           HLT_Ele17_SW_TightEleIdIsol_L1R_v1;
   Int_t           HLT_Ele17_SW_TightEleIdIsol_L1R_v1_Prescl;
   Int_t           HLT_Ele17_SW_TighterEleIdIsol_L1R_v1;
   Int_t           HLT_Ele17_SW_TighterEleIdIsol_L1R_v1_Prescl;
   Int_t           HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1;
   Int_t           HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1_Prescl;
   Int_t           HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1;
   Int_t           HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1_Prescl;
   Int_t           HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1;
   Int_t           HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1_Prescl;
   Int_t           HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1;
   Int_t           HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1_Prescl;
   Int_t           HLT_DoubleEle4_SW_eeRes_L1R;
   Int_t           HLT_DoubleEle4_SW_eeRes_L1R_Prescl;
   Int_t           HLT_DoubleEle15_SW_L1R_v1;
   Int_t           HLT_DoubleEle15_SW_L1R_v1_Prescl;
   Int_t           HLT_Photon10_Cleaned_L1R;
   Int_t           HLT_Photon10_Cleaned_L1R_Prescl;
   Int_t           HLT_Photon15_Cleaned_L1R;
   Int_t           HLT_Photon15_Cleaned_L1R_Prescl;
   Int_t           HLT_Photon17_SC17HE_L1R_v1;
   Int_t           HLT_Photon17_SC17HE_L1R_v1_Prescl;
   Int_t           HLT_Photon20_NoHE_L1R;
   Int_t           HLT_Photon20_NoHE_L1R_Prescl;
   Int_t           HLT_Photon20_Cleaned_L1R;
   Int_t           HLT_Photon20_Cleaned_L1R_Prescl;
   Int_t           HLT_Photon30_Cleaned_L1R;
   Int_t           HLT_Photon30_Cleaned_L1R_Prescl;
   Int_t           HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1;
   Int_t           HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_Prescl;
   Int_t           HLT_Photon35_Isol_Cleaned_L1R_v1;
   Int_t           HLT_Photon35_Isol_Cleaned_L1R_v1_Prescl;
   Int_t           HLT_Photon50_Cleaned_L1R_v1;
   Int_t           HLT_Photon50_Cleaned_L1R_v1_Prescl;
   Int_t           HLT_Photon50_NoHE_L1R;
   Int_t           HLT_Photon50_NoHE_L1R_Prescl;
   Int_t           HLT_Photon70_NoHE_Cleaned_L1R_v1;
   Int_t           HLT_Photon70_NoHE_Cleaned_L1R_v1_Prescl;
   Int_t           HLT_Photon100_NoHE_Cleaned_L1R_v1;
   Int_t           HLT_Photon100_NoHE_Cleaned_L1R_v1_Prescl;
   Int_t           HLT_DoublePhoton5_CEP_L1R;
   Int_t           HLT_DoublePhoton5_CEP_L1R_Prescl;
   Int_t           HLT_DoublePhoton17_L1R;
   Int_t           HLT_DoublePhoton17_L1R_Prescl;
   Int_t           HLT_SingleIsoTau20_Trk5_MET20;
   Int_t           HLT_SingleIsoTau20_Trk5_MET20_Prescl;
   Int_t           HLT_SingleIsoTau20_Trk15_MET20;
   Int_t           HLT_SingleIsoTau20_Trk15_MET20_Prescl;
   Int_t           HLT_SingleIsoTau30_Trk5_MET20;
   Int_t           HLT_SingleIsoTau30_Trk5_MET20_Prescl;
   Int_t           HLT_SingleIsoTau30_Trk5_v2;
   Int_t           HLT_SingleIsoTau30_Trk5_v2_Prescl;
   Int_t           HLT_DoubleIsoTau15_OneLeg_Trk5;
   Int_t           HLT_DoubleIsoTau15_OneLeg_Trk5_Prescl;
   Int_t           HLT_DoubleIsoTau15_Trk5;
   Int_t           HLT_DoubleIsoTau15_Trk5_Prescl;
   Int_t           HLT_BTagMu_DiJet10U_v1;
   Int_t           HLT_BTagMu_DiJet10U_v1_Prescl;
   Int_t           HLT_BTagMu_DiJet20U_v1;
   Int_t           HLT_BTagMu_DiJet20U_v1_Prescl;
   Int_t           HLT_BTagMu_DiJet20U_Mu5_v1;
   Int_t           HLT_BTagMu_DiJet20U_Mu5_v1_Prescl;
   Int_t           HLT_StoppedHSCP20_v3;
   Int_t           HLT_StoppedHSCP20_v3_Prescl;
   Int_t           HLT_StoppedHSCP35_v3;
   Int_t           HLT_StoppedHSCP35_v3_Prescl;
   Int_t           HLT_Mu5_Photon11_Cleaned_L1R_v1;
   Int_t           HLT_Mu5_Photon11_Cleaned_L1R_v1_Prescl;
   Int_t           HLT_Mu5_Ele5_v1;
   Int_t           HLT_Mu5_Ele5_v1_Prescl;
   Int_t           HLT_Mu5_Ele9_v1;
   Int_t           HLT_Mu5_Ele9_v1_Prescl;
   Int_t           HLT_Mu5_Jet35U_v1;
   Int_t           HLT_Mu5_Jet35U_v1_Prescl;
   Int_t           HLT_Mu5_Jet50U_v2;
   Int_t           HLT_Mu5_Jet50U_v2_Prescl;
   Int_t           HLT_Mu5_MET45_v1;
   Int_t           HLT_Mu5_MET45_v1_Prescl;
   Int_t           HLT_Mu5_HT50U_v1;
   Int_t           HLT_Mu5_HT50U_v1_Prescl;
   Int_t           HLT_Mu5_HT70U_v1;
   Int_t           HLT_Mu5_HT70U_v1_Prescl;
   Int_t           HLT_Ele10_MET45_v1;
   Int_t           HLT_Ele10_MET45_v1_Prescl;
   Int_t           HLT_ZeroBias;
   Int_t           HLT_ZeroBias_Prescl;
   Int_t           HLT_ZeroBiasPixel_SingleTrack;
   Int_t           HLT_ZeroBiasPixel_SingleTrack_Prescl;
   Int_t           HLT_MinBiasPixel_SingleTrack;
   Int_t           HLT_MinBiasPixel_SingleTrack_Prescl;
   Int_t           HLT_MultiVertex6;
   Int_t           HLT_MultiVertex6_Prescl;
   Int_t           HLT_MultiVertex8_L1ETT60;
   Int_t           HLT_MultiVertex8_L1ETT60_Prescl;
   Int_t           HLT_L1_BptxXOR_BscMinBiasOR;
   Int_t           HLT_L1_BptxXOR_BscMinBiasOR_Prescl;
   Int_t           HLT_L1Tech_BSC_minBias_OR;
   Int_t           HLT_L1Tech_BSC_minBias_OR_Prescl;
   Int_t           HLT_L1Tech_BSC_minBias;
   Int_t           HLT_L1Tech_BSC_minBias_Prescl;
   Int_t           HLT_L1Tech_BSC_halo;
   Int_t           HLT_L1Tech_BSC_halo_Prescl;
   Int_t           HLT_L1Tech_BSC_halo_forPhysicsBackground;
   Int_t           HLT_L1Tech_BSC_halo_forPhysicsBackground_Prescl;
   Int_t           HLT_L1Tech_BSC_HighMultiplicity;
   Int_t           HLT_L1Tech_BSC_HighMultiplicity_Prescl;
   Int_t           HLT_L1Tech_RPC_TTU_RBst1_collisions;
   Int_t           HLT_L1Tech_RPC_TTU_RBst1_collisions_Prescl;
   Int_t           HLT_L1Tech_HCAL_HF;
   Int_t           HLT_L1Tech_HCAL_HF_Prescl;
   Int_t           HLT_TrackerCosmics;
   Int_t           HLT_TrackerCosmics_Prescl;
   Int_t           HLT_IsoTrackHB_v2;
   Int_t           HLT_IsoTrackHB_v2_Prescl;
   Int_t           HLT_IsoTrackHE_v2;
   Int_t           HLT_IsoTrackHE_v2_Prescl;
   Int_t           HLT_RPCBarrelCosmics;
   Int_t           HLT_RPCBarrelCosmics_Prescl;
   Int_t           HLT_HcalPhiSym;
   Int_t           HLT_HcalPhiSym_Prescl;
   Int_t           HLT_HcalNZS;
   Int_t           HLT_HcalNZS_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity70;
   Int_t           HLT_PixelTracks_Multiplicity70_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity85;
   Int_t           HLT_PixelTracks_Multiplicity85_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity100;
   Int_t           HLT_PixelTracks_Multiplicity100_Prescl;
   Int_t           HLT_GlobalRunHPDNoise;
   Int_t           HLT_GlobalRunHPDNoise_Prescl;
   Int_t           HLT_TechTrigHCALNoise;
   Int_t           HLT_TechTrigHCALNoise_Prescl;
   Int_t           HLT_L1_BPTX;
   Int_t           HLT_L1_BPTX_Prescl;
   Int_t           HLT_L1_BPTX_MinusOnly;
   Int_t           HLT_L1_BPTX_MinusOnly_Prescl;
   Int_t           HLT_L1_BPTX_PlusOnly;
   Int_t           HLT_L1_BPTX_PlusOnly_Prescl;
   Int_t           HLT_DTErrors;
   Int_t           HLT_DTErrors_Prescl;
   Int_t           HLT_LogMonitor;
   Int_t           HLT_LogMonitor_Prescl;
   Int_t           HLT_Calibration;
   Int_t           HLT_Calibration_Prescl;
   Int_t           HLT_EcalCalibration;
   Int_t           HLT_EcalCalibration_Prescl;
   Int_t           HLT_HcalCalibration;
   Int_t           HLT_HcalCalibration_Prescl;
   Int_t           HLT_Random;
   Int_t           HLT_Random_Prescl;
   Int_t           AlCa_EcalPhiSym;
   Int_t           AlCa_EcalPhiSym_Prescl;
   Int_t           AlCa_EcalPi0;
   Int_t           AlCa_EcalPi0_Prescl;
   Int_t           AlCa_EcalEta;
   Int_t           AlCa_EcalEta_Prescl;
   Int_t           AlCa_RPCMuonNoHits;
   Int_t           AlCa_RPCMuonNoHits_Prescl;
   Int_t           AlCa_RPCMuonNoTriggers;
   Int_t           AlCa_RPCMuonNoTriggers_Prescl;
   Int_t           AlCa_RPCMuonNormalisation;
   Int_t           AlCa_RPCMuonNormalisation_Prescl;
   Int_t           DQM_FEDIntegrity;
   Int_t           DQM_FEDIntegrity_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           L1_BptxMinus;
   Int_t           L1_BptxMinus_Prescl;
   Int_t           L1_BptxMinus_5bx;
   Int_t           L1_BptxMinus_NotBptxPlus;
   Int_t           L1_BptxMinus_NotBptxPlus_Prescl;
   Int_t           L1_BptxMinus_NotBptxPlus_5bx;
   Int_t           L1_BptxPlus;
   Int_t           L1_BptxPlus_Prescl;
   Int_t           L1_BptxPlus_5bx;
   Int_t           L1_BptxPlusORMinus;
   Int_t           L1_BptxPlusORMinus_Prescl;
   Int_t           L1_BptxPlusORMinus_5bx;
   Int_t           L1_BptxPlus_NotBptxMinus;
   Int_t           L1_BptxPlus_NotBptxMinus_Prescl;
   Int_t           L1_BptxPlus_NotBptxMinus_5bx;
   Int_t           L1_BptxXOR_BscMinBiasOR;
   Int_t           L1_BptxXOR_BscMinBiasOR_Prescl;
   Int_t           L1_BptxXOR_BscMinBiasOR_5bx;
   Int_t           L1_Bsc2Minus_BptxMinus;
   Int_t           L1_Bsc2Minus_BptxMinus_Prescl;
   Int_t           L1_Bsc2Minus_BptxMinus_5bx;
   Int_t           L1_Bsc2Plus_BptxPlus;
   Int_t           L1_Bsc2Plus_BptxPlus_Prescl;
   Int_t           L1_Bsc2Plus_BptxPlus_5bx;
   Int_t           L1_BscHaloBeam1Inner;
   Int_t           L1_BscHaloBeam1Inner_Prescl;
   Int_t           L1_BscHaloBeam1Inner_5bx;
   Int_t           L1_BscHaloBeam1Outer;
   Int_t           L1_BscHaloBeam1Outer_Prescl;
   Int_t           L1_BscHaloBeam1Outer_5bx;
   Int_t           L1_BscHaloBeam2Inner;
   Int_t           L1_BscHaloBeam2Inner_Prescl;
   Int_t           L1_BscHaloBeam2Inner_5bx;
   Int_t           L1_BscHaloBeam2Outer;
   Int_t           L1_BscHaloBeam2Outer_Prescl;
   Int_t           L1_BscHaloBeam2Outer_5bx;
   Int_t           L1_BscHighMultiplicity;
   Int_t           L1_BscHighMultiplicity_Prescl;
   Int_t           L1_BscHighMultiplicity_5bx;
   Int_t           L1_BscMinBiasInnerThreshold1;
   Int_t           L1_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold1_5bx;
   Int_t           L1_BscMinBiasInnerThreshold2;
   Int_t           L1_BscMinBiasInnerThreshold2_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold2_5bx;
   Int_t           L1_BscMinBiasOR;
   Int_t           L1_BscMinBiasOR_Prescl;
   Int_t           L1_BscMinBiasOR_5bx;
   Int_t           L1_BscMinBiasOR_BptxPlusANDMinus;
   Int_t           L1_BscMinBiasOR_BptxPlusANDMinus_Prescl;
   Int_t           L1_BscMinBiasOR_BptxPlusANDMinus_5bx;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus_Prescl;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus_5bx;
   Int_t           L1_BscMinBiasThreshold1;
   Int_t           L1_BscMinBiasThreshold1_Prescl;
   Int_t           L1_BscMinBiasThreshold1_5bx;
   Int_t           L1_BscMinBiasThreshold2;
   Int_t           L1_BscMinBiasThreshold2_Prescl;
   Int_t           L1_BscMinBiasThreshold2_5bx;
   Int_t           L1_BscSplashBeam1;
   Int_t           L1_BscSplashBeam1_Prescl;
   Int_t           L1_BscSplashBeam1_5bx;
   Int_t           L1_BscSplashBeam2;
   Int_t           L1_BscSplashBeam2_Prescl;
   Int_t           L1_BscSplashBeam2_5bx;
   Int_t           L1_DoubleEG05_TopBottom;
   Int_t           L1_DoubleEG05_TopBottom_Prescl;
   Int_t           L1_DoubleEG05_TopBottom_5bx;
   Int_t           L1_DoubleEG2;
   Int_t           L1_DoubleEG2_Prescl;
   Int_t           L1_DoubleEG2_5bx;
   Int_t           L1_DoubleEG5;
   Int_t           L1_DoubleEG5_Prescl;
   Int_t           L1_DoubleEG5_5bx;
   Int_t           L1_DoubleForJet10_EtaOpp;
   Int_t           L1_DoubleForJet10_EtaOpp_Prescl;
   Int_t           L1_DoubleForJet10_EtaOpp_5bx;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1_Prescl;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1_5bx;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1_Prescl;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1_5bx;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200_5bx;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4_5bx;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200_5bx;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4_5bx;
   Int_t           L1_DoubleJet30;
   Int_t           L1_DoubleJet30_Prescl;
   Int_t           L1_DoubleJet30_5bx;
   Int_t           L1_DoubleMu3;
   Int_t           L1_DoubleMu3_Prescl;
   Int_t           L1_DoubleMu3_5bx;
   Int_t           L1_DoubleMuOpen;
   Int_t           L1_DoubleMuOpen_Prescl;
   Int_t           L1_DoubleMuOpen_5bx;
   Int_t           L1_DoubleMuTopBottom;
   Int_t           L1_DoubleMuTopBottom_Prescl;
   Int_t           L1_DoubleMuTopBottom_5bx;
   Int_t           L1_DoubleTauJet14;
   Int_t           L1_DoubleTauJet14_Prescl;
   Int_t           L1_DoubleTauJet14_5bx;
   Int_t           L1_ETM12;
   Int_t           L1_ETM12_Prescl;
   Int_t           L1_ETM12_5bx;
   Int_t           L1_ETM20;
   Int_t           L1_ETM20_Prescl;
   Int_t           L1_ETM20_5bx;
   Int_t           L1_ETM30;
   Int_t           L1_ETM30_Prescl;
   Int_t           L1_ETM30_5bx;
   Int_t           L1_ETM70;
   Int_t           L1_ETM70_Prescl;
   Int_t           L1_ETM70_5bx;
   Int_t           L1_ETT100;
   Int_t           L1_ETT100_Prescl;
   Int_t           L1_ETT100_5bx;
   Int_t           L1_ETT140;
   Int_t           L1_ETT140_Prescl;
   Int_t           L1_ETT140_5bx;
   Int_t           L1_ETT30;
   Int_t           L1_ETT30_Prescl;
   Int_t           L1_ETT30_5bx;
   Int_t           L1_ETT60;
   Int_t           L1_ETT60_Prescl;
   Int_t           L1_ETT60_5bx;
   Int_t           L1_HTM20;
   Int_t           L1_HTM20_Prescl;
   Int_t           L1_HTM20_5bx;
   Int_t           L1_HTM30;
   Int_t           L1_HTM30_Prescl;
   Int_t           L1_HTM30_5bx;
   Int_t           L1_HTT100;
   Int_t           L1_HTT100_Prescl;
   Int_t           L1_HTT100_5bx;
   Int_t           L1_HTT200;
   Int_t           L1_HTT200_Prescl;
   Int_t           L1_HTT200_5bx;
   Int_t           L1_HTT50;
   Int_t           L1_HTT50_Prescl;
   Int_t           L1_HTT50_5bx;
   Int_t           L1_IsoEG10_Jet6_ForJet6;
   Int_t           L1_IsoEG10_Jet6_ForJet6_Prescl;
   Int_t           L1_IsoEG10_Jet6_ForJet6_5bx;
   Int_t           L1_Mu3_EG5;
   Int_t           L1_Mu3_EG5_Prescl;
   Int_t           L1_Mu3_EG5_5bx;
   Int_t           L1_Mu3_Jet10;
   Int_t           L1_Mu3_Jet10_Prescl;
   Int_t           L1_Mu3_Jet10_5bx;
   Int_t           L1_Mu3_Jet6;
   Int_t           L1_Mu3_Jet6_Prescl;
   Int_t           L1_Mu3_Jet6_5bx;
   Int_t           L1_Mu5_Jet6;
   Int_t           L1_Mu5_Jet6_Prescl;
   Int_t           L1_Mu5_Jet6_5bx;
   Int_t           L1_QuadJet6;
   Int_t           L1_QuadJet6_Prescl;
   Int_t           L1_QuadJet6_5bx;
   Int_t           L1_QuadJet8;
   Int_t           L1_QuadJet8_Prescl;
   Int_t           L1_QuadJet8_5bx;
   Int_t           L1_SingleCenJet2;
   Int_t           L1_SingleCenJet2_Prescl;
   Int_t           L1_SingleCenJet2_5bx;
   Int_t           L1_SingleCenJet4;
   Int_t           L1_SingleCenJet4_Prescl;
   Int_t           L1_SingleCenJet4_5bx;
   Int_t           L1_SingleEG1;
   Int_t           L1_SingleEG1_Prescl;
   Int_t           L1_SingleEG1_5bx;
   Int_t           L1_SingleEG10;
   Int_t           L1_SingleEG10_Prescl;
   Int_t           L1_SingleEG10_5bx;
   Int_t           L1_SingleEG12;
   Int_t           L1_SingleEG12_Prescl;
   Int_t           L1_SingleEG12_5bx;
   Int_t           L1_SingleEG15;
   Int_t           L1_SingleEG15_Prescl;
   Int_t           L1_SingleEG15_5bx;
   Int_t           L1_SingleEG2;
   Int_t           L1_SingleEG2_Prescl;
   Int_t           L1_SingleEG2_5bx;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG20_5bx;
   Int_t           L1_SingleEG5;
   Int_t           L1_SingleEG5_Prescl;
   Int_t           L1_SingleEG5_5bx;
   Int_t           L1_SingleEG8;
   Int_t           L1_SingleEG8_Prescl;
   Int_t           L1_SingleEG8_5bx;
   Int_t           L1_SingleForJet2;
   Int_t           L1_SingleForJet2_Prescl;
   Int_t           L1_SingleForJet2_5bx;
   Int_t           L1_SingleForJet4;
   Int_t           L1_SingleForJet4_Prescl;
   Int_t           L1_SingleForJet4_5bx;
   Int_t           L1_SingleHfBitCountsRing1_1;
   Int_t           L1_SingleHfBitCountsRing1_1_Prescl;
   Int_t           L1_SingleHfBitCountsRing1_1_5bx;
   Int_t           L1_SingleHfBitCountsRing2_1;
   Int_t           L1_SingleHfBitCountsRing2_1_Prescl;
   Int_t           L1_SingleHfBitCountsRing2_1_5bx;
   Int_t           L1_SingleHfRingEtSumsRing1_200;
   Int_t           L1_SingleHfRingEtSumsRing1_200_Prescl;
   Int_t           L1_SingleHfRingEtSumsRing1_200_5bx;
   Int_t           L1_SingleHfRingEtSumsRing1_4;
   Int_t           L1_SingleHfRingEtSumsRing1_4_Prescl;
   Int_t           L1_SingleHfRingEtSumsRing1_4_5bx;
   Int_t           L1_SingleHfRingEtSumsRing2_200;
   Int_t           L1_SingleHfRingEtSumsRing2_200_Prescl;
   Int_t           L1_SingleHfRingEtSumsRing2_200_5bx;
   Int_t           L1_SingleHfRingEtSumsRing2_4;
   Int_t           L1_SingleHfRingEtSumsRing2_4_Prescl;
   Int_t           L1_SingleHfRingEtSumsRing2_4_5bx;
   Int_t           L1_SingleIsoEG10;
   Int_t           L1_SingleIsoEG10_Prescl;
   Int_t           L1_SingleIsoEG10_5bx;
   Int_t           L1_SingleIsoEG12;
   Int_t           L1_SingleIsoEG12_Prescl;
   Int_t           L1_SingleIsoEG12_5bx;
   Int_t           L1_SingleIsoEG15;
   Int_t           L1_SingleIsoEG15_Prescl;
   Int_t           L1_SingleIsoEG15_5bx;
   Int_t           L1_SingleIsoEG5;
   Int_t           L1_SingleIsoEG5_Prescl;
   Int_t           L1_SingleIsoEG5_5bx;
   Int_t           L1_SingleIsoEG8;
   Int_t           L1_SingleIsoEG8_Prescl;
   Int_t           L1_SingleIsoEG8_5bx;
   Int_t           L1_SingleJet10;
   Int_t           L1_SingleJet10_Prescl;
   Int_t           L1_SingleJet10_5bx;
   Int_t           L1_SingleJet10_NotBptxOR_Ext;
   Int_t           L1_SingleJet10_NotBptxOR_Ext_Prescl;
   Int_t           L1_SingleJet10_NotBptxOR_Ext_5bx;
   Int_t           L1_SingleJet20;
   Int_t           L1_SingleJet20_Prescl;
   Int_t           L1_SingleJet20_5bx;
   Int_t           L1_SingleJet30;
   Int_t           L1_SingleJet30_Prescl;
   Int_t           L1_SingleJet30_5bx;
   Int_t           L1_SingleJet40;
   Int_t           L1_SingleJet40_Prescl;
   Int_t           L1_SingleJet40_5bx;
   Int_t           L1_SingleJet50;
   Int_t           L1_SingleJet50_Prescl;
   Int_t           L1_SingleJet50_5bx;
   Int_t           L1_SingleJet6;
   Int_t           L1_SingleJet6_Prescl;
   Int_t           L1_SingleJet6_5bx;
   Int_t           L1_SingleJet60;
   Int_t           L1_SingleJet60_Prescl;
   Int_t           L1_SingleJet60_5bx;
   Int_t           L1_SingleMu0;
   Int_t           L1_SingleMu0_Prescl;
   Int_t           L1_SingleMu0_5bx;
   Int_t           L1_SingleMu10;
   Int_t           L1_SingleMu10_Prescl;
   Int_t           L1_SingleMu10_5bx;
   Int_t           L1_SingleMu14;
   Int_t           L1_SingleMu14_Prescl;
   Int_t           L1_SingleMu14_5bx;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu20_Prescl;
   Int_t           L1_SingleMu20_5bx;
   Int_t           L1_SingleMu3;
   Int_t           L1_SingleMu3_Prescl;
   Int_t           L1_SingleMu3_5bx;
   Int_t           L1_SingleMu5;
   Int_t           L1_SingleMu5_Prescl;
   Int_t           L1_SingleMu5_5bx;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu7_Prescl;
   Int_t           L1_SingleMu7_5bx;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuBeamHalo_5bx;
   Int_t           L1_SingleMuOpen;
   Int_t           L1_SingleMuOpen_Prescl;
   Int_t           L1_SingleMuOpen_5bx;
   Int_t           L1_SingleTauJet10;
   Int_t           L1_SingleTauJet10_Prescl;
   Int_t           L1_SingleTauJet10_5bx;
   Int_t           L1_SingleTauJet2;
   Int_t           L1_SingleTauJet2_Prescl;
   Int_t           L1_SingleTauJet2_5bx;
   Int_t           L1_SingleTauJet20;
   Int_t           L1_SingleTauJet20_Prescl;
   Int_t           L1_SingleTauJet20_5bx;
   Int_t           L1_SingleTauJet30;
   Int_t           L1_SingleTauJet30_Prescl;
   Int_t           L1_SingleTauJet30_5bx;
   Int_t           L1_SingleTauJet4;
   Int_t           L1_SingleTauJet4_Prescl;
   Int_t           L1_SingleTauJet4_5bx;
   Int_t           L1_SingleTauJet50;
   Int_t           L1_SingleTauJet50_Prescl;
   Int_t           L1_SingleTauJet50_5bx;
   Int_t           L1_TripleJet14;
   Int_t           L1_TripleJet14_Prescl;
   Int_t           L1_TripleJet14_5bx;
   Int_t           L1_ZdcLooseVertex;
   Int_t           L1_ZdcLooseVertex_Prescl;
   Int_t           L1_ZdcLooseVertex_5bx;
   Int_t           L1_ZdcMinusOverThreshold;
   Int_t           L1_ZdcMinusOverThreshold_Prescl;
   Int_t           L1_ZdcMinusOverThreshold_5bx;
   Int_t           L1_ZdcPlusOverThreshold;
   Int_t           L1_ZdcPlusOverThreshold_Prescl;
   Int_t           L1_ZdcPlusOverThreshold_5bx;
   Int_t           L1_ZdcTightVertex;
   Int_t           L1_ZdcTightVertex_Prescl;
   Int_t           L1_ZdcTightVertex_5bx;
   Int_t           L1_ZeroBias_Ext;
   Int_t           L1_ZeroBias_Ext_Prescl;
   Int_t           L1_ZeroBias_Ext_5bx;
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
   Int_t           L1Tech_ZDC_loose_vertex_v0;
   Int_t           L1Tech_ZDC_loose_vertex_v0_Prescl;
   Int_t           L1Tech_ZDC_loose_vertex_v0_5bx;
   Int_t           L1Tech_ZDC_minus_over_threshold_v0;
   Int_t           L1Tech_ZDC_minus_over_threshold_v0_Prescl;
   Int_t           L1Tech_ZDC_minus_over_threshold_v0_5bx;
   Int_t           L1Tech_ZDC_plus_over_threshold_v0;
   Int_t           L1Tech_ZDC_plus_over_threshold_v0_Prescl;
   Int_t           L1Tech_ZDC_plus_over_threshold_v0_5bx;
   Int_t           L1Tech_ZDC_tight_vertex_v0;
   Int_t           L1Tech_ZDC_tight_vertex_v0_Prescl;
   Int_t           L1Tech_ZDC_tight_vertex_v0_5bx;

   // List of branches
   TBranch        *b_NrecoJetCal;   //!
   TBranch        *b_NrecoJetGen;   //!
   TBranch        *b_NrecoTowCal;   //!
   TBranch        *b_recoJetCalPt;   //!
   TBranch        *b_recoJetCalPhi;   //!
   TBranch        *b_recoJetCalEta;   //!
   TBranch        *b_recoJetCalE;   //!
   TBranch        *b_recoJetCalEMF;   //!
   TBranch        *b_recoJetCalN90;   //!
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
   TBranch        *b_NrecoJetCorCal;   //!
   TBranch        *b_recoJetCorCalPt;   //!
   TBranch        *b_recoJetCorCalPhi;   //!
   TBranch        *b_recoJetCorCalEta;   //!
   TBranch        *b_recoJetCorCalE;   //!
   TBranch        *b_recoJetCorCalEMF;   //!
   TBranch        *b_recoJetCorCalN90;   //!
   TBranch        *b_NohTau;   //!
   TBranch        *b_ohTauEta;   //!
   TBranch        *b_ohTauPhi;   //!
   TBranch        *b_ohTauPt;   //!
   TBranch        *b_ohTauEiso;   //!
   TBranch        *b_ohTauL25Tpt;   //!
   TBranch        *b_ohTauL3Tiso;   //!
   TBranch        *b_NohpfTau;   //!
   TBranch        *b_ohpfTauPt;   //!
   TBranch        *b_ohpfTauEta;   //!
   TBranch        *b_ohpfTauPhi;   //!
   TBranch        *b_ohpfTauLeadTrackPt;   //!
   TBranch        *b_ohpfTauLeadPionPt;   //!
   TBranch        *b_ohpfTauTrkIso;   //!
   TBranch        *b_ohpfTauGammaIso;   //!
   TBranch        *b_ohpfTauJetPt;   //!
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
   TBranch        *b_pfMHT;   //!
   TBranch        *b_NohPFJet;   //!
   TBranch        *b_pfJetPt;   //!
   TBranch        *b_pfJetEta;   //!
   TBranch        *b_pfJetPhi;   //!
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
   TBranch        *b_ohBJetIPLooseL25Tag;   //!
   TBranch        *b_ohBJetIPLooseL3Tag;   //!
   TBranch        *b_ohBJetMuL25Tag;   //!
   TBranch        *b_ohBJetMuL3Tag;   //!
   TBranch        *b_ohBJetPerfL25Tag;   //!
   TBranch        *b_ohBJetPerfL3Tag;   //!
   TBranch        *b_NrecoElec;   //!
   TBranch        *b_recoElecPt;   //!
   TBranch        *b_recoElecPhi;   //!
   TBranch        *b_recoElecEta;   //!
   TBranch        *b_recoElecEt;   //!
   TBranch        *b_recoElecE;   //!
   TBranch        *b_recoElecEleID;   //!
   TBranch        *b_NrecoPhot;   //!
   TBranch        *b_recoPhotPt;   //!
   TBranch        *b_recoPhotPhi;   //!
   TBranch        *b_recoPhotEta;   //!
   TBranch        *b_recoPhotEt;   //!
   TBranch        *b_recoPhotE;   //!
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
   TBranch        *b_NohEle;   //!
   TBranch        *b_ohEleEt;   //!
   TBranch        *b_ohEleEta;   //!
   TBranch        *b_ohElePhi;   //!
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
   TBranch        *b_ohMuL2L1idx;   //!
   TBranch        *b_NohMuL3;   //!
   TBranch        *b_ohMuL3Pt;   //!
   TBranch        *b_ohMuL3Phi;   //!
   TBranch        *b_ohMuL3Eta;   //!
   TBranch        *b_ohMuL3Chg;   //!
   TBranch        *b_ohMuL3PtErr;   //!
   TBranch        *b_ohMuL3Iso;   //!
   TBranch        *b_ohMuL3Dr;   //!
   TBranch        *b_ohMuL3Dz;   //!
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
   TBranch        *b_ohMuL2NoVtxL1idx;   //!
   TBranch        *b_ohMuL2NoVtxNhits;   //!
   TBranch        *b_ohMuL2NoVtxNchambers;   //!
   TBranch        *b_ohHighestEnergyEERecHit;   //!
   TBranch        *b_ohHighestEnergyEBRecHit;   //!
   TBranch        *b_ohHighestEnergyHBHERecHit;   //!
   TBranch        *b_ohHighestEnergyHORecHit;   //!
   TBranch        *b_ohHighestEnergyHFRecHit;   //!
   TBranch        *b_Nalcapi0clusters;   //!
   TBranch        *b_ohAlcapi0ptClusAll;   //!
   TBranch        *b_ohAlcapi0etaClusAll;   //!
   TBranch        *b_ohAlcapi0phiClusAll;   //!
   TBranch        *b_ohAlcapi0s4s9ClusAll;   //!
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
   TBranch        *b_recoVrtX;   //!
   TBranch        *b_recoVrtY;   //!
   TBranch        *b_recoVrtZ;   //!
   TBranch        *b_recoVrtNtrk;   //!
   TBranch        *b_recoVrtChi2;   //!
   TBranch        *b_recoVrtNdof;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_HLT_Activity_CSC;   //!
   TBranch        *b_HLT_Activity_CSC_Prescl;   //!
   TBranch        *b_HLT_Activity_DT;   //!
   TBranch        *b_HLT_Activity_DT_Prescl;   //!
   TBranch        *b_HLT_Activity_DT_Tuned;   //!
   TBranch        *b_HLT_Activity_DT_Tuned_Prescl;   //!
   TBranch        *b_HLT_Activity_Ecal_SC7;   //!
   TBranch        *b_HLT_Activity_Ecal_SC7_Prescl;   //!
   TBranch        *b_HLT_Activity_Ecal_SC17;   //!
   TBranch        *b_HLT_Activity_Ecal_SC17_Prescl;   //!
   TBranch        *b_HLT_L1Jet6U;   //!
   TBranch        *b_HLT_L1Jet6U_Prescl;   //!
   TBranch        *b_HLT_L1Jet10U;   //!
   TBranch        *b_HLT_L1Jet10U_Prescl;   //!
   TBranch        *b_HLT_Jet15U;   //!
   TBranch        *b_HLT_Jet15U_Prescl;   //!
   TBranch        *b_HLT_Jet15U_HcalNoiseFiltered;   //!
   TBranch        *b_HLT_Jet15U_HcalNoiseFiltered_Prescl;   //!
   TBranch        *b_HLT_Jet30U;   //!
   TBranch        *b_HLT_Jet30U_Prescl;   //!
   TBranch        *b_HLT_Jet50U;   //!
   TBranch        *b_HLT_Jet50U_Prescl;   //!
   TBranch        *b_HLT_Jet70U_v2;   //!
   TBranch        *b_HLT_Jet70U_v2_Prescl;   //!
   TBranch        *b_HLT_Jet100U_v2;   //!
   TBranch        *b_HLT_Jet100U_v2_Prescl;   //!
   TBranch        *b_HLT_Jet140U_v1;   //!
   TBranch        *b_HLT_Jet140U_v1_Prescl;   //!
   TBranch        *b_HLT_DiJetAve15U;   //!
   TBranch        *b_HLT_DiJetAve15U_Prescl;   //!
   TBranch        *b_HLT_DiJetAve30U;   //!
   TBranch        *b_HLT_DiJetAve30U_Prescl;   //!
   TBranch        *b_HLT_DiJetAve50U;   //!
   TBranch        *b_HLT_DiJetAve50U_Prescl;   //!
   TBranch        *b_HLT_DiJetAve70U_v2;   //!
   TBranch        *b_HLT_DiJetAve70U_v2_Prescl;   //!
   TBranch        *b_HLT_DiJetAve100U_v1;   //!
   TBranch        *b_HLT_DiJetAve100U_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleJet15U_ForwardBackward;   //!
   TBranch        *b_HLT_DoubleJet15U_ForwardBackward_Prescl;   //!
   TBranch        *b_HLT_DoubleJet25U_ForwardBackward;   //!
   TBranch        *b_HLT_DoubleJet25U_ForwardBackward_Prescl;   //!
   TBranch        *b_HLT_ExclDiJet30U_HFAND_v1;   //!
   TBranch        *b_HLT_ExclDiJet30U_HFAND_v1_Prescl;   //!
   TBranch        *b_HLT_ExclDiJet30U_HFOR_v1;   //!
   TBranch        *b_HLT_ExclDiJet30U_HFOR_v1_Prescl;   //!
   TBranch        *b_HLT_QuadJet15U_v2;   //!
   TBranch        *b_HLT_QuadJet15U_v2_Prescl;   //!
   TBranch        *b_HLT_QuadJet20U_v2;   //!
   TBranch        *b_HLT_QuadJet20U_v2_Prescl;   //!
   TBranch        *b_HLT_QuadJet25U_v2;   //!
   TBranch        *b_HLT_QuadJet25U_v2_Prescl;   //!
   TBranch        *b_HLT_L1ETT100;   //!
   TBranch        *b_HLT_L1ETT100_Prescl;   //!
   TBranch        *b_HLT_L1ETT140_v1;   //!
   TBranch        *b_HLT_L1ETT140_v1_Prescl;   //!
   TBranch        *b_HLT_EcalOnly_SumEt160_v2;   //!
   TBranch        *b_HLT_EcalOnly_SumEt160_v2_Prescl;   //!
   TBranch        *b_HLT_L1MET20;   //!
   TBranch        *b_HLT_L1MET20_Prescl;   //!
   TBranch        *b_HLT_MET45;   //!
   TBranch        *b_HLT_MET45_Prescl;   //!
   TBranch        *b_HLT_MET45_HT100U_v1;   //!
   TBranch        *b_HLT_MET45_HT100U_v1_Prescl;   //!
   TBranch        *b_HLT_MET45_HT120U_v1;   //!
   TBranch        *b_HLT_MET45_HT120U_v1_Prescl;   //!
   TBranch        *b_HLT_MET65;   //!
   TBranch        *b_HLT_MET65_Prescl;   //!
   TBranch        *b_HLT_MET80_v1;   //!
   TBranch        *b_HLT_MET80_v1_Prescl;   //!
   TBranch        *b_HLT_MET100_v2;   //!
   TBranch        *b_HLT_MET100_v2_Prescl;   //!
   TBranch        *b_HLT_HT50U_v1;   //!
   TBranch        *b_HLT_HT50U_v1_Prescl;   //!
   TBranch        *b_HLT_HT100U;   //!
   TBranch        *b_HLT_HT100U_Prescl;   //!
   TBranch        *b_HLT_HT120U;   //!
   TBranch        *b_HLT_HT120U_Prescl;   //!
   TBranch        *b_HLT_HT140U;   //!
   TBranch        *b_HLT_HT140U_Prescl;   //!
   TBranch        *b_HLT_HT140U_Eta3_v1;   //!
   TBranch        *b_HLT_HT140U_Eta3_v1_Prescl;   //!
   TBranch        *b_HLT_HT160U_v1;   //!
   TBranch        *b_HLT_HT160U_v1_Prescl;   //!
   TBranch        *b_HLT_HT200U_v1;   //!
   TBranch        *b_HLT_HT200U_v1_Prescl;   //!
   TBranch        *b_HLT_L1MuOpen;   //!
   TBranch        *b_HLT_L1MuOpen_Prescl;   //!
   TBranch        *b_HLT_L1MuOpen_DT;   //!
   TBranch        *b_HLT_L1MuOpen_DT_Prescl;   //!
   TBranch        *b_HLT_L1MuOpen_AntiBPTX;   //!
   TBranch        *b_HLT_L1MuOpen_AntiBPTX_Prescl;   //!
   TBranch        *b_HLT_L1Mu7_v1;   //!
   TBranch        *b_HLT_L1Mu7_v1_Prescl;   //!
   TBranch        *b_HLT_L1Mu20;   //!
   TBranch        *b_HLT_L1Mu20_Prescl;   //!
   TBranch        *b_HLT_L2Mu0_NoVertex;   //!
   TBranch        *b_HLT_L2Mu0_NoVertex_Prescl;   //!
   TBranch        *b_HLT_L2Mu7_v1;   //!
   TBranch        *b_HLT_L2Mu7_v1_Prescl;   //!
   TBranch        *b_HLT_L2Mu30_v1;   //!
   TBranch        *b_HLT_L2Mu30_v1_Prescl;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu3_Prescl;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu5_Prescl;   //!
   TBranch        *b_HLT_Mu7;   //!
   TBranch        *b_HLT_Mu7_Prescl;   //!
   TBranch        *b_HLT_Mu9;   //!
   TBranch        *b_HLT_Mu9_Prescl;   //!
   TBranch        *b_HLT_Mu11;   //!
   TBranch        *b_HLT_Mu11_Prescl;   //!
   TBranch        *b_HLT_Mu13_v1;   //!
   TBranch        *b_HLT_Mu13_v1_Prescl;   //!
   TBranch        *b_HLT_Mu15_v1;   //!
   TBranch        *b_HLT_Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_IsoMu9;   //!
   TBranch        *b_HLT_IsoMu9_Prescl;   //!
   TBranch        *b_HLT_IsoMu11_v1;   //!
   TBranch        *b_HLT_IsoMu11_v1_Prescl;   //!
   TBranch        *b_HLT_Mu20_NoVertex;   //!
   TBranch        *b_HLT_Mu20_NoVertex_Prescl;   //!
   TBranch        *b_HLT_L1DoubleMuOpen;   //!
   TBranch        *b_HLT_L1DoubleMuOpen_Prescl;   //!
   TBranch        *b_HLT_L2DoubleMu0;   //!
   TBranch        *b_HLT_L2DoubleMu0_Prescl;   //!
   TBranch        *b_HLT_L2DoubleMu20_NoVertex_v1;   //!
   TBranch        *b_HLT_L2DoubleMu20_NoVertex_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_v1;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_LS_v1;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_LS_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleMu0;   //!
   TBranch        *b_HLT_DoubleMu0_Prescl;   //!
   TBranch        *b_HLT_DoubleMu3_v2;   //!
   TBranch        *b_HLT_DoubleMu3_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleMu5_v1;   //!
   TBranch        *b_HLT_DoubleMu5_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_L2Mu0;   //!
   TBranch        *b_HLT_Mu5_L2Mu0_Prescl;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_Prescl;   //!
   TBranch        *b_HLT_Mu3_Track5_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu3_Track5_Jpsi_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu5_Track0_Jpsi_Prescl;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Prescl;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1_Prescl;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Prescl;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Prescl;   //!
   TBranch        *b_HLT_L1SingleEG2;   //!
   TBranch        *b_HLT_L1SingleEG2_Prescl;   //!
   TBranch        *b_HLT_L1SingleEG8;   //!
   TBranch        *b_HLT_L1SingleEG8_Prescl;   //!
   TBranch        *b_HLT_Ele10_SW_L1R;   //!
   TBranch        *b_HLT_Ele10_SW_L1R_Prescl;   //!
   TBranch        *b_HLT_Ele12_SW_TightEleId_L1R;   //!
   TBranch        *b_HLT_Ele12_SW_TightEleId_L1R_Prescl;   //!
   TBranch        *b_HLT_Ele12_SW_TighterEleId_L1R_v1;   //!
   TBranch        *b_HLT_Ele12_SW_TighterEleId_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele12_SW_TighterEleIdIsol_L1R_v1;   //!
   TBranch        *b_HLT_Ele12_SW_TighterEleIdIsol_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_L1R;   //!
   TBranch        *b_HLT_Ele17_SW_L1R_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_TightEleId_L1R;   //!
   TBranch        *b_HLT_Ele17_SW_TightEleId_L1R_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_TighterEleId_L1R_v1;   //!
   TBranch        *b_HLT_Ele17_SW_TighterEleId_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_TightEleIdIsol_L1R_v1;   //!
   TBranch        *b_HLT_Ele17_SW_TightEleIdIsol_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_TighterEleIdIsol_L1R_v1;   //!
   TBranch        *b_HLT_Ele17_SW_TighterEleIdIsol_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1;   //!
   TBranch        *b_HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1;   //!
   TBranch        *b_HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1;   //!
   TBranch        *b_HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1;   //!
   TBranch        *b_HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_DoubleEle4_SW_eeRes_L1R;   //!
   TBranch        *b_HLT_DoubleEle4_SW_eeRes_L1R_Prescl;   //!
   TBranch        *b_HLT_DoubleEle15_SW_L1R_v1;   //!
   TBranch        *b_HLT_DoubleEle15_SW_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Photon10_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon10_Cleaned_L1R_Prescl;   //!
   TBranch        *b_HLT_Photon15_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon15_Cleaned_L1R_Prescl;   //!
   TBranch        *b_HLT_Photon17_SC17HE_L1R_v1;   //!
   TBranch        *b_HLT_Photon17_SC17HE_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Photon20_NoHE_L1R;   //!
   TBranch        *b_HLT_Photon20_NoHE_L1R_Prescl;   //!
   TBranch        *b_HLT_Photon20_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon20_Cleaned_L1R_Prescl;   //!
   TBranch        *b_HLT_Photon30_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon30_Cleaned_L1R_Prescl;   //!
   TBranch        *b_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Photon35_Isol_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon35_Isol_Cleaned_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Photon50_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon50_Cleaned_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Photon50_NoHE_L1R;   //!
   TBranch        *b_HLT_Photon50_NoHE_L1R_Prescl;   //!
   TBranch        *b_HLT_Photon70_NoHE_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon70_NoHE_Cleaned_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Photon100_NoHE_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon100_NoHE_Cleaned_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_DoublePhoton5_CEP_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_CEP_L1R_Prescl;   //!
   TBranch        *b_HLT_DoublePhoton17_L1R;   //!
   TBranch        *b_HLT_DoublePhoton17_L1R_Prescl;   //!
   TBranch        *b_HLT_SingleIsoTau20_Trk5_MET20;   //!
   TBranch        *b_HLT_SingleIsoTau20_Trk5_MET20_Prescl;   //!
   TBranch        *b_HLT_SingleIsoTau20_Trk15_MET20;   //!
   TBranch        *b_HLT_SingleIsoTau20_Trk15_MET20_Prescl;   //!
   TBranch        *b_HLT_SingleIsoTau30_Trk5_MET20;   //!
   TBranch        *b_HLT_SingleIsoTau30_Trk5_MET20_Prescl;   //!
   TBranch        *b_HLT_SingleIsoTau30_Trk5_v2;   //!
   TBranch        *b_HLT_SingleIsoTau30_Trk5_v2_Prescl;   //!
   TBranch        *b_HLT_DoubleIsoTau15_OneLeg_Trk5;   //!
   TBranch        *b_HLT_DoubleIsoTau15_OneLeg_Trk5_Prescl;   //!
   TBranch        *b_HLT_DoubleIsoTau15_Trk5;   //!
   TBranch        *b_HLT_DoubleIsoTau15_Trk5_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet10U_v1;   //!
   TBranch        *b_HLT_BTagMu_DiJet10U_v1_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet20U_v1;   //!
   TBranch        *b_HLT_BTagMu_DiJet20U_v1_Prescl;   //!
   TBranch        *b_HLT_BTagMu_DiJet20U_Mu5_v1;   //!
   TBranch        *b_HLT_BTagMu_DiJet20U_Mu5_v1_Prescl;   //!
   TBranch        *b_HLT_StoppedHSCP20_v3;   //!
   TBranch        *b_HLT_StoppedHSCP20_v3_Prescl;   //!
   TBranch        *b_HLT_StoppedHSCP35_v3;   //!
   TBranch        *b_HLT_StoppedHSCP35_v3_Prescl;   //!
   TBranch        *b_HLT_Mu5_Photon11_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Mu5_Photon11_Cleaned_L1R_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_Ele5_v1;   //!
   TBranch        *b_HLT_Mu5_Ele5_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_Ele9_v1;   //!
   TBranch        *b_HLT_Mu5_Ele9_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_Jet35U_v1;   //!
   TBranch        *b_HLT_Mu5_Jet35U_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_Jet50U_v2;   //!
   TBranch        *b_HLT_Mu5_Jet50U_v2_Prescl;   //!
   TBranch        *b_HLT_Mu5_MET45_v1;   //!
   TBranch        *b_HLT_Mu5_MET45_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_HT50U_v1;   //!
   TBranch        *b_HLT_Mu5_HT50U_v1_Prescl;   //!
   TBranch        *b_HLT_Mu5_HT70U_v1;   //!
   TBranch        *b_HLT_Mu5_HT70U_v1_Prescl;   //!
   TBranch        *b_HLT_Ele10_MET45_v1;   //!
   TBranch        *b_HLT_Ele10_MET45_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLT_ZeroBias_Prescl;   //!
   TBranch        *b_HLT_ZeroBiasPixel_SingleTrack;   //!
   TBranch        *b_HLT_ZeroBiasPixel_SingleTrack_Prescl;   //!
   TBranch        *b_HLT_MinBiasPixel_SingleTrack;   //!
   TBranch        *b_HLT_MinBiasPixel_SingleTrack_Prescl;   //!
   TBranch        *b_HLT_MultiVertex6;   //!
   TBranch        *b_HLT_MultiVertex6_Prescl;   //!
   TBranch        *b_HLT_MultiVertex8_L1ETT60;   //!
   TBranch        *b_HLT_MultiVertex8_L1ETT60_Prescl;   //!
   TBranch        *b_HLT_L1_BptxXOR_BscMinBiasOR;   //!
   TBranch        *b_HLT_L1_BptxXOR_BscMinBiasOR_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_minBias_OR;   //!
   TBranch        *b_HLT_L1Tech_BSC_minBias_OR_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_minBias;   //!
   TBranch        *b_HLT_L1Tech_BSC_minBias_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_halo;   //!
   TBranch        *b_HLT_L1Tech_BSC_halo_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_halo_forPhysicsBackground;   //!
   TBranch        *b_HLT_L1Tech_BSC_halo_forPhysicsBackground_Prescl;   //!
   TBranch        *b_HLT_L1Tech_BSC_HighMultiplicity;   //!
   TBranch        *b_HLT_L1Tech_BSC_HighMultiplicity_Prescl;   //!
   TBranch        *b_HLT_L1Tech_RPC_TTU_RBst1_collisions;   //!
   TBranch        *b_HLT_L1Tech_RPC_TTU_RBst1_collisions_Prescl;   //!
   TBranch        *b_HLT_L1Tech_HCAL_HF;   //!
   TBranch        *b_HLT_L1Tech_HCAL_HF_Prescl;   //!
   TBranch        *b_HLT_TrackerCosmics;   //!
   TBranch        *b_HLT_TrackerCosmics_Prescl;   //!
   TBranch        *b_HLT_IsoTrackHB_v2;   //!
   TBranch        *b_HLT_IsoTrackHB_v2_Prescl;   //!
   TBranch        *b_HLT_IsoTrackHE_v2;   //!
   TBranch        *b_HLT_IsoTrackHE_v2_Prescl;   //!
   TBranch        *b_HLT_RPCBarrelCosmics;   //!
   TBranch        *b_HLT_RPCBarrelCosmics_Prescl;   //!
   TBranch        *b_HLT_HcalPhiSym;   //!
   TBranch        *b_HLT_HcalPhiSym_Prescl;   //!
   TBranch        *b_HLT_HcalNZS;   //!
   TBranch        *b_HLT_HcalNZS_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity70;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity70_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity100;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity100_Prescl;   //!
   TBranch        *b_HLT_GlobalRunHPDNoise;   //!
   TBranch        *b_HLT_GlobalRunHPDNoise_Prescl;   //!
   TBranch        *b_HLT_TechTrigHCALNoise;   //!
   TBranch        *b_HLT_TechTrigHCALNoise_Prescl;   //!
   TBranch        *b_HLT_L1_BPTX;   //!
   TBranch        *b_HLT_L1_BPTX_Prescl;   //!
   TBranch        *b_HLT_L1_BPTX_MinusOnly;   //!
   TBranch        *b_HLT_L1_BPTX_MinusOnly_Prescl;   //!
   TBranch        *b_HLT_L1_BPTX_PlusOnly;   //!
   TBranch        *b_HLT_L1_BPTX_PlusOnly_Prescl;   //!
   TBranch        *b_HLT_DTErrors;   //!
   TBranch        *b_HLT_DTErrors_Prescl;   //!
   TBranch        *b_HLT_LogMonitor;   //!
   TBranch        *b_HLT_LogMonitor_Prescl;   //!
   TBranch        *b_HLT_Calibration;   //!
   TBranch        *b_HLT_Calibration_Prescl;   //!
   TBranch        *b_HLT_EcalCalibration;   //!
   TBranch        *b_HLT_EcalCalibration_Prescl;   //!
   TBranch        *b_HLT_HcalCalibration;   //!
   TBranch        *b_HLT_HcalCalibration_Prescl;   //!
   TBranch        *b_HLT_Random;   //!
   TBranch        *b_HLT_Random_Prescl;   //!
   TBranch        *b_AlCa_EcalPhiSym;   //!
   TBranch        *b_AlCa_EcalPhiSym_Prescl;   //!
   TBranch        *b_AlCa_EcalPi0;   //!
   TBranch        *b_AlCa_EcalPi0_Prescl;   //!
   TBranch        *b_AlCa_EcalEta;   //!
   TBranch        *b_AlCa_EcalEta_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoHits;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_Prescl;   //!
   TBranch        *b_DQM_FEDIntegrity;   //!
   TBranch        *b_DQM_FEDIntegrity_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxMinus_5bx;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_5bx;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxPlus_5bx;   //!
   TBranch        *b_L1_BptxPlusORMinus;   //!
   TBranch        *b_L1_BptxPlusORMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlusORMinus_5bx;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_5bx;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR_5bx;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus_Prescl;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus_5bx;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus_Prescl;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus_5bx;   //!
   TBranch        *b_L1_BscHaloBeam1Inner;   //!
   TBranch        *b_L1_BscHaloBeam1Inner_Prescl;   //!
   TBranch        *b_L1_BscHaloBeam1Inner_5bx;   //!
   TBranch        *b_L1_BscHaloBeam1Outer;   //!
   TBranch        *b_L1_BscHaloBeam1Outer_Prescl;   //!
   TBranch        *b_L1_BscHaloBeam1Outer_5bx;   //!
   TBranch        *b_L1_BscHaloBeam2Inner;   //!
   TBranch        *b_L1_BscHaloBeam2Inner_Prescl;   //!
   TBranch        *b_L1_BscHaloBeam2Inner_5bx;   //!
   TBranch        *b_L1_BscHaloBeam2Outer;   //!
   TBranch        *b_L1_BscHaloBeam2Outer_Prescl;   //!
   TBranch        *b_L1_BscHaloBeam2Outer_5bx;   //!
   TBranch        *b_L1_BscHighMultiplicity;   //!
   TBranch        *b_L1_BscHighMultiplicity_Prescl;   //!
   TBranch        *b_L1_BscHighMultiplicity_5bx;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_5bx;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2_5bx;   //!
   TBranch        *b_L1_BscMinBiasOR;   //!
   TBranch        *b_L1_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_5bx;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusANDMinus;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusANDMinus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusANDMinus_5bx;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus_5bx;   //!
   TBranch        *b_L1_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_5bx;   //!
   TBranch        *b_L1_BscMinBiasThreshold2;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_5bx;   //!
   TBranch        *b_L1_BscSplashBeam1;   //!
   TBranch        *b_L1_BscSplashBeam1_Prescl;   //!
   TBranch        *b_L1_BscSplashBeam1_5bx;   //!
   TBranch        *b_L1_BscSplashBeam2;   //!
   TBranch        *b_L1_BscSplashBeam2_Prescl;   //!
   TBranch        *b_L1_BscSplashBeam2_5bx;   //!
   TBranch        *b_L1_DoubleEG05_TopBottom;   //!
   TBranch        *b_L1_DoubleEG05_TopBottom_Prescl;   //!
   TBranch        *b_L1_DoubleEG05_TopBottom_5bx;   //!
   TBranch        *b_L1_DoubleEG2;   //!
   TBranch        *b_L1_DoubleEG2_Prescl;   //!
   TBranch        *b_L1_DoubleEG2_5bx;   //!
   TBranch        *b_L1_DoubleEG5;   //!
   TBranch        *b_L1_DoubleEG5_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_5bx;   //!
   TBranch        *b_L1_DoubleForJet10_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet10_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleForJet10_EtaOpp_5bx;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1_Prescl;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1_5bx;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1_Prescl;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1_5bx;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200_5bx;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4_5bx;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200_5bx;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4_5bx;   //!
   TBranch        *b_L1_DoubleJet30;   //!
   TBranch        *b_L1_DoubleJet30_Prescl;   //!
   TBranch        *b_L1_DoubleJet30_5bx;   //!
   TBranch        *b_L1_DoubleMu3;   //!
   TBranch        *b_L1_DoubleMu3_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_5bx;   //!
   TBranch        *b_L1_DoubleMuOpen;   //!
   TBranch        *b_L1_DoubleMuOpen_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_5bx;   //!
   TBranch        *b_L1_DoubleMuTopBottom;   //!
   TBranch        *b_L1_DoubleMuTopBottom_Prescl;   //!
   TBranch        *b_L1_DoubleMuTopBottom_5bx;   //!
   TBranch        *b_L1_DoubleTauJet14;   //!
   TBranch        *b_L1_DoubleTauJet14_Prescl;   //!
   TBranch        *b_L1_DoubleTauJet14_5bx;   //!
   TBranch        *b_L1_ETM12;   //!
   TBranch        *b_L1_ETM12_Prescl;   //!
   TBranch        *b_L1_ETM12_5bx;   //!
   TBranch        *b_L1_ETM20;   //!
   TBranch        *b_L1_ETM20_Prescl;   //!
   TBranch        *b_L1_ETM20_5bx;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM30_Prescl;   //!
   TBranch        *b_L1_ETM30_5bx;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM70_Prescl;   //!
   TBranch        *b_L1_ETM70_5bx;   //!
   TBranch        *b_L1_ETT100;   //!
   TBranch        *b_L1_ETT100_Prescl;   //!
   TBranch        *b_L1_ETT100_5bx;   //!
   TBranch        *b_L1_ETT140;   //!
   TBranch        *b_L1_ETT140_Prescl;   //!
   TBranch        *b_L1_ETT140_5bx;   //!
   TBranch        *b_L1_ETT30;   //!
   TBranch        *b_L1_ETT30_Prescl;   //!
   TBranch        *b_L1_ETT30_5bx;   //!
   TBranch        *b_L1_ETT60;   //!
   TBranch        *b_L1_ETT60_Prescl;   //!
   TBranch        *b_L1_ETT60_5bx;   //!
   TBranch        *b_L1_HTM20;   //!
   TBranch        *b_L1_HTM20_Prescl;   //!
   TBranch        *b_L1_HTM20_5bx;   //!
   TBranch        *b_L1_HTM30;   //!
   TBranch        *b_L1_HTM30_Prescl;   //!
   TBranch        *b_L1_HTM30_5bx;   //!
   TBranch        *b_L1_HTT100;   //!
   TBranch        *b_L1_HTT100_Prescl;   //!
   TBranch        *b_L1_HTT100_5bx;   //!
   TBranch        *b_L1_HTT200;   //!
   TBranch        *b_L1_HTT200_Prescl;   //!
   TBranch        *b_L1_HTT200_5bx;   //!
   TBranch        *b_L1_HTT50;   //!
   TBranch        *b_L1_HTT50_Prescl;   //!
   TBranch        *b_L1_HTT50_5bx;   //!
   TBranch        *b_L1_IsoEG10_Jet6_ForJet6;   //!
   TBranch        *b_L1_IsoEG10_Jet6_ForJet6_Prescl;   //!
   TBranch        *b_L1_IsoEG10_Jet6_ForJet6_5bx;   //!
   TBranch        *b_L1_Mu3_EG5;   //!
   TBranch        *b_L1_Mu3_EG5_Prescl;   //!
   TBranch        *b_L1_Mu3_EG5_5bx;   //!
   TBranch        *b_L1_Mu3_Jet10;   //!
   TBranch        *b_L1_Mu3_Jet10_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet10_5bx;   //!
   TBranch        *b_L1_Mu3_Jet6;   //!
   TBranch        *b_L1_Mu3_Jet6_Prescl;   //!
   TBranch        *b_L1_Mu3_Jet6_5bx;   //!
   TBranch        *b_L1_Mu5_Jet6;   //!
   TBranch        *b_L1_Mu5_Jet6_Prescl;   //!
   TBranch        *b_L1_Mu5_Jet6_5bx;   //!
   TBranch        *b_L1_QuadJet6;   //!
   TBranch        *b_L1_QuadJet6_Prescl;   //!
   TBranch        *b_L1_QuadJet6_5bx;   //!
   TBranch        *b_L1_QuadJet8;   //!
   TBranch        *b_L1_QuadJet8_Prescl;   //!
   TBranch        *b_L1_QuadJet8_5bx;   //!
   TBranch        *b_L1_SingleCenJet2;   //!
   TBranch        *b_L1_SingleCenJet2_Prescl;   //!
   TBranch        *b_L1_SingleCenJet2_5bx;   //!
   TBranch        *b_L1_SingleCenJet4;   //!
   TBranch        *b_L1_SingleCenJet4_Prescl;   //!
   TBranch        *b_L1_SingleCenJet4_5bx;   //!
   TBranch        *b_L1_SingleEG1;   //!
   TBranch        *b_L1_SingleEG1_Prescl;   //!
   TBranch        *b_L1_SingleEG1_5bx;   //!
   TBranch        *b_L1_SingleEG10;   //!
   TBranch        *b_L1_SingleEG10_Prescl;   //!
   TBranch        *b_L1_SingleEG10_5bx;   //!
   TBranch        *b_L1_SingleEG12;   //!
   TBranch        *b_L1_SingleEG12_Prescl;   //!
   TBranch        *b_L1_SingleEG12_5bx;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG15_Prescl;   //!
   TBranch        *b_L1_SingleEG15_5bx;   //!
   TBranch        *b_L1_SingleEG2;   //!
   TBranch        *b_L1_SingleEG2_Prescl;   //!
   TBranch        *b_L1_SingleEG2_5bx;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG20_5bx;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG5_Prescl;   //!
   TBranch        *b_L1_SingleEG5_5bx;   //!
   TBranch        *b_L1_SingleEG8;   //!
   TBranch        *b_L1_SingleEG8_Prescl;   //!
   TBranch        *b_L1_SingleEG8_5bx;   //!
   TBranch        *b_L1_SingleForJet2;   //!
   TBranch        *b_L1_SingleForJet2_Prescl;   //!
   TBranch        *b_L1_SingleForJet2_5bx;   //!
   TBranch        *b_L1_SingleForJet4;   //!
   TBranch        *b_L1_SingleForJet4_Prescl;   //!
   TBranch        *b_L1_SingleForJet4_5bx;   //!
   TBranch        *b_L1_SingleHfBitCountsRing1_1;   //!
   TBranch        *b_L1_SingleHfBitCountsRing1_1_Prescl;   //!
   TBranch        *b_L1_SingleHfBitCountsRing1_1_5bx;   //!
   TBranch        *b_L1_SingleHfBitCountsRing2_1;   //!
   TBranch        *b_L1_SingleHfBitCountsRing2_1_Prescl;   //!
   TBranch        *b_L1_SingleHfBitCountsRing2_1_5bx;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing1_200;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing1_200_Prescl;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing1_200_5bx;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing1_4;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing1_4_Prescl;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing1_4_5bx;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing2_200;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing2_200_Prescl;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing2_200_5bx;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing2_4;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing2_4_Prescl;   //!
   TBranch        *b_L1_SingleHfRingEtSumsRing2_4_5bx;   //!
   TBranch        *b_L1_SingleIsoEG10;   //!
   TBranch        *b_L1_SingleIsoEG10_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG10_5bx;   //!
   TBranch        *b_L1_SingleIsoEG12;   //!
   TBranch        *b_L1_SingleIsoEG12_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG12_5bx;   //!
   TBranch        *b_L1_SingleIsoEG15;   //!
   TBranch        *b_L1_SingleIsoEG15_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG15_5bx;   //!
   TBranch        *b_L1_SingleIsoEG5;   //!
   TBranch        *b_L1_SingleIsoEG5_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG5_5bx;   //!
   TBranch        *b_L1_SingleIsoEG8;   //!
   TBranch        *b_L1_SingleIsoEG8_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG8_5bx;   //!
   TBranch        *b_L1_SingleJet10;   //!
   TBranch        *b_L1_SingleJet10_Prescl;   //!
   TBranch        *b_L1_SingleJet10_5bx;   //!
   TBranch        *b_L1_SingleJet10_NotBptxOR_Ext;   //!
   TBranch        *b_L1_SingleJet10_NotBptxOR_Ext_Prescl;   //!
   TBranch        *b_L1_SingleJet10_NotBptxOR_Ext_5bx;   //!
   TBranch        *b_L1_SingleJet20;   //!
   TBranch        *b_L1_SingleJet20_Prescl;   //!
   TBranch        *b_L1_SingleJet20_5bx;   //!
   TBranch        *b_L1_SingleJet30;   //!
   TBranch        *b_L1_SingleJet30_Prescl;   //!
   TBranch        *b_L1_SingleJet30_5bx;   //!
   TBranch        *b_L1_SingleJet40;   //!
   TBranch        *b_L1_SingleJet40_Prescl;   //!
   TBranch        *b_L1_SingleJet40_5bx;   //!
   TBranch        *b_L1_SingleJet50;   //!
   TBranch        *b_L1_SingleJet50_Prescl;   //!
   TBranch        *b_L1_SingleJet50_5bx;   //!
   TBranch        *b_L1_SingleJet6;   //!
   TBranch        *b_L1_SingleJet6_Prescl;   //!
   TBranch        *b_L1_SingleJet6_5bx;   //!
   TBranch        *b_L1_SingleJet60;   //!
   TBranch        *b_L1_SingleJet60_Prescl;   //!
   TBranch        *b_L1_SingleJet60_5bx;   //!
   TBranch        *b_L1_SingleMu0;   //!
   TBranch        *b_L1_SingleMu0_Prescl;   //!
   TBranch        *b_L1_SingleMu0_5bx;   //!
   TBranch        *b_L1_SingleMu10;   //!
   TBranch        *b_L1_SingleMu10_Prescl;   //!
   TBranch        *b_L1_SingleMu10_5bx;   //!
   TBranch        *b_L1_SingleMu14;   //!
   TBranch        *b_L1_SingleMu14_Prescl;   //!
   TBranch        *b_L1_SingleMu14_5bx;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20_Prescl;   //!
   TBranch        *b_L1_SingleMu20_5bx;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu3_Prescl;   //!
   TBranch        *b_L1_SingleMu3_5bx;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu5_Prescl;   //!
   TBranch        *b_L1_SingleMu5_5bx;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_Prescl;   //!
   TBranch        *b_L1_SingleMu7_5bx;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo_5bx;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_5bx;   //!
   TBranch        *b_L1_SingleTauJet10;   //!
   TBranch        *b_L1_SingleTauJet10_Prescl;   //!
   TBranch        *b_L1_SingleTauJet10_5bx;   //!
   TBranch        *b_L1_SingleTauJet2;   //!
   TBranch        *b_L1_SingleTauJet2_Prescl;   //!
   TBranch        *b_L1_SingleTauJet2_5bx;   //!
   TBranch        *b_L1_SingleTauJet20;   //!
   TBranch        *b_L1_SingleTauJet20_Prescl;   //!
   TBranch        *b_L1_SingleTauJet20_5bx;   //!
   TBranch        *b_L1_SingleTauJet30;   //!
   TBranch        *b_L1_SingleTauJet30_Prescl;   //!
   TBranch        *b_L1_SingleTauJet30_5bx;   //!
   TBranch        *b_L1_SingleTauJet4;   //!
   TBranch        *b_L1_SingleTauJet4_Prescl;   //!
   TBranch        *b_L1_SingleTauJet4_5bx;   //!
   TBranch        *b_L1_SingleTauJet50;   //!
   TBranch        *b_L1_SingleTauJet50_Prescl;   //!
   TBranch        *b_L1_SingleTauJet50_5bx;   //!
   TBranch        *b_L1_TripleJet14;   //!
   TBranch        *b_L1_TripleJet14_Prescl;   //!
   TBranch        *b_L1_TripleJet14_5bx;   //!
   TBranch        *b_L1_ZdcLooseVertex;   //!
   TBranch        *b_L1_ZdcLooseVertex_Prescl;   //!
   TBranch        *b_L1_ZdcLooseVertex_5bx;   //!
   TBranch        *b_L1_ZdcMinusOverThreshold;   //!
   TBranch        *b_L1_ZdcMinusOverThreshold_Prescl;   //!
   TBranch        *b_L1_ZdcMinusOverThreshold_5bx;   //!
   TBranch        *b_L1_ZdcPlusOverThreshold;   //!
   TBranch        *b_L1_ZdcPlusOverThreshold_Prescl;   //!
   TBranch        *b_L1_ZdcPlusOverThreshold_5bx;   //!
   TBranch        *b_L1_ZdcTightVertex;   //!
   TBranch        *b_L1_ZdcTightVertex_Prescl;   //!
   TBranch        *b_L1_ZdcTightVertex_5bx;   //!
   TBranch        *b_L1_ZeroBias_Ext;   //!
   TBranch        *b_L1_ZeroBias_Ext_Prescl;   //!
   TBranch        *b_L1_ZeroBias_Ext_5bx;   //!
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
   TBranch        *b_L1Tech_ZDC_loose_vertex_v0;   //!
   TBranch        *b_L1Tech_ZDC_loose_vertex_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_loose_vertex_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_minus_over_threshold_v0;   //!
   TBranch        *b_L1Tech_ZDC_minus_over_threshold_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_minus_over_threshold_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_plus_over_threshold_v0;   //!
   TBranch        *b_L1Tech_ZDC_plus_over_threshold_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_plus_over_threshold_v0_5bx;   //!
   TBranch        *b_L1Tech_ZDC_tight_vertex_v0;   //!
   TBranch        *b_L1Tech_ZDC_tight_vertex_v0_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_tight_vertex_v0_5bx;   //!

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
   void fillNameArray(std::string * nameArray);

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

   void saveHistograms(const TString & name);
   void saveHistogram(TH1F * histo);
   TLorentzVector fromPtEtaPhiToPxPyPz( const double & pt, const double & eta, const double & phi );

   std::map<TString, TH1*> histoMap_;
   TString dir_;
};

#endif

#ifdef checkOpenHLT_cxx
checkOpenHLT::checkOpenHLT(TTree *tree) :
  dir_("plots/")
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

   fChain->SetBranchAddress("NrecoJetCal", &NrecoJetCal, &b_NrecoJetCal);
   fChain->SetBranchAddress("NrecoJetGen", &NrecoJetGen, &b_NrecoJetGen);
   fChain->SetBranchAddress("NrecoTowCal", &NrecoTowCal, &b_NrecoTowCal);
   fChain->SetBranchAddress("recoJetCalPt", recoJetCalPt, &b_recoJetCalPt);
   fChain->SetBranchAddress("recoJetCalPhi", recoJetCalPhi, &b_recoJetCalPhi);
   fChain->SetBranchAddress("recoJetCalEta", recoJetCalEta, &b_recoJetCalEta);
   fChain->SetBranchAddress("recoJetCalE", recoJetCalE, &b_recoJetCalE);
   fChain->SetBranchAddress("recoJetCalEMF", recoJetCalEMF, &b_recoJetCalEMF);
   fChain->SetBranchAddress("recoJetCalN90", recoJetCalN90, &b_recoJetCalN90);
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
   fChain->SetBranchAddress("NrecoJetCorCal", &NrecoJetCorCal, &b_NrecoJetCorCal);
   fChain->SetBranchAddress("recoJetCorCalPt", recoJetCorCalPt, &b_recoJetCorCalPt);
   fChain->SetBranchAddress("recoJetCorCalPhi", recoJetCorCalPhi, &b_recoJetCorCalPhi);
   fChain->SetBranchAddress("recoJetCorCalEta", recoJetCorCalEta, &b_recoJetCorCalEta);
   fChain->SetBranchAddress("recoJetCorCalE", recoJetCorCalE, &b_recoJetCorCalE);
   fChain->SetBranchAddress("recoJetCorCalEMF", recoJetCorCalEMF, &b_recoJetCorCalEMF);
   fChain->SetBranchAddress("recoJetCorCalN90", recoJetCorCalN90, &b_recoJetCorCalN90);
   fChain->SetBranchAddress("NohTau", &NohTau, &b_NohTau);
   fChain->SetBranchAddress("ohTauEta", ohTauEta, &b_ohTauEta);
   fChain->SetBranchAddress("ohTauPhi", ohTauPhi, &b_ohTauPhi);
   fChain->SetBranchAddress("ohTauPt", ohTauPt, &b_ohTauPt);
   fChain->SetBranchAddress("ohTauEiso", ohTauEiso, &b_ohTauEiso);
   fChain->SetBranchAddress("ohTauL25Tpt", ohTauL25Tpt, &b_ohTauL25Tpt);
   fChain->SetBranchAddress("ohTauL3Tiso", ohTauL3Tiso, &b_ohTauL3Tiso);
   fChain->SetBranchAddress("NohpfTau", &NohpfTau, &b_NohpfTau);
   fChain->SetBranchAddress("ohpfTauPt", ohpfTauPt, &b_ohpfTauPt);
   fChain->SetBranchAddress("ohpfTauEta", ohpfTauEta, &b_ohpfTauEta);
   fChain->SetBranchAddress("ohpfTauPhi", ohpfTauPhi, &b_ohpfTauPhi);
   fChain->SetBranchAddress("ohpfTauLeadTrackPt", ohpfTauLeadTrackPt, &b_ohpfTauLeadTrackPt);
   fChain->SetBranchAddress("ohpfTauLeadPionPt", ohpfTauLeadPionPt, &b_ohpfTauLeadPionPt);
   fChain->SetBranchAddress("ohpfTauTrkIso", ohpfTauTrkIso, &b_ohpfTauTrkIso);
   fChain->SetBranchAddress("ohpfTauGammaIso", ohpfTauGammaIso, &b_ohpfTauGammaIso);
   fChain->SetBranchAddress("ohpfTauJetPt", ohpfTauJetPt, &b_ohpfTauJetPt);
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
   fChain->SetBranchAddress("pfMHT", &pfMHT, &b_pfMHT);
   fChain->SetBranchAddress("NohPFJet", &NohPFJet, &b_NohPFJet);
   fChain->SetBranchAddress("pfJetPt", pfJetPt, &b_pfJetPt);
   fChain->SetBranchAddress("pfJetEta", pfJetEta, &b_pfJetEta);
   fChain->SetBranchAddress("pfJetPhi", pfJetPhi, &b_pfJetPhi);
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
   fChain->SetBranchAddress("ohBJetIPLooseL25Tag", ohBJetIPLooseL25Tag, &b_ohBJetIPLooseL25Tag);
   fChain->SetBranchAddress("ohBJetIPLooseL3Tag", ohBJetIPLooseL3Tag, &b_ohBJetIPLooseL3Tag);
   fChain->SetBranchAddress("ohBJetMuL25Tag", ohBJetMuL25Tag, &b_ohBJetMuL25Tag);
   fChain->SetBranchAddress("ohBJetMuL3Tag", ohBJetMuL3Tag, &b_ohBJetMuL3Tag);
   fChain->SetBranchAddress("ohBJetPerfL25Tag", ohBJetPerfL25Tag, &b_ohBJetPerfL25Tag);
   fChain->SetBranchAddress("ohBJetPerfL3Tag", ohBJetPerfL3Tag, &b_ohBJetPerfL3Tag);
   fChain->SetBranchAddress("NrecoElec", &NrecoElec, &b_NrecoElec);
   fChain->SetBranchAddress("recoElecPt", &recoElecPt, &b_recoElecPt);
   fChain->SetBranchAddress("recoElecPhi", &recoElecPhi, &b_recoElecPhi);
   fChain->SetBranchAddress("recoElecEta", &recoElecEta, &b_recoElecEta);
   fChain->SetBranchAddress("recoElecEt", &recoElecEt, &b_recoElecEt);
   fChain->SetBranchAddress("recoElecE", &recoElecE, &b_recoElecE);
   fChain->SetBranchAddress("recoElecEleID", &recoElecEleID, &b_recoElecEleID);
   fChain->SetBranchAddress("NrecoPhot", &NrecoPhot, &b_NrecoPhot);
   fChain->SetBranchAddress("recoPhotPt", &recoPhotPt, &b_recoPhotPt);
   fChain->SetBranchAddress("recoPhotPhi", &recoPhotPhi, &b_recoPhotPhi);
   fChain->SetBranchAddress("recoPhotEta", &recoPhotEta, &b_recoPhotEta);
   fChain->SetBranchAddress("recoPhotEt", &recoPhotEt, &b_recoPhotEt);
   fChain->SetBranchAddress("recoPhotE", &recoPhotE, &b_recoPhotE);
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
   fChain->SetBranchAddress("NohEle", &NohEle, &b_NohEle);
   fChain->SetBranchAddress("ohEleEt", ohEleEt, &b_ohEleEt);
   fChain->SetBranchAddress("ohEleEta", ohEleEta, &b_ohEleEta);
   fChain->SetBranchAddress("ohElePhi", ohElePhi, &b_ohElePhi);
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
   fChain->SetBranchAddress("ohMuL2L1idx", ohMuL2L1idx, &b_ohMuL2L1idx);
   fChain->SetBranchAddress("NohMuL3", &NohMuL3, &b_NohMuL3);
   fChain->SetBranchAddress("ohMuL3Pt", ohMuL3Pt, &b_ohMuL3Pt);
   fChain->SetBranchAddress("ohMuL3Phi", ohMuL3Phi, &b_ohMuL3Phi);
   fChain->SetBranchAddress("ohMuL3Eta", ohMuL3Eta, &b_ohMuL3Eta);
   fChain->SetBranchAddress("ohMuL3Chg", ohMuL3Chg, &b_ohMuL3Chg);
   fChain->SetBranchAddress("ohMuL3PtErr", ohMuL3PtErr, &b_ohMuL3PtErr);
   fChain->SetBranchAddress("ohMuL3Iso", ohMuL3Iso, &b_ohMuL3Iso);
   fChain->SetBranchAddress("ohMuL3Dr", ohMuL3Dr, &b_ohMuL3Dr);
   fChain->SetBranchAddress("ohMuL3Dz", ohMuL3Dz, &b_ohMuL3Dz);
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
   fChain->SetBranchAddress("ohMuL2NoVtxL1idx", ohMuL2NoVtxL1idx, &b_ohMuL2NoVtxL1idx);
   fChain->SetBranchAddress("ohMuL2NoVtxNhits", ohMuL2NoVtxNhits, &b_ohMuL2NoVtxNhits);
   fChain->SetBranchAddress("ohMuL2NoVtxNchambers", ohMuL2NoVtxNchambers, &b_ohMuL2NoVtxNchambers);
   fChain->SetBranchAddress("ohHighestEnergyEERecHit", &ohHighestEnergyEERecHit, &b_ohHighestEnergyEERecHit);
   fChain->SetBranchAddress("ohHighestEnergyEBRecHit", &ohHighestEnergyEBRecHit, &b_ohHighestEnergyEBRecHit);
   fChain->SetBranchAddress("ohHighestEnergyHBHERecHit", &ohHighestEnergyHBHERecHit, &b_ohHighestEnergyHBHERecHit);
   fChain->SetBranchAddress("ohHighestEnergyHORecHit", &ohHighestEnergyHORecHit, &b_ohHighestEnergyHORecHit);
   fChain->SetBranchAddress("ohHighestEnergyHFRecHit", &ohHighestEnergyHFRecHit, &b_ohHighestEnergyHFRecHit);
   fChain->SetBranchAddress("Nalcapi0clusters", &Nalcapi0clusters, &b_Nalcapi0clusters);
   fChain->SetBranchAddress("ohAlcapi0ptClusAll", ohAlcapi0ptClusAll, &b_ohAlcapi0ptClusAll);
   fChain->SetBranchAddress("ohAlcapi0etaClusAll", ohAlcapi0etaClusAll, &b_ohAlcapi0etaClusAll);
   fChain->SetBranchAddress("ohAlcapi0phiClusAll", ohAlcapi0phiClusAll, &b_ohAlcapi0phiClusAll);
   fChain->SetBranchAddress("ohAlcapi0s4s9ClusAll", ohAlcapi0s4s9ClusAll, &b_ohAlcapi0s4s9ClusAll);
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
   fChain->SetBranchAddress("ohPixelTracksL3Pt", ohPixelTracksL3Pt, &b_ohPixelTracksL3Pt);
   fChain->SetBranchAddress("ohPixelTracksL3Eta", ohPixelTracksL3Eta, &b_ohPixelTracksL3Eta);
   fChain->SetBranchAddress("ohPixelTracksL3Phi", ohPixelTracksL3Phi, &b_ohPixelTracksL3Phi);
   fChain->SetBranchAddress("ohPixelTracksL3Vz", ohPixelTracksL3Vz, &b_ohPixelTracksL3Vz);
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
   fChain->SetBranchAddress("recoNVrt", &recoNVrt, &b_NVrtx);
   fChain->SetBranchAddress("recoVrtX", recoVrtX, &b_recoVrtX);
   fChain->SetBranchAddress("recoVrtY", recoVrtY, &b_recoVrtY);
   fChain->SetBranchAddress("recoVrtZ", recoVrtZ, &b_recoVrtZ);
   fChain->SetBranchAddress("recoVrtNtrk", recoVrtNtrk, &b_recoVrtNtrk);
   fChain->SetBranchAddress("recoVrtChi2", recoVrtChi2, &b_recoVrtChi2);
   fChain->SetBranchAddress("recoVrtNdof", recoVrtNdof, &b_recoVrtNdof);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Bx", &Bx, &b_Bx);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("HLT_Activity_CSC", &HLT_Activity_CSC, &b_HLT_Activity_CSC);
   fChain->SetBranchAddress("HLT_Activity_CSC_Prescl", &HLT_Activity_CSC_Prescl, &b_HLT_Activity_CSC_Prescl);
   fChain->SetBranchAddress("HLT_Activity_DT", &HLT_Activity_DT, &b_HLT_Activity_DT);
   fChain->SetBranchAddress("HLT_Activity_DT_Prescl", &HLT_Activity_DT_Prescl, &b_HLT_Activity_DT_Prescl);
   fChain->SetBranchAddress("HLT_Activity_DT_Tuned", &HLT_Activity_DT_Tuned, &b_HLT_Activity_DT_Tuned);
   fChain->SetBranchAddress("HLT_Activity_DT_Tuned_Prescl", &HLT_Activity_DT_Tuned_Prescl, &b_HLT_Activity_DT_Tuned_Prescl);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC7", &HLT_Activity_Ecal_SC7, &b_HLT_Activity_Ecal_SC7);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC7_Prescl", &HLT_Activity_Ecal_SC7_Prescl, &b_HLT_Activity_Ecal_SC7_Prescl);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC17", &HLT_Activity_Ecal_SC17, &b_HLT_Activity_Ecal_SC17);
   fChain->SetBranchAddress("HLT_Activity_Ecal_SC17_Prescl", &HLT_Activity_Ecal_SC17_Prescl, &b_HLT_Activity_Ecal_SC17_Prescl);
   fChain->SetBranchAddress("HLT_L1Jet6U", &HLT_L1Jet6U, &b_HLT_L1Jet6U);
   fChain->SetBranchAddress("HLT_L1Jet6U_Prescl", &HLT_L1Jet6U_Prescl, &b_HLT_L1Jet6U_Prescl);
   fChain->SetBranchAddress("HLT_L1Jet10U", &HLT_L1Jet10U, &b_HLT_L1Jet10U);
   fChain->SetBranchAddress("HLT_L1Jet10U_Prescl", &HLT_L1Jet10U_Prescl, &b_HLT_L1Jet10U_Prescl);
   fChain->SetBranchAddress("HLT_Jet15U", &HLT_Jet15U, &b_HLT_Jet15U);
   fChain->SetBranchAddress("HLT_Jet15U_Prescl", &HLT_Jet15U_Prescl, &b_HLT_Jet15U_Prescl);
   fChain->SetBranchAddress("HLT_Jet15U_HcalNoiseFiltered", &HLT_Jet15U_HcalNoiseFiltered, &b_HLT_Jet15U_HcalNoiseFiltered);
   fChain->SetBranchAddress("HLT_Jet15U_HcalNoiseFiltered_Prescl", &HLT_Jet15U_HcalNoiseFiltered_Prescl, &b_HLT_Jet15U_HcalNoiseFiltered_Prescl);
   fChain->SetBranchAddress("HLT_Jet30U", &HLT_Jet30U, &b_HLT_Jet30U);
   fChain->SetBranchAddress("HLT_Jet30U_Prescl", &HLT_Jet30U_Prescl, &b_HLT_Jet30U_Prescl);
   fChain->SetBranchAddress("HLT_Jet50U", &HLT_Jet50U, &b_HLT_Jet50U);
   fChain->SetBranchAddress("HLT_Jet50U_Prescl", &HLT_Jet50U_Prescl, &b_HLT_Jet50U_Prescl);
   fChain->SetBranchAddress("HLT_Jet70U_v2", &HLT_Jet70U_v2, &b_HLT_Jet70U_v2);
   fChain->SetBranchAddress("HLT_Jet70U_v2_Prescl", &HLT_Jet70U_v2_Prescl, &b_HLT_Jet70U_v2_Prescl);
   fChain->SetBranchAddress("HLT_Jet100U_v2", &HLT_Jet100U_v2, &b_HLT_Jet100U_v2);
   fChain->SetBranchAddress("HLT_Jet100U_v2_Prescl", &HLT_Jet100U_v2_Prescl, &b_HLT_Jet100U_v2_Prescl);
   fChain->SetBranchAddress("HLT_Jet140U_v1", &HLT_Jet140U_v1, &b_HLT_Jet140U_v1);
   fChain->SetBranchAddress("HLT_Jet140U_v1_Prescl", &HLT_Jet140U_v1_Prescl, &b_HLT_Jet140U_v1_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve15U", &HLT_DiJetAve15U, &b_HLT_DiJetAve15U);
   fChain->SetBranchAddress("HLT_DiJetAve15U_Prescl", &HLT_DiJetAve15U_Prescl, &b_HLT_DiJetAve15U_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve30U", &HLT_DiJetAve30U, &b_HLT_DiJetAve30U);
   fChain->SetBranchAddress("HLT_DiJetAve30U_Prescl", &HLT_DiJetAve30U_Prescl, &b_HLT_DiJetAve30U_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve50U", &HLT_DiJetAve50U, &b_HLT_DiJetAve50U);
   fChain->SetBranchAddress("HLT_DiJetAve50U_Prescl", &HLT_DiJetAve50U_Prescl, &b_HLT_DiJetAve50U_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve70U_v2", &HLT_DiJetAve70U_v2, &b_HLT_DiJetAve70U_v2);
   fChain->SetBranchAddress("HLT_DiJetAve70U_v2_Prescl", &HLT_DiJetAve70U_v2_Prescl, &b_HLT_DiJetAve70U_v2_Prescl);
   fChain->SetBranchAddress("HLT_DiJetAve100U_v1", &HLT_DiJetAve100U_v1, &b_HLT_DiJetAve100U_v1);
   fChain->SetBranchAddress("HLT_DiJetAve100U_v1_Prescl", &HLT_DiJetAve100U_v1_Prescl, &b_HLT_DiJetAve100U_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleJet15U_ForwardBackward", &HLT_DoubleJet15U_ForwardBackward, &b_HLT_DoubleJet15U_ForwardBackward);
   fChain->SetBranchAddress("HLT_DoubleJet15U_ForwardBackward_Prescl", &HLT_DoubleJet15U_ForwardBackward_Prescl, &b_HLT_DoubleJet15U_ForwardBackward_Prescl);
   fChain->SetBranchAddress("HLT_DoubleJet25U_ForwardBackward", &HLT_DoubleJet25U_ForwardBackward, &b_HLT_DoubleJet25U_ForwardBackward);
   fChain->SetBranchAddress("HLT_DoubleJet25U_ForwardBackward_Prescl", &HLT_DoubleJet25U_ForwardBackward_Prescl, &b_HLT_DoubleJet25U_ForwardBackward_Prescl);
   fChain->SetBranchAddress("HLT_ExclDiJet30U_HFAND_v1", &HLT_ExclDiJet30U_HFAND_v1, &b_HLT_ExclDiJet30U_HFAND_v1);
   fChain->SetBranchAddress("HLT_ExclDiJet30U_HFAND_v1_Prescl", &HLT_ExclDiJet30U_HFAND_v1_Prescl, &b_HLT_ExclDiJet30U_HFAND_v1_Prescl);
   fChain->SetBranchAddress("HLT_ExclDiJet30U_HFOR_v1", &HLT_ExclDiJet30U_HFOR_v1, &b_HLT_ExclDiJet30U_HFOR_v1);
   fChain->SetBranchAddress("HLT_ExclDiJet30U_HFOR_v1_Prescl", &HLT_ExclDiJet30U_HFOR_v1_Prescl, &b_HLT_ExclDiJet30U_HFOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet15U_v2", &HLT_QuadJet15U_v2, &b_HLT_QuadJet15U_v2);
   fChain->SetBranchAddress("HLT_QuadJet15U_v2_Prescl", &HLT_QuadJet15U_v2_Prescl, &b_HLT_QuadJet15U_v2_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet20U_v2", &HLT_QuadJet20U_v2, &b_HLT_QuadJet20U_v2);
   fChain->SetBranchAddress("HLT_QuadJet20U_v2_Prescl", &HLT_QuadJet20U_v2_Prescl, &b_HLT_QuadJet20U_v2_Prescl);
   fChain->SetBranchAddress("HLT_QuadJet25U_v2", &HLT_QuadJet25U_v2, &b_HLT_QuadJet25U_v2);
   fChain->SetBranchAddress("HLT_QuadJet25U_v2_Prescl", &HLT_QuadJet25U_v2_Prescl, &b_HLT_QuadJet25U_v2_Prescl);
   fChain->SetBranchAddress("HLT_L1ETT100", &HLT_L1ETT100, &b_HLT_L1ETT100);
   fChain->SetBranchAddress("HLT_L1ETT100_Prescl", &HLT_L1ETT100_Prescl, &b_HLT_L1ETT100_Prescl);
   fChain->SetBranchAddress("HLT_L1ETT140_v1", &HLT_L1ETT140_v1, &b_HLT_L1ETT140_v1);
   fChain->SetBranchAddress("HLT_L1ETT140_v1_Prescl", &HLT_L1ETT140_v1_Prescl, &b_HLT_L1ETT140_v1_Prescl);
   fChain->SetBranchAddress("HLT_EcalOnly_SumEt160_v2", &HLT_EcalOnly_SumEt160_v2, &b_HLT_EcalOnly_SumEt160_v2);
   fChain->SetBranchAddress("HLT_EcalOnly_SumEt160_v2_Prescl", &HLT_EcalOnly_SumEt160_v2_Prescl, &b_HLT_EcalOnly_SumEt160_v2_Prescl);
   fChain->SetBranchAddress("HLT_L1MET20", &HLT_L1MET20, &b_HLT_L1MET20);
   fChain->SetBranchAddress("HLT_L1MET20_Prescl", &HLT_L1MET20_Prescl, &b_HLT_L1MET20_Prescl);
   fChain->SetBranchAddress("HLT_MET45", &HLT_MET45, &b_HLT_MET45);
   fChain->SetBranchAddress("HLT_MET45_Prescl", &HLT_MET45_Prescl, &b_HLT_MET45_Prescl);
   fChain->SetBranchAddress("HLT_MET45_HT100U_v1", &HLT_MET45_HT100U_v1, &b_HLT_MET45_HT100U_v1);
   fChain->SetBranchAddress("HLT_MET45_HT100U_v1_Prescl", &HLT_MET45_HT100U_v1_Prescl, &b_HLT_MET45_HT100U_v1_Prescl);
   fChain->SetBranchAddress("HLT_MET45_HT120U_v1", &HLT_MET45_HT120U_v1, &b_HLT_MET45_HT120U_v1);
   fChain->SetBranchAddress("HLT_MET45_HT120U_v1_Prescl", &HLT_MET45_HT120U_v1_Prescl, &b_HLT_MET45_HT120U_v1_Prescl);
   fChain->SetBranchAddress("HLT_MET65", &HLT_MET65, &b_HLT_MET65);
   fChain->SetBranchAddress("HLT_MET65_Prescl", &HLT_MET65_Prescl, &b_HLT_MET65_Prescl);
   fChain->SetBranchAddress("HLT_MET80_v1", &HLT_MET80_v1, &b_HLT_MET80_v1);
   fChain->SetBranchAddress("HLT_MET80_v1_Prescl", &HLT_MET80_v1_Prescl, &b_HLT_MET80_v1_Prescl);
   fChain->SetBranchAddress("HLT_MET100_v2", &HLT_MET100_v2, &b_HLT_MET100_v2);
   fChain->SetBranchAddress("HLT_MET100_v2_Prescl", &HLT_MET100_v2_Prescl, &b_HLT_MET100_v2_Prescl);
   fChain->SetBranchAddress("HLT_HT50U_v1", &HLT_HT50U_v1, &b_HLT_HT50U_v1);
   fChain->SetBranchAddress("HLT_HT50U_v1_Prescl", &HLT_HT50U_v1_Prescl, &b_HLT_HT50U_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT100U", &HLT_HT100U, &b_HLT_HT100U);
   fChain->SetBranchAddress("HLT_HT100U_Prescl", &HLT_HT100U_Prescl, &b_HLT_HT100U_Prescl);
   fChain->SetBranchAddress("HLT_HT120U", &HLT_HT120U, &b_HLT_HT120U);
   fChain->SetBranchAddress("HLT_HT120U_Prescl", &HLT_HT120U_Prescl, &b_HLT_HT120U_Prescl);
   fChain->SetBranchAddress("HLT_HT140U", &HLT_HT140U, &b_HLT_HT140U);
   fChain->SetBranchAddress("HLT_HT140U_Prescl", &HLT_HT140U_Prescl, &b_HLT_HT140U_Prescl);
   fChain->SetBranchAddress("HLT_HT140U_Eta3_v1", &HLT_HT140U_Eta3_v1, &b_HLT_HT140U_Eta3_v1);
   fChain->SetBranchAddress("HLT_HT140U_Eta3_v1_Prescl", &HLT_HT140U_Eta3_v1_Prescl, &b_HLT_HT140U_Eta3_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT160U_v1", &HLT_HT160U_v1, &b_HLT_HT160U_v1);
   fChain->SetBranchAddress("HLT_HT160U_v1_Prescl", &HLT_HT160U_v1_Prescl, &b_HLT_HT160U_v1_Prescl);
   fChain->SetBranchAddress("HLT_HT200U_v1", &HLT_HT200U_v1, &b_HLT_HT200U_v1);
   fChain->SetBranchAddress("HLT_HT200U_v1_Prescl", &HLT_HT200U_v1_Prescl, &b_HLT_HT200U_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MuOpen", &HLT_L1MuOpen, &b_HLT_L1MuOpen);
   fChain->SetBranchAddress("HLT_L1MuOpen_Prescl", &HLT_L1MuOpen_Prescl, &b_HLT_L1MuOpen_Prescl);
   fChain->SetBranchAddress("HLT_L1MuOpen_DT", &HLT_L1MuOpen_DT, &b_HLT_L1MuOpen_DT);
   fChain->SetBranchAddress("HLT_L1MuOpen_DT_Prescl", &HLT_L1MuOpen_DT_Prescl, &b_HLT_L1MuOpen_DT_Prescl);
   fChain->SetBranchAddress("HLT_L1MuOpen_AntiBPTX", &HLT_L1MuOpen_AntiBPTX, &b_HLT_L1MuOpen_AntiBPTX);
   fChain->SetBranchAddress("HLT_L1MuOpen_AntiBPTX_Prescl", &HLT_L1MuOpen_AntiBPTX_Prescl, &b_HLT_L1MuOpen_AntiBPTX_Prescl);
   fChain->SetBranchAddress("HLT_L1Mu7_v1", &HLT_L1Mu7_v1, &b_HLT_L1Mu7_v1);
   fChain->SetBranchAddress("HLT_L1Mu7_v1_Prescl", &HLT_L1Mu7_v1_Prescl, &b_HLT_L1Mu7_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Mu20", &HLT_L1Mu20, &b_HLT_L1Mu20);
   fChain->SetBranchAddress("HLT_L1Mu20_Prescl", &HLT_L1Mu20_Prescl, &b_HLT_L1Mu20_Prescl);
   fChain->SetBranchAddress("HLT_L2Mu0_NoVertex", &HLT_L2Mu0_NoVertex, &b_HLT_L2Mu0_NoVertex);
   fChain->SetBranchAddress("HLT_L2Mu0_NoVertex_Prescl", &HLT_L2Mu0_NoVertex_Prescl, &b_HLT_L2Mu0_NoVertex_Prescl);
   fChain->SetBranchAddress("HLT_L2Mu7_v1", &HLT_L2Mu7_v1, &b_HLT_L2Mu7_v1);
   fChain->SetBranchAddress("HLT_L2Mu7_v1_Prescl", &HLT_L2Mu7_v1_Prescl, &b_HLT_L2Mu7_v1_Prescl);
   fChain->SetBranchAddress("HLT_L2Mu30_v1", &HLT_L2Mu30_v1, &b_HLT_L2Mu30_v1);
   fChain->SetBranchAddress("HLT_L2Mu30_v1_Prescl", &HLT_L2Mu30_v1_Prescl, &b_HLT_L2Mu30_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu3_Prescl", &HLT_Mu3_Prescl, &b_HLT_Mu3_Prescl);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu5_Prescl", &HLT_Mu5_Prescl, &b_HLT_Mu5_Prescl);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu7_Prescl", &HLT_Mu7_Prescl, &b_HLT_Mu7_Prescl);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu9_Prescl", &HLT_Mu9_Prescl, &b_HLT_Mu9_Prescl);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   fChain->SetBranchAddress("HLT_Mu11_Prescl", &HLT_Mu11_Prescl, &b_HLT_Mu11_Prescl);
   fChain->SetBranchAddress("HLT_Mu13_v1", &HLT_Mu13_v1, &b_HLT_Mu13_v1);
   fChain->SetBranchAddress("HLT_Mu13_v1_Prescl", &HLT_Mu13_v1_Prescl, &b_HLT_Mu13_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu15_v1", &HLT_Mu15_v1, &b_HLT_Mu15_v1);
   fChain->SetBranchAddress("HLT_Mu15_v1_Prescl", &HLT_Mu15_v1_Prescl, &b_HLT_Mu15_v1_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu9", &HLT_IsoMu9, &b_HLT_IsoMu9);
   fChain->SetBranchAddress("HLT_IsoMu9_Prescl", &HLT_IsoMu9_Prescl, &b_HLT_IsoMu9_Prescl);
   fChain->SetBranchAddress("HLT_IsoMu11_v1", &HLT_IsoMu11_v1, &b_HLT_IsoMu11_v1);
   fChain->SetBranchAddress("HLT_IsoMu11_v1_Prescl", &HLT_IsoMu11_v1_Prescl, &b_HLT_IsoMu11_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu20_NoVertex", &HLT_Mu20_NoVertex, &b_HLT_Mu20_NoVertex);
   fChain->SetBranchAddress("HLT_Mu20_NoVertex_Prescl", &HLT_Mu20_NoVertex_Prescl, &b_HLT_Mu20_NoVertex_Prescl);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, &b_HLT_L1DoubleMuOpen);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen_Prescl", &HLT_L1DoubleMuOpen_Prescl, &b_HLT_L1DoubleMuOpen_Prescl);
   fChain->SetBranchAddress("HLT_L2DoubleMu0", &HLT_L2DoubleMu0, &b_HLT_L2DoubleMu0);
   fChain->SetBranchAddress("HLT_L2DoubleMu0_Prescl", &HLT_L2DoubleMu0_Prescl, &b_HLT_L2DoubleMu0_Prescl);
   fChain->SetBranchAddress("HLT_L2DoubleMu20_NoVertex_v1", &HLT_L2DoubleMu20_NoVertex_v1, &b_HLT_L2DoubleMu20_NoVertex_v1);
   fChain->SetBranchAddress("HLT_L2DoubleMu20_NoVertex_v1_Prescl", &HLT_L2DoubleMu20_NoVertex_v1_Prescl, &b_HLT_L2DoubleMu20_NoVertex_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &HLT_DoubleMu0_Quarkonium_v1, &b_HLT_DoubleMu0_Quarkonium_v1);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1_Prescl", &HLT_DoubleMu0_Quarkonium_v1_Prescl, &b_HLT_DoubleMu0_Quarkonium_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_LS_v1", &HLT_DoubleMu0_Quarkonium_LS_v1, &b_HLT_DoubleMu0_Quarkonium_LS_v1);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_LS_v1_Prescl", &HLT_DoubleMu0_Quarkonium_LS_v1_Prescl, &b_HLT_DoubleMu0_Quarkonium_LS_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
   fChain->SetBranchAddress("HLT_DoubleMu0_Prescl", &HLT_DoubleMu0_Prescl, &b_HLT_DoubleMu0_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu3_v2", &HLT_DoubleMu3_v2, &b_HLT_DoubleMu3_v2);
   fChain->SetBranchAddress("HLT_DoubleMu3_v2_Prescl", &HLT_DoubleMu3_v2_Prescl, &b_HLT_DoubleMu3_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleMu5_v1", &HLT_DoubleMu5_v1, &b_HLT_DoubleMu5_v1);
   fChain->SetBranchAddress("HLT_DoubleMu5_v1_Prescl", &HLT_DoubleMu5_v1_Prescl, &b_HLT_DoubleMu5_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu0", &HLT_Mu5_L2Mu0, &b_HLT_Mu5_L2Mu0);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu0_Prescl", &HLT_Mu5_L2Mu0_Prescl, &b_HLT_Mu5_L2Mu0_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi", &HLT_Mu3_Track3_Jpsi, &b_HLT_Mu3_Track3_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_Prescl", &HLT_Mu3_Track3_Jpsi_Prescl, &b_HLT_Mu3_Track3_Jpsi_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v1", &HLT_Mu3_Track5_Jpsi_v1, &b_HLT_Mu3_Track5_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v1_Prescl", &HLT_Mu3_Track5_Jpsi_v1_Prescl, &b_HLT_Mu3_Track5_Jpsi_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Track0_Jpsi", &HLT_Mu5_Track0_Jpsi, &b_HLT_Mu5_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_Track0_Jpsi_Prescl", &HLT_Mu5_Track0_Jpsi_Prescl, &b_HLT_Mu5_Track0_Jpsi_Prescl);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi", &HLT_Mu0_TkMu0_OST_Jpsi, &b_HLT_Mu0_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Prescl", &HLT_Mu0_TkMu0_OST_Jpsi_Prescl, &b_HLT_Mu0_TkMu0_OST_Jpsi_Prescl);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1_Prescl", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1_Prescl, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi", &HLT_Mu3_TkMu0_OST_Jpsi, &b_HLT_Mu3_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Prescl", &HLT_Mu3_TkMu0_OST_Jpsi_Prescl, &b_HLT_Mu3_TkMu0_OST_Jpsi_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi", &HLT_Mu5_TkMu0_OST_Jpsi, &b_HLT_Mu5_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Prescl", &HLT_Mu5_TkMu0_OST_Jpsi_Prescl, &b_HLT_Mu5_TkMu0_OST_Jpsi_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleEG2", &HLT_L1SingleEG2, &b_HLT_L1SingleEG2);
   fChain->SetBranchAddress("HLT_L1SingleEG2_Prescl", &HLT_L1SingleEG2_Prescl, &b_HLT_L1SingleEG2_Prescl);
   fChain->SetBranchAddress("HLT_L1SingleEG8", &HLT_L1SingleEG8, &b_HLT_L1SingleEG8);
   fChain->SetBranchAddress("HLT_L1SingleEG8_Prescl", &HLT_L1SingleEG8_Prescl, &b_HLT_L1SingleEG8_Prescl);
   fChain->SetBranchAddress("HLT_Ele10_SW_L1R", &HLT_Ele10_SW_L1R, &b_HLT_Ele10_SW_L1R);
   fChain->SetBranchAddress("HLT_Ele10_SW_L1R_Prescl", &HLT_Ele10_SW_L1R_Prescl, &b_HLT_Ele10_SW_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Ele12_SW_TightEleId_L1R", &HLT_Ele12_SW_TightEleId_L1R, &b_HLT_Ele12_SW_TightEleId_L1R);
   fChain->SetBranchAddress("HLT_Ele12_SW_TightEleId_L1R_Prescl", &HLT_Ele12_SW_TightEleId_L1R_Prescl, &b_HLT_Ele12_SW_TightEleId_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Ele12_SW_TighterEleId_L1R_v1", &HLT_Ele12_SW_TighterEleId_L1R_v1, &b_HLT_Ele12_SW_TighterEleId_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele12_SW_TighterEleId_L1R_v1_Prescl", &HLT_Ele12_SW_TighterEleId_L1R_v1_Prescl, &b_HLT_Ele12_SW_TighterEleId_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele12_SW_TighterEleIdIsol_L1R_v1", &HLT_Ele12_SW_TighterEleIdIsol_L1R_v1, &b_HLT_Ele12_SW_TighterEleIdIsol_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele12_SW_TighterEleIdIsol_L1R_v1_Prescl", &HLT_Ele12_SW_TighterEleIdIsol_L1R_v1_Prescl, &b_HLT_Ele12_SW_TighterEleIdIsol_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_L1R", &HLT_Ele17_SW_L1R, &b_HLT_Ele17_SW_L1R);
   fChain->SetBranchAddress("HLT_Ele17_SW_L1R_Prescl", &HLT_Ele17_SW_L1R_Prescl, &b_HLT_Ele17_SW_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightEleId_L1R", &HLT_Ele17_SW_TightEleId_L1R, &b_HLT_Ele17_SW_TightEleId_L1R);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightEleId_L1R_Prescl", &HLT_Ele17_SW_TightEleId_L1R_Prescl, &b_HLT_Ele17_SW_TightEleId_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_TighterEleId_L1R_v1", &HLT_Ele17_SW_TighterEleId_L1R_v1, &b_HLT_Ele17_SW_TighterEleId_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele17_SW_TighterEleId_L1R_v1_Prescl", &HLT_Ele17_SW_TighterEleId_L1R_v1_Prescl, &b_HLT_Ele17_SW_TighterEleId_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightEleIdIsol_L1R_v1", &HLT_Ele17_SW_TightEleIdIsol_L1R_v1, &b_HLT_Ele17_SW_TightEleIdIsol_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightEleIdIsol_L1R_v1_Prescl", &HLT_Ele17_SW_TightEleIdIsol_L1R_v1_Prescl, &b_HLT_Ele17_SW_TightEleIdIsol_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1", &HLT_Ele17_SW_TighterEleIdIsol_L1R_v1, &b_HLT_Ele17_SW_TighterEleIdIsol_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1_Prescl", &HLT_Ele17_SW_TighterEleIdIsol_L1R_v1_Prescl, &b_HLT_Ele17_SW_TighterEleIdIsol_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1", &HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1, &b_HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1_Prescl", &HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1_Prescl, &b_HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1", &HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1, &b_HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1_Prescl", &HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1_Prescl, &b_HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1", &HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1, &b_HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1_Prescl", &HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1_Prescl, &b_HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1", &HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1, &b_HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1);
   fChain->SetBranchAddress("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1_Prescl", &HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1_Prescl, &b_HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoubleEle4_SW_eeRes_L1R", &HLT_DoubleEle4_SW_eeRes_L1R, &b_HLT_DoubleEle4_SW_eeRes_L1R);
   fChain->SetBranchAddress("HLT_DoubleEle4_SW_eeRes_L1R_Prescl", &HLT_DoubleEle4_SW_eeRes_L1R_Prescl, &b_HLT_DoubleEle4_SW_eeRes_L1R_Prescl);
   fChain->SetBranchAddress("HLT_DoubleEle15_SW_L1R_v1", &HLT_DoubleEle15_SW_L1R_v1, &b_HLT_DoubleEle15_SW_L1R_v1);
   fChain->SetBranchAddress("HLT_DoubleEle15_SW_L1R_v1_Prescl", &HLT_DoubleEle15_SW_L1R_v1_Prescl, &b_HLT_DoubleEle15_SW_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon10_Cleaned_L1R", &HLT_Photon10_Cleaned_L1R, &b_HLT_Photon10_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon10_Cleaned_L1R_Prescl", &HLT_Photon10_Cleaned_L1R_Prescl, &b_HLT_Photon10_Cleaned_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Photon15_Cleaned_L1R", &HLT_Photon15_Cleaned_L1R, &b_HLT_Photon15_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon15_Cleaned_L1R_Prescl", &HLT_Photon15_Cleaned_L1R_Prescl, &b_HLT_Photon15_Cleaned_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Photon17_SC17HE_L1R_v1", &HLT_Photon17_SC17HE_L1R_v1, &b_HLT_Photon17_SC17HE_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon17_SC17HE_L1R_v1_Prescl", &HLT_Photon17_SC17HE_L1R_v1_Prescl, &b_HLT_Photon17_SC17HE_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon20_NoHE_L1R", &HLT_Photon20_NoHE_L1R, &b_HLT_Photon20_NoHE_L1R);
   fChain->SetBranchAddress("HLT_Photon20_NoHE_L1R_Prescl", &HLT_Photon20_NoHE_L1R_Prescl, &b_HLT_Photon20_NoHE_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Photon20_Cleaned_L1R", &HLT_Photon20_Cleaned_L1R, &b_HLT_Photon20_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon20_Cleaned_L1R_Prescl", &HLT_Photon20_Cleaned_L1R_Prescl, &b_HLT_Photon20_Cleaned_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Photon30_Cleaned_L1R", &HLT_Photon30_Cleaned_L1R, &b_HLT_Photon30_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon30_Cleaned_L1R_Prescl", &HLT_Photon30_Cleaned_L1R_Prescl, &b_HLT_Photon30_Cleaned_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1", &HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1, &b_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_Prescl", &HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_Prescl, &b_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon35_Isol_Cleaned_L1R_v1", &HLT_Photon35_Isol_Cleaned_L1R_v1, &b_HLT_Photon35_Isol_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon35_Isol_Cleaned_L1R_v1_Prescl", &HLT_Photon35_Isol_Cleaned_L1R_v1_Prescl, &b_HLT_Photon35_Isol_Cleaned_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon50_Cleaned_L1R_v1", &HLT_Photon50_Cleaned_L1R_v1, &b_HLT_Photon50_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon50_Cleaned_L1R_v1_Prescl", &HLT_Photon50_Cleaned_L1R_v1_Prescl, &b_HLT_Photon50_Cleaned_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon50_NoHE_L1R", &HLT_Photon50_NoHE_L1R, &b_HLT_Photon50_NoHE_L1R);
   fChain->SetBranchAddress("HLT_Photon50_NoHE_L1R_Prescl", &HLT_Photon50_NoHE_L1R_Prescl, &b_HLT_Photon50_NoHE_L1R_Prescl);
   fChain->SetBranchAddress("HLT_Photon70_NoHE_Cleaned_L1R_v1", &HLT_Photon70_NoHE_Cleaned_L1R_v1, &b_HLT_Photon70_NoHE_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon70_NoHE_Cleaned_L1R_v1_Prescl", &HLT_Photon70_NoHE_Cleaned_L1R_v1_Prescl, &b_HLT_Photon70_NoHE_Cleaned_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Photon100_NoHE_Cleaned_L1R_v1", &HLT_Photon100_NoHE_Cleaned_L1R_v1, &b_HLT_Photon100_NoHE_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon100_NoHE_Cleaned_L1R_v1_Prescl", &HLT_Photon100_NoHE_Cleaned_L1R_v1_Prescl, &b_HLT_Photon100_NoHE_Cleaned_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_DoublePhoton5_CEP_L1R", &HLT_DoublePhoton5_CEP_L1R, &b_HLT_DoublePhoton5_CEP_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton5_CEP_L1R_Prescl", &HLT_DoublePhoton5_CEP_L1R_Prescl, &b_HLT_DoublePhoton5_CEP_L1R_Prescl);
   fChain->SetBranchAddress("HLT_DoublePhoton17_L1R", &HLT_DoublePhoton17_L1R, &b_HLT_DoublePhoton17_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton17_L1R_Prescl", &HLT_DoublePhoton17_L1R_Prescl, &b_HLT_DoublePhoton17_L1R_Prescl);
   fChain->SetBranchAddress("HLT_SingleIsoTau20_Trk5_MET20", &HLT_SingleIsoTau20_Trk5_MET20, &b_HLT_SingleIsoTau20_Trk5_MET20);
   fChain->SetBranchAddress("HLT_SingleIsoTau20_Trk5_MET20_Prescl", &HLT_SingleIsoTau20_Trk5_MET20_Prescl, &b_HLT_SingleIsoTau20_Trk5_MET20_Prescl);
   fChain->SetBranchAddress("HLT_SingleIsoTau20_Trk15_MET20", &HLT_SingleIsoTau20_Trk15_MET20, &b_HLT_SingleIsoTau20_Trk15_MET20);
   fChain->SetBranchAddress("HLT_SingleIsoTau20_Trk15_MET20_Prescl", &HLT_SingleIsoTau20_Trk15_MET20_Prescl, &b_HLT_SingleIsoTau20_Trk15_MET20_Prescl);
   fChain->SetBranchAddress("HLT_SingleIsoTau30_Trk5_MET20", &HLT_SingleIsoTau30_Trk5_MET20, &b_HLT_SingleIsoTau30_Trk5_MET20);
   fChain->SetBranchAddress("HLT_SingleIsoTau30_Trk5_MET20_Prescl", &HLT_SingleIsoTau30_Trk5_MET20_Prescl, &b_HLT_SingleIsoTau30_Trk5_MET20_Prescl);
   fChain->SetBranchAddress("HLT_SingleIsoTau30_Trk5_v2", &HLT_SingleIsoTau30_Trk5_v2, &b_HLT_SingleIsoTau30_Trk5_v2);
   fChain->SetBranchAddress("HLT_SingleIsoTau30_Trk5_v2_Prescl", &HLT_SingleIsoTau30_Trk5_v2_Prescl, &b_HLT_SingleIsoTau30_Trk5_v2_Prescl);
   fChain->SetBranchAddress("HLT_DoubleIsoTau15_OneLeg_Trk5", &HLT_DoubleIsoTau15_OneLeg_Trk5, &b_HLT_DoubleIsoTau15_OneLeg_Trk5);
   fChain->SetBranchAddress("HLT_DoubleIsoTau15_OneLeg_Trk5_Prescl", &HLT_DoubleIsoTau15_OneLeg_Trk5_Prescl, &b_HLT_DoubleIsoTau15_OneLeg_Trk5_Prescl);
   fChain->SetBranchAddress("HLT_DoubleIsoTau15_Trk5", &HLT_DoubleIsoTau15_Trk5, &b_HLT_DoubleIsoTau15_Trk5);
   fChain->SetBranchAddress("HLT_DoubleIsoTau15_Trk5_Prescl", &HLT_DoubleIsoTau15_Trk5_Prescl, &b_HLT_DoubleIsoTau15_Trk5_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet10U_v1", &HLT_BTagMu_DiJet10U_v1, &b_HLT_BTagMu_DiJet10U_v1);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet10U_v1_Prescl", &HLT_BTagMu_DiJet10U_v1_Prescl, &b_HLT_BTagMu_DiJet10U_v1_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20U_v1", &HLT_BTagMu_DiJet20U_v1, &b_HLT_BTagMu_DiJet20U_v1);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20U_v1_Prescl", &HLT_BTagMu_DiJet20U_v1_Prescl, &b_HLT_BTagMu_DiJet20U_v1_Prescl);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20U_Mu5_v1", &HLT_BTagMu_DiJet20U_Mu5_v1, &b_HLT_BTagMu_DiJet20U_Mu5_v1);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20U_Mu5_v1_Prescl", &HLT_BTagMu_DiJet20U_Mu5_v1_Prescl, &b_HLT_BTagMu_DiJet20U_Mu5_v1_Prescl);
   fChain->SetBranchAddress("HLT_StoppedHSCP20_v3", &HLT_StoppedHSCP20_v3, &b_HLT_StoppedHSCP20_v3);
   fChain->SetBranchAddress("HLT_StoppedHSCP20_v3_Prescl", &HLT_StoppedHSCP20_v3_Prescl, &b_HLT_StoppedHSCP20_v3_Prescl);
   fChain->SetBranchAddress("HLT_StoppedHSCP35_v3", &HLT_StoppedHSCP35_v3, &b_HLT_StoppedHSCP35_v3);
   fChain->SetBranchAddress("HLT_StoppedHSCP35_v3_Prescl", &HLT_StoppedHSCP35_v3_Prescl, &b_HLT_StoppedHSCP35_v3_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Photon11_Cleaned_L1R_v1", &HLT_Mu5_Photon11_Cleaned_L1R_v1, &b_HLT_Mu5_Photon11_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Mu5_Photon11_Cleaned_L1R_v1_Prescl", &HLT_Mu5_Photon11_Cleaned_L1R_v1_Prescl, &b_HLT_Mu5_Photon11_Cleaned_L1R_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Ele5_v1", &HLT_Mu5_Ele5_v1, &b_HLT_Mu5_Ele5_v1);
   fChain->SetBranchAddress("HLT_Mu5_Ele5_v1_Prescl", &HLT_Mu5_Ele5_v1_Prescl, &b_HLT_Mu5_Ele5_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Ele9_v1", &HLT_Mu5_Ele9_v1, &b_HLT_Mu5_Ele9_v1);
   fChain->SetBranchAddress("HLT_Mu5_Ele9_v1_Prescl", &HLT_Mu5_Ele9_v1_Prescl, &b_HLT_Mu5_Ele9_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Jet35U_v1", &HLT_Mu5_Jet35U_v1, &b_HLT_Mu5_Jet35U_v1);
   fChain->SetBranchAddress("HLT_Mu5_Jet35U_v1_Prescl", &HLT_Mu5_Jet35U_v1_Prescl, &b_HLT_Mu5_Jet35U_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_Jet50U_v2", &HLT_Mu5_Jet50U_v2, &b_HLT_Mu5_Jet50U_v2);
   fChain->SetBranchAddress("HLT_Mu5_Jet50U_v2_Prescl", &HLT_Mu5_Jet50U_v2_Prescl, &b_HLT_Mu5_Jet50U_v2_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_MET45_v1", &HLT_Mu5_MET45_v1, &b_HLT_Mu5_MET45_v1);
   fChain->SetBranchAddress("HLT_Mu5_MET45_v1_Prescl", &HLT_Mu5_MET45_v1_Prescl, &b_HLT_Mu5_MET45_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_HT50U_v1", &HLT_Mu5_HT50U_v1, &b_HLT_Mu5_HT50U_v1);
   fChain->SetBranchAddress("HLT_Mu5_HT50U_v1_Prescl", &HLT_Mu5_HT50U_v1_Prescl, &b_HLT_Mu5_HT50U_v1_Prescl);
   fChain->SetBranchAddress("HLT_Mu5_HT70U_v1", &HLT_Mu5_HT70U_v1, &b_HLT_Mu5_HT70U_v1);
   fChain->SetBranchAddress("HLT_Mu5_HT70U_v1_Prescl", &HLT_Mu5_HT70U_v1_Prescl, &b_HLT_Mu5_HT70U_v1_Prescl);
   fChain->SetBranchAddress("HLT_Ele10_MET45_v1", &HLT_Ele10_MET45_v1, &b_HLT_Ele10_MET45_v1);
   fChain->SetBranchAddress("HLT_Ele10_MET45_v1_Prescl", &HLT_Ele10_MET45_v1_Prescl, &b_HLT_Ele10_MET45_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_ZeroBias_Prescl", &HLT_ZeroBias_Prescl, &b_HLT_ZeroBias_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBiasPixel_SingleTrack", &HLT_ZeroBiasPixel_SingleTrack, &b_HLT_ZeroBiasPixel_SingleTrack);
   fChain->SetBranchAddress("HLT_ZeroBiasPixel_SingleTrack_Prescl", &HLT_ZeroBiasPixel_SingleTrack_Prescl, &b_HLT_ZeroBiasPixel_SingleTrack_Prescl);
   fChain->SetBranchAddress("HLT_MinBiasPixel_SingleTrack", &HLT_MinBiasPixel_SingleTrack, &b_HLT_MinBiasPixel_SingleTrack);
   fChain->SetBranchAddress("HLT_MinBiasPixel_SingleTrack_Prescl", &HLT_MinBiasPixel_SingleTrack_Prescl, &b_HLT_MinBiasPixel_SingleTrack_Prescl);
   fChain->SetBranchAddress("HLT_MultiVertex6", &HLT_MultiVertex6, &b_HLT_MultiVertex6);
   fChain->SetBranchAddress("HLT_MultiVertex6_Prescl", &HLT_MultiVertex6_Prescl, &b_HLT_MultiVertex6_Prescl);
   fChain->SetBranchAddress("HLT_MultiVertex8_L1ETT60", &HLT_MultiVertex8_L1ETT60, &b_HLT_MultiVertex8_L1ETT60);
   fChain->SetBranchAddress("HLT_MultiVertex8_L1ETT60_Prescl", &HLT_MultiVertex8_L1ETT60_Prescl, &b_HLT_MultiVertex8_L1ETT60_Prescl);
   fChain->SetBranchAddress("HLT_L1_BptxXOR_BscMinBiasOR", &HLT_L1_BptxXOR_BscMinBiasOR, &b_HLT_L1_BptxXOR_BscMinBiasOR);
   fChain->SetBranchAddress("HLT_L1_BptxXOR_BscMinBiasOR_Prescl", &HLT_L1_BptxXOR_BscMinBiasOR_Prescl, &b_HLT_L1_BptxXOR_BscMinBiasOR_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_minBias_OR", &HLT_L1Tech_BSC_minBias_OR, &b_HLT_L1Tech_BSC_minBias_OR);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_minBias_OR_Prescl", &HLT_L1Tech_BSC_minBias_OR_Prescl, &b_HLT_L1Tech_BSC_minBias_OR_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_minBias", &HLT_L1Tech_BSC_minBias, &b_HLT_L1Tech_BSC_minBias);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_minBias_Prescl", &HLT_L1Tech_BSC_minBias_Prescl, &b_HLT_L1Tech_BSC_minBias_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_halo", &HLT_L1Tech_BSC_halo, &b_HLT_L1Tech_BSC_halo);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_halo_Prescl", &HLT_L1Tech_BSC_halo_Prescl, &b_HLT_L1Tech_BSC_halo_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_halo_forPhysicsBackground", &HLT_L1Tech_BSC_halo_forPhysicsBackground, &b_HLT_L1Tech_BSC_halo_forPhysicsBackground);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_halo_forPhysicsBackground_Prescl", &HLT_L1Tech_BSC_halo_forPhysicsBackground_Prescl, &b_HLT_L1Tech_BSC_halo_forPhysicsBackground_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_HighMultiplicity", &HLT_L1Tech_BSC_HighMultiplicity, &b_HLT_L1Tech_BSC_HighMultiplicity);
   fChain->SetBranchAddress("HLT_L1Tech_BSC_HighMultiplicity_Prescl", &HLT_L1Tech_BSC_HighMultiplicity_Prescl, &b_HLT_L1Tech_BSC_HighMultiplicity_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_RPC_TTU_RBst1_collisions", &HLT_L1Tech_RPC_TTU_RBst1_collisions, &b_HLT_L1Tech_RPC_TTU_RBst1_collisions);
   fChain->SetBranchAddress("HLT_L1Tech_RPC_TTU_RBst1_collisions_Prescl", &HLT_L1Tech_RPC_TTU_RBst1_collisions_Prescl, &b_HLT_L1Tech_RPC_TTU_RBst1_collisions_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech_HCAL_HF", &HLT_L1Tech_HCAL_HF, &b_HLT_L1Tech_HCAL_HF);
   fChain->SetBranchAddress("HLT_L1Tech_HCAL_HF_Prescl", &HLT_L1Tech_HCAL_HF_Prescl, &b_HLT_L1Tech_HCAL_HF_Prescl);
   fChain->SetBranchAddress("HLT_TrackerCosmics", &HLT_TrackerCosmics, &b_HLT_TrackerCosmics);
   fChain->SetBranchAddress("HLT_TrackerCosmics_Prescl", &HLT_TrackerCosmics_Prescl, &b_HLT_TrackerCosmics_Prescl);
   fChain->SetBranchAddress("HLT_IsoTrackHB_v2", &HLT_IsoTrackHB_v2, &b_HLT_IsoTrackHB_v2);
   fChain->SetBranchAddress("HLT_IsoTrackHB_v2_Prescl", &HLT_IsoTrackHB_v2_Prescl, &b_HLT_IsoTrackHB_v2_Prescl);
   fChain->SetBranchAddress("HLT_IsoTrackHE_v2", &HLT_IsoTrackHE_v2, &b_HLT_IsoTrackHE_v2);
   fChain->SetBranchAddress("HLT_IsoTrackHE_v2_Prescl", &HLT_IsoTrackHE_v2_Prescl, &b_HLT_IsoTrackHE_v2_Prescl);
   fChain->SetBranchAddress("HLT_RPCBarrelCosmics", &HLT_RPCBarrelCosmics, &b_HLT_RPCBarrelCosmics);
   fChain->SetBranchAddress("HLT_RPCBarrelCosmics_Prescl", &HLT_RPCBarrelCosmics_Prescl, &b_HLT_RPCBarrelCosmics_Prescl);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_HcalPhiSym_Prescl", &HLT_HcalPhiSym_Prescl, &b_HLT_HcalPhiSym_Prescl);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalNZS_Prescl", &HLT_HcalNZS_Prescl, &b_HLT_HcalNZS_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity70", &HLT_PixelTracks_Multiplicity70, &b_HLT_PixelTracks_Multiplicity70);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity70_Prescl", &HLT_PixelTracks_Multiplicity70_Prescl, &b_HLT_PixelTracks_Multiplicity70_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85", &HLT_PixelTracks_Multiplicity85, &b_HLT_PixelTracks_Multiplicity85);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_Prescl", &HLT_PixelTracks_Multiplicity85_Prescl, &b_HLT_PixelTracks_Multiplicity85_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity100", &HLT_PixelTracks_Multiplicity100, &b_HLT_PixelTracks_Multiplicity100);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity100_Prescl", &HLT_PixelTracks_Multiplicity100_Prescl, &b_HLT_PixelTracks_Multiplicity100_Prescl);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise", &HLT_GlobalRunHPDNoise, &b_HLT_GlobalRunHPDNoise);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise_Prescl", &HLT_GlobalRunHPDNoise_Prescl, &b_HLT_GlobalRunHPDNoise_Prescl);
   fChain->SetBranchAddress("HLT_TechTrigHCALNoise", &HLT_TechTrigHCALNoise, &b_HLT_TechTrigHCALNoise);
   fChain->SetBranchAddress("HLT_TechTrigHCALNoise_Prescl", &HLT_TechTrigHCALNoise_Prescl, &b_HLT_TechTrigHCALNoise_Prescl);
   fChain->SetBranchAddress("HLT_L1_BPTX", &HLT_L1_BPTX, &b_HLT_L1_BPTX);
   fChain->SetBranchAddress("HLT_L1_BPTX_Prescl", &HLT_L1_BPTX_Prescl, &b_HLT_L1_BPTX_Prescl);
   fChain->SetBranchAddress("HLT_L1_BPTX_MinusOnly", &HLT_L1_BPTX_MinusOnly, &b_HLT_L1_BPTX_MinusOnly);
   fChain->SetBranchAddress("HLT_L1_BPTX_MinusOnly_Prescl", &HLT_L1_BPTX_MinusOnly_Prescl, &b_HLT_L1_BPTX_MinusOnly_Prescl);
   fChain->SetBranchAddress("HLT_L1_BPTX_PlusOnly", &HLT_L1_BPTX_PlusOnly, &b_HLT_L1_BPTX_PlusOnly);
   fChain->SetBranchAddress("HLT_L1_BPTX_PlusOnly_Prescl", &HLT_L1_BPTX_PlusOnly_Prescl, &b_HLT_L1_BPTX_PlusOnly_Prescl);
   fChain->SetBranchAddress("HLT_DTErrors", &HLT_DTErrors, &b_HLT_DTErrors);
   fChain->SetBranchAddress("HLT_DTErrors_Prescl", &HLT_DTErrors_Prescl, &b_HLT_DTErrors_Prescl);
   fChain->SetBranchAddress("HLT_LogMonitor", &HLT_LogMonitor, &b_HLT_LogMonitor);
   fChain->SetBranchAddress("HLT_LogMonitor_Prescl", &HLT_LogMonitor_Prescl, &b_HLT_LogMonitor_Prescl);
   fChain->SetBranchAddress("HLT_Calibration", &HLT_Calibration, &b_HLT_Calibration);
   fChain->SetBranchAddress("HLT_Calibration_Prescl", &HLT_Calibration_Prescl, &b_HLT_Calibration_Prescl);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_EcalCalibration_Prescl", &HLT_EcalCalibration_Prescl, &b_HLT_EcalCalibration_Prescl);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration_Prescl", &HLT_HcalCalibration_Prescl, &b_HLT_HcalCalibration_Prescl);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_Random_Prescl", &HLT_Random_Prescl, &b_HLT_Random_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPhiSym", &AlCa_EcalPhiSym, &b_AlCa_EcalPhiSym);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_Prescl", &AlCa_EcalPhiSym_Prescl, &b_AlCa_EcalPhiSym_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPi0", &AlCa_EcalPi0, &b_AlCa_EcalPi0);
   fChain->SetBranchAddress("AlCa_EcalPi0_Prescl", &AlCa_EcalPi0_Prescl, &b_AlCa_EcalPi0_Prescl);
   fChain->SetBranchAddress("AlCa_EcalEta", &AlCa_EcalEta, &b_AlCa_EcalEta);
   fChain->SetBranchAddress("AlCa_EcalEta_Prescl", &AlCa_EcalEta_Prescl, &b_AlCa_EcalEta_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits", &AlCa_RPCMuonNoHits, &b_AlCa_RPCMuonNoHits);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_Prescl", &AlCa_RPCMuonNoHits_Prescl, &b_AlCa_RPCMuonNoHits_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers", &AlCa_RPCMuonNoTriggers, &b_AlCa_RPCMuonNoTriggers);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_Prescl", &AlCa_RPCMuonNoTriggers_Prescl, &b_AlCa_RPCMuonNoTriggers_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation", &AlCa_RPCMuonNormalisation, &b_AlCa_RPCMuonNormalisation);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_Prescl", &AlCa_RPCMuonNormalisation_Prescl, &b_AlCa_RPCMuonNormalisation_Prescl);
   fChain->SetBranchAddress("DQM_FEDIntegrity", &DQM_FEDIntegrity, &b_DQM_FEDIntegrity);
   fChain->SetBranchAddress("DQM_FEDIntegrity_Prescl", &DQM_FEDIntegrity_Prescl, &b_DQM_FEDIntegrity_Prescl);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFinalPath_Prescl", &HLTriggerFinalPath_Prescl, &b_HLTriggerFinalPath_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   fChain->SetBranchAddress("L1_BptxMinus_Prescl", &L1_BptxMinus_Prescl, &b_L1_BptxMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus_5bx", &L1_BptxMinus_5bx, &b_L1_BptxMinus_5bx);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus", &L1_BptxMinus_NotBptxPlus, &b_L1_BptxMinus_NotBptxPlus);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus_Prescl", &L1_BptxMinus_NotBptxPlus_Prescl, &b_L1_BptxMinus_NotBptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus_5bx", &L1_BptxMinus_NotBptxPlus_5bx, &b_L1_BptxMinus_NotBptxPlus_5bx);
   fChain->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   fChain->SetBranchAddress("L1_BptxPlus_Prescl", &L1_BptxPlus_Prescl, &b_L1_BptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlus_5bx", &L1_BptxPlus_5bx, &b_L1_BptxPlus_5bx);
   fChain->SetBranchAddress("L1_BptxPlusORMinus", &L1_BptxPlusORMinus, &b_L1_BptxPlusORMinus);
   fChain->SetBranchAddress("L1_BptxPlusORMinus_Prescl", &L1_BptxPlusORMinus_Prescl, &b_L1_BptxPlusORMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlusORMinus_5bx", &L1_BptxPlusORMinus_5bx, &b_L1_BptxPlusORMinus_5bx);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus", &L1_BptxPlus_NotBptxMinus, &b_L1_BptxPlus_NotBptxMinus);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus_Prescl", &L1_BptxPlus_NotBptxMinus_Prescl, &b_L1_BptxPlus_NotBptxMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus_5bx", &L1_BptxPlus_NotBptxMinus_5bx, &b_L1_BptxPlus_NotBptxMinus_5bx);
   fChain->SetBranchAddress("L1_BptxXOR_BscMinBiasOR", &L1_BptxXOR_BscMinBiasOR, &b_L1_BptxXOR_BscMinBiasOR);
   fChain->SetBranchAddress("L1_BptxXOR_BscMinBiasOR_Prescl", &L1_BptxXOR_BscMinBiasOR_Prescl, &b_L1_BptxXOR_BscMinBiasOR_Prescl);
   fChain->SetBranchAddress("L1_BptxXOR_BscMinBiasOR_5bx", &L1_BptxXOR_BscMinBiasOR_5bx, &b_L1_BptxXOR_BscMinBiasOR_5bx);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus", &L1_Bsc2Minus_BptxMinus, &b_L1_Bsc2Minus_BptxMinus);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus_Prescl", &L1_Bsc2Minus_BptxMinus_Prescl, &b_L1_Bsc2Minus_BptxMinus_Prescl);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus_5bx", &L1_Bsc2Minus_BptxMinus_5bx, &b_L1_Bsc2Minus_BptxMinus_5bx);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus", &L1_Bsc2Plus_BptxPlus, &b_L1_Bsc2Plus_BptxPlus);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus_Prescl", &L1_Bsc2Plus_BptxPlus_Prescl, &b_L1_Bsc2Plus_BptxPlus_Prescl);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus_5bx", &L1_Bsc2Plus_BptxPlus_5bx, &b_L1_Bsc2Plus_BptxPlus_5bx);
   fChain->SetBranchAddress("L1_BscHaloBeam1Inner", &L1_BscHaloBeam1Inner, &b_L1_BscHaloBeam1Inner);
   fChain->SetBranchAddress("L1_BscHaloBeam1Inner_Prescl", &L1_BscHaloBeam1Inner_Prescl, &b_L1_BscHaloBeam1Inner_Prescl);
   fChain->SetBranchAddress("L1_BscHaloBeam1Inner_5bx", &L1_BscHaloBeam1Inner_5bx, &b_L1_BscHaloBeam1Inner_5bx);
   fChain->SetBranchAddress("L1_BscHaloBeam1Outer", &L1_BscHaloBeam1Outer, &b_L1_BscHaloBeam1Outer);
   fChain->SetBranchAddress("L1_BscHaloBeam1Outer_Prescl", &L1_BscHaloBeam1Outer_Prescl, &b_L1_BscHaloBeam1Outer_Prescl);
   fChain->SetBranchAddress("L1_BscHaloBeam1Outer_5bx", &L1_BscHaloBeam1Outer_5bx, &b_L1_BscHaloBeam1Outer_5bx);
   fChain->SetBranchAddress("L1_BscHaloBeam2Inner", &L1_BscHaloBeam2Inner, &b_L1_BscHaloBeam2Inner);
   fChain->SetBranchAddress("L1_BscHaloBeam2Inner_Prescl", &L1_BscHaloBeam2Inner_Prescl, &b_L1_BscHaloBeam2Inner_Prescl);
   fChain->SetBranchAddress("L1_BscHaloBeam2Inner_5bx", &L1_BscHaloBeam2Inner_5bx, &b_L1_BscHaloBeam2Inner_5bx);
   fChain->SetBranchAddress("L1_BscHaloBeam2Outer", &L1_BscHaloBeam2Outer, &b_L1_BscHaloBeam2Outer);
   fChain->SetBranchAddress("L1_BscHaloBeam2Outer_Prescl", &L1_BscHaloBeam2Outer_Prescl, &b_L1_BscHaloBeam2Outer_Prescl);
   fChain->SetBranchAddress("L1_BscHaloBeam2Outer_5bx", &L1_BscHaloBeam2Outer_5bx, &b_L1_BscHaloBeam2Outer_5bx);
   fChain->SetBranchAddress("L1_BscHighMultiplicity", &L1_BscHighMultiplicity, &b_L1_BscHighMultiplicity);
   fChain->SetBranchAddress("L1_BscHighMultiplicity_Prescl", &L1_BscHighMultiplicity_Prescl, &b_L1_BscHighMultiplicity_Prescl);
   fChain->SetBranchAddress("L1_BscHighMultiplicity_5bx", &L1_BscHighMultiplicity_5bx, &b_L1_BscHighMultiplicity_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1", &L1_BscMinBiasInnerThreshold1, &b_L1_BscMinBiasInnerThreshold1);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1_Prescl", &L1_BscMinBiasInnerThreshold1_Prescl, &b_L1_BscMinBiasInnerThreshold1_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1_5bx", &L1_BscMinBiasInnerThreshold1_5bx, &b_L1_BscMinBiasInnerThreshold1_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold2", &L1_BscMinBiasInnerThreshold2, &b_L1_BscMinBiasInnerThreshold2);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold2_Prescl", &L1_BscMinBiasInnerThreshold2_Prescl, &b_L1_BscMinBiasInnerThreshold2_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold2_5bx", &L1_BscMinBiasInnerThreshold2_5bx, &b_L1_BscMinBiasInnerThreshold2_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasOR", &L1_BscMinBiasOR, &b_L1_BscMinBiasOR);
   fChain->SetBranchAddress("L1_BscMinBiasOR_Prescl", &L1_BscMinBiasOR_Prescl, &b_L1_BscMinBiasOR_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_5bx", &L1_BscMinBiasOR_5bx, &b_L1_BscMinBiasOR_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusANDMinus", &L1_BscMinBiasOR_BptxPlusANDMinus, &b_L1_BscMinBiasOR_BptxPlusANDMinus);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusANDMinus_Prescl", &L1_BscMinBiasOR_BptxPlusANDMinus_Prescl, &b_L1_BscMinBiasOR_BptxPlusANDMinus_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusANDMinus_5bx", &L1_BscMinBiasOR_BptxPlusANDMinus_5bx, &b_L1_BscMinBiasOR_BptxPlusANDMinus_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus", &L1_BscMinBiasOR_BptxPlusORMinus, &b_L1_BscMinBiasOR_BptxPlusORMinus);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus_Prescl", &L1_BscMinBiasOR_BptxPlusORMinus_Prescl, &b_L1_BscMinBiasOR_BptxPlusORMinus_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus_5bx", &L1_BscMinBiasOR_BptxPlusORMinus_5bx, &b_L1_BscMinBiasOR_BptxPlusORMinus_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1", &L1_BscMinBiasThreshold1, &b_L1_BscMinBiasThreshold1);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_Prescl", &L1_BscMinBiasThreshold1_Prescl, &b_L1_BscMinBiasThreshold1_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_5bx", &L1_BscMinBiasThreshold1_5bx, &b_L1_BscMinBiasThreshold1_5bx);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2", &L1_BscMinBiasThreshold2, &b_L1_BscMinBiasThreshold2);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2_Prescl", &L1_BscMinBiasThreshold2_Prescl, &b_L1_BscMinBiasThreshold2_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2_5bx", &L1_BscMinBiasThreshold2_5bx, &b_L1_BscMinBiasThreshold2_5bx);
   fChain->SetBranchAddress("L1_BscSplashBeam1", &L1_BscSplashBeam1, &b_L1_BscSplashBeam1);
   fChain->SetBranchAddress("L1_BscSplashBeam1_Prescl", &L1_BscSplashBeam1_Prescl, &b_L1_BscSplashBeam1_Prescl);
   fChain->SetBranchAddress("L1_BscSplashBeam1_5bx", &L1_BscSplashBeam1_5bx, &b_L1_BscSplashBeam1_5bx);
   fChain->SetBranchAddress("L1_BscSplashBeam2", &L1_BscSplashBeam2, &b_L1_BscSplashBeam2);
   fChain->SetBranchAddress("L1_BscSplashBeam2_Prescl", &L1_BscSplashBeam2_Prescl, &b_L1_BscSplashBeam2_Prescl);
   fChain->SetBranchAddress("L1_BscSplashBeam2_5bx", &L1_BscSplashBeam2_5bx, &b_L1_BscSplashBeam2_5bx);
   fChain->SetBranchAddress("L1_DoubleEG05_TopBottom", &L1_DoubleEG05_TopBottom, &b_L1_DoubleEG05_TopBottom);
   fChain->SetBranchAddress("L1_DoubleEG05_TopBottom_Prescl", &L1_DoubleEG05_TopBottom_Prescl, &b_L1_DoubleEG05_TopBottom_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG05_TopBottom_5bx", &L1_DoubleEG05_TopBottom_5bx, &b_L1_DoubleEG05_TopBottom_5bx);
   fChain->SetBranchAddress("L1_DoubleEG2", &L1_DoubleEG2, &b_L1_DoubleEG2);
   fChain->SetBranchAddress("L1_DoubleEG2_Prescl", &L1_DoubleEG2_Prescl, &b_L1_DoubleEG2_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG2_5bx", &L1_DoubleEG2_5bx, &b_L1_DoubleEG2_5bx);
   fChain->SetBranchAddress("L1_DoubleEG5", &L1_DoubleEG5, &b_L1_DoubleEG5);
   fChain->SetBranchAddress("L1_DoubleEG5_Prescl", &L1_DoubleEG5_Prescl, &b_L1_DoubleEG5_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5_5bx", &L1_DoubleEG5_5bx, &b_L1_DoubleEG5_5bx);
   fChain->SetBranchAddress("L1_DoubleForJet10_EtaOpp", &L1_DoubleForJet10_EtaOpp, &b_L1_DoubleForJet10_EtaOpp);
   fChain->SetBranchAddress("L1_DoubleForJet10_EtaOpp_Prescl", &L1_DoubleForJet10_EtaOpp_Prescl, &b_L1_DoubleForJet10_EtaOpp_Prescl);
   fChain->SetBranchAddress("L1_DoubleForJet10_EtaOpp_5bx", &L1_DoubleForJet10_EtaOpp_5bx, &b_L1_DoubleForJet10_EtaOpp_5bx);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1", &L1_DoubleHfBitCountsRing1_P1N1, &b_L1_DoubleHfBitCountsRing1_P1N1);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1_Prescl", &L1_DoubleHfBitCountsRing1_P1N1_Prescl, &b_L1_DoubleHfBitCountsRing1_P1N1_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1_5bx", &L1_DoubleHfBitCountsRing1_P1N1_5bx, &b_L1_DoubleHfBitCountsRing1_P1N1_5bx);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1", &L1_DoubleHfBitCountsRing2_P1N1, &b_L1_DoubleHfBitCountsRing2_P1N1);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1_Prescl", &L1_DoubleHfBitCountsRing2_P1N1_Prescl, &b_L1_DoubleHfBitCountsRing2_P1N1_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1_5bx", &L1_DoubleHfBitCountsRing2_P1N1_5bx, &b_L1_DoubleHfBitCountsRing2_P1N1_5bx);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200", &L1_DoubleHfRingEtSumsRing1_P200N200, &b_L1_DoubleHfRingEtSumsRing1_P200N200);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200_Prescl", &L1_DoubleHfRingEtSumsRing1_P200N200_Prescl, &b_L1_DoubleHfRingEtSumsRing1_P200N200_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200_5bx", &L1_DoubleHfRingEtSumsRing1_P200N200_5bx, &b_L1_DoubleHfRingEtSumsRing1_P200N200_5bx);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4", &L1_DoubleHfRingEtSumsRing1_P4N4, &b_L1_DoubleHfRingEtSumsRing1_P4N4);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4_Prescl", &L1_DoubleHfRingEtSumsRing1_P4N4_Prescl, &b_L1_DoubleHfRingEtSumsRing1_P4N4_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4_5bx", &L1_DoubleHfRingEtSumsRing1_P4N4_5bx, &b_L1_DoubleHfRingEtSumsRing1_P4N4_5bx);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200", &L1_DoubleHfRingEtSumsRing2_P200N200, &b_L1_DoubleHfRingEtSumsRing2_P200N200);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200_Prescl", &L1_DoubleHfRingEtSumsRing2_P200N200_Prescl, &b_L1_DoubleHfRingEtSumsRing2_P200N200_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200_5bx", &L1_DoubleHfRingEtSumsRing2_P200N200_5bx, &b_L1_DoubleHfRingEtSumsRing2_P200N200_5bx);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4", &L1_DoubleHfRingEtSumsRing2_P4N4, &b_L1_DoubleHfRingEtSumsRing2_P4N4);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4_Prescl", &L1_DoubleHfRingEtSumsRing2_P4N4_Prescl, &b_L1_DoubleHfRingEtSumsRing2_P4N4_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4_5bx", &L1_DoubleHfRingEtSumsRing2_P4N4_5bx, &b_L1_DoubleHfRingEtSumsRing2_P4N4_5bx);
   fChain->SetBranchAddress("L1_DoubleJet30", &L1_DoubleJet30, &b_L1_DoubleJet30);
   fChain->SetBranchAddress("L1_DoubleJet30_Prescl", &L1_DoubleJet30_Prescl, &b_L1_DoubleJet30_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet30_5bx", &L1_DoubleJet30_5bx, &b_L1_DoubleJet30_5bx);
   fChain->SetBranchAddress("L1_DoubleMu3", &L1_DoubleMu3, &b_L1_DoubleMu3);
   fChain->SetBranchAddress("L1_DoubleMu3_Prescl", &L1_DoubleMu3_Prescl, &b_L1_DoubleMu3_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3_5bx", &L1_DoubleMu3_5bx, &b_L1_DoubleMu3_5bx);
   fChain->SetBranchAddress("L1_DoubleMuOpen", &L1_DoubleMuOpen, &b_L1_DoubleMuOpen);
   fChain->SetBranchAddress("L1_DoubleMuOpen_Prescl", &L1_DoubleMuOpen_Prescl, &b_L1_DoubleMuOpen_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_5bx", &L1_DoubleMuOpen_5bx, &b_L1_DoubleMuOpen_5bx);
   fChain->SetBranchAddress("L1_DoubleMuTopBottom", &L1_DoubleMuTopBottom, &b_L1_DoubleMuTopBottom);
   fChain->SetBranchAddress("L1_DoubleMuTopBottom_Prescl", &L1_DoubleMuTopBottom_Prescl, &b_L1_DoubleMuTopBottom_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuTopBottom_5bx", &L1_DoubleMuTopBottom_5bx, &b_L1_DoubleMuTopBottom_5bx);
   fChain->SetBranchAddress("L1_DoubleTauJet14", &L1_DoubleTauJet14, &b_L1_DoubleTauJet14);
   fChain->SetBranchAddress("L1_DoubleTauJet14_Prescl", &L1_DoubleTauJet14_Prescl, &b_L1_DoubleTauJet14_Prescl);
   fChain->SetBranchAddress("L1_DoubleTauJet14_5bx", &L1_DoubleTauJet14_5bx, &b_L1_DoubleTauJet14_5bx);
   fChain->SetBranchAddress("L1_ETM12", &L1_ETM12, &b_L1_ETM12);
   fChain->SetBranchAddress("L1_ETM12_Prescl", &L1_ETM12_Prescl, &b_L1_ETM12_Prescl);
   fChain->SetBranchAddress("L1_ETM12_5bx", &L1_ETM12_5bx, &b_L1_ETM12_5bx);
   fChain->SetBranchAddress("L1_ETM20", &L1_ETM20, &b_L1_ETM20);
   fChain->SetBranchAddress("L1_ETM20_Prescl", &L1_ETM20_Prescl, &b_L1_ETM20_Prescl);
   fChain->SetBranchAddress("L1_ETM20_5bx", &L1_ETM20_5bx, &b_L1_ETM20_5bx);
   fChain->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   fChain->SetBranchAddress("L1_ETM30_Prescl", &L1_ETM30_Prescl, &b_L1_ETM30_Prescl);
   fChain->SetBranchAddress("L1_ETM30_5bx", &L1_ETM30_5bx, &b_L1_ETM30_5bx);
   fChain->SetBranchAddress("L1_ETM70", &L1_ETM70, &b_L1_ETM70);
   fChain->SetBranchAddress("L1_ETM70_Prescl", &L1_ETM70_Prescl, &b_L1_ETM70_Prescl);
   fChain->SetBranchAddress("L1_ETM70_5bx", &L1_ETM70_5bx, &b_L1_ETM70_5bx);
   fChain->SetBranchAddress("L1_ETT100", &L1_ETT100, &b_L1_ETT100);
   fChain->SetBranchAddress("L1_ETT100_Prescl", &L1_ETT100_Prescl, &b_L1_ETT100_Prescl);
   fChain->SetBranchAddress("L1_ETT100_5bx", &L1_ETT100_5bx, &b_L1_ETT100_5bx);
   fChain->SetBranchAddress("L1_ETT140", &L1_ETT140, &b_L1_ETT140);
   fChain->SetBranchAddress("L1_ETT140_Prescl", &L1_ETT140_Prescl, &b_L1_ETT140_Prescl);
   fChain->SetBranchAddress("L1_ETT140_5bx", &L1_ETT140_5bx, &b_L1_ETT140_5bx);
   fChain->SetBranchAddress("L1_ETT30", &L1_ETT30, &b_L1_ETT30);
   fChain->SetBranchAddress("L1_ETT30_Prescl", &L1_ETT30_Prescl, &b_L1_ETT30_Prescl);
   fChain->SetBranchAddress("L1_ETT30_5bx", &L1_ETT30_5bx, &b_L1_ETT30_5bx);
   fChain->SetBranchAddress("L1_ETT60", &L1_ETT60, &b_L1_ETT60);
   fChain->SetBranchAddress("L1_ETT60_Prescl", &L1_ETT60_Prescl, &b_L1_ETT60_Prescl);
   fChain->SetBranchAddress("L1_ETT60_5bx", &L1_ETT60_5bx, &b_L1_ETT60_5bx);
   fChain->SetBranchAddress("L1_HTM20", &L1_HTM20, &b_L1_HTM20);
   fChain->SetBranchAddress("L1_HTM20_Prescl", &L1_HTM20_Prescl, &b_L1_HTM20_Prescl);
   fChain->SetBranchAddress("L1_HTM20_5bx", &L1_HTM20_5bx, &b_L1_HTM20_5bx);
   fChain->SetBranchAddress("L1_HTM30", &L1_HTM30, &b_L1_HTM30);
   fChain->SetBranchAddress("L1_HTM30_Prescl", &L1_HTM30_Prescl, &b_L1_HTM30_Prescl);
   fChain->SetBranchAddress("L1_HTM30_5bx", &L1_HTM30_5bx, &b_L1_HTM30_5bx);
   fChain->SetBranchAddress("L1_HTT100", &L1_HTT100, &b_L1_HTT100);
   fChain->SetBranchAddress("L1_HTT100_Prescl", &L1_HTT100_Prescl, &b_L1_HTT100_Prescl);
   fChain->SetBranchAddress("L1_HTT100_5bx", &L1_HTT100_5bx, &b_L1_HTT100_5bx);
   fChain->SetBranchAddress("L1_HTT200", &L1_HTT200, &b_L1_HTT200);
   fChain->SetBranchAddress("L1_HTT200_Prescl", &L1_HTT200_Prescl, &b_L1_HTT200_Prescl);
   fChain->SetBranchAddress("L1_HTT200_5bx", &L1_HTT200_5bx, &b_L1_HTT200_5bx);
   fChain->SetBranchAddress("L1_HTT50", &L1_HTT50, &b_L1_HTT50);
   fChain->SetBranchAddress("L1_HTT50_Prescl", &L1_HTT50_Prescl, &b_L1_HTT50_Prescl);
   fChain->SetBranchAddress("L1_HTT50_5bx", &L1_HTT50_5bx, &b_L1_HTT50_5bx);
   fChain->SetBranchAddress("L1_IsoEG10_Jet6_ForJet6", &L1_IsoEG10_Jet6_ForJet6, &b_L1_IsoEG10_Jet6_ForJet6);
   fChain->SetBranchAddress("L1_IsoEG10_Jet6_ForJet6_Prescl", &L1_IsoEG10_Jet6_ForJet6_Prescl, &b_L1_IsoEG10_Jet6_ForJet6_Prescl);
   fChain->SetBranchAddress("L1_IsoEG10_Jet6_ForJet6_5bx", &L1_IsoEG10_Jet6_ForJet6_5bx, &b_L1_IsoEG10_Jet6_ForJet6_5bx);
   fChain->SetBranchAddress("L1_Mu3_EG5", &L1_Mu3_EG5, &b_L1_Mu3_EG5);
   fChain->SetBranchAddress("L1_Mu3_EG5_Prescl", &L1_Mu3_EG5_Prescl, &b_L1_Mu3_EG5_Prescl);
   fChain->SetBranchAddress("L1_Mu3_EG5_5bx", &L1_Mu3_EG5_5bx, &b_L1_Mu3_EG5_5bx);
   fChain->SetBranchAddress("L1_Mu3_Jet10", &L1_Mu3_Jet10, &b_L1_Mu3_Jet10);
   fChain->SetBranchAddress("L1_Mu3_Jet10_Prescl", &L1_Mu3_Jet10_Prescl, &b_L1_Mu3_Jet10_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet10_5bx", &L1_Mu3_Jet10_5bx, &b_L1_Mu3_Jet10_5bx);
   fChain->SetBranchAddress("L1_Mu3_Jet6", &L1_Mu3_Jet6, &b_L1_Mu3_Jet6);
   fChain->SetBranchAddress("L1_Mu3_Jet6_Prescl", &L1_Mu3_Jet6_Prescl, &b_L1_Mu3_Jet6_Prescl);
   fChain->SetBranchAddress("L1_Mu3_Jet6_5bx", &L1_Mu3_Jet6_5bx, &b_L1_Mu3_Jet6_5bx);
   fChain->SetBranchAddress("L1_Mu5_Jet6", &L1_Mu5_Jet6, &b_L1_Mu5_Jet6);
   fChain->SetBranchAddress("L1_Mu5_Jet6_Prescl", &L1_Mu5_Jet6_Prescl, &b_L1_Mu5_Jet6_Prescl);
   fChain->SetBranchAddress("L1_Mu5_Jet6_5bx", &L1_Mu5_Jet6_5bx, &b_L1_Mu5_Jet6_5bx);
   fChain->SetBranchAddress("L1_QuadJet6", &L1_QuadJet6, &b_L1_QuadJet6);
   fChain->SetBranchAddress("L1_QuadJet6_Prescl", &L1_QuadJet6_Prescl, &b_L1_QuadJet6_Prescl);
   fChain->SetBranchAddress("L1_QuadJet6_5bx", &L1_QuadJet6_5bx, &b_L1_QuadJet6_5bx);
   fChain->SetBranchAddress("L1_QuadJet8", &L1_QuadJet8, &b_L1_QuadJet8);
   fChain->SetBranchAddress("L1_QuadJet8_Prescl", &L1_QuadJet8_Prescl, &b_L1_QuadJet8_Prescl);
   fChain->SetBranchAddress("L1_QuadJet8_5bx", &L1_QuadJet8_5bx, &b_L1_QuadJet8_5bx);
   fChain->SetBranchAddress("L1_SingleCenJet2", &L1_SingleCenJet2, &b_L1_SingleCenJet2);
   fChain->SetBranchAddress("L1_SingleCenJet2_Prescl", &L1_SingleCenJet2_Prescl, &b_L1_SingleCenJet2_Prescl);
   fChain->SetBranchAddress("L1_SingleCenJet2_5bx", &L1_SingleCenJet2_5bx, &b_L1_SingleCenJet2_5bx);
   fChain->SetBranchAddress("L1_SingleCenJet4", &L1_SingleCenJet4, &b_L1_SingleCenJet4);
   fChain->SetBranchAddress("L1_SingleCenJet4_Prescl", &L1_SingleCenJet4_Prescl, &b_L1_SingleCenJet4_Prescl);
   fChain->SetBranchAddress("L1_SingleCenJet4_5bx", &L1_SingleCenJet4_5bx, &b_L1_SingleCenJet4_5bx);
   fChain->SetBranchAddress("L1_SingleEG1", &L1_SingleEG1, &b_L1_SingleEG1);
   fChain->SetBranchAddress("L1_SingleEG1_Prescl", &L1_SingleEG1_Prescl, &b_L1_SingleEG1_Prescl);
   fChain->SetBranchAddress("L1_SingleEG1_5bx", &L1_SingleEG1_5bx, &b_L1_SingleEG1_5bx);
   fChain->SetBranchAddress("L1_SingleEG10", &L1_SingleEG10, &b_L1_SingleEG10);
   fChain->SetBranchAddress("L1_SingleEG10_Prescl", &L1_SingleEG10_Prescl, &b_L1_SingleEG10_Prescl);
   fChain->SetBranchAddress("L1_SingleEG10_5bx", &L1_SingleEG10_5bx, &b_L1_SingleEG10_5bx);
   fChain->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12, &b_L1_SingleEG12);
   fChain->SetBranchAddress("L1_SingleEG12_Prescl", &L1_SingleEG12_Prescl, &b_L1_SingleEG12_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12_5bx", &L1_SingleEG12_5bx, &b_L1_SingleEG12_5bx);
   fChain->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   fChain->SetBranchAddress("L1_SingleEG15_Prescl", &L1_SingleEG15_Prescl, &b_L1_SingleEG15_Prescl);
   fChain->SetBranchAddress("L1_SingleEG15_5bx", &L1_SingleEG15_5bx, &b_L1_SingleEG15_5bx);
   fChain->SetBranchAddress("L1_SingleEG2", &L1_SingleEG2, &b_L1_SingleEG2);
   fChain->SetBranchAddress("L1_SingleEG2_Prescl", &L1_SingleEG2_Prescl, &b_L1_SingleEG2_Prescl);
   fChain->SetBranchAddress("L1_SingleEG2_5bx", &L1_SingleEG2_5bx, &b_L1_SingleEG2_5bx);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG20_Prescl", &L1_SingleEG20_Prescl, &b_L1_SingleEG20_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20_5bx", &L1_SingleEG20_5bx, &b_L1_SingleEG20_5bx);
   fChain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   fChain->SetBranchAddress("L1_SingleEG5_Prescl", &L1_SingleEG5_Prescl, &b_L1_SingleEG5_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5_5bx", &L1_SingleEG5_5bx, &b_L1_SingleEG5_5bx);
   fChain->SetBranchAddress("L1_SingleEG8", &L1_SingleEG8, &b_L1_SingleEG8);
   fChain->SetBranchAddress("L1_SingleEG8_Prescl", &L1_SingleEG8_Prescl, &b_L1_SingleEG8_Prescl);
   fChain->SetBranchAddress("L1_SingleEG8_5bx", &L1_SingleEG8_5bx, &b_L1_SingleEG8_5bx);
   fChain->SetBranchAddress("L1_SingleForJet2", &L1_SingleForJet2, &b_L1_SingleForJet2);
   fChain->SetBranchAddress("L1_SingleForJet2_Prescl", &L1_SingleForJet2_Prescl, &b_L1_SingleForJet2_Prescl);
   fChain->SetBranchAddress("L1_SingleForJet2_5bx", &L1_SingleForJet2_5bx, &b_L1_SingleForJet2_5bx);
   fChain->SetBranchAddress("L1_SingleForJet4", &L1_SingleForJet4, &b_L1_SingleForJet4);
   fChain->SetBranchAddress("L1_SingleForJet4_Prescl", &L1_SingleForJet4_Prescl, &b_L1_SingleForJet4_Prescl);
   fChain->SetBranchAddress("L1_SingleForJet4_5bx", &L1_SingleForJet4_5bx, &b_L1_SingleForJet4_5bx);
   fChain->SetBranchAddress("L1_SingleHfBitCountsRing1_1", &L1_SingleHfBitCountsRing1_1, &b_L1_SingleHfBitCountsRing1_1);
   fChain->SetBranchAddress("L1_SingleHfBitCountsRing1_1_Prescl", &L1_SingleHfBitCountsRing1_1_Prescl, &b_L1_SingleHfBitCountsRing1_1_Prescl);
   fChain->SetBranchAddress("L1_SingleHfBitCountsRing1_1_5bx", &L1_SingleHfBitCountsRing1_1_5bx, &b_L1_SingleHfBitCountsRing1_1_5bx);
   fChain->SetBranchAddress("L1_SingleHfBitCountsRing2_1", &L1_SingleHfBitCountsRing2_1, &b_L1_SingleHfBitCountsRing2_1);
   fChain->SetBranchAddress("L1_SingleHfBitCountsRing2_1_Prescl", &L1_SingleHfBitCountsRing2_1_Prescl, &b_L1_SingleHfBitCountsRing2_1_Prescl);
   fChain->SetBranchAddress("L1_SingleHfBitCountsRing2_1_5bx", &L1_SingleHfBitCountsRing2_1_5bx, &b_L1_SingleHfBitCountsRing2_1_5bx);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing1_200", &L1_SingleHfRingEtSumsRing1_200, &b_L1_SingleHfRingEtSumsRing1_200);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing1_200_Prescl", &L1_SingleHfRingEtSumsRing1_200_Prescl, &b_L1_SingleHfRingEtSumsRing1_200_Prescl);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing1_200_5bx", &L1_SingleHfRingEtSumsRing1_200_5bx, &b_L1_SingleHfRingEtSumsRing1_200_5bx);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing1_4", &L1_SingleHfRingEtSumsRing1_4, &b_L1_SingleHfRingEtSumsRing1_4);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing1_4_Prescl", &L1_SingleHfRingEtSumsRing1_4_Prescl, &b_L1_SingleHfRingEtSumsRing1_4_Prescl);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing1_4_5bx", &L1_SingleHfRingEtSumsRing1_4_5bx, &b_L1_SingleHfRingEtSumsRing1_4_5bx);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing2_200", &L1_SingleHfRingEtSumsRing2_200, &b_L1_SingleHfRingEtSumsRing2_200);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing2_200_Prescl", &L1_SingleHfRingEtSumsRing2_200_Prescl, &b_L1_SingleHfRingEtSumsRing2_200_Prescl);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing2_200_5bx", &L1_SingleHfRingEtSumsRing2_200_5bx, &b_L1_SingleHfRingEtSumsRing2_200_5bx);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing2_4", &L1_SingleHfRingEtSumsRing2_4, &b_L1_SingleHfRingEtSumsRing2_4);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing2_4_Prescl", &L1_SingleHfRingEtSumsRing2_4_Prescl, &b_L1_SingleHfRingEtSumsRing2_4_Prescl);
   fChain->SetBranchAddress("L1_SingleHfRingEtSumsRing2_4_5bx", &L1_SingleHfRingEtSumsRing2_4_5bx, &b_L1_SingleHfRingEtSumsRing2_4_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG10", &L1_SingleIsoEG10, &b_L1_SingleIsoEG10);
   fChain->SetBranchAddress("L1_SingleIsoEG10_Prescl", &L1_SingleIsoEG10_Prescl, &b_L1_SingleIsoEG10_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG10_5bx", &L1_SingleIsoEG10_5bx, &b_L1_SingleIsoEG10_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG12", &L1_SingleIsoEG12, &b_L1_SingleIsoEG12);
   fChain->SetBranchAddress("L1_SingleIsoEG12_Prescl", &L1_SingleIsoEG12_Prescl, &b_L1_SingleIsoEG12_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG12_5bx", &L1_SingleIsoEG12_5bx, &b_L1_SingleIsoEG12_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG15", &L1_SingleIsoEG15, &b_L1_SingleIsoEG15);
   fChain->SetBranchAddress("L1_SingleIsoEG15_Prescl", &L1_SingleIsoEG15_Prescl, &b_L1_SingleIsoEG15_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG15_5bx", &L1_SingleIsoEG15_5bx, &b_L1_SingleIsoEG15_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG5", &L1_SingleIsoEG5, &b_L1_SingleIsoEG5);
   fChain->SetBranchAddress("L1_SingleIsoEG5_Prescl", &L1_SingleIsoEG5_Prescl, &b_L1_SingleIsoEG5_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG5_5bx", &L1_SingleIsoEG5_5bx, &b_L1_SingleIsoEG5_5bx);
   fChain->SetBranchAddress("L1_SingleIsoEG8", &L1_SingleIsoEG8, &b_L1_SingleIsoEG8);
   fChain->SetBranchAddress("L1_SingleIsoEG8_Prescl", &L1_SingleIsoEG8_Prescl, &b_L1_SingleIsoEG8_Prescl);
   fChain->SetBranchAddress("L1_SingleIsoEG8_5bx", &L1_SingleIsoEG8_5bx, &b_L1_SingleIsoEG8_5bx);
   fChain->SetBranchAddress("L1_SingleJet10", &L1_SingleJet10, &b_L1_SingleJet10);
   fChain->SetBranchAddress("L1_SingleJet10_Prescl", &L1_SingleJet10_Prescl, &b_L1_SingleJet10_Prescl);
   fChain->SetBranchAddress("L1_SingleJet10_5bx", &L1_SingleJet10_5bx, &b_L1_SingleJet10_5bx);
   fChain->SetBranchAddress("L1_SingleJet10_NotBptxOR_Ext", &L1_SingleJet10_NotBptxOR_Ext, &b_L1_SingleJet10_NotBptxOR_Ext);
   fChain->SetBranchAddress("L1_SingleJet10_NotBptxOR_Ext_Prescl", &L1_SingleJet10_NotBptxOR_Ext_Prescl, &b_L1_SingleJet10_NotBptxOR_Ext_Prescl);
   fChain->SetBranchAddress("L1_SingleJet10_NotBptxOR_Ext_5bx", &L1_SingleJet10_NotBptxOR_Ext_5bx, &b_L1_SingleJet10_NotBptxOR_Ext_5bx);
   fChain->SetBranchAddress("L1_SingleJet20", &L1_SingleJet20, &b_L1_SingleJet20);
   fChain->SetBranchAddress("L1_SingleJet20_Prescl", &L1_SingleJet20_Prescl, &b_L1_SingleJet20_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_5bx", &L1_SingleJet20_5bx, &b_L1_SingleJet20_5bx);
   fChain->SetBranchAddress("L1_SingleJet30", &L1_SingleJet30, &b_L1_SingleJet30);
   fChain->SetBranchAddress("L1_SingleJet30_Prescl", &L1_SingleJet30_Prescl, &b_L1_SingleJet30_Prescl);
   fChain->SetBranchAddress("L1_SingleJet30_5bx", &L1_SingleJet30_5bx, &b_L1_SingleJet30_5bx);
   fChain->SetBranchAddress("L1_SingleJet40", &L1_SingleJet40, &b_L1_SingleJet40);
   fChain->SetBranchAddress("L1_SingleJet40_Prescl", &L1_SingleJet40_Prescl, &b_L1_SingleJet40_Prescl);
   fChain->SetBranchAddress("L1_SingleJet40_5bx", &L1_SingleJet40_5bx, &b_L1_SingleJet40_5bx);
   fChain->SetBranchAddress("L1_SingleJet50", &L1_SingleJet50, &b_L1_SingleJet50);
   fChain->SetBranchAddress("L1_SingleJet50_Prescl", &L1_SingleJet50_Prescl, &b_L1_SingleJet50_Prescl);
   fChain->SetBranchAddress("L1_SingleJet50_5bx", &L1_SingleJet50_5bx, &b_L1_SingleJet50_5bx);
   fChain->SetBranchAddress("L1_SingleJet6", &L1_SingleJet6, &b_L1_SingleJet6);
   fChain->SetBranchAddress("L1_SingleJet6_Prescl", &L1_SingleJet6_Prescl, &b_L1_SingleJet6_Prescl);
   fChain->SetBranchAddress("L1_SingleJet6_5bx", &L1_SingleJet6_5bx, &b_L1_SingleJet6_5bx);
   fChain->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60, &b_L1_SingleJet60);
   fChain->SetBranchAddress("L1_SingleJet60_Prescl", &L1_SingleJet60_Prescl, &b_L1_SingleJet60_Prescl);
   fChain->SetBranchAddress("L1_SingleJet60_5bx", &L1_SingleJet60_5bx, &b_L1_SingleJet60_5bx);
   fChain->SetBranchAddress("L1_SingleMu0", &L1_SingleMu0, &b_L1_SingleMu0);
   fChain->SetBranchAddress("L1_SingleMu0_Prescl", &L1_SingleMu0_Prescl, &b_L1_SingleMu0_Prescl);
   fChain->SetBranchAddress("L1_SingleMu0_5bx", &L1_SingleMu0_5bx, &b_L1_SingleMu0_5bx);
   fChain->SetBranchAddress("L1_SingleMu10", &L1_SingleMu10, &b_L1_SingleMu10);
   fChain->SetBranchAddress("L1_SingleMu10_Prescl", &L1_SingleMu10_Prescl, &b_L1_SingleMu10_Prescl);
   fChain->SetBranchAddress("L1_SingleMu10_5bx", &L1_SingleMu10_5bx, &b_L1_SingleMu10_5bx);
   fChain->SetBranchAddress("L1_SingleMu14", &L1_SingleMu14, &b_L1_SingleMu14);
   fChain->SetBranchAddress("L1_SingleMu14_Prescl", &L1_SingleMu14_Prescl, &b_L1_SingleMu14_Prescl);
   fChain->SetBranchAddress("L1_SingleMu14_5bx", &L1_SingleMu14_5bx, &b_L1_SingleMu14_5bx);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu20_Prescl", &L1_SingleMu20_Prescl, &b_L1_SingleMu20_Prescl);
   fChain->SetBranchAddress("L1_SingleMu20_5bx", &L1_SingleMu20_5bx, &b_L1_SingleMu20_5bx);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu3_Prescl", &L1_SingleMu3_Prescl, &b_L1_SingleMu3_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_5bx", &L1_SingleMu3_5bx, &b_L1_SingleMu3_5bx);
   fChain->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   fChain->SetBranchAddress("L1_SingleMu5_Prescl", &L1_SingleMu5_Prescl, &b_L1_SingleMu5_Prescl);
   fChain->SetBranchAddress("L1_SingleMu5_5bx", &L1_SingleMu5_5bx, &b_L1_SingleMu5_5bx);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu7_Prescl", &L1_SingleMu7_Prescl, &b_L1_SingleMu7_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7_5bx", &L1_SingleMu7_5bx, &b_L1_SingleMu7_5bx);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo", &L1_SingleMuBeamHalo, &b_L1_SingleMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &L1_SingleMuBeamHalo_Prescl, &b_L1_SingleMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_5bx", &L1_SingleMuBeamHalo_5bx, &b_L1_SingleMuBeamHalo_5bx);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_Prescl", &L1_SingleMuOpen_Prescl, &b_L1_SingleMuOpen_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_5bx", &L1_SingleMuOpen_5bx, &b_L1_SingleMuOpen_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet10", &L1_SingleTauJet10, &b_L1_SingleTauJet10);
   fChain->SetBranchAddress("L1_SingleTauJet10_Prescl", &L1_SingleTauJet10_Prescl, &b_L1_SingleTauJet10_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet10_5bx", &L1_SingleTauJet10_5bx, &b_L1_SingleTauJet10_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet2", &L1_SingleTauJet2, &b_L1_SingleTauJet2);
   fChain->SetBranchAddress("L1_SingleTauJet2_Prescl", &L1_SingleTauJet2_Prescl, &b_L1_SingleTauJet2_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet2_5bx", &L1_SingleTauJet2_5bx, &b_L1_SingleTauJet2_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet20", &L1_SingleTauJet20, &b_L1_SingleTauJet20);
   fChain->SetBranchAddress("L1_SingleTauJet20_Prescl", &L1_SingleTauJet20_Prescl, &b_L1_SingleTauJet20_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet20_5bx", &L1_SingleTauJet20_5bx, &b_L1_SingleTauJet20_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet30", &L1_SingleTauJet30, &b_L1_SingleTauJet30);
   fChain->SetBranchAddress("L1_SingleTauJet30_Prescl", &L1_SingleTauJet30_Prescl, &b_L1_SingleTauJet30_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet30_5bx", &L1_SingleTauJet30_5bx, &b_L1_SingleTauJet30_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet4", &L1_SingleTauJet4, &b_L1_SingleTauJet4);
   fChain->SetBranchAddress("L1_SingleTauJet4_Prescl", &L1_SingleTauJet4_Prescl, &b_L1_SingleTauJet4_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet4_5bx", &L1_SingleTauJet4_5bx, &b_L1_SingleTauJet4_5bx);
   fChain->SetBranchAddress("L1_SingleTauJet50", &L1_SingleTauJet50, &b_L1_SingleTauJet50);
   fChain->SetBranchAddress("L1_SingleTauJet50_Prescl", &L1_SingleTauJet50_Prescl, &b_L1_SingleTauJet50_Prescl);
   fChain->SetBranchAddress("L1_SingleTauJet50_5bx", &L1_SingleTauJet50_5bx, &b_L1_SingleTauJet50_5bx);
   fChain->SetBranchAddress("L1_TripleJet14", &L1_TripleJet14, &b_L1_TripleJet14);
   fChain->SetBranchAddress("L1_TripleJet14_Prescl", &L1_TripleJet14_Prescl, &b_L1_TripleJet14_Prescl);
   fChain->SetBranchAddress("L1_TripleJet14_5bx", &L1_TripleJet14_5bx, &b_L1_TripleJet14_5bx);
   fChain->SetBranchAddress("L1_ZdcLooseVertex", &L1_ZdcLooseVertex, &b_L1_ZdcLooseVertex);
   fChain->SetBranchAddress("L1_ZdcLooseVertex_Prescl", &L1_ZdcLooseVertex_Prescl, &b_L1_ZdcLooseVertex_Prescl);
   fChain->SetBranchAddress("L1_ZdcLooseVertex_5bx", &L1_ZdcLooseVertex_5bx, &b_L1_ZdcLooseVertex_5bx);
   fChain->SetBranchAddress("L1_ZdcMinusOverThreshold", &L1_ZdcMinusOverThreshold, &b_L1_ZdcMinusOverThreshold);
   fChain->SetBranchAddress("L1_ZdcMinusOverThreshold_Prescl", &L1_ZdcMinusOverThreshold_Prescl, &b_L1_ZdcMinusOverThreshold_Prescl);
   fChain->SetBranchAddress("L1_ZdcMinusOverThreshold_5bx", &L1_ZdcMinusOverThreshold_5bx, &b_L1_ZdcMinusOverThreshold_5bx);
   fChain->SetBranchAddress("L1_ZdcPlusOverThreshold", &L1_ZdcPlusOverThreshold, &b_L1_ZdcPlusOverThreshold);
   fChain->SetBranchAddress("L1_ZdcPlusOverThreshold_Prescl", &L1_ZdcPlusOverThreshold_Prescl, &b_L1_ZdcPlusOverThreshold_Prescl);
   fChain->SetBranchAddress("L1_ZdcPlusOverThreshold_5bx", &L1_ZdcPlusOverThreshold_5bx, &b_L1_ZdcPlusOverThreshold_5bx);
   fChain->SetBranchAddress("L1_ZdcTightVertex", &L1_ZdcTightVertex, &b_L1_ZdcTightVertex);
   fChain->SetBranchAddress("L1_ZdcTightVertex_Prescl", &L1_ZdcTightVertex_Prescl, &b_L1_ZdcTightVertex_Prescl);
   fChain->SetBranchAddress("L1_ZdcTightVertex_5bx", &L1_ZdcTightVertex_5bx, &b_L1_ZdcTightVertex_5bx);
   fChain->SetBranchAddress("L1_ZeroBias_Ext", &L1_ZeroBias_Ext, &b_L1_ZeroBias_Ext);
   fChain->SetBranchAddress("L1_ZeroBias_Ext_Prescl", &L1_ZeroBias_Ext_Prescl, &b_L1_ZeroBias_Ext_Prescl);
   fChain->SetBranchAddress("L1_ZeroBias_Ext_5bx", &L1_ZeroBias_Ext_5bx, &b_L1_ZeroBias_Ext_5bx);
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
   fChain->SetBranchAddress("L1Tech_ZDC_loose_vertex.v0", &L1Tech_ZDC_loose_vertex_v0, &b_L1Tech_ZDC_loose_vertex_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_loose_vertex.v0_Prescl", &L1Tech_ZDC_loose_vertex_v0_Prescl, &b_L1Tech_ZDC_loose_vertex_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_loose_vertex.v0_5bx", &L1Tech_ZDC_loose_vertex_v0_5bx, &b_L1Tech_ZDC_loose_vertex_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_minus_over_threshold.v0", &L1Tech_ZDC_minus_over_threshold_v0, &b_L1Tech_ZDC_minus_over_threshold_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_minus_over_threshold.v0_Prescl", &L1Tech_ZDC_minus_over_threshold_v0_Prescl, &b_L1Tech_ZDC_minus_over_threshold_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_minus_over_threshold.v0_5bx", &L1Tech_ZDC_minus_over_threshold_v0_5bx, &b_L1Tech_ZDC_minus_over_threshold_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_plus_over_threshold.v0", &L1Tech_ZDC_plus_over_threshold_v0, &b_L1Tech_ZDC_plus_over_threshold_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_plus_over_threshold.v0_Prescl", &L1Tech_ZDC_plus_over_threshold_v0_Prescl, &b_L1Tech_ZDC_plus_over_threshold_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_plus_over_threshold.v0_5bx", &L1Tech_ZDC_plus_over_threshold_v0_5bx, &b_L1Tech_ZDC_plus_over_threshold_v0_5bx);
   fChain->SetBranchAddress("L1Tech_ZDC_tight_vertex.v0", &L1Tech_ZDC_tight_vertex_v0, &b_L1Tech_ZDC_tight_vertex_v0);
   fChain->SetBranchAddress("L1Tech_ZDC_tight_vertex.v0_Prescl", &L1Tech_ZDC_tight_vertex_v0_Prescl, &b_L1Tech_ZDC_tight_vertex_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_tight_vertex.v0_5bx", &L1Tech_ZDC_tight_vertex_v0_5bx, &b_L1Tech_ZDC_tight_vertex_v0_5bx);
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
