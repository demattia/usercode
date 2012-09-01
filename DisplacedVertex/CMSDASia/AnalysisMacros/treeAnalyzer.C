#define treeAnalyzer_cxx
#include "treeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMath.h>

void treeAnalyzer::initializeCuts()
{
  if( electrons_ ) {
    //    ptCut_ = 38.;
    ptCut_ = 51.;
    d0SignificanceCut_ = 3.;
    dPhiCorrCut_ = 0.8;
    decayLengthSignificance2DCut_ = 8.;
  }
  else {
    ptCut_ = 26.;
    //     ptCut_ = 33.;
    d0SignificanceCut_ = 2.;
    dPhiCorrCut_ = 0.2;
    decayLengthSignificance2DCut_ = 5.;
  }
}

bool treeAnalyzer::acceptanceCuts( std::vector<TreeCandidate>::const_iterator cand )
{
  if( fabs(cand->leptonEtaH) < 2. &&
      fabs(cand->leptonEtaL) < 2. &&
      cand->leptonPtL > ptCut_ &&
      cand->leptonPtH > ptCut_
  ) return true;
  return false;
}

bool treeAnalyzer::trackSelectionCuts( std::vector<TreeCandidate>::const_iterator cand, const bool removeIsolationCut )
{
  if( // The PATtuples already select high purity tracks are high purity
      ( cand->leptonIsoL < 4 && cand->leptonIsoH < 4 ) || removeIsolationCut
  ) return true;
  return false;
}

bool treeAnalyzer::lifetimeRelatedCuts( std::vector<TreeCandidate>::const_iterator cand, const bool removeLifetimeRelatedCuts, const bool decayLengthSigniInverted)
{
  // std::cout<<removeLifetimeRelatedCuts<<","<<decayLengthSigniInverted<<std::endl;
  if(decayLengthSigniInverted==1){
    if((removeLifetimeRelatedCuts || cand->leptonAbsD0SignificanceH > d0SignificanceCut_) &&
        (removeLifetimeRelatedCuts || cand->leptonAbsD0SignificanceL > d0SignificanceCut_) &&
        (removeLifetimeRelatedCuts || cand->dPhiCorr < dPhiCorrCut_) && // this is the vertex colinearity cut
        //   (removeLifetimeRelatedCuts || decayLengthSignificance < decayLengthSignificance2DCut_)//inverted
        ( !removeLifetimeRelatedCuts||cand->decayLengthSignificance < decayLengthSignificance2DCut_)
    ) return true;
    return false;
  }
  else {
    if((removeLifetimeRelatedCuts || cand->leptonAbsD0SignificanceH > d0SignificanceCut_) &&
        (removeLifetimeRelatedCuts || cand->leptonAbsD0SignificanceL > d0SignificanceCut_) &&
        (removeLifetimeRelatedCuts || cand->dPhiCorr < dPhiCorrCut_) && // this is the vertex colinearity cut
        //     (removeLifetimeRelatedCuts || decayLengthSignificance > decayLengthSignificance2DCut_)
        (removeLifetimeRelatedCuts || cand->decayLengthSignificance > decayLengthSignificance2DCut_)
    ) return true;
    return false;
  }
}
bool treeAnalyzer::dileptonSelectionCuts( std::vector<TreeCandidate>::const_iterator cand )
{
  if((cand->leptonChargeL*cand->leptonChargeH < 0) &&
      cand->vertexChi2 < 5 &&
      cand->hitsBeforeVertexL <= 1 && cand->hitsBeforeVertexH <= 1 &&
      (electrons_ || cand->deltaR > 0.2 ) && // muons only
      (electrons_ || cand->cosine > -0.95 ) && // muons only
      !(cand->isStandAloneL && cand->isStandAloneH) && // we are not using standAloneMuons
      ( electrons_ || cand->corrDileptonMass > 15.)&& ( !electrons_ || cand->caloCorrMass > 15.)
  ) return true;
  return false;
}

bool treeAnalyzer::triggerMatching( std::vector<TreeCandidate>::const_iterator cand )
{
  if( cand->triggerMatchL != 0 && cand->triggerMatchH != 0 && cand->triggerMatchL != cand->triggerMatchH ) return true;
  return false;
}

bool treeAnalyzer::analysisCuts(std::vector<TreeCandidate>::const_iterator cand, const bool removeIsolationCut, const bool removeLifetimeRelatedCuts, const bool decayLengthSigniInverted )
{
  if( acceptanceCuts( cand ) &&
      trackSelectionCuts( cand, removeIsolationCut ) &&  // i.e. isolation cuts
      dileptonSelectionCuts( cand ) &&
      lifetimeRelatedCuts( cand, removeLifetimeRelatedCuts,decayLengthSigniInverted)
      &&   triggerMatching( cand )
  ) return true;
  return false;
}


void treeAnalyzer::Loop()
{ 

  //  if(electron_)
  //  std::cout<<"hello"<<std::endl;
  if (fChain == 0) return;

  initializeCuts();
  //   std::cout<<"hello"<<std::endl;
  Long64_t nentries = fChain->GetEntriesFast();

  TH1::SetDefaultSumw2();

  // nPV plots
  TH1F histoNPv_nolifetime_nodecayLength("nPV_nolifetime_nodecaylength","nPV_nolifetime_nodecaylength",60,0,60);
  TH1F histoNvtx_true_nolifetime_nodecayLength("nvtx_true_nolifetime_nodecayLength","nvtx_true_nolifetime_nodecayLength",60,0,60);

  double averagePUWeight = 0;
  unsigned int nPUWeights = 0;

  // Mass & deltaPhi histograms
  TH1F histo_nolifetime_inverted("Mass_nolifetime_inverted", "Mass_nolifetime_inverted", 75, 0, 500);
  TH1F histoDeltaPhi_nolifetime_inverted("deltaPhi_nolifetime_inverted", "deltaPhi_nolifetime_inverted", 100, 0, 3.2);

  TH1F histo_nolifetime_nodecaylength("Mass_nolifetime_nodecaylength", "Mass_nolifetime_nodecaylength", 100, 0, 500);
  TH1F histoDeltaPhi_nolifetime_nodecaylength("deltaPhi_nolifetime_nodecaylength", "deltaPhi_nolifetime_nodecaylength", 100, 0, 3.2);

  TH1F histo_lifetime_nodecaylength("Mass_lifetime_nodecaylength", "Mass_lifetime_nodecaylength", 100, 0, 500);
  TH1F histoDeltaPhi_lifetime_nodecaylength("deltaPhi_lifetime_nodecaylength", "deltaPhi_lifetime_nodecaylength", 100, 0, 3.2);

  TH1F histo_lifetime_notinverted("Mass_lifetime_notinverted", "Mass_lifetime_notinverted", 100, 0, 500);
  TH1F histoDeltaPhi_lifetime_notinverted("deltaPhi_lifetime_notinverted", "deltaPhi_lifetime_notinverted", 100, 0, 3.2); 


  // Decay Length Histograms
  double maxDecayLengthSignificance = 20;
  double maxDecayLength = 40;

  int nBins = 40;
  //  int nBins = 100;
  if( nBins%2 != 0 ) {
    std::cout << "Error: please set the binning of histoSignedLifetime to an even number." << std::endl;
  }
  TH1F histoSignedDecayLengthSignificance_nolifetime_inverted("signedDecayLengthSignificance_nolifetime_inverted", "signedDecayLengthSignificance_nolifetime_inverted", nBins, -maxDecayLengthSignificance, maxDecayLengthSignificance);
  TH1F histoSignedDecayLength_nolifetime_inverted("signedDecayLength_nolifetime_inverted", "signedDecayLength_nolifetime_inverted", nBins, -maxDecayLength, maxDecayLength);
  TH1F diff("diff", "diff", nBins/2, 0., maxDecayLengthSignificance);
  TH1F histoGenDecayLength_nolifetime_inverted("genDecayLength_nolifetime_inverted", "genDecayLength_nolifetime_inverted", 1000, -1, 1);


  TH1F histoSignedDecayLengthSignificance_nolifetime_nodecaylength("signedDecayLengthSignificance_nolifetime_nodecaylength", "signedDecayLengthSignificance_nolifetime_nodecaylength", nBins, -maxDecayLengthSignificance, maxDecayLengthSignificance);
  TH1F histoSignedDecayLength_nolifetime_nodecaylength("signedDecayLength_nolifetime_nodecaylength", "signedDecayLength_nolifetime_nodecaylength", nBins, -maxDecayLength, maxDecayLength);
  TH1F histoGenDecayLength_nolifetime_nodecaylength("genDecayLength_nolifetime_nodecaylength", "genDecayLength_nolifetime_nodecaylength", 1000, -1, 1);

  TH1F histoSignedDecayLengthSignificance_lifetime_notinverted("signedDecayLengthSignificance_lifetime_notinverted", "signedDecayLengthSignificance_lifetime_notinverted", nBins, -maxDecayLengthSignificance, maxDecayLengthSignificance);
  TH1F histoSignedDecayLength_lifetime_notinverted("signedDecayLength_lifetime_notinverted", "signedDecayLength_lifetime_notinverted", nBins, -maxDecayLength, maxDecayLength);
  TH1F histoGenDecayLength_lifetime_notinverted("genDecayLength_lifetime_notinverted", "genDecayLength_lifetime_notinverted", 1000, -1, 1);

  TH1F histoSignedDecayLengthSignificance_lifetime_nodecaylength("signedDecayLengthSignificance_lifetime_nodecaylength", "signedDecayLengthSignificance_lifetime_nodecaylength", nBins, -maxDecayLengthSignificance, maxDecayLengthSignificance);
  TH1F histoSignedDecayLength_lifetime_nodecaylength("signedDecayLength_lifetime_nodecaylength", "signedDecayLength_lifetime_nodecaylength", nBins, -maxDecayLength, maxDecayLength);
  TH1F histoGenDecayLength_lifetime_nodecaylength("genDecayLength_lifetime_nodecaylength", "genDecayLength_lifetime_nodecaylength", 1000, -1, 1);

  // Isolation histograms
  TH1F histoIsolationPtH_nolifetime_inverted("isolationPtH_nolifetime_inverted","IsolationPtH_nolifetime_inverted",50,0,25);
  TH1F histoIsolationPtL_nolifetime_inverted("isolationPtL_nolifetime_inverted","IsolationPtL_nolifetime_inverted",50,0,25);

  // Back-to-back
  TH1F histo_cosine_nolifetime_inverted("cosine_nolifetime_inverted","cosine_nolifetime_inverted", 44, -1.1, 1.1 );

  // Vertex Chi2
  TH1F histo_vertexChi2_nolifetime_inverted("vertexChi2_nolifetime_inverted","vertexChi2_nolifetime_inverted",40,0,20);

  // Lepton Pt
  TH1F histoPtH_nolifetime_inverted("ptH_nolifetime_inverted","ptH_nolifetime_inverted",100,0,200);
  TH1F histoPtL_nolifetime_inverted("ptL_nolifetime_inverted","ptL_nolifetime_inverted",100,0,200);

  int total = 0;
  // Loop over events
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    ++total;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Get PU weight for this event
    double totalWeight = 0;
    if ( dirName_.Contains("Data") ) {
      totalWeight = weight_;
    }
    else {
      double puweight = puweights_.weight( candidates->nvtx_true );

      //      double puweight = 1;
      totalWeight = puweight * weight_;
    }

    // Loop over candidates in event
    std::vector<TreeCandidate>::const_iterator cand = candidates->candidates_.begin();
    for( ; cand != candidates->candidates_.end(); ++cand ) {

      if( cand->decayLengthSignificance == -999 ) continue;

      // Apply analysis cuts.
      double Mass;
      if (electrons_) Mass = cand->caloCorrMass;
      else Mass = cand->corrDileptonMass;
      if( analysisCuts(cand,false,true,true) ) {

        // Pileup histos
        histoNPv_nolifetime_nodecayLength.Fill(candidates->numPV, totalWeight);
        histoNvtx_true_nolifetime_nodecayLength.Fill(candidates->nvtx_true, totalWeight);
        nPUWeights++;
        averagePUWeight += totalWeight / weight_;

        histo_nolifetime_inverted.Fill(Mass, totalWeight);
        histoDeltaPhi_nolifetime_inverted.Fill(cand->dPhiCorr, totalWeight);
        histo_cosine_nolifetime_inverted.Fill(cand->cosine,totalWeight);
        histo_vertexChi2_nolifetime_inverted.Fill(cand->vertexChi2,totalWeight);

        histoPtL_nolifetime_inverted.Fill(cand->leptonPtL, totalWeight);
        histoPtH_nolifetime_inverted.Fill(cand->leptonPtH, totalWeight);

        int sign = 1;
        if( cand->dPhiCorr > TMath::PiOver2() ) sign = -1;
        histoSignedDecayLengthSignificance_nolifetime_inverted.Fill(sign*cand->decayLengthSignificance, totalWeight);
        histoSignedDecayLength_nolifetime_inverted.Fill(sign*cand->decayLength, totalWeight);

        histoGenDecayLength_nolifetime_inverted.Fill(cand->genDecayLength2D, totalWeight);
      }

      if( analysisCuts(cand,false,true,false) ) {

        histo_nolifetime_nodecaylength.Fill(Mass, totalWeight);
        histoDeltaPhi_nolifetime_nodecaylength.Fill(cand->dPhiCorr, totalWeight);
        int sign = 1;
        if( cand->dPhiCorr > TMath::PiOver2() ) sign = -1;
        histoSignedDecayLengthSignificance_nolifetime_nodecaylength.Fill(sign*cand->decayLengthSignificance, totalWeight);
        histoSignedDecayLength_nolifetime_nodecaylength.Fill(sign*cand->decayLength, totalWeight);
        histoGenDecayLength_nolifetime_nodecaylength.Fill(cand->genDecayLength2D, totalWeight);
      }

      if( analysisCuts(cand,false,false,true) ) {

        histo_lifetime_nodecaylength.Fill(Mass, totalWeight);
        histoDeltaPhi_lifetime_nodecaylength.Fill(cand->dPhiCorr, totalWeight);
        int sign = 1;
        if( cand->dPhiCorr > TMath::PiOver2() ) sign = -1;
        histoSignedDecayLengthSignificance_lifetime_nodecaylength.Fill(sign*cand->decayLengthSignificance, totalWeight);
        histoSignedDecayLength_lifetime_nodecaylength.Fill(sign*cand->decayLength, totalWeight);

        histoGenDecayLength_lifetime_nodecaylength.Fill(cand->genDecayLength2D, totalWeight);
      }


      if( analysisCuts(cand,false,false,false) ) {

        histo_lifetime_notinverted.Fill(Mass, totalWeight);
        histoDeltaPhi_lifetime_notinverted.Fill(cand->dPhiCorr, totalWeight);
        int sign = 1;
        if( cand->dPhiCorr > TMath::PiOver2() ) sign = -1;
        histoSignedDecayLengthSignificance_lifetime_notinverted.Fill(sign*cand->decayLengthSignificance, totalWeight);
        histoSignedDecayLength_lifetime_notinverted.Fill(sign*cand->decayLength, totalWeight);
        histoGenDecayLength_lifetime_notinverted.Fill(cand->genDecayLength2D, totalWeight);

      }

      if ( analysisCuts(cand,true,true,true) ) {
        histoIsolationPtH_nolifetime_inverted.Fill(cand->leptonIsoH, totalWeight);
        histoIsolationPtL_nolifetime_inverted.Fill(cand->leptonIsoL, totalWeight);
      }



      //    std::vector<std::string>::const_iterator it = triggers->begin();
      //    for( ; it != triggers->end(); ++it ) {
      //	 std::cout << "trigger name = " << *it << std::endl;
      //   }


    }
  }

  std::cout << "total entries = " << total << std::endl;

  std::cout << "Average pu weight = " << averagePUWeight << " / " << nPUWeights << " = " << averagePUWeight/float(nPUWeights) << endl;

  TString type("_muons");
  if( electrons_ ) type = "_electrons";
  TFile outputFile("WeightedFiles/"+dirName_+"_weighted"+type+".root", "RECREATE");
  outputFile.cd();

  int overflowvalue = 1; 

  //  int  N= histo_nolifetime_inverted.GetNbinsX();
  //  histo_nolifetime_inverted.SetBinContent(N,histo_nolifetime_inverted.GetBinContent(N)+histo_nolifetime_inverted.GetBinContent(N+1));
  if (overflowvalue==1){
    overflow( & histo_nolifetime_inverted );
    overflow( & histoDeltaPhi_nolifetime_inverted );
    overflow( & histo_cosine_nolifetime_inverted );
    overflow( & histo_vertexChi2_nolifetime_inverted );
    overflow( & histoPtL_nolifetime_inverted );
    overflow( & histoPtH_nolifetime_inverted );
    overflow( & histoSignedDecayLengthSignificance_nolifetime_inverted );
    overflow( & histoSignedDecayLength_nolifetime_inverted );
    overflow( & histoGenDecayLength_nolifetime_inverted );
    overflow( & histo_nolifetime_nodecaylength );
    overflow( & histoNPv_nolifetime_nodecayLength );
    overflow( & histoNvtx_true_nolifetime_nodecayLength );
    overflow( & histoDeltaPhi_nolifetime_nodecaylength );
    overflow( & histoSignedDecayLengthSignificance_nolifetime_nodecaylength );
    overflow( & histoSignedDecayLength_nolifetime_nodecaylength );
    overflow( & histoGenDecayLength_nolifetime_nodecaylength );
    overflow( & histo_lifetime_nodecaylength );
    overflow( & histoDeltaPhi_lifetime_nodecaylength );
    overflow( & histoSignedDecayLengthSignificance_lifetime_nodecaylength );
    overflow( & histoSignedDecayLength_lifetime_nodecaylength );
    overflow( & histoGenDecayLength_lifetime_nodecaylength );
    overflow( & histo_lifetime_notinverted );
    overflow( & histoDeltaPhi_lifetime_notinverted );
    overflow( & histoSignedDecayLengthSignificance_lifetime_notinverted );
    overflow( & histoSignedDecayLength_lifetime_notinverted );
    overflow( & histoGenDecayLength_lifetime_notinverted );
    overflow( & histoIsolationPtH_nolifetime_inverted );
    overflow( & histoIsolationPtL_nolifetime_inverted );
  }


  histo_nolifetime_inverted.Write();
  histoDeltaPhi_nolifetime_inverted.Write();
  histo_cosine_nolifetime_inverted.Write();
  histo_vertexChi2_nolifetime_inverted.Write();
  histoPtL_nolifetime_inverted.Write();
  histoPtH_nolifetime_inverted.Write();
  histoSignedDecayLengthSignificance_nolifetime_inverted.Write();
  histoSignedDecayLength_nolifetime_inverted.Write();
  histoGenDecayLength_nolifetime_inverted.Write();

  histo_nolifetime_nodecaylength.Write();
  histoNPv_nolifetime_nodecayLength.Write();
  histoNvtx_true_nolifetime_nodecayLength.Write();
  histoDeltaPhi_nolifetime_nodecaylength.Write();
  histoSignedDecayLengthSignificance_nolifetime_nodecaylength.Write();
  histoSignedDecayLength_nolifetime_nodecaylength.Write();
  histoGenDecayLength_nolifetime_nodecaylength.Write();

  histo_lifetime_nodecaylength.Write();
  histoDeltaPhi_lifetime_nodecaylength.Write();
  histoSignedDecayLengthSignificance_lifetime_nodecaylength.Write();
  histoSignedDecayLength_lifetime_nodecaylength.Write();
  histoGenDecayLength_lifetime_nodecaylength.Write();

  histo_lifetime_notinverted.Write();
  histoDeltaPhi_lifetime_notinverted.Write();
  histoSignedDecayLengthSignificance_lifetime_notinverted.Write();
  histoSignedDecayLength_lifetime_notinverted.Write();
  histoGenDecayLength_lifetime_notinverted.Write();

  histoIsolationPtH_nolifetime_inverted.Write();
  histoIsolationPtL_nolifetime_inverted.Write();

  //  int  N= histo_nolifetime_inverted->GetNbinsX();
  //  histo_nolifetime_inverted->SetBinContent(N,GetBinContent(N)+GetBinContent(N+1));


  // Computing the difference between negative and positive decay legnth significance
  for(int i=1; i<=nBins/2; ++i) {
    // std::cout << "histoSignedDecayLengthSignificance->GetBinContent("<<i<<") = " << histoSignedDecayLengthSignificance->GetBinContent(i) << std::endl;
    // std::cout << "histoSignedDecayLengthSignificance->GetBinContent("<<nBins+1-i<<") = " << histoSignedDecayLengthSignificance->GetBinContent(nBins+1-i) << std::endl;
    diff.SetBinContent(i, histoSignedDecayLengthSignificance_nolifetime_inverted.GetBinContent(i)/histoSignedDecayLengthSignificance_nolifetime_inverted.GetBinContent(nBins+1-i));
  }
  diff.Write();

  outputFile.Write();
  outputFile.Close();

}
