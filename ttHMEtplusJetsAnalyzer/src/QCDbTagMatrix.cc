#ifndef QCDBTAGMATRIX_CC
#define QCDBTAGMATRIX_CC

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/QCDbTagMatrix.h"

QCDbTagMatrix::QCDbTagMatrix( const string & higgsFileName, const string & hadronicTopFileName, const string & qcdFileName, TString suffix, TFile * outputFile, bool fillHistograms, const string & qcdHistoFileName, const int bJetNumCut ) : EventVariables( higgsFileName, hadronicTopFileName, qcdFileName, suffix+"_tagMatrix", outputFile, fillHistograms ) {
  bJetNumCut_ = bJetNumCut;

  inputFileSignal_ = new TFile(qcdHistoFileName.c_str());
  if ( inputFileSignal_ == 0 ) {
    cout << "ERROR: Cannot open qcd histograms file" << endl;
    exit(1);
  }

  // taggedJetEt_         = dynamic_cast<TH1F*>(inputFileSignal_->Get("taggedJetEt"));
  // taggedJetEta_        = dynamic_cast<TH1F*>(inputFileSignal_->Get("taggedJetEta"));
  // taggedJetS1_         = dynamic_cast<TH1F*>(inputFileSignal_->Get("taggedJetS1"));
  taggedJetTagMass_    = dynamic_cast<TH1F*>(inputFileSignal_->Get("taggedJetTagMass"));
  taggedJetDiscriminatorHighEff_ = dynamic_cast<TH1F*>(inputFileSignal_->Get("taggedJetDiscriminatorHighEff"));
  // notTaggedJetEt_      = dynamic_cast<TH1F*>(inputFileSignal_->Get("notTaggedJetEt"));
  // notTaggedJetEta_     = dynamic_cast<TH1F*>(inputFileSignal_->Get("notTaggedJetEta"));
  // notTaggedJetS1_      = dynamic_cast<TH1F*>(inputFileSignal_->Get("notTaggedJetS1"));
  notTaggedJetTagMass_ = dynamic_cast<TH1F*>(inputFileSignal_->Get("notTaggedJetTagMass"));
  notTaggedJetDiscriminatorHighEff_ = dynamic_cast<TH1F*>(inputFileSignal_->Get("notTaggedJetDiscriminatorHighEff"));
}

QCDbTagMatrix::~QCDbTagMatrix() {
  inputFileSignal_->Close();
}

void QCDbTagMatrix::multiply( const vector<const OfflineJet *> jetCollection, const OfflineMEt * offlineMEt ) {

  double jetsNumber = jetCollection.size();

  if ( jetsNumber > 8 ) {
    cout << "Error: use <=8 jets, or the number of combinations will be high" << endl;
    exit(1);
  }


  // Extract the probability of each jet to be tagged or not and store it in arrays
  // ------------------------------------------------------------------------------
  double * probTagged = new double[int(jetsNumber)];
  double * probNotTagged = new double[int(jetsNumber)];
  for( int i=0; i<jetsNumber; ++i ) {
    probTagged[i] = 0;
    probNotTagged[i] = 0;
  }

  vector<const OfflineJet *>::const_iterator goodJetIt = jetCollection.begin();
  unsigned int goodJetNum = 0;
  for ( ; goodJetIt != jetCollection.end(); ++goodJetIt, ++goodJetNum ) {
    int numTkS1 = 0;
    edm::RefVector<std::vector<SimpleTrack> > tkRefVec((*goodJetIt)->tkRefVec());
    edm::RefVector<std::vector<SimpleTrack> >::const_iterator tkRefIt = tkRefVec.begin();
    for ( ; tkRefIt != tkRefVec.end(); ++tkRefIt ) {
      if ( (*tkRefIt)->ip2Dsignificance() >= 1. ) ++numTkS1;
    }

    double jetEt = (*goodJetIt)->et();
    double jetEta = (*goodJetIt)->eta();
    // Determine bin index
    int etId = int(jetEt/qcdBinSize_[0]);
    int etaId = int(fabs(jetEta)/qcdBinSize_[1]);
    int numTkS1Id = int(numTkS1/qcdBinSize_[2]);
    if ( etId < 0 || etaId < 0 || numTkS1Id < 0 ) cout << "Error: index < 0, will crash..." << endl;
    if ( etId>int(qcdBinNum_[0]-1) ) etId = qcdBinNum_[0]-1;
    if ( etaId>int(qcdBinNum_[1]-1) ) etaId = qcdBinNum_[1]-1;
    if ( numTkS1Id>int(qcdBinNum_[2]-1) ) numTkS1Id = qcdBinNum_[2]-1;

    probTagged[goodJetNum] = taggedJet_[etId][etaId][numTkS1Id];
    probNotTagged[goodJetNum] = notTaggedJet_[etId][etaId][numTkS1Id];
    cout << "goodJetNum = " << goodJetNum << endl;
    cout << "taggedJet_["<<etId<<"]["<<etaId<<"]["<<numTkS1Id<<"] = " << taggedJet_[etId][etaId][numTkS1Id] << endl;
    cout << "notTaggedJet_["<<etId<<"]["<<etaId<<"]["<<numTkS1Id<<"] = " << notTaggedJet_[etId][etaId][numTkS1Id] << endl;
  }

  // Copy the collection of jets in a modifiable one
  vector<OfflineJet> tempJets;
  vector<const OfflineJet *>::const_iterator jet = jetCollection.begin();
  for( ; jet != jetCollection.end(); ++jet ) {
    tempJets.push_back(**jet);
  }

  // Build all the combinations
  // --------------------------
  for( int comb=0; comb<pow(2, jetsNumber); ++comb ) {

    // To multiply by 2 or to move the bit on the left is the same in this case
    // for ( int check=1; check<=128 && check <= comb; check = check << 1) {
    // for ( int check=1; check<=128 && check <= comb; check*=2 ) {

    // Probability of this combination
    double prob = 1.;

    // New collection of goodJets to be passed to the fill method
    vector<const OfflineJet *> goodJets;
    // New collection of bJets to be passed to the fill method
    vector<const OfflineJet *> goodbJets;

    // Count the number of b-jets in the combination
    int bJetNum = 0;

    int jetNum = 0;
    vector<OfflineJet>::iterator jet = tempJets.begin();
    for( int check=1; jet != tempJets.end(); ++jet, check*=2, ++jetNum ) {

      // bitwise AND between the two integers
      // check acts as a mask to check for individual bits in the
      // comb int. check has always only one value as 1 and the
      // others are all 0.
      if ( comb & check ) {
        // this jet is b-tagged
//         cout << "value match for: comb = " << comb
//              << "; and check = " << check << endl; 
        // Extract new discriminant value from the b-tagged histogram
        jet->setDiscriminatorHighEff(taggedJetDiscriminatorHighEff_->GetRandom());
        jet->setbTagTkInvMass(taggedJetTagMass_->GetRandom());
        goodJets.push_back(&(*jet));
        goodbJets.push_back(&(*jet));
        ++bJetNum;
        prob *= probTagged[jetNum];
      }
      else {
        // Extract new discriminant value from the nob-tagged histogram
        jet->setDiscriminatorHighEff(notTaggedJetDiscriminatorHighEff_->GetRandom());
        jet->setbTagTkInvMass(notTaggedJetTagMass_->GetRandom());
        goodJets.push_back(&(*jet));
        prob *= probNotTagged[jetNum];
      }
    }
    // Call the fill method to store all the values of this pseudoevent.
    // Pass the probability as a weight.
    // Only fill events that pass the cut on the number of b-jets.
    cout << "bJetNum["<<comb<<"] = " << bJetNum << " with probability = " << prob << endl;
    if ( bJetNum >= bJetNumCut_ ) {
      fill( goodJets, goodbJets, offlineMEt, prob );
    }
  }


  delete[] probTagged;
  delete[] probNotTagged;

}

#endif // QCDBTAGMATRIX_CC
