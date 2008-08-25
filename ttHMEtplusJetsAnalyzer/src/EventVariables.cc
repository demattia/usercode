#ifndef EVENTVARIABLES_CC
#define EVENTVARIABLES_CC

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/EventVariables.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"

// Defined in a nameless namespace so that they only exist in this file
namespace {
  bool sortParticlesByProbability( const pair<double, Particle<const OfflineJet> > & a, const pair<double, Particle<const OfflineJet> > & b ) {
    return a.first < b.first;
  }
  bool sortOfflineJetPtr( const OfflineJet * a, const OfflineJet * b ) {
    return a->et() < b->et();
  }
  template <class T> class SortByMassDifference {
    double referenceMass_;
  public:
    SortByMassDifference( const double & referenceMass ) { referenceMass_ = referenceMass; }
    bool operator()( const T & p1, const T & p2 ) {
      return fabs( p1.mass() - referenceMass_ ) < fabs( p2.mass() - referenceMass_ );
    }
  };
}

EventVariables::EventVariables( const string & higgsFileName, const string & hadronicTopFileName, const string & qcdFileName, TString suffix, TFile * outputFile, bool fillHistograms ) {

  fillHistograms_ = fillHistograms;

  referenceHiggsMass_ = 120;
  referenceTopMass_ = 172.6;
  referenceWmass_ = 80.4;

  // Load the matrices for Higgs, hadronic top and QCD b-tag probabiliy from the files
  // ---------------------------------------------------------------------------------
  fillProbabilityMatrices( higgsFileName, higgsBinNum_, higgsBinSize_, trueH_, falseH_ );
  fillProbabilityMatrices( hadronicTopFileName, hadronicTopBinNum_, hadronicTopBinSize_, trueHadronicTop_, falseHadronicTop_ );
  fillProbabilityMatrices( qcdFileName, qcdBinNum_, qcdBinSize_, taggedJet_, notTaggedJet_ );

  outputFile_ = outputFile;
  TString dirName = ("EventVariables");

  if ( suffix != "" ) {
    suffix.Prepend("_");
    dirName.Append(suffix);
  }

  outputDir_ = outputFile_->mkdir(dirName);
  outputDir_->cd();

  if ( fillHistograms_ ) {

    higgsMass_       = new TH1D( "higgsMass" + suffix, "reconstructed higgs mass" + suffix, 50, 0, 300 );
    hadronicTopMass_ = new TH1D( "hadronicTopMass" + suffix, "reconstructed hadronic top mass" + suffix, 50, 0, 800 );
    hadronicWmass_   = new TH1D( "hadronicWmass" + suffix, "reconstructed W from the hadronic top" + suffix, 50, 0, 200 );
    chi2ofMasses_    = new TH1D( "chi2ofMasses" + suffix, "chi2 computed with Higgs, Top and W masses" + suffix, 50, 0, 200 );

    firstNjetsMass_[0]        = new TH1D( "first6jetsMass" + suffix, "Mass reconstructed with the first 6 jets" + suffix, 50, 0, 2000 );
    firstNjetsCentrality_[0]  = new TH1D( "first6jetsCentrality" + suffix, "Centrality reconstructed with the first 6 jets" + suffix, 50, 0, 10 );
    firstNjetsMass_[1]        = new TH1D( "first8jetsMass" + suffix, "Mass reconstructed with the first 8 jets" + suffix, 50, 0, 2000 );
    firstNjetsCentrality_[1]  = new TH1D( "first8jetsCentrality" + suffix, "Centrality reconstructed with the first 8 jets" + suffix, 50, 0, 10 );
    hadronicTopProjectionAlongHiggsDirection_ = new TH1D( "hadronicTopProjectionAlongHiggsDirection" + suffix, "Hadronic Top projection along Higgs direction" + suffix, 50, 0, 500 );
    deltaEtaHadronicTopHiggs_ = new TH1D( "deltaEtaHadronicTopHiggs" + suffix, "deltaEta hadronicTop-Higgs" + suffix, 50, -3., 3. );
    goodHt_                   = new TH1D( "goodHt" + suffix, "Ht evaluated with all selected jets" + suffix, 50, 0, 1000 );
    mEtSig_                   = new TH1D( "mEtSig" + suffix, "Missing Et significance" + suffix, 50, 0, 10 );
    deltaPhiMEtNthLeadingJet_[0] = new TH1D( "deltaPhiMEt1stLeadingJet" + suffix, "deltaPhi between MEt and 1st leading jet" + suffix, 50, 0., 3.14 );
    deltaPhiMEtNthLeadingJet_[1] = new TH1D( "deltaPhiMEt2ndLeadingJet" + suffix, "deltaPhi between MEt and 2nd leading jet" + suffix, 50, 0., 3.14 );
    deltaPhiMEtNthLeadingJet_[2] = new TH1D( "deltaPhiMEt3rdLeadingJet" + suffix, "deltaPhi between MEt and 3rd leading jet" + suffix, 50, 0., 3.14 );
    hadronicTopPlusHiggsMass_ = new TH1D( "hadronicTopPlusHiggsMass" + suffix, "Mass of the 5 jets from hadronic Top and Higgs" + suffix, 50, 0., 1000. );
    sumHighEffDiscriminantFirst4Jets_ = new TH1D( "sumHighEffDiscriminantFirst4Jets" + suffix, "Sum of the high efficiency discriminant of the 4 most energetic jets" + suffix, 50, 0., 100. );
    sumHighEffDiscriminantFirst6Jets_ = new TH1D( "sumHighEffDiscriminantFirst6Jets" + suffix, "Sum of the high efficiency discriminant of the 6 most energetic jets" + suffix, 50, 0., 100. );
    remainingJetsMass_ = new TH1D( "remainingJetsMass" + suffix, "Mass of the remaining jets after removal of Higgs and hadronic Top jets" + suffix, 50, 0., 1000. );
    sixthJetEt_ = new TH1D( "sixthJetEt" + suffix, "Et of the sixth jet" + suffix, 50, 0., 200 );
    bTagTkInvMass_ = new TH1D( "bTagTkInvMass" + suffix, "Sum of the mass of all b-tagged jets" + suffix, 50, 0., 200. );
  }
}

vector<double> EventVariables::fill( vector<const OfflineJet *> jetCollection, const vector<const OfflineJet *> & bTaggedJetCollection, const OfflineMEt * offlineMEt ) {

  // Empty the vector from any previous event variables
  eventVariablesVector_.clear();

  // First of all sort the collection in Et. The first is the most energetic
  // Do this here, not inside the method to do it only once.
  sort( jetCollection.rbegin(), jetCollection.rend(), sortOfflineJetPtr );

  // This also fills the corresponding histograms
  Particle<const OfflineJet> first6jets( firstNjetsParticle( jetCollection, 6 ) );
  Particle<const OfflineJet> first8jets( firstNjetsParticle( jetCollection, 8 ) );

  // Evaluates variables on all "good" jets and fills corresponding histograms
  allGoodJetsVariables( jetCollection, offlineMEt);

  // Higgs, hadronic Top and W determination
  // ---------------------------------------

  // Create pairs of b-jets and evaluate their probability to come from the Higgs decay
  //  vector<pair<true/false ratio, candidate> >
  vector<pair<double, Particle<const OfflineJet> > > bTaggedPairs;
  vector<const OfflineJet *>::const_iterator bTaggedJetIt = bTaggedJetCollection.begin();
  for ( ; bTaggedJetIt != bTaggedJetCollection.end(); ++bTaggedJetIt ) {
    vector<const OfflineJet *>::const_iterator subbTaggedJetIt = bTaggedJetIt+1;
    for ( ; subbTaggedJetIt != bTaggedJetCollection.end(); ++subbTaggedJetIt ) {
      // bTaggedPairs.push_back(pairStruct( *bTaggedJetIt, *subbTaggedJetIt ));
      Particle<const OfflineJet> higgsCandidate( *bTaggedJetIt );
      higgsCandidate.add( *subbTaggedJetIt );
      bTaggedPairs.push_back( make_pair( evalHiggsPairProbability(higgsCandidate), higgsCandidate ) ); 
      // cout << "higgsCandidate true/false ratio = " << bTaggedPairs.back().first << endl;
    }
  }
  sort( bTaggedPairs.rbegin(), bTaggedPairs.rend(), sortParticlesByProbability );
  Particle<const OfflineJet> * selectedHiggs = &(bTaggedPairs.front().second);

  eventVariablesVector_.push_back(selectedHiggs->mass());
  if ( fillHistograms_ ) higgsMass_->Fill(selectedHiggs->mass());

  // Remove the Higgs associated jets
  jetCollection.erase( find(jetCollection.begin(), jetCollection.end(), (selectedHiggs->components())[0]) );
  jetCollection.erase( find(jetCollection.begin(), jetCollection.end(), (selectedHiggs->components())[1]) );

  // If there are at least 3 jets, once those from the selected Higgs have been removed, search for hadronic top candidates
  if ( jetCollection.size() >= 3 ) {
    // Once the Higgs is selected use the remaining jets to select the hadronic top triplet
    vector<pair<double, Particle<const OfflineJet> > > hadronicTopTriplet;
    vector<const OfflineJet *>::const_iterator jetIt = jetCollection.begin();
    for ( ; jetIt != jetCollection.end(); ++jetIt ) {
      vector<const OfflineJet *>::const_iterator subJetIt = jetIt+1;
      for ( ; subJetIt != jetCollection.end(); ++subJetIt ) {
        vector<const OfflineJet *>::const_iterator subSubJetIt = subJetIt+1;
        for ( ; subSubJetIt != jetCollection.end(); ++subSubJetIt ) {
          Particle<const OfflineJet> hadronicTopCandidate( *jetIt );
          hadronicTopCandidate.add( *subJetIt );
          hadronicTopCandidate.add( *subSubJetIt );
          hadronicTopTriplet.push_back( make_pair( evalTopTripletProbability(hadronicTopCandidate, *selectedHiggs), hadronicTopCandidate ) ); 
          // cout << "hadronicTopCandidate true/false ratio = " << hadronicTopTriplet.back().first << endl;
        }
      }
    }
    sort( hadronicTopTriplet.rbegin(), hadronicTopTriplet.rend(), sortParticlesByProbability );
    Particle<const OfflineJet> * selectedHadronicTop = &(hadronicTopTriplet.front().second);

    // Determine the hadronic W in the reconstructed hadronic top
    Particle<const OfflineJet> selectedHadronicW( getWfromHadronicTop(*selectedHadronicTop) );

    double chi2ofMasses = pow((selectedHiggs->mass() - referenceHiggsMass_)/15.,2) +
      pow((selectedHadronicTop->mass() - referenceTopMass_)/25.,2) +
      pow((selectedHadronicW.mass() - referenceWmass_)/15.,2);

    // Selected hadronic Top momentum projection along Higgs direction.
    // Note that the object has the operator* overloaded and that this operator requires a pointer to a BaseParticle.
    double scalarProd = (*selectedHadronicTop*selectedHiggs)/sqrt(*selectedHiggs*selectedHiggs);
    double deltaEtaHadronicTopHiggs = selectedHadronicTop->eta() - selectedHiggs->eta();

    // Mass reconstructed with the 5 total jets from Higgs and hadronicTop
    Particle<const OfflineJet> hadronicTopPlusHiggs;
    hadronicTopPlusHiggs.add(selectedHiggs);
    hadronicTopPlusHiggs.add(selectedHadronicTop);

    // Now remove also the hadronic Top associated jets
    jetCollection.erase( find(jetCollection.begin(), jetCollection.end(), (selectedHadronicTop->components())[0]) );
    jetCollection.erase( find(jetCollection.begin(), jetCollection.end(), (selectedHadronicTop->components())[1]) );
    jetCollection.erase( find(jetCollection.begin(), jetCollection.end(), (selectedHadronicTop->components())[2]) );

    // Compute the mass of the remaining jets
    Particle<const OfflineJet> remainingJetsParticle;
    vector<const OfflineJet *>::const_iterator remainingJetsIter = jetCollection.begin();
    for( ; remainingJetsIter != jetCollection.end(); ++remainingJetsIter ) {
      remainingJetsParticle.add(*remainingJetsIter);
    }

    eventVariablesVector_.push_back( selectedHadronicTop->mass() );
    eventVariablesVector_.push_back( selectedHadronicW.mass() );
    eventVariablesVector_.push_back( chi2ofMasses );
    eventVariablesVector_.push_back( scalarProd );
    eventVariablesVector_.push_back( deltaEtaHadronicTopHiggs );
    eventVariablesVector_.push_back( hadronicTopPlusHiggs.mass() );
    eventVariablesVector_.push_back( remainingJetsParticle.mass() );

    // If asked for, fill all the histograms
    if ( fillHistograms_ ) {
      hadronicTopMass_->Fill(selectedHadronicTop->mass());
      hadronicWmass_->Fill( selectedHadronicW.mass() );
      // Evaluate chi2 with Higgs, hadronic Top and W masses
      chi2ofMasses_->Fill( chi2ofMasses );
      hadronicTopProjectionAlongHiggsDirection_->Fill( scalarProd );
      deltaEtaHadronicTopHiggs_->Fill( deltaEtaHadronicTopHiggs );
      hadronicTopPlusHiggsMass_->Fill( hadronicTopPlusHiggs.mass() );
      remainingJetsMass_->Fill( remainingJetsParticle.mass() );
    }

  } // end if jetCollection.size() >= 3

  return eventVariablesVector_;

}

EventVariables::~EventVariables() {

  outputDir_->cd();

  higgsMass_->Write();
  hadronicTopMass_->Write();
  hadronicWmass_->Write();
  chi2ofMasses_->Write();
  firstNjetsMass_[0]->Write();
  firstNjetsMass_[1]->Write();
  firstNjetsCentrality_[0]->Write();
  firstNjetsCentrality_[1]->Write();
  hadronicTopProjectionAlongHiggsDirection_->Write();
  deltaEtaHadronicTopHiggs_->Write();
  goodHt_->Write();
  mEtSig_->Write();
  for( int i=0; i<3; ++i ) deltaPhiMEtNthLeadingJet_[i]->Write();
  hadronicTopPlusHiggsMass_->Write();
  sumHighEffDiscriminantFirst4Jets_->Write();
  sumHighEffDiscriminantFirst6Jets_->Write();
  remainingJetsMass_->Write();
  sixthJetEt_->Write();
  bTagTkInvMass_->Write();

  // delete the multidimensional arrays
  for(unsigned int i=0; i != higgsBinNum_[0]; ++i) {
    for(unsigned int j=0; j != higgsBinNum_[1]; ++j) {
      delete[] trueH_[i][j];
      delete[] falseH_[i][j];
    }
    delete[] trueH_[i];
    delete[] falseH_[i];
  }

  for(unsigned int i=0; i != hadronicTopBinNum_[0]; ++i) {
    for(unsigned int j=0; j != hadronicTopBinNum_[1]; ++j) {
      delete[] trueHadronicTop_[i][j];
      delete[] falseHadronicTop_[i][j];
    }
    delete[] trueHadronicTop_[i];
    delete[] falseHadronicTop_[i];
  }

  for(unsigned int i=0; i != qcdBinNum_[0]; ++i) {
    for(unsigned int j=0; j != qcdBinNum_[1]; ++j) {
      delete[] taggedJet_[i][j];
      delete[] notTaggedJet_[i][j];
    }
    delete[] taggedJet_[i];
    delete[] notTaggedJet_[i];
  }
}

void EventVariables::fillProbabilityMatrices(const string & probabilityFileName, unsigned int * binNum, double * binSize, unsigned int ***& trueArray, unsigned int ***&falseArray ) {

  ifstream probabilityFile(probabilityFileName.c_str());
  if ( !probabilityFile.is_open() ) {
    cout << "file: " << probabilityFileName << "not found or unable to open" << endl;
    exit(1);
  }

  string line;
  // Read the first line with the bin numbers and sizes
  getline(probabilityFile, line);
  stringstream probabilityCounts(line);
  // This loop is based on the text file structure
  string skip;
  int wordCount=0;
  int binNumCount=0;
  int binSizeCount=0;
  while ( !probabilityCounts.eof() ) {
    ++wordCount;
    // Take every three words
    if ( wordCount%3 == 0 ) {
      // Every six words there is the size of the bins, in the other case is the number of bins
      if ( wordCount%6 != 0 ) {
        probabilityCounts >> binNum[binNumCount];
        // cout << "binNum_["<<binNumCount<<"] = " << binNum[binNumCount] << endl;
        ++binNumCount;
      }
      else {
        probabilityCounts >> binSize[binSizeCount];
        // cout << "binSize_["<<binSizeCount<<"] = " << binSize[binSizeCount] << endl;
        ++binSizeCount;
      }
    }
    // Skip the rest of the words
    else {
      probabilityCounts >> skip;
    }
  }

  // Create and fill the arrays
  trueArray = new unsigned int**[binNum[0]];
  falseArray = new unsigned int**[binNum[0]];
  for(unsigned int i=0; i < binNum[0]; ++i) {
    trueArray[i] = new unsigned int*[binNum[1]];
    falseArray[i] = new unsigned int*[binNum[1]];
    for(unsigned int j=0; j < binNum[1]; ++j) {
      trueArray[i][j] = new unsigned int[binNum[2]];
      falseArray[i][j] = new unsigned int[binNum[2]];
      for (unsigned int k=0; k < binNum[2]; ++k) {
        if ( probabilityFile.eof() ) { 
          cout << "ERROR: not enough lines in the file for the required bins" << endl;
          exit(1);
        }

        getline(probabilityFile, line);
        stringstream probabilityCounts(line);
        for ( int w=0; w<6; ++w ) {
          // The second word is for the true case probability
          if( w == 2 ) {
            probabilityCounts >> trueArray[i][j][k];
            // cout << "trueArray["<<i<<"]["<<j<<"]["<<k<<"] = " << trueArray[i][j][k];
            // The fourth word is for the false case probability
          }
          if( w == 4 ) {
            probabilityCounts >> falseArray[i][j][k];
            // cout << " trueArray["<<i<<"]["<<j<<"]["<<k<<"] = " << falseArray[i][j][k] << endl;;
          }
          // Skip the other words
          else probabilityCounts >> skip;
        }
      }

    }
  }
}

Particle<const OfflineJet> EventVariables::firstNjetsParticle( const vector<const OfflineJet *> & jetCollection, const int N ) {
  // Create a particle with the first N jets
  Particle<const OfflineJet> firstNjets;
  int jetCollectionSize = jetCollection.size();
  if ( jetCollectionSize >= N ) {
    int jetsN = 0;
    vector<const OfflineJet *>::const_iterator allJetIt = jetCollection.begin();
    for ( ; allJetIt != jetCollection.end(); ++allJetIt, ++jetsN ) {
      if ( jetsN <= N ) firstNjets.add(*allJetIt);
    }
  }

  unsigned int index = 100;
  if ( N == 6 ) index = 0;
  if ( N == 8 ) index = 1;

  if ( index != 100 ) {
    double firstNjetsMass = firstNjets.mass();
    double firstNjetsCentrality = 0.;
    if ( firstNjetsMass != 0 ) firstNjetsCentrality = firstNjets.e()/firstNjets.mass();
    eventVariablesVector_.push_back(firstNjetsMass);
    eventVariablesVector_.push_back(firstNjetsCentrality);
    if (fillHistograms_) {
      firstNjetsMass_[index]->Fill( firstNjetsMass );
      firstNjetsCentrality_[index]->Fill( firstNjetsCentrality );
    }
  }

  return firstNjets;
}

void EventVariables::allGoodJetsVariables( const vector<const OfflineJet *> & offlineJets, const OfflineMEt * offlineMEt ) {

  double offlineMEtPhi = offlineMEt->phi();

  double goodSumEt = 0.;
  double sumHed4 = 0.;
  double sumHed6 = 0.;
  double bTagTkInvMass = 0.;
  int nthLeadingJet = 0;
  vector<const OfflineJet *>::const_iterator jetPtr = offlineJets.begin();
  for ( ; jetPtr != offlineJets.end(); ++jetPtr, ++nthLeadingJet ) {
    goodSumEt += (*jetPtr)->et();
    // DeltaPhi between the MEt and the second leading jet. Evaluate at most only for the first 3 jets.
    if ( nthLeadingJet < 3 ) {
      double deltaPhiMEtJet = DeltaPhi( offlineMEtPhi, (*jetPtr)->phi() );
      eventVariablesVector_.push_back( deltaPhiMEtJet );
      if (fillHistograms_) deltaPhiMEtNthLeadingJet_[nthLeadingJet]->Fill( deltaPhiMEtJet );
    }

    // Compute the sum of the high efficiency b-tagger discriminant for the 4 and 6 most energetic jets
    double highEffDiscriminant = (*jetPtr)->discriminatorHighEff();
    if ( nthLeadingJet < 4 ) sumHed4 += highEffDiscriminant;
    if ( nthLeadingJet < 6 ) sumHed6 += highEffDiscriminant;
    if ( nthLeadingJet == 5 ) {
      eventVariablesVector_.push_back( (*jetPtr)->et() );
      if (fillHistograms_) sixthJetEt_->Fill((*jetPtr)->et());
    }
    // Medium discriminant, ATTENTION this strongly depends on the b-tagging algorithm
    if ( highEffDiscriminant >  5.3 ) bTagTkInvMass += (*jetPtr)->bTagTkInvMass();
  }
  double goodHt = goodSumEt + offlineMEt->corrL3et();
  eventVariablesVector_.push_back( goodHt );
  eventVariablesVector_.push_back( offlineMEt->corrL3mEtSig() );
  eventVariablesVector_.push_back( sumHed4 );
  eventVariablesVector_.push_back( sumHed6 );
  eventVariablesVector_.push_back( bTagTkInvMass );
  if (fillHistograms_) {
    goodHt_->Fill(goodHt);
    mEtSig_->Fill(offlineMEt->corrL3mEtSig());
    sumHighEffDiscriminantFirst4Jets_->Fill(sumHed4);
    sumHighEffDiscriminantFirst6Jets_->Fill(sumHed6);
    bTagTkInvMass_->Fill( bTagTkInvMass );
  }
}

double EventVariables::evalHiggsPairProbability(const Particle<const OfflineJet> & higgsCandidate) const {

  // Determine bin index
  int etId = int(higgsCandidate.pt()/higgsBinSize_[0]);
  int etaId = int(fabs(higgsCandidate.eta())/higgsBinSize_[1]);
  if ( higgsCandidate.components().size() != 2 ) {
    cout << "Error: higgsCandidate does not have 2 component jets" << endl;
    exit(1);
  }
  const OfflineJet * jet1 = (higgsCandidate.components())[0];
  const OfflineJet * jet2 = (higgsCandidate.components())[1];
  int dRId = int( DeltaR( jet1->eta(), jet1->phi(), jet2->eta(), jet2->phi())/higgsBinSize_[2] );
  if ( etId < 0 || etaId < 0 || dRId < 0 ) cout << "Error: index < 0, will crash..." << endl;
  if ( etId>int(higgsBinNum_[0]-1) ) etId = higgsBinNum_[0]-1;
  if ( etaId>int(higgsBinNum_[1]-1) ) etaId = higgsBinNum_[1]-1;
  if ( dRId>int(higgsBinNum_[2]-1) ) dRId = higgsBinNum_[2]-1;

//   cout << "etId = " << etId << endl;
//   cout << "etaId = " << etaId << endl;
//   cout << "dRId = " << dRId << endl;
//   cout << "trueH_["<<etId<<"]["<<etaId<<"]["<<dRId<<"] = " << trueH_[etId][etaId][dRId] << endl;
//   cout << "falseH_["<<etId<<"]["<<etaId<<"]["<<dRId<<"] = " << falseH_[etId][etaId][dRId] << endl;

  return (double(trueH_[etId][etaId][dRId])/double(falseH_[etId][etaId][dRId]));
}

double EventVariables::evalTopTripletProbability(const Particle<const OfflineJet> & hadronicTopCandidate, const Particle<const OfflineJet> & selectedHiggs) const {

  // Determine bin index
  int etId = int(hadronicTopCandidate.pt()/hadronicTopBinSize_[0]);
  int etaId = int(fabs(hadronicTopCandidate.eta())/hadronicTopBinSize_[1]);
  double dPhiHiggsHadronicTop = DeltaPhi(selectedHiggs.phi(), hadronicTopCandidate.phi());
  int dPhiHiggsHadronicTopId = int(dPhiHiggsHadronicTop/hadronicTopBinSize_[2]);
  if ( hadronicTopCandidate.components().size() != 3 ) {
    cout << "Error: hadronicTopCandidate does not have 3 component jets" << endl;
    exit(1);
  }
  if ( etId < 0 || etaId < 0 || dPhiHiggsHadronicTop < 0 ) cout << "Error: index < 0, will crash..." << endl;
  if ( etId>int(hadronicTopBinNum_[0]-1) ) etId = hadronicTopBinNum_[0]-1;
  if ( etaId>int(hadronicTopBinNum_[1]-1) ) etaId = hadronicTopBinNum_[1]-1;
  if ( dPhiHiggsHadronicTopId>int(hadronicTopBinNum_[2]-1) ) dPhiHiggsHadronicTopId = dPhiHiggsHadronicTopId-1;

//   cout << "etId = " << etId << endl;
//   cout << "etaId = " << etaId << endl;
//   cout << "dPhiHiggsHadronicTopId = " << dPhiHiggsHadronicTopId << endl;
//   cout << "trueHadronicTop_["<<etId<<"]["<<etaId<<"]["<<dPhiHiggsHadronicTopId<<"] = " << trueHadronicTop_[etId][etaId][dPhiHiggsHadronicTopId] << endl;
//   cout << "falseHadronicTop_["<<etId<<"]["<<etaId<<"]["<<dPhiHiggsHadronicTopId<<"] = " << falseHadronicTop_[etId][etaId][dPhiHiggsHadronicTopId] << endl;

  return (double(trueHadronicTop_[etId][etaId][dPhiHiggsHadronicTopId])/double(falseHadronicTop_[etId][etaId][dPhiHiggsHadronicTopId]));
}

Particle<const OfflineJet> EventVariables::getWfromHadronicTop(const Particle<const OfflineJet> & selectedHadronicTop ) const {

  vector<const OfflineJet *> hadronicTopDecayJet( selectedHadronicTop.components() );
  vector<const OfflineJet *>::const_iterator jet = hadronicTopDecayJet.begin();

  Particle<const OfflineJet> selectedW;

  if ( hadronicTopDecayJet.size() != 3 ) return ( selectedW );

  // Create the three possible combinations and take the one with mass closest to the W mass
  vector<Particle<const OfflineJet> > hadronicWcandidates;
  hadronicWcandidates.push_back( Particle<const OfflineJet>( hadronicTopDecayJet[0]) );
  hadronicWcandidates[0].add( hadronicTopDecayJet[1] );
  // cout << "mass[0] = " << hadronicWcandidates[0].mass() << endl;
  hadronicWcandidates.push_back( Particle<const OfflineJet>( hadronicTopDecayJet[0]) );
  hadronicWcandidates[1].add( hadronicTopDecayJet[2] );
  // cout << "mass[1] = " << hadronicWcandidates[1].mass() << endl;
  hadronicWcandidates.push_back( Particle<const OfflineJet>( hadronicTopDecayJet[1]) );
  hadronicWcandidates[2].add( hadronicTopDecayJet[2] );
  // cout << "mass[2] = " << hadronicWcandidates[2].mass() << endl;

  // SortByMassDifference<Particle<const OfflineJet> > sortFunctor(referenceWmass_);
  sort( hadronicWcandidates.begin(), hadronicWcandidates.end(), SortByMassDifference<Particle<const OfflineJet> >(referenceWmass_) );

  // cout << "selected W mass = " << hadronicWcandidates[0].mass() << endl;

  return hadronicWcandidates.front();
}

#endif // EVENTVARIABLES_CC
