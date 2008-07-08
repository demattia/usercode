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
}

EventVariables::EventVariables( const string & higgsFileName, const string & hadronicTopFileName, const string & qcdFileName, TString suffix, TFile * outputFile ) {

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

  outputDir_ = outputFile_->mkdir("EventVariables");
  outputDir_->cd();

  higgsMass_       = new TH1D( "higgsMass" + suffix, "reconstructed higgs mass" + suffix, 50, 0, 300 );
  hadronicTopMass_ = new TH1D( "hadronicTopMass" + suffix, "reconstructed hadronic top mass" + suffix, 50, 0, 800 );
}

void EventVariables::fill( vector<const OfflineJet *> jetCollection, const vector<const OfflineJet *> & bTaggedJetCollection ) {

  // First of all sort the collection in Et. The first is the most energetic
  // Do this here, not inside the method to do it only once.
  sort( jetCollection.rbegin(), jetCollection.rend(), sortOfflineJetPtr );

  Particle<const OfflineJet> first6jets( firstNjetsParticle( jetCollection, 6 ) );
  Particle<const OfflineJet> first8jets( firstNjetsParticle( jetCollection, 8 ) );


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

  higgsMass_->Fill(selectedHiggs->mass());

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
    hadronicTopMass_->Fill((hadronicTopTriplet.front().second).mass());
  } // end if jetCollection.size() >= 3

}

EventVariables::~EventVariables() {

  outputDir_->cd();

  higgsMass_->Write();
  hadronicTopMass_->Write();

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

Particle<const OfflineJet> EventVariables::firstNjetsParticle( const vector<const OfflineJet *> & jetCollection, const int N ) const {
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

  double firstNjetsMass = firstNjets.mass();
  double firstNjetsCentrality = 0.;
  if ( firstNjetsMass != 0 ) firstNjetsCentrality = firstNjets.e()/firstNjets.mass();
  cout << "mass of the first "<<N<< " jets = " << firstNjetsMass << endl;
  cout << "centrality of the first "<<N<<" jets = " << firstNjetsCentrality << endl;

  return firstNjets;
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

#endif // EVENTVARIABLES_CC
