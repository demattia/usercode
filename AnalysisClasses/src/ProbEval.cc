#ifndef PROBEVAL_CC
#define PROBEVAL_CC

#include "AnalysisExamples/AnalysisClasses/interface/ProbEval.h"

using namespace std;
using namespace anaobj;

void ProbEval::eval_(int iJet, double tempProb) {
  // Condition to check if the collection is finished
  //    cout<<"iJet: "<<iJet<<", jetNum_: "<<jetNum_<<endl;
  if ( iJet < jetNum_ ) {
    if(iJet >= jetNum_ - jetNumCut_){
      //      cout << "if: iJet = " << iJet <<endl;
      // Loop on the two possible configurations of a jet:
      // - associated to the closest vertex
      // - associated to the second closest vertex
      for ( int conf=0; conf<=1; ++conf ) {
	// call the function on the next jet passing the product of the
	// probability for this configuration with that of all the previous
	// jets.
	// Store also the index of the associated vertex in an array.
	// This is used in the esle part to evaluate the probability of the vertex configuration.
	// The values of this array are overwritten at each step, so the final configuration is
	// always updated correctly.
	if(conf == 0){
	  // cout<<"if conf == "<<conf << ", tempProb = " << tempProb + TMath::Log(jetArray_[iJet].probBestVtx()) <<endl;
	  jetVtxId_[iJet] = jetArray_[iJet].idBestVtx();
	  eval_( iJet+1, tempProb + TMath::Log(jetArray_[iJet].probBestVtx() ) );
	}else{
	  // cout<<"if conf == "<<conf<< ", tempProb = " << tempProb + TMath::Log(jetArray_[iJet].probSecondBestVtx()) <<endl;
	  jetVtxId_[iJet] = jetArray_[iJet].idSecondBestVtx();
	  eval_( iJet+1, tempProb + TMath::Log(jetArray_[iJet].probSecondBestVtx() ) );
	}
      }
    }
    //fixed to most probable vertex
    else{
      jetVtxId_[iJet] = jetArray_[iJet].idBestVtx();
      eval_( iJet+1, tempProb + TMath::Log(jetArray_[iJet].probBestVtx() ) );
    }
  }
  // if it is after the last jet, fill the vector with the probability
  else {
    //     cout << "else: iJet = " << iJet <<endl;
    // first evaluate the probability for the multiplicities of the verteces
    // the index in the vector is in biunivocal correspondence with the
    // configuration.
    // Use this information to extract the configuration and loop on the
    // jet collection to associate the jets to the verteces.
    // Loop then on the vertex collection and evalute the product of the
    // probabilities associated with that configuration.
    // Call this product probMultiplicity, then:

    //Returns the sum of the logs of probabilities
    double probMultiplicity = probCal_.computeEventProb( vtxCounter_.countVtx( vtxCounter_.countJet(jetVtxId_, jetNum_) ) );
    double finalProb = tempProb + probMultiplicity;

    //    cout << "finalProb["<<index_<<"] = " << tempProb << " + " << probMultiplicity << " = " << finalProb << endl;

    if( finalProb > maxProb_ ){
      maxProb_ = finalProb;
      maxProbId_ = index_;
    }

    //      probVec_.push_back(tempProb + probMultiplicity);
    //      cout << "Probability["<<probVec_.size()<<"] = " << tempProb << " + "  <<probMultiplicity<< " = " <<tempProb + probMultiplicity << endl;
    ++index_;
  }
}

map<int, int> ProbEval::evalProb(SimpleCaloJetCollection vec_CaloJet, int jetNumCut) {

  jetNumCut_ = jetNumCut;

  map<int, int> jetVtxMap;

  if ( !(vec_CaloJet.empty()) ) {
    // Store the jet collection in an array for faster random access
    jetNum_ = vec_CaloJet.size();
    jetArray_ = new SimpleCaloJet[jetNum_];
    jetVtxId_ = new int[jetNum_];
    int tempCounter = 0;
    //sorting in decreasing order with probability ratio
    sort(vec_CaloJet.rbegin(),vec_CaloJet.rend(),simpleCaloJetCollSort);
    SimpleCaloJetCollection::const_iterator  vec_CaloJetIt = vec_CaloJet.begin();
    for(; vec_CaloJetIt != vec_CaloJet.end() ; vec_CaloJetIt++, ++tempCounter){
      //      cout<<"p1/p2["<<tempCounter<<"]= "<<vec_CaloJetIt->probBestVtx()/vec_CaloJetIt->probSecondBestVtx()<<endl;
      jetArray_[tempCounter] = *vec_CaloJetIt;
      jetVtxId_[tempCounter] = 0;
    }
    //     cout<<"before eval_ (0.0)"<<endl;
    // call eval_ passing the index to the first jet and 1 as the starting
    // probability.
    eval_(0,0);
    // Possibily perform the association for the most probable configuration
    // and return it (verteces with number of jets associated to them).

    //      cout<<"probVec_.size(): "<<probVec_.size()<<endl;    

    
    SimpleCaloJetCollection::reverse_iterator  rItCaloJet = vec_CaloJet.rbegin();
    int tempId = maxProbId_;
    //    if(maxProbId_ == 0) cout<<"All best"<<endl;
    //    cout<<"maxProb: "<<maxProb_<<endl;
    for(; rItCaloJet != vec_CaloJet.rend(); ++rItCaloJet ){
      //If it is even it is the best vtx
      if(tempId%2 == 0) {
	jetVtxMap.insert(make_pair(rItCaloJet->idJet(), rItCaloJet->idBestVtx()));
      }else{
	jetVtxMap.insert(make_pair(rItCaloJet->idJet(), rItCaloJet->idSecondBestVtx()));
      }
      tempId = tempId/2;
    }

    delete[] jetArray_;
    delete[] jetVtxId_;
    maxProb_ = -10000;
    maxProbId_ = -1;
    index_ = 0;
    jetNumCut_ = 0;
  }
  return jetVtxMap;
}

#endif // PROBEVAL_CC
