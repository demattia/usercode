#ifndef PROBCAL_H
#define PROBCAL_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include "TMath.h"

using namespace std;  
using namespace anaobj;

/**
 *
 * class ProbCal
 * 
 * 
 * Return log of combinated probability for the event configuration
 * 
 *
 * Author R.Casagrande - 07/05/2008
 *
 */


class ProbCal {
 private:
  
  int nJetMax_;
  int nCollMax_;
  double pSing25_[9][10]; 
  
 public:

  ProbCal( const char * matrixFile = "MatrixProducerSing25.asc" ) {
    vector<double> pMatrixSing25;
    //--------------------------------------------------------------------------
    // Lettura da file .asc
    //--------------------------------------------------------------------------
    //matrice probabilita' singole
    //read from file
    //  ifstream b("MatrixProducerSing50.asc");
    ifstream b(matrixFile);
    string line1;
    while (b) {
      getline(b,line1);
      stringstream events;
      double number = 0.;
      if (line1 != "") {
	events << line1;
	events >> number;
      	pMatrixSing25.push_back(number);
      }
    }
    
    nJetMax_ = 9;
    nCollMax_ = 10;
    for(int nJet = 0; nJet < nJetMax_; nJet++){
      for(int nColl = 0; nColl < nCollMax_; nColl++){
	pSing25_[nJet][nColl] = 0.;
      }
    }

    //    if(pMatrixSing25.empty()) cout<<"Probability matrix in ProbCal is empty!"<<endl;
    
    if(!pMatrixSing25.empty()){
      //importo la matrice e la faccio diventare da vettore un array
      std::vector<double>::const_iterator pMatrixIt = pMatrixSing25.begin();
      for(int nJet = 0; nJet < nJetMax_; nJet++){
	for(int nColl = 0; nColl < nCollMax_; nColl++, pMatrixIt++){
	  if(pMatrixIt != pMatrixSing25.end()){
	    pSing25_[nJet][nColl] = *pMatrixIt;
	  }else{
	    //pongo le probabilita' superiori a quelle disponibili nell'ascfile pari all'ultima disponibile
	    //	    cout<<"pSing25 out of range ( ProbCal constructor )"<<endl;
	    pSing25_[nJet][nColl] = *pMatrixSing25.rbegin();
	  }
	}
      }
    }
    //    ProbCal( pMatrixSing25 );
  }

  // ProbCal ( pMatrixSing25 ){
  //  ProbCal (  vector<double> & pMatrix_  ){ 
  //   nJetMax_ = 9;
  //  nCollMax_ = 10;
  //  eventProbability_ = 0.;
  //  for(int nJet = 0; nJet < nJetMax_; nJet++){
  //   for(int nColl = 0; nColl < nCollMax_; nColl++){
  //	pSing25_[nJet][nColl] = 0.;
  //    }
  //  }

  //   if(pMatrix_.empty()) cout<<"Probability matrix in ProbCal is empty!"<<endl;

  //  if(!pMatrix_.empty()){
  //  //importo la matrice e la faccio diventare da vettore un array
  //     std::vector<double>::const_iterator pMatrixIt = pMatrix_.begin();
  //    for(int nJet = 0; nJet < nJetMax_; nJet++){
  //	for(int nColl = 0; nColl < nCollMax_; nColl++, pMatrixIt++){
  //	  if(pMatrixIt != pMatrix_.end()){
  //	    pSing25_[nJet][nColl] = *pMatrixIt;
  //	  }else{
  //	    //pongo le probabilita' superiori a quelle disponibili nell'ascfile pari all'ultima disponibile
  //	    cout<<"pSing25 out of range ( ProbCal constructor )"<<endl;
  //	    pSing25_[nJet][nColl] = *pMatrix_.rbegin();
  //	  }
  //	}
  //      }
  //    }


  //  }
  
  double computeEventProb(std::map<int,int>  nVtxMap){
  
    double eventProbability = 0.;

    //    if(nVtxMap.empty()) cout<<"Map passed to computeEventProb empty!"<<endl;
    
    if(!nVtxMap.empty()){

      //njet max 8
      //loop fino a njet 20, perche' siano consistenti le probabilita' anche se le mappe sono di size diverse 
      //20 e' un valore ragionevolemente alto perche' non si abbiano vertici con un numero superiore di jet
      //se considero la probabilita' superiore a 8 jet come uguale a quella di averne 8 potrei arrivare a 20
      //e per considerare anche le probabilita' di aver 0 vertici con un dato numero di jet
      //nColl calcolato fino a 9, possibilita' aumentare ricalcolando la matrice con MatrixProducer
      ////loop nJet fino a 8 jet, perche' non mi e' possibile calcolare oltre  -> calcolo matrice da NOPU sample

      for(int nVtxMapIt = 0 ; nVtxMapIt < 9 ; nVtxMapIt++ ){
	
	int nColl = nVtxMap[nVtxMapIt];
	if(pSing25_[nVtxMapIt][nColl] != 0.){
	  //    cout<<nColl<<endl;
	  if(nColl < nCollMax_){
	    //	cout<<"in loop"<<endl;
	    eventProbability +=  TMath::Log(pSing25_[nVtxMapIt][nColl]);
	  }
	  if(nColl >= nCollMax_){
	    //	    cout<<"Warning: nColl out of range"<<endl;
	    eventProbability +=  TMath::Log( pSing25_[nVtxMapIt][nCollMax_ - 1] );
	  }
	}
	//	cout<<eventProbability<<endl;
      }
    }
  
    return eventProbability;
  }
  
};


#endif
   
   
