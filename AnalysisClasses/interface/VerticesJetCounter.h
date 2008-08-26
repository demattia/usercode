#ifndef VERTICESJETCOUNTER_H
#define VERTICESJETCOUNTER_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "AnalysisExamples/AnalysisObjects/interface/AssociationSummary.h"


using namespace std;  
using namespace anaobj;

/**
 *
 * class VerticesJetCounter
 *
 * Used to count the # of jet assigned to every vertices
 * 
 * Using class AssociationSummary
 * 
 *
 *
 * Author R.Casagrande - 06/05/2008
 *
 */


class VerticesJetCounter {
 private:
  
 public:
  VerticesJetCounter (){
    
  }
  
  std::map<int,int> countJet(anaobj::AssociationSummaryCollection & jetMap_ ){
    
    // mappa di numero jet associati a ogni vertice
    std::map<int,int> nJetMap;

    //    if( jetMap_.empty() ) cout<<"jetMap passed to VerticesJetCounter.countJet() is empty!"<<endl;
    
    if( !jetMap_.empty() ){
      // conto il numero di jet associati a ogni vertice
      anaobj::AssociationSummaryCollection::const_iterator  jetMap_It = jetMap_.begin();
      for(;jetMap_It != jetMap_.end(); ++jetMap_It ){
	int idVtx1 = jetMap_It->idVtx();
	bool assignedJet = jetMap_It->Associated();
	
	//se il jet e' stato assegnato conto il numero di jet assegnati a ciascun vertice
	if (assignedJet){
	  pair<map<int, int>::iterator, bool> mapInsert(nJetMap.insert( make_pair( idVtx1 , 1 ) ) );
	  
	  //se il vertice non e' stato inserito c'era gia', allora incremento il contatore corrispondente
	  if ( !(mapInsert.second) ) {
	    ++((mapInsert.first)->second);
	  }
	}  
      }//end loop
    }
    //    if(nJetMap.empty()) cout<<"njetMap returned from VerticesJetCounter.countJet(1) is empty!"<<endl;
    return nJetMap;
  }
  // Il metodo prende una mappa <id jet , id vertice> e restituisce una mappa <id vertice, # jet associati>
  std::map<int,int> countJet(std::map<int,int> & idRJAlgMap){
    
    std::map<int,int> nJetMap;

    //    if( idRJAlgMap.empty() ) cout<<"idRJAlgMap passed to VerticesJetCounter.countJet() is empty!"<<endl;

    if( !idRJAlgMap.empty() ){
      // conto il numero di jet associati a ogni vertice
      std::map<int,int>::const_iterator  idRJAlgMapIt = idRJAlgMap.begin();
      for(;idRJAlgMapIt != idRJAlgMap.end(); ++idRJAlgMapIt ){
	int idVtx1 = idRJAlgMapIt->second;
	
	pair<map<int, int>::iterator, bool> mapInsert(nJetMap.insert( make_pair( idVtx1 , 1 ) ) );
	
	//se il vertice non e' stato inserito c'era gia', allora incremento il contatore corrispondente
	if ( !(mapInsert.second) ) {
	  ++((mapInsert.first)->second);
	}
      }//end loop
    }
    //    if(nJetMap.empty()) cout<<"njetMap returned from VerticesJetCounter.countJet(2) is empty!"<<endl;
    return nJetMap;
  }

 // Il metodo prende una mappa <id jet , id vertice> e restituisce una mappa <id vertice, # jet associati>
  std::map<int,int> countJet(int* idRJAlgMap, int size){
  
    std::map<int,int> nJetMap;

    if( size != 0 ){
      // conto il numero di jet associati a ogni vertice

      for( int i = 0 ;i < size ; ++i ){
	int idVtx1 = idRJAlgMap[i];
	
	pair<map<int, int>::iterator, bool> mapInsert(nJetMap.insert( make_pair( idVtx1 , 1 ) ) );
	
	//se il vertice non e' stato inserito c'era gia', allora incremento il contatore corrispondente
	if ( !(mapInsert.second) ) {
	  ++((mapInsert.first)->second);
	}
      }//end loop
    }
    return nJetMap;
  }
  
  // Il metodo prende una mappa <id vertice,# jet associati> e restituisce una mappa < # jet associati, # vertici con quel # jet>
  std::map<int,int> countVtx(const std::map<int,int> & mapNjet_){

    std::map<int,int> nVtxMap;

    //    if( mapNjet_.empty()) cout<<"mapNjet passed to VerticesJetCounter.countVtx() is empty!"<<endl;
 
    if( !mapNjet_.empty()){
      std::map<int,int>::const_iterator mapNjet_It = mapNjet_.begin();
      for(;mapNjet_It != mapNjet_.end(); ++mapNjet_It){
	
	int nJet = mapNjet_It->second;
	
	pair<map<int, int>::iterator, bool> mapInsert(nVtxMap.insert( make_pair( nJet , 1 ) ) );
	if ( !(mapInsert.second) ) {
	  ++((mapInsert.first)->second);
	}
      }//end loop
    }
    //    if(nVtxMap.empty()) cout<<"nVtxMap returned from VerticesJetCounter.countVtx() is empty!"<<endl;
    return nVtxMap;
  }
  
};

   
#endif
   
   
