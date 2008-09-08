#ifndef JETVERTEXASSOCIATOR_CC
#define JETVERTEXASSOCIATOR_CC

#include "AnalysisExamples/AnalysisClasses/interface/JetVertexAssociator.h"

JetVertexAssociator::JetVertexAssociator( double jetEtCut, double jetEtaCut, int sigmaCut, double sigmaVtx ) {
  jetEtCut_ = jetEtCut;
  jetEtaCut_ = jetEtaCut;
  sigmaCut_ = sigmaCut;
  sigmaVtx_ = sigmaVtx;

  //--------------------------------------------------------------------------
  // Lettura da file .asc
  //--------------------------------------------------------------------------
  //matrice probabilita' singole
  //read from file
  ifstream b("MatrixProducerSing25.asc");
  string line1;
  while (b) {
    getline(b,line1);
    stringstream events;
    double number = 0.;
    if (line1 != "") {
      events << line1;
      events >> number;
      pMatrixSing25_.push_back(number);
    }
  }
  if(pMatrixSing25_.size() == 0) cout<<"WARNING: no probability matrix loaded!"<<endl;
  //matrice probabilita' combinate
  //read from file
  ifstream a("MatrixProducer25.asc"); 
  string line;
  while (a) {
    getline(a,line);
    stringstream events;
    double number = 0.;
    if (line != "") {
      events << line;
      events >> number;
      pMatrix25_.push_back(number);
    }
  }
   if(pMatrix25_.size() == 0) cout<<"WARNING: no probability matrix loaded!"<<endl;
}

OfflineJetCollection JetVertexAssociator::associate( const OfflineJetCollection & caloJets, const BaseVertexCollection & recoVertexes ) {

  OfflineJetCollection vtxAssocJet;
  vector<const OfflineJet *> vec_JetWithNoTracks;

  //--------------------------------
  //Reco tracks and Offline calojet
  //--------------------------------

  //Et, eta, phi Offline caloJet
  double Et_Jet = 0.;
  double eta_Jet = 0.;
  double phi_Jet = 0.;

  SimpleCaloJetCollection vec_CaloJet; //Collection of simpleCalojet passing al cuts

  int Njet = 0 ; //# jet

  for ( OfflineJetCollection::const_iterator allJetIt = caloJets.begin(); allJetIt != caloJets.end(); ++allJetIt ) {

    //vector of reco selected tracks' Pt, eta, phi, Z
    SimpleTrackCollection vec_RecoTk;
	  
    //Et, eta, phi Offline caloJet
    Et_Jet = allJetIt->et();
    eta_Jet = allJetIt->eta();
    phi_Jet =  allJetIt->phi();
    
    //vettore di reference alle selected track dell'offline jet
    SimpleTrackRefVector assocTk( allJetIt->tkRefVec() ) ;
	  
    SimpleTrackRefVector::const_iterator TkColl_it = assocTk.begin();
    for (; TkColl_it != assocTk.end(); ++TkColl_it){ 
	    
      //2D Significance
      double Tk_IpS2D = 0.;
      Tk_IpS2D =(*TkColl_it)->ip2Dsignificance();

      //Reco Track Z
      double Tk_dz = 0.;
      Tk_dz = (*TkColl_it)->z();
	    
      //Reco Track Z error
      double Tk_dz_error = 0.;
      Tk_dz_error= (*TkColl_it)->zError();
	    
      //Reco Track eta
      double Tk_eta = 0.;
      Tk_eta = (*TkColl_it)->eta();
	    
      //Reco Track phi
      double Tk_phi = 0.;
      Tk_phi = (*TkColl_it)->phi(); 
	    
      //Reco Track Pt
      double Tk_Pt = 0.;
      Tk_Pt = (*TkColl_it)->pt();
    
      //cuts in filling reco SimpleTrackCollection
      //-----------------------------------

      if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_  && fabs(Tk_IpS2D) < 3.) {
	vec_RecoTk.push_back(SimpleTrack(Tk_Pt,
					 Tk_eta,
					 Tk_phi,
					 Tk_dz, 
					 Tk_dz_error,
					 Tk_IpS2D)
			     );
	
      }//end cut Et(jet)>25, IP significance < 3. , eta (jet) < 3.   
    }// end loop on selected tracks
	  
    //----------------------------------
    //Taglio Et(jet)>25 , eta(jet) < 3
    //----------------------------------
    if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_ ) {

      //----------------------------------------------
      //Algoritmo di reiezione tracce ricostruite nel jet
      //Z jet calcolata con la media pesata, scartando 
      //le tracce piu' distanti di sigmaCut * sigma dal valore 
      //della media pesata con algoritmo iterativo
      //----------------------------------------------

      //chiamo l'algoritmo di reiezione (e' una classe con il metodo eval(# sigma taglio)) sulla collezione di tracce ricostruite
      //agisce se il numero di trace nel jet e' > 2
      RejectionAlg tkRejection(vec_RecoTk);
      SimpleTrackCollection recoTkClosest(tkRejection.eval(sigmaCut_));

      //chiamo la funzione che fa la media pesata wAvenger, che restituisce una pair(# tracce, pair(value, error))
      pair< int, pair< double, double > > zJetRejectionAlg( wAverager(recoTkClosest));
      //z jet e suo errore
      double zJet = zJetRejectionAlg.second.first;
      double zJetError = zJetRejectionAlg.second.second;

      int NTksRJAlg_S3 = zJetRejectionAlg.first ; //contatore per il numero di tracce con significanza minore di 3	
      
      //sono nello stesso numero di dei jet prima della reiezione, perche' quelli con 2 o meno tracce skippano l'algoritmo
      //e gli altri hanno tracce comunque, cambiera' il contenuto in tracce dei jet

      //---------------------------------
      // Delta Z minimo tra jet e vertici 
      // con Z media pesata tracce
      //---------------------------------
	
      //Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con la media pesata dopo la reiezione e i vertici
      //la prima coppia contiene l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
      //la funzione vuole in input una collezione di vertici e un double
      //----------------------------------------------------------------------------------------------------------------------
      //if # of selected tracks with IP Significance < 3. != 0
      if (NTksRJAlg_S3 != 0){
	pair<pair<int, double>, pair<int, double> > dzAfterRejection;

	dzAfterRejection = dzAssociator(recoVertexes,zJetRejectionAlg.second.first ) ;

	int  idRJAlg = dzAfterRejection.first.first; //id best vtx
	int  id2RJAlg = dzAfterRejection.second.first; //id second best vtx
	double dzRJAlg = dzAfterRejection.first.second; //minimum dz jet vtx
	double dz2RJAlg = dzAfterRejection.second.second; //second best minimum dz jet vtx

	//------------------
	// Calcolo P1/P2
	//------------------
	//Uso invece dell'errore sui vertici dato dalla collezione, l'errore sulla DeltaZ jet vertice dopo l'algoritmo di 
	//reiezione calcolato per il campione QCD170-230 = 0.03556
	
	//Ho ricavato l'errore sul vertice dalla "deconvoluzione" con l'errore medio sui jet, 
	//ora lo combino con l'errore per ogni jet per ottenere l'errore sulla Dz del jet
	
	double dzJetError = sqrt(pow(zJetError,2)+ pow(sigmaVtx_ ,2));
	double dz2JetError = sqrt(pow(zJetError,2)+ pow(sigmaVtx_,2));
	
	//calcolo la probabilita' per il primo e per il secondo vertice
	double p1 = TMath::Erfc(fabs(dzRJAlg/dzJetError));
	double p2 = TMath::Erfc(fabs(dz2RJAlg/dz2JetError));
	
	//se e' inferiore a 10e-30 per limite del calcolatore la pongo al valore minimo 10e-30, per evitare eccezioni nell'algoritmo
	if(p1 < 10e-30) p1 = 10e-30;
	if(p2 < 10e-30) p2 = 10e-30;

	//riempio una collezione di simpleCaloJet che passano le richieste su Et ed Eta del jet e hanno tracce 
	vec_CaloJet.push_back(SimpleCaloJet( Et_Jet,
					     eta_Jet,
					     phi_Jet,
					     zJet,
					     zJetError,
					     Njet,
					     idRJAlg,
					     id2RJAlg,
					     dzRJAlg,
					     dz2RJAlg,
					     p1,
					     p2,
					     &(*allJetIt) ) );
	Njet++;
      }//end cut NTksRJAlg_S3 != 0
      else{
	//fill a vector of offlineJet with jet with non tracks
	vec_JetWithNoTracks.push_back(&(*allJetIt));

      }

    }// end cut Et(jet>25) eta<3
  }//end loop calojet

  //if jetCollection after cuts !empty 
  if( !vec_CaloJet.empty()){
    //----------------------
    //PROBEVAL
    //----------------------

    ProbEval ValueProb_;
    //mappa di associazione (idJet, idVtx)
    map<int,int> mapProbEval(ValueProb_.evalProb(vec_CaloJet));

    VerticesJetCounter JetCountProbEval;
    //Overloaded method conto quanti jet ha ogni vertice e salvo una mappa con idVtx e # jet
    std::map<int,int> nJetMapProbEval(JetCountProbEval.countJet(mapProbEval));

    //cerco il vertice con max # jet assegnati	
    std::map<int , int> VtxMaxProbEval; //mappa con key il # di jet e value l'id del vertice

    //riempimento mappa con key il # di jet e value l'id del vertice
    std::map<int,int>::const_iterator nJetMapProbEvalIt = nJetMapProbEval.begin();
    for(; nJetMapProbEvalIt != nJetMapProbEval.end() ; nJetMapProbEvalIt++){
      VtxMaxProbEval.insert( make_pair( nJetMapProbEvalIt->second, nJetMapProbEvalIt->first ) );
    }
      
    //id vertice con max num jet e sua Z (la mappa e' ordinata in valori crescenti 
    //della key che per queste mappe e' scelta esser il # di jet)
    int  idPVProbEval = VtxMaxProbEval.rbegin()->second;
 
    //loop on mapProbEval (idJet , idVtx) to select the PV jet
    std::map<int,int>::const_iterator mapProbEvalIt = mapProbEval.begin();
    for(; mapProbEvalIt != mapProbEval.end() ; ++mapProbEvalIt ){
      //if the id vtx is the id of the primary vertex push_back a vector of jet
      if (idPVProbEval == (*mapProbEvalIt).second) {
	int idJet =(*mapProbEvalIt).first;
	vtxAssocJet.push_back(*((vec_CaloJet[idJet]).offlineJetPtr()));
      }
    }
  } //end if !caloJet.empty()  

  //add to the PV jet collection jets with no tracks
  vector<const OfflineJet *>::const_iterator vec_JetWithNoTracksIt = vec_JetWithNoTracks.begin();
  for(; vec_JetWithNoTracksIt != vec_JetWithNoTracks.end(); ++vec_JetWithNoTracksIt){
    vtxAssocJet.push_back(**vec_JetWithNoTracksIt);
  }

  return vtxAssocJet;

}



#endif
