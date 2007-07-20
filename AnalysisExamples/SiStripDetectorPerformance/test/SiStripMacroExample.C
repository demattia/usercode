/**

FWLite Macro Example to retrieve from the event plots on following info:
Trigger Bits
Tracks
Clusters

HowTo Run


root [0] .L SiStripMacroExample.C
root [1] SiStripMacroExample(<runNumber>,<root_file>)

eg: under lxcmsf1
root [1] SiStripMacroExample(2501,"rfio:/data/noeding/pass2_subset/2501/reco_full_2501.root")


**/

gROOT->Reset();

int DEBUG=0;

TRFIOFile *file;
TFile *outputfile;
TTree *events;

TCanvas myCanvas("c","c");

TObjArray Hlist(0);

std::vector<reco::Track> vTracks;
TBranch* TracksB;

std::vector<edm::DetSet<SiStripClusterInfo> > dsvSiStripClusterInfo;
TBranch* SiStripClusterInfoB;

std::vector<LTCDigi> ltcdigis;
TBranch *LTCDigiB; 

int totTrigs[6];

TString outfilename;

//-------------------------------

void book_histos(){
  Hlist.Add(new TH1F("TriggerBits","TriggerBits",6,-0.5,5.5));
  Hlist.Add(new TH1F("nTracks","nTracks",10,-0.5,9.5));
  Hlist.Add(new TH1F("nRecHits","nRecHits",10,-0.5,9.5));
  Hlist.Add(new TH1F("cCharge","cCharge",50,-0.5,200.5));
  Hlist.Add(new TH1F("cNoise","cNoise",16,-0.5,15.5));
  Hlist.Add(new TH1F("cStoN","cStoN",12,-0.5,60.5));
  Hlist.Add(new TH2F("cEta","cEta",5,-2.5,2.5,5,-2.5,2.5));
  Hlist.Add(new TH1F("cWidth","cWidth",8,-0.5,7.5));
}

//-------------------------------

void save_histos(TString outfilename){
  
  TString filename=outfilename+TString(".root");
  std::cout << "........ Opening file "<< filename << std::endl;
  outputfile = new TFile(filename,"RECREATE");
  Hlist.Write();
  outputfile->Close();  
  std::cout << "........ Closed"<< std::endl;

  filename=outfilename+TString(".ps");
  std::cout << "........ Opening file "<< filename << std::endl;
  TPostScript ps(filename);

  for (int ih=0; ih<Hlist.GetEntries();ih++){
    Hlist[ih]->Draw();
    myCanvas.Update();
    ps.NewPage();
  }
  
  ps.Close();
  std::cout << "........ Closed"<< std::endl;
  std::cout << "to see the file please do \n gv " << filename << std::endl;

}

//-------------------------------

void list_of_branches(){
  
  std::cout << "\nList of branches in the Event tree" << std::endl;
  TObjArray * branches = Events->GetListOfBranches();
  int n = branches->GetEntries();
  for( int i = 0; i < n; ++i ) {
    TBranch * branch = (TBranch*) branches->At( i );
    
    std::cout << ">>> " << branch->GetName() << std::endl;
  }
}

//-------------------------------

void set_branches(){
  events = (TTree*)file->Get("Events");
  //Needed for SetAddress to work right
  events->GetEntry();

  events->SetBranchAddress("recoTracks_cosmictrackfinder__SISTRIPRECODIGISCLUSTERS.obj",&vTracks,&TracksB);

  events->SetBranchAddress("SiStripClusterInfoedmDetSetVector_siStripClusterInfoProducer__SISTRIPRECODIGISCLUSTERS.obj._sets",&dsvSiStripClusterInfo,&SiStripClusterInfoB);

  events->SetBranchAddress("LTCDigis_ltcraw2digi__FU.obj",&ltcdigis,&LTCDigiB);
}

//-------------------------------

void cluster_study(int idetset){
  std::vector<SiStripClusterInfo> vSiStripClusterInfo = dsvSiStripClusterInfo[idetset].data;

  for (int icluster=0; icluster!=vSiStripClusterInfo.size(); icluster++){

    SiStripClusterInfo siStripClusterInfo=vSiStripClusterInfo[icluster];

    ((TH1F*) Hlist.FindObject("cCharge"))->Fill(siStripClusterInfo.charge());
    ((TH1F*) Hlist.FindObject("cNoise"))->Fill(siStripClusterInfo.noise());
    if (siStripClusterInfo.noise())
      ((TH1F*) Hlist.FindObject("cStoN"))->Fill(siStripClusterInfo.charge()/siStripClusterInfo.noise());
    ((TH1F*) Hlist.FindObject("cWidth"))->Fill(siStripClusterInfo.width());
    ((TH2F*) Hlist.FindObject("cEta"))->Fill(
					     (siStripClusterInfo.chargeL()-siStripClusterInfo.chargeR())/siStripClusterInfo.charge(),
					     (siStripClusterInfo.chargeL()+siStripClusterInfo.chargeR())/siStripClusterInfo.charge()
					     );

    if (DEBUG>2){
      std::cout << "\t\tcluster first strip " << siStripClusterInfo.firstStrip() << std::endl;
      std::cout << "\t\tcluster charge " << siStripClusterInfo.charge()    << std::endl;
      std::cout << "\t\tcluster noise " <<   siStripClusterInfo.noise()     << std::endl;
      std::cout << "\t\tcluster position " <<     siStripClusterInfo.position()  << std::endl;
      std::cout << "\t\tcluster width " <<    siStripClusterInfo.width()     << std::endl;
      std::cout << "\t\tcluster maxCharge " <<    siStripClusterInfo.maxCharge() << std::endl;
      std::cout << "\t\tcluster maxPos " <<   siStripClusterInfo.maxPos()  << std::endl;  
      std::cout << "\t\tcluster chargeL " <<   siStripClusterInfo.chargeL() << std::endl;  
      std::cout << "\t\tcluster chargeR " <<    siStripClusterInfo.chargeR() << std::endl;  

      std::cout << "\t\tcluster detid " <<   siStripClusterInfo.geographicalId() << std::endl;  

      const std::vector<uint16_t>&  stripAmplitudes_ = siStripClusterInfo.stripAmplitudes();
      for (int i=0;i<stripAmplitudes_.size();i++){
	std::cout << "\t\t\t strip amplitudes " << i << " " << stripAmplitudes_[i]<< std::endl;  
      }
      const std::vector<uint16_t>&  rawdigiAmplitudesL_ = siStripClusterInfo.rawdigiAmplitudesL();
      for (int i=0;i<rawdigiAmplitudesL_.size();i++){
	std::cout << "\t\t\t rawdigi amplitudes L " << i << " " << rawdigiAmplitudesL_[i]<< std::endl;  
      }

      const std::vector<uint16_t>&  rawdigiAmplitudesR_ = siStripClusterInfo.rawdigiAmplitudesR();
      for (int i=0;i<rawdigiAmplitudesR_.size();i++){
	std::cout << "\t\t\t rawdigi amplitudes R " << i << " " << rawdigiAmplitudesR_[i]<< std::endl;  
      }
      const std::vector<float>&  stripNoises_ = siStripClusterInfo.stripNoises();
      for (int i=0;i<stripNoises_.size();i++){
	std::cout << "\t\t\t strip noise " << i << " " << stripNoises_[i]<< std::endl;  
      }
    }
  }  
}

//-------------------------------

void find_cluster(int detid=0){
  if (DEBUG>1){
    if (detid)
      std::cout << "[find_cluster] detid " << detid << std::endl;
    else
      std::cout << "cluster size " <<    dsvSiStripClusterInfo.size() << std::endl;
  }

  for (int idetset=0; idetset<dsvSiStripClusterInfo.size();idetset++){
    
    int thedetid = dsvSiStripClusterInfo[idetset].id;

    if ( detid==0 || detid == thedetid ){
      std::cout << "detid " << idetset << " " << thedetid<< std::endl;
      cluster_study(idetset);
    }
  }
}

//-------------------------------

void track_study(){

  int nTracks=vTracks.size();
  std::cout << "Number of tracks = " << nTracks << std::endl;
  ((TH1F*) Hlist.FindObject("nTracks"))->Fill(nTracks);

  for ( unsigned int iTrack = 0;  iTrack < nTracks; iTrack++ ) {
    const reco::Track& track = vTracks[iTrack];

    std::cout << "Track number "<< iTrack << std::endl ;
    if (DEBUG>2){
      std::cout<< " pt = " << track.pt() << std::endl;    
      std::cout << "\tmomentum: " << track.momentum().x()<< " " << track.momentum().y()<< " " << track.momentum().z()<< std::endl;
      std::cout << "\tPT: " << track.pt()<< std::endl;
      //std::cout << "\tvertex: " << track.vertex().x() << " " << track.vertex.y() << " " << track.vertex.z()<< std::endl;
      std::cout << "\timpact parameter: " << track.d0()<< std::endl;
      std::cout << "\tcharge: " << track.charge()<< std::endl;
      //std::cout << "\tnormalizedChi2: " << track.normalizedChi2()<< std::endl;      
      cout<<"\tFrom EXTRA : "<<endl;
      cout<<"\t\touter PT "<< track.outerPt()<<endl;
    }
    
    trackingRecHit_iterator dummy = track.recHitsBegin(); //!!! I don't know why, but without this the following crashes !!!!
    int recHitsSize=track.recHitsSize();
    cout <<"\t\tNumber of RecHits "<<recHitsSize<<endl;
    ((TH1F*) Hlist.FindObject("nRecHits"))->Fill(recHitsSize);
  
    for (size_t i = 0; i < recHitsSize; i++){
      TrackingRecHitRef recHit=track.recHit(i);
      if (recHit->isValid()){
	if (DEBUG>1){
	  std::cout <<"\t\t\tRecHit on det "<< recHit->geographicalId().rawId()<<std::endl;
	  std::cout <<"\t\t\tRecHit in LP "<<recHit->localPosition().x() << " " << recHit->localPosition().y () << " " << recHit->localPosition().z()<<std::endl;
	  std::cout <<"\t\t\tRecHit position error "<< recHit->localPositionError().xx() << " \t " << recHit->localPositionError().xy() << " \t " << recHit->localPositionError().yy()<<std::endl;
	}
	
	find_cluster(recHit->geographicalId().rawId());
      }
    }
  }
}

//-------------------------------

void get_trigger_bits(){
  if (DEBUG>2)
    std::cout<<"[LTCDigiCollection]: size "<<ltcdigis.size()<<std::endl;
  for(unsigned int iltc =0; iltc < ltcdigis.size(); iltc++)
    {
      LTCDigi ltcdigi = ltcdigis[iltc];
      for (int i = 0; i < 6; i++){
	if (DEBUG>1)
	  cout<<"[LTCDigi]: bit "<<i <<" has Trigger  "<<ltcdigi.HasTriggered(i)<<endl;
  
	if(ltcdigi.HasTriggered(i)){ 
	  totTrigs[i]++;
	  ((TH1F*) Hlist.FindObject("TriggerBits"))->Fill(i);
	}
      }    
    }
}

//-------------------------------

void loop(){
  std::cout << "\nLoop on events" << std::endl;
  for ( unsigned int ievent = 0; ievent < events->GetEntries(); ievent++ ) {
    if(ievent==0){
      for(int i=0;i<6;i++) totTrigs[i]=0;
    }
    std::cout << ".... Event: " << ievent << std::endl;


    LTCDigiB->GetEntry(ievent);
    TracksB->GetEntry(ievent);
    SiStripClusterInfoB->GetEntry(ievent);
    events->GetEntry(ievent,0);

    TracksB->SetAddress(&vTracks);
    SiStripClusterInfoB->SetAddress(&dsvSiStripClusterInfo);
    LTCDigiB->SetAddress(&ltcdigis);
    get_trigger_bits();
    track_study();
    //find_cluster(); //Now called by track_study
  }//close ievent loop
}

//-------------------------------

SiStripMacroExample(TString afilename, TString bfilename="pippo"){
  
  outfilename=bfilename;

  book_histos();

  file= new TRFIOFile(afilename);

  list_of_branches();

  set_branches();

  loop();

  cout << "\n\nFinal Trigger Summary events read = " << events->GetEntries()<< endl;
  for (int i = 0; i < 6; i++)
    cout<<"[LTCDigi]: bit "<< i <<" has Trigger "<<totTrigs[i] << " times"<<endl;

  save_histos(outfilename);
}



