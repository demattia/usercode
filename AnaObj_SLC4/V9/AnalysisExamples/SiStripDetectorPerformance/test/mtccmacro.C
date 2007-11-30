#include "langaus.C"

TFile	    *poFileIn;
TCanvas     *poCanvas;
TPaveLabel  *title;

TTree	    *poTrackTree;
TTree	    *poLATree;

char cMaxPadsX = 2;
char cMaxPadsY = 2;
char cCurPad   = 0;
bool bCanvasSaved = false;
std::string oPageTitle = "";

void createPagePre() {
  if( 0 == cCurPad) {
    poCanvas->Clear();
    poCanvas->Divide( cMaxPadsX, cMaxPadsY);
    title->SetLabel( oPageTitle.c_str());
    title->Draw();
  }
  
  ++cCurPad;
  poCanvas->cd( cCurPad);

  bCanvasSaved = false;
}

void createPagePost( const char *pcFILE_OUT) {
  poCanvas->Update();

  if( cMaxPadsX * cMaxPadsY == cCurPad) {
    poCanvas->Print( pcFILE_OUT);
    bCanvasSaved = true;
    // Divide Canvas
    cCurPad = 0;
  }
}

void insertPageBreak( const char *pcFILE_OUT) {
  cCurPad = cMaxPadsX * cMaxPadsY;
  createPagePost( pcFILE_OUT);
}

void setPadsLayout( const char *pcFILE_OUT,
		    const char cMAX_PADSX,
		    const char cMAX_PADSY) {
  insertPageBreak( pcFILE_OUT);
  cMaxPadsX = cMAX_PADSX;
  cMaxPadsY = cMAX_PADSY;
  cCurPad = 0;
}

void setPageTitle( const char *pcPAGE_TITLE) {
  oPageTitle = pcPAGE_TITLE;
}

void createPage( TTree *poTree,
		 const char *pcFILE_OUT,
		 const std::string &roFUNC,
		 const char cREBIN = 0) {

  createPagePre();
  poTree->Draw( roFUNC.c_str());
  createPagePost( pcFILE_OUT);
}

void createPage( TTree *poTree,
		 const char *pcFILE_OUT,
		 const std::string &roFUNC,
		 const std::string &roARGS,
		 const char cREBIN = 0) {

  createPagePre();
  if( cREBIN) {
    poTree->Draw( roFUNC.c_str(), roARGS.c_str());
    TH1F *htemp = dynamic_cast<TH1F *>( gPad->GetPrimitive("htemp"));
    TH1F *hnew =  dynamic_cast<TH1F*>( htemp->Rebin(3,"hnew"));
    hnew->Draw();
  } else {
    poTree->Draw( roFUNC.c_str(), roARGS.c_str());
  }
  createPagePost( pcFILE_OUT);
}

void closeFile( char *pcFILE_OUT) {
  if( !bCanvasSaved) {
    poCanvas->Print( pcFILE_OUT);
  }

  {
    std::string oFile( pcFILE_OUT);
    oFile += "]";
    poCanvas->Print( oFile.c_str());
  }

}

void mtccmacro( const char *pcFILE_IN,
		const char *pcFILE_OUT = "out.ps") {
  SetStyle();
  poFileIn  = new TRFIOFile( pcFILE_IN);

  poCanvas = new TCanvas();
  poCanvas->Draw();

  title = new TPaveLabel( 0.01, 0.01, 0.9, 0.04, oPageTitle.c_str());
  title->SetFillColor(0);
  title->SetTextColor(1);
  title->SetTextFont(52);
  title->SetBorderSize(0);
  title->Draw();

  poTrackTree = dynamic_cast<TTree *>( poFileIn->Get( "TrackTree"));
  poLATree    = dynamic_cast<TTree *>( poFileIn->Get( "MTCCNtupleMakerTree"));

  {
    std::string oFile( pcFILE_OUT);
    oFile += "[";
    poCanvas->Print( oFile.c_str());
  }

  setPageTitle( "Track: pt and eta");
  createPage( poTrackTree, pcFILE_OUT, "pt", "pt<100");
  createPage( poTrackTree, pcFILE_OUT, "pt", "pt<100 && hitspertrack>3");
  createPage( poTrackTree, pcFILE_OUT, "eta", "pt<100");
  createPage( poTrackTree, pcFILE_OUT, "eta", "pt<100 && hitspertrack>3");

  setPageTitle("Tracks: phi, hits per track");
  createPage( poTrackTree, pcFILE_OUT, "phi", "pt<100");
  createPage( poTrackTree, pcFILE_OUT, "hitspertrack", "pt<100");	
  createPage( poTrackTree, pcFILE_OUT, "hitspertrack:eta", "pt<100");
  createPage( poTrackTree, pcFILE_OUT, "hitspertrack:phi", "pt<100");

  setPageTitle("Tracks: chi2");
  createPage( poTrackTree, pcFILE_OUT, "hitspertrack:chi2", "pt<100 && chi2<50");
  createPage( poTrackTree, pcFILE_OUT, "pt:chi2", "pt<100 && chi2<50");
  createPage( poTrackTree, pcFILE_OUT, "eta:chi2", "pt<100 && chi2<50");
  createPage( poTrackTree, pcFILE_OUT, "phi:chi2", "pt<100 && chi2<50");

  // cluster leaves: angle, bTrack, bwfw, charge, chi2, clusterchg, clusterchgl,
  // clusterchgr, clustereta, clustermaxchg, clusternoise, clusterpos, eta, event, 
  // extint, hitspertrack, layer, localmagfield, module, momentum, monostereo, 
  // ndof, normchi2, phi, pt, rod, run, sign, size, stereocorrection, string, type, 
  // bTriggerDT, bTriggerCSC, BTriggerRBC1, bTriggerRBC2, bTRiggerRPC, eventcounter
  // wheel

  // hists: # of clusters per layer, # of track clusters / total # of clusters
  // # of hits per track, phi, cluster size vs angle, angle per layer, cluster charge,
  // cluster charge vs angle
 
  setPageTitle("track phi for layers");
  createPage( poLATree, pcFILE_OUT, "phi", "pt<100 && rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "phi", "pt<100 && rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "phi", "pt<100 && rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "phi", "pt<100 && rod>0 && layer==4 && bTrack==1");

  setPageTitle("track eta for layers");
  createPage( poLATree, pcFILE_OUT, "eta", "pt<100 && rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "eta", "pt<100 && rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "eta", "pt<100 && rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "eta", "pt<100 && rod>0 && layer==4 && bTrack==1");
  
  setPageTitle("hitspertrack for layers");
  createPage( poLATree, pcFILE_OUT, "hitspertrack", "pt<100 && rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "hitspertrack", "pt<100 && rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "hitspertrack", "pt<100 && rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "hitspertrack", "pt<100 && rod>0 && layer==4 && bTrack==1");

  setPageTitle("chi2 for layers");
  createPage( poLATree, pcFILE_OUT, "chi2", "pt<100 && rod<0 && layer==2 && bTrack==1 && chi2<50");
  createPage( poLATree, pcFILE_OUT, "chi2", "pt<100 && rod<0 && layer==3 && bTrack==1 && chi2<50");
  createPage( poLATree, pcFILE_OUT, "chi2", "pt<100 && rod>0 && layer==3 && bTrack==1 && chi2<50");
  createPage( poLATree, pcFILE_OUT, "chi2", "pt<100 && rod>0 && layer==4 && bTrack==1 && chi2<50");
  
  setPageTitle("track pt vs eta for layers");
  createPage( poLATree, pcFILE_OUT, "pt:eta", "pt<50 && rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "pt:eta", "pt<50 && rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "pt:eta", "pt<50 && rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "pt:eta", "pt<50 && rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster position for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod<0 && layer==2");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod<0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod>0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod>0 && layer==4");

  setPageTitle("cluster position for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusterpos", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge vs position for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod<0 && layer==2 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod<0 && layer==3 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod>0 && layer==3 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod>0 && layer==4 && clusterchg<500");
  
  setPageTitle("cluster charge vs position for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod<0 && layer==2 && bTrack==1 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod<0 && layer==3 && bTrack==1 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod>0 && layer==3 && bTrack==1 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:clusterpos", "rod>0 && layer==4 && bTrack==1 && clusterchg<500");

  setPageTitle("cluster noise vs position for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod<0 && layer==2");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod<0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod>0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod>0 && layer==4");

  setPageTitle("cluster noise vs position for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:clusterpos", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge vs chi2 for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusterchg:chi2", "rod<0 && layer==2 && bTrack==1 && clusterchg<500 && chi2<50");
  createPage( poLATree, pcFILE_OUT, "clusterchg:chi2", "rod<0 && layer==3 && bTrack==1 && clusterchg<500 && chi2<50");
  createPage( poLATree, pcFILE_OUT, "clusterchg:chi2", "rod>0 && layer==3 && bTrack==1 && clusterchg<500 && chi2<50");
  createPage( poLATree, pcFILE_OUT, "clusterchg:chi2", "rod>0 && layer==4 && bTrack==1 && clusterchg<500 && chi2<50");

  setPageTitle("cluster size for all clusters");
  createPage( poLATree, pcFILE_OUT, "size", "rod<0 && layer==2 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "size", "rod<0 && layer==3 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "size", "rod>0 && layer==3 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "size", "rod>0 && layer==4 && clusterchg<500");

  setPageTitle("cluster size for track clusters");
  createPage( poLATree, pcFILE_OUT, "size", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "size", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "size", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "size", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge-to-noise for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod<0 && layer==2 && clusterchg<500 && module!= 369197581");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod<0 && layer==3 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod>0 && layer==3 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod>0 && layer==4 && clusterchg<500 && module != 436371976 && module != 436371984 && module != 436371716");

  setPageTitle("cluster charge-to-noise for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod<0 && layer==2 && bTrack==1 && clusterchg<500 && module!= 369197581");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod<0 && layer==3 && bTrack==1 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod>0 && layer==3 && bTrack==1 && clusterchg<500");
  createPage( poLATree, pcFILE_OUT, "clusterchg/clusternoise", "rod>0 && layer==4 && bTrack==1 && clusterchg<500 && module != 436371976 && module != 436371984 && module != 436371716");

  setPageTitle("track angle");
  createPage( poLATree, pcFILE_OUT, "angle", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "angle", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "angle", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "angle", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster noise vs size number for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod<0 && layer==2 && size<10");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod<0 && layer==3 && size<10");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod>0 && layer==3 && size<10");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod>0 && layer==4 && size<10");

  setPageTitle("cluster noise vs size number for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod<0 && layer==2 && bTrack==1 && size<10");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod<0 && layer==3 && bTrack==1 && size<10");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod>0 && layer==3 && bTrack==1 && size<10");
  createPage( poLATree, pcFILE_OUT, "clusternoise:size", "rod>0 && layer==4 && bTrack==1 && size<10");
  
  setPageTitle("module number");
  createPage( poLATree, pcFILE_OUT, "module-369000000", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "module-369000000", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "module-436000000", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "module-436000000", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge vs module number");
  createPage( poLATree, pcFILE_OUT, "clusterchg:(module-369000000)", "clusterchg < 500 && rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusterchg:(module-369000000)", "clusterchg < 500 && rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusterchg:(module-436000000)", "clusterchg < 500 && rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusterchg:(module-436000000)", "clusterchg < 500 && rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster noise vs module number for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-369000000)", "rod<0 && layer==2");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-369000000)", "rod<0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-436000000)", "rod>0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-436000000)", "rod>0 && layer==4");

  setPageTitle("cluster noise vs module number for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-369000000)", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-369000000)", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-436000000)", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:(module-436000000)", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster noise vs run number for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod<0 && layer==2");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod<0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod>0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod>0 && layer==4");

  setPageTitle("cluster noise vs run number for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod<0 && layer==2 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod<0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod>0 && layer==3 && bTrack==1");
  createPage( poLATree, pcFILE_OUT, "clusternoise:run", "rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge vs eventcounter number for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusterchg:eventcounter", "rod<0 && layer==2 && bTrack==1 && clusterchg < 500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:eventcounter", "rod<0 && layer==3 && bTrack==1 && clusterchg < 500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:eventcounter", "rod>0 && layer==3 && bTrack==1 && clusterchg < 500");
  createPage( poLATree, pcFILE_OUT, "clusterchg:eventcounter", "rod>0 && layer==4 && bTrack==1 && clusterchg < 500");

  setPageTitle("cluster noise vs eventcounter for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise:eventcounter", "rod<0 && layer==2");
  createPage( poLATree, pcFILE_OUT, "clusternoise:eventcounter", "rod<0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:eventcounter", "rod>0 && layer==3");
  createPage( poLATree, pcFILE_OUT, "clusternoise:eventcounter", "rod>0 && layer==4");

  setPageTitle("cluster eta for all clusters");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod<0 && layer==2 && clustereta>0");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod<0 && layer==3 && clustereta>0");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod>0 && layer==3 && clustereta>0");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod>0 && layer==4 && clustereta>0");

  setPageTitle("cluster eta for track clusters");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod<0 && layer==2 && bTrack==1 && clustereta>0");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod<0 && layer==3 && bTrack==1 && clustereta>0");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod>0 && layer==3 && bTrack==1 && clustereta>0");
  createPage( poLATree, pcFILE_OUT, "clustereta", "rod>0 && layer==4 && bTrack==1 && clustereta>0");

  setPageTitle("cluster noise for all clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod<0 && layer==2", 1);
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod<0 && layer==3", 1);
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod>0 && layer==3", 1);
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod>0 && layer==4", 1);

  setPageTitle("cluster noise for track clusters");
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod<0 && layer==2 && bTrack==1", 1);
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod<0 && layer==3 && bTrack==1", 1);
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod>0 && layer==3 && bTrack==1", 1);
  createPage( poLATree, pcFILE_OUT, "clusternoise", "rod>0 && layer==4 && bTrack==1", 1);


  setPageTitle( "Track: various parameters");
  createPage( poTrackTree, pcFILE_OUT, "phi", "pt<200");
  createPage( poTrackTree, pcFILE_OUT, "hitspertrack");	
  createPage( poTrackTree, pcFILE_OUT, "hitspertrack:eta");
  insertPageBreak( pcFILE_OUT);

  // Fits
  setPageTitle("cluster charge for all clusters");
  fitCharge( pcFILE_OUT, "TIB L2", "clusterchg", "clusterchg < 500 && rod<0 && layer==2");
  fitCharge( pcFILE_OUT, "TIB L3", "clusterchg", "clusterchg < 500 && rod<0 && layer==3");
  fitCharge( pcFILE_OUT, "TOB L1", "clusterchg", "clusterchg < 500 && rod>0 && layer==3");
  fitCharge( pcFILE_OUT, "TOB L5", "clusterchg", "clusterchg < 500 && rod>0 && layer==4");

  setPageTitle("cluster charge for track clusters");
  fitCharge( pcFILE_OUT, "TIB L2", "clusterchg", "clusterchg < 500 && rod<0 && layer==2 && bTrack==1");
  fitCharge( pcFILE_OUT, "TIB L3", "clusterchg", "clusterchg < 500 && rod<0 && layer==3 && bTrack==1");
  fitCharge( pcFILE_OUT, "TOB L1", "clusterchg", "clusterchg < 500 && rod>0 && layer==3 && bTrack==1");
  fitCharge( pcFILE_OUT, "TOB L5", "clusterchg", "clusterchg < 500 && rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge-to-noise for all clusters");
  fitChargeN( pcFILE_OUT, "TIB L2", "clusterchg/clusternoise", "rod<0 && layer==2 && clusterchg<500 && module!= 369197581");
  fitChargeN( pcFILE_OUT, "TIB L3", "clusterchg/clusternoise", "rod<0 && layer==3 && clusterchg<500");
  fitChargeN( pcFILE_OUT, "TOB L1", "clusterchg/clusternoise", "rod>0 && layer==3 && clusterchg<500");
  fitChargeN( pcFILE_OUT, "TOB L5", "clusterchg/clusternoise", "rod>0 && layer==4 && clusterchg<500 && module != 436371976 && module != 436371984 && module != 436371716");

  setPageTitle("cluster charge-to-noise for track clusters");
  fitChargeN( pcFILE_OUT, "TIB L2", "clusterchg/clusternoise", "rod<0 && layer==2 && bTrack==1 && clusterchg<500 && module!= 369197581");
  fitChargeN( pcFILE_OUT, "TIB L3", "clusterchg/clusternoise", "rod<0 && layer==3 && bTrack==1 && clusterchg<500");
  fitChargeN( pcFILE_OUT, "TOB L1", "clusterchg/clusternoise", "rod>0 && layer==3 && bTrack==1 && clusterchg<500");
  fitChargeN( pcFILE_OUT, "TOB L5", "clusterchg/clusternoise", "rod>0 && layer==4 && bTrack==1 && clusterchg<500 && module != 436371976 && module != 436371984 && module != 436371716");

  // Normalized clusterCharge values
  setPageTitle("cluster charge for track clusters (normalized)");
  fitCharge( pcFILE_OUT, "TIB L2", "clusterchg * cos( 1.0 * angle / 180 * 3.14)", "clusterchg * cos( 1.0 * angle / 180 * 3.14) < 500 && rod<0 && layer==2 && bTrack==1");
  fitCharge( pcFILE_OUT, "TIB L3", "clusterchg * cos( 1.0 * angle / 180 * 3.14)", "clusterchg * cos( 1.0 * angle / 180 * 3.14) < 500 && rod<0 && layer==3 && bTrack==1");
  fitCharge( pcFILE_OUT, "TOB L1", "clusterchg * cos( 1.0 * angle / 180 * 3.14)", "clusterchg * cos( 1.0 * angle / 180 * 3.14) < 500 && rod>0 && layer==3 && bTrack==1");
  fitCharge( pcFILE_OUT, "TOB L5", "clusterchg * cos( 1.0 * angle / 180 * 3.14)", "clusterchg * cos( 1.0 * angle / 180 * 3.14) < 500 && rod>0 && layer==4 && bTrack==1");

  setPageTitle("cluster charge-to-noise for track clusters (normalized)");
  fitChargeN( pcFILE_OUT, "TIB L2", "clusterchg * cos( 1.0 * angle / 180 * 3.14)/clusternoise", "rod<0 && layer==2 && bTrack==1 && clusterchg<500 && module!= 369197581");
  fitChargeN( pcFILE_OUT, "TIB L3", "clusterchg * cos( 1.0 * angle / 180 * 3.14)/clusternoise", "rod<0 && layer==3 && bTrack==1 && clusterchg<500");
  fitChargeN( pcFILE_OUT, "TOB L1", "clusterchg * cos( 1.0 * angle / 180 * 3.14)/clusternoise", "rod>0 && layer==3 && bTrack==1 && clusterchg<500");
  fitChargeN( pcFILE_OUT, "TOB L5", "clusterchg * cos( 1.0 * angle / 180 * 3.14)/clusternoise", "rod>0 && layer==4 && bTrack==1 && clusterchg<500 && module != 436371976 && module != 436371984 && module != 436371716");


  closeFile( pcFILE_OUT);

  poFileIn->Close();

  delete poCanvas;
  delete poFileIn;

  poTrackTree = 0;
  poLATree    = 0;
}

void fitCharge( const char *pcFILE_OUT,
		std::string oHistName,
		std::string oLeafName,
		std::string oCuts) {
  gStyle->SetStatH(0.2);  
  gStyle->SetStatW(0.2);
  createPagePre();

  TH1F* chrg = new TH1F( "chrg", oHistName.c_str(),  50, 0., 500.);
  oLeafName += ">>chrg";
  poLATree->Draw( oLeafName.c_str(), oCuts.c_str());

  langaus( chrg);

  createPagePost( pcFILE_OUT);
  gStyle->SetStatH(0.15);  
  gStyle->SetStatW(0.15);
}

void fitChargeN( const char *pcFILE_OUT,
		std::string oHistName,
		std::string oLeafName,
		std::string oCuts) {
  gStyle->SetStatH(0.2);  
  gStyle->SetStatW(0.2);
  createPagePre();

  TH1F* chrg = new TH1F( "chrg", oHistName.c_str(),  50, 0., 100.);
  oLeafName += ">>chrg";
  poLATree->Draw( oLeafName.c_str(), oCuts.c_str());

  langausN( chrg);

  createPagePost( pcFILE_OUT);
  gStyle->SetStatH(0.15);  
  gStyle->SetStatW(0.15);
}

void SetStyle() {
  TStyle *CMSStyle = new TStyle("CMS-Style","The Perfect Style for Plots ;-)");
  gStyle = CMSStyle;

  // Canvas
  CMSStyle->SetCanvasColor     (0);
  CMSStyle->SetCanvasBorderSize(0);
  CMSStyle->SetCanvasBorderMode(0);
  CMSStyle->SetCanvasDefH      (700);
  CMSStyle->SetCanvasDefW      (700);
  CMSStyle->SetCanvasDefX      (100);
  CMSStyle->SetCanvasDefY      (100);

  // Pads
  CMSStyle->SetPadColor       (0);
  CMSStyle->SetPadBorderSize  (0);
  CMSStyle->SetPadBorderMode  (0);
  CMSStyle->SetPadBottomMargin(0.18);
  CMSStyle->SetPadTopMargin   (0.08);
  CMSStyle->SetPadLeftMargin  (0.16);
  CMSStyle->SetPadRightMargin (0.05);
  CMSStyle->SetPadGridX       (0);
  CMSStyle->SetPadGridY       (0);
  CMSStyle->SetPadTickX       (0);
  CMSStyle->SetPadTickY       (0);

  // Frames
  CMSStyle->SetFrameFillStyle ( 0);
  CMSStyle->SetFrameFillColor ( 0);
  CMSStyle->SetFrameLineColor ( 1);
  CMSStyle->SetFrameLineStyle ( 0);
  CMSStyle->SetFrameLineWidth ( 1);
  CMSStyle->SetFrameBorderSize(10);
  CMSStyle->SetFrameBorderMode( 0);

  // Histograms
  CMSStyle->SetHistFillColor(2);
  CMSStyle->SetHistFillStyle(0);
  CMSStyle->SetHistLineColor(1);
  CMSStyle->SetHistLineStyle(0);
  CMSStyle->SetHistLineWidth(2);
  CMSStyle->SetNdivisions(505);

  // Functions
  CMSStyle->SetFuncColor(1);
  CMSStyle->SetFuncStyle(0);
  CMSStyle->SetFuncWidth(2);

  // Various
  CMSStyle->SetTitleSize  (0.050,"X");
  CMSStyle->SetTitleOffset(1.000,"X");
  CMSStyle->SetTitleFillColor (0);
  CMSStyle->SetLabelOffset(0.003,"X");
  CMSStyle->SetLabelSize  (0.050,"X");
  CMSStyle->SetLabelFont  (42   ,"X");

  CMSStyle->SetStripDecimals(kFALSE);

  CMSStyle->SetTitleSize  (0.050,"Y");
  CMSStyle->SetTitleOffset(1.500,"Y");
  CMSStyle->SetLabelOffset(0.008,"Y");
  CMSStyle->SetLabelSize  (0.050,"Y");
  CMSStyle->SetLabelFont  (42   ,"Y");

  CMSStyle->SetOptStat("emr");
  CMSStyle->SetOptFit();
  CMSStyle->SetStatColor(0);
  CMSStyle->SetStatBorderSize( 1);
  
  CMSStyle->SetTitleFont  (42);
  CMSStyle->SetTitleBorderSize( 1);

  CMSStyle->SetHistFillColor(kYellow);
}
