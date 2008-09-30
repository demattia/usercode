{
  gROOT->ProcessLine(".L tdrstyle.C");
  gROOT->ProcessLine(".L plot.C");
  gROOT->ProcessLine(".L plot_backup.C");

  gROOT->Reset();

  // Set the appropiate style
  setTDRStyle();

  TFile data       = TFile("hdata_ZS_wzero.root", "read");
  TFile simulation = TFile("hsim_ConfC_3CC.root", "read");

  TCanvas * c = 0;

  TString histoName = "eta_s_all_TIB_on";
  // Read and normalize histograms
  TH1F * histoData       = ((TH1F*)(data.Get(histoName)));
  TH1F * histoSimulation = ((TH1F*)(simulation.Get(histoName)));

  //plotSimData( histoData, "Data ZS", histoSimulation, "MC", "TIB" );

  TString histoName = "eta_s_all_TOB_on";
  // Read and normalize histograms
  histoData       = ((TH1F*)(data.Get(histoName)));
  histoSimulation = ((TH1F*)(simulation.Get(histoName)));

  //plotSimData( histoData, "Data ZS", histoSimulation, "MC", "TOB" );

  /*
  // Comparison data for tracks and data for no-tracks clusters
  // ----------------------------------------------------------
  TFile vrData = TFile("hdata_VR.root", "read");
  histoData           = ((TH1F*)(data.Get("size_TIB_on")));
  TH1F * histoDataCut = ((TH1F*)(data.Get("size_cut_TIB_on")));
  c = plot( histoData, "All tracks", histoDataCut, "Perpendicular tracks", "SizeTIB", 0.6,
            "Cluster Width[strips]", "Number of clusters", "", true, "TIB", 0.7, 0.4 );
  c->Print("ClusterWidthTrack_TIB.pdf");

  histoData    = ((TH1F*)(data.Get("size_TOB_on")));
  histoDataCut = ((TH1F*)(data.Get("size_cut_TOB_on")));
  c = plot( histoData, "All tracks", histoDataCut, "Perpendicular tracks", "SizeTOB", 0.6,
            "Cluster Width[strips]", "Number of clusters", "", true, "TOB", 0.7, 0.4 );
  c->Print("ClusterWidthTrack_TOB.pdf");

  */

  // Data-Simulation  comparison for cluster size
  // --------------------------------------------
  histoData       = ((TH1F*)(data.Get("size_cut_TIB_on")));
  histoSimulation = ((TH1F*)(simulation.Get("size_cut_TIB_on")));
  c = plotSimData( histoData, "Data ZS", histoSimulation, "MC", "sizePerpTIB", 0.6,
                   "Cluster Width[strips]", "Number of clusters", "", true, "TIB", 0.7, 0.4, 0, true );
  c->Print("ClusterWidthTrackPerp_TIB.pdf");

  histoData       = ((TH1F*)(data.Get("size_cut_TOB_on")));
  histoSimulation = ((TH1F*)(simulation.Get("size_cut_TOB_on")));
  c = plotSimData( histoData, "Data ZS", histoSimulation, "MC", "sizePerpTOB", 0.6,
                   "Cluster Width[strips]", "Number of clusters", "", true, "TOB", 0.7, 0.4, 0, true );
  c->Print("ClusterWidthTrackPerp_TOB.pdf");

  /*
  // XZ angle
  // --------
  histoData       = ((TH1F*)(data.Get("angleXZ_TIB_on")));
  histoSimulation = ((TH1F*)(simulation.Get("angleXZ_TIB_on")));
  c = plotSimData( histoData, "Data ZS", histoSimulation, "MC", "XZangleTIB", 0.07,
                   "XZ Angle (deg)", "", "", false, "TIB", 0.7, 0.7, 5 );
  c->Print("Angle_XZ_TIB.pdf");

  histoData       = ((TH1F*)(data.Get("angleXZ_TOB_on")));
  histoSimulation = ((TH1F*)(simulation.Get("angleXZ_TOB_on")));
  c = plotSimData( histoData, "Data ZS", histoSimulation, "MC", "XZangleTOB", 0.07,
                   "XZ Angle (deg)", "", "", false, "TOB", 0.7, 0.7, 5 );
  c->Print("Angle_XZ_TOB.pdf");
  */

}
