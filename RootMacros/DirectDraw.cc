{
  // A simple macro usefull for quick looks at files
  // Shows how to take a TChain and create the aliases for the branches

  TChain file("Events");
  file.Add("../QCD_30-50_1.root");
  file.Add("../QCD_30-50_2.root");
  file.Add("../QCD_30-50_11.root");
  file.Add("../QCD_30-50_12.root");
  file.Add("../QCD_30-50_13.root");
  file.Add("../QCD_30-50_14.root");
  file.Add("../QCD_30-50_15.root");
  file.Add("../QCD_30-50_16.root");
  file.Add("../QCD_30-50_17.root");
  file.Add("../QCD_30-50_18.root");
  file.Add("../QCD_30-50_19.root");
  file.Add("../QCD_30-50_1.root");
  file.Add("../QCD_30-50_20.root");
  file.Add("../QCD_30-50_21.root");
  file.Add("../QCD_30-50_22.root");
  file.Add("../QCD_30-50_23.root");
  file.Add("../QCD_30-50_24.root");
  file.Add("../QCD_30-50_25.root");
  file.Add("../QCD_30-50_26.root");
  file.Add("../QCD_30-50_27.root");
  file.Add("../QCD_30-50_28.root");
  file.Add("../QCD_30-50_29.root");
  file.Add("../QCD_30-50_2.root");
  file.Add("../QCD_30-50_30.root");
  file.Add("../QCD_30-50_31.root");
  file.Add("../QCD_30-50_32.root");
  file.Add("../QCD_30-50_33.root");
  file.Add("../QCD_30-50_34.root");
  file.Add("../QCD_30-50_35.root");
  file.Add("../QCD_30-50_36.root");
  file.Add("../QCD_30-50_37.root");
  file.Add("../QCD_30-50_3.root");
  file.Add("../QCD_30-50_40.root");
  file.Add("../QCD_30-50_41.root");
  file.Add("../QCD_30-50_43.root");
  file.Add("../QCD_30-50_44.root");
  file.Add("../QCD_30-50_46.root");
  file.Add("../QCD_30-50_47.root");
  file.Add("../QCD_30-50_48.root");
  file.Add("../QCD_30-50_49.root");
  file.Add("../QCD_30-50_4.root");
  file.Add("../QCD_30-50_50.root");
  file.Add("../QCD_30-50_51.root");
  file.Add("../QCD_30-50_52.root");
  file.Add("../QCD_30-50_53.root");
  file.Add("../QCD_30-50_54.root");
  file.Add("../QCD_30-50_57.root");
  file.Add("../QCD_30-50_58.root");
  file.Add("../QCD_30-50_5.root");
  file.Add("../QCD_30-50_60.root");
  file.Add("../QCD_30-50_7.root");
  file.Add("../QCD_30-50_8.root");
  file.Add("../QCD_30-50_9.root");

  Events->SetAlias("cenjets", "anaobjBaseJets_offlineProd_cenJets_PROD.obj"); 
  Events->SetAlias("taujets", "anaobjBaseJets_offlineProd_tauJets_PROD.obj"); 
  Events->SetAlias("forjets", "anaobjBaseJets_offlineProd_forJets_PROD.obj"); 
  Events->SetAlias("met", "anaobjBaseMEt_offlineProd_l1MEt_PROD.obj"); 

  // Draw with some conditions
  Events->Draw("met.et_", "cenjets.et_>100 || taujets.et_>100"); 
}
