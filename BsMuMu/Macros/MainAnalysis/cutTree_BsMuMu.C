void cutTree_BsMuMu( TString eb = "0", TString index = "0", const TString & inFile="tmva-trees-3000.root", TString outFile="data_afterCuts_" )
{
  cout << "eb = " << eb << endl;
  cout << "index = " << index << endl;
  outFile += eb+"_"+index+".root";
  TFile infile(inFile);
  gDirectory->Cd("sidebandChan"+eb+"Events"+index);
  TTree* inTree = (TTree*)gROOT->FindObject("events");
  if(!inTree){
    cout<<"Could not access yield tree!"<<endl;
    return;
  }

  TFile outfile(outFile,"recreate");

  TString lifetimeCuts("pvip < 0.008 && pvips < 2.000");
  TString trackCuts("iso > 0.8 && docatrk > 0.015 && closetrk < 2");
  TString barrelCuts("m1eta < 1.4 && m1eta > -1.4 && m2eta < 1.4 && m2eta > -1.4 && ((m1pt > m2pt && m1pt > 4.5 && m2pt > 4.0)||(m1pt < m2pt && m1pt > 4.0 && m2pt > 4.5))");
  TString endcapsCuts("!(m1eta < 1.4 && m1eta > -1.4 && m2eta < 1.4 && m2eta > -1.4) && ((m1pt > m2pt && m1pt > 4.5 && m2pt > 4.2)||(m1pt < m2pt && m1pt > 4.2 && m2pt > 4.5))");
  TString candidateBarrelCuts("pt > 6.5 && alpha < 0.050 && fls3d > 13.0 && chi2/dof < 2.2");
  TString candidateEndcapsCuts("pt > 8.5 && alpha < 0.030 && fls3d > 15.0 && chi2/dof < 1.8");


  TTree* outTree = (TTree*)inTree->CopyTree(lifetimeCuts+"&&"+trackCuts+"&& (("+barrelCuts+"&&"+candidateBarrelCuts+")"+"||"+"("+endcapsCuts+"&&"+candidateEndcapsCuts+"))");


  //close everything
  outTree->Write();
  outfile.Close();
  infile.Close();
}
