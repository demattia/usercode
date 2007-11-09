#ifndef HIVARIABLES_C
#define HIVARIABLES_C

#include "AnalysisExamples/AnalysisClasses/interface/HiVariables.h"

HiVariables::HiVariables(const char * Histo_Name)
{

  // This copies only the first 10 bytes
  strncpy(_histo_name, Histo_Name, 10);
  // Add the null termination character to the string
  _histo_name[10] = '\0';

  char SumEtName[30];
  char CentralityName[30];
  char AplanarityName[30];
  char SphericityName[30];
  char Y_SName[30];
  char ThrustName[30];
  char HeavyJetMassName[30];
  char LightJetMassName[30];
  char JetMassDifferenceName[30];
  char MajorBroadeningName[30];
  char MinorBroadeningName[30];
  char BroadeningDifferenceName[30];
  sprintf(SumEtName, "%s_SumEt", _histo_name);
  sprintf(CentralityName, "%s_Centrality", _histo_name);
  sprintf(AplanarityName, "%s_Aplanarity", _histo_name);
  sprintf(SphericityName, "%s_Sphericity", _histo_name);
  sprintf(Y_SName, "%s_Y_SName", _histo_name);
  sprintf(ThrustName, "%s_Thrust", _histo_name);
  sprintf(HeavyJetMassName, "%s_HeavyJetMass", _histo_name);
  sprintf(LightJetMassName, "%s_LightJetMass", _histo_name);
  sprintf(JetMassDifferenceName, "%s_JetMassDifference", _histo_name);
  sprintf(MajorBroadeningName, "%s_MajorBroadening", _histo_name);
  sprintf(MinorBroadeningName, "%s_MinorBroadening", _histo_name);
  sprintf(BroadeningDifferenceName, "%s_BroadeningDifference", _histo_name);

  SumEt = new TH1F ( SumEtName, SumEtName, 200, 0, 2000);
  Centrality = new TH1F ( CentralityName, CentralityName, 200, 0, 1);
  Aplanarity = new TH1F ( AplanarityName, AplanarityName, 200, 0, 0.5);
  Sphericity = new TH1F ( SphericityName, SphericityName, 200, 0, 1);
  Y_S = new TH2F (Y_SName, Y_SName, 100, -0.1, 1, 100, -0.6, 0.2);
  Thrust = new TH1F (ThrustName, ThrustName, 200, 0.5, 1.2);
  HeavyJetMass = new TH1F (HeavyJetMassName, HeavyJetMassName, 200, 0, 0.3);
  LightJetMass = new TH1F (LightJetMassName, LightJetMassName, 200, 0, 0.3);
  JetMassDifference = new TH1F (JetMassDifferenceName, JetMassDifferenceName, 200, 0, 0.1);
  MajorBroadening = new TH1F (MajorBroadeningName, MajorBroadeningName, 200, 0, 100.);
  MinorBroadening = new TH1F (MinorBroadeningName, MinorBroadeningName, 200, 0, 100.);
  BroadeningDifference = new TH1F (BroadeningDifferenceName, BroadeningDifferenceName, 200, 0, 100.);

  _thrust[0] = 0;
  _thrust[1] = 0;
  _thrust[2] = 0;
  _SumEt = 0;
  _Centrality = 0;
  _Sphericity = 0;
  _Aplanarity = 0;
  _Y = 0;
  _Thrust = 0;
  _HjetMass = 0;
  _LjetMass = 0;
  _MjetBroadening = 0;
  _mjetBroadening = 0;
}

void HiVariables::Fill( vector<SimpleJet> & vec_jet ) {
  // Order the jets in E (not Et) for the thrust
//   sort(vec_jet.begin(), vec_jet.end());

//   SimpleJet temp_jet;
//   for (unsigned int i=0; i<vec_jet.size(); ++i) {
//     for (int j=vec_jet.size()-2; j>-1; j=j-1) {
//       if (vec_jet.at(j+1).E() > vec_jet.at(j).E()) {
// 	temp_jet = vec_jet.at(j);
// 	vec_jet.at(j) = vec_jet.at(j+1);
// 	vec_jet.at(j+1) = temp_jet;
//       }
//     }
//   }

  sort( vec_jet.begin(), vec_jet.end() );
  reverse( vec_jet.begin(),vec_jet.end() );

  float ** momJB;
  float sum_E = 0;
  float sum_Et = 0;
  float sum_Pz = 0;
  for (unsigned int i=0; i<vec_jet.size(); ++i) {
    _tempJet = &(vec_jet.at(i));
    _temp_e = _tempJet->E();
    _temp_pt = _tempJet->pt();
    _temp_pz = _tempJet->pz();
    sum_E += _temp_e;
    sum_Et += _temp_pt;
    sum_Pz += _temp_pz;
  }
   SumEt->Fill(sum_Et);
//  _HiVar_SumEt_HistoSampler->Fill(sum_Et, fC, TOTEVENTS);
  _SumEt = sum_Et;
  _Centrality = sum_Et/sqrt(sum_E*sum_E - sum_Pz*sum_Pz);
   Centrality->Fill(_Centrality);
//  _HiVar_Centrality_HistoSampler->Fill(_Centrality, fC, TOTEVENTS);

  // Matrix evaluation
  // -----------------
  // boost vector
  momJB = new float*[vec_jet.size()];
  if (sum_E != 0.) {
    _bst.SetZ(-sum_Pz/sum_E);
    for (unsigned int i=0; i<vec_jet.size(); ++i) {
      _tempJet = &(vec_jet.at(i));
      _momJ.SetPx(_tempJet->px());
      _momJ.SetPy(_tempJet->py());
      _momJ.SetPz(_tempJet->pz());
      _momJ.SetE(_tempJet->E());
      // boost the vector
      _momJ.Boost(_bst);

      momJB[i] = new float[4];

      momJB[i][0] = _momJ.Px();
      momJB[i][1] = _momJ.Py();
      momJB[i][2] = _momJ.Pz();
      momJB[i][3] = _momJ.E();
    }
  }
  TMatrixD mom_ten_2(3,3);
  mom_ten_2(0,0) = 0;
  mom_ten_2(1,1) = 0;
  mom_ten_2(2,2) = 0;
  mom_ten_2(0,1) = 0;
  mom_ten_2(0,2) = 0;
  mom_ten_2(1,2) = 0;
  float sump2 = 0;
  for (unsigned int i=0; i<vec_jet.size(); ++i) {
    mom_ten_2(0,0) += momJB[i][0]*momJB[i][0];
    mom_ten_2(1,1) += momJB[i][1]*momJB[i][1];
    mom_ten_2(2,2) += momJB[i][2]*momJB[i][2];
    mom_ten_2(0,1) += momJB[i][0]*momJB[i][1];
    mom_ten_2(0,2) += momJB[i][0]*momJB[i][2];
    mom_ten_2(1,2) += momJB[i][1]*momJB[i][2];
    sump2 += momJB[i][3]*momJB[i][3];
  }
  if (sump2 != 0) {
    mom_ten_2(0,0) = mom_ten_2(0,0)/sump2;
    mom_ten_2(1,1) = mom_ten_2(1,1)/sump2;
    mom_ten_2(2,2) = mom_ten_2(2,2)/sump2;
    mom_ten_2(0,1) = mom_ten_2(0,1)/sump2;
    mom_ten_2(0,2) = mom_ten_2(0,2)/sump2;
    mom_ten_2(1,2) = mom_ten_2(1,2)/sump2;
    mom_ten_2(1,0) = mom_ten_2(0,1);
    mom_ten_2(2,0) = mom_ten_2(0,2);
    mom_ten_2(2,1) = mom_ten_2(1,2);
  }
  else {
    mom_ten_2(0,0) = 0;
    mom_ten_2(1,1) = 0;
    mom_ten_2(2,2) = 0;
    mom_ten_2(0,1) = mom_ten_2(1,0) = 0;
    mom_ten_2(0,2) = mom_ten_2(2,0) = 0;
    mom_ten_2(1,2) = mom_ten_2(2,1) = 0;
  }

  // Eigenvalues and Eigenvectors
  // ----------------------------
  TVectorD eigval(3);
  TMatrixD eigvec(3,3);
  eigvec = mom_ten_2.EigenVectors(eigval);

  _Aplanarity = 3./2.*(eigval(2));
  _Sphericity = 3./2.*(eigval(2)+eigval(1));
  Aplanarity->Fill(_Aplanarity);
//  _HiVar_Aplanarity_HistoSampler->Fill(_Aplanarity, fC, TOTEVENTS);
  Sphericity->Fill(_Sphericity);
//  _HiVar_Sphericity_HistoSampler->Fill(_Sphericity, fC, TOTEVENTS);
  _Y = sqrt(3.)/2.*(eigval(2)-eigval(1));
  Y_S->Fill(_Sphericity,_Y);

  // Evaluate the Thrust and Heavy and Light jet mass and Broadening
  EvalThrust(momJB, vec_jet.size());
  float jet_temp[3] = {0,0,0};
  float sum_jet_e = 0;
  float M1_vec[3] = {0,0,0};
  float M2_vec[3] = {0,0,0};
  float M1 = 0;
  float M2 = 0;
  float B1 = 0;
  float B2 = 0;
  for (unsigned int i=0; i<vec_jet.size(); ++i) {
    jet_temp[0] = momJB[i][0];
    jet_temp[1] = momJB[i][1];
    jet_temp[2] = momJB[i][2];
    _Thrust += fabs(_thrust[0]*jet_temp[0]+_thrust[1]*jet_temp[1]+_thrust[2]*jet_temp[2]);
    sum_jet_e += sqrt(jet_temp[0]*jet_temp[0]+jet_temp[1]*jet_temp[1]+jet_temp[2]*jet_temp[2]);
    // Evaluate Heavy and Light jet mass and Broadening
    if ((_thrust[0]*jet_temp[0]+_thrust[1]*jet_temp[1]+_thrust[2]*jet_temp[2])>0) {
      M1_vec[0] += jet_temp[0];
      M1_vec[1] += jet_temp[1];
      M1_vec[2] += jet_temp[2];
      B1 += pow(_thrust[1]*jet_temp[2]-_thrust[2]*jet_temp[1],2)+
	pow(_thrust[2]*jet_temp[0]-_thrust[0]*jet_temp[2],2)+
	pow(_thrust[0]*jet_temp[1]+_thrust[1]*jet_temp[0],2);
    }
    else {
      M2_vec[0] += jet_temp[0];
      M2_vec[1] += jet_temp[1];
      M2_vec[2] += jet_temp[2];
      B2 += pow(_thrust[1]*jet_temp[2]-_thrust[2]*jet_temp[1],2)+
	pow(_thrust[2]*jet_temp[0]-_thrust[0]*jet_temp[2],2)+
	pow(_thrust[0]*jet_temp[1]-_thrust[1]*jet_temp[0],2);
    }
  }
  if (sum_jet_e != 0) {
    _Thrust = _Thrust/sum_jet_e;

    M1 = (M1_vec[0]*M1_vec[0]+M1_vec[1]*M1_vec[1]+M1_vec[2]*M1_vec[2])/(sum_jet_e*sum_jet_e);
    M2 = (M2_vec[0]*M2_vec[0]+M2_vec[1]*M2_vec[1]+M2_vec[2]*M2_vec[2])/(sum_jet_e*sum_jet_e);
    B1 = B1/sum_jet_e;
    B2 = B2/sum_jet_e;
  }

  Thrust->Fill(_Thrust);
//  _HiVar_Thrust_HistoSampler->Fill(_Thrust, fC, TOTEVENTS);

  if (M1>M2) {
    _HjetMass = M1;
    _LjetMass = M2;
    HeavyJetMass->Fill(M1);
//    _HiVar_HeavyJetMass_HistoSampler->Fill(M1, fC, TOTEVENTS);
    LightJetMass->Fill(M2);
//    _HiVar_LightJetMass_HistoSampler->Fill(M2, fC, TOTEVENTS);
    JetMassDifference->Fill(M1-M2);
//    _HiVar_JetMassDifference_HistoSampler->Fill(M1-M2, fC, TOTEVENTS);
  }
  else {
    _HjetMass = M2;
    _LjetMass = M1;
    HeavyJetMass->Fill(M2);
//    _HiVar_HeavyJetMass_HistoSampler->Fill(M2, fC, TOTEVENTS);
    LightJetMass->Fill(M1);
//    _HiVar_LightJetMass_HistoSampler->Fill(M1, fC, TOTEVENTS);
    JetMassDifference->Fill(M2-M1);
//    _HiVar_JetMassDifference_HistoSampler->Fill(M2-M1, fC, TOTEVENTS);
  }

  if (B1>B2) {
    _MjetBroadening = B1;
    _mjetBroadening = B2;
    MajorBroadening->Fill(B1);
    MinorBroadening->Fill(B2);
    BroadeningDifference->Fill(B1-B2);
//    _HiVar_MajorBroadening_HistoSampler->Fill(B1, fC, TOTEVENTS);
//    _HiVar_MinorBroadening_HistoSampler->Fill(B2, fC, TOTEVENTS);
//    _HiVar_BroadeningDifference_HistoSampler->Fill(B1-B2, fC, TOTEVENTS);
  }
  else {
    _MjetBroadening = B2;
    _mjetBroadening = B1;
    MajorBroadening->Fill(B2);
    MinorBroadening->Fill(B1);
    BroadeningDifference->Fill(B2-B1);
//    _HiVar_MajorBroadening_HistoSampler->Fill(B2, fC, TOTEVENTS);
//    _HiVar_MinorBroadening_HistoSampler->Fill(B1, fC, TOTEVENTS);
//    _HiVar_BroadeningDifference_HistoSampler->Fill(B2-B1, fC, TOTEVENTS);
  }

  for (unsigned int i=0; i<vec_jet.size(); ++i) {
    delete [] momJB[i];
  }
  delete [] momJB;
}

void HiVariables::EvalThrust (float ** momJB, int vec_jet_size) {
  // Evaluate starting value for the iteration
  float N0[3] = {0,0,0};
  float N1[3] = {0,0,0};
  float N2[3] = {0,0,0};
  float N3[3] = {0,0,0};
  bool N0_good = false;
  bool N1_good = false;
  bool N2_good = false;
  bool N3_good = false;
  float n0[3] = {0,0,0};
  float n0_j[3] = {0,0,0};
  float jet_temp[3] = {0,0,0};
  float Psi_temp = 0;
  float diff = 2;
  float n0_mod = 0;
  float n0_j_mod = 0;
  // first iteration: n0 = versor of jet with k-highest E
  // loops 4 times and writes the N0, N1, N2 and N3.
  int kmax = 0;
  if (vec_jet_size < 4) kmax = vec_jet_size;
  else kmax = 4;
  for (int k=0; k<kmax; ++k) {
    n0[0] = momJB[k][0];
    n0[1] = momJB[k][1];
    n0[2] = momJB[k][2];
    n0_mod = sqrt(n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2]);
    if (n0_mod != 0) {
      n0[0] = n0[0]/n0_mod;
      n0[1] = n0[1]/n0_mod;
      n0[2] = n0[2]/n0_mod;
    }

    // iteration
//     for (int j=0; diff < 0.05; ++j) {
    int limit = 0;
    do {
      for (int i=0; i<vec_jet_size; ++i) {
	jet_temp[0] = momJB[i][0];
	jet_temp[1] = momJB[i][1];
	jet_temp[2] = momJB[i][2];
	Psi_temp = Psi(n0,jet_temp);
	n0_j[0] += Psi_temp*jet_temp[0];
	n0_j[1] += Psi_temp*jet_temp[1];
	n0_j[2] += Psi_temp*jet_temp[2];
      }
      n0_j_mod = sqrt(n0_j[0]*n0_j[0]+n0_j[1]*n0_j[1]+n0_j[2]*n0_j[2]);
      if (n0_j_mod != 0) {
	n0_j[0] = n0_j[0]/n0_j_mod;
	n0_j[1] = n0_j[1]/n0_j_mod;
	n0_j[2] = n0_j[2]/n0_j_mod;
      }
      diff = sqrt((n0_j[0]-n0[0])*(n0_j[0]-n0[0])+
		  (n0_j[1]-n0[1])*(n0_j[1]-n0[1])+
		  (n0_j[2]-n0[2])*(n0_j[2]-n0[2]));
      // store n0_j in n0 to prepare for the next iteration
      n0[0] = n0_j[0];
      n0[1] = n0_j[1];
      n0[2] = n0_j[2];
      n0_j[0] = 0;
      n0_j[1] = 0;
      n0_j[2] = 0;
      ++limit;
      if (limit == 10) {
	break;
	std::cout << "break" << std::endl; 
      }
    } while (diff > 0.05);
    // Store the thrust
    if (k == 0) {
      N0[0] = n0[0];
      N0[1] = n0[1];
      N0[2] = n0[2];
      N0_good = true;
    }
    if (k == 1) {
      N1[0] = n0[0];
      N1[1] = n0[1];
      N1[2] = n0[2];
      N1_good = true;
    }
    if (k == 2) {
      N2[0] = n0[0];
      N2[1] = n0[1];
      N2[2] = n0[2];
      N2_good = true;
    }
    if (k == 3) {
      N3[0] = n0[0];
      N3[1] = n0[1];
      N3[2] = n0[2];
      N3_good = true;
    }
  }
  // Conpare the thrust values
  float Diff[6];
  if (N0_good && N1_good)
    Diff[0] = sqrt((N0[0]-N1[0])*(N0[0]-N1[0])+
		   (N0[1]-N1[1])*(N0[1]-N1[1])+
		   (N0[2]-N1[2])*(N0[2]-N1[2]));
  else Diff[0] = 1;
  if (N0_good && N2_good)
    Diff[1] = sqrt((N0[0]-N2[0])*(N0[0]-N2[0])+
		   (N0[1]-N2[1])*(N0[1]-N2[1])+
		   (N0[2]-N2[2])*(N0[2]-N2[2]));
  else Diff[1] = 1;
  if (N0_good && N3_good)
    Diff[2] = sqrt((N0[0]-N3[0])*(N0[0]-N3[0])+
		   (N0[1]-N3[1])*(N0[1]-N3[1])+
		   (N0[2]-N3[2])*(N0[2]-N3[2]));
  else Diff[2] = 1;
  if (N1_good && N2_good)
    Diff[3] = sqrt((N1[0]-N2[0])*(N1[0]-N2[0])+
		   (N1[1]-N2[1])*(N1[1]-N2[1])+
		   (N1[2]-N2[2])*(N1[2]-N2[2]));
  else Diff[3] = 1;
  if (N1_good && N3_good)
    Diff[4] = sqrt((N1[0]-N3[0])*(N1[0]-N3[0])+
		   (N1[1]-N3[1])*(N1[1]-N3[1])+
		   (N1[2]-N3[2])*(N1[2]-N3[2]));
  else Diff[4] = 1;
  if (N2_good && N3_good)
    Diff[5] = sqrt((N2[0]-N3[0])*(N2[0]-N3[0])+
		   (N2[1]-N3[1])*(N2[1]-N3[1])+
		   (N2[2]-N3[2])*(N2[2]-N3[2]));
  else Diff[5] = 1;
  // check compatibility of the thrust axises
  int Diff_check[6];
  Diff_check[0] = (Diff[0]<0.1) ? 1 : 0;
  Diff_check[1] = (Diff[1]<0.1) ? 1 : 0;
  Diff_check[2] = (Diff[2]<0.1) ? 1 : 0;
  Diff_check[3] = (Diff[3]<0.1) ? 1 : 0;
  Diff_check[4] = (Diff[4]<0.1) ? 1 : 0;
  Diff_check[5] = (Diff[5]<0.1) ? 1 : 0;
  if (Diff_check[0]+Diff_check[1]+Diff_check[2] >= 2) {
    _thrust[0] = N0[0];
    _thrust[1] = N0[1];
    _thrust[2] = N0[2];
  }
  if (Diff_check[0]+Diff_check[3]+Diff_check[4] >= 2) {
    _thrust[0] = N1[0];
    _thrust[1] = N1[1];
    _thrust[2] = N1[2];
  }
  if (Diff_check[1]+Diff_check[3]+Diff_check[5] >= 2) {
    _thrust[0] = N2[0];
    _thrust[1] = N2[1];
    _thrust[2] = N2[2];
  }
  if (Diff_check[2]+Diff_check[4]+Diff_check[5] >= 2) {
    _thrust[0] = N3[0];
    _thrust[1] = N3[1];
    _thrust[2] = N3[2];
  } 
}

void HiVariables::Plot() {
  //  char HiVarName[30];
//  std::string OutputFileName( _histo_name );
//  std::string CanvasName( "Canvas_" );
  //  sprintf(HiVarName, "%s_HiVarPlots", _histo_name);
  //  TFile * outfile = new TFile (HiVarName,"RECREATE");
  //  TCanvas * HiVar = new TCanvas( HiVarName, HiVarName, 500, 500 );
//   OutputFileName += ".root";
//  CanvasName += _histo_name;
//  const char * canvasname = CanvasName.c_str();
//  TFile * outfile = new TFile ( OutputFileName.c_str(),"RECREATE" );
//   TCanvas * HiVar = new TCanvas( canvasname, canvasname, 500, 500 );
//   HiVar->SetFillColor(kNone);
//   HiVar->SetBorderSize(2);
//   HiVar->SetBorderMode(0);
//   HiVar->Divide(2,6);

//   HiVar->cd(1);
//   _HiVar_SumEt_HistoSampler->Draw();
  SumEtPlot()->Write();
//   HiVar->cd(2);
//  _HiVar_Centrality_HistoSampler->Draw();
  CentralityPlot()->Write();
//   HiVar->cd(3);
//   _HiVar_Aplanarity_HistoSampler->Draw();
  AplanarityPlot()->Write();
//   HiVar->cd(4);
//   _HiVar_Sphericity_HistoSampler->Draw();
  SphericityPlot()->Write();
//   HiVar->cd(5);
  Y_SPlot()->Write();
//   HiVar->cd(6);
//  _HiVar_Thrust_HistoSampler->Draw();
  ThrustPlot()->Write();
//   HiVar->cd(7);
//  _HiVar_HeavyJetMass_HistoSampler->Draw();
  HeavyJetMassPlot()->Write();
//   HiVar->cd(8);
//  _HiVar_LightJetMass_HistoSampler->Draw();
  LightJetMassPlot()->Write();
//   HiVar->cd(9);
//  _HiVar_JetMassDifference_HistoSampler->Draw();
  JetMassDifferencePlot()->Write();
//   HiVar->cd(10);
//  _HiVar_MajorBroadening_HistoSampler->Draw();
  MajorBroadeningPlot()->Write();
//   HiVar->cd(11);
//  _HiVar_MinorBroadening_HistoSampler->Draw();
  MinorBroadeningPlot()->Write();
//   HiVar->cd(12);
//  _HiVar_BroadeningDifference_HistoSampler->Draw();
  BroadeningDifferencePlot()->Write();

//   char HiVarFile[30];
//   sprintf(HiVarFile, "%s_HiVarPlots.eps", _histo_name);
//   HiVar->Print(HiVarFile);

//   SumEtPlot()->Write();
//   AplanarityPlot()->Write();
//   SphericityPlot()->Write();
//   Y_SPlot()->Write();
//   ThrustPlot()->Write();
//   HeavyJetMassPlot()->Write();
//   LightJetMassPlot()->Write();
//   JetMassDifferencePlot()->Write();
//   MajorBroadeningPlot()->Write();
//   MinorBroadeningPlot()->Write();
//   BroadeningDifferencePlot()->Write();
//  outfile->Close();
}

#endif // HIVARIABLES_C
