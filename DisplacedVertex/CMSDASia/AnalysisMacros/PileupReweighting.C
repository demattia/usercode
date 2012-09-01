#define PileupReweighting_cxx

#include "PileupReweighting.h"


PileupReweighting::PileupReweighting( TString dataFileName):
dataFileName_(dataFileName)
{

  bool drawThings = false;


  // Get data distribution
  dataFile_ = new TFile( dataFileName_.Data() );
  Data_distr_ = (TH1F*) dataFile_->Get("pileup")->Clone();

  static const unsigned int nBinsX = Data_distr_->GetNbinsX();

//  std::cout << "N bins in data : " << Data_distr_->GetNbinsX() << endl;
//  std::cout << "Low edge : " << Data_distr_->GetBinLowEdge(1) << endl;
//  std::cout << "High edge : " << Data_distr_->GetBinLowEdge(nBinsX) + Data_distr_->GetBinWidth(nBinsX) << endl;
  weights_ = 0;

  // MC distribution, read from array below
  MC_distr_ = new TH1F("PU_generated","Generated pileup distribution (i.e., MC)",60,-0.5,59.5);


  // Summer 12 MC distribution

  Double_t Summer2012[60] = {
      2.344E-05,
      2.344E-05,
      2.344E-05,
      2.344E-05,
      4.687E-04,
      4.687E-04,
      7.032E-04,
      9.414E-04,
      1.234E-03,
      1.603E-03,
      2.464E-03,
      3.250E-03,
      5.021E-03,
      6.644E-03,
      8.502E-03,
      1.121E-02,
      1.518E-02,
      2.033E-02,
      2.608E-02,
      3.171E-02,
      3.667E-02,
      4.060E-02,
      4.338E-02,
      4.520E-02,
      4.641E-02,
      4.735E-02,
      4.816E-02,
      4.881E-02,
      4.917E-02,
      4.909E-02,
      4.842E-02,
      4.707E-02,
      4.501E-02,
      4.228E-02,
      3.896E-02,
      3.521E-02,
      3.118E-02,
      2.702E-02,
      2.287E-02,
      1.885E-02,
      1.508E-02,
      1.166E-02,
      8.673E-03,
      6.190E-03,
      4.222E-03,
      2.746E-03,
      1.698E-03,
      9.971E-04,
      5.549E-04,
      2.924E-04,
      1.457E-04,
      6.864E-05,
      3.054E-05,
      1.282E-05,
      5.081E-06,
      1.898E-06,
      6.688E-07,
      2.221E-07,
      6.947E-08,
      2.047E-08
  };

  // Fill array with nBinsX
  Double_t TrueData[nBinsX];

  for (unsigned int i=1;i<nBinsX+1;i++) {
    TrueData[i-1]=Data_distr_->GetBinContent(i);
//    std::cout << "Bin : " << i << " Content : " << TrueData[i-1] << std::endl;
  }

  Data_distr_->Reset();
  Data_distr_ = new TH1F("PU_intended","Intended pileup distribution (i.e., Data)",nBinsX,-0.5,59.5);

  //Fill histograms with one generated, one desired distribution

  for (unsigned int i=1;i<61;i++) {
    MC_distr_->SetBinContent(i,Summer2012[i-1]);
  }
  for (unsigned int i=1;i<nBinsX+1;i++) {
    Data_distr_->SetBinContent(i,TrueData[i-1]);
//    std::cout << "Setting bin content : " << i << " as " << TrueData[i-1] << endl;
  }

//  double integral_data = Data_distr_->Integral();
//  double integral_mc = MC_distr_->Integral();
//
//  std::cout << "Integral data : " << integral_data << " MC : " << integral_mc << std::endl;

  TCanvas *can1 = 0;

  if ( drawThings ) {
    can1 = new TCanvas("can1,","can1",1200,800);
    can1->Divide(2,2);

    can1->cd(1);
    MC_distr_->SetLineColor(kRed);
    Data_distr_->SetLineColor(kBlue);

    Data_distr_->DrawClone();
    MC_distr_->DrawClone("Same");

    can1->cd(2);
  }
  // Normalize histograms
  Data_distr_->Scale(nBinsX/Data_distr_->Integral());
  MC_distr_->Scale(MC_distr_->GetNbinsX()/MC_distr_->Integral());

//  integral_data = Data_distr_->Integral();
//  integral_mc = MC_distr_->Integral();
//
//  std::cout << "Integral data : " << integral_data << " MC : " << integral_mc << std::endl;

  if ( drawThings ) {
    Data_distr_->Draw();
    MC_distr_->Draw("Same");

    can1->cd(3);
  }

  // Calculate weights
  weights_ = (TH1F*) Data_distr_->Clone();
  for (unsigned int i=1;i<nBinsX+1;i++) {

    double dataBinContent = weights_->GetBinContent( i );

    // Find value of x axis in this bin
    double nPV_true_data = weights_->GetXaxis()->GetBinCenter( i );

    // Find corresponding bin in MC distribution
    int binInMC = MC_distr_->GetXaxis()->FindBin( nPV_true_data );

    // Bin content in MC
    double mcBinContent = MC_distr_->GetBinContent( binInMC );

    // Set weight Content
    weights_->SetBinContent( i, dataBinContent / mcBinContent );
  }
  //  weights_->Divide( MC_distr_ );

  if ( drawThings ) {
    weights_->SetLineColor(kGreen);
    weights_->Draw();
  }

}

PileupReweighting::~PileupReweighting() {
  //  dataFile_->Close();
}



double PileupReweighting::weight( float nPV ) {
  // Which bin?
  int bin = weights_->GetXaxis()->FindBin( nPV );
  // Return the weight for this value of nPV
  //  double bin = weights_->GetBinContent( bin );
  return weights_->GetBinContent( bin );
}

