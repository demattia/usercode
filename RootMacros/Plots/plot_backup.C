/**
 * This macro draws two input histograms superimposed. If bands is passed as true,
 * it will draw the simulation histogram 3 times, emulating error bands.
 */

TCanvas * plotSimData (TH1F* histoDataIn, TString legendData, TH1F* histoSimulationIn, TString legendSimulation,
                       TString & canvasName, Float_t maximum = 0.15, TString xAxisTitle = "#eta",
                       TString yAxisTitle = "Number of Clusters", TString error = "", bool useLegend = true,
                       TString text = "", Float_t textX = 0.7, Float_t textY = 0.4, Float_t rebin = 0, bool bands = false ) {

  TH1F * histoData = (TH1F*)histoDataIn->Clone();
  TH1F * histoSimulation = (TH1F*)histoSimulationIn->Clone();

  int fillColor = 17;

  if ( rebin != 0 ) {
    histoSimulation->Rebin(rebin);
    histoData->Rebin(rebin);
  }

  TH1F * histoSimulationOrig = histoSimulation->Clone();
  histoSimulationOrig->Scale(1/(histoSimulation->Integral()));  
  if (bands) {
    histoSimulation->Sumw2();
  }
  histoData->Sumw2();
  histoData->Scale(1/(histoData->Integral()));
  histoSimulation->Scale(1/(histoSimulation->Integral()));

  // Create also the legend and add the histograms
  TLegend * legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );
  legend->AddEntry( histoData, legendData );
  legend->AddEntry( histoSimulation, legendSimulation, "F" );

  cout << "histoData = " << histoData << endl;
  cout << "histoSimulation = " << histoSimulation << endl;

  TCanvas * c = new TCanvas( canvasName, canvasName, 1000, 800 );
  c->Draw();

  histoSimulation->SetMinimum(0);
  histoSimulation->SetMaximum(maximum);

  //gStyle->SetErrorX(0.5);

  histoSimulation->SetFillColor(fillColor);
  //histoSimulation->SetLineWidth(2);
  histoSimulation->SetLineColor(fillColor);
  histoSimulation->SetXTitle(xAxisTitle);
  histoSimulation->SetYTitle(yAxisTitle);
  histoSimulation->SetTitleOffset(1.6,"Y");
  histoSimulation->SetTitle();

  histoData->SetMarkerStyle(21);

  if (bands) drawBands(histoSimulationOrig, histoSimulation, fillColor);
  else histoSimulation->Draw(error);

  histoData->Draw("sameep");

  legend->SetFillColor(kWhite);
  if (useLegend) legend->Draw("same");

  if ( text != "" ) {
    TPaveText * pt = new TPaveText(textX, textY, textX+0.2, textY+0.17, "NDC" ); // "NDC" option sets coords relative to pad dimensions
    pt->SetFillColor(0); // text is black on white
    pt->SetTextSize(0.08); 
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->AddText(text);
    pt->Draw("same");       //to draw your text object
  }

  return c;
};

drawBands(TH1F * histoOrig, TH1F * histo, int fillColor = 17) {

  TH1F * lowerErrorHisto = histo->Clone();
  lowerErrorHisto->Reset();
  TH1F * upperErrorHisto = histo->Clone();
  upperErrorHisto->Reset();

  int binsNumber = histo->GetNbinsX();
  // the 0 bin is the underflow and the binsNumber+1 is the overflow
  for( int iBin=1; iBin<binsNumber; ++iBin ) {
    float binValue = histo->GetBinContent(iBin);
    float binError = histo->GetBinError(iBin);
    lowerErrorHisto->SetBinContent(iBin, binValue-binError);
    upperErrorHisto->SetBinContent(iBin, binValue+binError);
  }

  lowerErrorHisto->SetFillColor(fillColor);
  upperErrorHisto->SetFillColor(fillColor+1);
  histo->SetFillColor(fillColor-1);

  upperErrorHisto->Draw();
  //histo->Draw("same");
  lowerErrorHisto->Draw("same");
  histoOrig->Draw("same");
}
