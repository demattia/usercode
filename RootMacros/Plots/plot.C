/**
 * This macro draws two input histograms superimposed. It does not
 * draw errors on the data histogram. To have errors, use plot_backup.C
 */

TCanvas * plot (TH1F* histoDataIn, TString legendData, TH1F* histoSimulationIn, TString legendSimulation,
                TString & canvasName, Float_t maximum = 0.15, TString xAxisTitle = "#eta",
                TString yAxisTitle = "Number of Clusters", TString error = "", bool useLegend = true,
                TString text = "", Float_t textX = 0.7, Float_t textY = 0.4, Float_t rebin = 0 ) {

  TH1F * histoData = (TH1F*)histoDataIn->Clone();
  TH1F * histoSimulation = (TH1F*)histoSimulationIn->Clone();

  // histoData->Sumw2();
  histoData->Scale(1/(histoData->Integral()));
  histoSimulation->Scale(1/(histoSimulation->Integral()));

  // Create also the legend and add the histograms
  TLegend * legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );
  legend->AddEntry( histoData, xAxisTitle );
  legend->AddEntry( histoSimulation, yAxisTitle, "F" );

  cout << "histoData = " << histoData << endl;
  cout << "histoSimulation = " << histoSimulation << endl;

  TCanvas * c = new TCanvas( canvasName, canvasName, 1000, 800 );
  c->Draw();

  histoSimulation->SetMaximum(maximum);

  histoSimulation->SetFillColor(kRed);
  // histoSimulation->SetLineWidth(0);
  histoSimulation->SetLineColor(kRed);
  histoSimulation->SetXTitle(xAxisTitle);
  histoSimulation->SetYTitle(yAxisTitle);
  histoSimulation->SetTitleOffset(1.6,"Y");
  histoSimulation->SetTitle();

  histoData->SetLineStyle(1);
  histoData->SetLineWidth(2.5);

  histoSimulation->Draw();
  histoData->Draw("same");

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
