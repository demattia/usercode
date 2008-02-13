#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

int main()
{

  TFile * file = new TFile("file.root", "RECREATE");

  TH1F histo( "histo", "histo", 100, -3., 3. );
  histo.FillRandom( "gaus", 10000 );
  histo.SetDirectory(0);

  TCanvas canvas( "canvas", "canvas", 1000, 600 );
  canvas.cd();

  histo.Draw();

  canvas.Write();

  file->Write();

  printf("Howdy world!\n");
}
