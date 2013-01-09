// ROOT includes
#include <TROOT.h>
#include <TFile.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>

int main(int argc, char* argv[]) {

  char inputFileName[150];
  if ( argc < 2 ){
    std::cout << "missing argument: insert inputFile with list of root files" << std::endl; 
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  char Buffer[500];
  char MyRootFile[2000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);

  /// WEIGHTS AND SCALE HERE
  // float weights[] = {0.987,0.0100,0.0136,0.0151,0.0032,0.0858,0.0727};
  // int scaleStats = 2;
  float weights[] = {1.0,1.0};
  int scaleStats = 1;

  ///
  /// CATEGORIES HERE
  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);
  ///
  
  RooRealVar* MCweight = new RooRealVar("MCweight","Monte Carlo Weight",0.,1000.);

  TFile *theFile;
  TFile theOutput("totalDataSet.root","RECREATE");

  RooDataSet* thisData;
  RooDataSet* newData;
  RooDataSet* finalData;

  int nfiles=0;
  int nevents=0;

  while( !(inputFile->eof()) ){

    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))) { 
      sscanf(Buffer,"%s",MyRootFile);
      std::cout << "Merging " << MyRootFile << std::endl;

      theFile = TFile::Open(MyRootFile);
      thisData = (RooDataSet*)theFile->Get("data");

      const RooArgSet* thisRow = thisData->get(0); 
      RooArgSet* newRow = new RooArgSet(*thisRow);
      newRow->add(MCType);   newRow->add(*MCweight);

      newData = new RooDataSet("data","new data",*newRow,RooFit::WeightVar(*MCweight));

      for (Int_t iSamp = 0; iSamp < thisData->numEntries(); iSamp++)
	{
          nevents++;
          if (nevents%scaleStats == 0) {
	    thisRow = thisData->get(iSamp);
	    
	    RooCategory* myMatched = (RooCategory*)thisRow->find("matchType");
	    int isMatched = (int)(myMatched->getIndex());
	    int theMCType = 2;
	    if (isMatched && strstr(MyRootFile,"prompt")) theMCType = 0; 
	    if (isMatched && strstr(MyRootFile,"b")) theMCType = 1;
	    
	    MCType.setIndex(theMCType);
	    MCweight->setVal(weights[nfiles]);          
	    
	    RooArgSet* tempRow = new RooArgSet(*thisRow);
	    tempRow->add(MCType);   tempRow->add(*MCweight);
	    newData->add(*tempRow);
	   
	  }
	}

      if (nfiles == 0) {
	finalData = new RooDataSet(*newData);
      } else {
	finalData->append(*newData);
      }

      nfiles++;
    }
  }

  theOutput.cd();
  // finalData->setWeightVar(*MCweight);
  finalData->Write();
  theOutput.Close();
  inputFile->close();
  delete inputFile;
  
  return 0;

}
