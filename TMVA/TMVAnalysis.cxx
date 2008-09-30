// @(#)root/tmva $Id: TMVAnalysis.cxx,v 1.106 2008/08/06 20:41:56 andreas.hoecker Exp $
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAnalysis                                                        *
 *                                                                                *
 * This executable provides examples for the training and testing of the          *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans.        *
 *                                                                                *
 * Compile and run the example with the following commands                        *
 *                                                                                *
 *    make                                                                        *
 *    ./TMVAnalysis <Methods>                                                     *
 *                                                                                *
 * where: <Methods> = "method1 method2"                                           *
 *        are the TMVA classifier names                                           *
 *                                                                                *
 * example:                                                                       *
 *    ./TMVAnalysis Fisher LikelihoodPCA BDT                                      *
 *                                                                                *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <../macros/macro.C>), which can be conveniently    *
 * invoked through a GUI launched by the command                                  *
 *                                                                                *
 *    root -l ../macros/TMVAGui.C                                                 *
 **********************************************************************************/

#include <iostream> // Stream declarations
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/MethodTMlpANN.h"
#include "TMVA/Tools.h"

// read input data file with ascii format (otherwise ROOT) ?
Bool_t ReadDataFromAsciiIFormat = kFALSE;

int main( int argc, char** argv ) 
{
   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<TString,Bool_t> Use;

   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // ---
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   // ---
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDERSkNN"]        = 0; // depreciated until further notice
   Use["KNN"]             = 0;
   // ---
   Use["HMatrix"]         = 0;
   Use["Fisher"]          = 0;
   // ---
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   // ---
   Use["MLP"]             = 0; // this is the recommended ANN
   Use["CFMlpANN"]        = 0; 
   Use["TMlpANN"]         = 0; 
   // ---
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;
   // ---
   Use["BDT"]             = 1;
   Use["BDTD"]            = 0;
   // ---
   Use["RuleFit"]         = 0;
   Use["RuleFitJF"]       = 0;
   // ---------------------------------------------------------------

   std::cout << "Start TMVAnalysis" << std::endl
             << "=================" << std::endl
             << std::endl;
   std::cout << "Running the following methods" << std::endl;
   if (argc>1) {
      for (std::map<TString,Bool_t>::iterator it = Use.begin(); it != Use.end(); it++) (*it).second=kFALSE;
   }
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if (Use.find(regMethod) == Use.end()) {
         std::cout << "Method " << regMethod << " not known in TMVA under this name. Please try one of:" << std::endl;
         for (std::map<TString,Bool_t>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << (*it).first << " ";
         std::cout << std::endl;
         return 1;
      }
      Use[regMethod] = kTRUE;
   }

   // Create a new root output file.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAnalysis", outputFile, 
                                               Form("!V:!Silent:%sColor", gROOT->IsBatch()?"!":"") );

   // If you wish to modify default settings 
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
//    factory->AddVariable( "var1+var2", 'F' );
//    factory->AddVariable( "var1-var2", 'F' );
//    factory->AddVariable( "var3",      'F' );
//    factory->AddVariable( "var4",      'F' );

   // read training and test data
   if (ReadDataFromAsciiIFormat) {
      // load the signal and background event samples from ascii files
      // format in file must be:
      // var1/F:var2/F:var3/F:var4/F
      // 0.04551   0.59923   0.32400   -0.19170
      // ...

      TString datFileS = "tmva_example_sig.dat";
      TString datFileB = "tmva_example_bkg.dat";
      factory->SetInputTrees( datFileS, datFileB );
   }
   else {
      // load the signal and background event samples from ROOT trees
      // TFile *input(0);
      // TString fname = "./tmva_example.root";
      // if (!gSystem->AccessPathName( fname )) {
      //    // first we try to find tmva_example.root in the local directory
      //    std::cout << "--- TMVAnalysis    : Accessing " << fname << std::endl;
      //    input = TFile::Open( fname );
      // }
//       else { 
//          // second we try accessing the file via the web from
//          // http://root.cern.ch/files/tmva_example.root
//          std::cout << "--- TMVAnalysis    : Accessing tmva_example.root file from http://root.cern.ch/files" << std::endl;
//          std::cout << "--- TMVAnalysis    : For faster startup you may consider downloading it into you local directory" << std::endl;
//          input = TFile::Open( "http://root.cern.ch/files/tmva_example.root" );
//       }

//       if (!input) {
//          std::cout << "ERROR: could not open data file" << std::endl;
//          exit(1);
//       }




     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------
     TFile * inputSig              = TFile::Open( "/data/demattia/test/tmva_2tagsTTH_120.root" );
     TFile * inputBkg_qcd_30_50    = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_30-50.root" );
     TFile * inputBkg_qcd_50_80    = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_50-80.root" );
     TFile * inputBkg_qcd_80_120   = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_80-120.root" );
     TFile * inputBkg_qcd_120_170  = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_120-170.root" );
     TFile * inputBkg_qcd_170_230  = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_170-230.root" );
     TFile * inputBkg_qcd_230_300  = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_230-300.root" );
     TFile * inputBkg_qcd_300_380  = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_300-380.root" );
     TFile * inputBkg_qcd_380_incl = TFile::Open( "/data/demattia/test/tmva_2tagsQCD_380-incl.root" );
     TFile * inputBkg_tt_0jets     = TFile::Open( "/data/demattia/test/tmva_2tagsTT_0JETS.root" );
     TFile * inputBkg_tt_1jets     = TFile::Open( "/data/demattia/test/tmva_2tagsTT_1JETS.root" );
     TFile * inputBkg_tt_2jets     = TFile::Open( "/data/demattia/test/tmva_2tagsTT_2JETS.root" );
     TFile * inputBkg_tt_3jets     = TFile::Open( "/data/demattia/test/tmva_2tagsTT_3JETS.root" );
     TFile * inputBkg_tt_4jets     = TFile::Open( "/data/demattia/test/tmva_2tagsTT_4JETS.root" );
     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------
     
     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------
     TTree *signal                  = (TTree*)inputSig->Get("tree__2tagsTTH_120");
     TTree *background_qcd_30_50    = (TTree*)inputBkg_qcd_30_50->Get("tree__2tagsQCD_30-50");
     TTree *background_qcd_50_80    = (TTree*)inputBkg_qcd_50_80->Get("tree__2tagsQCD_50-80");
     TTree *background_qcd_80_120   = (TTree*)inputBkg_qcd_80_120->Get("tree__2tagsQCD_80-120");
     TTree *background_qcd_120_170  = (TTree*)inputBkg_qcd_120_170->Get("tree__2tagsQCD_120-170");
     TTree *background_qcd_170_230  = (TTree*)inputBkg_qcd_170_230->Get("tree__2tagsQCD_170-230");
     TTree *background_qcd_230_300  = (TTree*)inputBkg_qcd_230_300->Get("tree__2tagsQCD_230-300");
     TTree *background_qcd_300_380  = (TTree*)inputBkg_qcd_300_380->Get("tree__2tagsQCD_300-380");
     TTree *background_qcd_380_incl = (TTree*)inputBkg_qcd_380_incl->Get("tree__2tagsQCD_380-incl");
     TTree *background_tt_0jets     = (TTree*)inputBkg_tt_0jets->Get("tree__2tagsTT_0JETS");
     TTree *background_tt_1jets     = (TTree*)inputBkg_tt_1jets->Get("tree__2tagsTT_1JETS");
     TTree *background_tt_2jets     = (TTree*)inputBkg_tt_2jets->Get("tree__2tagsTT_2JETS");
     TTree *background_tt_3jets     = (TTree*)inputBkg_tt_3jets->Get("tree__2tagsTT_3JETS");
     TTree *background_tt_4jets     = (TTree*)inputBkg_tt_4jets->Get("tree__2tagsTT_4JETS");
     
     Double_t signalWeight                  = 0.667 / 76516;
     Double_t backgroundWeight_qcd_30_50    = 163000000. / 80709;
     Double_t backgroundWeight_qcd_50_80    = 21600000. / 79347;
     Double_t backgroundWeight_qcd_80_120   = 3080000. / 86011;
     Double_t backgroundWeight_qcd_120_170  = 494000. / 81308;
     Double_t backgroundWeight_qcd_170_230  = 101000. / 83469;
     Double_t backgroundWeight_qcd_230_300  = 24500. / 79026;
     Double_t backgroundWeight_qcd_300_380  = 6240. / 81077;
     Double_t backgroundWeight_qcd_380_incl = 2860. / 81478;
     Double_t backgroundWeight_tt_0jets     = 619. / 7388;
     Double_t backgroundWeight_tt_1jets     = 176. / 18171;
     Double_t backgroundWeight_tt_2jets     = 34. / 2000;
     Double_t backgroundWeight_tt_3jets     = 6. / 3000;
     Double_t backgroundWeight_tt_4jets     = 1.5 / 535;
      
     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------
     factory->AddSignalTree    ( signal,                  signalWeight );
     factory->AddBackgroundTree( background_qcd_30_50,    backgroundWeight_qcd_30_50 );
     factory->AddBackgroundTree( background_qcd_50_80,    backgroundWeight_qcd_50_80 );
     factory->AddBackgroundTree( background_qcd_80_120,   backgroundWeight_qcd_80_120 );
     factory->AddBackgroundTree( background_qcd_120_170,  backgroundWeight_qcd_120_170 );
     factory->AddBackgroundTree( background_qcd_170_230,  backgroundWeight_qcd_170_230 );
     factory->AddBackgroundTree( background_qcd_230_300,  backgroundWeight_qcd_230_300 );
     factory->AddBackgroundTree( background_qcd_300_380,  backgroundWeight_qcd_300_380 );
     factory->AddBackgroundTree( background_qcd_380_incl, backgroundWeight_qcd_380_incl );
     factory->AddBackgroundTree( background_tt_0jets,     backgroundWeight_tt_0jets );
     // This file seems to give problems
     factory->AddBackgroundTree( background_tt_1jets,     backgroundWeight_tt_1jets );
     //
     factory->AddBackgroundTree( background_tt_2jets,     backgroundWeight_tt_2jets );
     factory->AddBackgroundTree( background_tt_3jets,     backgroundWeight_tt_3jets );
     factory->AddBackgroundTree( background_tt_4jets,     backgroundWeight_tt_4jets );
     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------








//       TTree *signal     = (TTree*)input->Get("TreeS");
//       TTree *background = (TTree*)input->Get("TreeB");

//       // global event weights per tree (see below for setting event-wise weights)
//       Double_t signalWeight     = 1.0;
//       Double_t backgroundWeight = 1.0;

      // ====== register trees ====================================================
      //
      // the following method is the prefered one:
      // you can add an arbitrary number of signal or background trees
//      factory->AddSignalTree    ( signal,     signalWeight     );
//      factory->AddBackgroundTree( background, backgroundWeight );

      // To give different trees for training and testing, do as follows:
      //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
      //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

      // Use the following code instead of the above two or four lines to add signal and background 
      // training and test events "by hand"
      // NOTE that in this case one should not give expressions (such as "var1+var2") in the input 
      //      variable definition, but simply compute the expression before adding the event
      // 
      //    // --- begin ----------------------------------------------------------
      //    std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
      //    Float_t  treevars[4];
      //    for (Int_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
      //    for (Int_t i=0; i<signal->GetEntries(); i++) {
      //       signal->GetEntry(i);
      //       for (Int_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
      //       // add training and test events; here: first half is training, second is testing
      //       // note that the weight can also be event-wise	
      //       if (i < signal->GetEntries()/2) factory->AddSignalTrainingEvent( vars, signalWeight ); 
      //       else                            factory->AddSignalTestEvent    ( vars, signalWeight ); 
      //    }
      //
      //    for (Int_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
      //    for (Int_t i=0; i<background->GetEntries(); i++) {
      //       background->GetEntry(i); 
      //       for (Int_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
      //       // add training and test events; here: first half is training, second is testing
      //       // note that the weight can also be event-wise	
      //       if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight ); 
      //       else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight ); 
      //    }
      //    // --- end ------------------------------------------------------------
      //
      // ====== end of register trees ==============================================
   }
   






   // Adding variables

   factory->AddVariable("bTagTkInvMass", 'D');
   factory->AddVariable("chi2ofMasses", 'D');
   factory->AddVariable("deltaEtaHadronicTopHiggs", 'D');
   factory->AddVariable("deltaPhiMEtNthLeadingJet_1", 'D');
   factory->AddVariable("deltaPhiMEtNthLeadingJet_2", 'D');
   factory->AddVariable("deltaPhiMEtNthLeadingJet_3", 'D');
   factory->AddVariable("firstNjetsCentrality_1", 'D');
   factory->AddVariable("firstNjetsCentrality_2", 'D');
   factory->AddVariable("firstNjetsMass_1", 'D');
   factory->AddVariable("firstNjetsMass_2", 'D');
   factory->AddVariable("goodHt", 'D');
   factory->AddVariable("hadronicTopMass", 'D');
   factory->AddVariable("hadronicTopPlusHiggsMass", 'D');
   factory->AddVariable("hadronicTopProjectionAlongHiggsDirection", 'D');
   factory->AddVariable("hadronicWmass", 'D');
   factory->AddVariable("higgsMass", 'D');
   factory->AddVariable("mEtSig", 'D');
   // factory->AddVariable("remainingJetsMass", 'D');
   factory->AddVariable("sixthJetEt", 'D');
   factory->AddVariable("sumHighEffDiscriminantFirst4Jets", 'D');
   factory->AddVariable("sumHighEffDiscriminantFirst6Jets", 'D');


   // This would set individual event weights (the variables defined in the 
   // expression need to exist in the original TTree)
   // factory->SetWeightExpression("weight1*weight2");
   factory->SetSignalWeightExpression("__WEIGHT__");
   factory->SetBackgroundWeightExpression("__WEIGHT__");





   // This would set individual event weights (the variables defined in the 
   // expression need to exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // tell the factory to use all remaining events in the trees after training for testing:
   // factory->PrepareTrainingAndTestTree( mycuts, mycutb, 
   //                                      "NSigTrain=3000:NBkgTrain=3000:SplitMode=Random:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used for training, and 
   // the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  
   // To also specify the number of testing events, use:
   // factory->PrepareTrainingAndTestTree( mycuts, 
   //                                      "NSigTrain=1500:NBkgTrain=10000:NBkgTest=10000:SplitMode=Random:!V" );  

   factory->PrepareTrainingAndTestTree( mycuts, 
                                        "SplitMode=Random:!V" );  

   // ---- Book MVA methods
   //
   // please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc.
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"]) 
     factory->BookMethod( TMVA::Types::kCuts, "Cuts", 
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"]) 
     factory->BookMethod( TMVA::Types::kCuts, "CutsD", 
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"]) 
     factory->BookMethod( TMVA::Types::kCuts, "CutsPCA", 
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
     factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                         "H:!V:FitMethod=GA:EffSel:Steps=30:Cycles=3:PopSize=100:SC_steps=10:SC_rate=5:SC_factor=0.95:VarProp=FSmart" );
   
   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemperature=IncreasingAdaptive:InitialTemperature=1e+6:MinTemperature=1e-6:Eps=1e-10:UseDefaultScale" );
   
   // Likelihood
   if (Use["Likelihood"]) 
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", 
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=10:NSmoothBkg[0]=10:NSmoothBkg[1]=10:NSmooth=10:NAvEvtPerBin=50" ); 

   // test the decorrelated likelihood
   if (Use["LikelihoodD"]) 
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD", 
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=10:NSmoothBkg[0]=10:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ); 

   if (Use["LikelihoodPCA"]) 
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA", 
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=10:NSmoothBkg[0]=10:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 
 
   // test the new kernel density estimator
   if (Use["LikelihoodKDE"]) 
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE", 
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // test the mixed splines and kernel density estimator (depending on which variable)
   if (Use["LikelihoodMIX"]) 
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX", 
                           "!H:!V:!TransformOutput:PDFInterpol[0]=KDE:PDFInterpol[1]=KDE:PDFInterpol[2]=Spline2:PDFInterpol[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // PDE - RS method
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS", 
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
   // And the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   

   if (Use["PDERSkNN"]) // depreciated until further notice
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSkNN", 
                           "!H:!V:VolumeRangeMode=kNN:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"]) 
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD", 
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"]) 
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA", 
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
  
   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN", 
                           "nkNN=400:TreeOptDepth=6:ScaleFrac=0.8:!UseKernel:!Trim" );  

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" ); 

   // Fisher discriminant
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", 
                           "H:!V:!Normalise:CreateMVAPdfs:Fisher:NbinsMVAPdf=50:NsmoothMVAPdf=1" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
   
   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options)
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=3:Steps=20:Trim=True:SaveBestGen=0" );

   if (Use["FDA_MT"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"]) 
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:!Normalise:NeuronType=tanh:NCycles=200:HiddenLayers=N+1,N:TestRate=5" );
      // factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:!Normalise:NeuronType=tanh:NCycles=200:HiddenLayers=6:TestRate=5" );

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=500:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  
  
   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...
  
   // Support Vector Machines using three different Kernel types (Gauss, polynomial and linear)
   if (Use["SVM_Gauss"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM_Gauss", "Sigma=2:C=1:Tol=0.001:Kernel=Gauss" );

   if (Use["SVM_Poly"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM_Poly", "Order=4:Theta=1:C=0.1:Tol=0.001:Kernel=Polynomial" );

   if (Use["SVM_Lin"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM_Lin", "!H:!V:Kernel=Linear:C=1:Tol=0.001" );  

   // Boosted Decision Trees (second one with decorrelation)
   if (Use["BDT"])
      factory->BookMethod( TMVA::Types::kBDT, "BDT", 
                           "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=1.5" );
   if (Use["BDTD"])
      factory->BookMethod( TMVA::Types::kBDT, "BDTD", 
                           "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=1.5:VarTransform=Decorrelate" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // Friedman's RuleFit method, implementation by J. Friedman
   if (Use["RuleFitJF"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFitJF",
                           "!V:RuleFitModule=RFFriedman:Model=ModRuleLinear:GDStep=0.01:GDNSteps=10000:GDErrScale=1.1:RFNendnodes=4" );
 
   // --------------------------------------------------------------------------------------------------

   // For an example how to use the ROOT plugin mechanism, please see TMVA/macros/TMVAnalysis.C
   // ...
   
   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAnalysis is done!" << std::endl;      

   std::cout << std::endl;
   std::cout << "==> Too view the results, launch the GUI: \"root -l ../macros/TMVAGui.C\"" << std::endl;
   std::cout << std::endl;

   // Clean up
   delete factory;
}

