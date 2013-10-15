/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include "setdirs.h"

#include "TMVAGui.C"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

//typedef float mytype;

using namespace TMVA;

Int_t TMVAClassificationApplication_main( const TString inputFileName = 
					  //"test_t3.root",
					  //"trees_main/data_afterCuts_0_0.root", 
					  "/home/data/MainAnalysis/data_afterCuts_0_0.root",
					 //"rootfiles/Barrel_preselection.root", 
					  const TString outputFileName = "test_mvaApp.root", 
					  const TString weightfile = //"weights_main/TMVA-2-Events0_BDT.weights.xml", 
					  "/home/nuno/TMVA-2-Events0_BDT.weights.xml",
					  const float cutValue=0.35, TString myMethodList = "BDT", TString region = "barrel", 
					 //TString tree_name = "probe_tree" 
					 TString tree_name = "events", TString anaType = "main" 
					 ) {   
  
  bool barrel = false;
  //  if(outputFileName.Contains("arrel"))
  if(region=="barrel")
    barrel = true;
  bool main = false;
  if(anaType=="main")
    main = true;

  /*
  bool mainData = false;
  if(inputFileName.Contains("trees_main"))
    mainData = true;
  cout << "mainData:" << mainData << endl;
  */

  //const TString & weightfile = "weights/barrelWeights/TMVAClassification_BDT.weights.xml",
  //weightfile = "weights_main/TMVA-0-Events0_BDT.weights.xml"

  //theTree->Show();

  //TString treefilename = outputFileName;
  //treefilename.ReplaceAll(".root","_tree.root");
  //cout << treefilename << endl;


   //  return 1;


  //TCanvas c;
  //theTree->Draw("m");
  //c.SaveAs("plots/tmp.gif");
  
  // this redefinition is necessary as the different trees are formed with floats or doubles
  //  reading a float tree with a double variable (or vice-versa) gives garbish!!
  // the workaround works with interpreter root, but does not compile
  // ... however there is a limitation of the TMVA::Reader which accpepts only floats, not doubles!!!!
  /*
  if(!mainData) {
    typedef float mytype;
  } else {
    typedef double mytype;
  }
  */



  /*
  /// test of double vs float -- use of typedef
  mytype m;
  theTree->SetBranchAddress( "m",                            &m );
  for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    cout << "mass:"  << m << endl;
  }
  return;
  */





  //  cout << "aaaaaa\n" << endl;

#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   bool useBDT = myMethodList=="BDT"?1:0;
   bool useMLP = myMethodList=="MLP"?1:0;
   bool useCutsSA = myMethodList=="CutsSA"?1:0;

   cout << "running TMVAClassificationApplication_main for method:" << myMethodList << " BDT:" << useBDT << " MLP:" << useMLP << " Cuts:" << useCutsSA << " withy cut value:" << cutValue << " for " << region << endl;


   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;//useCutsSA;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0;//useMLP; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1;//useBDT; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return 10;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
   // TMVA::ReaderMuonID *readerMuonID = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used

   //mytype
     //  float 
     //  double 
     //m, fls3d, alpha, pvips, iso, m1iso, m2iso, chi2dof, eta, pt, maxdoca, docatrk, pvip;
   Float_t m_, fls3d_, alpha_, pvips_, iso_, m1iso_, m2iso_, chi2dof_, eta_, pt_, maxdoca_, docatrk_, pvip_; // MVA reader requires float!!

     //float m_(m), fls3d_(fls3d), alpha_(alpha), pvips_(pvips), iso_(iso), m1iso_(m1iso), m2iso_(m2iso), chi2dof_(chi2dof), eta_(eta), pt_(pt), maxdoca_(maxdoca), docatrk_(docatrk), pvip_(pvip);

   if(barrel) {
     reader->AddVariable( "fls3d",                     (&fls3d_) );
     reader->AddVariable( "alpha",     		       (&alpha_) );
     reader->AddVariable( "pvips",     		       (&pvips_) );
     reader->AddVariable( "iso",       		       (&iso_) );
     reader->AddVariable( "m1iso",     		       (&m1iso_) );
     reader->AddVariable( "m2iso",     		       (&m2iso_) );
     reader->AddVariable( "chi2dof",   		       (&chi2dof_) );
     reader->AddVariable( "eta",       		       (&eta_) );
     reader->AddVariable( "maxdoca",                   (&maxdoca_) );
     reader->AddVariable( "docatrk",   		       (&docatrk_) );
   } else {
     reader->AddVariable( "fls3d",                     (&fls3d_) );
     reader->AddVariable( "alpha",     		       (&alpha_) );
     reader->AddVariable( "pvips",     		       (&pvips_) );
     reader->AddVariable( "iso",       		       (&iso_) );
     reader->AddVariable( "m1iso",     		       (&m1iso_) );
     reader->AddVariable( "m2iso",     		       (&m2iso_) );
     reader->AddVariable( "chi2dof",   		       (&chi2dof_) );
     reader->AddVariable( "pt",                        (&pt_) );
     reader->AddVariable( "pvip",      		       (&pvip_) );
     reader->AddVariable( "docatrk",   		       (&docatrk_) );
   }

   // Spectator variables declared in the training have to be added to the reader, too
   // Float_t spec1,spec2;
   // reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   // reader->AddSpectator( "spec2 := var1*3",   &spec2 );

   reader->AddSpectator( "m", (&m_) );
   //reader->AddSpectator( "m", static_cast<float*> (&m) );


   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }

   // --- Book the MVA methods

   // TString dir    = "endcapsWeights/";
   // TString dir    = "barrelWeights/";
   TString prefix = "TMVAClassification";

   // TString dirMuonID    = "tmvaMuonID/";
   // TString prefixMuonID = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         // TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         //TString weightfile = weightDir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 

         // TString weightfileMuonID = dirMuonID + prefixMuonID + TString("_") + TString(it->first) + TString(".weights.xml");
         // readerMuonID->BookMVA( methodName, weightfileMuonID ); 
      }
   }
   
   /// END READER

   /// START INPUT TREE

   TFile *inputS = TFile::Open( inputFileName, "READ" );
   std::cout << "--- TMVAClassificationApplication       : Using input file: " << inputS->GetName() << std::endl;
   TTree *theTree = (TTree*)inputS->Get(tree_name); 


   mytype m, fls3d, alpha, pvips, iso, m1iso, m2iso, chi2dof, eta, pt, maxdoca, docatrk, pvip;
 
   //theTree->SetBranchAddress( "m",                            static_cast<float*> (&m) );
   theTree->SetBranchAddress( "m",                            &m );

   //   if(barrel) {
     theTree->SetBranchAddress( "fls3d",                           &fls3d );
     theTree->SetBranchAddress( "alpha",                           &alpha );
     theTree->SetBranchAddress( "pvips",                           &pvips );
     theTree->SetBranchAddress( "iso",                             &iso );
     theTree->SetBranchAddress( "m1iso",                           &m1iso );
     theTree->SetBranchAddress( "m2iso",                           &m2iso );
     theTree->SetBranchAddress( "chi2dof",                         &chi2dof );
     theTree->SetBranchAddress( "eta",                             &eta );
     theTree->SetBranchAddress( "maxdoca",                         &maxdoca );
     theTree->SetBranchAddress( "docatrk",                         &docatrk );
//   } else {
//     theTree->SetBranchAddress( "fls3d",                           &fls3d );
//     theTree->SetBranchAddress( "alpha",                           &alpha );
//     theTree->SetBranchAddress( "pvips",                           &pvips );
//     theTree->SetBranchAddress( "iso",                             &iso );
//     theTree->SetBranchAddress( "m1iso",                           &m1iso );
//     theTree->SetBranchAddress( "m2iso",                          &m2iso );
//     theTree->SetBranchAddress( "chi2dof",                         &chi2dof );
     theTree->SetBranchAddress( "pt",                              &pt );
     theTree->SetBranchAddress( "pvip",                            &pvip );
//     theTree->SetBranchAddress( "docatrk",                         &docatrk );
//   }
  

   /// should try to confirm preselection cuts, to avoid inconsistencies

   // read extra variables

     mytype me, pvw8, m1pt, m2pt, m1eta,m2eta, fl3d, pvlip, pvlips;
     UInt_t evtnum_xcheck; mytype closetrk_xcheck;  //xcheck
	 Long64_t evtnum; Int_t closetrk;
     //lxysig, flsxy, .. these names differ for main and xcheck

     if (main) {
       theTree->SetBranchAddress( "evt",                            &evtnum ); //main
       theTree->SetBranchAddress( "pvlips",                              &pvlips ); //main
       theTree->SetBranchAddress( "closetrk",                              &closetrk ); //main
     }
     else {
       theTree->SetBranchAddress( "event",                            &evtnum_xcheck ); //xcheck
       theTree->SetBranchAddress( "pvlipErr",                              &pvlips ); //xcheck
       theTree->SetBranchAddress( "closetrk",                              &closetrk_xcheck ); //xcheck
     }     
     theTree->SetBranchAddress( "me",                            &me );
     theTree->SetBranchAddress( "pvw8",                              &pvw8 );
     theTree->SetBranchAddress( "m1pt",                              &m1pt );
     theTree->SetBranchAddress( "m2pt",                              &m2pt );
     theTree->SetBranchAddress( "m1eta",                              &m1eta );
     theTree->SetBranchAddress( "m2eta",                              &m2eta );
     theTree->SetBranchAddress( "fl3d",                              &fl3d );
     theTree->SetBranchAddress( "pvlip",                              &pvlip );
      
   /// END INPUT TREE

   /// START OUTPUT TREE

   TFile *resTreeFile = new TFile( outputFileName,"RECREATE" );
   TTree *resTree = (TTree*) theTree->CloneTree(0);
   //TTree* resTree = (TTree*)theTree->CopyTree("");
   //if (resTree == 0) return 0;
   //resTree->Reset();
   
   Float_t bdt_v;
   TBranch* bdtOutput = resTree->Branch("bdt_v", &bdt_v, "bdt_v/F");

   Float_t new_m;
   TBranch *newBranch = resTree->Branch("new_m", &new_m, "new_m/F");


   // Book output histograms
   UInt_t nbin = 100;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );



   TH1F *histMass = new TH1F("mass", "mass", 44, 4.85, 5.95);
   //TH1F *histMass = new TH1F("mass", "mass", 40, 4.9, 5.9);
   TH1F *histPt = new TH1F("pt", "pt", 50, 0., 50.);
   TH1F *histEta = new TH1F("eta", "eta", 24, -2.4, 2.4);
   TH1F *histFls3d = new TH1F("fls3d", "fls3d", 100, 0., 100.);
   TH1F *histAlpha = new TH1F("alpha", "alpha", 30, 0., 0.3);
   TH1F *histMaxdoca = new TH1F("maxdoca", "maxdoca", 50, 0., 0.05);
   TH1F *histPvip = new TH1F("pvip", "pvip", 50, 0., 0.05);
   TH1F *histPvips = new TH1F("pvips", "pvips", 50, 0., 5.);
   TH1F *histIso = new TH1F("iso", "iso", 100, 0., 1.);
   TH1F *histDocatrk = new TH1F("docatrk", "docatrk", 20, 0., 0.2);
   TH1F *histClosetrk = new TH1F("closetrk", "closetrk", 21, 0, 21);
   TH1F *histChi2dof = new TH1F("chi2dof", "chi2dof", 100, 0., 10.);

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   // TFile *input(0);
   // TString fname = "";

   // TString fnameS = "BsMC12_barrel_preselection.root";
   // TString fnameS = "Endcaps_preselection_unblinded.root";
   //TString fnameS = "Barrel_preselection_unblinded.root";

   // if (!gSystem->AccessPathName( fname ))
   //    input = TFile::Open( fname ); // check if file in local directory exists
   // else
   //    input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
   // 
   // if (!input) {
   //    std::cout << "ERROR: could not open data file" << std::endl;
   //    exit(1);
   // }
   // std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // TFile *inputS = TFile::Open( fnameS );
   // TFile *inputB = TFile::Open( fnameB );


  
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   // std::cout << "--- Select signal sample" << std::endl;
   // TTree* theTree = (TTree*)input->Get("TreeS");
   // Float_t userVar1, userVar2;
   // theTree->SetBranchAddress( "var1", &userVar1 );
   // theTree->SetBranchAddress( "var2", &userVar2 );
   // theTree->SetBranchAddress( "var3", &var3 );
   // theTree->SetBranchAddress( "var4", &var4 );

 
   // Efficiency calculator for cut method
   Int_t    nSelCutsSA = 0;
   Double_t effS       = cutValue;


   bool preSelPass = false;
   bool maskDuplicateEvents = false;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      //if(ievt>4000) continue;

      theTree->GetEntry(ievt);
      if (!main) {  
        evtnum = (Long64_t)evtnum_xcheck;
        closetrk = (Int_t)closetrk_xcheck;
      }
      /*
      new_m = m;  
      bdt_v = 10.2;
      printf("==: %f %f\n",m,bdt_v);

      histMass->Fill(m);
      bdtOutput->Fill(); //this does not work for some reason XXXXX
      newBranch->Fill();
      resTree->Fill();

      continue;
      */


      //TString massCut = "m > 4.9 && m < 5.9";
      //      massCut = massCut + "&& ( m < 5.2 || m > 5.45)";
      bool massCutPass = m > 4.9 && m < 5.9;

      preSelPass =  me < 0.2 && pt > 5. && pt < 9999. && m1pt > 4. && m2pt > 4. && m1eta<2.0 && m1eta>-2.0 && m2eta<2.0 && m2eta>-2.0 && fl3d < 2. && fls3d > 0. && fls3d < 200. && chi2dof < 20. && pvip < 0.1 && pvips < 5. && maxdoca < 0.1 && alpha < 1. && closetrk < 21 && docatrk < 2.5 && iso > 0. && pvlip < 1.0 && pvlip>-1.0 && pvlips < 5.0 && pvlips>-5.0 && pvw8>0.7;  // && (lxysig > 3.) 

      //maskDuplicateEvents = evtnum == 659534454 || evtnum == 1127720136 || evtnum == 262785584 || evtnum == 249102695 || (evtnum == 491167455 && 18.5-fls3d<0.1) || (evtnum==997510994 && 52.5-fls3d<0.1);

   preSelPass = preSelPass && !maskDuplicateEvents && massCutPass;

      //preSelPass = false;
      //preSelPass =  m > 4.9 && m < 5.9 && me < 0.2 && m1eta<2.0 && m1eta>-2.0 && m2eta<2.0 && m2eta>-2.0;// &&  m1pt>4.0 && m2pt>4.0;
      //cout << "pass:" << preSelPass;
      //preSelPass = preSelPass && pvw8>0.7;  
      //cout << " " << preSelPass;
      //preSelPass = preSelPass && fl3d < 2. && fls3d > 0. && fls3d < 200. && chi2dof < 20. && pvip < 0.1 && pvips < 5. && maxdoca < 0.1 && alpha < 1. && closetrk < 21 && docatrk < 2.5 && iso > 0. && pvlip < 1.0 && pvlip>-1.0 && pvlips < 5.0 && pvlips>-5.0;
      //cout << " " << preSelPass << endl;
   // && (lxysig > 3.) 
   
      if(!preSelPass) {
	printf("%f(4.9-5.9) %f(<0.2) %f(<2) %f(<2) %f(>4) %f(>4) %f(>0.7)\t",m,me,m1eta,m2eta,m1pt,m2pt,pvw8);
	printf("%f(<1) %f(<5) %f(>0) %d(<21) %f(<0.1) %f(<2.5) %f(<1)\n",pvlip,pvlips,iso,closetrk,maxdoca,docatrk,alpha);
	cout << "failed pre selection!" << endl;
	continue;
      }
   /*
   TString preselection("(m > 4.9 && m < 5.9 && me < 0.2 && pt > 5. && pt < 9999. && m1pt > 4. && m1pt < 999. && m2pt > 4. && m2pt < 999. && fl3d < 2. && fls3d > 0. && fls3d < 200. && chi2dof < 20. && pvip < 0.1 && pvips < 5. && maxdoca < 0.1 && alpha < 1. && closetrk < 21 && docatrk < 2.5 && iso > 0. && (mu1_charge*mu2_charge == -1) && (lxysig > 3.) && (abs(pvlip) < 1.0) && (abs(pvlips) < 5.0) && pvw8>0.7");
   TString barrelCuts("(m1eta < 1.4 && m1eta > -1.4 && m2eta < 1.4 && m2eta > -1.4)");
   if( endcaps ) 
   TString endcapsCuts("!"+barrelCuts+" && m1eta < 2.0 && m2eta < 2.0 && m1eta > -2.0 && m2eta > -2.0");
   */




      //assign variables for mva reader (need to eb floats)

      m_       =(Float_t)m       ; 
      fls3d_   =(Float_t)fls3d 	 ;
      alpha_   =(Float_t)alpha 	 ;
      pvips_   =(Float_t)pvips 	 ;
      iso_     =(Float_t)iso 	 ;
      m1iso_   =(Float_t)m1iso 	 ;
      m2iso_   =(Float_t)m2iso 	 ;
      chi2dof_ =(Float_t)chi2dof ;
      eta_     =(Float_t)eta 	 ;
      pt_      =(Float_t)pt 	 ;
      maxdoca_ =(Float_t)maxdoca ;
      docatrk_ =(Float_t)docatrk ;
      pvip_    =(Float_t)pvip    ;
      

      //cout << " xxx m:" << m << " m_:" << m_ << endl;

      // --- Return the MVA outputs and fill into histograms

      if (Use["CutsSA"]) {
	assert(1); //should not pass here
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsSA method", effS );
         if (passed) { 
//	   nSelCutsSA++;
//	   histMass->Fill(m);
//	   histPt->Fill(pt);
//	   histEta->Fill(eta);
//	   histFls3d->Fill(fls3d);
//	   histAlpha->Fill(alpha);
//	   histMaxdoca->Fill(maxdoca);
//	   histPvip->Fill(pvip);
//	   histPvips->Fill(pvips);
//	   histIso->Fill(iso);
//	   histDocatrk->Fill(docatrk);
//	   //	   histClosetrk->Fill(ntrk);
//	   histChi2dof->Fill(chi2dof);
	 }
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   {
	assert(1); //should not pass here!!
	histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
	if (reader->EvaluateMVA( "MLP method"           ) > cutValue) {
//	  histMass->Fill(m);
//	  histPt->Fill(pt);
//	  histEta->Fill(eta);
//	  histFls3d->Fill(fls3d);
//	  histAlpha->Fill(alpha);
//	  histMaxdoca->Fill(maxdoca);
//	  histPvip->Fill(pvip);
//	  histPvips->Fill(pvips);
//	  histIso->Fill(iso);
//	  histDocatrk->Fill(docatrk);
//	  //	  histClosetrk->Fill(ntrk);
//	  histChi2dof->Fill(chi2dof);
	}

      }
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["BDT"          ]) {
	//should always pass here!
	//	histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
	if (reader->EvaluateMVA( "BDT method"           ) > cutValue) {


	  new_m = m;  
	  bdt_v = (reader->EvaluateMVA( "BDT method"));
	  printf("==: %f %f\n",m,bdt_v);
	  
	  //	  histMass->Fill(m);
	  histBdt ->Fill(bdt_v); //this works

	  //bdtOutput->Fill(); //this does not work for some reason XXXXX
	  //newBranch->Fill();
	  resTree->Fill();
	  
	  //bdt_v = (mytype)(reader->EvaluateMVA( "BDT method"));
	  //resTree->Fill(); 
	  //bdtOutput->Fill(); //this does not work for some reason XXXXX

	  cout << "BDT value:" << reader->EvaluateMVA( "BDT method") << endl;
	  printf("%s %f %f\n",region.Data(),m,bdt_v);

	  //note: when running either barrel or endcap, pt or eta will have no sense as they are only read in turn for one of the two cases
	  printf("-- m:%3.2f pt:%3.2f eta:%3.2f fls3d:%3.2f alpha:%3.2f maxdoca:%3.2f pvip:%3.2f pvips:%3.2f iso:%3.2f docatrk:%3.2f chi2dof:%3.2f\n",m, pt, eta, fls3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof);

	  //	  continue;


	  histMass->Fill(m);
	  histPt->Fill(pt);
	  histEta->Fill(eta);
	  histFls3d->Fill(fls3d);
	  histAlpha->Fill(alpha);
	  histMaxdoca->Fill(maxdoca);
	  histPvip->Fill(pvip);
	  histPvips->Fill(pvips);
	  histIso->Fill(iso);
	  histDocatrk->Fill(docatrk);
	  //	  histClosetrk->Fill(ntrk);
	  histChi2dof->Fill(chi2dof);
	}
	// if (reader->EvaluateMVA( "BDT method"           ) > 0.1361) histMass->Fill(mass);
	// if (reader->EvaluateMVA( "BDT method"           ) > 0.2163) histMass->Fill(mass);
	// if (reader->EvaluateMVA( "BDT method"           ) > 0.1547) histMass->Fill(mass);
      }
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsSA"]) std::cout << "--- Efficiency for CutsSA method: " << double(nSelCutsSA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsSA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsSA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( effS, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of " << effS << " from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }


   // --- Write filtered tree   

   resTree->Write();  
   resTree->Write("", TObject::kOverwrite);

   // --- Write histograms

   //TFile *target  = new TFile( outputFileName,"RECREATE" );
   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   histMass->Write();
   histPt->Write();
   histEta->Write();
   histFls3d->Write();
   histAlpha->Write();
   histMaxdoca->Write();
   histPvip->Write();
   histPvips->Write();
   histIso->Write();
   histDocatrk->Write();
   histClosetrk->Write();
   histChi2dof->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }


   std::cout << "--- Created root file: \"" << outputFileName << "\" containing the MVA output histograms" << std::endl;
  
   delete reader;

   // histMass->Draw();
   // histPt->Draw();
   // histEta->Draw();
   // histFls3d->Draw();
   // histAlpha->Draw();
   // histMaxdoca->Draw();
   // histPvip->Draw();
   // histPvips->Draw();
   // histIso->Draw();
   // histDocatrk->Draw();
   // histClosetrk->Draw();
   // histChi2dof->Draw();


   resTreeFile->Close();
   //target->Close();
   inputS->Close();


  
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
