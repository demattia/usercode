#ifndef EVENTSREADER_H
#define EVENTSREADER_H

#include <vector>

using namespace std;

class EventsReader
{
public:
  EventsReader( const vector<TString> & eventVariablesNamesVector, const TString & fileName, const TString & treeName, const int variablesNumber )
  {
    //
    // create the Reader object
    //
    reader_ = new TMVA::Reader("!Color");    

    variables_ = new Float_t[variablesNumber];

    // create a set of variables and declare them to the reader
    // - the variable names must corresponds in name and type to 
    // those given in the weight file(s) that you use

    // All minus the last one which is the event number

    for( int iVar = 0; iVar < variablesNumber - 1; ++iVar ) {
      reader_->AddVariable( eventVariablesNamesVector[iVar], &(variables_[iVar]) );
    }

    //
    // book the MVA methods
    //
    string dir    = "weights/";
    string prefix = "TMVAnalysis";

    reader_->BookMVA( "BDT method", dir + prefix + "_BDT.weights.txt" );

    // book output histograms
    UInt_t nbin = 100;
    TString name("MVA_BDT");
    histBdt_ = new TH1F( name+treeName, name+treeName, nbin, -0.8, 0.8 );

    // book examsple histogram for probability (the other methods are done similarly)
    name = "PROBA_MVA_BDT";
    probHistBdt_   = new TH1F( name+treeName, name+treeName, nbin, 0, 1 );
    name = "RARITY_MVA_BDT";
    rarityHistBdt_ = new TH1F( name+treeName, name+treeName, nbin, 0, 1 );

    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //

    TFile *input = TFile::Open( fileName );

    if (!input) {
      cout << "ERROR: could not open data file: " << fileName << endl;
      exit(1);
    }

    //
    // prepare the tree
    // - here the variable names have to corresponds to your tree
    // - you can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //

    theTree_ = (TTree*)input->Get(treeName);
    cout << "--- Select signal sample" << endl;

    vector<TString>::const_iterator varName = eventVariablesNamesVector.begin();
    int iVar = 0;
    for( ; varName != eventVariablesNamesVector.end(); ++varName, ++iVar ) {
      theTree_->SetBranchAddress( *varName, &(variables_[iVar]) );
    }
  }

  /// return the discriminant value
  double discriminant()
  {
    // ATTENTION: if the same method is applyied two times, the second time
    // does not work.
    // Also the variables_ array after application is filled with random values.
    double discriminator = reader_->EvaluateMVA( "BDT method" );
    histBdt_->Fill( discriminator );
    return discriminator;
  }

  /// retrieve probability instead of MVA output
  void probability()
  {
    probHistBdt_->Fill( reader_->GetProba ( "BDT method" ) );
  }
  /// retrieve rarity instead of MVA output
  void rarity()
  {
    rarityHistBdt_->Fill( reader_->GetRarity( "BDT method" ) );
  }
  /// Write histograms
  void write()
  {
    histBdt_->Write();
    //     probHistBdt_->Write();
    //     rarityHistBdt_->Write();
  }
  Float_t * variables()
  {
    return variables_;
  }
  /// Returns the corresponding variable as if accessing directly the variables[] array
  double operator[]( const int iVar ) const
  {
    return variables_[iVar];
  }
  /// Get the pointer to the TTree
  TTree * tree() const
  {
    return theTree_;
  }
protected:
  TMVA::Reader *reader_;
  TTree* theTree_;
  TH1F *histBdt_;
  TH1F *probHistBdt_, *rarityHistBdt_;
  Float_t * variables_;
};

#endif // EVENTSREADER_H
