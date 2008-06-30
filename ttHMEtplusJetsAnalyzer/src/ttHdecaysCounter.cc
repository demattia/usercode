#ifndef TTHDECAYSCOUNTER_CC
#define TTHDECAYSCOUNTER_CC

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHdecaysCounter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

ttHdecaysCounter::ttHdecaysCounter(const TString & OUTFILENAME) {
  higgsDecayNumbers_ = 5;
  ttDecayNumbers_ = 10;
  // Initializing the counters to 0
  for( int i=0; i<higgsDecayNumbers_; ++i ) {
    for( int j=0; j<ttDecayNumbers_; ++j ) {
        decayTypes_[i][j] = 0;
    }
  }
  eventCounter_ = 0;
  outFileName_ = OUTFILENAME;
}

void ttHdecaysCounter::countDecays(const MCParticleCollection & mcParticles) {

  ++eventCounter_;

  int foundw=0;
  int iparton=0;

  // Array with the pointers to the partons from the decay of the H and the ttbar for quicker access later.
  // At most 8, but not all are necessarily used.
  vector<const MCParticle *> decayPartons;
  vector<const MCParticle *> hadronicDecayPartons;

  // NB we know that the block is filled with top and antitop first, then Higgs
  // So we need not worry about messing up W's from top and from H ...
  // --------------------------------------------------------------------------
  int HdecayType = -1;
  int WplusDecayType = 0;
  int WminusDecayType = 0;

  for ( MCParticleCollection::const_iterator MCp = mcParticles.begin(); MCp != mcParticles.end(); ++MCp ) {

    int mPid = MCp->mPid();
    int pid = MCp->pid();
    int absMpid = abs(mPid);
    int absPid = abs(pid);

    // See if this is a daughter of the W
    if ( absMpid==24 && foundw<4 ) { // first two W's in list are always from t, tbar
      foundw++;
      // For hadronic W decay from top, store info on partons
      if ( iparton<8 ) {
        decayPartons.push_back(&(*MCp));
        // The decay types we distinguish for W+ and W- are 4: jets, enu, munu, taunu
        // For the W+
        if ( mPid == 24 ) {
          if ( absPid < 6 ) WplusDecayType = 1;      // (W+)->jetjet
          if ( absPid == 11 ) WplusDecayType = 10;   // (W+)->enu
          if ( absPid == 13 ) WplusDecayType = 100;  // (W+)->munu
          if ( absPid == 15 ) WplusDecayType = 1000; // (W+)->taunu
        }
        // For the W-
        if ( mPid == -24 ) {
          if ( absPid < 6 ) WminusDecayType = 1;      // (W-)->jetjet
          if ( absPid == 11 ) WminusDecayType = 10;   // (W-)->enu
          if ( absPid == 13 ) WminusDecayType = 100;  // (W-)->munu
          if ( absPid == 15 ) WminusDecayType = 1000; // (W-)->taunu
        }
        // Store the hadronic partons from W decays, if any
        if (absPid<6) {
          hadronicDecayPartons.push_back(&(*MCp));
        }
        iparton++;
      }
    }
    // Store information on b partons from t, tbar
    if ( absMpid==6 && absPid==5 && iparton<8 ) {
      decayPartons.push_back(&(*MCp));
      hadronicDecayPartons.push_back(&(*MCp));
      iparton++;
    }
    // Store the partons from the Higgs decay
    if ( mPid==25 && iparton<8 ) {
      decayPartons.push_back(&(*MCp));
      if (absPid<6) {
        hadronicDecayPartons.push_back(&(*MCp));
      }
      // Determine Higgs decay type
      if ( absPid == 4 ) HdecayType = 0;       // H->cc
      else if ( absPid == 5 ) HdecayType = 1;  // H->bb
      else if ( absPid == 15 ) HdecayType = 2; // H->tautau
      else if ( absPid == 24 ) HdecayType = 3; // H->WW
      else HdecayType = 4;                     // H->others

      iparton++;
    }
  }

  // If no Higgs was found we set HdecayType to 4 (H->others)
  // so that this will still work for ttbar events.
  // Since only events with two Ws from top quarks are considered,
  // this will not count anything for the backgrounds (QCD, W+jets, ...).
  // --------------------------------------------------------------------
  if ( HdecayType == -1 ) HdecayType = 4;

  // Parton vectors filled. Check the type of the event and increase the counter.
  if ( WplusDecayType + WminusDecayType == 2 )    decayTypes_[HdecayType][0]++; // tt->6jets
  if ( WplusDecayType + WminusDecayType == 11 )   decayTypes_[HdecayType][1]++; // tt->enu4jets
  if ( WplusDecayType + WminusDecayType == 101 )  decayTypes_[HdecayType][2]++; // tt->munu4jets
  if ( WplusDecayType + WminusDecayType == 1001 ) decayTypes_[HdecayType][3]++; // tt->taunu4jets
  if ( WplusDecayType + WminusDecayType == 20 )   decayTypes_[HdecayType][4]++; // tt->enuenu2jets
  if ( WplusDecayType + WminusDecayType == 110 )  decayTypes_[HdecayType][5]++; // tt->enumunu2jets
  if ( WplusDecayType + WminusDecayType == 1010 ) decayTypes_[HdecayType][6]++; // tt->enutaunu2jets
  if ( WplusDecayType + WminusDecayType == 200 )  decayTypes_[HdecayType][7]++; // tt->munumunu2jets
  if ( WplusDecayType + WminusDecayType == 1100 ) decayTypes_[HdecayType][8]++; // tt->munutaunu2jets
  if ( WplusDecayType + WminusDecayType == 2000 ) decayTypes_[HdecayType][9]++; // tt->taunutaunu2jets
}

void ttHdecaysCounter::writeDecays(bool append) {

  cout << "Total events            = " << eventCounter_ << endl;
  cout << "ttH->bb+6jets           = " << decayTypes_[1][0] << "; fraction = " << float(decayTypes_[1][0])/float(eventCounter_) << endl;
  cout << "ttH->bb+enu4jets        = " << decayTypes_[1][1] << "; fraction = " << float(decayTypes_[1][1])/float(eventCounter_) << endl;
  cout << "ttH->bb+munu4jets       = " << decayTypes_[1][2] << "; fraction = " << float(decayTypes_[1][2])/float(eventCounter_) << endl;
  cout << "ttH->bb+taunu4jets      = " << decayTypes_[1][3] << "; fraction = " << float(decayTypes_[1][3])/float(eventCounter_) << endl;
  cout << "ttH->bb+enuenu2jets     = " << decayTypes_[1][4] << "; fraction = " << float(decayTypes_[1][4])/float(eventCounter_) << endl;
  cout << "ttH->bb+enumunu2jets    = " << decayTypes_[1][5] << "; fraction = " << float(decayTypes_[1][5])/float(eventCounter_) << endl;
  cout << "ttH->bb+enutaunu2jets   = " << decayTypes_[1][6] << "; fraction = " << float(decayTypes_[1][6])/float(eventCounter_) << endl;
  cout << "ttH->bb+munumunu2jets   = " << decayTypes_[1][7] << "; fraction = " << float(decayTypes_[1][7])/float(eventCounter_) << endl;
  cout << "ttH->bb+munutaunu2jets  = " << decayTypes_[1][8] << "; fraction = " << float(decayTypes_[1][8])/float(eventCounter_) << endl;
  cout << "ttH->bb+taunutaunu2jets = " << decayTypes_[1][9] << "; fraction = " << float(decayTypes_[1][9])/float(eventCounter_) << endl;

  TString higgsDecaysNames[5] = { "H->cc", "H->bb", "H->tautau", "H->WW", "H->others" };
  for (int i=0; i<higgsDecayNumbers_; ++i) {
    int higgsDecays = 0;
    for (int j=0; j<ttDecayNumbers_; ++j) {
      higgsDecays += decayTypes_[i][j];
    }
    cout << higgsDecaysNames[i] + " number = " << higgsDecays << "; fraction = " << float(higgsDecays)/float(eventCounter_) << endl;
  }



  // Write the counts in the outputfile
  // ----------------------------------

  TString ttDecayNames[10] = { "tt->jjbjjb", "tt->evbjjb", "tt->mvbjjb", "tt->tvbjjb", 
                               "tt->evbenb", "tt->evbmvb", "tt->evbtvb", "tt->mnbmnb", 
                               "tt->mvbtvb", "tt->tnbtnb" };
  // File with HEP table
  // -------------------
  ofstream * decayFile = 0;
  // Create a new output file or append to an existing one
  if (append) decayFile = new ofstream(outFileName_, ios_base::app);
  else decayFile = new ofstream(outFileName_);

  for ( int h=0; h<5; h++ ) {
    *decayFile << "--------" << endl << higgsDecaysNames[h] << endl << "--------" << endl;
    for ( int tt=0; tt<10; tt++ ) {
      float eventsFrac = float(decayTypes_[h][tt])/eventCounter_;
      float eventsFracError = sqrt(eventsFrac*(1-eventsFrac)/eventCounter_);
      *decayFile << ttDecayNames[tt] << " = " << setprecision(8) << eventsFrac
                 << " +- "  << setprecision(8) << eventsFracError << endl;
    }
  }
  *decayFile << "----------------------------------------------------------------------" << endl;
  *decayFile << "Total events = " << eventCounter_ << endl;
  *decayFile << "----------------------------------------------------------------------" << endl;
  *decayFile << "----------------------------------------------------------------------" << endl;

  decayFile->close();
  delete decayFile;
}

#endif // TTHDECAYSCOUNTER_CC
