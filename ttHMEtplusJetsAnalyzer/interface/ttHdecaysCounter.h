#ifndef TTHDECAYSCOUNTER_H
#define TTHDECAYSCOUNTER_H

/**
 * Counts the different decay types in ttH events.
 * Receives a MCParticleCollection and fills two
 * vectors of MCParticles with all the hadronic partons
 * from ttH decay and all the partons from ttH decay.
 * It writes a text file with all the decays numbers and the full event count.
 * It also writes the decays for the ttbar pair, so it can be used for ttbar events.
 *
 * Author: M. De Mattia
 * Date: 29/6/2008
 */

#include <vector>
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "TString.h"

using namespace std;
using namespace anaobj;

class ttHdecaysCounter {
public:
  ttHdecaysCounter(const TString & outFileName = TString("ttHdecays.txt"));
  // ~ttHdecaysCounter();
  /// uses the MCParticleCollection to determine the type of event and adds the count in the matrix
  void countDecays(const MCParticleCollection & mcParticles);
  /** write the decay type counts in the text file. The bool it receives is used to decide if to append
   * to an existing file or to create a new one (deleting any previously existing file with the same name).
   * By default it creates a new file (append = false).
   */
  void writeDecays(bool append = false);

private:
  // Notation:
  // the first 5 are the considered decay modes of the Higgs: cc, bb, tautau, WW, others
  // W -> jetjet, enu, munu, taunu.
  // The second number is for the 10 possible decay modes of the tt pair (assuming t->bW).
  int decayTypes_[5][10];
  TString outFileName_;
  int eventCounter_;
  int higgsDecayNumbers_;
  int ttDecayNumbers_;
};

#endif // TTHDECAYSCOUNTER_H
