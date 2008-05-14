#ifndef XSECLO_H
#define XSECLO_H

namespace anaobj {

  // production cross-Section @ LO
  // https://twiki.cern.ch/twiki/bin/view/CMS/GeneratorProduction2007CSA07Signal
  // https://twiki.cern.ch/twiki/bin/view/Main/AlpgenSummer07

  // alpgen:
  // -------

  // ZZ+n jets
  const double xSec_ZZ0Jets_ = 1.369; // in pb
  const double xSec_ZZ1Jets_ = 0.637; // in pb
  const double xSec_ZZ2Jets_ = 0.247; // in pb

  // Zbb+n jets
  const double xSec_Zbb0Jets_ = 1.66; // in pb
  const double xSec_Zbb1Jets_ = 0.29; // in pb
  const double xSec_Zbb2Jets_ = 0.05; // in pb

}

#endif //XSECLO_H
