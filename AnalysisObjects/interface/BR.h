#ifndef BR_H
#define BR_H

namespace anaobj {

  // Higgs branching ratio for m_top = 175
  // as function as its mass value
  // -------------------------------------
  const double BR_H160ZZ_ = 0.4334;  // m_H = 160 GeV
  const double BR_H200ZZ_ = 0.2613;  // m_H = 200 GeV
  const double BR_H400ZZ_ = 0.2742;  // m_H = 400 GeV
  const double BR_H800ZZ_ = 0.2928;  // m_H = 800 GeV

  // Z branching ratio
  // -----------------
  const double BR_Zee_     = 0.03363;
  const double BR_Zmumu_   = 0.03366;
  const double BR_Ztautau_ = 0.03370;
  const double BR_Znunu_   = 0.2010;
  const double BR_Zuu_     = 0.1117;
  const double BR_Zdd_     = 0.1582;
  const double BR_Zcc_     = 0.1203;
  const double BR_Zss_     = 0.1582;
  const double BR_Ztt_     = 0.;
  const double BR_Zbb_     = 0.1512;

  const double BR_Zll_ = 3*0.033658;
  const double BR_Zqq_ = 0.698;
}

#endif //BR_H
