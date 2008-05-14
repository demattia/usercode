#ifndef XSECNLO_H
#define XSECNLO_H

namespace anaobj {

  // Higgs production cross-Section @ NLO
  // as function as its mass value
  // xSec errors are due to the scale variation
  //  (in per cent and different above and below)
  // - the errors for the VBF are unprecise,
  //  since they are of the order of the integration errors in many cases
  // - the scale dependence of VBF is of the order of +-2%
  // https://twiki.cern.ch/twiki/bin/view/CMS/GeneratorProduction2007CSA07Signal

  // m_H = 160 GeV
  // -------------
  const double xSec_ggH160_  = 21.0; // in pb	  
  const double xSecPosUncert_ggH160_ =  18.; // %   
  const double xSecNegUncert_ggH160_ = -14.; // %	  
  const double xSec_VBFH160_ = 3.32; // in pb	  
  const double xSecPosUncert_VBFH160_ =  0.5; // % 
  const double xSecNegUncert_VBFH160_ = -1.0; // % 

  // m_H = 200 GeV
  // -------------
  const double xSec_ggH200_  = 14.8; // in pb	  
  const double xSecPosUncert_ggH200_ =  17.; // %   
  const double xSecNegUncert_ggH200_ = -14.; // %	  
  const double xSec_VBFH200_ = 2.53; // in pb	  
  const double xSecPosUncert_VBFH200_ =  0.4; // % 
  const double xSecNegUncert_VBFH200_ = -1.1; // % 

  // m_H = 400 GeV
  // -------------
  const double xSec_ggH400_  = 7.88; // in pb	  
  const double xSecPosUncert_ggH400_ =  17.; // %   
  const double xSecNegUncert_ggH400_ = -14.; // %	  
  const double xSec_VBFH400_ = 0.869; // in pb	  
  const double xSecPosUncert_VBFH400_ =  0.2; // % 
  const double xSecNegUncert_VBFH400_ = -1.1; // % 

  // m_H = 800 GeV
  // -------------
  const double xSec_ggH800_  = 0.397; // in pb	  
  const double xSecPosUncert_ggH800_ =  16.; // %   
  const double xSecNegUncert_ggH800_ = -14.; // %	  
  const double xSec_VBFH800_ = 0.196; // in pb	  
  const double xSecPosUncert_VBFH800_ =  0.7; // % 
  const double xSecNegUncert_VBFH800_ = -1.5; // % 
}

#endif //XSECNLO_H
