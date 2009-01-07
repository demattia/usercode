#include "CondCore/PluginSystem/interface/registration_macros.h"
#include "MuonAnalysis/MomentumScaleCalibrationObjects/interface/MuScleFitLikelihoodPdf.h"
#include "MuonAnalysis/MomentumScaleCalibrationObjects/interface/MuScleFitLikelihoodPdfRcd.h"

DEFINE_SEAL_MODULE();
REGISTER_PLUGIN(MuScleFitLikelihoodPdfRcd,MuScleFitLikelihoodPdf);

