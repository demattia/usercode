#include "CondCore/PluginSystem/interface/registration_macros.h"
#include "CondFormats/MomentumScaleCalibrationObjects/interface/MuScleFitLikelihoodPdf.h"
#include "CondFormats/DataRecord/interface/MuScleFitLikelihoodPdfRcd.h"

DEFINE_SEAL_MODULE();
REGISTER_PLUGIN(MuScleFitLikelihoodPdfRcd,MuScleFitLikelihoodPdf);

