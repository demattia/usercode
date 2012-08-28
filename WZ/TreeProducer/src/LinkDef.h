#ifndef LINKDEF_H
#define LINKDEF_H

#include "TreeProducer/TreeProducer/interface/Candidates.h"
#include "TreeProducer/TreeProducer/interface/TreeCandidate.h"
#include "TreeProducer/TreeProducer/interface/TreeLepton.h"

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class Candidates+;
// #pragma link C++ class std::vector<std::string>+;
#pragma link C++ class TreeCandidate+;
#pragma link C++ class TreeLepton+;
#endif

#endif // LINKDEF_H
