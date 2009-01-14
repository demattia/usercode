#ifndef MuScleFitDBobject_h
#define MuScleFitDBobject_h

#include <vector>

struct MuScleFitDBobject
{
  std::vector<int> identifiers;
  std::vector<std::vector<double> > parameters;
};

#endif // MuScleFitDBobject
