#ifndef SERIALIZABLEEFFICIENCY_H
#define SERIALIZABLEEFFICIENCY_H

#include <vector>
#include <utility>
#include <TObject.h>

struct SerializableEfficiency : public TObject
{
  SerializableEfficiency() {}

  // Data memebers

  std::vector<unsigned int> vSizes;
  std::vector<double> vBinSizes;
  std::vector<double> vMin;
  std::vector<double> vMax;
  std::vector<unsigned int> vIndexes;
  // Linear representation of the N-dimensional matrix
  std::vector<std::pair<unsigned int, unsigned int> > values;

  ClassDef(SerializableEfficiency, 1)
};
ClassImp(SerializableEfficiency)

#endif // SERIALIZABLEEFFICIENCY_H
