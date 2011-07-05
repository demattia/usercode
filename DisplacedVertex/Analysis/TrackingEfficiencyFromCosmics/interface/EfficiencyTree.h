#ifndef EFFICIENCYTREE_H
#define EFFICIENCYTREE_H

#include <iostream>

#include <TFile.h>

#include <Analysis/TrackingEfficiencyFromCosmics/interface/SerializableEfficiency.h>
#include <Analysis/TrackingEfficiencyFromCosmics/interface/Efficiency.h>
#include <stdlib.h>

/**
 * This class allows to write the efficiency objects to file and read them back.
 */

class EfficiencyTree
{
public:
  void writeTree( const TString & fileName, const Efficiency * eff)
  {
    TFile * f1 = new TFile(fileName, "RECREATE");
    SerializableEfficiency sEff;

    // Set all the values

//    sEff.S = eff->linearSize();
//    sEff.N = eff->getN();
    unsigned int N = eff->getN();
    for( unsigned int i=0; i<N; ++i ) {
      sEff.vIndexes.push_back(eff->getVIndexes(i));
      sEff.vSizes.push_back(eff->bins(i));
      sEff.vBinSizes.push_back(eff->binsSize(i));
      sEff.vMin.push_back(eff->min(i));
      sEff.vMax.push_back(eff->max(i));
      sEff.names.push_back(eff->getName(i));
    }
    unsigned int S = eff->getLinearSize();
    for( unsigned int i=0; i<S; ++i ) {
      sEff.values.push_back(eff->getValues(i));
    }

    sEff.Write("Efficiency");
    f1->Write();
    f1->Close();
  }

  void readTree( const TString & fileName, Efficiency * eff )
  {
    TFile * file = new TFile(fileName, "READ");
    if( file->IsOpen() ) {
      SerializableEfficiency * sEff = (SerializableEfficiency*)file->Get("Efficiency");
      unsigned int N = sEff->vIndexes.size();
      unsigned int S = sEff->values.size();
      boost::shared_array<unsigned int> vIndexes(new unsigned int[N]);
      boost::shared_array<unsigned int> vSizes(new unsigned int[N]);
      boost::shared_array<double> vBinSizes(new double[N]);
      boost::shared_array<double> vMin(new double[N]);
      boost::shared_array<double> vMax(new double[N]);
      for( unsigned int i=0; i<N; ++i ) {
        vIndexes[i] = sEff->vIndexes[i];
        vSizes[i] = sEff->vSizes[i];
        vBinSizes[i] = sEff->vBinSizes[i];
        vMin[i] = sEff->vBinSizes[i];
        vMax[i] = sEff->vBinSizes[i];
      }
      boost::shared_array<std::pair<unsigned int, unsigned int> > values(new std::pair<unsigned int, unsigned int>[S]);
      for( unsigned int i=0; i<S; ++i ) {
        values[i] = sEff->values[i];
      }
      // Set the values in the efficiency object
      eff->setN(N);
      eff->setLinearSize(S);
      eff->setVIndexes(vIndexes);
      eff->setVSizes(vSizes);
      eff->setVBinSizes(vBinSizes);
      eff->setVMin(vMin);
      eff->setVMax(vMax);
      eff->setValues(values);
      eff->setNames(sEff->names);
    }
    file->Close();
  }
};

#endif
