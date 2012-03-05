#ifndef ROOTTREEHANDLER_H
#define ROOTTREEHANDLER_H

#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include <Analysis/RootTreeProducers/interface/GenParticle.h>
#include <Analysis/RootTreeProducers/interface/Track.h>
#include <TH1F.h>
#include <stdlib.h>
#include <map>

/**
 * This class can be used to save an arbitrary number of vectors<type> in a root tree. <br>
 * The constructor opens the file and creates the tree <br>
 * The method addBranch(name, type) (where type must inherit from TObject) adds a branch of
 * given name with a vector<type>. A map stores the pointers to the objects and the keys are the names.
 */

class RootTreeHandler
{
public:
  RootTreeHandler(const TString & fileName, const TString & treeName);
  // ~RootTreeHandler();

  /** The method must be called specifying the type. The type will be used to create a pointer in the
  * map and to add a branch to the tree.
  */
  template <class T>
  void addBranch(const TString & name, const TString & collectionType, T* pointer)
  {
//    tree_->Branch("tracks", vectorName+object.GetName()+">", &objects_[name]);
//    T object;
//    if( !object.InheritsFrom("TObject") ) {
//      std::cout << "Error: object type does not inherit from TObject" << std::endl;
//      throw;
//    }
//    TString vectorName("std::vector<");
//    objects_.insert(std::pair<TString, std::vector<T>*>(name, 0));
//    std::cout << "collection name = " << vectorName+object.GetName()+">" << std::endl;
//    tree_->Branch("tracks", vectorName+object.GetName()+">", &objects_[name]);
    tree_->Branch(name, collectionType, pointer);
  }

//  template <class T>
//  void fill( const TString & name, const std::vector<T> & collection )
//  {
//    objects_[name] = &collection;
//  }


  void saveToTree( const unsigned int eventNumber, const unsigned int runNumber );
  void writeTree();

protected:
  TFile * file_;
  TTree * tree_;
  unsigned int eventNumber_;
  unsigned int runNumber_;
  // std::map<TString, std::vector<TObject>* > objects_;
};

#endif // ROOTTREEHANDLER_H
