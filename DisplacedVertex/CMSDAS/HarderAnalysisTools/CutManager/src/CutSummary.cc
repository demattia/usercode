#include <iostream>
#include "HarderAnalysisTools/CutManager/interface/CutSummary.h"


CutSummary::CutSummary(std::string title) : numCutNames_(0), bigTree_(NULL) {

  cutNameMap_.clear();

  histos_ = new HistMap(title);
  trees_.clear();
  treeVariables_.clear();
}


void CutSummary::addEntry(CutFlow& cuts,
			  std::map<std::string,double>& weightMap,
			  const int numDecays,
			  const double ctau1,
			  const double ctau2,
			  const double leptonD01,
			  const double leptonD02,
			  const double mass) {

  // temporary: only dump pile-up weights
  Float_t weight=1;
  Float_t weight_up=1;
  Float_t weight_down=1;
  if (weightMap["lumiWeight_central"]) {
    weight = weightMap["lumiWeight_central"];
  }
  if (weightMap["lumiWeight_high"]) {
    weight_up = weightMap["lumiWeight_high"];
  }
  if (weightMap["lumiWeight_low"]) {
    weight_down = weightMap["lumiWeight_low"];
  }

  // count how many candidates we looked at and how many pass all cuts
  histos_->fill("passAllCuts",(cuts.passAll() ? 1 : 0),1,
		      "does candidate pass all cuts",2,0,2);

  // the big tree for all variables: book it if not yet done
  if (!bigTree_) {
    bigTree_=histos_->makeTree("bigTree");
    bigTreeValues_.clear();
    bigTreeBools_.clear();
    for (unsigned i=0; i<cuts.values_.size(); i++) {
      bigTreeValues_.push_back(0);
      bigTreeBools_.push_back(false);
    }
    // tree entries for all variables
    unsigned idx=0;
    for (std::map<std::string,double>::iterator it=cuts.values_.begin();
	 it!=cuts.values_.end(); it++) {
      const std::string& cutName = it->first;
      bigTreeNames_.push_back(cutName);
      bigTree_->Branch(cutName.c_str(),&(bigTreeValues_[idx]),
		       (cutName+"/F").c_str());
      bigTree_->Branch((cutName+"_pass").c_str(),
		       &(bigTreeBools_[idx]),
		       (cutName+"_pass/F").c_str());
      idx++;
    }
    // additional entries for reweighting and diagnostics
    bigTree_->Branch("_mass",&bigTreeVars_.mass,"_mass/F");
    bigTree_->Branch("_weight",&bigTreeVars_.weight,"_weight/F");
    bigTree_->Branch("_weight_up",&bigTreeVars_.weight_up,"_weight_up/F");
    bigTree_->Branch("_weight_down",&bigTreeVars_.weight_down,"_weight_down/F");
    bigTree_->Branch("_numDecays",&bigTreeVars_.numDecays,"_numDecays/F");
    bigTree_->Branch("_ctau1",&bigTreeVars_.ctau1,"_ctau1/F");
    bigTree_->Branch("_ctau2",&bigTreeVars_.ctau2,"_ctau2/F");
    bigTree_->Branch("_leptonD01",&bigTreeVars_.leptonD01,"_leptonD01/F");
    bigTree_->Branch("_leptonD02",&bigTreeVars_.leptonD02,"_leptonD02/F");
    bigTree_->Branch("passesAllCuts",&bigTreeVars_.passesAllCuts,
		     "passesAllCuts/O");
    bigTree_->Branch("passesAllCutsIgnoreLifetime",
		     &bigTreeVars_.passesAllCutsIgnoreLifetime,
		     "passesAllCutsIgnoreLifetime/O");

  } else if (cuts.values_.size()!=bigTreeValues_.size()) {
    std::cout << "ERROR: CutSummary::AddEntry different tree size - was"
	      << bigTreeValues_.size() << ", is " << cuts.values_.size()
	      << std::endl;
    exit(999);
  }

  unsigned idx=0;
  for (std::map<std::string,double>::iterator it=cuts.values_.begin();
       it!=cuts.values_.end(); it++) {
    if (bigTreeNames_[idx]!=it->first) {
      std::cout << "ERROR: bigTree variables in wrong order: " << it->first
		<< " vs " << bigTreeNames_[idx] << " at index " << idx
		<< std::endl;
      exit(999);
    }
    bigTreeValues_[idx]=it->second;
    bigTreeBools_[idx]=cuts.passCut(it->first);
    std::cout << "KHDEBUG: " << it->first << "=" << bigTreeValues_[idx]
	      << ", pass=" << bigTreeBools_[idx] << std::endl;
    idx++;
  }
  bigTreeVars_.mass=mass;
  bigTreeVars_.weight=weight;
  bigTreeVars_.weight_up=weight_up;
  bigTreeVars_.weight_down=weight_down;
  bigTreeVars_.numDecays=numDecays;
  bigTreeVars_.ctau1=ctau1;
  bigTreeVars_.ctau2=ctau2;
  bigTreeVars_.leptonD01=leptonD01;
  bigTreeVars_.leptonD02=leptonD02;
  bigTreeVars_.passesAllCuts=(cuts.passAll());
  bigTreeVars_.passesAllCutsIgnoreLifetime=(cuts.passAllIgnoreLifetime());
  bigTree_->Fill();

  // look at all cut results for this candidate
  for (std::map<std::string,double>::const_iterator it=cuts.values_.begin();
       it!=cuts.values_.end(); it++) {
    const std::string& cutName = it->first;

    // make sure we take this kind of cut into account
    // (not all candidates are subjected to all cuts,
    // so there might be later additions)
    if (cutNameMap_.find(cutName)==cutNameMap_.end()) {
      cutNameMap_[cutName]=numCutNames_;
      ++numCutNames_;
      TTree* newTree = histos_->makeTree(cutName);
      tree_type* newTreeVars = new tree_type();
      newTree->Branch("value",&newTreeVars->value,"value/F");
      newTree->Branch("mass",&newTreeVars->mass,"mass/F");
      newTree->Branch("weight",&newTreeVars->weight,"weight/F");
      newTree->Branch("weight_up",&newTreeVars->weight_up,"weight_up/F");
      newTree->Branch("weight_down",&newTreeVars->weight_down,"weight_down/F");
      newTree->Branch("numDecays",&newTreeVars->numDecays,"numDecays/F");
      newTree->Branch("ctau1",&newTreeVars->ctau1,"ctau1/F");
      newTree->Branch("ctau2",&newTreeVars->ctau2,"ctau2/F");
      newTree->Branch("leptonD01",&newTreeVars->leptonD01,"leptonD01/F");
      newTree->Branch("leptonD02",&newTreeVars->leptonD02,"leptonD02/F");
      newTree->Branch("passesThisCut",&newTreeVars->passesThisCut,
		      "passesThisCut/B");
      newTree->Branch("passesAllOtherCuts",&newTreeVars->passesAllOtherCuts,
		      "passesAllOtherCuts/B");
      newTree->Branch("passesAllOtherCutsIgnoreLifetime",
		      &newTreeVars->passesAllOtherCutsIgnoreLifetime,
		      "passesAllOtherCutsIgnoreLifetime/B");
      trees_.push_back(newTree);
      treeVariables_.push_back(newTreeVars);
    }
    unsigned cutNameIndex=cutNameMap_[cutName];

    // fill ROOT tree
    treeVariables_[cutNameIndex]->value=cuts.values_[cutName];
    treeVariables_[cutNameIndex]->mass=mass;
    treeVariables_[cutNameIndex]->weight=weight;
    treeVariables_[cutNameIndex]->weight_up=weight_up;
    treeVariables_[cutNameIndex]->weight_down=weight_down;
    treeVariables_[cutNameIndex]->numDecays=numDecays;
    treeVariables_[cutNameIndex]->ctau1=ctau1;
    treeVariables_[cutNameIndex]->ctau2=ctau2;
    treeVariables_[cutNameIndex]->leptonD01=leptonD01;
    treeVariables_[cutNameIndex]->leptonD02=leptonD02;
    treeVariables_[cutNameIndex]->passesThisCut=(cuts.passCut(cutName));
    treeVariables_[cutNameIndex]->passesAllOtherCuts=(cuts.passAllIgnoreOne(cutName));
    treeVariables_[cutNameIndex]->passesAllOtherCutsIgnoreLifetime=(cuts.passAllIgnoreOneIgnoreLifetime(cutName));
    trees_[cutNameIndex]->Fill();
  }
}
