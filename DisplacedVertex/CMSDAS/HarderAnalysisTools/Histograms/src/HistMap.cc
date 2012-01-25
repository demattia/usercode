#include "HarderAnalysisTools/Histograms/interface/HistMap.h"
#include "FWCore/ServiceRegistry/interface/ServiceRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

HistMap::HistMap(std::string folderName) {

  edm::Service<TFileService> fs;
  if (folderName!="") {
    _histDir = new TFileDirectory(fs->mkdir(folderName));
  } else {
    _histDir = &edm::ServiceRegistry::instance().get<TFileService>();
  }

  _histmap1D.clear();
  _histmap2D.clear();
}

HistMap::~HistMap() {
}

void HistMap::fill(std::string histname,double val, double weight,
                   std::string title, int nbins, double min, double max,
                   std::string axistitle) {

  if (!_histmap1D[histname]) {
    if (nbins==0) {
      //streamlog_out(ERROR) << "must define histogram "+histname+
      //" properties on first call of fillhisto(...)";
      exit(1);
    }
    _histmap1D[histname] = _histDir->make<TH1F>(histname.c_str(),title.c_str(),
                                                nbins,min,max);
    if (!_histmap1D[histname]) {
      //streamlog_out(ERROR) << "failed to create histogram "
      //+_dirName+"/"+histname+" (out of memory?)";
      exit(1);
    }
    _histmap1D[histname]->Sumw2();
    _histmap1D_xnan[histname]=0;
  }

  // fill histogram and record any unusual entries (under-/overflow, NaN)
  if ((val>=0) || (val<0)) {
    if (_histmap1D[histname]->GetEntries()==0) {
      _histmap1D_xmin[histname] = val;
      _histmap1D_xmax[histname] = val;
    } else {
      _histmap1D_xmin[histname] = std::min(val,_histmap1D_xmin[histname]);
      _histmap1D_xmax[histname] = std::max(val,_histmap1D_xmax[histname]);
    }
    _histmap1D[histname]->Fill(val,weight);
  } else {
    ++_histmap1D_xnan[histname];
  }
}


void HistMap::fill(std::string histname,double xval, double yval,
                   double weight,
                   std::string title, int nbinsx, double xmin, double xmax,
                   int nbinsy, double ymin, double ymax,
                   std::string xaxistitle, std::string yaxistitle) {

  
  if (!_histmap2D[histname]) {
    if (nbinsx==0) {
      //streamlog_out(ERROR) << "must define histogram "+histname
      //+" properties on first call of fillhisto(...)";
      exit(1);
    }

    _histmap2D[histname] = _histDir->make<TH2F>(histname.c_str(),title.c_str(),
                                                nbinsx, xmin, xmax,
                                                nbinsy, ymin, ymax ); 
    if (!_histmap2D[histname]) {
      //streamlog_out(ERROR) << "failed to create histogram "
      //+_dirName+"/"+histname+" (out of memory?)";
      exit(1);
    }
    _histmap2D_xnan[histname]=0;
  }
  _histmap2D[histname]->Fill(xval,yval,weight);
}

TH1* HistMap::getHist1D(const std::string& histname)
{
  HistMap1D::iterator result = _histmap1D.find(histname);
  if(result != _histmap1D.end()) { return result->second; }
  // throw cms::Exception("InvalidHistogram") << "No 1D histogram with name '" << histname << "' exists!" << std::endl;
  return 0;
}

TH2* HistMap::getHist2D(const std::string& histname)
{
  HistMap2D::iterator result = _histmap2D.find(histname);
  if(result != _histmap2D.end()) { return result->second; }
  //throw cms::Exception("InvalidHistogram") << "No 2D histogram with name '" << histname << "' exists!" << std::endl;
  return 0;
}


bool HistMap::makeEfficiencyHistogram(const std::string& selected,
				      const std::string& all,
				      const std::string& effi) {
  
  HistMap1D::iterator selhist = _histmap1D.find(selected);
  if (selhist == _histmap1D.end()) return false;
  HistMap1D::iterator refhist = _histmap1D.find(all);
  if (refhist == _histmap1D.end()) return false;

  _histmap1D[effi]= _histDir->make<TH1F>(effi.c_str(),selhist->second->GetTitle(),
					 selhist->second->GetNbinsX(),
					 selhist->second->GetXaxis()->GetXmin(),
					 selhist->second->GetXaxis()->GetXmax());
  _histmap1D[effi]->Divide(selhist->second,refhist->second,1,1,"B");
  _histmap1D_xnan[effi]=0;
  _histmap1D_xmin[effi] = _histmap1D_xmin[selected];
  _histmap1D_xmax[effi] = _histmap1D_xmax[selected];
  return true;
}


bool HistMap::makeEfficiencyHistogram(const std::string& selected,
				      const TH1* all,
				      const std::string& effi) {
  
  HistMap1D::iterator selhist = _histmap1D.find(selected);
  if (selhist == _histmap1D.end()) return false;
  if (!all) return false;

  _histmap1D[effi]= _histDir->make<TH1F>(effi.c_str(),selhist->second->GetTitle(),
					 selhist->second->GetNbinsX(),
					 selhist->second->GetXaxis()->GetXmin(),
					 selhist->second->GetXaxis()->GetXmax());
  _histmap1D[effi]->Divide(selhist->second,all,1,1,"B");
  _histmap1D_xnan[effi]=0;
  _histmap1D_xmin[effi] = _histmap1D_xmin[selected];
  _histmap1D_xmax[effi] = _histmap1D_xmax[selected];
  return true;
}


bool HistMap::scaleHistogram(const std::string& histName, double scale)
{
  TH1* plot = getHist1D(histName);
  if(!plot) { return false; }
  plot->Scale(scale);
  return true;
}


TTree* HistMap::makeTree(const std::string& title) {

  return _histDir->make<TTree>(title.c_str(),title.c_str());
}
