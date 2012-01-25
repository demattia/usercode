#ifndef HistMap_h
#define HistMap_h 1

#include <string>
#include <map>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class HistMap {
  
 public:
 
  HistMap(std::string folderName="");
  ~HistMap();
    
  void fill(std::string histname,double val, double weight=1,
	    std::string title="", int nbins=0, double min=0, double max=0,
	    std::string axistitle="");
  void fill(std::string histname,double xval, double yval, double weight=1,
	    std::string title="", int nbinsx=0, double xmin=0, double xmax=0,
	    int nbinsy=0, double ymin=0, double ymax=0,
	    std::string xaxistitle="", std::string yaxistitle=""); 

  TH1* getHist1D(const std::string& histname);  ///< Returns a null pointer if histogram doesn't exist.
  TH2* getHist2D(const std::string& histname);  ///< Returns a null pointer if histogram doesn't exist.
  bool makeEfficiencyHistogram(const std::string& selected,
			       const std::string& all,
			       const std::string& effi);
  bool makeEfficiencyHistogram(const std::string& selected,
			       const TH1* all,
			       const std::string& effi);

  /// Multiplies the histogram by the given scale factor.
  /*! False is returned on error. */
  bool scaleHistogram(const std::string& histName, double scale);

  TTree* makeTree(const std::string& title);

 private:

  typedef std::map<std::string, TH1*> HistMap1D;
  typedef std::map<std::string, TH2*> HistMap2D;

  TFileDirectory* _histDir;

  HistMap1D _histmap1D;
  HistMap2D _histmap2D;

  std::map<std::string, double> _histmap1D_xmin;
  std::map<std::string, double> _histmap2D_xmin;
  std::map<std::string, double> _histmap2D_ymin;
  std::map<std::string, double> _histmap1D_xmax;
  std::map<std::string, double> _histmap2D_xmax;
  std::map<std::string, double> _histmap2D_ymax;
  std::map<std::string, bool> _histmap1D_xnan;
  std::map<std::string, bool> _histmap2D_xnan;
  std::map<std::string, bool> _histmap2D_ynan;

} ;

#endif
