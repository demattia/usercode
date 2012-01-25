#ifndef CutFlow_h
#define CutFlow_h 1

#include <string>
#include <map>

class CutFlow {
  
 public:
 
  CutFlow();
  ~CutFlow() {};

  // register a specific cut result.
  void applyCut(std::string cutName, double value, double minValue, double maxValue);
  void applyCut(std::string cutName, bool passThisCut);

  // cuts registered by the following function are treated separately
  void applyLifetimeCut(std::string cutName, double value, double minValue, double maxValue);
  void applyLifetimeCut(std::string cutName, bool passThisCut);

  // check whether all cuts were passed
  bool passAll();

  // check whether a particular cut was passed
  bool passCut(std::string cutName);

  // return value
  double getValue(std::string cutName);

  // check whether all cuts except a particular one were passed.
  // two options: either require requested cut to be failed, or ignore result of this cut.
  bool passAllFailOne(std::string cutName);
  bool passAllIgnoreOne(std::string cutName);
  // third option: ignore result of this cut and of all lifetime cuts
  bool passAllIgnoreOneIgnoreLifetime(std::string cutName);
  // or how about this one:
  bool passAllIgnoreLifetime();

 private:

  friend class CutSummary;

  std::map<std::string,double> values_;
  std::map<std::string,bool> cutResults_;
  std::map<std::string,bool> lifetimeCutResults_;
  bool passAll_;
  bool passAllIgnoreLifetime_;

} ;

#endif
