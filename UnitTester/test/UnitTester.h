#ifndef UnitTester_h
#define UnitTester_h
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "RecoLocalTracker/SiStripClusterizer/test/TestObjectWrappers.h"

#include <string>

using namespace std;

class UnitTester : public edm::EDAnalyzer {

  typedef edm::ParameterSet PSet;
  typedef std::vector<PSet> VPSet;
  typedef VPSet::const_iterator iter_t;
  
 public:
  
  UnitTester(const PSet& conf) : testGroups(conf.getParameter<VPSet>("TestGroups")) {}
  virtual ~UnitTester() {}

 protected:
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
  baseObjectWrapper * baseObjectFactory(const PSet & test, const edm::EventSetup & es ) {
    try {
      return objectFactory(test, es);
    }
    catch( const string & a ) {
      if( a == "unidentifiedType" ) {
        string type = test.getParameter<string>("Type");
//         if( type == "int" ) return new objectWrapper<int>( test.getParameter<uint32_t>("Value") );
//         else if( type == "double" ) return new objectWrapper<double>( test.getParameter<double>("Value") );
//         else if( type == "string" ) return new objectWrapper<string>( test.getParameter<string>("Value") );
//         else {
          cout << "Error: unrecognized type = " << type << endl;
          abort();
//         }
      }
      else throw( a );
    }
    return 0;
  }

  /**
   * This method reads the parameters and creates the objectWrappers. </br>
   * It needs to be implemented in the derived class.
   */
  virtual baseObjectWrapper * objectFactory(const PSet& test, const edm::EventSetup& es) = 0;

  void initializeTheGroup(const PSet &, const edm::EventSetup &);
  void testTheGroup(const PSet &, const edm::EventSetup &);
  void runTheTest(const PSet &, const edm::EventSetup &);

  VPSet testGroups;

  // Test methods
  // ------------
  /// Continues on error
  void check(const boost::shared_ptr<baseObjectWrapper> & inputObject, const boost::shared_ptr<baseObjectWrapper> expectedOutputObject);
  /// throws on error
  void require(const boost::shared_ptr<baseObjectWrapper> & inputObject, const boost::shared_ptr<baseObjectWrapper> expectedOutputObject);

  map<string, baseObjectWrapper * > map_;
};

#endif
