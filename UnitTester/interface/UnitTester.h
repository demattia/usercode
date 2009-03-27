#ifndef UnitTester_h
#define UnitTester_h
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CondTools/UnitTester/interface/TestObjectWrappers.h"

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
  
  BaseObjectWrapper * BaseObjectFactory(const PSet & test, const edm::EventSetup & es ) {
    try {
      return objectFactory(test, es);
    }
    catch( const string & a ) {
      if( a == "unidentifiedType" ) {
        string type = test.getParameter<string>("Type");
        //         if( type == "bool" ) return new BoolWrapper( test.getParameter<bool>("Value") );
        //         else if( type == "double" ) return new ObjectWrapper<double>( test.getParameter<double>("Value") );
        //         else if( type == "string" ) return new ObjectWrapper<string>( test.getParameter<string>("Value") );
        cout << "Error: unrecognized type = " << type << endl;
        abort();
      }
      else throw( a );
    }
    return 0;
  }

  /**
   * This method reads the parameters and creates the ObjectWrappers. </br>
   * It needs to be implemented in the derived class.
   */
  virtual BaseObjectWrapper * objectFactory(const PSet& test, const edm::EventSetup& es) = 0;

  void initializeTheGroup(const PSet &, const edm::EventSetup &);
  void testTheGroup(const PSet &, const edm::EventSetup &);
  void runTheTest(const PSet &, const edm::EventSetup &);

  VPSet testGroups;

  // Test methods
  // ------------
  /// Continues on error
  void check(const boost::shared_ptr<BaseObjectWrapper> & inputObject, const boost::shared_ptr<BaseObjectWrapper> expectedOutputObject);
  /// throws on error
  void require(const boost::shared_ptr<BaseObjectWrapper> & inputObject, const boost::shared_ptr<BaseObjectWrapper> expectedOutputObject);

  map<string, BaseObjectWrapper * > map_;
};

#endif
