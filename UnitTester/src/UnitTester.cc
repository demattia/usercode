#include "CondTools/UnitTester/interface/UnitTester.h"

#include <functional>
#include <numeric>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

using namespace std;

void UnitTester::analyze(const edm::Event &, const edm::EventSetup & es)
{
  for(iter_t group = testGroups.begin(); group < testGroups.end(); ++group) {
    initializeTheGroup(*group, es);
    testTheGroup(*group, es);
  }
}

void UnitTester::initializeTheGroup(const PSet& group, const edm::EventSetup& es)
{
  // Create any fixtures
  VPSet fixtures = group.getParameter<VPSet>("Fixtures");
  for(iter_t fixture = fixtures.begin();  fixture < fixtures.end(); ++fixture) {
    // Use the internal label as name of the fixture as the external one is the type of the object.
    // cout << "fixture->getParameter<string>(\"Label\") = " << fixture->getParameter<string>("FixtureName") << endl;
    map_.insert( make_pair(fixture->getParameter<string>("FixtureName"), objectFactory(*fixture, es)) );
  }
}

void UnitTester::testTheGroup(const PSet & group, const edm::EventSetup & es)
{
  std::string label = group.getParameter<std::string>("Label");
  // PSet params = group.getParameter<PSet>("Parameters");
  VPSet tests = group.getParameter<VPSet>("Tests");

  std::cout << "\nTesting group: \"" << label << "\"\n               " << endl;
  // << params << std::endl;
  for(iter_t test = tests.begin();  test < tests.end();  test++) {
    runTheTest(*test, es);
  }
}

void UnitTester::runTheTest(const PSet& test, const edm::EventSetup & es)
{
  std::string label =  test.getParameter<std::string>("Label");

  boost::shared_ptr<BaseObjectWrapper> inputObject;
  boost::shared_ptr<BaseObjectWrapper> expectedOutputObject;

  // Get it from a map of fixtures
  string inputFixture = test.getParameter<string>("InputFixture");
  PSet inputParameters = test.getParameter<PSet>("InputParameters");
  if( inputFixture != "" ) {
    inputObject.reset(map_.find(inputFixture)->second->eval(inputParameters));
  }
  else inputObject.reset(objectFactory(inputParameters, es));

  string expectedOutputFixture = test.getParameter<string>("ExpectedOutputFixture");
  PSet expectedOutputParameters = test.getParameter<PSet>("ExpectedOutputParameters");
  // cout << "test.getParameter<string>(\"ExpectedOutputFixture\") = " << test.getParameter<string>("ExpectedOutputFixture") << endl;
  if( test.getParameter<string>("ExpectedOutputFixture") != "" ) {
    // cout << "before expectedOutput: inside if" << endl;
    expectedOutputObject.reset(map_.find(expectedOutputFixture)->second->eval(expectedOutputParameters));
  }
  else {
    // cout << "before expectedOutput: inside else" << endl;
    expectedOutputObject.reset(objectFactory(expectedOutputParameters, es));
  }
  // std::cout << "Testing: \"" << label << "\"\n";
  std::cout << "Testing: \"" << label << "\"\t\t ";

  // Select the type of tests based on input parameters
  string testType = test.getParameter<string>("TestType");
  if( testType == "Check" ) check(inputObject, expectedOutputObject);
  else if ( testType == "Require" ) {
    try {
      require(inputObject, expectedOutputObject);
    }
    catch(string a) {
      if( a == "Error" ) cout << "Require error" << endl;
      else throw(a);
    }
  }
  else {
    cout << "Unrecognized test type = " << testType << endl;
    abort();
  }
}

void UnitTester::check(const boost::shared_ptr<BaseObjectWrapper> & inputObject, const boost::shared_ptr<BaseObjectWrapper> expectedOutputObject)
{
  if( *inputObject != *expectedOutputObject ) cout << "Error: check failed" << endl;
  else cout << endl;
}

void UnitTester::require(const boost::shared_ptr<BaseObjectWrapper> & inputObject, const boost::shared_ptr<BaseObjectWrapper> expectedOutputObject)
{
  if( *inputObject != *expectedOutputObject ) {
    cout << "Error: check failed" << endl;
    throw(string("Error"));
  }
  else cout << endl;
}
