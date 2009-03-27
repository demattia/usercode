import FWCore.ParameterSet.Config as cms    

def BuildTestGroup(label, parameters, tests, fixtures = []) :
    return cms.PSet(
        Label = cms.string(label),
        Parameters = cms.string(parameters),
        Tests = cms.VPSet() + tests,
        Fixtures = cms.VPSet() + fixtures
        )

def BuildTest(label, inputFixture, inputParameters, expectedOutputFixture, expectedOutputParameters, testType) :
    return cms.PSet(
        Label = cms.string(label),
        InputFixture = cms.string(inputFixture),
        InputParameters = inputParameters,
        ExpectedOutputFixture = cms.string(expectedOutputFixture),
        ExpectedOutputParameters = expectedOutputParameters,
        TestType = cms.string(testType)
        )

def BuildTestElement(label, parameters) :
    return cms.PSet(
        Label = cms.string(label),
        Parameters = parameters
        )

def BuildFixture(label, fixtureName, parameters) :
    return cms.PSet(
        Label = cms.string(label),
        FixtureName = cms.string(fixtureName),
        Parameters = parameters
        )

def BuildBool(bool) :
    return cms.PSet(
        Label = cms.string("bool"),
        Value = cms.bool(bool)
        )
