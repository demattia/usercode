import FWCore.ParameterSet.Config as cms    

def BuildTestGroup(label, parameters, tests, fixtures) :
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

# ----------------------------------------------------- #
# End of fixed functions to be moved to a separate file #
# ----------------------------------------------------- #

def BuildClusterSet(detId, cluster) :
    return cms.PSet(
        DetId = cms.uint32(detId),
        ClusterSet = cms.VPSet() + cluster
        )

def BuildCluster(firstStrip, amplitudes) :
    return cms.PSet(
        FirstStrip = cms.uint32(firstStrip),
        Amplitudes = cms.vuint32(amplitudes),
        )

def BuildAlgo(label, firstStrip, amplitudes) :
    return cms.PSet(
        Label = cms.string(label),
        Parameters = cms.PSet(),
        )

Tests = [
    BuildTest( "First Test",
               "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(0, [0])])),
               "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(0, [0])])),
               "Require"
               ),
    BuildTest( "Second Test",
               "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(0, [0])])),
               "", BuildTestElement("Cluster", BuildClusterSet(1, [BuildCluster(0, [0])])),
               "Require"
               ),
    ]


AllFixtures = []
#AllFixtures = [
#    BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(0, [0])]))
#    ]

TestGroup = BuildTestGroup( "First Group of Test",
                            "test1",
                            Tests,
                            AllFixtures
                            )
