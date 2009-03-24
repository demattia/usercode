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

def BuildFixture(label, fixtureName, parameters) :
    return cms.PSet(
        Label = cms.string(label),
        FixtureName = cms.string(fixtureName),
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

def BuildOldAlgo(algorithm, channelThreshold, seedThreshold, clusterThreshold, maxSequentialHoles, qualityLabel) :
    #return cms.PSet(
    return cms.PSet( Algorithm = cms.string(algorithm),
                     ChannelThreshold = cms.double(channelThreshold),
                     SeedThreshold    = cms.double(seedThreshold),
                     ClusterThreshold = cms.double(clusterThreshold),
                     MaxSequentialHoles = cms.uint32(maxSequentialHoles),
                     QualityLabel = cms.string(qualityLabel)
                     # ),
                     )

def BuildNewAlgo(algorithm, channelThreshold, seedThreshold,
                 clusterThreshold, maxSequentialHoles, maxSequentialBad,
                 maxAdjacentBad, qualityLabel) :
    return cms.PSet( Algorithm = cms.string(algorithm),
                     ChannelThreshold = cms.double(channelThreshold),
                     SeedThreshold = cms.double(seedThreshold),
                     ClusterThreshold = cms.double(clusterThreshold),
                     MaxSequentialHoles = cms.uint32(maxSequentialHoles),
                     MaxSequentialBad = cms.uint32(maxSequentialBad),
                     MaxAdjacentBad = cms.uint32(maxAdjacentBad),
                     QualityLabel = cms.string(qualityLabel)
                     )

def BuildDigiSet(detId, digi) :
    return cms.PSet(
        DetId = cms.uint32(detId),
        DigiSet = cms.VPSet() + digi
        )

def BuildDigi(strip, adc, noise, gain, quality) :
    return cms.PSet( Strip = cms.uint32(strip),
                     ADC   = cms.uint32(adc),
                     Noise = cms.uint32(noise),
                     Gain  = cms.uint32(gain),
                     Quality = cms.uint32(quality)
                     )

# Tests are made of: Label, InputFixture, InputParameters, ExpectedOutputFixture, ExpectedOutputParameters, TestType
# The parameters should be built with the BuildTestElement which requires a Label and a PSet.
def BuildAllTests(fixtureName) :
    return [
        BuildTest( "(sanity check) : Cluster = Cluster",
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(0, [0])])),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(0, [0])])),
                   "Require"
                   ),
        BuildTest( "[] = []",
                   fixtureName, BuildDigiSet(0, []),
                   "", BuildTestElement("Cluster", BuildClusterSet(0,[])),
                   "Require"
                   ),
        BuildTest( "(4/1) = []",
                   # In this case the label is specified as the fixture label, so it is sufficient to create a PSet.
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 4, 1, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [])),
                   "Require"
                   ),
        BuildTest( "(5/1) = [5]",
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 5, 1, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(10, [5])])),
                   "Check"
                   ),
        BuildTest( "(110/1) = [110]",
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 110, 1, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(10, [110])])),
                   "Require"
                   ),
        BuildTest( "(24/5) = []",
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 24, 5, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [])),
                   "Require"
                   ),
        BuildTest( "(24/5) = [5]",
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 24, 5, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(10, [24])])),
                   "Require"
                   ),
        BuildTest( "(25/5) = [5]",
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 25, 5, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(10, [25])])),
                   "Require"
                   ),
        BuildTest( "(111/5) = [111]",
                   fixtureName, BuildDigiSet(0, [BuildDigi(10, 111, 5, 1, True)]),
                   "", BuildTestElement("Cluster", BuildClusterSet(0, [BuildCluster(10, [111])])),
                   "Require"
                   ),
        ]

# Fixtures are made of a testElement. They should be built with the BuildTestElement function which requires a Label
# and a PSet.
OldAlgoFixtures = [
    BuildFixture("Algo", "OldThreeThresholdAlgorithm", BuildOldAlgo( "OldThreeThresholdAlgorithm", 2, 3, 5, 0, "" ) ),
    ]

NewAlgoFixtures = [
    BuildFixture("Algo", "ThreeThresholdAlgorithm", BuildNewAlgo( "ThreeThresholdAlgorithm", 2, 3, 5, 0, 0, 1, "" ) )
    ]

OldAlgoTests = BuildAllTests( "OldThreeThresholdAlgorithm" )
NewAlgoTests = BuildAllTests( "ThreeThresholdAlgorithm" )


OldAlgoTestGroup = BuildTestGroup( "Old Algo group of Test",
                            "test1",
                            OldAlgoTests,
                            OldAlgoFixtures
                            )

NewAlgoTestGroup = BuildTestGroup( "New Algo group of Test",
                            "test2",
                            NewAlgoTests,
                            NewAlgoFixtures
                            )
