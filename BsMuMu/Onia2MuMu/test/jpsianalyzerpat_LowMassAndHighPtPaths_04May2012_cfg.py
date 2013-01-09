import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
filename = ""
if len(args) > 0: filename = args[0]


process = cms.Process("ANAPAT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_52_V7::All" #Run2012A

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
			# "file:onia2MuMuPAT.root"
      '%s' % filename, # INPUTFILE PASSED BY ARGUMENT: cmsRun cfg_name file:filename.root [or rfio:/filepath_and_filename_on_castor.root]
#			INPUTFILE
			)
)

process.hltMuF = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_DoubleMu0"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)


# filter on good vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
                                           )

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 500

#read the muon scale correction parameters from the DB
process.poolDBESSource = cms.ESSource("PoolDBESSource",
     BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
     DBParameters = cms.PSet(
          messageLevel = cms.untracked.int32(2)
#          authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
      ),
      timetype = cms.untracked.string('runnumber'),
#      connect = cms.string('oracle://cms_orcoff_prod/CMS_COND_31X_PHYSICSTOOLS'), #needs to point to the correct file
      connect = cms.string('frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS'),
      toGet = cms.VPSet(cms.PSet(
          record = cms.string('MuScleFitDBobjectRcd'),
          tag = cms.string('MuScleFit_Scale_JPsi_19_invPb_innerTrack') #adjust the tag once a new one is available
      ))
  )


process.demo = cms.EDAnalyzer('JPsiAnalyzerPAT',

    src = cms.InputTag("onia2MuMuPatTrkTrk"),
    srcWithCaloMuons = cms.InputTag("onia2MuMuPatGlbCal"),

    writeTree = cms.bool(True),
    treeFileName = cms.string("onia2MuMu_tree.root"),
    #treeFileName = cms.string("OUTPUTFILE1"),

    writeDataSet = cms.bool(True),                 
    dataSetName = cms.string("dataSet.root"),
    #dataSetName = cms.string("OUTPUTFILE2"),
    triggersForDataset = cms.vstring("HLT_Dimuon10_Jpsi_Barrel_v5"),

    massMin = cms.double(2.6),
    massMax = cms.double(3.5),
    pTBinRanges = cms.vdouble(0.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 40.0),
    etaBinRanges = cms.vdouble(0.0, 1.3, 2.5),
    onlyTheBest = cms.bool(True),		
    applyCuts = cms.bool(True),
    applyExpHitCuts = cms.untracked.bool(False),
    applyDiMuonCuts = cms.untracked.bool(False), #H:                         
    useBeamSpot = cms.bool(False),
    useCaloMuons = cms.untracked.bool(False),
    removeSignalEvents = cms.untracked.bool(False),
    removeTrueMuons = cms.untracked.bool(False),
    writeOutCandidates = cms.untracked.bool(False),
    massCorrectionMode = cms.int32(3),    # mode 0 no correction,
                                          # mode 1 constant corr,
                                          # mode 2 pt dependent corr,
                                          # mode 3 pt and eta dependent corr
    oniaPDG = cms.int32(443),
    genParticles = cms.InputTag("genMuons"),
    isMC = cms.untracked.bool(False),
    storeAllMCEvents = cms.untracked.bool(True),
    isPromptMC = cms.untracked.bool(True),
                              
    # Configuration for the extrapolation at the muon system 
    propagatorStation1 = cms.PSet(
        useStation2 = cms.bool(False), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # for AOD
        useSimpleGeometry = cms.bool(True), 
    ),
    propagatorStation2 = cms.PSet(
        useStation2 = cms.bool(True), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # for AOD
        useSimpleGeometry = cms.bool(True), 
    ),

    # Configuration of trigger matching                           
    triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
                              
    ####### 2011, 5E32 ########
                              
    HLTBitNames_SingleMu = cms.vstring(),
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_SingleMu = cms.vstring(),

    HLTBitNames_DoubleMu = cms.vstring(                                      
                                      # PARKED TRIGGERS
                                      "HLT_DoubleMu3p5_LowMassNonResonant_Displaced_v3",
                                      "HLT_DoubleMu3p5_LowMass_Displaced_v3",
                                      "HLT_Mu13_Mu8_v16" 
                                       ),
   # ONE FILTER NAME PER PATH    
   HLTLastFilterNames_DoubleMu = cms.vstring(
                                             "hltDisplacedmumuFilterDoubleMu3p5LowMassNonResonant", #HLT_DoubleMu3p5_LowMassNonResonant_Displaced_v3
                                             "hltDisplacedmumuFilterDoubleMu3p5LowMass", #HLT_DoubleMu3p5_LowMass_Displaced_v3
                                             "hltDiMuonMu13Mu8DzFiltered0p2", #HLT_Mu13_Mu8_v16
                                             ),

    HLTBitNames_MuL2Mu = cms.vstring(),
    # TWO FILTER NAMES PER PATH (FIRST is L3, SECOND is L2)                           
    HLTLastFilterNames_MuL2Mu = cms.vstring(),

    HLTBitNames_MuTrack = cms.vstring(), 
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_MuTrack = cms.vstring(),
                              
    HLTBitNames_MuTkMu = cms.vstring(),
    # ONE FILTER NAME PER PATH
    HLTLastFilterNames_MuTkMu = cms.vstring(),

)

## no filter
# process.p = cms.Path(process.demo)

## filter on vertex
process.p = cms.Path(process.primaryVertexFilter*process.demo)

#import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
myLumis = LumiList.LumiList(filename = 'Cert_190456-191859_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)
## filter on vertex and HLT
# process.p = cms.Path(process.primaryVertexFilter*process.hltMuF*process.demo)
