import FWCore.ParameterSet.Config as cms

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tagAndProbe"),
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "2.8", "3.5", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        phi = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        pair_drM2 = cms.vstring("dR in 2nd MS", "-1e6", "1e6", ""),
        run = cms.vstring("Run", "0", "1000000", ""),
        L1dR = cms.vstring("L1dR", "-1e6", "1e6", ""),
        tag_L1dR = cms.vstring("tag_L1dR", "-1e6", "1e6", ""),
    ),
    Categories = cms.PSet(
        # Muon selectors
        TMLSAT = cms.vstring("TMLSAT", "dummy[pass=1,fail=0]"),
        TMA    = cms.vstring("TMA", "dummy[pass=1,fail=0]"),
        TM     = cms.vstring("TM", "dummy[pass=1,fail=0]"),
        TM2    = cms.vstring("TM2", "dummy[pass=1,fail=0]"),                
        TMOST  = cms.vstring("TMOST", "dummy[pass=1,fail=0]"),
        Glb    = cms.vstring("Glb", "dummy[pass=1,fail=0]"),
        GlbPT  = cms.vstring("GlbPT", "dummy[pass=1,fail=0]"),
        Calo   = cms.vstring("Calo", "dummy[pass=1,fail=0]"),

        # Track Quality Cuts
        TQ = cms.vstring("TQ", "dummy[pass=1,fail=0]"),
        TQ1= cms.vstring("TQ1", "dummy[pass=1,fail=0]"),

        # Trigger for the tag muon
        tag_Mu5_Track0_Jpsi_MU = cms.vstring("tag_Mu3_Track0_Jpsi_MU", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track3_Jpsi_MU = cms.vstring("tag_Mu3_Track3_Jpsi_MU", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track5_Jpsi_MU = cms.vstring("tag_Mu3_Track5_Jpsi_MU", "dummy[pass=1,fail=0]"),
        tag_MuX_TrackY_Jpsi_MU = cms.vstring("tag_Mu3_OR_Mu5_Track0_Jpsi_MU", "dummy[pass=1,fail=0]"),
        tag_Mu5L2Mu0 = cms.vstring("tag_Mu5L2Mu0", "dummy[pass=1,fail=0]"),
        tag_Trigger                   = cms.vstring("tag_Trigger", "dummy[pass=1,fail=0]"),
        tag_Mu3                   = cms.vstring("tag_Mu3", "dummy[pass=1,fail=0]"),                
        tag_L1DoubleMuOpen_Tight = cms.vstring("tag_L1DoubleMuOpen_Tight", "dummy[pass=1,fail=0]"),
        tag_L1DoubleMuOpen = cms.vstring("tag_L1DoubleMuOpen", "dummy[pass=1,fail=0]"),
        tag_L1MuOpen = cms.vstring("tag_L1MuOpen", "dummy[pass=1,fail=0]"),
        tag_L2DoubleMu0 = cms.vstring("tag_L2DoubleMu0", "dummy[pass=1,fail=0]"),

        # Trigger for the probes
        L1SingleMuOpen = cms.vstring("L1SingleMuOpen", "dummy[pass=1,fail=0]"),
        L1SingleMu0 = cms.vstring("L1SingleMu0", "dummy[pass=1,fail=0]"),
        L1SingleMu3 = cms.vstring("L1SingleMu3", "dummy[pass=1,fail=0]"),
        L1SingleMu5 = cms.vstring("L1SingleMu5", "dummy[pass=1,fail=0]"),
        L1SingleMu7 = cms.vstring("L1SingleMu7", "dummy[pass=1,fail=0]"),
        L1SingleMu20 = cms.vstring("L1SingleMu20", "dummy[pass=1,fail=0]"),
        L1DoubleMuOpen       = cms.vstring("L1DoubleMuOpen", "dummy[pass=1,fail=0]"),
        L1DoubleMuOpen_Tight = cms.vstring("L1DoubleMuOpen_Tight", "dummy[pass=1,fail=0]"),
        L1DoubleMuOpen_OR    = cms.vstring("L1DoubleMuOpen_OR", "dummy[pass=1,fail=0]"),                                
        L2DoubleMu0       = cms.vstring("L2DoubleMu0", "dummy[pass=1,fail=0]"),
        DoubleMu0         = cms.vstring("DoubleMu0", "dummy[pass=1,fail=0]"),
        DoubleMu0Quarkonium = cms.vstring("DoubleMu0Quarkonium", "dummy[pass=1,fail=0]"),
        DoubleMu3         = cms.vstring("DoubleMu3", "dummy[pass=1,fail=0]"),
        Mu0         = cms.vstring("Mu0", "dummy[pass=1,fail=0]"),
        pair_badCowboy = cms.vstring("badCowboy", "dummy[pass=1,fail=0]"),
    ),
    Cuts = cms.PSet(
        matched = cms.vstring("Matched", "L1dR", "0.3"),
    ),
    ## The PDF is a Gaussian  (peak) + Exponential (background) 
    PDFs = cms.PSet(
        gaussPlusExpo = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0.,1.]",
            "signalFractionInPassing[0.9]"
        )
    ),
    Efficiencies = cms.PSet(),
    binsForMassPlots = cms.uint32(35),
)

