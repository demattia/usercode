#!/bin/bash

# Muon Specific
###############

# QCD mu
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu15_analysis_20120831/histograms.root", 2738580.0, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu20_analysis_20120831/histograms.root", 1865500.0, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu30_analysis_20120831/histograms.root", 806298.0, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu50_analysis_20120831/histograms.root", 176187.6, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu80_analysis_20120831/histograms.root", 40448.0, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu120_analysis_20120831/histograms.root", 7463.94, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu170_analysis_20120831/histograms.root", 2299.752, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu300_analysis_20120831/histograms.root", 151.8048, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu470_analysis_20120831/histograms.root", 11.79648, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu600_analysis_20120831/histograms.root", 2.690196, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu800_analysis_20120831/histograms.root", 0.368781, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDmu1000_analysis_20120831/histograms.root", 0.0849078, false)'

# TTbar
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/TTJets_analysis_20120831/histograms.root", 225.2, false)'

# DY JETS
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/DYJets10_analysis_20120831/histograms.root", 915, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/DYJets50_analysis_20120831/histograms.root", 3503.7, false)'

# Di-boson
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/ZZ_analysis_20120831/histograms.root", 8.3, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/WZ_analysis_20120831/histograms.root", 22, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/WW_analysis_20120831/histograms.root", 57.1, false)'

# Data
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Mu_Run2012A1_analysis_20120831/histograms.root", 1, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Mu_Run2012B1_analysis_20120831/histograms.root", 1, false)'
# root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Mu_Run2012C1_analysis_20120831/histograms.root", 1, false)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Mu_Run2012C2_analysis_20120903/histograms.root", 1, false)'

# Electron specific
###################

# QCD em
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDem20_analysis_20120831/histograms.root", 2914860.0, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDem30_analysis_20120831/histograms.root", 4615893.0, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDem80_analysis_20120831/histograms.root", 183294.9, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDem170_analysis_20120831/histograms.root", 4586.52, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDem250_analysis_20120831/histograms.root", 556.75, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/QCDem350_analysis_20120831/histograms.root", 89.1, true)'

# Di-boson
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/ZZ_analysis_20120831/histograms.root", 8.3, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/WZ_analysis_20120831/histograms.root", 22, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/WW_analysis_20120831/histograms.root", 57.1, true)'

# TTbar
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/TTJets_analysis_20120831/histograms.root", 225.2, true)'

# DY JETS
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/DYJets10_analysis_20120831/histograms.root", 915, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/DYJets50_analysis_20120831/histograms.root", 3503.7, true)'

# Data
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Photon_Run2012A1_analysis_20120831/histograms.root", 1, true)'
root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Photon_Run2012B1_analysis_20120831/histograms.root", 1, true)'
# root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Photon_Run2012C1_analysis_20120831/histograms.root", 1, true)'
# root -l -b -q 'treeAnalyzer.C+("/home/cmsdas/store/user/cmsdas/2012/ExoticaDisplacedVertices/2012/workdirs/Data_Photon_Run2012C2_analysis_20120831/histograms.root", 1, true)'
