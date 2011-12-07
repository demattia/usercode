#!/usr/bin/python

import os

os.system("cmsRun MuonAnalyzerTreeWriter_cfg.py")
os.system("mv standAloneMuons.root standAloneMuons_reco.root")
os.system("mv cleanedstandAloneMuons.root cleanedStandAloneMuons_reco.root")

types=["Refit", "RefitOneIteration", "RefitOneIterationSegmentBased", "SegmentBased", "RefitSegmentBased"]
# types=["RefitOneIteration", "RefitOneIterationSegmentBased", "SegmentBased", "RefitSegmentBased"]

for type in types:
    print "Running for", type
    os.system("cat MuonAnalyzerTreeWriter_cfg.py | sed s@/Reco/@/"+type+"/@ > MuonAnalyzerTreeWriter_"+type+"_cfg.py")
    os.system("cmsRun MuonAnalyzerTreeWriter_"+type+"_cfg.py")
    os.system("mv standAloneMuons.root standAloneMuons_"+type+".root")
    os.system("mv cleanedstandAloneMuons.root cleanedStandAloneMuons_"+type+".root")

