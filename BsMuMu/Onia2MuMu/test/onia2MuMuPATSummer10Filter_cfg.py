# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms

process = cms.Process("Onia2MuMuPAT")

from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import *

onia2MuMuPAT(process, GlobalTag="START38_V8::All", MC=True, HLT="REDIGI36X", Filter=True)

process.source.fileNames = cms.untracked.vstring(
       '/store/mc/Summer10/BpToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/AODSIM/START36_V9_S09-v1/0020/945DD77B-C181-DF11-B02B-003048D15D04.root',
       '/store/mc/Summer10/BpToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/AODSIM/START36_V9_S09-v1/0017/5803E159-A781-DF11-9DE5-00304867BFAA.root',
       '/store/mc/Summer10/BpToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/AODSIM/START36_V9_S09-v1/0014/629C78F8-9B81-DF11-9882-001BFCDBD176.root',
       '/store/mc/Summer10/BpToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/AODSIM/START36_V9_S09-v1/0014/18C906C4-9B81-DF11-81B7-002354EF3BDF.root',
       '/store/mc/Summer10/BpToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/AODSIM/START36_V9_S09-v1/0011/06DC9018-8081-DF11-9D24-002354EF3BE6.root',
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

