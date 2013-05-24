#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
cd /uscms_data/d3/demattia/BsMuMu/CMSSW_5_3_8/src/
cmsenv
cd -
cmsRun makeSelectionBs_cfg.py
