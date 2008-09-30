#!/bin/csh
cd cmsswDir
eval `scramv1 runtime -csh`
cd -
cmsRun cfgFile
rfcp outFileName outDir
