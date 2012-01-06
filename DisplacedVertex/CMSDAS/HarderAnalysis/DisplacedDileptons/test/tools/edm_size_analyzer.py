#!/bin/env python
import os,sys

#
# check file name argument
#
if len(sys.argv)!=2:
    print "usage:",sys.argv[0],"edmFileName.root"
    sys.exit(1)
edmfile=sys.argv[1]
if not os.path.exists(edmfile):
    print "file not found"
    sys.exit(1)
fullsize=os.path.getsize(edmfile)
print "original file size:",fullsize

#
# loop over all collections and drop each one in turn
#
dummylabel="dummy_dummy_dummy_dummy"
result=os.popen("edmDumpEventContent --name "+edmfile).readlines()
result=[dummylabel]+result
nprocessed=0
size_reference=0
for line in result:
    branchname=line.strip("\n").replace("__","_*_")

    # create cmsRun config file that drops this particular collection
    configfile=open("edm_size_analyzer_cfg.py","w")
    configfile.write("import FWCore.ParameterSet.Config as cms\n")
    configfile.write("process = cms.Process(\"Test\")\n")
    configfile.write("process.source = cms.Source(\"PoolSource\",\n")
    configfile.write("  duplicateCheckMode = cms.untracked.string(\"noDuplicateCheck\"),\n")
    configfile.write("  fileNames = cms.untracked.vstring(\"file:"+edmfile+"\"))\n")
    configfile.write("process.output = cms.OutputModule(\"PoolOutputModule\",\n")
    configfile.write("  outputCommands = cms.untracked.vstring(\"keep *\",\n")
    configfile.write("     \"drop "+branchname+"\"),\n")
    configfile.write("  fileName = cms.untracked.string(\"edm_size_analyzer.root\"))\n")
    configfile.write("process.outpath = cms.EndPath(process.output)\n")
    configfile.close()

    # run cmsRun
    if os.path.exists("edm_size_analyzer.root"):
        os.remove("edm_size_analyzer.root")
    os.system("cmsRun edm_size_analyzer_cfg.py &> /dev/null")
    if not os.path.exists("edm_size_analyzer.root"):
        print "cmsRun did not create the output root file"
        sys.exit(4)
    newsize=os.path.getsize("edm_size_analyzer.root")
    if branchname==dummylabel:
        size_reference=newsize
    else:
        spacer=""
        for i in range(60-len(branchname)): spacer+=" "
        print branchname+spacer+": %3i%%"%int((1.005-float(newsize)\
                                               /float(size_reference))*100.)
    
    os.remove("edm_size_analyzer_cfg.py")
    os.remove("edm_size_analyzer.root")
    nprocessed+=1



if nprocessed!=len(result):
    print "ERROR: only",nprocessed,"out of",len(result),"collections dealt with!"
    sys.exit(3)

