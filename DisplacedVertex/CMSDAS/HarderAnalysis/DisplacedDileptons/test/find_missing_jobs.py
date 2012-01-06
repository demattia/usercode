#!/bin/env python

import os,sys,ROOT

workdir=sys.argv[1]
print "cd "+workdir
filenames=os.listdir(workdir)
for filename in filenames:
    try:
        jobnum=int(filename.replace("job","").replace(".sh",""))
    except:
        continue
    file_ok=1
    if "histograms_%i.root"%jobnum in filenames:
        # output file is there, but need to check whether it is corrupted
        histfile=ROOT.TFile.Open(workdir+"/histograms_%i.root"%jobnum)
        try:
            num_events=histfile.Get("prefilterPassTwo/numSignal").GetEntries()
        except:
            file_ok=0
            pass
        pass
    else:
        # file is missing alltogether
        file_ok=0
        pass

    if not file_ok:
        print "qsub -q prod "+filename
        pass
    pass
print "cd $CMSSW_BASE/src"
       
