#!/usr/bin/env python

import os

for dir in os.listdir("./"):
    if os.path.isdir(dir):
        print dir
        os.system("cp mergeFiles.py "+dir)
        os.system("cd "+dir+"; ./mergeFiles.py")
        os.system("cd "+dir+"; rm condor*")
        os.system("cd "+dir+"; rm job*")
        os.system("cd "+dir+"; rm run*")
        # os.system("cd "+dir+"; rm histograms_*")
