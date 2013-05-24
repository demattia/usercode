#!/usr/env python

import os
import sys

dirName = sys.argv[1]
outDirName = sys.argv[2]
print dirName

def all_files(directory):
    for path, dirs, files in os.walk(directory):
        for f in files:
            yield os.path.join(path, f)
            
root_files = [f for f in all_files(dirName)
             if f.endswith('.root')]


outputFile = open("makeSelectionBs_cfg.py", "w")

templateFile = open("makeSelectionBs_template_cfg.py")
for line in templateFile:
    if line.find("FILELIST") != -1:
        i = 1
        for file in root_files:
            lineToSave = file.replace('/pnfs/cms/WAX/11', '"').replace('.root', '.root",\n')
            if i < len(root_files):
                outputFile.write(lineToSave)
            else:
                outputFile.write(lineToSave.replace('",\n', '"\n'))
            i += 1
    else:
        outputFile.write(line)
outputFile.close()

os.system('mkdir -p '+outDirName+'; mv makeSelectionBs_cfg.py '+outDirName)
os.system('cp condor '+outDirName+'; cp Run.csh '+outDirName)

