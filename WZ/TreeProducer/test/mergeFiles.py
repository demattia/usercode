#!/usr/bin/env python

import os

command = "hadd histograms.root"
for file in os.listdir("./"):
    if file.find("histograms") != -1:
        command += " "+file

os.system(command)
# os.system("mv histograms.root ../")
# os.system("rm *")
# os.system("mv ../histograms.root .")

# print command
