#!/bin/env python

import os,sys

printit = True

def get_samples(mdir):

    ldir = sorted(os.listdir(mdir),reverse=True)

    thelist = []
    for dir in ldir:
        if os.path.isfile(dir):
            continue
        if dir.find("analysis")== -1:
            continue
        skip = 0
        for x in thelist:
            y = x.split("_analysis_")[0]
            z = dir.split("_analysis_")[0]
            if y == z:
                skip = 1
        if not skip:
            thelist.append(dir)

    if printit:
        print "processing sample list: ", thelist
            
    return thelist


def checkPrefilterFiles(workdirs):
    for workdir in workdirs:
        if not os.path.exists(workdir+"/prefilter.root"):
            # Try copying it from dCache
            print "ERROR: missing file",workdir+"/prefilter.root, trying to copy it from dcache"
            print "srmcp -2 \"srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/longlived/CMSSW_4_2_7/"+workdir.split("/")[-1].split("_analysis")[0]+"/prefilter.root\" \"file:///"+workdir+"/prefilter.root\""
            os.system("srmcp -2 \"srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/longlived/CMSSW_4_2_7/"+workdir.split("/")[-1].split("_analysis")[0]+"/prefilter.root\" \"file:///"+workdir+"/prefilter.root\"")
            if not os.path.exists(workdir+"/prefilter.root"):
                print "ERROR: prefilter.root copy from dCache failed"
                return -1

def mergeHistogramsFiles(workdirs):
    for dir in workdirs:
        if os.path.isfile(dir+"/histograms.root"):
            continue
        command = "cd "+dir+"; hadd histograms.root"
        listDir = os.listdir(dir)
        for fileName in listDir:
            # print fileName
            if fileName.startswith("histograms_") > 0:
                # print "found file: "+fileName
                command += " "+fileName
        if printit:
            print "merging histograms files in "+dir
        os.system(command)


def print_s(sample,name):
    print name, ":"
    for s in sample:
        print "\t", s


class Workdirs:
    workdirs_data_mu = []
    workdirs_data_e = []
    workdirs_background_mu = []
    workdirs_background_e = []
    workdirs_signal = []
    workdirs_other = []
    workdirs_benchmark_mu = []
    workdirs_benchmark_e = []

def getWorkdirs(mydir):

    samples = get_samples(mydir)
    wdir = Workdirs()

    for s in samples:
        if   s.find("Signal") != -1:
            wdir.workdirs_signal.append(s)        
        elif s.find("Data_Mu") != -1:
            wdir.workdirs_data_mu.append(s)        
        elif s.find("Data_Photon") != -1:
            wdir.workdirs_data_e.append(s)        
        elif s.find("DisplacedMu") != -1:
            wdir.workdirs_benchmark_mu.append(s)        
        elif s.find("DisplacedE") != -1:
            wdir.workdirs_benchmark_e.append(s)        
        elif s.find("mu") != -1:
            wdir.workdirs_background_mu.append(s)        
        elif (s.find("em") != -1) or (s.find("ee") != -1):
            wdir.workdirs_background_e.append(s)
        else:
            wdir.workdirs_other.append(s)

    wdir.workdirs_background_mu += wdir.workdirs_other    
    wdir.workdirs_background_e  += wdir.workdirs_other    

    if printit:
        print_s(wdir.workdirs_data_mu,"workdirs_data_mu")
        print_s(wdir.workdirs_data_e,"workdirs_data_e")
        print_s(wdir.workdirs_background_mu,"workdirs_background_mu")
        print_s(wdir.workdirs_background_e,"workdirs_background_e")
        print_s(wdir.workdirs_signal,"workdirs_signal")
        print_s(wdir.workdirs_benchmark_mu,"workdirs_benchmark_mu")
        print_s(wdir.workdirs_benchmark_e,"workdirs_benchmark_e")

    return wdir


#myDir = os.environ["LOCALRT"]+"/src/workdirs/"
#myDir = os.environ["PWD"]
#print myDir

#getWorkdirs(myDir)

#mergeHistogramsFiles(samples)
#checkPrefilterFiles(samples)

