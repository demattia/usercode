#!/bin/env python

import os, sys
import subprocess

# if len(sys.argv) != 1:
#     print "Please, provide the dir to check"


def checkAndResubmit(crabDir):
    print "Checking status and resubmitting failed jobs for:", crabDir
    # p = subprocess.Popen(["crab", "-status -getoutput -c "+crabDir], stdout=subprocess.PIPE)
    # out, err = p.communicate()
    # print out
    os.system("crab -status -getoutput -c "+crabDir+" > monitorTmpLog.txt")

    p = subprocess.Popen(["crab", "-status", "-c", crabDir], stdout=subprocess.PIPE)
    out, err = p.communicate()
    # print out

    totalJobs = 0
    jobsWithExitCodeZero = 0
    jobsToResubmit = []
    takeNext = False
    created = False
    createdList = ""
    running = False
    runningList = ""
    cancelled = False
    # cancelledList = ""
    takeNextToNext = 0
    for line in out.split('\n'):
        if line.find("Total") != -1:
            totalJobs = int(line.split()[1])
        # print line
        if takeNextToNext == 1:
            takeNext = True
        if takeNextToNext > 0:
            takeNextToNext -= 1
        if created:
            print line
            createdList += line
            created = False
        if running:
            print line
            runningList += line
            running = False
        if cancelled:
            print line
            os.system("crab -kill "+line.split(':')[1].strip()+" -c "+crabDir)
            # cancelledList += line
            cancelled = False
        if takeNext:
            # print line
            # if jobsToResubmit != "": jobsToResubmit += ','
            # jobsToResubmit += line.split(':')[1].strip()
            jobsToResubmit += [line.split(':')[1].strip()]
            takeNext = False
        if line.find('>>>>') != -1 and line.find('Retrieved') == -1 and line.find('Submitted') == -1 and line.find('Done') == -1:
            if line.find('Created') != -1:
                print "Warning: the following jobs are in status created"
                created = True
            elif line.find('Running') != -1:
                print "The following jobs are in status running"
                running = True
            elif line.find('Cancelled') != -1:
                print "the following jobs are Cancelled and will be resubmitted"
                takeNext = True
                cancelled = True
            elif line.find('Aborted') != -1:
                print "the following jobs are Aborted and will be resubmitted"
                takeNextToNext = 2
            else:
                print line.split(':')
                exitCode = line.split(':')[1].strip()
                if exitCode != "0":
                    print "The following jobs have exitCode =", exitCode, "and will be resubmitted"
                    takeNext = True
                else:
                    jobsWithExitCodeZero = int(line.split()[1])

    if len(jobsToResubmit) == 0:
        print "No job to resubmit at this time"
    else:
        print "\nResubmitting all jobs:", jobsToResubmit

        print jobsToResubmit

    for element in jobsToResubmit:
        print "Resubmitting jobs:", element
        # p = subprocess.Popen(["crab", "-resubmit "+element+" -c "+crabDir], stdout=subprocess.PIPE)
        # out, err = p.communicate()
        # print "out: ", out
        os.system("crab -resubmit "+element+" -c "+crabDir)

    print ""
    print "Summary for "+crabDir+":"
    print "--------\n"
    print "Total number of jobs:", totalJobs
    print "Number of jobs in status done (exit code 0):", jobsWithExitCodeZero
    print "Jobs in status created:"
    print createdList
    print "Jobs in status running:"
    print runningList
    print ""
    print "Percentange of successfully completed jobs: "+str(int(float(jobsWithExitCodeZero)/float(totalJobs)*100))+"%"
    if totalJobs == jobsWithExitCodeZero:
        print "All jobs completed sucessfully!"
    # Leave an empty line at the end
    print ""


# Running
while True:
    # checkAndResubmit("workdirs/DYJets10_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/DYJets50_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDem170_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDem20_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDem250_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDem30_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDem350_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDem80_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDmu15_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/QCDmu20_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/TTJets_FullLept_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/WW_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/WZ_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/ZZ_pat_20130203/crabDir")
    # checkAndResubmit("workdirs/WJetsToLNu_pat_20130220/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH1000_MFF150_CTau10To1000_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH125_MFF50_CTau50To5000_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH1000_MFF20_CTau1p5To150_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH200_MFF20_CTau7To700_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH200_MFF50_CTau20To2000_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH400_MFF150_CTau40To4000_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH400_MFF20_CTau4To400_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH400_MFF50_CTau8To800_pat_20140211/crabDir")
    # checkAndResubmit("workdirs/Chi0ToNuLL_MSquark120_MChi48_pat_20140218/crabDir")
    # checkAndResubmit("workdirs/Chi0ToNuLL_MSquark1000_MChi148_pat_20140218/crabDir")
    # checkAndResubmit("workdirs/Chi0ToNuLL_MSquark1500_MChi494_pat_20140218/crabDir")
    # checkAndResubmit("workdirs/Chi0ToNuLL_MSquark350_MChi148_pat_20140218/crabDir")
    # checkAndResubmit("workdirs/HTo2LongLivedTo4F_MH125_MFF20_CTau13To1300_pat_20140219/crabDir")
    # checkAndResubmit("workdirs/Data_Mu_Run2012A22Jan_pat_20140218/crabDir")
    # checkAndResubmit("workdirs/Data_Mu_Run2012B22Jan_pat_20140219/crabDir")
    # checkAndResubmit("workdirs/Data_Mu_Run2012C22Jan_pat_20140219/crabDir")
    checkAndResubmit("workdirs/Data_Mu_Run2012D22Jan_pat_20140219/crabDir")
